//
//  GLVariableOperations.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 4/21/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import "GLVariableOperations.h"
#import "GLDimension.h"
#import "GLVariable.h"
#import "GLMemoryPool.h"

/************************************************/
/*		GLVariableOperation						*/
/************************************************/

@interface GLVariableOperation ()

+ (void) setVariableOperationPrototype: (variableOperation) opPrototype forSubclass: (Class) subclass;
+ (variableOperation) variableOperationPrototypeForSubclass: (Class) subclass;

+ (void) setOperandPrototype: (NSArray *) operandPrototype forSubclass: (Class) subclass;
+ (NSArray *) operandPrototypeForSubclass: (Class) subclass;

+ (void) setResultPrototype: (NSArray *) resultPrototype forSubclass: (Class) subclass;
+ (NSArray *) resultrototypeForSubclass: (Class) subclass;

+ (void) setBufferLengthsPrototype: (NSArray *) bufferLengthsPrototype forSubclass: (Class) subclass;
+ (NSArray *) bufferLengthsPrototypeForSubclass: (Class) subclass;

+ (void) setPreOperationPrototype: (dispatch_block_t) preOpPrototype forSubclass: (Class) subclass;
+ (dispatch_block_t) preOperationPrototypeForSubclass: (Class) subclass;

+ (void) setPostOperationPrototype: (dispatch_block_t) postOperationPrototype forSubclass: (Class) subclass;
+ (dispatch_block_t) postOperationPrototypeForSubclass: (Class) subclass;

// To be called after the result and operands are set.
- (void) setupDependencies;

// To be called after the operation is complete.
- (void) tearDownDependencies;

@end

@implementation GLVariableOperation : NSOperation

static NSMapTable *classVariableOperationTable = nil;
static NSMapTable *classOperandTable = nil;
static NSMapTable *classResultTable = nil;
static NSMapTable *classBufferLengthTable = nil;
static NSMapTable *classPreOperationTable = nil;
static NSMapTable *classPostOperationTable = nil;

- (id) init
{
    if ((self=[super init])) {
        self.graphvisDescription = @"unnamed";
    }
    return self;
}

@synthesize isEnqueued;

- (BOOL) canOperateInPlace {
	return NO;
}

- (BOOL) isEqualToOperation: (id) otherOperation {
    if ( ![[self class] isSubclassOfClass: [otherOperation class]] )  {
        return NO;
    }
	
	GLVariableOperation * op = otherOperation;
	if (self.result.count != op.result.count) {
        return NO;
    } else if (self.operand.count != op.operand.count) {
        return NO;
    } else if (self.buffer.count != op.buffer.count) {
        return NO;
    }
	
	for ( NSUInteger i=0; i<self.result.count; i++) {
        if ( self.result[i] != op.result[i]) {
            return NO;
        }
    }
	
	for ( NSUInteger i=0; i<self.operand.count; i++) {
        if ( self.operand[i] != op.operand[i]) {
            return NO;
        }
    }
	
    return YES;
}

- (id) initWithResult: (NSArray *) result operand: (NSArray *) operand buffers: (NSArray *) buffers operation: (variableOperation) op
{
	if (( self = [super init] ))
	{
		self.result = result;
		self.operand = operand;
		self.buffer = buffers;
		self.operation = op;
		
		[self setupDependencies];
	}
	
    return self;
}

- (id) initWithResult: (NSArray *) result operand: (NSArray *) operand
{
	if (( self = [super init] ))
	{
		self.result = result;
		self.operand = operand;
		
		[self setupDependencies];
	}
	
    return self;
}

- (id) initWithOperand: (NSArray *) operand
{
	if (( self = [super init] ))
	{
		NSMutableArray *array = [[NSMutableArray alloc] initWithCapacity: operand.count];
		for (GLVariable *variable in operand) {
			[array addObject: [[variable class] variableOfType: variable.dataFormat withDimensions: variable.dimensions forEquation: variable.equation]];
		}
		self.result = array;
		self.operand = operand;
		
		[self setupDependencies];
	}
	
    return self;
}

- (void) main
{
	if (self.preOperation) {
		self.preOperation();
	}
	
	NSMutableArray *resultBuffer = [[NSMutableArray alloc] initWithCapacity: self.result.count];
	NSMutableArray *operandBuffer = [[NSMutableArray alloc] initWithCapacity: self.operand.count];
	for (GLVariable *variable in self.result) {
		[resultBuffer addObject: variable.data];
	}
	for (GLVariable *variable in self.operand) {
		[operandBuffer addObject: variable.data];
	}
	
	if (self.buffer.count) {
		NSMutableArray *dataBuffers = [[NSMutableArray alloc] initWithCapacity: self.buffer.count];
		for (GLBuffer *aBuffer in self.buffer) {
			[dataBuffers addObject: [[GLMemoryPool sharedMemoryPool] dataWithLength: aBuffer.numBytes]];
		}
		self.operation( resultBuffer, operandBuffer, dataBuffers );
		for (NSMutableData *data in dataBuffers) {
			[[GLMemoryPool sharedMemoryPool] returnData: data];
		}
	} else {
		self.operation( resultBuffer, operandBuffer, nil );
	}
	
	if (self.postOperation) {
		self.postOperation();
	}
	
	[self tearDownDependencies];
}

- (void) setupDependencies
{
	for (GLVariable *variable in self.operand) {
		if (variable.lastOperation && ![self.dependencies containsObject:variable.lastOperation]) {
			[self addDependency: variable.lastOperation];
		}
	}
	for (GLVariable *variable in self.result) {
		if (variable.lastOperation && ![self.dependencies containsObject:variable.lastOperation]) {
			[self addDependency: variable.lastOperation];
		}
	}
	for (GLVariable *variable in self.result) {
		[variable addOperation: self];
	}
}

- (void) tearDownDependencies
{
	for (GLVariable *variable in self.result) {
		[variable removeOperation: self];
	}
	self.result = nil;
	self.operand = nil;
}

@end













