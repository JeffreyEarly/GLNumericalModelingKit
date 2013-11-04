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

/************************************************/
/*		GLVariableOperation						*/
/************************************************/

@interface GLVariableOperation ()

// To be called after the result and operands are set.
- (void) setupDependencies;

// To be called after the operation is complete.
- (void) tearDownDependencies;

@end

@implementation GLVariableOperation : NSOperation

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
    } else if (self.buffers.count != op.buffers.count) {
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
		self.buffers = buffers;
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
	NSMutableArray *resultBuffer = [[NSMutableArray alloc] initWithCapacity: self.result.count];
	NSMutableArray *operandBuffer = [[NSMutableArray alloc] initWithCapacity: self.operand.count];
	for (GLVariable *variable in self.result) {
		[resultBuffer addObject: variable.data];
	}
	for (GLVariable *variable in self.operand) {
		[operandBuffer addObject: variable.data];
	}
	self.operation( resultBuffer, operandBuffer, self.buffers );
	
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













