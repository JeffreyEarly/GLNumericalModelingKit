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
#import "GLOperationOptimizer.h"
#import "GLLinearTransform.h"
#import <objc/runtime.h>

/************************************************/
/*		GLOptimizedVariableOperation Interface	*/
/************************************************/

#pragma mark -
#pragma mark GLOptimizedVariableOperation Interface
#pragma mark

@interface GLOptimizedVariableOperation : GLVariableOperation

+ (NSUInteger) totalSubclasses;

+ (void) setVariableOperationPrototype: (variableOperation) opPrototype forSubclassWithName: (NSString *) subclassName;
+ (variableOperation) variableOperationPrototypeForSubclassWithName: (NSString *) subclassName;

+ (void) setOperandPrototype: (NSArray *) operandPrototype forSubclassWithName: (NSString *) subclassName;
+ (NSArray *) operandPrototypeForSubclassWithName: (NSString *) subclassName;

+ (void) setResultPrototype: (NSArray *) resultPrototype forSubclassWithName: (NSString *) subclassName;
+ (NSArray *) resultPrototypeForSubclassWithName: (NSString *) subclassName;

+ (void) setBufferPrototype: (NSArray *) bufferLengthsPrototype forSubclassWithName: (NSString *) subclassName;
+ (NSArray *) bufferPrototypeForSubclassWithName: (NSString *) subclassName;

+ (void) setPreOperationPrototype: (dispatch_block_t) preOpPrototype forSubclassWithName: (NSString *) subclassName;
+ (dispatch_block_t) preOperationPrototypeForSubclassWithName: (NSString *) subclassName;

+ (void) setPostOperationPrototype: (dispatch_block_t) postOperationPrototype forSubclassWithName: (NSString *) subclassName;
+ (dispatch_block_t) postOperationPrototypeForSubclassWithName: (NSString *) subclassName;

@end

/************************************************/
/*		GLVariableOperation						*/
/************************************************/

#pragma mark -
#pragma mark GLVariableOperation
#pragma mark

@interface GLVariableOperation ()

// To be called after the result and operands are set.
- (void) setupDependencies;

// To be called after the operation is complete.
- (void) tearDownDependencies;

@end

@implementation GLVariableOperation : NSOperation

+ (Class) variableOperationSubclassWithOperand: (NSArray *) operandVariables result: (NSArray *) resultVariables
{
    GLOperationOptimizer *optimizer = [[GLOperationOptimizer alloc] initWithTopVariables: operandVariables bottomVariables: resultVariables];
    if (!optimizer.operationBlock) return nil;
    
    NSUInteger suffix = [GLOptimizedVariableOperation totalSubclasses];
    NSString *subclassName = [NSString stringWithFormat: @"GLOptimizedVariableOperation_%2lu",suffix ];
    Class subclass = NSClassFromString(subclassName);
    if (subclass) {
        NSLog(@"Whoops! A subclass already exists with this name: %@", subclassName);
    } else {
        subclass = objc_allocateClassPair([GLOptimizedVariableOperation class], [subclassName UTF8String], 0);
        objc_registerClassPair(subclass);
    }
    
    [GLOptimizedVariableOperation setVariableOperationPrototype: optimizer.operationBlock forSubclassWithName: subclassName];
    [GLOptimizedVariableOperation setOperandPrototype: optimizer.operandVariablePrototypes forSubclassWithName: subclassName];
    [GLOptimizedVariableOperation setResultPrototype: optimizer.resultVariablePrototypes forSubclassWithName: subclassName];
    [GLOptimizedVariableOperation setBufferPrototype: optimizer.internalDataBuffers forSubclassWithName: subclassName];
    
    return subclass;
}

/************************************************/
/*		Primary Interface						*/
/************************************************/



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
		for (GLTensor *variable in operand) {
			if (variable.rank == 0) {
				GLScalar *scalar = (GLScalar *) variable;
				[array addObject: [[GLScalar alloc] initWithType: scalar.dataFormat forEquation:scalar.equation]];
			} else if (variable.rank == 1) {
				GLVariable *function = (GLVariable *) variable;
				[array addObject: [[function class] variableOfType: function.dataFormat withDimensions: function.dimensions forEquation: function.equation]];
			}  else if (variable.rank == 2) {
				GLLinearTransform *matrix = (GLLinearTransform *) variable;
				[array addObject: [GLLinearTransform transformOfType: matrix.dataFormat withFromDimensions: matrix.fromDimensions toDimensions: matrix.toDimensions inFormat: matrix.matrixFormats forEquation:matrix.equation matrix:nil]];
			}
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
		if (variable.lastOperation && ![self.dependencies containsObject:variable.lastOperation] && variable.lastOperation !=self) {
			[self addDependency: variable.lastOperation];
		}
	}
	for (GLVariable *variable in self.result) {
		if (variable.lastOperation && ![self.dependencies containsObject:variable.lastOperation] && variable.lastOperation !=self) {
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

@implementation GLOptimizedVariableOperation

static NSMapTable *classVariableOperationTable = nil;
static NSMapTable *classOperandTable = nil;
static NSMapTable *classResultTable = nil;
static NSMapTable *classBufferTable = nil;
static NSMapTable *classPreOperationTable = nil;
static NSMapTable *classPostOperationTable = nil;

+ (NSUInteger) totalSubclasses {
    return  !classVariableOperationTable ? 0 : classVariableOperationTable.count;
}

+ (void) setVariableOperationPrototype: (variableOperation) opPrototype forSubclassWithName: (NSString *) subclassName {
    if (!classVariableOperationTable) {
        classVariableOperationTable = [NSMapTable strongToStrongObjectsMapTable];
    }
    [classVariableOperationTable setObject: opPrototype forKey: subclassName];
}

+ (variableOperation) variableOperationPrototypeForSubclassWithName: (NSString *) subclassName {
    return [classVariableOperationTable objectForKey: subclassName];
}

+ (void) setOperandPrototype: (NSArray *) operandPrototype forSubclassWithName: (NSString *) subclassName {
    if (!classOperandTable) {
        classOperandTable = [NSMapTable strongToStrongObjectsMapTable];
    }
    [classOperandTable setObject: operandPrototype forKey: subclassName];
}
+ (NSArray *) operandPrototypeForSubclassWithName: (NSString *) subclassName {
    return [classOperandTable objectForKey: subclassName];
}

+ (void) setResultPrototype: (NSArray *) resultPrototype forSubclassWithName: (NSString *) subclassName {
    if (!classResultTable) {
        classResultTable = [NSMapTable strongToStrongObjectsMapTable];
    }
    [classResultTable setObject: resultPrototype forKey: subclassName];
}
+ (NSArray *) resultPrototypeForSubclassWithName: (NSString *) subclassName {
    return [classResultTable objectForKey: subclassName];
}

+ (void) setBufferPrototype: (NSArray *) bufferPrototype forSubclassWithName: (NSString *) subclassName {
    if (!classBufferTable) {
        classBufferTable = [NSMapTable strongToStrongObjectsMapTable];
    }
    [classBufferTable setObject: bufferPrototype forKey: subclassName];
}
+ (NSArray *) bufferPrototypeForSubclassWithName: (NSString *) subclassName {
    return [classBufferTable objectForKey: subclassName];
}

+ (void) setPreOperationPrototype: (dispatch_block_t) preOpPrototype forSubclassWithName: (NSString *) subclassName {
    if (!classPreOperationTable) {
        classPreOperationTable = [NSMapTable strongToStrongObjectsMapTable];
    }
    [classPreOperationTable setObject: preOpPrototype forKey: subclassName];
}
+ (dispatch_block_t) preOperationPrototypeForSubclassWithName: (NSString *) subclassName {
    return [classPreOperationTable objectForKey: subclassName];
}

+ (void) setPostOperationPrototype: (dispatch_block_t) postOperationPrototype forSubclassWithName: (NSString *) subclassName {
    if (!classPostOperationTable) {
        classPostOperationTable = [NSMapTable strongToStrongObjectsMapTable];
    }
    [classPostOperationTable setObject: postOperationPrototype forKey: subclassName];
}
+ (dispatch_block_t) postOperationPrototypeForSubclassWithName: (NSString *) subclassName {
    return [classPostOperationTable objectForKey: subclassName];
}

- (id) initWithResult: (NSArray *) result operand: (NSArray *) operand buffers: (NSArray *) buffers operation: (variableOperation) op
{
    [NSException raise: @"InvalidInitializer" format: @"An optimized operation can only be initialized with -initWithOperand:"];
    return nil;
}

- (id) initWithResult: (NSArray *) result operand: (NSArray *) operand
{
    [NSException raise: @"InvalidInitializer" format: @"An optimized operation can only be initialized with -initWithOperand:"];
    return nil;
}

- (id) initWithOperand: (NSArray *) operand
{
    NSArray *operandPrototypes = [GLOptimizedVariableOperation operandPrototypeForSubclassWithName: NSStringFromClass(self.class)];
    if (!operandPrototypes) {
        [NSException raise: @"InvalidInitializer" format: @"Unable to find the implementation of this optimized subclass: %@", NSStringFromClass(self.class)];
    }
    if (operandPrototypes.count != operand.count) {
        [NSException raise: @"InvalidOperands" format: @"Invalid number of operands. Expected %lu, found %lu", operandPrototypes.count, operand.count];
    }
    for (NSUInteger i=0; i<operand.count; i++) {
        if ( ![[operandPrototypes[i] dimensions] isEqualToArray: [operand[i] dimensions]] ) {
            [NSException raise: @"DimensionsNotEqualException" format: @"The dimensions are not equal."];
        }
        if ([operandPrototypes[i] dataFormat] != [operand[i] dataFormat]) {
            [NSException raise: @"FormatsNotEqualException" format: @"The formats are not equal."];
        }
    }
    
	if (( self = [super init] ))
	{
        self.operand = operand;
        
        NSArray *resultPrototypes = [GLOptimizedVariableOperation resultPrototypeForSubclassWithName: NSStringFromClass(self.class)];
		NSMutableArray *array = [NSMutableArray array];
		for (GLVariable *variable in resultPrototypes) {
			[array addObject: [[variable class] variableOfType: variable.dataFormat withDimensions: variable.dimensions forEquation: variable.equation]];
		}
		self.result = array;
        
        NSArray *bufferPrototypes = [GLOptimizedVariableOperation bufferPrototypeForSubclassWithName: NSStringFromClass(self.class)];
		array = [NSMutableArray array];
		for (GLBuffer *buffer in bufferPrototypes) {
			[array addObject: [[GLBuffer alloc] initWithLength: buffer.numBytes]];
		}
		self.buffer = array;
		
        self.operation = [GLOptimizedVariableOperation variableOperationPrototypeForSubclassWithName: NSStringFromClass(self.class)];
		
		[self setupDependencies];
	}
	
    return self;
}

@end











