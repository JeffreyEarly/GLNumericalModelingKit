//
//  GLVariableOperations.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 4/21/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>

#import <GLNumericalModelingKit/Precision.h>
#import <GLNumericalModelingKit/GLVariable.h>


// GLVariableOperations operate on the data associated GLVariables.
//
// You *must* add as dependencies the operations that the other variables depend on.
// Note that if a mutable variable is being modified, the *result* variable will need
// its dependencies added as well.
//
// The instance variables that hold on to the dependency should NOT be weak.
// However, after the subclass is done computing, it SHOULD set the variable to nil
// and allow the dependency to be released.

// The first argument is always the result, subsequent arguments are the operands.
// The same is true of vector operations. Vector operations need not have the same
// sized arrays---that's entirely up to the implementation of each operation.

// The arguments are arrays of buffers (assumed to be NSMutableData objects.
// 1. output, 2. input, 3. internal usage
typedef void (^variableOperation)(NSArray *, NSArray *, NSArray *);

/************************************************/
/*		GLVariableOperation						*/
/************************************************/

@interface GLVariableOperation : NSOperation

@property BOOL isEnqueued;

// Returns YES if the result buffer can be the same as an operand buffer.
@property(readonly) BOOL canOperateInPlace;

- (BOOL) isEqualToOperation: (id) otherOperation;

@property(readwrite, copy, nonatomic) NSString *graphvisDescription;


- (id) initWithResult: (NSArray *) result operand: (NSArray *) operand buffers: (NSArray *) buffers operation: (variableOperation) op;

// This assumes that there are no internal buffers necessary, and that the operation will be set later.
- (id) initWithResult: (NSArray *) result operand: (NSArray *) operand;

// This assumes that the result variables are the same type as the operand variables.
- (id) initWithOperand: (NSArray *) operand;

@property(readwrite, strong, nonatomic) NSArray *result;
@property(readwrite, strong, nonatomic) NSArray *operand;
@property(readwrite, strong, nonatomic) NSArray *buffers;

@property(copy) variableOperation operation;

@end




