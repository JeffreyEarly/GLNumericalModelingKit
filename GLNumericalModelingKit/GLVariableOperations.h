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
#import <GLNumericalModelingKit/GLBuffer.h>

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

/** Creates a new variable operation class optimized for the operation graph defined by the given operand and result variables.
 @discussion Instances of this class should be created with -initWithOperand.
 @param operandVariables An array of GLVariable objects that define the 'top' of the operation graph.
 @param resultVariables An array of GLVariable objects that define the 'bottom' of the operation graph.
 @returns A new GLVariableOperation subclass..
 */
+ (Class) variableOperationSubclassWithOperand: (NSArray *) operandVariables result: (NSArray *) resultVariables;

/** Creates a new variable operation objects.
 @discussion Do not call this method on a fully implemented subclass or you will overwrite the internal behavior.
 @param result An array of GLVariable objects that store the results from the operation.
 @param operand An array of GLVariable objects that store the operands of the operation.
 @param buffers An array of NSMutableData objects used a buffers during the operations.
 @param op The actual operation that performs the calculation on the raw buffers.
 @returns A new variable operation object.
 */
- (instancetype) initWithResult: (NSArray *) result operand: (NSArray *) operand buffers: (NSArray *) buffers operation: (variableOperation) op;

/** Creates a new variable operation objects.
 @discussion Do not call this method on a fully implemented subclass or you will overwrite the internal behavior.
 @discussion This assumes that there are no internal buffers necessary, and that the operation will be set later.
 @param result An array of GLVariable objects that store the results from the operation.
 @param operand An array of GLVariable objects that store the operands of the operation.
 @returns A new variable operation object.
 */
- (instancetype) initWithResult: (NSArray *) result operand: (NSArray *) operand;

/** Creates a new variable operation objects.
 @discussion This assumes that there are no internal buffers necessary, and that the operation will be set later.
 @discussion This also assumes that the result variables are the same type as the operand variables.
 @param operand An array of GLVariable objects that store the operands of the operation.
 @returns A new variable operation object.
 */
- (instancetype) initWithOperand: (NSArray *) operand;

/// An array of GLVariable objects that store the results from the operation.
@property(readwrite, strong, nonatomic) NSArray *result;

/// An array of GLVariable objects that store the operands of the operation.
@property(readwrite, strong, nonatomic) NSArray *operand;

/// An array of GLBuffer objects indicating the length, in bytes, of the buffers required during the operations.
@property(readwrite, strong, nonatomic) NSArray *buffer;

/// The operation that performs the calculation on the raw buffers
@property(copy) variableOperation operation;

/// If set, this block will be run sometime before the variableOperation is executed. This would be a logical place create a specific buffer to be used a runtime.
@property(copy) dispatch_block_t preOperation;

/// If set, this block will be run sometime after the variableOperation is executed. Useful for deallocating/cleaning up.
@property(copy) dispatch_block_t postOperation;

/// Whether or not this operation is already on the queue, slated to be run.
@property BOOL isEnqueued;

/// Returns YES if the result buffer can be the same as an operand buffer.
@property(readonly) BOOL canOperateInPlace;

/// Returns yes if another operation is of the same class, and has the same result and operand variables.
- (BOOL) isEqualToOperation: (id) otherOperation;

/// A description of this operation appropriate for output to graphvis.
@property(readwrite, copy, nonatomic) NSString *graphvisDescription;

@end




