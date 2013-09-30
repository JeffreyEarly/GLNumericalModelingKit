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

enum {
	kGLUndefinedOperation = 0,
    kGLUnaryOperation = 1,
	kGLBinaryOperation = 2,
	kGLUnaryVectorOperation = 3,
	kGLBinaryVectorOperation = 4,
	kGLNullaryOperation = 5
};
typedef NSUInteger GLOperationType;

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
typedef void (^nullaryOperation)(NSMutableData *);
typedef void (^unaryOperation)(NSMutableData *, NSData *);
typedef void (^binaryOperation)(NSMutableData *, NSData *, NSData *);
typedef void (^unaryVectorOperation)(NSArray *, NSArray *);
typedef void (^binaryVectorOperation)(NSArray *, NSArray *, NSArray *);

/************************************************/
/*		GLVariableOperation						*/
/************************************************/

@interface GLVariableOperation : NSOperation

@property GLOperationType operationType;

//@property(readwrite, strong, nonatomic) GLVariable *result;

@property BOOL isEnqueued;

// Returns YES if the result buffer can be the same as an operand buffer.
@property(readonly) BOOL canOperateInPlace;

- (BOOL) isEqualToOperation: (id) otherOperation;

@property(readwrite, copy, nonatomic) NSString *graphvisDescription;

@end

/************************************************/
/*		GLNullaryOperation						*/
/************************************************/

@interface GLNullaryOperation : GLVariableOperation

// Init with both the target variable and the operand. In-place operation is therefore possible.
- (id) initWithResult: (GLVariable *) result;

@property(readwrite, strong, nonatomic) GLVariable *result;

@property(copy) nullaryOperation blockOperation;

@end

/************************************************/
/*		GLUnaryOperation						*/
/************************************************/

@interface GLUnaryOperation : GLVariableOperation

// Init with the operand. The result variable is automatically created with the same properties as the operand.
- (id) initWithOperand: (GLVariable *) operand;

// Init with both the target variable and the operand. In-place operation is therefore possible.
- (id) initWithResult: (GLVariable *) result operand: (GLVariable *) operand;

@property(readwrite, strong, nonatomic) GLVariable *result;
@property(readwrite, strong, nonatomic) GLVariable *operand;

@property(copy) unaryOperation blockOperation;

@end


/************************************************/
/*		GLBinaryOperation						*/
/************************************************/

@interface GLBinaryOperation : GLVariableOperation
{
    GLVariable *_result;
	GLVariable *_firstOperand;
	GLVariable *_secondOperand;
}

// Init with only the operand to automatically create the target variable.
- (id) initWithFirstOperand: (GLVariable *) fOperand secondOperand: (GLVariable *) sOperand;

// Init with both the target variable and the operand. In-place negation is therefore possible.
- (id) initWithResult: (GLVariable *) result firstOperand: (GLVariable *) fOperand secondOperand: (GLVariable *) sOperand;

@property(readwrite, strong, nonatomic) GLVariable *result;
@property(readwrite, strong, nonatomic) GLVariable *firstOperand;
@property(readwrite, strong, nonatomic) GLVariable *secondOperand;

@property(copy) binaryOperation blockOperation;

@end


/************************************************/
/*		GLUnaryVectorOperation					*/
/************************************************/

@interface GLUnaryVectorOperation : GLVariableOperation

// Init with the operand. The result variable is automatically created with the same properties as the operand.
- (id) initWithOperand: (NSArray *) operand;

// Init with both the target variable and the operand. In-place operation is therefore possible.
- (id) initWithResult: (NSArray *) result operand: (NSArray *) operand;

@property(readwrite, strong, nonatomic) NSArray *result;
@property(readwrite, strong, nonatomic) NSArray *operand;

@property(copy) unaryVectorOperation blockOperation;

@end

/************************************************/
/*		GLBinaryVectorOperation					*/
/************************************************/

@interface GLBinaryVectorOperation : GLVariableOperation

// Init with the operand. The result variable is automatically created with the same properties as the operand.
- (id) initWithFirstOperand: (NSArray *) fOperand secondOperand: (NSArray *) sOperand;

// Init with both the target variable and the operand. In-place operation is therefore possible.
- (id) initWithResult: (NSArray *) result firstOperand: (NSArray *) fOperand secondOperand: (NSArray *) sOperand;

@property(readwrite, strong, nonatomic) NSArray *result;
@property(readwrite, strong, nonatomic) NSArray *firstOperand;
@property(readwrite, strong, nonatomic) NSArray *secondOperand;

@property(copy) binaryVectorOperation blockOperation;

@end





