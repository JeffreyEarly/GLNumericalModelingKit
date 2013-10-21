//
//  GLVectorVectorOperations.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 1/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import <GLNumericalModelingKit/GLVariableOperations.h>

/************************************************/
/*		GLAdditionOperation						*/
/************************************************/
// variable = leftVariable + rightVariable
/** C = A + B. This operation takes scalar-scalar, scalar-vector, vector-vector, scalar-matrix and matrix-matrix arguments. At the moment, the matrices must be in the same format.
 @param A An input scalar, vector or matrix.
 @param B An input scalar, vector or matrix.
 @returns The result C.
 */
@interface GLAdditionOperation : GLVariableOperation
- (id) initWithFirstOperand: (GLTensor *) A secondOperand: (GLTensor *) B;
@end


/************************************************/
/*		GLSubtractionOperation					*/
/************************************************/
// variable = leftVariable - rightVariable
@interface GLSubtractionOperation : GLBinaryOperation
@end


/************************************************/
/*		GLMultiplicationOperation				*/
/************************************************/
// variable = leftVariable * rightVariable
@interface GLMultiplicationOperation : GLBinaryOperation
@end

/************************************************/
/*		GLAbsoluteLargestOperation				*/
/************************************************/
// variable = max( abs(leftVariable), abs(rightVariable ) element-wise
@interface GLAbsoluteLargestOperation : GLBinaryOperation
@end

/************************************************/
/*		GLDivisionOperation						*/
/************************************************/
// variable = leftVariable / rightVariable
@interface GLDivisionOperation : GLBinaryOperation

// You can optionally set useComplexDivision to NO, and it will simply do an element-wise divide, ignoring complex math.
- (id) initWithFirstOperand: (GLVariable *) fOperand secondOperand: (GLVariable *) sOperand shouldUseComplexArithmetic: (BOOL) useComplexArithmetic;
@property BOOL useComplexArithmetic;

@end

/************************************************/
/*		GLDotProductOperation					*/
/************************************************/

// variable = leftVariable dot rightVariable
@interface GLDotProductOperation : GLBinaryOperation
@end

/************************************************/
/*		GLSetVariableValueOperation				*/
/************************************************/

@interface GLSetVariableValueOperation : GLBinaryOperation
- (id) initWithVectorOperand: (GLVariable *) variable scalarVariableOperand: (GLVariable *) aScalarVariable indexString: (NSString *) indexString;
@property(copy) NSString * indexString;
@end