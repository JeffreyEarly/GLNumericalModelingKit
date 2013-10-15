//
//  GLLinearTransformationOperations.h
//  GLNumericalModelingKitCopy
//
//  Created by Jeffrey J. Early on 1/18/13.
//
//

#import <GLNumericalModelingKit/GLVariableOperations.h>

/************************************************/
/*		GLTriadiagonalOperation					*/
/************************************************/
// Solve a_i x_{i-1} + b_i x_{i} + c_i x_{i+1} = d_i for x_i.
// Input a tridiagonal matrix (essentially a,b,c) and a vector (d), the function returns x.
@interface GLTriadiagonalOperation : GLBinaryOperation
@end


/************************************************/
/*		GLTriadiagonalTransformOperation		*/
/************************************************/
// Solve a_i x_{i-1} + b_i x_{i} + c_i x_{i+1} = d_i for x_i.
// Input a tridiagonal matrix (essentially a,b,c) and a vector (d), the function returns x.
@interface GLTriadiagonalTransformOperation : GLBinaryOperation
@end

/************************************************/
/*		GLDenseMatrixTransformOperation         */
/************************************************/

@interface GLDenseMatrixTransformOperation : GLBinaryOperation
@end

/************************************************/
/*		GLMatrixMatrixMultiplicationOperation   */
/************************************************/

@interface GLMatrixMatrixMultiplicationOperation : GLBinaryOperation
@end

/************************************************/
/*		GLMatrixInversionOperation              */
/************************************************/

@interface GLMatrixInversionOperation : GLUnaryOperation
@end

/************************************************/
/*		GLLinearTransformAdditionOperation		*/
/************************************************/

@interface GLLinearTransformAdditionOperation : GLBinaryOperation
@end

/************************************************/
/*		GLDenseMatrixSolver						*/
/************************************************/
@interface GLDenseMatrixSolver : GLBinaryOperation
@end

/************************************************/
/*		GLSingleDiagonalTransformOperation		*/
/************************************************/
// Find b in A x = b for matrices A that have exactly one diagonal
// in each dimension, be a proper-, super-, or sub-diagonal.
@interface GLSingleDiagonalTransformOperation : GLVariableOperation
- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform function: (GLVariable *) function;
@end