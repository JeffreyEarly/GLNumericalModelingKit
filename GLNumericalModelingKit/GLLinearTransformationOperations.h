//
//  GLLinearTransformationOperations.h
//  GLNumericalModelingKitCopy
//
//  Created by Jeffrey J. Early on 1/18/13.
//
//

#import <GLNumericalModelingKit/GLVariableOperations.h>

#pragma mark -
#pragma mark Transforms
#pragma mark

/************************************************/
/*		GLSingleDiagonalTransformOperation		*/
/************************************************/

/** Find b in A x = b for matrices A that have exactly one diagonal in each dimension, be it a proper-, super-, or sub-diagonal.
 @param linearTransform The linear transform representing matrix A.
 @param function The function representing vector x.
 @returns The transformed function, b.
 */
@interface GLSingleDiagonalTransformOperation : GLVariableOperation
- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform function: (GLVariable *) function;
@end

/************************************************/
/*		GLTriadiagonalTransformOperation		*/
/************************************************/
/** Find b in A x = b for matrices A that are in tridiagonal format.
 @param linearTransform The linear transform representing matrix A.
 @param function The function representing vector x.
 @returns The solution function, b.
 */
@interface GLTriadiagonalTransformOperation : GLVariableOperation
- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform function: (GLVariable *) function;
@end

/************************************************/
/*		GLDenseMatrixTransformOperation         */
/************************************************/
/** Find b in A x = b for matrices A that are in dense format.
 @param linearTransform The linear transform representing matrix A.
 @param function The function representing vector x.
 @returns The solution function, b.
 */
@interface GLDenseMatrixTransformOperation : GLBinaryOperation
@end



#pragma mark -
#pragma mark Solvers
#pragma mark

/************************************************/
/*		GLTriadiagonalSolverOperation			*/
/************************************************/
/** Find x in A x = b for matrices A that are in tridiagonal format. In particular, this solve a_i x_{i-1} + b_i x_{i} + c_i x_{i+1} = d_i for x_i. Input a tridiagonal matrix (essentially a,b,c) and a vector (d), the function returns x.
 @param linearTransform The linear transform representing matrix A.
 @param function The function representing vector b.
 @returns The solution function, x.
 */
@interface GLTriadiagonalSolverOperation : GLVariableOperation
- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform function: (GLVariable *) function;
@end


/************************************************/
/*		GLDenseMatrixSolver						*/
/************************************************/
/** Find x in A x = b for matrices A that are in dense format.
 @param linearTransform The linear transform representing matrix A.
 @param function The function representing vector b.
 @returns The solution function, x.
 */
@interface GLDenseMatrixSolver : GLVariableOperation
- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform function: (GLVariable *) function;
@end


#pragma mark -
#pragma mark Manipulation
#pragma mark

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


