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
- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform function: (GLFunction *) function;
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
- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform function: (GLFunction *) function;
@end

/************************************************/
/*		GLDenseMatrixTransformOperation         */
/************************************************/
/** Find b in A x = b for matrices A that are in dense format.
 @param linearTransform The linear transform representing matrix A.
 @param function The function representing vector x.
 @returns The solution function, b.
 */
@interface GLDenseMatrixTransformOperation : GLVariableOperation
- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform function: (GLFunction *) function;
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
- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform function: (GLFunction *) function;
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
- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform function: (GLFunction *) function;
@end


#pragma mark -
#pragma mark Manipulation
#pragma mark

/************************************************/
/*		GLMatrixMatrixMultiplicationOperation   */
/************************************************/

@interface GLMatrixMatrixMultiplicationOperation : GLVariableOperation
- (id) initWithFirstOperand: (GLLinearTransform *) A secondOperand: (GLLinearTransform *) B;
@end

/************************************************/
/*		GLMatrixInversionOperation              */
/************************************************/

@interface GLMatrixInversionOperation : GLVariableOperation
- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform;
@end

/************************************************/
/*		GLMatrixEigensystemOperation             */
/************************************************/
/** Finds the eigenvectors and eigenvalues of the matrix A.
 @discussion The eigenvalues are returned as the linear transformation S which takes vectors from the eigenbasis to the original basis.
 @discussion The matrix A must be an endomorphism, meaning that its fromBasis must be the same as its toBasis.
 @param linearTransform The linear transform representing matrix A.
 @returns An NSArray containing the function v of the eigenvalues and the linear transformation S.
 */
@interface GLMatrixEigensystemOperation : GLVariableOperation
- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform;
@end


