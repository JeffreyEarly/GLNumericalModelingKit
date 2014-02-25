//
//  GLLinearTransformationOperations.h
//  GLNumericalModelingKitCopy
//
//  Created by Jeffrey J. Early on 1/18/13.
//
//

#import <GLNumericalModelingKit/GLVariableOperations.h>

/*!
 * @function apply_matrix_loop
 *
 * @abstract
 * Loops over other matrix dimensions, allowing you to treat one dimension in isolation from the others. Given an equation of the form A x = b.
 */
void apply_matrix_vector_loop( GLMatrixDescription *matrixDescription, GLMatrixDescription *vectorDescription, NSUInteger loopIndex, dispatch_queue_t queue, void (^block)(NSUInteger, NSUInteger, NSUInteger ));

NSUInteger compute_total_matrix_vector_loops( GLMatrixDescription *matrixDescription, GLMatrixDescription *vectorDescription, NSUInteger loopIndex );

void apply_matrix_matrix_loop( GLMatrixDescription *matrixA, GLMatrixDescription *matrixB, GLMatrixDescription *matrixC, NSUInteger loopIndex, dispatch_queue_t queue, void (^block)(NSUInteger, NSUInteger, NSUInteger ));


/************************************************/
/*                                              */
/*		Creation                                */
/*                                              */
/************************************************/

#pragma mark -
#pragma mark Creation
#pragma mark

/************************************************/
/*		GLDiagonalTransformCreationOperation	*/
/************************************************/

/** Places the values of the function along the diagonal of matrix
 @param function The function to be placed along the diagonal.
 @returns The transformed function with fromDimensions and toDimensions matching the dimension of the function.
 */
@interface GLDiagonalTransformCreationOperation : GLVariableOperation
- (id) initWithFunction: (GLFunction *) function;
@end

/************************************************/
/*		GLExpandMatrixDimensionsOperation		*/
/************************************************/

/** Add new dimension to a transformation, by assuming the identity transform for the new dimensions.
 @param linearTransform Existing linear transformation.
 @param fromDims New fromDimensions
 @param toDims New toDimensions
 @returns The more-or-less the same linear transformation, but with new fromDimensions and toDimensions.
 */
@interface GLExpandMatrixDimensionsOperation : GLVariableOperation
- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform fromDimensions: (NSArray *) fromDims toDimensions: (NSArray *) toDims;
@end

/************************************************/
/*		GLReduceMatrixDimensionsOperation		*/
/************************************************/

/** Reduce the dimensions of the linear transformation.
 @param linearTransform Existing linear transformation.
 @param fromString A string, in matlab format, describing the subset of indices requested.
 @param toString A string, in matlab format, describing the subset of indices requested.
 @returns The more-or-less the same linear transformation, but with new fromDimensions and toDimensions.
 */

@interface GLReduceMatrixDimensionsOperation : GLVariableOperation
- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform fromDimensionsIndexString: (NSString *) fromString toDimensionsIndexString: (NSString *) toString;
@end

/************************************************/
/*                                              */
/*		Vector Multiplication (Transforms)      */
/*                                              */
/************************************************/

#pragma mark -
#pragma mark Vector Multiplication (Transforms)
#pragma mark

// Operations that find b, in A x = b

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
 @discussion The matrices must have exactly one tridiagonal dimension. The remaining dimensions can be diagonal or identity.
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
  @discussion The matrices must have exactly one dense dimension. The remaining dimensions can be diagonal or identity.
 @param linearTransform The linear transform representing matrix A.
 @param function The function representing vector x.
 @returns The solution function, b.
 */
@interface GLDenseMatrixTransformOperation : GLVariableOperation
- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform function: (GLFunction *) function;
@end



/************************************************/
/*                                              */
/*		Matrix Multiplication                   */
/*                                              */
/************************************************/

#pragma mark -
#pragma mark Matrix Multiplication
#pragma mark

// Operations that find C, in A.B = C

/********************************************************/
/*		GLDiagonalMatrixMatrixMultiplicationOperation   */
/********************************************************/

@interface GLDiagonalMatrixMatrixMultiplicationOperation : GLVariableOperation
- (id) initWithFirstOperand: (GLLinearTransform *) A secondOperand: (GLLinearTransform *) B;
@end

/************************************************/
/*		GLMatrixMatrixMultiplicationOperation   */
/************************************************/

@interface GLMatrixMatrixMultiplicationOperation : GLVariableOperation
- (id) initWithFirstOperand: (GLLinearTransform *) A secondOperand: (GLLinearTransform *) B;
@end

/************************************************/
/*		GLMatrixMatrixMultiplicationOperation   */
/************************************************/

@interface GLMatrixMatrixDiagonalDenseMultiplicationOperation : GLVariableOperation
- (id) initWithFirstOperand: (GLLinearTransform *) A secondOperand: (GLLinearTransform *) B;
@end

/************************************************/
/*                                              */
/*		Solvers                              */
/*                                              */
/************************************************/

#pragma mark -
#pragma mark Solvers
#pragma mark

// Operations that find x in A x = b

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



/************************************************/
/*                                              */
/*		Inversion                               */
/*                                              */
/************************************************/

#pragma mark -
#pragma mark Inversion
#pragma mark

// Operations that find B in B.A=I

/************************************************/
/*		GLMatrixInversionOperation              */
/************************************************/

@interface GLMatrixInversionOperation : GLVariableOperation
- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform;
@end



/************************************************/
/*                                              */
/*		Other                                   */
/*                                              */
/************************************************/

#pragma mark -
#pragma mark Other
#pragma mark

/************************************************/
/*		GLMatrixNormalizationOperation          */
/************************************************/

#pragma mark -
#pragma mark GLMatrixNormalizationOperation
#pragma mark

@interface GLMatrixNormalizationOperation : GLVariableOperation
- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform normalizationConstant: (GLFloat) aConst dimensionIndex: (NSUInteger) index;
- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform normalizationFunction: (GLFunction *) aFunction;
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
- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform sort: (NSComparisonResult) sortOrder;
@end

/************************************************/
/*		GLGeneralizedMatrixEigensystemOperation */
/************************************************/
/** Finds the eigenvectors and eigenvalues of the matrices A and B.
 @discussion The eigenvalues are returned as the linear transformation S which takes vectors from the eigenbasis to the original basis.
 @discussion The matrices A and B must be an endomorphism, meaning that its fromBasis must be the same as its toBasis.
 @param A The linear transform representing matrix A.
 @param B The linear transform representing matrix B.
 @returns An NSArray containing the function v of the eigenvalues and the linear transformation S.
 */
@interface GLGeneralizedMatrixEigensystemOperation : GLVariableOperation
- (id) initWithFirstOperand: (GLLinearTransform *) A secondOperand: (GLLinearTransform *) B;
- (id) initWithFirstOperand: (GLLinearTransform *) A secondOperand: (GLLinearTransform *) B sort: (NSComparisonResult) sortOrder;
@end

/************************************************/
/*		GLTensorProductOperation				*/
/************************************************/
/** Takes the tensor product (outer product) of an array linear transformations.
 @param linearTransformations An array of linear transformations.
 @returns A linear transformation with fromDimensions and toDimensions of all the linear transformations in the array.
 */
@interface GLTensorProductOperation : GLVariableOperation
- (id) initWithLinearTransformations: (NSArray *) linearTransformations;
@end

/************************************************/
/*		GLFormatShiftOperation				*/
/************************************************/
/** Copy the linear transformation with the data formatted according to new parameters.
 @discussion This is not implemented for speed (and therefore shouldn't be used in a loop where speed is needed), but it should be able to convert to any format.
 @param linearTransformation A linear transformation
 @param dataFormat Specify the desired data format.
 @param matrixFormats An array of GLMatrixFormat types specifying the necessary storage requires that should be allocated for a particular dimension pair.
 @param ordering Whether the dense matrix indices should be column and row-major ordered.
 @returns A new created GLLinearTransform instance in the requested format..
 */
@interface GLFormatShiftOperation : GLVariableOperation
- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform dataType: (GLDataFormat) dataFormat matrixFormat: (NSArray *) matrixFormats ordering: (GLMatrixOrder) ordering;
@end


