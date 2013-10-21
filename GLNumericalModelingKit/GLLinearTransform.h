//
//  GLLinearTransform.h
//  GLNumericalModelingKitCopy
//
//  Created by Jeffrey J. Early on 1/22/13.
//
//

#import <GLNumericalModelingKit/GLNumericalModelingKit.h>
#import <GLNumericalModelingKit/GLMatrixDescription.h>
#import <complex.h>

// To convert from an element to an index: vec[(i*ny+j)*nz+k] = matrix[i][j][k]

// Simple example:
// If we have a vector along the diagonal in one dimension, and the identity in the others, the diagonal vector is simply repeated in the other dimensions along the diagonal.
//

enum {
	kGLFourierTransform,
    kGLFourierInverseTransform,
	kGLCosineTransform,
    kGLCosineInverseTransform,
    kGLSineTransform,
    kGLSineInverseTransform,
    kGLChebyshevTransform,
    kGLChebyshevInverseTransform,
    kGLDifferentiationTransform
};
typedef NSUInteger GLLinearTransformType;

// This block type is how you define a matrix, independent of how it is stored.
// The first argument indicates the row/destination indices, the second argument indicates the column/starting indices.
typedef GLFloatComplex (^transformMatrix)(NSUInteger *, NSUInteger *);

@interface GLLinearTransform : GLTensor

/************************************************/
/*		Initialization							*/
/************************************************/

#pragma mark -
#pragma mark Initialization
#pragma mark

/** Create a new, full specified, GLLinearTransform.
@param dataFormat Specify whether or not the transformation is real or complex.
@param fromDims An array of dimensions specifying the basis to transform from.
@param toDims An array of dimensions specifying the basis to transform to.
@param matrixFormats An array of GLMatrixFormat types specifying the necessary storage requires that should be allocated for a particular dimension pair.
@param equation The GLEquation object being used.
@param matrix A transformMatrix block capable of populating the complete matrix.
@returns A new created GLLinearTransform instance.
*/
+ (GLLinearTransform *) transformOfType: (GLDataFormat) dataFormat withFromDimensions: (NSArray *) fromDims toDimensions: (NSArray *) toDims inFormat: (NSArray *) matrixFormats forEquation: (GLEquation *) equation matrix:(GLFloatComplex (^)(NSUInteger *, NSUInteger *)) matrix;

/** Create a new GLLinearTransform. This initializer will automatically determine the type (real or complex) and the necessary storage format (dense, diagonal, etc.), although at the cost of greater initialization time. If initialization time is an issue, use the fully specified format above.
 @param fromDims An array of dimensions specifying the basis to transform from.
 @param toDims An array of dimensions specifying the basis to transform to.
 @param equation The GLEquation object being used.
 @param matrix A transformMatrix block capable of populating the complete matrix.
 @returns A new created GLLinearTransform instance.
 */
+ (GLLinearTransform *) transformWithFromDimensions: (NSArray *) fromDims toDimensions: (NSArray *) toDims forEquation: (GLEquation *) equation matrix:(GLFloatComplex (^)(NSUInteger *, NSUInteger *)) matrix;

/** Create a new, full specified, GLLinearTransform.
 @param dataFormat Specify whether or not the transformation is real or complex.
 @param fromDims An array of dimensions specifying the basis to transform from.
 @param toDims An array of dimensions specifying the basis to transform to.
 @param matrixFormats An array of GLMatrixFormat types specifying the necessary storage requires that should be allocated for a particular dimension pair.
 @param equation The GLEquation object being used.
 @param matrix A transformMatrix block capable of populating the complete matrix.
 @returns A new created GLLinearTransform instance.
*/
- (GLLinearTransform *) initTransformOfType: (GLDataFormat) dataFormat withFromDimensions: (NSArray *) fromDims toDimensions: (NSArray *) toDims inFormat: (NSArray *) matrixFormats forEquation: (GLEquation *) theEquation matrix:(GLFloatComplex (^)(NSUInteger *, NSUInteger *)) matrix;

/************************************************/
/*		Pre-defined transformations             */
/************************************************/

#pragma mark -
#pragma mark Pre-defined transformations
#pragma mark

/** Create a discrete transform that can act on 1-dimensional functions in the given dimension and transform them to the new basis.
 @param aDimension The dimension (and therefore basis) that we are transforming from.
 @param aBasis The basis that we are transforming to.
 @param equation The GLEquation object being used.
 @returns A GLLinearTransform with fromDimensions that match aDimension, and toDimensions in the new basis.
 */
+ (GLLinearTransform *) discreteTransformFromDimension: (GLDimension *) aDimension toBasis: (GLBasisFunction) aBasis forEquation: (GLEquation *) equation;

/** Create a differential operator of arbitrary order that can act on 1-dimensional functions in the given dimension.
 @param numDerivs Order of differentiation.
 @param aDimension The dimension (and therefore basis) which should be used to take the derivative.
 @param equation The GLEquation object being used.
 @returns A GLLinearTransform with fromDimensions that match aDimension, and toDimensions which may be different.
 */
+ (GLLinearTransform *) differentialOperatorOfOrder: (NSUInteger) numDerivs fromDimension: (GLDimension *) aDimension forEquation: (GLEquation *) equation;

/** Create a differential operator of arbitrary order that can act on multi-dimensional functions in the given dimensions.
 @param numDerivs An array of NSNumbers indicating the order of differentiation for each dimension.
 @param dimensions Ordered dimensions (and therefore basis) which should be used to take the derivative.
 @param equation The GLEquation object being used.
 @returns A GLLinearTransform with fromDimensions that match dimensions, and toDimensions which may be different.
 */
+ (GLLinearTransform *) differentialOperatorWithDerivatives: (NSArray *) numDerivs fromDimensions: (NSArray *) dimensions forEquation: (GLEquation *) equation;

// Assuming the linear transform was initialized with fromDims that match the diagonalVariable and a format of diag in each dimension
- (void) setVariableAlongDiagonal: (GLVariable *) diagonalVariable;

// Starting with the subdiagonal, diagonal, superdiagonal.
- (void) setVariablesAlongTridiagonal: (NSArray *) tridiagonalVariables;

/************************************************/
/*		Dimensionality							*/
/************************************************/

#pragma mark -
#pragma mark Dimensionality
#pragma mark

/// At the moment the from dimensions need to be the same rank as the to dimensions.
@property(readwrite, strong) NSArray *fromDimensions;
@property(readwrite, strong) NSArray *toDimensions;

/// Does this even make sense? I don't think any of it does here.

/// The symmetry (none, even, or odd) of each dimension.
@property(readwrite, assign, nonatomic) NSMutableArray *realSymmetry;
@property(readwrite, assign, nonatomic) NSMutableArray *imaginarySymmetry;

/// Returns YES if the variable is Hermitian,  H(-f)=H†(f).
@property(readonly, assign, nonatomic) BOOL isHermitian;

// Return the dimension restricted to positive values that can be used to recover
// the negative values. Returns nil if the variable is not hermitian.
@property(readonly, strong, nonatomic) GLDimension *hermitianDimension;


/************************************************/
/*		Data									*/
/************************************************/

#pragma mark -
#pragma mark Data
#pragma mark


// row-major or column-major
@property(readwrite) GLMatrixOrder matrixOrder;

// The end-all-be-all description of how this matrix is stored in memory.
@property(readwrite, strong) GLMatrixDescription *matrixDescription;

// A subset of the above, this simply indicates the type of matrix in each dimension.
@property(readwrite, strong) NSArray *matrixFormats;

// A block that can be used to create the matrix.
@property(copy) transformMatrix matrixBlock;


/************************************************/
/*		Operations								*/
/************************************************/

#pragma mark -
#pragma mark Operations
#pragma mark

// Returns b in the equation A x = b, where A is this linear transformation
- (GLVariable *) transform: (GLVariable *) x;

// Returns x in the equation A x = b, where A is this linear transformation
- (GLVariable *) solve: (GLVariable *) b;

// Returns x in the equation A x = b, where A is this linear transformation
- (GLVariable *) tridiagonalSolveWithVector: (GLVariable *)	b;

- (GLLinearTransform *) inverse;

- (GLLinearTransform *) plus: (GLLinearTransform *) otherVariable;

@end
