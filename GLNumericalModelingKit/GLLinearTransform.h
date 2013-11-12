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

@interface GLLinearTransform : GLVariable

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

/**  Returns the linear transformation of the given name, if it can be found or created.
 @discussion You may request derivatives in the form 'xxx' or known operators such as 'svv' or 'harmonicOperator'.
 @param name The name of the linear transformation.
 @param dimensions The dimensions that we are transforming from.
 @param equation The GLEquation object being used.
 @returns A GLLinearTransform with fromDimensions that match dimensions, and toDimensions that may be different.
 */
+ (GLLinearTransform *) linearTransformWithName: (NSString *) name forDimensions: (NSArray *) dimensions equation: (GLEquation *) equation;

/**  Add an operator to the pool so that it can be used later.
 @param transform The linear transform that you want to save.
 @param name The name of the linear transformation.
*/
+ (void) setLinearTransform: (GLLinearTransform *) transform withName: (NSString *) name;

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

/** Create a harmonic operator of arbitrary order that can act on multi-dimensional functions in the given dimensions.
 @param order Order n, is given by \nabla^{2n}
 @param dimensions Ordered dimensions (and therefore basis) which should be used to take the derivative.
 @param equation The GLEquation object being used.
 @returns A GLLinearTransform with fromDimensions that match dimensions, and toDimensions which may be different.
 */
+ (GLLinearTransform *) harmonicOperatorOfOrder: (NSUInteger) order fromDimensions: (NSArray *) dimensions forEquation: (GLEquation *) equation;

/** Create a harmonic operator of order 1 that can act on multi-dimensional functions in the given dimensions.
 @param dimensions Ordered dimensions (and therefore basis) which should be used to take the derivative.
 @param equation The GLEquation object being used.
 @returns A GLLinearTransform with fromDimensions that match dimensions, and toDimensions which may be different.
 */
+ (GLLinearTransform *) harmonicOperatorFromDimensions: (NSArray *) dimensions forEquation: (GLEquation *) equation;

/** Create the spectral vanishing viscosity operator for a set of (spectral) dimensions.
 @discussion This linear transformation is intended to act as a filter on some isotropic damping operator.
 @param dimensions Ordered dimensions that, at the moment, must all be spectral.
 @param isAntialiasing Whether or not we should assume the Nyquist is decreased by 2/3 for anti-aliasing.
 @param equation The GLEquation object being used.
 @returns A GLLinearTransform with fromDimensions and toDimensions that match dimensions.
 */
+ (GLLinearTransform *) spectralVanishingViscosityFilterWithDimensions: (NSArray *) dimensions scaledForAntialiasing: (BOOL) isAntialiasing forEquation: (GLEquation *) equation;

// Assuming the linear transform was initialized with fromDims that match the diagonalVariable and a format of diag in each dimension
- (void) setVariableAlongDiagonal: (GLFunction *) diagonalVariable;

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

/** Returns b in the equation A x = b, where A is this linear transformation
 @param x A function with dimensions matching the fromDimensions of this linear transformation
 @returns A transformed function with dimensions matching the toDimensions of this linear transformation
 */
- (GLFunction *) transform: (GLFunction *) x;

/** Returns x in the equation A x = b, where A is this linear transformation
 @param x A function with dimensions matching the toDimensions of this linear transformation
 @returns A solution function with dimensions matching the fromDimensions of this linear transformation
 */
- (GLFunction *) solve: (GLFunction *) b;

- (GLLinearTransform *) inverse;

- (GLLinearTransform *) matrixMultiply: (GLLinearTransform *) otherVariable;
@end
