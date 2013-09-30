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

// At the moment the from dimensions need to be the same rank as the to dimensions.
@property(readwrite, strong) NSArray *fromDimensions;
@property(readwrite, strong) NSArray *toDimensions;

// row-major or column-major
@property(readwrite) GLMatrixOrder matrixOrder;

// The end-all-be-all description of how this matrix is stored in memory.
@property(readwrite, strong) GLMatrixDescription *matrixDescription;

// A subset of the above, this simply indicates the type of matrix in each dimension.
@property(readwrite, strong) NSArray *matrixFormats;

// A block that can be used to create the matrix.
@property(copy) transformMatrix matrixBlock;

+ (id) transformOfType: (GLDataFormat) dataFormat withFromDimensions: (NSArray *) fromDims toDimensions: (NSArray *) toDims inFormat: (NSArray *) matrixFormats forEquation: (GLEquation *) equation;

- (id) initTransformOfType: (GLDataFormat) dataFormat withFromDimensions: (NSArray *) fromDims toDimensions: (NSArray *) toDims inFormat: (NSArray *) matrixFormats forEquation: (GLEquation *) theEquation;

+ (id) dftMatrixFromDimension: (GLDimension *) aDimension forEquation: (GLEquation *) equation;
+ (id) idftMatrixFromDimension: (GLDimension *) aDimension forEquation: (GLEquation *) equation;

+ (id) cosineTransformMatrixFromDimension: (GLDimension *) aDimension forEquation: (GLEquation *) equation;
+ (id) inverseCosineTransformMatrixFromDimension: (GLDimension *) aDimension forEquation: (GLEquation *) equation;

+ (id) differentiationMatrixFromDimension: (GLDimension *) aDimension forEquation: (GLEquation *) equation;

// Assuming the linear transform was initialized with fromDims that match the diagonalVariable and a format of diag in each dimension
- (void) setVariableAlongDiagonal: (GLVariable *) diagonalVariable;

// Starting with the subdiagonal, diagonal, superdiagonal.
- (void) setVariablesAlongTridiagonal: (NSArray *) tridiagonalVariables;

// Returns x in the equation A x = b, where A is this linear transformation
- (GLVariable *) solve: (GLVariable *) b;

// Returns x in the equation A x = b, where A is this linear transformation
- (GLVariable *) tridiagonalSolveWithVector: (GLVariable *)	b;

// Returns b in the equation A x = b, where A is this linear transformation
- (GLVariable *) transform: (GLVariable *) x;

- (GLLinearTransform *) inverse;

- (GLLinearTransform *) plus: (GLLinearTransform *) otherVariable;

@end
