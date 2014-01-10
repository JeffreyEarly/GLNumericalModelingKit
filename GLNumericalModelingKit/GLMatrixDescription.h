//
//  GLMatrixDescription.h
//  GLNumericalModelingKitCopy
//
//  Created by Jeffrey Early on 1/24/13.
//
//

#import <Foundation/Foundation.h>

// The format must be specified for each set of dimensions.
// The identity matrix format assumes an identity transformation, and therefore has no storage requirements.
// The dense matrix format assume all rows and columns must be specified for that dimension.
// The diagonal matrix format assume only the diagonal must be specified.
typedef NS_ENUM(NSUInteger, GLMatrixFormat) {
	kGLIdentityMatrixFormat = 0,
	kGLDenseMatrixFormat = 1,
    kGLDiagonalMatrixFormat = 2,
	kGLSubdiagonalMatrixFormat = 3,
	kGLSuperdiagonalMatrixFormat = 4,
	kGLTridiagonalMatrixFormat = 5
};

typedef NS_ENUM(NSUInteger, GLMatrixOrder) {
	kGLRowMatrixOrder = 0,
	kGLColumnMatrixOrder = 1
};

// This specifies how the data is organized in the memory buffer.
// kGLRealDataFormat means that there is no memory allocated for the imaginary part.
// kGLSplitComplexDataFormat means that the imaginary nPoints follow the real nPoints in the buffer.
// kGLInterleavedComplexFormat means that the imaginary part of a point immediately follows the real part. Strides are set to 2.
typedef NS_ENUM(NSUInteger, GLDataFormat) {
	kGLRealDataFormat = 0,
    kGLSplitComplexDataFormat = 1,
    kGLInterleavedComplexDataFormat = 2
};

// Note all strides (distances) are measure in sizeof(GLFloat).
// The complex stride will be 0 for real numbers, 1 for complex interleaved numbers and nPoints (total point in all dims) for split formats.
typedef struct {
    GLMatrixFormat matrixFormat;
	GLDataFormat dataFormat;
    
	NSUInteger nPoints;         // total number of points
	NSUInteger nRows;           // total number of rows
	NSUInteger nColumns;        // total number of columns
    NSUInteger nDiagonals;      // total number of diagonal (mutually exclusive with rows and columns).
    NSUInteger nDiagonalPoints; // total number of points along the diagonal (mutually exclusive with rows and columns).
    
	NSUInteger stride;          // distance to the next real element
	NSUInteger complexStride;	// distance to the imaginary part of the element
    
	NSUInteger rowStride;       // distance to the next row (may be the same as stride)
	NSUInteger columnStride;    // distance to the next column  (may be the same as stride)
    NSUInteger diagonalStride;  // distance to the next diagonal
} GLDataStride;

@class GLLinearTransform, GLFunction;
@interface GLMatrixDescription : NSObject

/** Creates a new matrix description based on the dimensions and formatting of the function.
 @discussion Functions are treated as "dense" column vectors.
 @discussion Row-major and column-major ordering are meaningless for a function.
 @param aFunction A function object.
 @returns A new GLMatrixDescription object.
 */
- (GLMatrixDescription *) initWithFunction: (GLFunction *) aFunction;

/** Creates a new matrix description based on the dimensions and formatting of the linear transformation.
 @param aLinearTransform A linear transformation object.
 @returns A new GLMatrixDescription object.
 */
- (GLMatrixDescription *) initWithLinearTransform: (GLLinearTransform *) aLinearTransform;

/// Total number of dimensions that are represented. In the case of linear transformation, this is the total number of toDimensions (or fromDimensions).
@property NSUInteger nDimensions;

/// Total number of points *stored*.
@property NSUInteger nPoints;

/// Total number of elements (double the number of points for complex formats)
@property NSUInteger nElements;

/// Total number of bytes required to store the matrix.
@property NSUInteger nBytes;

/// The stride/distance to the imaginary part of the point.
@property NSUInteger complexStride;

/// The stride/distance to the next point.
@property NSUInteger elementStride;

/// The format being used to store the data
@property GLDataFormat dataFormat;

/// An array of size nDimensions
@property GLDataStride *strides;

/// Returns yes if the matrices are in the same format.
- (BOOL) isEqualToMatrixDescription: (GLMatrixDescription *) otherMatrixDescription;

/// Returns the smallest format required to store both the left and right format.
+ (NSArray *) commonFormatsFromLeft: (NSArray *) A right: (NSArray *) B;

@end
