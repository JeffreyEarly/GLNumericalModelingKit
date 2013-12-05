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

typedef struct {
    GLMatrixFormat format;
    
	NSUInteger nPoints;         // total number of points
	NSUInteger nRows;           // total number of rows
	NSUInteger nColumns;        // total number of columns
    NSUInteger nDiagonals;      // total number of diagonal (mutually exclusive with rows and columns).
    NSUInteger nDiagonalPoints; // total number of points along the diagonal (mutually exclusive with rows and columns).
    
	NSUInteger stride;          // distance to the next element
    
	NSUInteger rowStride;       // distance to the next row (may be the same as stride)
	NSUInteger columnStride;    // distance to the next column  (may be the same as stride)
    NSUInteger diagonalStride;  // distance to the next diagonal
} GLDataStride;

@class GLLinearTransform, GLFunction;
@interface GLMatrixDescription : NSObject

/** Creates a new matrix description based on the dimensions and formatting of the function.
 @discussion Functions are treated as "dense" column vectors.
 @param aFunction A function object.
 @returns A new GLMatrixDescription object.
 */
- (GLMatrixDescription *) initWithVariable: (GLFunction *) aFunction;


- (GLMatrixDescription *) initWithLinearTransform: (GLLinearTransform *) aLinearTransform;

@property NSUInteger nDimensions;

@property GLDataStride *strides;

// Returns yes if the matrices are in the same format.
- (BOOL) isEqualToMatrixDescription: (GLMatrixDescription *) otherMatrixDescription;

@end
