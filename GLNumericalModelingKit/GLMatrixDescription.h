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
enum {
	kGLIdentityMatrixFormat = 0,
	kGLDenseMatrixFormat = 1,
    kGLDiagonalMatrixFormat = 2,
	kGLSubdiagonalMatrixFormat = 3,
	kGLSuperdiagonalMatrixFormat = 4,
	kGLTridiagonalMatrixFormat = 5
};
typedef NSUInteger GLMatrixFormat;

enum {
	kGLRowMatrixOrder = 0,
	kGLColumnMatrixOrder = 1
};
typedef NSUInteger GLMatrixOrder;

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

@class GLLinearTransform, GLVariable;
@interface GLMatrixDescription : NSObject

- (GLMatrixDescription *) initWithVariable: (GLVariable *) variable;
- (GLMatrixDescription *) initWithLinearTransform: (GLLinearTransform *) aLinearTransform;

@property NSUInteger nDimensions;

@property GLDataStride *strides;

@end
