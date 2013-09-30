//
//  GLMatrixDescription.m
//  GLNumericalModelingKitCopy
//
//  Created by Jeffrey Early on 1/24/13.
//
//

#import "GLMatrixDescription.h"
#import "GLLinearTransform.h"

@interface GLMatrixDescription ()
@property(readwrite, strong) NSMutableData *stridesData;
@end

@implementation GLMatrixDescription

- (GLMatrixDescription *) initWithVariable: (GLVariable *) variable
{
    if ((self=[super init])) {
        self.nDimensions = variable.dimensions.count;
        self.stridesData = [NSMutableData dataWithLength: self.nDimensions*sizeof(GLDataStride)];
        self.strides = self.stridesData.mutableBytes;
        
        for (NSInteger iDim=variable.dimensions.count-1; iDim >= 0; iDim--)
		{
			GLDimension *dim = variable.dimensions[iDim];
            
            self.strides[iDim].format = kGLDenseMatrixFormat;
            
            self.strides[iDim].nPoints = dim.nPoints;
            self.strides[iDim].nRows = dim.nPoints;
            self.strides[iDim].nColumns = 1;
			
            if (iDim == variable.dimensions.count-1) {
                self.strides[iDim].stride = 1;
            } else {
                self.strides[iDim].stride = self.strides[iDim+1].nPoints * self.strides[iDim+1].stride;
            }
  			
			// Finally, we need to determine how to get to the next row or column within this matrix.
            // Thsi shouldn't matter for column vectors. 
			if ( 1 ) {
                self.strides[iDim].rowStride = self.strides[iDim].stride * self.strides[iDim].nColumns;
                self.strides[iDim].columnStride = self.strides[iDim].stride;
			} else if (kGLRowMatrixOrder == kGLColumnMatrixOrder) {
                self.strides[iDim].rowStride = self.strides[iDim].stride;
                self.strides[iDim].columnStride = self.strides[iDim].stride * self.strides[iDim].nRows;
			}
		}
    }
    return self;
}

- (GLMatrixDescription *) initWithLinearTransform: (GLLinearTransform *) linearTransform
{
    if ((self=[super init])) {
        self.nDimensions = linearTransform.fromDimensions.count;
        self.stridesData = [NSMutableData dataWithLength: self.nDimensions*sizeof(GLDataStride)];
        self.strides = self.stridesData.mutableBytes;
        
        for (NSInteger iDim=linearTransform.fromDimensions.count-1; iDim >= 0; iDim--)
		{
			GLDimension *fromDim = linearTransform.fromDimensions[iDim];
			GLDimension *toDim = linearTransform.toDimensions[iDim];
			GLMatrixFormat matrixFormat = [linearTransform.matrixFormats[iDim] unsignedIntegerValue];
            
            self.strides[iDim].format = matrixFormat;
			
			// First, note the number of points for this dimension.
			if ( matrixFormat == kGLIdentityMatrixFormat) {
				self.strides[iDim].nPoints = 0;
				self.strides[iDim].nRows = 0;
				self.strides[iDim].nColumns = 0;
			} else if ( matrixFormat == kGLDenseMatrixFormat) {
				self.strides[iDim].nPoints = fromDim.nPoints * toDim.nPoints;
				self.strides[iDim].nRows = toDim.nPoints;
				self.strides[iDim].nColumns = fromDim.nPoints;
			} else if ( matrixFormat == kGLDiagonalMatrixFormat) {
				self.strides[iDim].nPoints = toDim.nPoints;
                self.strides[iDim].nDiagonals = 1;
                self.strides[iDim].nDiagonalPoints = toDim.nPoints;
			} else if ( matrixFormat == kGLTridiagonalMatrixFormat) {
                self.strides[iDim].nPoints = 3*toDim.nPoints;
                self.strides[iDim].nDiagonals = 3;
                self.strides[iDim].nDiagonalPoints = toDim.nPoints;
			}
			
			// Second, determine the stride between elements in this dimension.
			if ( matrixFormat == kGLIdentityMatrixFormat )
			{	// The identity matrix has no data, and therefore no stride.
				self.strides[iDim].stride = 0;
				self.strides[iDim].rowStride = 0;
				self.strides[iDim].columnStride = 0;
			}
			else
			{	// Otherwise, this is the previous (nontrivial) dimension's number of points multiplied by its stride.
				BOOL foundNontrivialDimension = NO;
				for (NSInteger jDim = iDim+1; jDim < linearTransform.fromDimensions.count; jDim++)
				{
					if ( self.strides[jDim].nPoints )
					{
						self.strides[iDim].stride = self.strides[jDim].nPoints * self.strides[jDim].stride;
						
						foundNontrivialDimension = YES;
						break;
					}
				}
				if ( !foundNontrivialDimension) {
					self.strides[iDim].stride = 1;
				}
			}
			
			// Finally, we need to determine how to get to the next row or column within this matrix.
			if (linearTransform.matrixOrder == kGLRowMatrixOrder) {
                if ( !self.strides[iDim].nDiagonals ) {
                    self.strides[iDim].rowStride = self.strides[iDim].stride * self.strides[iDim].nColumns;
                    self.strides[iDim].columnStride = self.strides[iDim].stride;
                } else {
                    self.strides[iDim].diagonalStride = self.strides[iDim].stride * self.strides[iDim].nDiagonalPoints;
                }
			} else if (linearTransform.matrixOrder == kGLColumnMatrixOrder) {
                if ( !self.strides[iDim].nDiagonals ) {
                    self.strides[iDim].rowStride = self.strides[iDim].stride;
                    self.strides[iDim].columnStride = self.strides[iDim].stride * self.strides[iDim].nRows;
                } else {
                    self.strides[iDim].diagonalStride = self.strides[iDim].stride * self.strides[iDim].nDiagonalPoints;
                }
			}
		}
    }
    return self;
}


@end
