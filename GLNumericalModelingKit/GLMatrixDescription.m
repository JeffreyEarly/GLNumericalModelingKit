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

static NSString *GLMatrixDescriptionNDimensionsKey = @"GLMatrixDescriptionNDimensionsKey";
static NSString *GLMatrixDescriptionNPointsKey = @"GLMatrixDescriptionNPointsKey";
static NSString *GLMatrixDescriptionNElementsKey = @"GLMatrixDescriptionNElementsKey";
static NSString *GLMatrixDescriptionNBytesKey = @"GLMatrixDescriptionNBytesKey";
static NSString *GLMatrixDescriptionComplexStrideKey = @"GLMatrixDescriptionComplexStrideKey";
static NSString *GLMatrixDescriptionElementStrideKey = @"GLMatrixDescriptionElementStrideKey";
static NSString *GLMatrixDescriptionDataFormatKey = @"GLMatrixDescriptionDataFormatKey";
static NSString *GLMatrixDescriptionStridesDataKey = @"GLMatrixDescriptionStridesDataKey";

@implementation GLMatrixDescription

- (void)encodeWithCoder:(NSCoder *)coder
{
    [coder encodeObject: @(self.nDimensions) forKey:GLMatrixDescriptionNDimensionsKey];
    [coder encodeObject: @(self.nPoints) forKey:GLMatrixDescriptionNPointsKey];
    [coder encodeObject: @(self.nElements) forKey:GLMatrixDescriptionNElementsKey];
    [coder encodeObject: @(self.nBytes) forKey:GLMatrixDescriptionNBytesKey];
    [coder encodeObject: @(self.complexStride) forKey:GLMatrixDescriptionComplexStrideKey];
    [coder encodeObject: @(self.elementStride) forKey:GLMatrixDescriptionElementStrideKey];
    [coder encodeObject: @(self.dataFormat) forKey:GLMatrixDescriptionDataFormatKey];
    [coder encodeObject: self.stridesData forKey:GLMatrixDescriptionStridesDataKey];
}

- (id)initWithCoder:(NSCoder *)decoder
{
    if ((self=[super init])) {
        _nDimensions = [[decoder decodeObjectForKey: GLMatrixDescriptionNDimensionsKey] unsignedIntegerValue];
        _nPoints = [[decoder decodeObjectForKey: GLMatrixDescriptionNPointsKey] unsignedIntegerValue];
        _nElements = [[decoder decodeObjectForKey: GLMatrixDescriptionNElementsKey] unsignedIntegerValue];
        _nBytes = [[decoder decodeObjectForKey: GLMatrixDescriptionNBytesKey] unsignedIntegerValue];
        _complexStride = [[decoder decodeObjectForKey: GLMatrixDescriptionComplexStrideKey] unsignedIntegerValue];
        _elementStride = [[decoder decodeObjectForKey: GLMatrixDescriptionElementStrideKey] unsignedIntegerValue];
        _dataFormat = [[decoder decodeObjectForKey: GLMatrixDescriptionDataFormatKey] unsignedIntegerValue];
        _stridesData = [decoder decodeObjectForKey: GLMatrixDescriptionStridesDataKey];
        self.strides = self.stridesData.mutableBytes;
    }
    return self;
}

- (GLMatrixDescription *) initWithFunction: (GLFunction *) variable
{
    if ((self=[super init])) {
        self.nDimensions = variable.dimensions.count;
        self.stridesData = [NSMutableData dataWithLength: self.nDimensions*sizeof(GLDataStride)];
        self.strides = self.stridesData.mutableBytes;
        self.dataFormat = variable.dataFormat;
		
		self.nPoints = 1;
        for (NSInteger iDim=variable.dimensions.count-1; iDim >= 0; iDim--)
		{
			GLDimension *dim = variable.dimensions[iDim];
            
			self.nPoints *= dim.nPoints;
			
            self.strides[iDim].matrixFormat = kGLDenseMatrixFormat;
            
            self.strides[iDim].nPoints = dim.nPoints;
            self.strides[iDim].nRows = dim.nPoints;
            self.strides[iDim].nColumns = 1;
			
            if (iDim == variable.dimensions.count-1) {
                self.strides[iDim].stride = variable.dataFormat == kGLInterleavedComplexDataFormat ? 2 : 1;
                self.elementStride = self.strides[iDim].stride;
            } else {
                self.strides[iDim].stride = self.strides[iDim+1].nPoints * self.strides[iDim+1].stride;
            }
  			
			// Finally, we need to determine how to get to the next row or column within this matrix.
            // Thsi shouldn't matter for column vectors.
			self.strides[iDim].rowStride = self.strides[iDim].stride * self.strides[iDim].nColumns;
			self.strides[iDim].columnStride = self.strides[iDim].stride;
		}
		
		if (variable.dataFormat == kGLRealDataFormat) {
			self.complexStride = 0;
			self.nElements = self.nPoints;
		} else if (variable.dataFormat == kGLInterleavedComplexDataFormat) {
			self.complexStride = 1;
			self.nElements = 2*self.nPoints;
		} else {
			self.complexStride = self.nPoints;
			self.nElements = 2*self.nPoints;
		}
		
		self.nBytes = self.nElements*sizeof(GLFloat);
		
		for (NSUInteger iDim=0; iDim<self.nDimensions; iDim++) {
			self.strides[iDim].complexStride = self.complexStride;
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
        self.dataFormat = linearTransform.dataFormat;
		
		self.nPoints = 0;
        for (NSInteger iDim=linearTransform.fromDimensions.count-1; iDim >= 0; iDim--)
		{
			GLDimension *fromDim = linearTransform.fromDimensions[iDim];
			GLDimension *toDim = linearTransform.toDimensions[iDim];
			GLMatrixFormat matrixFormat = [linearTransform.matrixFormats[iDim] unsignedIntegerValue];
            
            self.strides[iDim].matrixFormat = matrixFormat;
			
			if (self.nPoints == 0 && matrixFormat != kGLIdentityMatrixFormat) {
				self.nPoints = 1;
			}
			
			// First, note the number of points for this dimension.
			if ( matrixFormat == kGLIdentityMatrixFormat) {
				self.nPoints *= 1;
				self.strides[iDim].nPoints = 0;
				self.strides[iDim].nRows = 0;
				self.strides[iDim].nColumns = 0;
			} else if ( matrixFormat == kGLDenseMatrixFormat) {
				self.nPoints *= fromDim.nPoints * toDim.nPoints;
				self.strides[iDim].nPoints = fromDim.nPoints * toDim.nPoints;
				self.strides[iDim].nRows = toDim.nPoints;
				self.strides[iDim].nColumns = fromDim.nPoints;
			} else if ( matrixFormat == kGLDiagonalMatrixFormat || matrixFormat == kGLSubdiagonalMatrixFormat || matrixFormat == kGLSuperdiagonalMatrixFormat) {
				self.nPoints *= toDim.nPoints;
				self.strides[iDim].nPoints = toDim.nPoints;
                self.strides[iDim].nDiagonals = 1;
                self.strides[iDim].nDiagonalPoints = toDim.nPoints;
			} else if ( matrixFormat == kGLTridiagonalMatrixFormat) {
				self.nPoints *= 3*toDim.nPoints;
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
					self.strides[iDim].stride = linearTransform.dataFormat == kGLInterleavedComplexDataFormat ? 2 : 1;
                    self.elementStride = self.strides[iDim].stride;
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
		
		if (linearTransform.dataFormat == kGLRealDataFormat) {
			self.complexStride = 0;
			self.nElements = self.nPoints;
		} else if (linearTransform.dataFormat == kGLInterleavedComplexDataFormat) {
			self.complexStride = 1;
			self.nElements = 2*self.nPoints;
		} else {
			self.complexStride = self.nPoints;
			self.nElements = 2*self.nPoints;
		}
		
		self.nBytes = self.nElements*sizeof(GLFloat);
		
		for (NSUInteger iDim=0; iDim<self.nDimensions; iDim++) {
			self.strides[iDim].complexStride = self.complexStride;
		}
    }
    return self;
}

- (BOOL) isEqualToMatrixDescription: (GLMatrixDescription *) otherMatrixDescription
{
	if (self.nDimensions != otherMatrixDescription.nDimensions) return NO;
	if (self.nPoints != otherMatrixDescription.nPoints) return NO;
//	if (self.nElements != otherMatrixDescription.nElements) return NO;
//	if (self.nBytes != otherMatrixDescription.nBytes) return NO;
//	if (self.complexStride != otherMatrixDescription.complexStride) return NO;
//	if (self.dataFormat != otherMatrixDescription.dataFormat) return NO;
	
	BOOL flag = YES;
	for (NSUInteger iDim=0; iDim<self.nDimensions; iDim++) {
		flag &= self.strides[iDim].matrixFormat == otherMatrixDescription.strides[iDim].matrixFormat;
		flag &= self.strides[iDim].nPoints == otherMatrixDescription.strides[iDim].nPoints;
		flag &= self.strides[iDim].nRows == otherMatrixDescription.strides[iDim].nRows;
		flag &= self.strides[iDim].nColumns == otherMatrixDescription.strides[iDim].nColumns;
		flag &= self.strides[iDim].nDiagonals == otherMatrixDescription.strides[iDim].nDiagonals;
		flag &= self.strides[iDim].nDiagonalPoints == otherMatrixDescription.strides[iDim].nDiagonalPoints;
		flag &= self.strides[iDim].stride == otherMatrixDescription.strides[iDim].stride;
		flag &= self.strides[iDim].rowStride == otherMatrixDescription.strides[iDim].rowStride;
		flag &= self.strides[iDim].columnStride == otherMatrixDescription.strides[iDim].columnStride;
		flag &= self.strides[iDim].diagonalStride == otherMatrixDescription.strides[iDim].diagonalStride;
	}
	
	return flag;
}

+ (NSArray *) commonFormatsFromLeft: (NSArray *) left right: (NSArray *) right
{
	NSMutableArray *result = [NSMutableArray array];
	for (NSUInteger i=0; i<left.count; i++) {
		GLMatrixFormat a = [left[i] unsignedIntegerValue];
		GLMatrixFormat b = [right[i] unsignedIntegerValue];
		GLMatrixFormat c;
		
		if (a==b) {
			c = a;
		} else if (a == kGLIdentityMatrixFormat) {
			if (b == kGLSubdiagonalMatrixFormat || b == kGLSuperdiagonalMatrixFormat) {
				c = kGLTridiagonalMatrixFormat;
			} else {
				c = b;
			}
		} else if (a == kGLDenseMatrixFormat) {
			c = kGLDenseMatrixFormat;
		} else if (a == kGLDiagonalMatrixFormat) {
			if (b == kGLIdentityMatrixFormat) {
				c = kGLDiagonalMatrixFormat;
			} else if (b == kGLDenseMatrixFormat) {
				c = kGLDenseMatrixFormat;
			} else {
				c = kGLTridiagonalMatrixFormat;
			}
			
		} else {
			if (b == kGLDenseMatrixFormat) {
				c = kGLDenseMatrixFormat;
			} else {
				c = kGLTridiagonalMatrixFormat;
			}
		}
		
		result[i] = @(c);
	}
	
	return result;
}

@end
