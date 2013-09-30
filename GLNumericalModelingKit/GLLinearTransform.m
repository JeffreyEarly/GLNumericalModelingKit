//
//  GLLinearTransform.m
//  GLNumericalModelingKitCopy
//
//  Created by Jeffrey J. Early on 1/22/13.
//
//

#import "GLLinearTransform.h"
#import "GLEquation.h"
#import "GLDimension.h"

#import "GLLinearTransformationOperations.h"

#include <mach/mach_time.h>

@implementation GLLinearTransform

+ (id) transformOfType: (GLDataFormat) dataFormat withFromDimensions: (NSArray *) fromDims toDimensions: (NSArray *) toDims inFormat: (NSArray *) matrixFormats forEquation: (GLEquation *) equation
{	
	return [[GLLinearTransform alloc] initTransformOfType: dataFormat withFromDimensions: fromDims toDimensions: toDims inFormat: matrixFormats forEquation: equation];
}

- (void) populateDataFromMatrixBlock
{
    if (self.matrixBlock) {
        transformMatrix theMatrixBlock = self.matrixBlock;
        
        // First we need a block to iterate over the inner most matrix.
        NSInteger iDim = self.matrixDescription.nDimensions-1;
        GLDataStride stride = self.matrixDescription.strides[iDim];
        
        void (^innerLoop)( NSUInteger *, NSUInteger *, GLFloat *, NSUInteger );
        
        if (stride.format == kGLDenseMatrixFormat) {
            innerLoop = ^( NSUInteger *row, NSUInteger *col, GLFloat *f, NSUInteger index ) {
                for (NSUInteger i=0; i<stride.nRows; i++) {
                    for (NSUInteger j=0; j<stride.nColumns; j++) {
                        row[iDim] = i;
                        col[iDim] = j;
                        
                        NSUInteger memIndex = index + i*stride.rowStride + j*stride.columnStride;
                        GLFloatComplex value = theMatrixBlock(row, col);
                        
                        if (self.dataFormat == kGLRealDataFormat) {
                            f[memIndex] = creal(value);
                        } else if (self.dataFormat == kGLInterleavedComplexDataFormat) {
                            f[2*memIndex] = creal(value);
                            f[2*memIndex+1] = cimag(value);
                        } else if (self.dataFormat == kGLSplitComplexDataFormat) {
                            f[memIndex] = creal(value);
                            f[memIndex+self.nDataPoints] = cimag(value);
                        }
                    }
                }
            }; 
        }
        
        for (iDim = self.matrixDescription.nDimensions-2; iDim >= 0; iDim--)
        {
            stride = self.matrixDescription.strides[iDim];
            void (^outerLoop)( NSUInteger *, NSUInteger *, GLFloat *, NSUInteger );
            
            if (stride.format == kGLDenseMatrixFormat) {
                outerLoop = ^( NSUInteger *row, NSUInteger *col, GLFloat *f, NSUInteger memIndex ) {
                    for (NSUInteger i=0; i<stride.nRows; i++) {
                        for (NSUInteger j=0; j<stride.nColumns; j++) {
                            row[iDim] = i;
                            col[iDim] = j;
                            
                            memIndex += i*stride.rowStride + j*stride.columnStride;
                            
                            innerLoop( row, col, f, memIndex);
                        }
                    }
                };
                
            }
            
            innerLoop = outerLoop;
        }
        
        // And finally, we can now excute the block.
        NSUInteger *rows = malloc(self.matrixDescription.nDimensions * sizeof(NSUInteger));
        NSUInteger *cols = malloc(self.matrixDescription.nDimensions * sizeof(NSUInteger));
        innerLoop( rows, cols, self.pointerValue, 0);
        free(rows);
        free(cols);
    }
}

+ (id) dftMatrixFromDimension: (GLDimension *) x forEquation: (GLEquation *) equation
{
	GLDimension *k = [[GLDimension alloc] initAsDimension: x transformedToBasis: kGLExponentialBasis strictlyPositive: NO];
	GLLinearTransform *dft = [self transformOfType: kGLSplitComplexDataFormat withFromDimensions: @[x] toDimensions: @[k] inFormat: @[@(kGLDenseMatrixFormat)] forEquation: equation];
	
//	GLFloat *kVal = (GLFloat *) k.data.bytes;
//	GLFloat *xVal = (GLFloat *) x.data.bytes;
//	GLFloat N = x.nPoints;
//	
//	GLSplitComplex f = dft.splitComplex;
//	for (NSUInteger i=0; i<k.nPoints; i++) {
//		for (NSUInteger j=0; j<x.nPoints; j++) {
//			f.realp[i*x.nPoints+j] = cos(2*M_PI*kVal[i]*xVal[j])/N;
//			f.imagp[i*x.nPoints+j] = -sin(2*M_PI*kVal[i]*xVal[j])/N;
//		}
//	}
    
    transformMatrix matrix = ^( NSUInteger *row, NSUInteger *col ) {
        GLFloat *kVal = (GLFloat *) k.data.bytes;
        GLFloat *xVal = (GLFloat *) x.data.bytes;
        GLFloat N = x.nPoints;
        
        GLFloatComplex value = cos(2*M_PI*kVal[row[0]]*xVal[col[0]])/N - I*(sin(2*M_PI*kVal[row[0]]*xVal[col[0]])/N);
        
        return value;
    };
    
    dft.matrixBlock = matrix;
    [dft populateDataFromMatrixBlock];
	
	return dft;
}

+ (id) idftMatrixFromDimension: (GLDimension *) k forEquation: (GLEquation *) equation
{
	GLDimension *x = [[GLDimension alloc] initAsDimension: k transformedToBasis: kGLDeltaBasis strictlyPositive: NO];
	GLLinearTransform *dft = [self transformOfType: kGLSplitComplexDataFormat withFromDimensions: @[k] toDimensions: @[x] inFormat: @[@(kGLDenseMatrixFormat)] forEquation: equation];
	
	GLFloat *kVal = (GLFloat *) k.data.bytes;
	GLFloat *xVal = (GLFloat *) x.data.bytes;
	
	GLSplitComplex f = dft.splitComplex;
	
	for (NSUInteger i=0; i<x.nPoints; i++) {
		for (NSUInteger j=0; j<k.nPoints; j++) {
			f.realp[i*k.nPoints+j] = cos(2*M_PI*kVal[i]*xVal[j]);
			f.imagp[i*k.nPoints+j] = sin(2*M_PI*kVal[i]*xVal[j]);
		}
	}
	
	return dft;
}

+ (id) cosineTransformMatrixFromDimension: (GLDimension *) x forEquation: (GLEquation *) equation
{
	GLDimension *k = [[GLDimension alloc] initAsDimension: x transformedToBasis: kGLCosineHalfShiftBasis strictlyPositive: YES];
	GLLinearTransform *dft = [self transformOfType: kGLRealDataFormat withFromDimensions: @[x] toDimensions: @[k] inFormat: @[@(kGLDenseMatrixFormat)] forEquation: equation];
	
	GLFloat N = x.nPoints;
	
	GLFloat *f = dft.pointerValue;
	for (NSUInteger i=0; i<k.nPoints; i++) {
		for (NSUInteger j=0; j<x.nPoints; j++) {
			f[i*x.nPoints+j] = cos(M_PI*((GLFloat)j+0.5)*((GLFloat)i)/N)/N;
		}
	}
	
	return dft;
}

+ (id) inverseCosineTransformMatrixFromDimension: (GLDimension *) someDimension forEquation: (GLEquation *) equation
{
    GLDimension *x, *k;
    if (someDimension.basisFunction == kGLCosineHalfShiftBasis) {
        k = someDimension;
        x = [[GLDimension alloc] initAsDimension: k transformedToBasis: kGLDeltaBasis strictlyPositive: NO];
    } else if (someDimension.basisFunction == kGLDeltaBasis) {
        x = someDimension;
        k = [[GLDimension alloc] initAsDimension: x transformedToBasis: kGLCosineHalfShiftBasis strictlyPositive: YES];
    } else {
        [NSException exceptionWithName: @"BadDimension" reason:@"I don't understand the words that are coming out of your mouth." userInfo:nil];
    }
    
	GLLinearTransform *dft = [self transformOfType: kGLRealDataFormat withFromDimensions: @[k] toDimensions: @[x] inFormat: @[@(kGLDenseMatrixFormat)] forEquation: equation];
	
	GLFloat N = x.nPoints;
	
	GLFloat *f = dft.pointerValue;
	for (NSUInteger i=0; i<x.nPoints; i++) {
		for (NSUInteger j=0; j<k.nPoints; j++) {
			if (j==0) {
				f[i*k.nPoints+j] = 1.0;
			} else {
				f[i*k.nPoints+j] = 2.0*cos(M_PI*((GLFloat)j)*((GLFloat)i+0.5)/N);
			}
		}
	}
	
	return dft;
}

+ (id) differentiationMatrixFromDimension: (GLDimension *) aDimension forEquation: (GLEquation *) equation
{
	GLLinearTransform *diff;
	if (aDimension.basisFunction == kGLExponentialBasis)
    {
		diff = [GLLinearTransform transformOfType: kGLSplitComplexDataFormat withFromDimensions: @[aDimension] toDimensions: @[aDimension] inFormat: @[@(kGLDenseMatrixFormat)] forEquation: equation];
		
		NSUInteger N = aDimension.nPoints;
		GLSplitComplex f = diff.splitComplex;
		GLFloat *kVal = (GLFloat *) aDimension.data.bytes;
		
		for (NSUInteger i=0; i<N; i++) {
			for (NSUInteger j=0; j<N; j++) {
				if ( i == j ) {
					f.realp[i*N+j] = 0.0;
					f.imagp[i*N+j] = 2*M_PI*kVal[i];
				} else {
					f.realp[i*N+j] = 0.0;
					f.imagp[i*N+j] = 0.0;
				}
			}
		}
	}
	else if (aDimension.basisFunction == kGLCosineHalfShiftBasis)
    {
        GLDimension *transformedDimension = [[GLDimension alloc] initAsDimension: aDimension transformedToBasis: kGLSineHalfShiftBasis strictlyPositive: YES];
		diff = [GLLinearTransform transformOfType: kGLRealDataFormat withFromDimensions: @[aDimension] toDimensions: @[transformedDimension] inFormat: @[@(kGLDenseMatrixFormat)] forEquation: equation];
		
		NSUInteger N = aDimension.nPoints;
		GLFloat *f = diff.pointerValue;
		GLFloat *kVal = (GLFloat *) aDimension.data.bytes;
		
		for (NSInteger i=0; i<N; i++) {
			for (NSInteger j=0; j<N; j++) { // j loops through the 'from' dimension
				if ( i+1 == j ) {
					f[i*N+j] = -2*M_PI*kVal[j];
				} else {
					f[i*N+j] = 0.0;
				}
			}
		}
	}
    else if (aDimension.basisFunction == kGLSineHalfShiftBasis)
    {
        GLDimension *transformedDimension = [[GLDimension alloc] initAsDimension: aDimension transformedToBasis: kGLCosineHalfShiftBasis strictlyPositive: YES];
		diff = [GLLinearTransform transformOfType: kGLRealDataFormat withFromDimensions: @[aDimension] toDimensions: @[transformedDimension] inFormat: @[@(kGLDenseMatrixFormat)] forEquation: equation];
		
		NSUInteger N = aDimension.nPoints;
		GLFloat *f = diff.pointerValue;
		GLFloat *kVal = (GLFloat *) aDimension.data.bytes;
		
		for (NSInteger i=0; i<N; i++) {
			for (NSInteger j=0; j<N; j++) {
				if ( i == j+1 ) {
					f[i*N+j] = 2*M_PI*kVal[j];
				} else {
					f[i*N+j] = 0.0;
				}
			}
		}
	}
    else if (aDimension.basisFunction == kGLChebyshevBasis)
    {
		diff = [GLLinearTransform transformOfType: kGLRealDataFormat withFromDimensions: @[aDimension] toDimensions: @[aDimension] inFormat: @[@(kGLDenseMatrixFormat)] forEquation: equation];
		
		NSUInteger N = aDimension.nPoints;
		GLFloat *f = diff.pointerValue;
		GLFloat *kVal = (GLFloat *) aDimension.data.bytes;
		
		for (NSInteger i=0; i<N; i++) {
			for (NSInteger j=0; j<N; j++) {
				if ( j >= i+1 && (i+j)%2==1 ) {
                    f[i*N+j] = 2*kVal[j];
				} else {
					f[i*N+j] = 0.0;
				}
			}
		}
	}
	
	return diff;
}

- (id) init
{
	[NSException raise: @"BadInitialization" format: @"Cannot initialize GLLinearTransfor with -init method."];
	
	return self;
}

- (id) initTransformOfType: (GLDataFormat) dataFormat withFromDimensions: (NSArray *) fromDims toDimensions: (NSArray *) toDims inFormat: (NSArray *) matrixFormats forEquation: (GLEquation *) theEquation
{
	if (!theEquation || fromDims.count != toDims.count || fromDims.count != matrixFormats.count) {
		NSLog(@"Attempted to initialize GLLinearTransform without an equation or consistent set of dimensions!!!");
		return nil;
	}
	
	if ((self = [super init])) {
		_variableDifferentialOperationMaps = [[NSMutableSet alloc] init];
        _existingOperations = [[NSMutableArray alloc] init];
		
		_pendingOperations = [[NSMutableArray alloc] init];
		_equation = theEquation;
		_uniqueID = mach_absolute_time();
		
		_isFrequencyDomain = 0;
		_isComplex = dataFormat != kGLRealDataFormat;
		_isImaginaryPartZero = dataFormat == kGLRealDataFormat;
		_dataFormat = dataFormat;
        self.matrixFormats = matrixFormats;
		
		self.toDimensions = [NSArray arrayWithArray: toDims];
		self.fromDimensions = [NSArray arrayWithArray: fromDims];
		
		// We loop through the dimensions and allocate enough memory for the variable
		// defined on each dimension.
		_nDataPoints = 0;
		_nDataElements = 0;
		
		for (NSUInteger iDim=0; iDim < fromDims.count; iDim++)
		{
			GLDimension *fromDim = fromDims[iDim];
			GLDimension *toDim = toDims[iDim];
			GLMatrixFormat matrixFormat = [matrixFormats[iDim] unsignedIntegerValue];
			
			if (_nDataPoints == 0 && matrixFormat != kGLIdentityMatrixFormat) {
				_nDataPoints = 1;
				_nDataElements = 1;
			}
			
			if ( matrixFormat == kGLIdentityMatrixFormat) {
				_nDataPoints *= 1;
				_nDataElements *= 1;
			} else if ( matrixFormat == kGLDenseMatrixFormat) {
				_nDataPoints *= fromDim.nPoints * toDim.nPoints;
				_nDataElements *= fromDim.nPoints * toDim.nPoints;
			} else if ( matrixFormat == kGLDiagonalMatrixFormat) {
				_nDataPoints *= toDim.nPoints;
				_nDataElements *= toDim.nPoints;
			} else if ( matrixFormat == kGLTridiagonalMatrixFormat) {
				_nDataPoints *= 3*toDim.nPoints;
				_nDataElements *= 3*toDim.nPoints;
			}
			
			if (fromDim.basisFunction == kGLDeltaBasis) {
                self.realSymmetry[iDim] = @(kGLNoSymmetry);
                self.imaginarySymmetry[iDim] = (dataFormat == kGLRealDataFormat ? @(kGLZeroSymmetry) : @(kGLNoSymmetry));
            } else if (fromDim.basisFunction == kGLCosineBasis || fromDim.basisFunction == kGLCosineHalfShiftBasis) {
                self.realSymmetry[iDim] = @(kGLEvenSymmetry);
                self.imaginarySymmetry[iDim] = @(kGLZeroSymmetry);
            } else if (fromDim.basisFunction == kGLSineBasis || fromDim.basisFunction == kGLSineHalfShiftBasis) {
                self.realSymmetry[iDim] = @(kGLOddSymmetry);
                self.imaginarySymmetry[iDim] = @(kGLZeroSymmetry);
            } else if (fromDim.basisFunction == kGLExponentialBasis ) {
                self.realSymmetry[iDim] = @(kGLNoSymmetry);
                self.imaginarySymmetry[iDim] = (dataFormat == kGLRealDataFormat ? @(kGLZeroSymmetry) : @(kGLNoSymmetry));
            }            
		}
        
        if (dataFormat == kGLSplitComplexDataFormat || dataFormat == kGLInterleavedComplexDataFormat) {
            _nDataElements *= 2;
        }
		
		_dataBytes = _nDataElements*sizeof(GLFloat);
        
        self.matrixDescription = [[GLMatrixDescription alloc] initWithLinearTransform: self];
	}
	
	return self;
}

- (void) setVariableAlongDiagonal: (GLVariable *) diagonalVariable
{
    if (self.matrixDescription.nDimensions == 1)
    {
        if (self.matrixDescription.strides[0].format == kGLDenseMatrixFormat) {
            
            GLFloat *f = self.pointerValue;
            GLFloat *a = diagonalVariable.pointerValue;
            
            NSUInteger rowStride = self.matrixDescription.strides[0].rowStride;
            NSUInteger columnStride = self.matrixDescription.strides[0].columnStride;
            
            for ( NSUInteger i=0; i<self.matrixDescription.strides[0].nRows; i++) {
                for ( NSUInteger j=0; j<self.matrixDescription.strides[0].nColumns; j++) {
                    if (i==j) {
                        f[i*rowStride + j*columnStride] = a[j];
                    } else {
                        f[i*rowStride + j*columnStride] = 0.0;
                    }
                }
            }
            
        }
    }
}

- (void) setVariablesAlongTridiagonal: (NSArray *) tridiagonalVariables
{
    NSUInteger triIndex = NSNotFound;
	NSUInteger firstNonTriIndex = NSNotFound;
	NSUInteger numTriIndices = 0;
    for ( NSNumber *num in self.matrixFormats) {
        if ([num unsignedIntegerValue] == kGLTridiagonalMatrixFormat) {
            triIndex = [self.matrixFormats indexOfObject: num];
			numTriIndices++;
        } else if ([num unsignedIntegerValue] != kGLIdentityMatrixFormat && firstNonTriIndex == NSNotFound ) {
			firstNonTriIndex = [self.matrixFormats indexOfObject: num];
		}
    }
    
    GLFloat *a = [tridiagonalVariables[0] pointerValue];
    GLFloat *b = [tridiagonalVariables[1] pointerValue];
    GLFloat *c = [tridiagonalVariables[2] pointerValue];
    
    GLFloat *d = self.pointerValue;
    NSUInteger elementStride = self.matrixDescription.strides[triIndex].stride;
    NSUInteger diagonalStride = self.matrixDescription.strides[triIndex].diagonalStride;
    NSUInteger m,n;
    for ( NSUInteger i=0; i<[tridiagonalVariables[0] nDataPoints]; i++) {
        m = i%diagonalStride;
        n = i/diagonalStride;
        d[ (3*n+0)*diagonalStride + m*elementStride] = a[i];
        d[ (3*n+1)*diagonalStride + m*elementStride] = b[i];
        d[ (3*n+2)*diagonalStride + m*elementStride] = c[i];
    }
}



- (GLVariable *) tridiagonalSolveWithVector: (GLVariable *)	b
{
	GLTriadiagonalOperation *operation = [[GLTriadiagonalOperation alloc] initWithFirstOperand: self secondOperand: b];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result;
}

- (GLVariable *) solve: (GLVariable *) b
{
	NSUInteger numIdentityIndices = 0;
	NSUInteger numDiagonalIndices = 0;
	NSUInteger numTriIndices = 0;
	NSUInteger numDenseIndices = 0;
	for ( NSNumber *num in self.matrixFormats ) {
        if ([num unsignedIntegerValue] == kGLIdentityMatrixFormat) {
			numIdentityIndices++;
        } else if ([num unsignedIntegerValue] == kGLDiagonalMatrixFormat) {
			numDiagonalIndices++;
        } else if ([num unsignedIntegerValue] == kGLTridiagonalMatrixFormat) {
			numTriIndices++;
        } else if ([num unsignedIntegerValue] == kGLDenseMatrixFormat) {
			numDenseIndices++;
        }
    }
	
	if ( numIdentityIndices && !numDiagonalIndices && !numTriIndices && !numDenseIndices )
	{	// Trivially solution.
		return b;
	}
	else if ( numDiagonalIndices && !numTriIndices && !numDenseIndices )
	{	// Diagonal only
		
	}
	else if ( numTriIndices == 1 && !numDenseIndices )
	{	// A single tridiagonal dimension
		GLTriadiagonalOperation *operation = [[GLTriadiagonalOperation alloc] initWithFirstOperand: self secondOperand: b];
		operation = [self replaceWithExistingOperation: operation];
		return operation.result;
	}
	else if ( !numTriIndices && numDenseIndices == 1 )
	{	// A single dense dimension
		GLDenseMatrixSolver *operation = [[GLDenseMatrixSolver alloc] initWithFirstOperand: self secondOperand: b];
		operation = [self replaceWithExistingOperation: operation];
		return operation.result;
	}
	
	NSString *descrip = [NSString stringWithFormat: @"We can only solve a matrix with one nondiagonal dimension. This matrix contains %lu tridiagonal dimensions and %lu dense dimensions.", numTriIndices, numDenseIndices];
	[NSException exceptionWithName: @"BadFormat" reason:descrip userInfo:nil];
	
	return nil;
}

- (GLVariable *) transform: (GLVariable *) x
{
	NSUInteger numTriIndices = 0;
	for ( NSNumber *num in self.matrixFormats ) {
        if ([num unsignedIntegerValue] == kGLTridiagonalMatrixFormat) {
			numTriIndices++;
        }
    }
		
	if (numTriIndices) {
		GLTriadiagonalTransformOperation *operation = [[GLTriadiagonalTransformOperation alloc] initWithFirstOperand: self secondOperand: x];
		operation = [self replaceWithExistingOperation: operation];
		return operation.result;
	} else {
		GLDenseMatrixTransformOperation *operation = [[GLDenseMatrixTransformOperation alloc] initWithFirstOperand: self secondOperand: x];
		operation = [self replaceWithExistingOperation: operation];
		return operation.result;
	}

}

- (GLVariable *) times: (GLVariable *) otherVariable
{
	return [self multiply: otherVariable];
}

- (GLVariable *) multiply: (GLVariable *) otherVariable
{
	GLMatrixMatrixMultiplicationOperation *operation = [[GLMatrixMatrixMultiplicationOperation alloc] initWithFirstOperand: self secondOperand: otherVariable];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result;
}

- (GLLinearTransform *) inverse
{
    GLMatrixInversionOperation *operation = [[GLMatrixInversionOperation alloc] initWithOperand: self];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result;
}

- (GLLinearTransform *) plus: (GLLinearTransform *) otherVariable
{
	GLLinearTransformAdditionOperation *operation = [[GLLinearTransformAdditionOperation alloc] initWithFirstOperand: self secondOperand: otherVariable];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result;
}

- (NSString *) matrixDescriptionString
{
	NSUInteger n = self.matrixDescription.strides[0].rowStride;
	n = n==0?1:n;
	NSMutableString *descrip = [NSMutableString string];
	
	GLFloat max, min;
	vGL_maxv( self.data.mutableBytes, 1, &max, self.nDataElements);
	vGL_minv( self.data.mutableBytes, 1, &min, self.nDataElements);
	
	if ( fabs(min) > max) {
		max = fabs(min);
	}
	
	GLFloat divisor = pow(10, floor(log10(max)));
	if ( divisor == 0.0) divisor = 1;
	
	if (0 && self.dataFormat == kGLSplitComplexDataFormat)
	{
		GLSplitComplex splitComplex = self.splitComplex;
		[descrip appendFormat: @"%f * ", divisor];
        //		for (NSUInteger i=0; i<self.nDataPoints; i++)
        //		{
        //			if ( i % n == 0 ) {
        //				[descrip appendFormat: @"\n"];
        //			}
        //			[descrip appendFormat: @"%1.1f ", sqrt(fabs(splitComplex.realp[i] * splitComplex.realp[i] - splitComplex.imagp[i] * splitComplex.imagp[i]))/divisor];
        //		}
		
		for (NSUInteger i=0; i<self.nDataPoints; i++)
		{
			if ( i % n == 0 ) {
				[descrip appendFormat: @"\n"];
			}
			[descrip appendFormat: @"%+1.1f ", splitComplex.realp[i]/divisor];
		}
		
		[descrip appendFormat: @" imagp \n"];
		for (NSUInteger i=0; i<self.nDataPoints; i++)
		{
			if ( i % n == 0 ) {
				[descrip appendFormat: @"\n"];
			}
			[descrip appendFormat: @"%+1.1f ", splitComplex.imagp[i]/divisor];
		}
	}
    if ( self.dimensions.count == 3)
    {
        NSUInteger m = [self.dimensions[2] nPoints] * [self.dimensions[1] nPoints];
        GLFloat *f = self.pointerValue;
		[descrip appendFormat: @"%g * ", divisor];
		for (NSUInteger i=0; i<self.nDataElements; i++)
		{
			if ( i % m == 0 ) {
				[descrip appendFormat: @"\n"];
			}
            if ( i % n == 0 ) {
				[descrip appendFormat: @"\n"];
			}
            
			[descrip appendFormat: @"%+1.1f ", f[i]/divisor];
		}
    }
	else
	{
		GLFloat *f = self.pointerValue;
		[descrip appendFormat: @"%g * ", divisor];
		for (NSUInteger i=0; i<self.nDataElements; i++)
		{
			if ( i % n == 0 ) {
				[descrip appendFormat: @"\n"];
			}
			[descrip appendFormat: @"%+1.1f ", f[i]/divisor];
		}
	}
	
	return descrip;
}

@end
