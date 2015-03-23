//
//  GLBasisTransformationOperations.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 6/14/12.
//
//

#import "GLBasisTransformationOperations.h"
#import <fftw3.h>
#import "GLMemoryPool.h"

/************************************************/
/*		GLBasisTransformOperation               */
/************************************************/

#pragma mark -
#pragma mark GLBasisTransformOperation
#pragma mark

// We need to partition the transformation in

// We support three situations
// 1. Spatial domain, complex functions.
// 2. Spatial domain, real functions.
// 3. Frequency domain, complex functions.
// 4. Frequency domain, real functions.
//
// 1. complex-to-complex
// 2. real-to-real, followed by real-to-complex
// 3. if hermitian {complex-to-real, real-to-real}, else {complex-to-complex}
// 4. real-to-real

@implementation GLBasisTransformOperation

static BOOL _shouldAntiAlias = NO;

+ (NSString *) wisdomFilePath
{
    NSFileManager *fileManager = [[NSFileManager alloc] init];
    NSArray *urls = [fileManager URLsForDirectory: NSApplicationSupportDirectory inDomains: NSUserDomainMask];
    
    NSURL *path;
    if (urls.count) {
        path = [urls.lastObject URLByAppendingPathComponent:@"com.earlyinnovations.GLNumericalModelingKit"];
        [fileManager createDirectoryAtURL: path withIntermediateDirectories: YES attributes: nil error: nil];
        path = [path URLByAppendingPathComponent: @"wisdom"];
    }
    return path.path;
}

+ (void) readWisdom
{
    NSUInteger cores = [[NSProcessInfo processInfo] activeProcessorCount]/2;
    if (cores < 1) cores=1;
    vGL_fftw_plan_with_nthreads( (int) cores);
    
    NSString *path = [self wisdomFilePath];
    
     NSFileManager *fileManager = [[NSFileManager alloc] init];
    if ([fileManager fileExistsAtPath: path])
    {
        if (path) {
            vGL_fftw_import_wisdom_from_filename( [path cStringUsingEncoding: NSASCIIStringEncoding] );
        }
    }
}

+ (void) saveWisdom
{
    NSString *path = [self wisdomFilePath];
    if (path) {
        if (!vGL_fftw_export_wisdom_to_filename( [path cStringUsingEncoding: NSASCIIStringEncoding] )) {
            NSLog(@"Failed to save wisdom");
        }
    }
}

+ (void) setShouldAntialias: (BOOL) flag
{
	_shouldAntiAlias = flag;
}
+ (BOOL) shouldAntialias {
	return _shouldAntiAlias;
}

// antialiasingFilter = 1 if k < (2/3)K_max, 0 otherwise
+ (GLFunction *) antialiasingFilterFromDimensions: (NSArray *) dimensions forEquation: (GLEquation *) equation
{
	GLFloat minNyquist;
	
	GLFunction *bigK;
	for (GLDimension *dim in dimensions) {
		if (dim.basisFunction == kGLDeltaBasis) {
			continue;
		}
		GLFunction *k = [GLFunction functionOfRealTypeFromDimension: dim withDimensions: dimensions forEquation: equation];
		GLFloat nyquist = dim.domainMin + dim.domainLength;
		if (!bigK) {
			bigK = [k multiply: k];
			minNyquist = nyquist;
		} else {
			bigK = [bigK plus: [k multiply: k]];
			if (nyquist < minNyquist) {
				minNyquist = nyquist;
			}
		}
	}
	if (!bigK) {
		[NSException raise: @"Dimension Issue" format:@"Requesting an antialiasing filter for a variable with no wavenumber dimensions"];
	}
	
	bigK = [bigK sqrt];
	
	GLFloat wavenumberCutoff = 2*minNyquist/3;
	
	
	// expf( -alpha * powf( (k-max)/(k-cutoff), p) );
	GLFunction *filter = [GLFunction functionOfRealTypeWithDimensions: dimensions forEquation: equation];
	
	[bigK solve];
	
	GLFloat *kk = [bigK pointerValue];
	GLFloat *f = [filter pointerValue];
	for (NSUInteger i=0; i<filter.nDataPoints; i++) {
		if ( kk[i] < wavenumberCutoff ) {
			f[i] = 1.0;
		} else {
			f[i] = 0.0;
		}
	}
	
	return filter;
}

+ (id) basisTransformationWithOperand: (GLFunction *) variable destinationBasis: (NSArray *) toBasis
{
    // If only one basis was specified, assume it was intended for all dimensions
    if ( toBasis.count == 1 && variable.dimensions.count > 1) {
        NSMutableArray *array = [NSMutableArray array];
        for (NSUInteger i=0; i < variable.dimensions.count; i++) {
            [array addObject: toBasis.lastObject];
        }
        toBasis = array;
    }
    
	if (variable.dimensions.count != toBasis.count) {
		[NSException raise: @"DimensionsNotEqualException" format: @"The destination basis must be specified for all dimensions."];
	}
    
    // Basis array of the variable in its current state.
    NSMutableArray *fromBasis = [NSMutableArray array];
    for (GLDimension *aDim in variable.dimensions) {
		[fromBasis addObject: @(aDim.basisFunction)];
    }
	
	// Check to see if a transformation is even necessary
	BOOL finishedTransformation = YES;
	for (NSInteger i=0; i < toBasis.count; i++) {
		finishedTransformation &= [toBasis[i] isEqualTo: fromBasis[i]];
	}
	if (finishedTransformation) {
		// No transformation necessary
		return nil;
	}
	
    //
    //  1. Look for real cosine and sine transformations forward
    //
    NSUInteger transformationRank = 0;
    NSMutableArray *transformationBasis = [fromBasis mutableCopy];
    if (variable.isPurelyReal)
    {   // Any forward r2r transforms?
        for (GLDimension *aDim in variable.dimensions) {
            NSUInteger idx = [variable.dimensions indexOfObject: aDim];
            GLBasisFunction destinationBasis = [toBasis[idx] unsignedIntegerValue];
            
            if (aDim.basisFunction == kGLDeltaBasis) {
                if (destinationBasis == kGLChebyshevBasis || destinationBasis == kGLDiscreteCosineTransformIBasis || destinationBasis == kGLCosineBasis || destinationBasis == kGLDiscreteSineTransformIBasis || destinationBasis == kGLSineBasis) {
                    transformationRank++;
                    transformationBasis[idx] = @(destinationBasis);
                }
            }
        }
    }
    
    if (transformationRank) {
        GLRealToRealTransformOperation *operation = [[GLRealToRealTransformOperation alloc] initWithOperand: variable destinationBasis: transformationBasis];
        operation = [variable replaceWithExistingOperation: operation];
        BOOL finishedTransformation = YES;
        for (NSInteger i=0; i < toBasis.count; i++) {
            finishedTransformation &= [toBasis[i] isEqualTo: transformationBasis[i]];
        }
        if (finishedTransformation) {
            return operation;
        } else {
            variable = operation.result[0];
        }
    }
    
    //
    //  2. Look for forward real exponential transformations
    //
    fromBasis = [transformationBasis mutableCopy];
    transformationRank = 0;
    if (variable.isPurelyReal)
    { 
        for (GLDimension *aDim in variable.dimensions) {
            NSUInteger idx = [variable.dimensions indexOfObject: aDim];
            GLBasisFunction destinationBasis = [toBasis[idx] unsignedIntegerValue];
            
            if (aDim.basisFunction == kGLDeltaBasis) {
                if (destinationBasis == kGLExponentialBasis) {
                    transformationRank++;
                    transformationBasis[idx] = @(destinationBasis);
                }
            }
        }
    }
    
    if (transformationRank) {
        GLRealToComplexTransformOperation *operation = [[GLRealToComplexTransformOperation alloc] initWithOperand: variable destinationBasis: transformationBasis];
        operation = [variable replaceWithExistingOperation: operation];
        BOOL finishedTransformation = YES;
        for (NSInteger i=0; i < toBasis.count; i++) {
            finishedTransformation &= [toBasis[i] isEqualTo: transformationBasis[i]];
        }
        if (finishedTransformation) {
            return operation;
        } else {
            variable = operation.result[0];
        }
    }
    
    //
    //  3. Look forward complex transformations
    //
    fromBasis = [transformationBasis mutableCopy];
    transformationRank = 0;
    for (GLDimension *aDim in variable.dimensions) {
        NSUInteger idx = [variable.dimensions indexOfObject: aDim];
        GLBasisFunction destinationBasis = [toBasis[idx] unsignedIntegerValue];
        
        if (aDim.basisFunction == kGLDeltaBasis) {
            if (destinationBasis == kGLExponentialBasis) {
                transformationRank++;
                transformationBasis[idx] = @(destinationBasis);
            }
        }
    }
    
    if (transformationRank) {
        GLComplexToComplexTransformOperation *operation = [[GLComplexToComplexTransformOperation alloc] initWithOperand: variable destinationBasis: transformationBasis];
        operation = [variable replaceWithExistingOperation: operation];
        BOOL finishedTransformation = YES;
        for (NSInteger i=0; i < toBasis.count; i++) {
            finishedTransformation &= [toBasis[i] isEqualTo: transformationBasis[i]];
        }
        if (finishedTransformation) {
            return operation;
        } else {
            variable = operation.result[0];
        }
    }
    
	//
	//	At this point is appears as if we're transforming back to the grid, so we need to anti-alias
	//
	if ([GLBasisTransformOperation shouldAntialias]) {
		GLFunction *aa = [self antialiasingFilterFromDimensions: variable.dimensions forEquation: variable.equation];
		variable = [variable multiply: aa];
	}
	
    //
    //  4. Look for inverse real exponential transformations
    //
    fromBasis = [transformationBasis mutableCopy];
    transformationRank = 0;
    GLDimension *lastDim;
    for (GLDimension *aDim in variable.dimensions) {
        NSUInteger idx = [variable.dimensions indexOfObject: aDim];
        GLBasisFunction destinationBasis = [toBasis[idx] unsignedIntegerValue];
        
        if (aDim.basisFunction == kGLExponentialBasis) {
            if (destinationBasis == kGLDeltaBasis) {
                transformationRank++;
                transformationBasis[idx] = @(destinationBasis);
                lastDim = aDim;
            }
        }
    }
    
    // Actually need to test if the variable is hermitian
    if (transformationRank && lastDim.isStrictlyPositive) {
        GLRealToComplexTransformOperation *operation = [[GLRealToComplexTransformOperation alloc] initWithOperand: variable destinationBasis: transformationBasis];
        operation = [variable replaceWithExistingOperation: operation];
        BOOL finishedTransformation = YES;
        for (NSInteger i=0; i < toBasis.count; i++) {
            finishedTransformation &= [toBasis[i] isEqualTo: transformationBasis[i]];
        }
        if (finishedTransformation) {
            return operation;
        } else {
            variable = operation.result[0];
        }
    }
    
	//
    //  5. Look for inverse cosine and sine transformations
    //
    transformationRank = 0;
    fromBasis = [transformationBasis mutableCopy];
    if (variable.isPurelyReal)
    {   // Any forward r2r transforms?
        for (GLDimension *aDim in variable.dimensions) {
            NSUInteger idx = [variable.dimensions indexOfObject: aDim];
            GLBasisFunction destinationBasis = [toBasis[idx] unsignedIntegerValue];
            
            if (aDim.basisFunction == kGLChebyshevBasis || aDim.basisFunction == kGLDiscreteCosineTransformIBasis || aDim.basisFunction == kGLCosineBasis || aDim.basisFunction == kGLDiscreteSineTransformIBasis || aDim.basisFunction == kGLSineBasis){
                if (destinationBasis == kGLDeltaBasis) {
                    transformationRank++;
                    transformationBasis[idx] = @(destinationBasis);
                }
            }
        }
    }
    
    if (transformationRank) {
        GLRealToRealTransformOperation *operation = [[GLRealToRealTransformOperation alloc] initWithOperand: variable destinationBasis: transformationBasis];
        operation = [variable replaceWithExistingOperation: operation];
        BOOL finishedTransformation = YES;
        for (NSInteger i=0; i < toBasis.count; i++) {
            finishedTransformation &= [toBasis[i] isEqualTo: transformationBasis[i]];
        }
        if (finishedTransformation) {
            return operation;
        } else {
            variable = operation.result[0];
        }
    }
    
    //
    //  6. Look for inverse complex transformations
    //
    fromBasis = [transformationBasis mutableCopy];
    transformationRank = 0;
    for (GLDimension *aDim in variable.dimensions) {
        NSUInteger idx = [variable.dimensions indexOfObject: aDim];
        GLBasisFunction destinationBasis = [toBasis[idx] unsignedIntegerValue];
        
        if (aDim.basisFunction == kGLExponentialBasis) {
            if (destinationBasis == kGLDeltaBasis) {
                transformationRank++;
                transformationBasis[idx] = @(destinationBasis);
            }
        }
    }
    
    if (transformationRank) {
        GLComplexToComplexTransformOperation *operation = [[GLComplexToComplexTransformOperation alloc] initWithOperand: variable destinationBasis: transformationBasis];
        operation = [variable replaceWithExistingOperation: operation];
        BOOL finishedTransformation = YES;
        for (NSInteger i=0; i < toBasis.count; i++) {
            finishedTransformation &= [toBasis[i] isEqualTo: transformationBasis[i]];
        }
        if (finishedTransformation) {
            return operation;
        } else {
            variable = operation.result[0];
        }
    }
    
	[NSException raise: @"TransformationFailure" format: @"Unable to find the appropriate combination of transformations."];
    return nil;
}

- (BOOL) isEqualToOperation: (id) otherOperation {
    if ( ![super isEqualToOperation: otherOperation] )  {
        return NO;
    }
    
    GLBasisTransformOperation * op = otherOperation;
    if (self.toBasis.count != op.toBasis.count) {
        return NO;
    }
    if (self.fromBasis.count != op.fromBasis.count) {
        return NO;
    }
    for (NSUInteger i=0; i<self.fromBasis.count; i++) {
        if ( [self.fromBasis[i] unsignedIntegerValue] != [op.fromBasis[i] unsignedIntegerValue] ) {
            return NO;
        }
    }
    for (NSUInteger i=0; i<self.toBasis.count; i++) {
        if ( [self.toBasis[i] unsignedIntegerValue] != [op.toBasis[i] unsignedIntegerValue] ) {
            return NO;
        }
    }
    
    return YES;
}

- (NSString *) graphvisDescription {
    NSMutableString *basisDescription = [NSMutableString stringWithFormat:@"("];
    for ( NSNumber *basis in self.fromBasis) {
        if (basis.unsignedIntegerValue == kGLDeltaBasis) {
            [basisDescription appendString: @"delta"];
        } else if (basis.unsignedIntegerValue == kGLExponentialBasis) {
            [basisDescription appendString: @"exponential"];
        } else if (basis.unsignedIntegerValue == kGLDiscreteCosineTransformIBasis) {
            [basisDescription appendString: @"cosine"];
        } else if (basis.unsignedIntegerValue == kGLDiscreteSineTransformIBasis) {
            [basisDescription appendString: @"sine"];
        } else if (basis.unsignedIntegerValue == kGLCosineBasis) {
            [basisDescription appendString: @"cosine"];
        } else if (basis.unsignedIntegerValue == kGLSineBasis) {
            [basisDescription appendString: @"sine"];
        } else if (basis.unsignedIntegerValue == kGLChebyshevBasis) {
            [basisDescription appendString: @"chebyshev"];
        }
        [basisDescription appendString: @", "];
    }
    [basisDescription deleteCharactersInRange: NSMakeRange(basisDescription.length-2, 2)];
    [basisDescription appendString: @") -> ("];
    for ( NSNumber *basis in self.toBasis) {
        if (basis.unsignedIntegerValue == kGLDeltaBasis) {
            [basisDescription appendString: @"delta"];
        } else if (basis.unsignedIntegerValue == kGLExponentialBasis) {
            [basisDescription appendString: @"exponential"];
        } else if (basis.unsignedIntegerValue == kGLDiscreteCosineTransformIBasis) {
            [basisDescription appendString: @"cosine"];
        } else if (basis.unsignedIntegerValue == kGLDiscreteSineTransformIBasis) {
            [basisDescription appendString: @"sine"];
        } else if (basis.unsignedIntegerValue == kGLCosineBasis) {
            [basisDescription appendString: @"cosine"];
        } else if (basis.unsignedIntegerValue == kGLSineBasis) {
            [basisDescription appendString: @"sine"];
        } else if (basis.unsignedIntegerValue == kGLChebyshevBasis) {
            [basisDescription appendString: @"chebyshev"];
        }
        [basisDescription appendString: @", "];
    }
    [basisDescription deleteCharactersInRange: NSMakeRange(basisDescription.length-2, 2)];
    [basisDescription appendString: @")"];
    
    return basisDescription;
}

- (NSString *) description
{
    NSMutableString *basisDescription = [NSMutableString stringWithFormat:@"("];
    for ( NSNumber *basis in self.fromBasis) {
        if (basis.unsignedIntegerValue == kGLDeltaBasis) {
            [basisDescription appendString: @"delta"];
        } else if (basis.unsignedIntegerValue == kGLExponentialBasis) {
            [basisDescription appendString: @"exponential"];
        } else if (basis.unsignedIntegerValue == kGLDiscreteCosineTransformIBasis) {
            [basisDescription appendString: @"cosine"];
        } else if (basis.unsignedIntegerValue == kGLDiscreteSineTransformIBasis) {
            [basisDescription appendString: @"sine"];
        } else if (basis.unsignedIntegerValue == kGLCosineBasis) {
            [basisDescription appendString: @"cosine"];
        } else if (basis.unsignedIntegerValue == kGLSineBasis) {
            [basisDescription appendString: @"sine"];
        } else if (basis.unsignedIntegerValue == kGLChebyshevBasis) {
            [basisDescription appendString: @"chebyshev"];
        }
        [basisDescription appendString: @", "];
    }
    [basisDescription deleteCharactersInRange: NSMakeRange(basisDescription.length-2, 2)];
    [basisDescription appendString: @") -> ("];
    for ( NSNumber *basis in self.toBasis) {
        if (basis.unsignedIntegerValue == kGLDeltaBasis) {
            [basisDescription appendString: @"delta"];
        } else if (basis.unsignedIntegerValue == kGLExponentialBasis) {
            [basisDescription appendString: @"exponential"];
        } else if (basis.unsignedIntegerValue == kGLDiscreteCosineTransformIBasis) {
            [basisDescription appendString: @"cosine"];
        } else if (basis.unsignedIntegerValue == kGLDiscreteSineTransformIBasis) {
            [basisDescription appendString: @"sine"];
        } else if (basis.unsignedIntegerValue == kGLCosineBasis) {
            [basisDescription appendString: @"cosine"];
        } else if (basis.unsignedIntegerValue == kGLSineBasis) {
            [basisDescription appendString: @"sine"];
        } else if (basis.unsignedIntegerValue == kGLChebyshevBasis) {
            [basisDescription appendString: @"chebyshev"];
        }
        [basisDescription appendString: @", "];
    }
    [basisDescription deleteCharactersInRange: NSMakeRange(basisDescription.length-2, 2)];
    [basisDescription appendString: @")"];
    
    if (((GLFunction *)self.operand[0]).name) {
        return [NSString stringWithFormat: @"%@\t<0x%lx> (%@) %@", NSStringFromClass([self class]), (NSUInteger) self, ((GLFunction *)self.operand[0]).name, basisDescription];
    } else {
        return [NSString stringWithFormat: @"%@\t<0x%lx> %@", NSStringFromClass([self class]), (NSUInteger) self, basisDescription];
    }
}

@end

/************************************************/
/*		GLRealToRealTransformOperation          */
/************************************************/

#pragma mark -
#pragma mark GLRealToRealTransformOperation
#pragma mark

@implementation GLRealToRealTransformOperation

- (id) initWithOperand: (GLFunction *) variable destinationBasis: (NSArray *) toBasis
{
    if (!variable.isPurelyReal) {
		[NSException raise: @"VariableNotRealException" format: @"A real-to-real transformation can only act on a real variable."];
	}
    
	if (variable.dimensions.count != toBasis.count) {
		[NSException raise: @"DimensionsNotEqualException" format: @"The destination basis must be specified for all dimensions."];
	}
    
    // Basis array of the variable in its current state.
    NSMutableArray *initialBasis = [NSMutableArray array];
    for (GLDimension *aDim in variable.dimensions) {
		[initialBasis addObject: @(aDim.basisFunction)];
    }
	
    // Look for cosine and sine transformations, *only*
    NSUInteger finalRank = 0;
    NSMutableArray *finalBasis = [initialBasis mutableCopy];
    NSMutableArray *finalDimensions = [NSMutableArray array];
    for (GLDimension *aDim in variable.dimensions) {
        NSUInteger idx = [variable.dimensions indexOfObject: aDim];
        GLBasisFunction destinationBasis = [toBasis[idx] unsignedIntegerValue];
        
        if (aDim.basisFunction == destinationBasis)
        {   // This dimension isn't getting transformed
            [finalDimensions addObject: aDim];
        } else if (aDim.basisFunction == kGLDeltaBasis && (destinationBasis == kGLChebyshevBasis || destinationBasis == kGLDiscreteCosineTransformIBasis || destinationBasis == kGLCosineBasis || destinationBasis == kGLDiscreteSineTransformIBasis || destinationBasis == kGLSineBasis)) {
            finalRank++;
            finalBasis[idx] = @(destinationBasis);
            [finalDimensions addObject: [[GLDimension alloc] initAsDimension: aDim transformedToBasis: destinationBasis strictlyPositive: YES]];
        } else if (destinationBasis == kGLDeltaBasis && (aDim.basisFunction == kGLChebyshevBasis || aDim.basisFunction == kGLDiscreteCosineTransformIBasis || aDim.basisFunction == kGLCosineBasis || aDim.basisFunction == kGLDiscreteSineTransformIBasis || aDim.basisFunction == kGLSineBasis)) {
            finalRank++;
            finalBasis[idx] = @(destinationBasis);
            [finalDimensions addObject: [[GLDimension alloc] initAsDimension: aDim transformedToBasis: destinationBasis strictlyPositive: YES]];
        } else {
            [NSException raise: @"InvalidDestinationBasis" format: @"Invalid destination basis for a real-to-real transform."];
        }
    }
    
    // Okay, we're good to go!
    GLFunction *resultVariable = [GLFunction functionOfRealTypeWithDimensions: finalDimensions forEquation: variable.equation];
    resultVariable.name = variable.name;
	if (( self = [super initWithResult: @[resultVariable] operand: @[variable]] ))
	{
        self.fromBasis = initialBasis;
        self.toBasis = finalBasis;
        
        int rank = 0;
        int howmany_rank = 0;
        
        // First compute the strides between dimensions
        vGL_fftw_iodim *alldims = malloc(finalDimensions.count*sizeof(vGL_fftw_iodim));
        for (NSInteger i=finalDimensions.count-1; i >= 0; i--)
        {
            GLDimension *aDim = [finalDimensions objectAtIndex: i];
            
            if (i == finalDimensions.count-1 ) {
                alldims[i].n = (int) aDim.nPoints;
                alldims[i].is = 1;
                alldims[i].os = 1;
            } else {
                alldims[i].n = (int) aDim.nPoints;
                alldims[i].is = alldims[i+1].is*alldims[i+1].n;
                alldims[i].os = alldims[i+1].os*alldims[i+1].n;
            }
            
            if ( [initialBasis[i] isEqualTo: finalBasis[i]] ) {
                howmany_rank++;
            } else {
                rank++;
            }
        }
        
        // Now sort them into transformed and untransformed dimensions
        vGL_fftw_r2r_kind *transformKind = malloc(rank*sizeof(vGL_fftw_r2r_kind));
        vGL_fftw_iodim *dims = malloc(rank*sizeof(vGL_fftw_iodim));
        vGL_fftw_iodim *howmany_dims = malloc(howmany_rank*sizeof(vGL_fftw_iodim));
        
        NSUInteger dataPoints = 1;
        NSUInteger transformedDataPoints = 1;
        NSUInteger iDims = 0;
        NSUInteger iHowManyDims = 0;
        GLFloat scaleFactor = 1.0;
        for (NSUInteger i=0; i<finalDimensions.count; i++)
        {
            GLDimension *aDim = [finalDimensions objectAtIndex: i];
            dataPoints *= aDim.nPoints;
            
            if ( [initialBasis[i] isEqualTo: finalBasis[i]] )
            {   // This dimension does NOT get transformed
                howmany_dims[iHowManyDims].n =alldims[i].n;
                howmany_dims[iHowManyDims].is =alldims[i].is;
                howmany_dims[iHowManyDims].os =alldims[i].os;
                iHowManyDims++;
            }
            else
            {   // This dimension DOES get transformed
                dims[iDims].n =alldims[i].n;
                dims[iDims].is =alldims[i].is;
                dims[iDims].os =alldims[i].os;
                
                GLBasisFunction fromBasis = [[initialBasis objectAtIndex: i] unsignedIntegerValue];
                GLBasisFunction toBasis = [[finalBasis objectAtIndex: i] unsignedIntegerValue];
                
                if (fromBasis == kGLDeltaBasis && toBasis == kGLExponentialBasis) {
                    (transformKind)[iDims] = FFTW_R2HC;
                    scaleFactor *= 1.0 / ( (GLFloat) alldims[i].n );
                } else if (fromBasis == kGLDeltaBasis && toBasis == kGLDiscreteCosineTransformIBasis) {
                    (transformKind)[iDims] = FFTW_REDFT00;
                    scaleFactor *= 0.5 / ( (GLFloat) alldims[i].n - 1 );
                } else if (fromBasis == kGLDeltaBasis && toBasis == kGLCosineBasis) {
                    (transformKind)[iDims] = FFTW_REDFT10;
                    scaleFactor *= 0.5 / ( (GLFloat) alldims[i].n );
                } else if (fromBasis == kGLDeltaBasis && toBasis == kGLDiscreteSineTransformIBasis) {
                    (transformKind)[iDims] = FFTW_RODFT00;
                    scaleFactor *= 0.5 / ( (GLFloat) alldims[i].n + 1 );
                } else if (fromBasis == kGLDeltaBasis && toBasis == kGLSineBasis) {
                    (transformKind)[iDims] = FFTW_RODFT10;
                    scaleFactor *= 0.5 / ( (GLFloat) alldims[i].n );
                } else if (fromBasis == kGLDeltaBasis && toBasis == kGLChebyshevBasis) {
                    (transformKind)[iDims] = FFTW_REDFT00;
                    scaleFactor *= 1.0 / ( (GLFloat) alldims[i].n - 1 );
                }  else if (fromBasis == kGLExponentialBasis && toBasis == kGLDeltaBasis) {
                    (transformKind)[iDims] = FFTW_HC2R;
                } else if (fromBasis == kGLDiscreteCosineTransformIBasis && toBasis == kGLDeltaBasis) {
                    (transformKind)[iDims] = FFTW_REDFT00;
                } else if (fromBasis == kGLCosineBasis && toBasis == kGLDeltaBasis) {
                    (transformKind)[iDims] = FFTW_REDFT01;
                } else if (fromBasis == kGLDiscreteSineTransformIBasis && toBasis == kGLDeltaBasis) {
                    (transformKind)[iDims] = FFTW_RODFT00;
                } else if (fromBasis == kGLSineBasis && toBasis == kGLDeltaBasis) {
                    (transformKind)[iDims] = FFTW_RODFT01;
                } else if (fromBasis == kGLChebyshevBasis && toBasis == kGLDeltaBasis) {
                    (transformKind)[iDims] = FFTW_REDFT00;
					scaleFactor *= 0.5;
                }
                
                iDims++;
                
                transformedDataPoints *= aDim.nPoints;
                
            }
        }
        
        NSMutableData *tempf = [[GLMemoryPool sharedMemoryPool] dataWithLength: variable.nDataElements*sizeof(GLFloat)];
        NSMutableData *tempfbar = [[GLMemoryPool sharedMemoryPool] dataWithLength: resultVariable.nDataElements*sizeof(GLFloat)];
//        fftwf_plan plan = fftwf_plan_guru_r2r(rank, dims, howmany_rank, howmany_dims, tempf.mutableBytes, tempfbar.mutableBytes, transformKind, FFTW_PATIENT);
        
        [GLBasisTransformOperation readWisdom];
        vGL_fftw_plan plan = vGL_fftw_plan_guru_r2r(rank, dims, howmany_rank, howmany_dims, tempf.mutableBytes, tempfbar.mutableBytes, transformKind, FFTW_PATIENT);
        [GLBasisTransformOperation saveWisdom];
        
        // Clean-up
        [[GLMemoryPool sharedMemoryPool] returnData: tempf];
        [[GLMemoryPool sharedMemoryPool] returnData: tempfbar];
        free(alldims);
        free(transformKind);
        free(dims);
        free(howmany_dims);
        
        // Note that we leak the plan.
        if (scaleFactor == 1.0) {
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				NSMutableData *result = resultArray[0];
				NSMutableData *operand = operandArray[0];
				
                GLFloat *f = (void *) operand.bytes;
                GLFloat *fbar = (void *) result.mutableBytes;
                
                vGL_fftw_execute_r2r( plan, f, fbar );
            };
        } else {
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				NSMutableData *result = resultArray[0];
				NSMutableData *operand = operandArray[0];
				
                GLFloat *f = (void *) operand.bytes;
                GLFloat *fbar = (void *) result.mutableBytes;
                
                vGL_fftw_execute_r2r( plan, f, fbar );
                
                GLFloat inversesize = scaleFactor;
                vGL_vsmul( fbar, 1, &inversesize, fbar, 1, dataPoints );			
            };
        }
		
		

    }
    
    return self;
}

@end

/************************************************/
/*		GLRealToComplexTransformOperation       */
/************************************************/

#pragma mark -
#pragma mark GLRealToComplexTransformOperation
#pragma mark

@implementation GLRealToComplexTransformOperation

- (id) initWithOperand: (GLFunction *) variable destinationBasis: (NSArray *) toBasis
{   
	if (variable.dimensions.count != toBasis.count) {
		[NSException raise: @"DimensionsNotEqualException" format: @"The destination basis must be specified for all dimensions."];
	}
    
    // Basis array of the variable in its current state.
    NSMutableArray *initialBasis = [NSMutableArray array];
    for (GLDimension *aDim in variable.dimensions) {
		[initialBasis addObject: @(aDim.basisFunction)];
    }
	
    // Look for exponential transforms, *only*
    NSUInteger forwardRank = 0;
    NSUInteger inverseRank = 0;
    GLDimension *lastDim;
    for (GLDimension *aDim in variable.dimensions) {
        NSUInteger idx = [variable.dimensions indexOfObject: aDim];
        GLBasisFunction destinationBasis = [toBasis[idx] unsignedIntegerValue];
        
        if (aDim.basisFunction == kGLDeltaBasis && destinationBasis == kGLExponentialBasis) {
            forwardRank++;
            lastDim = aDim;
        } else if (aDim.basisFunction == kGLExponentialBasis && destinationBasis == kGLDeltaBasis) {
            inverseRank++;
            lastDim = aDim;
        } else if (aDim.basisFunction != destinationBasis) {
            [NSException raise: @"InvalidDestinationBasis" format: @"Invalid destination basis for a DFT."];
        }
    }
    
    if (forwardRank && inverseRank) {
        [NSException raise: @"InvalidDestinationBasis" format: @"You can't do a forward and inverse transformation simultaneously."];
    }
    
    if (!forwardRank && !inverseRank) {
        [NSException raise: @"InvalidDestinationBasis" format: @"No dimensions to transform!"];
    }
    
    if (forwardRank && !variable.isPurelyReal) {
		[NSException raise: @"VariableNotRealException" format: @"A forward real-to-complex transformation can only act on a real variable."];
	}
    
    if (inverseRank) {
        // Now check that the variable really is hermitian for the last dimension, and that the last dimension has only positive points.
    }
    
    NSMutableArray *finalBasis = [initialBasis mutableCopy];
    NSMutableArray *finalDimensions = [NSMutableArray array];
    for (GLDimension *aDim in variable.dimensions) {
        NSUInteger idx = [variable.dimensions indexOfObject: aDim];
        GLBasisFunction destinationBasis = [toBasis[idx] unsignedIntegerValue];
        
        if (aDim.basisFunction == destinationBasis) {
            [finalDimensions addObject: aDim];
        } else if (aDim.basisFunction == kGLDeltaBasis && destinationBasis == kGLExponentialBasis) {
            finalBasis[idx] = @(destinationBasis);
            // Note--the last dimension of the transform has only ~half the number of points.
            [finalDimensions addObject: [[GLDimension alloc] initAsDimension: aDim transformedToBasis: destinationBasis strictlyPositive: aDim == lastDim]];
        } else if (aDim.basisFunction == kGLExponentialBasis && destinationBasis == kGLDeltaBasis) {
            finalBasis[idx] = @(destinationBasis);
            // All the transformations, being the inverse, should ignore the strictly positive flag
            [finalDimensions addObject: [[GLDimension alloc] initAsDimension: aDim transformedToBasis: destinationBasis strictlyPositive: NO]];
        }
    }
    
    // Okay, we're good to go!
    GLFunction *resultVariable;
    if (forwardRank) {
        resultVariable = [GLFunction functionOfComplexTypeWithDimensions: finalDimensions forEquation: variable.equation];
    } else {
        resultVariable = [GLFunction functionOfRealTypeWithDimensions: finalDimensions forEquation: variable.equation];
    }
    
    resultVariable.name = variable.name;
	if (( self = [super initWithResult: @[resultVariable] operand: @[variable]] ))
	{
        self.fromBasis = initialBasis;
        self.toBasis = finalBasis;
        
        int rank = 0;
        int howmany_rank = 0;
        
        // First compute the strides between dimensions
        // It's important to use the operand's dimensions to get the strides correct.
        vGL_fftw_iodim *alldims = malloc(variable.dimensions.count*sizeof(vGL_fftw_iodim));
        NSArray *indim = variable.dimensions;
        NSArray *outdim = resultVariable.dimensions;
        for (NSInteger i=variable.dimensions.count-1; i >= 0; i--)
        {
//            GLDimension *inDim = [variable.dimensions objectAtIndex: i];
//            GLDimension *outDim = [resultVariable.dimensions objectAtIndex: i];
            
            if (i == variable.dimensions.count-1 ) {
                alldims[i].n = (int) [indim[i] nPoints];
                alldims[i].is = 1;
                alldims[i].os = 1;
            } else {
                alldims[i].n = (int) [indim[i] nPoints];
                alldims[i].is = alldims[i+1].is*((int)[indim[i+1] nPoints]);
                alldims[i].os = alldims[i+1].os*((int)[outdim[i+1] nPoints]);
            }
            
            if ( [initialBasis[i] isEqualTo: finalBasis[i]] ) {
                howmany_rank++;
            } else {
                rank++;
            }
        }
        
        // Now sort them into transformed and untransformed dimensions
        vGL_fftw_iodim *dims = malloc(rank*sizeof(vGL_fftw_iodim));
        vGL_fftw_iodim *howmany_dims = malloc(howmany_rank*sizeof(vGL_fftw_iodim));
        
        // From the FFTW manual:
        // The last dimension of dims is interpreted specially: that dimension of the real array has size dims[rank-1].n, but that dimension of the complex array has size dims[rank-1].n/2+1 (division rounded down). The strides, on the other hand, are taken to be exactly as specified.
        // This means that our above calculation needs to be corrected for the inverse transforms
        NSUInteger transformedDataPoints = 1;
        NSUInteger iDims = 0;
        NSUInteger iHowManyDims = 0;
        GLFloat scaleFactor = 1.0;
        for (NSUInteger i=0; i<variable.dimensions.count; i++)
        {
            GLDimension *aDim = [variable.dimensions objectAtIndex: i];
            
            if ( [initialBasis[i] isEqualTo: finalBasis[i]] )
            {   // This dimension does NOT get transformed
                howmany_dims[iHowManyDims].n =alldims[i].n;
                howmany_dims[iHowManyDims].is =alldims[i].is;
                howmany_dims[iHowManyDims].os =alldims[i].os;
                iHowManyDims++;
            }
            else
            {   // This dimension DOES get transformed
                if ( inverseRank && iDims == rank-1 ) {
                    dims[iDims].n = 2*(alldims[i].n-1); // Here's our funny fix
                } else {
                    dims[iDims].n =alldims[i].n;
                }
                dims[iDims].is =alldims[i].is;
                dims[iDims].os =alldims[i].os;
                
                if (forwardRank) {
                    scaleFactor *= 1.0 / ( (GLFloat) dims[iDims].n );
                }
                
                iDims++;
                
                transformedDataPoints *= aDim.nPoints;
            }
        }
        
        NSMutableData *tempf = [[GLMemoryPool sharedMemoryPool] dataWithLength: variable.nDataElements*sizeof(GLFloat)];
        NSMutableData *tempfbar = [[GLMemoryPool sharedMemoryPool] dataWithLength: resultVariable.nDataElements*sizeof(GLFloat)];
        
        [GLBasisTransformOperation readWisdom];
        
        if (forwardRank)
        {
            GLSplitComplex split = splitComplexFromData(tempfbar);
            vGL_fftw_plan plan = vGL_fftw_plan_guru_split_dft_r2c(rank, dims, howmany_rank, howmany_dims, tempf.mutableBytes, split.realp, split.imagp, FFTW_PATIENT);
            
            // Note that we leak the plan.
            NSUInteger dataPoints = resultVariable.nDataPoints;
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				NSMutableData *result = resultArray[0];
				NSMutableData *operand = operandArray[0];
				
                GLFloat *f = (void *) operand.bytes;
                GLSplitComplex fbar = splitComplexFromData(result);
                
                vGL_fftw_execute_split_dft_r2c( plan, f, fbar.realp, fbar.imagp );
                
                GLFloat inversesize = scaleFactor;
                vGL_vsmul( fbar.realp, 1, &inversesize, fbar.realp, 1, dataPoints );
                vGL_vsmul( fbar.imagp, 1, &inversesize, fbar.imagp, 1, dataPoints );
            };
        }
        else
        {
            GLSplitComplex split = splitComplexFromData(tempf);
            vGL_fftw_plan plan = vGL_fftw_plan_guru_split_dft_c2r(rank, dims, howmany_rank, howmany_dims, split.realp, split.imagp, tempfbar.mutableBytes, FFTW_PATIENT);

            // Note that we leak the plan.
//			self.blockOperation = ^(NSMutableData *result, NSData *operand) {
//                GLFloat *f = (void *) result.mutableBytes;
//                GLSplitComplex fbar = splitComplexFromData(operand);
//                
//                fftwf_execute_split_dft_c2r( plan, fbar.realp, fbar.imagp, f );
//            };
			
			// We are inefficiently copying to a buffer because the input is destroyed.
			NSUInteger numBytes = variable.nDataElements*sizeof(GLFloat);
			self.buffer = @[[[GLBuffer alloc] initWithLength: numBytes]];
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				NSMutableData *result = resultArray[0];
				NSMutableData *operand = operandArray[0];
				NSMutableData *buffer = bufferArray[0];
				
                GLFloat *f = (void *) result.mutableBytes;
				memcpy(buffer.mutableBytes, operand.bytes, numBytes);
                GLSplitComplex fbar = splitComplexFromData(buffer);
                
                vGL_fftw_execute_split_dft_c2r( plan, fbar.realp, fbar.imagp, f );
            };
        }
        
        [GLBasisTransformOperation saveWisdom];
        
        // Clean-up
        [[GLMemoryPool sharedMemoryPool] returnData: tempf];
        [[GLMemoryPool sharedMemoryPool] returnData: tempfbar];
        free(alldims);
        free(dims);
        free(howmany_dims);
    }
    
    return self;
}

@end

/************************************************/
/*		GLComplexToComplexTransformOperation    */
/************************************************/

#pragma mark -
#pragma mark GLComplexToComplexTransformOperation
#pragma mark

@implementation GLComplexToComplexTransformOperation

- (id) initWithOperand: (GLFunction *) variable destinationBasis: (NSArray *) toBasis
{   
	if (variable.dimensions.count != toBasis.count) {
		[NSException raise: @"DimensionsNotEqualException" format: @"The destination basis must be specified for all dimensions."];
	}
    
    if (!variable.isComplex) {
		[NSException raise: @"VariableNotComplexException" format: @"This DFT operation can only act on complex variables."];
	}
    
    // Basis array of the variable in its current state.
    NSMutableArray *initialBasis = [NSMutableArray array];
    for (GLDimension *aDim in variable.dimensions) {
		[initialBasis addObject: @(aDim.basisFunction)];
    }
	
    // Look for exponential transforms, *only*
    NSUInteger forwardRank = 0;
    NSUInteger inverseRank = 0;
    NSMutableArray *finalBasis = [initialBasis mutableCopy];
    NSMutableArray *finalDimensions = [NSMutableArray array];
    for (GLDimension *aDim in variable.dimensions) {
        NSUInteger idx = [variable.dimensions indexOfObject: aDim];
        GLBasisFunction destinationBasis = [toBasis[idx] unsignedIntegerValue];
        
        if (aDim.basisFunction == destinationBasis) {
            [finalDimensions addObject: aDim];
        } else if (aDim.basisFunction == kGLDeltaBasis && destinationBasis == kGLExponentialBasis) {
            forwardRank++;
            [finalDimensions addObject: [[GLDimension alloc] initAsDimension: aDim transformedToBasis: destinationBasis strictlyPositive: NO]];
        } else if (aDim.basisFunction == kGLExponentialBasis && destinationBasis == kGLDeltaBasis) {
            inverseRank++;
            [finalDimensions addObject: [[GLDimension alloc] initAsDimension: aDim transformedToBasis: destinationBasis strictlyPositive: NO]];
        } else {
            [NSException raise: @"InvalidDestinationBasis" format: @"Invalid destination basis for a DFT."];
        }
    }
    
    if (forwardRank && inverseRank) {
        [NSException raise: @"InvalidDestinationBasis" format: @"You can't do a forward and inverse transformation simultaneously."];
    }
    
    if (!forwardRank && !inverseRank) {
        [NSException raise: @"InvalidDestinationBasis" format: @"No dimensions to transform!"];
    }
    
    // Okay, we're good to go!
    GLFunction *resultVariable = [GLFunction functionOfComplexTypeWithDimensions: finalDimensions forEquation: variable.equation];
    resultVariable.name = variable.name;
	if (( self = [super initWithResult: @[resultVariable] operand: @[variable]] ))
	{
        self.fromBasis = initialBasis;
        self.toBasis = finalBasis;
        
        int rank = 0;
        int howmany_rank = 0;
        
        // First compute the strides between dimensions
        vGL_fftw_iodim *alldims = malloc(finalDimensions.count*sizeof(vGL_fftw_iodim));
        for (NSInteger i=finalDimensions.count-1; i >= 0; i--)
        {
            GLDimension *aDim = [finalDimensions objectAtIndex: i];
            
            if (i == finalDimensions.count-1 ) {
                alldims[i].n = (int) aDim.nPoints;
                alldims[i].is = 1;
                alldims[i].os = 1;
            } else {
                alldims[i].n = (int) aDim.nPoints;
                alldims[i].is = alldims[i+1].is*alldims[i+1].n;
                alldims[i].os = alldims[i+1].os*alldims[i+1].n;
            }
            
            if ( [initialBasis[i] isEqualTo: finalBasis[i]] ) {
                howmany_rank++;
            } else {
                rank++;
            }
        }
        
        // Now sort them into transformed and untransformed dimensions
        vGL_fftw_iodim *dims = malloc(rank*sizeof(vGL_fftw_iodim));
        vGL_fftw_iodim *howmany_dims = malloc(howmany_rank*sizeof(vGL_fftw_iodim));
        
        NSUInteger dataPoints = 1;
        NSUInteger transformedDataPoints = 1;
        NSUInteger iDims = 0;
        NSUInteger iHowManyDims = 0;
        for (NSUInteger i=0; i<finalDimensions.count; i++)
        {
            GLDimension *aDim = [finalDimensions objectAtIndex: i];
            dataPoints *= aDim.nPoints;
            
            if ( [initialBasis[i] isEqualTo: finalBasis[i]] )
            {   // This dimension does NOT get transformed
                howmany_dims[iHowManyDims].n =alldims[i].n;
                howmany_dims[iHowManyDims].is =alldims[i].is;
                howmany_dims[iHowManyDims].os =alldims[i].os;
                iHowManyDims++;
            }
            else
            {   // This dimension DOES get transformed
                dims[iDims].n =alldims[i].n;
                dims[iDims].is =alldims[i].is;
                dims[iDims].os =alldims[i].os;
                iDims++;
                
                transformedDataPoints *= aDim.nPoints;
            }
        }
        
        NSMutableData *tempf = [[GLMemoryPool sharedMemoryPool] dataWithLength: variable.nDataElements*sizeof(GLFloat)];
        NSMutableData *tempfbar = [[GLMemoryPool sharedMemoryPool] dataWithLength: resultVariable.nDataElements*sizeof(GLFloat)];
        GLSplitComplex split = splitComplexFromData(tempf);
        GLSplitComplex splitbar = splitComplexFromData(tempfbar);
        [GLBasisTransformOperation readWisdom];
        vGL_fftw_plan plan = vGL_fftw_plan_guru_split_dft(rank, dims, howmany_rank, howmany_dims, split.realp, split.imagp, splitbar.realp, splitbar.imagp, FFTW_PATIENT);
        [GLBasisTransformOperation saveWisdom];
        [[GLMemoryPool sharedMemoryPool] returnData: tempf];
        [[GLMemoryPool sharedMemoryPool] returnData: tempfbar];
        
        if (forwardRank)
        {
            // Compute the scale factor we need after the transformation.
            GLFloat scaleFactor = 1.0;
            for (NSUInteger i=0; i<initialBasis.count; i++)
            {
                GLBasisFunction fromBasis = [[initialBasis objectAtIndex: i] unsignedIntegerValue];
                GLBasisFunction toBasis = [[finalBasis objectAtIndex: i] unsignedIntegerValue];
                
                if (fromBasis == kGLDeltaBasis && toBasis == kGLExponentialBasis) {
                    scaleFactor *= 1.0 / ( (GLFloat) transformedDataPoints );
                } else if (fromBasis == kGLExponentialBasis && toBasis == kGLDeltaBasis) {
                    scaleFactor *= 1.0 / ( (GLFloat) transformedDataPoints );
                }
            }
            
            // Note that we leak the plan.
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				NSMutableData *result = resultArray[0];
				NSMutableData *operand = operandArray[0];
				
                GLSplitComplex f = splitComplexFromData(operand);
                GLSplitComplex fbar = splitComplexFromData(result);
                
                vGL_fftw_execute_split_dft( plan, f.realp, f.imagp, fbar.realp, fbar.imagp );
                
                GLFloat inversesize = scaleFactor;
                vGL_vsmul( fbar.realp, 1, &inversesize, fbar.realp, 1, dataPoints );
                vGL_vsmul( fbar.imagp, 1, &inversesize, fbar.imagp, 1, dataPoints );
            };
        }
        else
        {
            // Note that we leak the plan.
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				NSMutableData *result = resultArray[0];
				NSMutableData *operand = operandArray[0];
				
                GLSplitComplex fbar = splitComplexFromData(operand);
                GLSplitComplex f = splitComplexFromData(result);
                
                // Swap the real and imaginary parts for an inverse transform
                vGL_fftw_execute_split_dft( plan, fbar.imagp, fbar.realp, f.imagp, f.realp );
            };
        }
        // Clean-up
        
        free(alldims);
        free(dims);
        free(howmany_dims);
    }
    
    return self;
}

@end

@implementation GLMatrixFFTTransformOperation

+ (id) basisTransformationWithOperand: (GLFunction *) variable destinationBasis: (NSArray *) toBasis
{
    // If only one basis was specified, assume it was intended for all dimensions
    if ( toBasis.count == 1 && variable.dimensions.count > 1) {
        NSMutableArray *array = [NSMutableArray array];
        for (NSUInteger i=0; i < variable.dimensions.count; i++) {
            [array addObject: toBasis.lastObject];
        }
        toBasis = array;
    }
    
	if (variable.dimensions.count != toBasis.count) {
		[NSException raise: @"DimensionsNotEqualException" format: @"The destination basis must be specified for all dimensions."];
	}
    
    // Basis array of the variable in its current state.
    NSMutableArray *fromBasis = [NSMutableArray array];
    for (GLDimension *aDim in variable.dimensions) {
		[fromBasis addObject: @(aDim.basisFunction)];
    }
	
	// Check to see if a transformation is even necessary
	BOOL finishedTransformation = YES;
	for (NSInteger i=0; i < toBasis.count; i++) {
		finishedTransformation &= [toBasis[i] isEqualTo: fromBasis[i]];
	}
	if (finishedTransformation) {
		// No transformation necessary
		return nil;
	}
	
//    GLFourierTransformOperation *operation = [[GLFourierTransformOperation alloc] initWithOperand: @[variable]];
//    return operation;
	return nil;
}

@end

