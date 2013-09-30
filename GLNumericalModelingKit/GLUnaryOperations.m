//
//  GLUnaryOperations.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 1/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import "GLUnaryOperations.h"

#import "GLDimension.h"
#import "GLFourierTransformPool.h"
#import "GLFourierTransform.h"
#import "GLMemoryPool.h"

/************************************************/
/*		GLNegationOperation                     */
/************************************************/

@implementation GLNegationOperation

- (id) initWithOperand: (GLVariable *) variable;
{
	if (( self = [super initWithOperand: variable] ))
	{
		NSUInteger numElements = self.result.nDataElements;
		self.result.isPurelyReal = self.operand.isPurelyReal;
		self.result.isPurelyImaginary = self.operand.isPurelyImaginary;
		self.blockOperation = ^(NSMutableData *result, NSData *operand) {
			vGL_vneg( (void *) operand.bytes, 1, result.mutableBytes, 1, numElements );
		};
        self.graphvisDescription = @"negate";
	}
	
    return self;
}

- (BOOL) canOperateInPlace {
	return YES;
}

@end

/************************************************/
/*		GLAbsoluteValueOperation                */
/************************************************/

@implementation GLAbsoluteValueOperation

- (id) initWithOperand: (GLVariable *) variable {
	return [self initWithOperand: variable shouldUseComplexArithmetic: YES];
}

- (id) initWithOperand: (GLVariable *) variable shouldUseComplexArithmetic: (BOOL) useComplexArithmetic
{
	GLVariable *resultVar = [[variable class] variableOfRealTypeWithDimensions: variable.dimensions forEquation: variable.equation];
	
	if (( self = [super initWithResult: resultVar operand: variable] ))
	{
		self.result.isPurelyReal = YES;
        self.useComplexArithmetic = useComplexArithmetic;
		
		if (self.operand.isComplex && useComplexArithmetic == YES) {
			NSUInteger numPoints = self.result.nDataPoints;
			self.blockOperation = ^(NSMutableData *result, NSData *operand) {
				GLSplitComplex fromSplit = splitComplexFromData(operand);
				vGL_zvabs( &fromSplit, 1, result.mutableBytes, 1, numPoints );
			};
            self.graphvisDescription = @"abs (complex)";
		} else {
			NSUInteger numElements = self.result.nDataElements;
			self.blockOperation = ^(NSMutableData *result, NSData *operand) {
				vGL_vabs( (void *) operand.bytes, 1, result.mutableBytes, 1, numElements );
			};
            self.graphvisDescription = @"abs (real)";
		}
	}
	
    return self;
}

- (BOOL) canOperateInPlace {
	return YES;
}

- (BOOL) isEqualToOperation: (id) otherOperation {
    if ( ![super isEqualToOperation: otherOperation] )  {
        return NO;
    }
    
    GLAbsoluteValueOperation * op = otherOperation;
    if (self.useComplexArithmetic != op.useComplexArithmetic) {
        return NO;
    }
    
    return YES;
}

@end

/************************************************/
/*		GLExponentialOperation                  */
/************************************************/

@implementation GLExponentialOperation

- (id) initWithOperand: (GLVariable *) variable;
{
	// If the operand is purely real, we don't need a complex number
	GLDataFormat format = variable.isPurelyReal ? kGLRealDataFormat : kGLSplitComplexDataFormat;
	GLVariable *resultVar = [[variable class] variableOfType: format withDimensions: variable.dimensions forEquation: variable.equation];
	if (( self = [super initWithResult: resultVar operand: variable] ))
	{
		const int numElements = (int) self.result.nDataElements;
		const int numPoints = (int) self.result.nDataPoints;
		self.result.isPurelyReal = self.operand.isPurelyReal;
		
		NSMutableData *buffer = [[GLMemoryPool sharedMemoryPool] dataWithLength: numPoints*sizeof(GLFloat)];
		
		if (self.operand.isComplex) {
			if (self.operand.isPurelyReal) {
				self.blockOperation = ^(NSMutableData *result, NSData *operand) {
					vGL_vvexp( result.mutableBytes, operand.bytes, &numElements );
				};
                self.graphvisDescription = @"exp (complex, purely real)";
			} else if (self.operand.isPurelyImaginary) {
				self.blockOperation = ^(NSMutableData *result, NSData *operand) {
					GLSplitComplex fromSplit = splitComplexFromData(operand);
					GLSplitComplex toSplit = splitComplexFromData(result);
					
					vGL_vvsincos( toSplit.imagp, toSplit.realp, fromSplit.imagp, &numPoints);
				};
                self.graphvisDescription = @"exp (complex, purely imaginary)";
			} else {
				self.blockOperation = ^(NSMutableData *result, NSData *operand) {
					GLSplitComplex fromSplit = splitComplexFromData(operand);
					GLSplitComplex toSplit = splitComplexFromData(result);
					
					vGL_vvsincos( toSplit.imagp, toSplit.realp, fromSplit.imagp, &numPoints);
					vGL_vvexp( buffer.mutableBytes, fromSplit.realp, &numPoints );
					
					vGL_vmul( toSplit.realp, 1, buffer.mutableBytes, 1, toSplit.realp, 1, numPoints);
					vGL_vmul( toSplit.imagp, 1, buffer.mutableBytes, 1, toSplit.imagp, 1, numPoints);
				};
                self.graphvisDescription = @"exp (complex)";
			}
		} else {
			self.blockOperation = ^(NSMutableData *result, NSData *operand) {
				vGL_vvexp( result.mutableBytes, operand.bytes, &numElements );
			};
            self.graphvisDescription = @"exp (real)";
		}
	}
	
    return self;
}

- (BOOL) canOperateInPlace {
	return YES;
}

@end

/************************************************/
/*		GLSineOperation							*/
/************************************************/

@implementation GLSineOperation

- (id) initWithOperand: (GLVariable *) variable;
{
	if (( self = [super initWithOperand: variable] ))
	{
		const int numElements = (int) self.result.nDataElements;
		self.result.isPurelyReal = self.operand.isPurelyReal;
		self.result.isPurelyImaginary = self.operand.isPurelyImaginary;
		self.blockOperation = ^(NSMutableData *result, NSData *operand) {
			vGL_vvsin( result.mutableBytes, operand.bytes, &numElements );
		};
        self.graphvisDescription = @"sine";
	}
	
    return self;
}

- (BOOL) canOperateInPlace {
	return YES;
}

@end


/************************************************/
/*		GLCosineOperation							*/
/************************************************/

@implementation GLCosineOperation

- (id) initWithOperand: (GLVariable *) variable;
{
	if (( self = [super initWithOperand: variable] ))
	{
		const int numElements = (int) self.result.nDataElements;
		self.result.isPurelyReal = self.operand.isPurelyReal;
		self.result.isPurelyImaginary = self.operand.isPurelyImaginary;
		self.blockOperation = ^(NSMutableData *result, NSData *operand) {
			vGL_vvcos( result.mutableBytes, operand.bytes, &numElements );
		};
        self.graphvisDescription = @"cosine";
	}
	
    return self;
}

- (BOOL) canOperateInPlace {
	return YES;
}

@end

/************************************************/
/*		GLInverseTangentOperation				*/
/************************************************/

@implementation GLInverseTangentOperation

- (id) initWithOperand: (GLVariable *) variable;
{
	if (( self = [super initWithOperand: variable] ))
	{
		const int numElements = (int) self.result.nDataElements;
		self.result.isPurelyReal = self.operand.isPurelyReal;
		self.result.isPurelyImaginary = self.operand.isPurelyImaginary;
		self.blockOperation = ^(NSMutableData *result, NSData *operand) {
			vGL_vvatan( result.mutableBytes, operand.bytes, &numElements );
		};
        self.graphvisDescription = @"atan";
	}
	
    return self;
}

- (BOOL) canOperateInPlace {
	return YES;
}

@end

/************************************************/
/*		GLSquareRootOperation						*/
/************************************************/
// variable = sqrt( variable )
@implementation GLSquareRootOperation
- (id) initWithOperand: (GLVariable *) variable;
{
	if (( self = [super initWithOperand: variable] ))
	{
		const int numElements = (int) self.result.nDataElements;
		self.result.isPurelyReal = self.operand.isPurelyReal;
		self.result.isPurelyImaginary = self.operand.isPurelyImaginary;
		self.blockOperation = ^(NSMutableData *result, NSData *operand) {
			vGL_vvsqrt( result.mutableBytes, operand.bytes, &numElements );
		};
        self.graphvisDescription = @"sqrt";
	}
	
    return self;
}

- (BOOL) canOperateInPlace {
	return YES;
}
@end


/************************************************/
/*		GLFourierTransformOperation             */
/************************************************/

@implementation  GLFourierTransformOperation

- (id) initWithOperand: (GLVariable *) variable;
{
	NSMutableArray *transformedDimensions = [NSMutableArray array];
	for (GLDimension *aDim in variable.dimensions) {
		[transformedDimensions addObject: [aDim fourierTransformedDimension]];
	}
	
	// TODO - for the moment I'm assuming that anything NOT in the frequency domain is real.
	GLDataFormat format = !variable.isFrequencyDomain ? kGLSplitComplexDataFormat : kGLRealDataFormat;
	GLVariable *resultVariable = [[variable class] variableOfType: format withDimensions: transformedDimensions forEquation: variable.equation];
	
	if (( self = [super initWithResult: resultVariable operand: variable] ))
	{
		NSNumber *key = [[GLFourierTransformPool sharedFourierTransformPool] keyForDimensions: self.result.dimensions];
		if (self.result.isFrequencyDomain) {
			self.blockOperation = ^(NSMutableData *result, NSData *operand) {
				GLFourierTransform *transform = [[GLFourierTransformPool sharedFourierTransformPool] transformForKey: key];
				GLSplitComplex fbar = splitComplexFromData(result);
				GLFloat *f = (void *) operand.bytes;
				
				[transform transform: f forward: &fbar];
				
				[[GLFourierTransformPool sharedFourierTransformPool] returnTransform: transform];
			};
		} else  {
			self.blockOperation = ^(NSMutableData *result, NSData *operand) {
				GLFourierTransform *transform = [[GLFourierTransformPool sharedFourierTransformPool] transformForKey: key];
				GLSplitComplex fbar = splitComplexFromData(operand);
				GLFloat *f = (void *) result.bytes;
				
				[transform transform: &fbar inverse: f];
				
				[[GLFourierTransformPool sharedFourierTransformPool] returnTransform: transform];
			};
		}
	}
	
    return self;
}

// At the moment the fourier transforms are assumed to be out-of-place.
- (BOOL) canOperateInPlace {
	return NO;
}

@end

/************************************************/
/*		GLSwapComplexOperation                  */
/************************************************/

@implementation GLSwapComplexOperation

- (id) initWithOperand: (GLVariable *) variable
{
	GLVariable *resultVariable = [[variable class] variableOfComplexTypeWithDimensions: variable.dimensions forEquation: variable.equation];
	resultVariable.isPurelyImaginary = variable.isPurelyReal;
	resultVariable.isPurelyReal = variable.isPurelyImaginary;
	
	if (( self = [super initWithResult: resultVariable operand: variable] ))
	{		
		NSUInteger numBytes = self.result.nDataPoints*sizeof(GLFloat);
		if (self.operand.isComplex) {
			self.blockOperation = ^(NSMutableData *result, NSData *operand) {
				GLSplitComplex toSplit = splitComplexFromData(result);
				GLSplitComplex fromSplit = splitComplexFromData(operand);
				
				memcpy( toSplit.realp, fromSplit.imagp, numBytes );
				memcpy( toSplit.imagp, fromSplit.realp, numBytes );	
			};
            self.graphvisDescription = @"swap complex (complex)";
		} else {
			NSUInteger numPoints = self.result.nDataPoints;
			self.blockOperation = ^(NSMutableData *result, NSData *operand) {
				GLSplitComplex toSplit = splitComplexFromData(result);
				
				vGL_vclr( toSplit.realp, 1, numPoints);
				memcpy( toSplit.imagp, operand.bytes, numBytes );
			};
            self.graphvisDescription = @"swap complex (real)";
		}		
	}
	
    return self;
}

@end

/************************************************/
/*		GLCopyVariableOperation                 */
/************************************************/

@implementation GLCopyVariableOperation

- (id) initWithOperand: (GLVariable *) variable;
{
	if (( self = [super initWithOperand: variable] ))
	{
        self.result.name = variable.name;
		NSUInteger numBytes = self.result.nDataElements*sizeof(GLFloat);
		self.result.isPurelyReal = self.operand.isPurelyReal;
		self.result.isPurelyImaginary = self.operand.isPurelyImaginary;
		self.blockOperation = ^(NSMutableData *result, NSData *operand) {
			memcpy( result.mutableBytes, operand.bytes,  numBytes );
        
		};
        self.graphvisDescription = @"copy";
	}
	
    return self;
}

- (BOOL) canOperateInPlace {
	return NO;
}

@end

/************************************************/
/*		GLMaxOperation							*/
/************************************************/

@implementation GLMaxOperation

- (id) initWithOperand: (GLVariable *) variable
{
	
	GLVariable *resultVariable = [GLVariable variableOfRealTypeWithDimensions: [NSArray array] forEquation: variable.equation];
	
	if (( self = [super initWithResult: resultVariable operand: variable] ))
	{		
		NSUInteger nDataElements = self.operand.nDataElements;
		self.blockOperation = ^(NSMutableData *result, NSData *operand) {
			vGL_maxv( (GLFloat *) operand.bytes, 1, result.mutableBytes, nDataElements);
		};
        self.graphvisDescription = @"max";
	}
	
    return self;
}

@end

/************************************************/
/*		GLAverageOperation						*/
/************************************************/

@implementation GLAverageOperation

- (id) initWithOperand: (GLVariable *) variable dimensionIndex: (NSUInteger) index
{
	NSMutableArray *newDimensions = [variable.dimensions mutableCopy];
	[newDimensions removeObjectAtIndex: index];
	
	GLVariable *resultVariable = [GLVariable variableOfType: variable.dataFormat withDimensions: newDimensions forEquation: variable.equation];
	
	if (( self = [super initWithResult: resultVariable operand: variable] ))
	{
        self.dimIndex = index;
        
		// index = i*ny*nz + j*nz + k
		// For example: for a given (i,k), we want to sum over all j's. The stride between j elements, is nz and there are clearly ny of them.
		
		// The spacing between the elements that we want to sum over.
		// If we're collapsing the last index
		NSUInteger summingStride=1;
		for ( NSUInteger i=variable.dimensions.count-1; i>index; i--) {
			summingStride *= [[variable.dimensions objectAtIndex: i] nPoints];
		}
		// The number of elements that we want to sum over
		NSUInteger summingPoints = [[newDimensions objectAtIndex: index] nPoints];
		
        // The tricky part is covering all combinations of (i,k) --- of which there are nx*nz. In this example, if we iterate m=0, m<nx*nz,
		// then when m==nz, we actually want to skip to ny*nz, and then increment by ones again. (i/nz)*ny*nz + (i%nz)
        // Collapse nx:
        //  k = m%nz
        //  j = m/nz
        //  index = (m/nz)*nz + (m%nz)*1
        // Collapse ny:
        //  k = m%nz
        //  i = m/nz
        //  index = (m/nz)*ny*nz + (m%nz)*1
        // Collapse nz:
        //  j = m%ny
        //  i = m/ny
        //  index = (m/ny)*ny*nz + (m%ny)*nz
        //
        // Collapse nx (alt):
        //  index = (m/ny*nz)*0 + (m%ny*nz)*1
        // Collapse ny (alt):
        //  index = (m/nz)*ny*nz + (m%nz)*1
        // Collapse nz (alt):
        //  index = (m/1)*nz + (m%1)*0
        
        // This takes the more general form of: index = (m/divisor)*a + (m%divisor)*b
        // The divisor is the stride of the dimension being collapsed.
        // The coefficient a is the stride of the next slower dimension than the one being collapsed (or 0)
        // The coefficient b is 1 if a faster dimension remains, zero otherwise.
        
        // The number of points we'll return in the end.
		NSUInteger nDataPoints = resultVariable.nDataPoints;
        NSUInteger divisor = [GLDimension strideOfDimensionAtIndex: index inArray: variable.dimensions];
        NSUInteger a = index > 0 ? [GLDimension strideOfDimensionAtIndex: index-1 inArray: variable.dimensions] : 0;
        NSUInteger b = index < variable.dimensions.count-1 ? 1 : 0;
		
		if (resultVariable.dataFormat == kGLRealDataFormat)
		{
			self.blockOperation = ^(NSMutableData *result, NSData *operand) {
				
				dispatch_apply( nDataPoints, dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0), ^(size_t iteration) {
                    NSUInteger index = (iteration/divisor)*a + (iteration%divisor)*b;
					GLFloat *op = (GLFloat *) operand.bytes;
					GLFloat *res = (GLFloat *) result.mutableBytes;
					vGL_meanv(&op[index], summingStride, &res[iteration], summingPoints);
				});
				
			};
            self.graphvisDescription = [NSString stringWithFormat: @"average dimension %lu (real formatting)", index] ;
		}
        else if (resultVariable.dataFormat == kGLSplitComplexDataFormat)
		{
			self.blockOperation = ^(NSMutableData *result, NSData *operand) {
				
				dispatch_apply( nDataPoints, dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0), ^(size_t iteration) {
                    NSUInteger index = (iteration/divisor)*a + (iteration%divisor)*b;
                    GLSplitComplex toSplit = splitComplexFromData(result);
                    GLSplitComplex fromSplit = splitComplexFromData(operand);
					vGL_meanv(fromSplit.realp + index*sizeof(GLFloat), summingStride, toSplit.realp + index*sizeof(GLFloat), summingPoints);
                    vGL_meanv(fromSplit.imagp + index*sizeof(GLFloat), summingStride, toSplit.imagp + index*sizeof(GLFloat), summingPoints);

				});
				
			};
            self.graphvisDescription = [NSString stringWithFormat: @"average dimension %lu (split complex formatting)", index] ;
		}

	}
	
    return self;
}

- (BOOL) isEqualToOperation: (id) otherOperation {
    if ( ![super isEqualToOperation: otherOperation] )  {
        return NO;
    }
    
    GLAverageOperation * op = otherOperation;
    if (self.dimIndex != op.dimIndex) {
        return NO;
    }
    
    return YES;
}

@end

