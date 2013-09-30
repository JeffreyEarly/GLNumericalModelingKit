//
//  GLVectorVectorOperations.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 1/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import "GLVectorVectorOperations.h"
#import "GLMatrixDescription.h"

/************************************************/
/*		GLAdditionOperation						*/
/************************************************/

@implementation GLAdditionOperation

- (id) initWithFirstOperand: (GLVariable *) fOperand secondOperand: (GLVariable *) sOperand
{
    // We order the operands so that scalars are always in the first position.
    // We can do this in this case because order doesn't matter for addition.
    GLDataFormat format = (fOperand.isPurelyReal && sOperand.isPurelyReal) ? kGLRealDataFormat : kGLSplitComplexDataFormat;
    GLVariable *op1 = (sOperand.dimensions.count < fOperand.dimensions.count) ? sOperand : fOperand;
    GLVariable *op2 = (sOperand.dimensions.count < fOperand.dimensions.count) ? fOperand : sOperand;
    GLVariable *result = [[op2 class] variableOfType:format withDimensions: op2.dimensions forEquation: fOperand.equation];
    
	if (( self = [super initWithResult: result firstOperand: op1 secondOperand: op2] )) {
		self.result.isPurelyReal = self.firstOperand.isPurelyReal && self.secondOperand.isPurelyReal;
		self.result.isPurelyImaginary = self.firstOperand.isPurelyImaginary && self.secondOperand.isPurelyImaginary;
		
		if (!self.firstOperand.dimensions.count && !self.secondOperand.dimensions.count)
		{ // scalar-scalar addition
			if (!self.firstOperand.isComplex && !self.secondOperand.isComplex) {
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					GLFloat *a = (GLFloat *) fOperand.bytes;
					GLFloat *b = (GLFloat *) sOperand.bytes;
					GLFloat *c = (GLFloat *) result.mutableBytes;
					*c = (*a) + (*b);
				};
                self.graphvisDescription = [NSString stringWithFormat: @"add (real scalar, real scalar)"];
			} else {
				[NSException raise: @"MethodNotImplemented" format: @"This case has not yet been implemented."];
			}
		}
        else if (!self.firstOperand.dimensions.count && self.secondOperand.dimensions.count)
        {   // scalar-vector addition
            if (!self.firstOperand.isComplex && !self.secondOperand.isComplex) {
                NSUInteger numPoints = self.result.nDataPoints;
                self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
                    vGL_vsadd( (void *) sOperand.bytes, 1, (void *) fOperand.bytes, result.mutableBytes, 1, numPoints);
                };
                self.graphvisDescription = [NSString stringWithFormat: @"add (scalar, vector)"];
            } else {
				[NSException raise: @"MethodNotImplemented" format: @"This case has not yet been implemented."];
			}
        }
        else if (self.firstOperand.dimensions.count == self.secondOperand.dimensions.count)
        {
			if (!fOperand.isComplex && sOperand.isComplex) {
				NSUInteger nDataPoints = self.result.nDataPoints;
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					GLSplitComplex leftComplex = splitComplexFromData( fOperand );
					GLSplitComplex rightComplex = splitComplexFromData( sOperand );
					GLSplitComplex destComplex = splitComplexFromData( result );
					
					vGL_vadd( leftComplex.realp, 1, rightComplex.realp, 1, destComplex.realp, 1, nDataPoints);      
					vGL_mmov( rightComplex.imagp, destComplex.imagp, nDataPoints, 1, nDataPoints, nDataPoints);
				};
                self.graphvisDescription = [NSString stringWithFormat: @"add (real,complex)"];
			} else if (fOperand.isComplex && !sOperand.isComplex) {
				NSUInteger nDataPoints = self.result.nDataPoints;
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					GLSplitComplex leftComplex = splitComplexFromData( fOperand );
					GLSplitComplex rightComplex = splitComplexFromData( sOperand );
					GLSplitComplex destComplex = splitComplexFromData( result );
					
					vGL_vadd( leftComplex.realp, 1, rightComplex.realp, 1, destComplex.realp, 1, nDataPoints);      
					vGL_mmov( leftComplex.imagp, destComplex.imagp, nDataPoints, 1, nDataPoints, nDataPoints);
				};
                self.graphvisDescription = [NSString stringWithFormat: @"add (complex,real)"];
			} else {
				NSUInteger nDataElements = self.result.nDataElements;
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					vGL_vadd( fOperand.bytes, 1, sOperand.bytes, 1, result.mutableBytes, 1, nDataElements);
				};
                self.graphvisDescription = [NSString stringWithFormat: @"add"];
			}
		} else {
            [NSException raise: @"MethodNotImplemented" format: @"This case has not yet been implemented."];
        }
    }
    return self;
}

- (BOOL) canOperateInPlace {
	return YES;
}

@end

/************************************************/
/*		GLSubtractionOperation					*/
/************************************************/

// Output = leftVariable - rightVariable

@implementation GLSubtractionOperation

- (id) initWithFirstOperand: (GLVariable *) fOperand secondOperand: (GLVariable *) sOperand
{
    // We order the operands so that scalars are always in the first position.
    GLDataFormat format = (fOperand.isPurelyReal && sOperand.isPurelyReal) ? kGLRealDataFormat : kGLSplitComplexDataFormat;
    NSArray *dims = fOperand.dimensions.count ? fOperand.dimensions : sOperand.dimensions;
    Class aClass = fOperand.dimensions.count ? [fOperand class] : [sOperand class];
    GLVariable *result = [[aClass class] variableOfType:format withDimensions: dims forEquation: fOperand.equation];
    
	if (( self = [super initWithResult: result firstOperand:fOperand secondOperand:sOperand] )) {
		self.result.isPurelyReal = self.firstOperand.isPurelyReal && self.secondOperand.isPurelyReal;
		self.result.isPurelyImaginary = self.firstOperand.isPurelyImaginary && self.secondOperand.isPurelyImaginary;
		
		if (!self.firstOperand.dimensions.count && !self.secondOperand.dimensions.count)
		{ // scalar-scalar subtraction
			if (!self.firstOperand.isComplex && !self.secondOperand.isComplex) {
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					GLFloat *a = (GLFloat *) fOperand.bytes;
					GLFloat *b = (GLFloat *) sOperand.bytes;
					GLFloat *c = (GLFloat *) result.mutableBytes;
					*c = (*a) - (*b);
				};
                self.graphvisDescription = [NSString stringWithFormat: @"subtract (real scalar, real scalar)"];
			} else {
				[NSException raise: @"MethodNotImplemented" format: @"This case has not yet been implemented."];
			}
		}
        else if ((!self.firstOperand.dimensions.count && self.secondOperand.dimensions.count) || (self.firstOperand.dimensions.count && !self.secondOperand.dimensions.count))
        {   // scalar-vector subtraction
            if (!self.firstOperand.isComplex && !self.secondOperand.isComplex) {
                NSUInteger numElements = self.result.nDataElements;
                if (!self.firstOperand.dimensions.count) {
                    // result = scalar - vector
                    self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
                        vGL_vneg( (void *) sOperand.bytes, 1, result.mutableBytes, 1, numElements );
                        vGL_vsadd( (void *) result.bytes, 1, (void *) fOperand.bytes, result.mutableBytes, 1, numElements);
                    };
                } else {
                    // result = vector - scalar
                    self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
                        GLFloat *a = (GLFloat *)sOperand.bytes;
                        GLFloat b = - (*a);
                        vGL_vsadd( (void *) fOperand.bytes, 1, &b, result.mutableBytes, 1, numElements);
                    };
                }
                self.graphvisDescription = [NSString stringWithFormat: @"subtract (scalar, vector)"];
            } else {
				[NSException raise: @"MethodNotImplemented" format: @"This case has not yet been implemented."];
			}
        }
        else if (self.firstOperand.dimensions.count == self.secondOperand.dimensions.count)
        {
			if (!fOperand.isComplex && sOperand.isComplex) {
				NSUInteger nDataPoints = self.result.nDataPoints;
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					GLSplitComplex leftComplex = splitComplexFromData( fOperand );
					GLSplitComplex rightComplex = splitComplexFromData( sOperand );
					GLSplitComplex destComplex = splitComplexFromData( result );
					
					vGL_vsub( rightComplex.realp, 1, leftComplex.realp, 1, destComplex.realp, 1, nDataPoints);
					vGL_vneg( rightComplex.imagp, 1, destComplex.imagp, 1, nDataPoints);
				};
                self.graphvisDescription = [NSString stringWithFormat: @"subtract (real,complex)"];
			} else if (fOperand.isComplex && !sOperand.isComplex) {
				NSUInteger nDataPoints = self.result.nDataPoints;
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					GLSplitComplex leftComplex = splitComplexFromData( fOperand );
					GLSplitComplex rightComplex = splitComplexFromData( sOperand );
					GLSplitComplex destComplex = splitComplexFromData( result );
					
					vGL_vsub( rightComplex.realp, 1, leftComplex.realp, 1, destComplex.realp, 1, nDataPoints);      
					vGL_mmov( leftComplex.imagp, destComplex.imagp, nDataPoints, 1, nDataPoints, nDataPoints);
				};
                self.graphvisDescription = [NSString stringWithFormat: @"subtract (complex,real)"];
			} else {
				NSUInteger nDataElements = self.result.nDataElements;
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					// Note that vsub does: C = B - A
					vGL_vsub( sOperand.bytes, 1, fOperand.bytes, 1, result.mutableBytes, 1, nDataElements);
				};
                self.graphvisDescription = [NSString stringWithFormat: @"subtract"];
			}
		} else {
            [NSException raise: @"MethodNotImplemented" format: @"This case has not yet been implemented."];
        }
    }
    return self;
}

- (BOOL) canOperateInPlace {
	return YES;
}

@end

/************************************************/
/*		GLAbsoluteLargestOperation				*/
/************************************************/
// variable = max( abs(leftVariable), abs(rightVariable ) element-wise

@implementation GLAbsoluteLargestOperation

- (id) initWithFirstOperand: (GLVariable *) fOperand secondOperand: (GLVariable *) sOperand
{
	if (( self = [super initWithFirstOperand:fOperand secondOperand:sOperand] )) {
		self.result.isPurelyReal = self.firstOperand.isPurelyReal && self.secondOperand.isPurelyReal;
		self.result.isPurelyImaginary = self.firstOperand.isPurelyImaginary && self.secondOperand.isPurelyImaginary;
		
		if (!self.firstOperand.dimensions.count && !self.secondOperand.dimensions.count)
		{ // scalar-scalar multiplication
			[NSException raise: @"MethodNotImplemented" format: @"This case has not yet been implemented."];
		}
		else if (!self.firstOperand.dimensions.count && self.secondOperand.dimensions.count)
		{ // scalar-vector multiplication
			if (!self.firstOperand.isComplex) {
				NSUInteger nDataElements = self.secondOperand.nDataElements;
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					vGL_vthr( (GLFloat *) sOperand.bytes, 1, (GLFloat *) fOperand.bytes, result.mutableBytes, 1, nDataElements);
				};
                self.graphvisDescription = [NSString stringWithFormat: @"element-wise max (real scalar, other)"];
			}
		}
		else if (self.firstOperand.dimensions.count && !self.secondOperand.dimensions.count)
		{ // vector-scalar multiplication
			if (!self.secondOperand.isComplex) {
				NSUInteger nDataElements = self.firstOperand.nDataElements;
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					vGL_vthr( (GLFloat *) fOperand.bytes, 1, (GLFloat *) sOperand.bytes, result.mutableBytes, 1, nDataElements);
				};
                self.graphvisDescription = [NSString stringWithFormat: @"element-wise max (other, real scalar)"];
			}
		}
		else
		{
			NSUInteger nDataElements = self.result.nDataElements;
			self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
				vGL_vmaxmg( (GLFloat *) sOperand.bytes, 1, (GLFloat *) fOperand.bytes, 1, result.mutableBytes, 1, nDataElements);
			};
            self.graphvisDescription = [NSString stringWithFormat: @"element-wise max"];
		}
    }
    return self;
}

- (BOOL) canOperateInPlace {
	return YES;
}

@end

/************************************************/
/*		GLMultiplicationOperation				*/
/************************************************/

@interface GLMultiplicationOperation ()
@property(readwrite) BOOL canOperateInPlace;
@end

@implementation GLMultiplicationOperation

@synthesize canOperateInPlace;

- (id) initWithFirstOperand: (GLVariable *) fOperand secondOperand: (GLVariable *) sOperand
{
	/********************************************************************/
	/*		vector - lower dimensional vector multiplication			*/
	/********************************************************************/
	if (fOperand.dimensions.count && sOperand.dimensions.count && (fOperand.dimensions.count < sOperand.dimensions.count || sOperand.dimensions.count < fOperand.dimensions.count))
	{
		GLVariable *lowerDimVariable;
		GLVariable *higherDimVariable;
		BOOL flipOperands = NO;
		if (fOperand.dimensions.count < sOperand.dimensions.count) {
			lowerDimVariable = fOperand;
			higherDimVariable = sOperand;
		} else {
			lowerDimVariable = sOperand;
			higherDimVariable = fOperand;
			flipOperands = YES;
		}
		
		// In this scenario the first-operand has fewer dimensions than the second operand, but those dimensions are in the same order.
		// So, this would look like h(x,y) = f(x)*g(x,y)
		// Note that we are not allowing h(x,y,z) = f(x,z)*g(x,y,z) because the indexing is trickier
		NSInteger lastIndex = NSNotFound;
		for (GLDimension *dim in lowerDimVariable.dimensions) {
			NSUInteger index = [higherDimVariable.dimensions indexOfObject: dim];
			if (index == NSNotFound) {
				[NSException raise: @"Dimensional mismatch" format: @"The lower dimensional variable must have a subset of dimensions from the higher dimensional variable. It does not!"];
			} else if ( lastIndex != NSNotFound && lastIndex+1 != index) {
				[NSException raise: @"Dimensional mismatch" format: @"The lower dimensional variable must have a subset of dimensions *in the same order (with no gaps)* as the higher dimensional variable. It does not!"];
			}
			lastIndex = index;
		}
		
		// Now we know enough to build the result.
		BOOL isPurelyReal = (fOperand.isPurelyReal && sOperand.isPurelyReal) || (fOperand.isPurelyImaginary && sOperand.isPurelyImaginary);
		GLVariable *result = [[higherDimVariable class] variableOfType: isPurelyReal ? kGLRealDataFormat : kGLSplitComplexDataFormat withDimensions: higherDimVariable.dimensions forEquation: higherDimVariable.equation];
		result.isPurelyReal = isPurelyReal;
		result.isPurelyImaginary= (fOperand.isPurelyReal && sOperand.isPurelyImaginary) || (fOperand.isPurelyImaginary && sOperand.isPurelyReal);
		
		if (( self = [super init] )) {
			self.firstOperand = lowerDimVariable;
			self.secondOperand = higherDimVariable;
			self.result = result;
			[self performSelector:@selector(setupDependencies)];
			
			// Example: h(x,y,z) = f(y)*g(x,y,z)
			// where g_ijk is indexed with i*ny*nz + j*nz + k
			// and therefore we want to loop over i*ny*nz+k while striding by nz
			//
			// multiplicationStride = nz
			// multiplicationLength = ny
			//
			// innerLoopStride = 1
			// outerLoopStride = ny*nz
			// innerLoopSize = nz
			// innerLoopLength = nx*nz
			
			GLMatrixDescription *matrixDescription = higherDimVariable.matrixDescription;
			NSUInteger lastMultiplicationIndex = [higherDimVariable.dimensions indexOfObject: lowerDimVariable.dimensions.lastObject];
			
			NSUInteger multiplicationStride = matrixDescription.strides[lastMultiplicationIndex].stride;
			NSUInteger multiplicationLength = lowerDimVariable.nDataPoints;
			
			NSMutableArray *missingDimensions = [NSMutableArray arrayWithArray: higherDimVariable.dimensions];
			[missingDimensions removeObjectsInArray: lowerDimVariable.dimensions];
			NSUInteger lastMissingDimIndex = [higherDimVariable.dimensions indexOfObject: missingDimensions.lastObject];
			
			NSUInteger innerLoopStride = matrixDescription.strides[lastMissingDimIndex].stride;
			NSUInteger innerLoopSize = matrixDescription.strides[lastMissingDimIndex].nPoints;
			
			NSUInteger outerLoopStride = 0;
			NSUInteger outerLoopSize = 1;
			
			if ( missingDimensions.count > 1) {
				NSUInteger firstMissingDimIndex = [higherDimVariable.dimensions indexOfObject: missingDimensions[0]];
				outerLoopStride = matrixDescription.strides[firstMissingDimIndex].stride;
				outerLoopSize = matrixDescription.strides[firstMissingDimIndex].nPoints;
			}
			
			if ( !self.firstOperand.isComplex && !self.secondOperand.isComplex )
			{
				dispatch_queue_t globalQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
				NSUInteger totalLoops = outerLoopSize*innerLoopSize;
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					dispatch_apply(totalLoops, globalQueue, ^(size_t iteration) {
						
						const GLFloat *lowerF = fOperand.bytes;
						const GLFloat *higherF = sOperand.bytes;
						GLFloat *higherOut = result.mutableBytes;
						
						NSUInteger index = (iteration/innerLoopSize)*outerLoopStride + (iteration%innerLoopSize)*innerLoopStride;
						
						vGL_vmul( lowerF, 1, &(higherF[index]), multiplicationStride, &(higherOut[index]), multiplicationStride, multiplicationLength);
					});
				};
				self.canOperateInPlace = YES;
                self.graphvisDescription = [NSString stringWithFormat: @"multiplication (real %lu dim, real %lu dim)", (unsigned long)self.firstOperand.dimensions.count, (unsigned long)self.secondOperand.dimensions.count];
			} else if ( self.firstOperand.isComplex && !self.secondOperand.isComplex )
			{
				dispatch_queue_t globalQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
				NSUInteger totalLoops = outerLoopSize*innerLoopSize;
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					dispatch_apply(totalLoops, globalQueue, ^(size_t iteration) {
						const GLSplitComplex lowerF = splitComplexFromData(fOperand);
						const GLFloat *higherF = sOperand.bytes;
						const GLSplitComplex higherOut = splitComplexFromData(result);
						
						NSUInteger index = (iteration/innerLoopSize)*outerLoopStride + (iteration%innerLoopSize)*innerLoopStride;
						// (a + i b)*(x) = a x + i b x
						vGL_vmul( lowerF.realp, 1, &(higherF[index]), multiplicationStride, &(higherOut.realp[index]), multiplicationStride, multiplicationLength);
						vGL_vmul( lowerF.imagp, 1, &(higherF[index]), multiplicationStride, &(higherOut.imagp[index]), multiplicationStride, multiplicationLength);
					});
				};
				self.canOperateInPlace = YES;
                self.graphvisDescription = [NSString stringWithFormat: @"multiplication (complex %lu dim, real %lu dim)", (unsigned long)self.firstOperand.dimensions.count, (unsigned long)self.secondOperand.dimensions.count];
			} else {
				[NSException raise: @"MethodNotImplemented" format: @"This case has not yet been implemented."];
			}

		}
		return self;
	}
	
	if (( self = [super initWithFirstOperand:fOperand secondOperand:sOperand] )) {
		self.result.isPurelyReal = (self.firstOperand.isPurelyReal && self.secondOperand.isPurelyReal) || (self.firstOperand.isPurelyImaginary && self.secondOperand.isPurelyImaginary);
		self.result.isPurelyImaginary= (self.firstOperand.isPurelyReal && self.secondOperand.isPurelyImaginary) || (self.firstOperand.isPurelyImaginary && self.secondOperand.isPurelyReal);
		
        if (sOperand.name && fOperand.name) {
            self.result.name = [NSString stringWithFormat: @"%@_%@", sOperand.name, fOperand.name];
        }
		
		/********************************************************************/
		/*		scalar - scalar multiplication								*/
		/********************************************************************/
		if (!self.firstOperand.dimensions.count && !self.secondOperand.dimensions.count)
		{ // scalar-scalar multiplication
			if (!self.firstOperand.isComplex && !self.secondOperand.isComplex) {
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					GLFloat *a = (GLFloat *) fOperand.bytes;
					GLFloat *b = (GLFloat *) sOperand.bytes;
					GLFloat *c = (GLFloat *) result.mutableBytes;
					*c = (*a) * (*b);
				};
                self.graphvisDescription = @"multiplication (real scalar, real scalar)";
			} else {
				[NSException raise: @"MethodNotImplemented" format: @"This case has not yet been implemented."];
			}
		}
		/********************************************************************/
		/*		scalar - vector multiplication								*/
		/********************************************************************/
		else if (!self.firstOperand.dimensions.count && self.secondOperand.dimensions.count)
		{ // scalar-vector multiplication
			if (!self.firstOperand.isComplex) {
				NSUInteger nDataElements = self.secondOperand.nDataElements;
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					vGL_vsmul( sOperand.bytes, 1, fOperand.bytes, result.mutableBytes, 1, nDataElements);
				};
                self.graphvisDescription = @"multiplication (real scalar, vector)";
			} else {
				[NSException raise: @"MethodNotImplemented" format: @"This case has not yet been implemented."];
			}
		}
		else if (self.firstOperand.dimensions.count && !self.secondOperand.dimensions.count)
		{ // vector-scalar multiplication
			if (!self.secondOperand.isComplex) {
				NSUInteger nDataElements = self.firstOperand.nDataElements;
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					vGL_vsmul( fOperand.bytes, 1, sOperand.bytes, result.mutableBytes, 1, nDataElements);
				};
                self.graphvisDescription = @"multiplication (vector, real scalar)";
			} else {
				[NSException raise: @"MethodNotImplemented" format: @"This case has not yet been implemented."];
			}
		}
		/********************************************************************/
		/*		vector - vector multiplication								*/
		/********************************************************************/
		else
		{
			if ( !self.firstOperand.isComplex && !self.secondOperand.isComplex )
			{
				NSUInteger nDataElements = self.result.nDataElements;
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					vGL_vmul( fOperand.bytes, 1, sOperand.bytes, 1, result.mutableBytes, 1, nDataElements);
				};
				self.canOperateInPlace = YES;
                self.graphvisDescription = @"multiplication (real, real)";
			}
			else if ( self.firstOperand.isComplex && !self.secondOperand.isComplex )
			{
				NSUInteger nDataPoints = self.result.nDataPoints;
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					GLSplitComplex leftComplex = splitComplexFromData( fOperand );
					GLSplitComplex destComplex = splitComplexFromData( result );
					// (a + i b)*(x) = a x + i b x
					vGL_vmul( leftComplex.realp, 1, sOperand.bytes, 1, destComplex.realp, 1, nDataPoints);
					vGL_vmul( leftComplex.imagp, 1, sOperand.bytes, 1, destComplex.imagp, 1, nDataPoints);
				};
				self.canOperateInPlace = YES;
                self.graphvisDescription = @"multiplication (complex, real)";
			}
			else if ( !self.firstOperand.isComplex && self.secondOperand.isComplex )
			{
				NSUInteger nDataPoints = self.result.nDataPoints;
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					GLSplitComplex rightComplex = splitComplexFromData( sOperand );
					GLSplitComplex destComplex = splitComplexFromData( result );
					// a *(x + i y) = a x + i a y
					vGL_vmul( rightComplex.realp, 1, fOperand.bytes, 1, destComplex.realp, 1, nDataPoints);
					vGL_vmul( rightComplex.imagp, 1, fOperand.bytes, 1, destComplex.imagp, 1, nDataPoints);
				};
				self.canOperateInPlace = YES;
                self.graphvisDescription = @"multiplication (real, complex)";
			}
			else {
				NSUInteger nDataPoints = self.result.nDataPoints;
				
				if ((!self.firstOperand.isPurelyReal && !self.firstOperand.isPurelyImaginary) && self.secondOperand.isPurelyReal) {
					self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
						
						GLSplitComplex leftComplex = splitComplexFromData( fOperand );
						GLSplitComplex rightComplex = splitComplexFromData( sOperand );
						GLSplitComplex destComplex = splitComplexFromData( result );
						
						// (a + i b)*(x) = a x + i b x
						vGL_vmul( leftComplex.realp, 1, rightComplex.realp, 1, destComplex.realp, 1, nDataPoints);
						vGL_vmul( leftComplex.imagp, 1, rightComplex.realp, 1, destComplex.imagp, 1, nDataPoints);					
					};
					self.canOperateInPlace = NO;
                    self.graphvisDescription = @"multiplication (complex, complex purely real)";
				} else if ((!self.firstOperand.isPurelyReal && !self.firstOperand.isPurelyImaginary) && self.secondOperand.isPurelyImaginary) {
					self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
						
						GLSplitComplex leftComplex = splitComplexFromData( fOperand );
						GLSplitComplex rightComplex = splitComplexFromData( sOperand );
						GLSplitComplex destComplex = splitComplexFromData( result );
						
						// (a + i b)*(i y) = - b y + i a y
						vGL_vmul( leftComplex.imagp, 1, rightComplex.imagp, 1, destComplex.realp, 1, nDataPoints);
						vGL_vneg( destComplex.realp, 1, destComplex.realp, 1, nDataPoints );
						vGL_vmul( leftComplex.realp, 1, rightComplex.imagp, 1, destComplex.imagp, 1, nDataPoints);					
					};
					self.canOperateInPlace = NO;
                    self.graphvisDescription = @"multiplication (complex, complex purely imaginary)";
				} else if (self.firstOperand.isPurelyReal && (!self.secondOperand.isPurelyReal && !self.secondOperand.isPurelyImaginary)) {
					self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
						
						GLSplitComplex leftComplex = splitComplexFromData( fOperand );
						GLSplitComplex rightComplex = splitComplexFromData( sOperand );
						GLSplitComplex destComplex = splitComplexFromData( result );
						
						// a *(x + i y) = a x + i a y
						vGL_vmul( leftComplex.realp, 1, rightComplex.realp, 1, destComplex.realp, 1, nDataPoints);
						vGL_vmul( leftComplex.realp, 1, rightComplex.imagp, 1, destComplex.imagp, 1, nDataPoints);					
					};
					self.canOperateInPlace = NO;
                    self.graphvisDescription = @"multiplication (complex purely real, complex)";
				} else if (self.firstOperand.isPurelyImaginary && (!self.secondOperand.isPurelyReal && !self.secondOperand.isPurelyImaginary)) {
					self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
						
						GLSplitComplex leftComplex = splitComplexFromData( fOperand );
						GLSplitComplex rightComplex = splitComplexFromData( sOperand );
						GLSplitComplex destComplex = splitComplexFromData( result );
						
						// i b *(x + i y) = - b y + i b x
						vGL_vmul( leftComplex.imagp, 1, rightComplex.imagp, 1, destComplex.realp, 1, nDataPoints);
						vGL_vneg( destComplex.realp, 1, destComplex.realp, 1, nDataPoints );
						vGL_vmul( leftComplex.imagp, 1, rightComplex.realp, 1, destComplex.imagp, 1, nDataPoints);					
					};
					self.canOperateInPlace = NO;
                    self.graphvisDescription = @"multiplication (complex purely imaginary, complex)";
				} else {
					self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
						
						GLSplitComplex leftComplex = splitComplexFromData( fOperand );
						GLSplitComplex rightComplex = splitComplexFromData( sOperand );
						GLSplitComplex destComplex = splitComplexFromData( result );
						
						// (a + i b)*(x + i y) = (a x - b y) + i (b x + ay)
						vGL_vmmsb( leftComplex.realp, 1, rightComplex.realp, 1, leftComplex.imagp, 1, rightComplex.imagp, 1, destComplex.realp, 1, nDataPoints);
						vGL_vmma( leftComplex.imagp, 1, rightComplex.realp, 1, leftComplex.realp, 1, rightComplex.imagp, 1, destComplex.imagp, 1, nDataPoints);		
					};
					self.canOperateInPlace = NO;
                    self.graphvisDescription = @"multiplication (complex, complex)";
				}
	//#warning overriding this for testing purposes.
	//			self.canOperateInPlace = NO;
			}
		}
    }
    return self;
}

@end

/************************************************/
/*		GLDivisionOperation						*/
/************************************************/
// variable = leftVariable / rightVariable
@implementation GLDivisionOperation

- (id) initWithFirstOperand: (GLVariable *) fOperand secondOperand: (GLVariable *) sOperand {
	return [self initWithFirstOperand: fOperand secondOperand: sOperand shouldUseComplexArithmetic: YES];
}

- (id) initWithFirstOperand: (GLVariable *) fOperand secondOperand: (GLVariable *) sOperand shouldUseComplexArithmetic: (BOOL) useComplexArithmetic
{
	BOOL complexResult = fOperand.isComplex || sOperand.isComplex;
	GLDataFormat format = complexResult ? kGLSplitComplexDataFormat : kGLRealDataFormat;
	GLVariable *resultVariable = [[fOperand class] variableOfType: format withDimensions: fOperand.dimensions.count ? fOperand.dimensions : sOperand.dimensions forEquation: fOperand.equation];
	
	if (( self = [super initWithResult: resultVariable firstOperand:fOperand secondOperand: sOperand] ))
	{
		self.result.isPurelyReal = self.firstOperand.isPurelyReal && self.secondOperand.isPurelyReal;
        self.useComplexArithmetic = useComplexArithmetic;
		
		if (!self.firstOperand.dimensions.count && !self.secondOperand.dimensions.count)
		{ // scalar-scalar division
			if (!self.firstOperand.isComplex && !self.secondOperand.isComplex) {
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					GLFloat *a = (GLFloat *) fOperand.bytes;
					GLFloat *b = (GLFloat *) sOperand.bytes;
					GLFloat *c = (GLFloat *) result.mutableBytes;
					*c = (*a) / (*b);
				};
                self.graphvisDescription = @"division (real scalar, real scalar)";
			} else {
				[NSException raise: @"MethodNotImplemented" format: @"This case has not yet been implemented."];
			}
		}
		else if (!self.firstOperand.dimensions.count && self.secondOperand.dimensions.count)
		{ // scalar-vector multiplication
			if (!self.firstOperand.isComplex && !self.secondOperand.isComplex)  {
				NSUInteger nDataElements = self.secondOperand.nDataElements;
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					vGL_svdiv( (GLFloat *)fOperand.bytes, (GLFloat *)sOperand.bytes, 1, result.mutableBytes, 1, nDataElements);
				};
                self.graphvisDescription = @"division (real scalar, real vector)";
			} else {
				[NSException raise: @"MethodNotImplemented" format: @"This case has not yet been implemented."];
			}
		}
		else if (self.firstOperand.dimensions.count && !self.secondOperand.dimensions.count)
		{ // vector-scalar multiplication
			if (!self.firstOperand.isComplex && !self.secondOperand.isComplex)  {
				NSUInteger nDataElements = self.firstOperand.nDataElements;
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					vGL_svdiv( (GLFloat *)sOperand.bytes, (GLFloat *)fOperand.bytes, 1, result.mutableBytes, 1, nDataElements);
				};
                self.graphvisDescription = @"division (real vector, real scalar)";
			} else {
				[NSException raise: @"MethodNotImplemented" format: @"This case has not yet been implemented."];
			}
		}
		else
		{
			if (!self.firstOperand.isComplex && !self.secondOperand.isComplex) {
				NSUInteger nDataElements = self.secondOperand.nDataElements;
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					vGL_vdiv( (GLFloat *) sOperand.bytes, 1, (GLFloat *) fOperand.bytes, 1, result.mutableBytes, 1, nDataElements);
				};
                self.graphvisDescription = @"division (real, real)";
			} else if (self.firstOperand.isComplex && self.secondOperand.isComplex && useComplexArithmetic == NO) {
				NSUInteger nDataPoints = self.result.nDataPoints;
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					
					GLSplitComplex leftComplex = splitComplexFromData( fOperand );
					GLSplitComplex rightComplex = splitComplexFromData( sOperand );
					GLSplitComplex destComplex = splitComplexFromData( result );
					
					vGL_vdiv( rightComplex.realp, 1, leftComplex.realp, 1, destComplex.realp, 1, nDataPoints);
					vGL_vdiv( rightComplex.imagp, 1, leftComplex.imagp, 1, destComplex.imagp, 1, nDataPoints);
				};
                self.graphvisDescription = @"element-wise division (complex, complex)";
			} else if (self.firstOperand.isComplex && !self.secondOperand.isComplex && useComplexArithmetic == YES ) {
				NSUInteger nDataPoints = self.result.nDataPoints;
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					
					GLSplitComplex leftComplex = splitComplexFromData( fOperand );
					GLFloat *right = (GLFloat *)sOperand.bytes;
					GLSplitComplex destComplex = splitComplexFromData( result );
					
					// C = B/A
					vGL_vdiv( right, 1, leftComplex.realp, 1, destComplex.realp, 1, nDataPoints);
					vGL_vdiv( right, 1, leftComplex.imagp, 1, destComplex.imagp, 1, nDataPoints);
				};
                self.graphvisDescription = @"division (complex, complex)";
			}
			else {
				[NSException raise: @"MethodNotImplemented" format: @"This case has not yet been implemented."];
			}
		}
		
		
		
	}

	return self;
}

- (BOOL) isEqualToOperation: (id) otherOperation {
    if ( ![super isEqualToOperation: otherOperation] )  {
        return NO;
    }
    
    GLDivisionOperation * op = otherOperation;
    if (self.useComplexArithmetic != op.useComplexArithmetic) {
        return NO;
    }
    
    return YES;
}

@end

/************************************************/
/*		GLDotProductOperation					*/
/************************************************/

@implementation GLDotProductOperation

- (id) initWithFirstOperand: (GLVariable *) fOperand secondOperand: (GLVariable *) sOperand {
	
	BOOL complexResult = fOperand.isComplex || sOperand.isComplex;
	GLDataFormat format = complexResult ? kGLSplitComplexDataFormat : kGLRealDataFormat;
	GLVariable *resultVariable = [GLVariable variableOfType: format withDimensions: [NSArray array] forEquation: fOperand.equation];
	
	if (( self = [super initWithResult: resultVariable firstOperand:fOperand secondOperand: sOperand] ))
	{
		self.result.isPurelyReal = self.firstOperand.isPurelyReal && self.secondOperand.isPurelyReal;
		NSUInteger nDataPoints = self.firstOperand.nDataPoints;
		
		if (!self.firstOperand.isComplex && !self.secondOperand.isComplex) {
			self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
				vGL_dotpr( fOperand.bytes, 1, sOperand.bytes, 1, result.mutableBytes, nDataPoints);
			};
            self.graphvisDescription = @"dot (real, real)";
		} else if (self.firstOperand.isComplex && self.secondOperand.isComplex) {
			
			self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
				
				GLSplitComplex leftComplex = splitComplexFromData( fOperand );
				GLSplitComplex rightComplex = splitComplexFromData( sOperand );
				GLSplitComplex destComplex = splitComplexFromData( result );
				
				vGL_zdotpr( &leftComplex, 1, &rightComplex, 1, &destComplex, nDataPoints );
			};
            self.graphvisDescription = @"dot (complex, complex)";
		} else {
			[NSException raise: @"MethodNotImplemented" format: @"This case has not yet been implemented."];
		}
		
	}
	
	return self;
}

@end

/************************************************/
/*		GLSetVariableValueOperation				*/
/************************************************/

@implementation GLSetVariableValueOperation

- (id) initWithVectorOperand: (GLVariable *) variable scalarVariableOperand: (GLVariable *) aScalarVariable indexString: (NSString *) indexString
{
    NSArray *ranges = [GLDimension rangesFromIndexString: indexString usingDimensions: variable.dimensions];
	
	BOOL complexResult = variable.isComplex || aScalarVariable.isComplex;
	if (complexResult) [NSException raise: @"MethodNotImplemented" format: @"Complex numbers not implemented here."];
	
	GLDataFormat format = complexResult ? kGLSplitComplexDataFormat : kGLRealDataFormat;
	GLVariable *resultVariable = [GLVariable variableOfType: format withDimensions: variable.dimensions forEquation: variable.equation];
	
    if (( self = [super initWithResult: resultVariable firstOperand: variable secondOperand: aScalarVariable]))
	{
        self.indexString = indexString;
        
        NSUInteger numBytes = self.result.nDataElements*sizeof(GLFloat);
        
        self.result.name = variable.name;
		self.result.isPurelyReal = self.firstOperand.isPurelyReal;
		self.result.isPurelyImaginary = self.firstOperand.isPurelyImaginary;
        
        if ( variable.dimensions.count == 1 )
		{
            NSRange fastRange = [ranges[0] rangeValue];
            
            NSUInteger startIndex = fastRange.location;
            NSUInteger fastIndexLength = fastRange.length;
            
            self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand)  {
                GLFloat *toData = (GLFloat *) result.mutableBytes;
                // first copy the data
                memcpy( result.mutableBytes, fOperand.bytes,  numBytes );
                // then replace the value at the desired indices
                vGL_vfill( (GLFloat *) sOperand.bytes, &toData[startIndex], 1, fastIndexLength);
            };
            self.graphvisDescription = [NSString stringWithFormat: @"set leftVar=rightVar (1 dim)"];
		}
		else if ( variable.dimensions.count == 2 )
		{
			NSRange fastRange = [ranges[1] rangeValue];
            NSUInteger fastDimLength = [variable.dimensions[1] nPoints];
            
            NSRange slowRange = [ranges[0] rangeValue];
			
            self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
				// first copy the data
				memcpy( result.mutableBytes, fOperand.bytes,  numBytes );
				
                GLFloat *toData = (GLFloat *) result.mutableBytes;
                dispatch_apply(slowRange.length, dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0), ^(size_t iteration) {
                    // then replace the value at the desired indices
                    vGL_vfill( (GLFloat *)sOperand.bytes, &(toData[(slowRange.location + iteration)*fastDimLength + fastRange.location]), 1, fastRange.length);
                    
                });
            };
            self.graphvisDescription = [NSString stringWithFormat: @"set=rightVar (1 dim)"];
        }
        else if ( variable.dimensions.count == 3 )
        {
            NSUInteger ny = [variable.dimensions[1] nPoints];
            NSUInteger nz = [variable.dimensions[2] nPoints];
            
            NSRange xrange = [ranges[0] rangeValue];
            NSRange yrange = [ranges[1] rangeValue];
            NSRange zrange = [ranges[2] rangeValue];
            
            self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
				// first copy the data
				memcpy( result.mutableBytes, fOperand.bytes,  numBytes );
				GLFloat *aScalar = (GLFloat *) sOperand.bytes;
                GLFloat *toData = (GLFloat *) result.mutableBytes;
                for (NSUInteger i=xrange.location; i<xrange.location+xrange.length; i++) {
                    for (NSUInteger j=yrange.location; j<yrange.location+yrange.length; j++) {
                        for (NSUInteger k=zrange.location; k<zrange.location+zrange.length; k++) {
                            toData[(i*ny+j)*nz+k] = *aScalar;
                        }
                    }
                }
            };
        }
    }
    return self;
}

- (BOOL) isEqualToOperation: (id) otherOperation {
    if ( ![super isEqualToOperation: otherOperation] )  {
        return NO;
    }
    
    GLSetVariableValueOperation * op = otherOperation;
    if ( [self.indexString isEqualToString: op.indexString]) {
        return NO;
    }
    
    return YES;
}

@end


