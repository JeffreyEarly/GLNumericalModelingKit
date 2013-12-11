//
//  GLLinearTransformationOperations.m
//  GLNumericalModelingKitCopy
//
//  Created by Jeffrey J. Early on 1/18/13.
//
//

#import "GLLinearTransformationOperations.h"
#import "GLLinearTransform.h"

#import "GLMemoryPool.h"

#import <Accelerate/Accelerate.h>

/************************************************/
/*		GLSingleDiagonalTransformOperation		*/
/************************************************/

#pragma mark -
#pragma mark GLSingleDiagonalTransformOperation
#pragma mark

@implementation GLSingleDiagonalTransformOperation

- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform function: (GLFunction *) function
{
	BOOL isComplex = linearTransform.isComplex || function.isComplex;
	GLDataFormat format = isComplex ? kGLSplitComplexDataFormat : kGLRealDataFormat;
	GLFunction *result = [GLFunction functionOfType: format withDimensions: linearTransform.toDimensions forEquation: linearTransform.equation];
	
    if (linearTransform.name && function.name) {
        result.name = [NSString stringWithFormat: @"%@_%@", function.name, linearTransform.name];
    }
    
	if (( self = [super initWithResult: @[result] operand: @[linearTransform, function]] )) {
		result.isPurelyReal = (linearTransform.isPurelyReal && function.isPurelyReal) || (linearTransform.isPurelyImaginary && function.isPurelyImaginary);
		result.isPurelyImaginary= (linearTransform.isPurelyReal && function.isPurelyImaginary) || (linearTransform.isPurelyImaginary && function.isPurelyReal);
		
		// First we determine how many elements over we need to shift the result.
		// This is necessary because shifting between a sine and cosine basis results
		// in a different set of wavenumbers available. A positive shift indicates that the
		// result array shifts forward, a negative shift indicates that the operator shifts forward.
		
		// Going from a cosine to sine should be a superdiagonal, meaning that we reference the first element of the operand, and the second element of transformation.
		// In otherwards, we shift the operator forward.
		
		NSInteger totalShift = 0;
		GLMatrixDescription *matrix = linearTransform.matrixDescription;
		for (NSUInteger i=0; i<matrix.nDimensions; i++) {
			if (matrix.strides[i].matrixFormat == kGLSuperdiagonalMatrixFormat) {
				totalShift -= matrix.strides[i].stride;
			} else if (matrix.strides[i].matrixFormat == kGLSubdiagonalMatrixFormat) {
				totalShift += matrix.strides[i].stride;
			}  else if (matrix.strides[i].matrixFormat == kGLDiagonalMatrixFormat) {
				totalShift += 0;
			} else {
				[NSException raise: @"BadFormat" format: @"This operation type can only transform with matrices in a diagonal format."];
			}
		}
		
		BOOL shiftResult = totalShift > 0 ? YES : NO; // Yes if we shift the result and transform, no if we shift the operand
		totalShift = labs(totalShift); // total number of elements to shift by
		NSInteger totalShiftBytes = totalShift*sizeof(GLFloat); // total number of bytes to shift by
		NSInteger nDataElements = result.nDataElements - totalShift;
		NSInteger nDataPoints = result.nDataPoints - totalShift;
		NSInteger totalEndShiftBytes = nDataPoints*sizeof(GLFloat);
		
		if ( !linearTransform.isComplex && !function.isComplex )
		{
			if (shiftResult) {
				self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
					NSMutableData *result = resultArray[0];
					NSMutableData *transform = operandArray[0];
					NSMutableData *function = operandArray[1];
					vGL_vmul( transform.bytes+totalShiftBytes, 1, function.bytes, 1, result.mutableBytes+totalShiftBytes, 1, nDataElements);
					vGL_vclr( result.mutableBytes, 1, totalShift );
				};
			} else {
				self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
					NSMutableData *result = resultArray[0];
					NSMutableData *transform = operandArray[0];
					NSMutableData *function = operandArray[1];
					vGL_vmul( transform.bytes, 1, function.bytes+totalShiftBytes, 1, result.mutableBytes, 1, nDataElements);
					vGL_vclr( result.mutableBytes+totalEndShiftBytes, 1, totalShift );
				};
			}
		}
		else if ( !linearTransform.isComplex && function.isComplex )
		{
			if (shiftResult) {
				self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
					GLSplitComplex result = splitComplexFromData( resultArray[0] );
					NSMutableData *transform = operandArray[0];
					GLSplitComplex function = splitComplexFromData( operandArray[1] );
					
					vGL_vmul( transform.bytes+totalShiftBytes, 1, function.realp, 1, result.realp+totalShiftBytes, 1, nDataPoints);
					vGL_vclr( result.realp, 1, totalShift );
					
					vGL_vmul( transform.bytes+totalShiftBytes, 1, function.imagp, 1, result.imagp+totalShiftBytes, 1, nDataPoints);
					vGL_vclr( result.realp, 1, totalShift );
				};
			} else {
				self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
					GLSplitComplex result = splitComplexFromData( resultArray[0] );
					NSMutableData *transform = operandArray[0];
					GLSplitComplex function = splitComplexFromData( operandArray[1] );
					
					vGL_vmul( transform.bytes, 1, function.realp+totalShiftBytes, 1, result.realp, 1, nDataPoints);
					vGL_vclr( result.realp+totalEndShiftBytes, 1, totalShift );
					
					vGL_vmul( transform.bytes, 1, function.imagp+totalShiftBytes, 1, result.imagp, 1, nDataPoints);
					vGL_vclr( result.realp+totalEndShiftBytes, 1, totalShift );
				};
			}
		} else {
			[NSException raise: @"MethodNotImplemented" format: @"MethodNotImplemented"];
		}
	}
    return self;
}

- (BOOL) canOperateInPlace {
	return YES;
}

@end

/************************************************/
/*		GLTriadiagonalTransformOperation		*/
/************************************************/

#pragma mark -
#pragma mark GLTriadiagonalTransformOperation
#pragma mark

@implementation GLTriadiagonalTransformOperation

- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform function: (GLFunction *) function
{
    if ( ![linearTransform.fromDimensions isEqualToArray: function.dimensions] ) {
        [NSException raise: @"DimensionsNotEqualException" format: @"From dimensions of the linear operator must equal the operand vector."];
    }
    
    NSUInteger triIndex = NSNotFound;
	NSUInteger lastNonTrivialNonTriIndex = NSNotFound;
	NSUInteger lastTrivialIndex = NSNotFound;
	NSUInteger totalTrivialPoints = 1;
	NSUInteger numTriIndices = 0;
    for ( NSUInteger index=0; index < linearTransform.matrixFormats.count; index++) {
        NSNumber *num = linearTransform.matrixFormats[index];
		
        if ([num unsignedIntegerValue] == kGLTridiagonalMatrixFormat) {
            triIndex = index;
			numTriIndices++;
        } else if ([num unsignedIntegerValue] == kGLIdentityMatrixFormat) {
            lastTrivialIndex = index;
			totalTrivialPoints *= [linearTransform.fromDimensions[index] nPoints];
        }
        
        if ([num unsignedIntegerValue] != kGLIdentityMatrixFormat && [num unsignedIntegerValue] != kGLTridiagonalMatrixFormat ) {
			lastNonTrivialNonTriIndex = index;
		}
		
		
    }
	
    if ( numTriIndices != 1 ) {
        [NSException raise: @"TridiagonalIndexNotFound" format: @"Unable to find a tridiagonal index."];
    }
    
	BOOL isComplex = linearTransform.isComplex || function.isComplex;
	GLDataFormat format = isComplex ? kGLSplitComplexDataFormat : kGLRealDataFormat;
	GLFunction *result = [GLFunction functionOfType: format withDimensions: linearTransform.toDimensions forEquation: linearTransform.equation];
	
    if (linearTransform.name && function.name) {
        result.name = [NSString stringWithFormat: @"%@_%@", function.name, linearTransform.name];
    }
    
	if (( self = [super initWithResult: @[result] operand: @[linearTransform, function]] )) {
        GLMatrixDescription *matrixDescription = linearTransform.matrixDescription;
        
        dispatch_queue_t globalQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
        
        NSUInteger inNDiagonalPoints = matrixDescription.strides[triIndex].nDiagonalPoints;
        NSUInteger inElementStride = matrixDescription.strides[triIndex].stride;
        NSUInteger inDiagonalStride = matrixDescription.strides[triIndex].diagonalStride;
		
		// How many 'inner loop' steps we need to take depends on how many other nontrivial dimensions there are.
        NSUInteger inEquationStride = lastNonTrivialNonTriIndex == NSNotFound ? 0 : matrixDescription.strides[lastNonTrivialNonTriIndex].stride;
		NSUInteger totalEquations = linearTransform.nDataPoints / matrixDescription.strides[triIndex].nPoints;
		
		// Now we need the strides to match up to the inner loop
        NSUInteger outElementStride = result.matrixDescription.strides[triIndex].stride;
        NSUInteger outEquationStride = lastNonTrivialNonTriIndex == NSNotFound ? 0 : result.matrixDescription.strides[lastNonTrivialNonTriIndex].stride;
        
		// Finally, we need the outer loop strides and totals
		NSUInteger totalOuterLoops = totalTrivialPoints;
		NSUInteger outerLoopStride = lastTrivialIndex == NSNotFound ? 0 : result.matrixDescription.strides[lastTrivialIndex].stride;
		
		// This operation has two loops.
		// The inner loop walks over non-trivial elements of the linear operator, while
		// the out loop walks over trivial elements of the linear operator, really just to new x and b elements
        
        self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
            
			dispatch_apply(totalOuterLoops, globalQueue, ^(size_t outerIteration) {
				
				NSInteger outerOffset = outerIteration*outerLoopStride;
				
				dispatch_apply(totalEquations, globalQueue, ^(size_t iteration) {
					
					NSInteger inEquationPos = iteration*inEquationStride;
					NSInteger outEquationPos = outerOffset + iteration*outEquationStride;
					
					GLFloat *f = (GLFloat *) [operandArray[0] bytes];
					
					// plus the offset!!!!
					GLFloat *a = &(f[0*inDiagonalStride]);
					GLFloat *b = &(f[1*inDiagonalStride]);
					GLFloat *c = &(f[2*inDiagonalStride]);
					
					GLFloat *x = (GLFloat *) [operandArray[1] bytes];
					GLFloat *d = (GLFloat *) [resultArray[0] bytes];
					
					NSUInteger iIn0 = inEquationPos + 0*inElementStride;
					NSUInteger iOut0 = outEquationPos + 0*outElementStride;
					NSUInteger iIn1 = inEquationPos + 1*inElementStride;
					NSUInteger iOut1 = outEquationPos + 1*outElementStride;
					NSUInteger iOut2 = outEquationPos + 2*outElementStride;
					
					// The first multiplication doesn't involve the subdiagonal.
					// d(0) = b(0)*x(0) + c(0)*x(1)
					d[iOut0] = b[iIn0]*x[iOut0] + c[iIn0]*x[iOut1];
					
					// d(i) = a(i)*x(i-1) + b(i)*x(i)
					vGL_vmma( &(a[iIn1]), inElementStride, &(x[iOut0]), outElementStride, &(b[iIn1]), inElementStride, &(x[iOut1]), outElementStride, &(d[iOut1]), outElementStride, inNDiagonalPoints-2);
					
					// d(i) = c(i)*x(i+1) + d(i)
					vGL_vma( &(c[iIn1]), inElementStride, &(x[iOut2]), outElementStride, &(d[iOut1]), outElementStride, &(d[iOut1]), outElementStride, inNDiagonalPoints-2);
					
					NSUInteger iInN1 = inEquationPos + (inNDiagonalPoints-1)*inElementStride;
					NSUInteger iOutN1 = outEquationPos + (inNDiagonalPoints-1)*outElementStride;
					NSUInteger iOutN2 = outEquationPos + (inNDiagonalPoints-2)*outElementStride;
										
					d[iOutN1] = a[iInN1]*x[iOutN2] + b[iInN1]*x[iOutN1];
				});
			});
        };
        
    }
    return self;
}

@end

/************************************************/
/*		GLDenseMatrixTransformOperation		*/
/************************************************/

#pragma mark -
#pragma mark GLDenseMatrixTransformOperation
#pragma mark

@implementation GLDenseMatrixTransformOperation

- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform function: (GLFunction *) function
{
    if ( ![linearTransform.fromDimensions isEqualToArray: function.dimensions] ) {
        [NSException raise: @"DimensionsNotEqualException" format: @"From dimensions of the linear operator must equal the operand vector."];
    }
    
    for (NSNumber *format in linearTransform.matrixFormats) {
        if (format.unsignedIntegerValue != kGLDenseMatrixFormat) {
            [NSException raise: @"MatrixWrongFormat" format: @"This operation can only be performed with a dense matrix."];
        }
    }
    
    if (linearTransform.matrixDescription.nDimensions != 1) {
        [NSException raise: @"MatrixWrongFormat" format: @"We can only do one dimensional matrices at the moment."];
    }

	BOOL isComplex = linearTransform.isComplex || function.isComplex;
	GLDataFormat format = isComplex ? kGLSplitComplexDataFormat : kGLRealDataFormat;
	GLFunction *result = [GLFunction functionOfType: format withDimensions: linearTransform.toDimensions forEquation: linearTransform.equation];
	
    if (linearTransform.name && function.name) {
        result.name = [NSString stringWithFormat: @"%@_%@", function.name, linearTransform.name];
    }
    
	if (( self = [super initWithResult: @[result] operand: @[linearTransform, function]] )) {
        GLMatrixDescription *matrixDescription = linearTransform.matrixDescription;
        
//        int M = (int) matrixDescription.strides[0].nRows;
//        int N = (int) matrixDescription.strides[0].nColumns;
//        self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
//            cblas_sgemv( CblasRowMajor,  CblasNoTrans, M, N, 1.0, fOperand.bytes, M, sOperand.bytes, 1.0, 1.0, result.mutableBytes, 1);
//        };
		
		int M = (int) matrixDescription.strides[0].nRows;
        int N = (int) 1;
		int K = (int) matrixDescription.strides[0].nColumns;
		if ( !linearTransform.isComplex && !function.isComplex)
		{	// C = A.X
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *A = (GLFloat *) [operandArray[0] bytes];
				GLFloat *B = (GLFloat *) [operandArray[1] bytes];
				GLFloat *C = (GLFloat *) [resultArray[0] bytes];
				vDSP_mmul( A, 1, B, 1, C, 1, M, N, K);
			};
		}
		else if ( linearTransform.isComplex && !function.isComplex)
		{	// (A+iB).(X) = A.X + iB.X
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLSplitComplex A = splitComplexFromData( operandArray[0] );
				GLFloat *B = (GLFloat *) [operandArray[1] bytes];
				GLSplitComplex C = splitComplexFromData( resultArray[0] );
				
				vDSP_mmul( A.realp, 1, B, 1, C.realp, 1, M, N, K);
				vDSP_mmul( A.imagp, 1, B, 1, C.imagp, 1, M, N, K);
			};
		}
		else if ( !linearTransform.isComplex && function.isComplex)
		{	// A.(X+iY) = A.X + iA.Y
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *A = (GLFloat *) [operandArray[0] bytes];
				GLSplitComplex B = splitComplexFromData( operandArray[1] );
				GLSplitComplex C = splitComplexFromData( resultArray[0] );
				
				vDSP_mmul( A, 1, B.realp, 1, C.realp, 1, M, N, K);
				vDSP_mmul( A, 1, B.imagp, 1, C.imagp, 1, M, N, K);
			};
		}
		else if ( linearTransform.isComplex && function.isComplex)
		{
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLSplitComplex A = splitComplexFromData( operandArray[0] );
				GLSplitComplex B = splitComplexFromData( operandArray[1] );
				GLSplitComplex C = splitComplexFromData( resultArray[0] );
				
				vDSP_zmmul( &A, 1, &B, 1, &C, 1, M, N, K);
			};
		}
    }
    return self;
}

@end

/************************************************/
/*		GLTriadiagonalSolverOperation			*/
/************************************************/

#pragma mark -
#pragma mark GLTriadiagonalSolverOperation
#pragma mark

@implementation GLTriadiagonalSolverOperation

- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform function: (GLFunction *) function
{
    if ( ![linearTransform.toDimensions isEqualToArray: function.dimensions] ) {
        [NSException raise: @"DimensionsNotEqualException" format: @"Destination dimensions of the linear operator must equal the resultant vector."];
    }
    
    NSUInteger triIndex = NSNotFound;
	NSUInteger lastNonTriIndex = NSNotFound;
	NSUInteger numTriIndices = 0;
    for ( NSUInteger index=0; index < linearTransform.matrixFormats.count; index++) {
        NSNumber *num = linearTransform.matrixFormats[index];
        if ([num unsignedIntegerValue] == kGLTridiagonalMatrixFormat) {
            triIndex = index;
			numTriIndices++;
        }
        
        if ([num unsignedIntegerValue] != kGLIdentityMatrixFormat && [num unsignedIntegerValue] != kGLTridiagonalMatrixFormat ) {
			lastNonTriIndex = index;
		}
    }
	
    if ( numTriIndices != 1 ) {
        [NSException raise: @"TridiagonalIndexNotFound" format: @"Unable to find a tridiagonal index."];
    }
	
	BOOL isComplex = linearTransform.isComplex || function.isComplex;
	GLDataFormat format = isComplex ? kGLSplitComplexDataFormat : kGLRealDataFormat;
	GLFunction *result = [GLFunction functionOfType: format withDimensions: linearTransform.fromDimensions forEquation: linearTransform.equation];
    
	if (( self = [super initWithResult: @[result] operand: @[linearTransform, function]] )) {
        
		
        GLMatrixDescription *matrixDescription = linearTransform.matrixDescription;
        
#warning this needs the same outer loop as the transform below.
        
        // This buffer is 3 times as large as it needs to be, but it makes bookkeeping easier.
        NSMutableData *buffer = [[GLMemoryPool sharedMemoryPool] dataWithLength: result.nDataPoints*sizeof(GLFloat)];
		
        dispatch_queue_t globalQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
        
        NSUInteger inNDiagonalPoints = matrixDescription.strides[triIndex].nDiagonalPoints;
        NSUInteger inElementStride = matrixDescription.strides[triIndex].stride;
        NSUInteger inDiagonalStride = matrixDescription.strides[triIndex].diagonalStride;
        NSUInteger inEquationStride = lastNonTriIndex == NSNotFound ? 0 : matrixDescription.strides[lastNonTriIndex].stride;
        NSUInteger outElementStride = result.matrixDescription.strides[triIndex].stride;
        NSUInteger outEquationStride = lastNonTriIndex == NSNotFound ? 0 : result.matrixDescription.strides[lastNonTriIndex].stride;
        NSUInteger totalEquations = linearTransform.nDataPoints / matrixDescription.strides[triIndex].nPoints;
        
        self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
            
            dispatch_apply(totalEquations, globalQueue, ^(size_t iteration) {
				
                NSInteger inEquationPos = iteration*inEquationStride;
                NSInteger outEquationPos = iteration*outEquationStride;
                
                GLFloat *f = (GLFloat *) [operandArray[0] bytes];
                
                // plus the offset!!!!
                GLFloat *a = &(f[0*inDiagonalStride]);
                GLFloat *b = &(f[1*inDiagonalStride]);
                GLFloat *c = &(f[2*inDiagonalStride]);
                
                GLFloat *cprime = (GLFloat *) buffer.mutableBytes;
                
                GLFloat *d = (GLFloat *) [operandArray[1] bytes];
                GLFloat *x = (GLFloat *) [resultArray[0] bytes];
                
                NSUInteger iIn = inEquationPos + 0*inElementStride;
                NSUInteger iOut = outEquationPos + 0*outElementStride;
                
                cprime[iOut] = c[iIn] / b[iIn];
                x[iOut] = d[iOut] / b[iIn];
                
                for (NSInteger i=1; i<inNDiagonalPoints; i++) {
                    iIn = inEquationPos + i*inElementStride;
                    iOut = outEquationPos + i*outElementStride;
                    NSUInteger iOutMinus1 = outEquationPos + (i-1)*outElementStride;
                    
                    GLFloat m = 1.0 / ( b[iIn] - a[iIn] * cprime[iOutMinus1]);
                    cprime[iOut] = c[iIn] * m;
                    x[iOut] = (d[iOut] - a[i]*x[iOutMinus1])*m;
                }
                
                for (NSInteger i=inNDiagonalPoints-2; i >= 0; i-- ) {
                    iOut = outEquationPos + i*outElementStride;
                    NSUInteger iOutPlus1 = outEquationPos + (i+1)*outElementStride;
                    
                    x[iOut] = x[iOut] - cprime[iOut] * x[iOutPlus1];
                }
                
				//                if (iteration == 2) { for (NSInteger i=0; i<inNDiagonalPoints; i++) {
				//                    iIn = inEquationPos + i*inElementStride;
				//                    iOut = outEquationPos + i*outElementStride;
				//
				//                    printf("i=%ld (a,b,c,d,x,c')=(%f,%f,%f,%f, %f, %f)\n", i, a[iIn], b[iIn], c[iIn], d[iOut], x[iOut], cprime[iOut]);
				//                }}
            });
        };
        
    }
    return self;
}

@end


/************************************************/
/*		GLDenseMatrixSolver						*/
/************************************************/

#pragma mark -
#pragma mark GLDenseMatrixSolver
#pragma mark

@implementation GLDenseMatrixSolver

- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform function: (GLFunction *) function
{
    if ( ![linearTransform.toDimensions isEqualToArray: function.dimensions] ) {
        [NSException raise: @"DimensionsNotEqualException" format: @"Destination dimensions of the linear operator must equal the resultant vector."];
    }
    
    NSUInteger denseIndex = NSNotFound;
	NSUInteger lastNonTrivialNonDenseIndex = NSNotFound;
	NSUInteger lastTrivialIndex = NSNotFound;
	NSUInteger totalTrivialPoints = 1;
	NSUInteger numDenseIndices = 0;
    for ( NSUInteger index=0; index < linearTransform.matrixFormats.count; index++) {
        NSNumber *num = linearTransform.matrixFormats[index];
		
        if ([num unsignedIntegerValue] == kGLDenseMatrixFormat) {
            denseIndex = index;
			numDenseIndices++;
        } else if ([num unsignedIntegerValue] == kGLIdentityMatrixFormat) {
            lastTrivialIndex = index;
			totalTrivialPoints *= [linearTransform.fromDimensions[index] nPoints];
        }
        if ([num unsignedIntegerValue] != kGLIdentityMatrixFormat && [num unsignedIntegerValue] != kGLDenseMatrixFormat ) {
			lastNonTrivialNonDenseIndex = index;
		}
    }
	
    if ( numDenseIndices != 1 ) {
        [NSException raise: @"BadInputFormat" format: @"The GLDenseMatrixSolver can only take exactly one dense dimension."];
    }
	
    
	BOOL isComplex = linearTransform.isComplex || function.isComplex;
	GLDataFormat format = isComplex ? kGLSplitComplexDataFormat : kGLRealDataFormat;
	GLFunction *result = [GLFunction functionOfType: format withDimensions: linearTransform.fromDimensions forEquation: linearTransform.equation];
    
	if (( self = [super initWithResult: @[result] operand: @[linearTransform, function]] )) {
        
        GLMatrixDescription *matrixDescription = linearTransform.matrixDescription;
        		
        dispatch_queue_t globalQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
		
		// How many 'inner loop' steps we need to take depends on how many other nontrivial dimensions there are.
        NSUInteger inEquationStride = lastNonTrivialNonDenseIndex == NSNotFound ? 0 : matrixDescription.strides[lastNonTrivialNonDenseIndex].stride;
		NSUInteger totalEquations = linearTransform.nDataPoints / matrixDescription.strides[denseIndex].nPoints;
		
		// Now we need the strides to match up to the inner loop
        NSUInteger outElementStride = result.matrixDescription.strides[denseIndex].stride;
        NSUInteger outEquationStride = lastNonTrivialNonDenseIndex == NSNotFound ? 0 : result.matrixDescription.strides[lastNonTrivialNonDenseIndex].stride;
        
		// Finally, we need the outer loop strides and totals
		// The totalTrivialPoints is the number of times we need to repeat a point
		// trivialPointStride is the distance between those trivial points
		// totalNonTrivialPoints is the number of points we need to copy
		NSUInteger totalNonTrivialPoints = result.nDataPoints / totalTrivialPoints;
		NSUInteger trivialPointStride = lastTrivialIndex == NSNotFound ? 0 : result.matrixDescription.strides[lastTrivialIndex].stride;
		
		// This operation has two loops.
		// The inner loop walks over non-trivial elements of the linear operator, while
		// the out loop walks over trivial elements of the linear operator, really just to new x and b elements
        
		NSUInteger N = matrixDescription.strides[denseIndex].nRows;
		
		if ( lastNonTrivialNonDenseIndex != NSNotFound && lastNonTrivialNonDenseIndex > denseIndex) {
			NSLog(@"Warning! The dense dimension is not the last non-trivial dimensions, which results in slower solutions times.");
		}
		// sgesv is only capable of handling input arrays with a stride of 1. For us, this means that if that lastNonTrivialNonDenseIndex is greater
		// than the denseIndex, we have to copy the elements to a temporary, contiguous buffer.
		
		// sgesv also overwrites the b vector (in Mx=b) with the output x. So we need to copy the output as well.
		
        self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			
			dispatch_apply(totalEquations, globalQueue, ^(size_t iteration) {
				
				NSInteger inEquationPos = iteration*inEquationStride;
				NSInteger outEquationPos = iteration*outEquationStride;
				
				NSMutableData *ipiv = [[GLMemoryPool sharedMemoryPool] dataWithLength: N*sizeof(__CLPK_integer)];
				__CLPK_integer n = (__CLPK_integer) N;
				__CLPK_integer nrhs = 1;
				__CLPK_integer info;
				
				GLFloat *MData = (GLFloat *) [operandArray[0] bytes];
				GLFloat *bData = (GLFloat *) [operandArray[1] bytes];
				GLFloat *xData = (GLFloat *) [resultArray[0] bytes];
				
				GLFloat *M;
				GLFloat *b;
				GLFloat *x;
				if ( lastNonTrivialNonDenseIndex != NSNotFound && lastNonTrivialNonDenseIndex > denseIndex) {
					[NSException raise: @"NotYetImplemented" format: @"This case is not yet implemented."];
					return;
				} else {
					M = &(MData[inEquationPos]);
					b = &(bData[outEquationPos]);
					x = &(xData[outEquationPos]);
					memcpy( x, b, n*n*sizeof(GLFloat));
				}
				
				sgesv_( &n, &nrhs, M, &n, ipiv.mutableBytes, x, &n, (__CLPK_integer *) &info );
				
				if (info != 0) {
					printf("sgesv failed with error code %d\n", (int)info);
				}
				
				[[GLMemoryPool sharedMemoryPool] returnData: ipiv];
				
				if (totalTrivialPoints) {
					dispatch_apply(totalNonTrivialPoints, globalQueue, ^(size_t outIteration) {
						GLFloat *theOutput = (GLFloat *) [resultArray[0] bytes];
						vGL_vfill( &theOutput[ outIteration*outElementStride], &theOutput[outIteration*outElementStride], trivialPointStride, totalTrivialPoints);
					});
					
				}
				
			});
        };
        
    }
    return self;
}

@end

/************************************************/
/*		GLMatrixMatrixMultiplicationOperation   */
/************************************************/

#pragma mark -
#pragma mark GLMatrixMatrixMultiplicationOperation
#pragma mark

@implementation GLMatrixMatrixMultiplicationOperation

// This is copy and pasted from the superclass, needs to be properly retooled.
- (id) initWithFirstOperand: (GLLinearTransform *) A secondOperand: (GLLinearTransform *) B;
{
    
    if ( ![A.fromDimensions isEqualToArray: B.toDimensions] ) {
        [NSException raise: @"DimensionsNotEqualException" format: @"fromDimensions of A, must equal the toDimensions of B."];
    }
    
    for ( NSUInteger index=0; index < A.matrixFormats.count; index++) {
        NSNumber *format = A.matrixFormats[index];
        if (format.unsignedIntegerValue != kGLDenseMatrixFormat) {
            [NSException raise: @"MatrixWrongFormat" format: @"This operation can only be performed with a dense matrix."];
        }
    }
    
    if (A.matrixDescription.nDimensions != 1) {
        [NSException raise: @"MatrixWrongFormat" format: @"We can only do one dimensional matrices at the moment."];
    }
	
	BOOL isComplex = A.isComplex || B.isComplex;
	GLDataFormat format = isComplex ? kGLSplitComplexDataFormat : kGLRealDataFormat;
	GLLinearTransform *result = [GLLinearTransform transformOfType: format withFromDimensions: B.fromDimensions toDimensions:A.toDimensions inFormat:B.matrixFormats forEquation:B.equation matrix: nil];
    
	if (( self = [super initWithResult: @[result] operand: @[A, B]] )) {
        
		//GLLinearTransform *C = result;
		        
        int M = (int) A.matrixDescription.strides[0].nRows;
        int N = (int) B.matrixDescription.strides[0].nColumns;
		int K = (int) A.matrixDescription.strides[0].nColumns;
		
//		int lda = (int) A.matrixDescription.strides[0].rowStride;
//        int ldb = (int) B.matrixDescription.strides[0].rowStride;
//		int ldc = (int) C.matrixDescription.strides[0].rowStride;
//		
//        self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
//			cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0, fOperand.bytes, lda, sOperand.bytes, ldb, 0.0, result.mutableBytes, ldc);
//        };
		
		if ( !A.isComplex && !B.isComplex)
		{	// C = A.X
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *A = (GLFloat *) [operandArray[0] bytes];
				GLFloat *B = (GLFloat *) [operandArray[1] bytes];
				GLFloat *C = (GLFloat *) [resultArray[0] bytes];
				vDSP_mmul( A, 1, B, 1, C, 1, M, N, K);
			};
		}
		else if ( A.isComplex && !B.isComplex)
		{	// (A+iB).(X) = A.X + iB.X
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLSplitComplex A = splitComplexFromData( operandArray[0] );
				GLFloat *B = (GLFloat *) [operandArray[1] bytes];
				GLSplitComplex C = splitComplexFromData( resultArray[0] );
				
				vDSP_mmul( A.realp, 1, B, 1, C.realp, 1, M, N, K);
				vDSP_mmul( A.imagp, 1, B, 1, C.imagp, 1, M, N, K);
			};
		}
		else if ( !A.isComplex && B.isComplex)
		{	// A.(X+iY) = A.X + iA.Y
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *A = (GLFloat *) [operandArray[0] bytes];
				GLSplitComplex B = splitComplexFromData( operandArray[1] );
				GLSplitComplex C = splitComplexFromData( resultArray[0] );
				
				vDSP_mmul( A, 1, B.realp, 1, C.realp, 1, M, N, K);
				vDSP_mmul( A, 1, B.imagp, 1, C.imagp, 1, M, N, K);
			};
		}
		else if ( A.isComplex && B.isComplex)
		{
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLSplitComplex A = splitComplexFromData(operandArray[0]);
				GLSplitComplex B = splitComplexFromData(operandArray[1]);
				GLSplitComplex C = splitComplexFromData(resultArray[0]);
				
				vDSP_zmmul( &A, 1, &B, 1, &C, 1, M, N, K);
			};
		}
		
    }
    return self;
}

@end

/************************************************/
/*		GLMatrixInversionOperation              */
/************************************************/

#pragma mark -
#pragma mark GLMatrixInversionOperation
#pragma mark

@implementation GLMatrixInversionOperation

- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform;
{
    if (linearTransform.matrixDescription.nDimensions != 1) {
        [NSException raise: @"MatrixWrongFormat" format: @"We can only do one dimensional matrices at the moment."];
    }
    
    GLLinearTransform *A = (GLLinearTransform *) linearTransform;
    if (A.fromDimensions.count != A.toDimensions.count ) {
        [NSException raise: @"MatrixWrongFormat" format: @"We can only do square matrices."];
    }
	
	GLDataFormat format = A.isComplex ? kGLSplitComplexDataFormat : kGLRealDataFormat;
	NSUInteger N = [A.fromDimensions[0] nPoints];
	GLLinearTransform *result = [GLLinearTransform transformOfType: format withFromDimensions: A.toDimensions toDimensions:A.fromDimensions inFormat:A.matrixFormats forEquation:A.equation matrix: nil];
    GLBuffer *buffer1 = [[GLBuffer alloc] initWithLength: N*sizeof(__CLPK_integer)];
	GLBuffer *buffer2 = [[GLBuffer alloc] initWithLength: N*N*sizeof(GLFloat)];
	NSArray *buffers = @[buffer1, buffer2];
	
	variableOperation op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
		GLFloat *A = (GLFloat *) [operandArray[0] bytes];
		GLFloat *C = (GLFloat *) [resultArray[0] bytes];
		NSMutableData *ipiv = bufferArray[0];
		NSMutableData *work = bufferArray[1];
		
		// clapack takes matrices in column-major format.
		// However, the transpose of the inverse is the inverse of the transpose---so we don't need to worry here.
		
		__CLPK_integer n = (__CLPK_integer) N;
		__CLPK_integer info;
		memcpy( C, A, n*n*sizeof(GLFloat));
		sgetrf_(&n, &n, C, &n, ipiv.mutableBytes, (__CLPK_integer *)&info);
		
		if (info != 0) {
			printf("sgetrf failed with error code %d\n", (int)info);
		}
		
		__CLPK_integer lwork = n*n;
		sgetri_(&n, C, &n, ipiv.mutableBytes, work.mutableBytes, &lwork, (__CLPK_integer *)&info);
		
		if (info != 0) {
			printf("sgetri failed with error code %d\n", (int)info);
		}
	};
	
	if (( self = [super initWithResult: @[result] operand: @[A] buffers: buffers operation: op] )) {
        
        
    }
    return self;
}

@end

/************************************************/
/*		GLMatrixEigensystemOperation            */
/************************************************/

#pragma mark -
#pragma mark GLMatrixEigensystemOperation
#pragma mark

@implementation GLMatrixEigensystemOperation

- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform;
{
    if (linearTransform.matrixDescription.nDimensions != 1) {
        [NSException raise: @"MatrixWrongFormat" format: @"We can only do one dimensional matrices at the moment."];
    }
    
    GLLinearTransform *A = (GLLinearTransform *) linearTransform;
    if (A.fromDimensions.count != A.toDimensions.count ) {
        [NSException raise: @"MatrixWrongFormat" format: @"We can only do square matrices."];
    }
	
	for (NSUInteger i=0; i<A.fromDimensions.count; i++) {
		if ( ![A.fromDimensions[i] isEqualToDimension: A.toDimensions[i]] ) {
			[NSException raise: @"MatrixWrongFormat" format: @"By assumption, a linear transformation must be an endomorphism to compute the eigensystem."];
		}
	}
	
	// We need to construct a *new* eigenbasis.
	// I'm not quite sure the right definitions to use.
	NSMutableArray *eigenbasis = [NSMutableArray array];
	for (GLDimension *dim in A.fromDimensions) {
		GLDimension *newDim = [[GLDimension alloc] initDimensionWithGrid: dim.gridType nPoints: dim.nPoints domainMin:dim.domainMin length: dim.domainLength];
		if (dim.name && ![dim.name isEqualToString: @""]) {
			newDim.name = [NSString stringWithFormat: @"%@_eigen", dim.name];
		} else {
			newDim.name = @"eigen";
		}
		[eigenbasis addObject: newDim];
	}
	
	// Should check whether or not the matrix is symmetric.
	GLVariable *eigenvalues = [GLFunction functionOfComplexTypeWithDimensions: eigenbasis forEquation: A.equation];
	GLLinearTransform *eigentransform = [GLLinearTransform transformOfType: kGLSplitComplexDataFormat withFromDimensions: eigenbasis toDimensions:eigenbasis inFormat:A.matrixFormats forEquation:A.equation matrix: nil];
    NSArray *results = @[eigenvalues, eigentransform];
	
	// http://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/sgeev_ex.c.htm
	
	NSUInteger N = [A.fromDimensions[0] nPoints];
	// first buffer will be used to store the transpose (which will be overwritten)
	GLBuffer *buffer1 = [[GLBuffer alloc] initWithLength: A.nDataElements*sizeof(GLFloat)];
	// second buffer will store the annoyingly formatted output
	GLBuffer *buffer2 = [[GLBuffer alloc] initWithLength: A.nDataElements*sizeof(GLFloat)];
	// third buffer is the lapack work buffer
	GLBuffer *buffer3 = [[GLBuffer alloc] initWithLength: 8*N*sizeof(GLFloat)];
	NSArray *buffers = @[buffer1, buffer2, buffer3];
	
	variableOperation op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
		GLFloat *A = (GLFloat *) [operandArray[0] bytes];
		
		GLFloat *B = (GLFloat *) [bufferArray[0] bytes];
		GLFloat *output = (GLFloat *) [bufferArray[1] bytes];
		NSMutableData *work = bufferArray[2];
		
		GLSplitComplex v = splitComplexFromData(resultArray[0]);
		GLSplitComplex C = splitComplexFromData(resultArray[1]);
		
		// clapack takes matrices in column-major format.
		// To be clever, we could use the left-eigenvectors instead of the right-eigenvectors and just take the conjugate, but we need to prevent the input from being overwritten anyway.
		vGL_mtrans(A, 1, B, 1, N, N);
		
		char JOBVL ='N';
		char JOBVR ='V';
		__CLPK_integer n = (__CLPK_integer) N;
		__CLPK_integer lwork = 8*n;
		__CLPK_integer info;
		
		sgeev_(&JOBVL, &JOBVR, &n, B, &n, v.realp, v.imagp, NULL, &n, output, &n, work.mutableBytes, &lwork, (__CLPK_integer *)&info);
		
		if (info != 0) {
			printf("sgeev failed with error code %d\n", (int)info);
		}
		
		// Now we have to get the eigenvectors in the proper format.
		// If the j-th eigenvalue is real, then v(j) = VR(:,j), the j-th column of VR.
		// If the j-th and (j+1)-st eigenvalues form a complex conjugate pair,
		// then v(j) = VR(:,j) + i*VR(:,j+1) and v(j+1) = VR(:,j) - i*VR(:,j+1).
		// And, don't forget, we need to fix the transpose.
		
		vGL_vclr( C.imagp, 1,  N*N);
		
		NSUInteger j=0;
		for( NSUInteger i = 0; i < N; i++ ) { // i indicates which eigenvector we're copying
			if ( v.imagp[i] == (GLFloat)0.0 ) {
				for( NSUInteger k = 0; k < N; k++ ) { // k walks down the column
					C.realp[k*N+i] = output[j*N+k];
				}
				j++;
			} else {
				for( NSUInteger k = 0; k < N; k++ ) {
					C.realp[k*N+i] = output[j*N+k];
					C.imagp[k*N+i] = output[(j+1)*N+k];
					
					C.realp[k*N+(i+1)] = output[j*N+k];
					C.imagp[k*N+(i+1)] = -output[(j+1)*N+k];
				}
				j+=2;
				i++;
			}
		}
	};
	
	if (( self = [super initWithResult: results operand: @[A] buffers: buffers operation: op] )) {
        
    }
    return self;
}

@end

/************************************************/
/*		GLGeneralizedMatrixEigensystemOperation */
/************************************************/

#pragma mark -
#pragma mark GLGeneralizedMatrixEigensystemOperation
#pragma mark

@implementation GLGeneralizedMatrixEigensystemOperation

- (id) initWithFirstOperand: (GLLinearTransform *) A secondOperand: (GLLinearTransform *) B
{
    if (A.matrixDescription.nDimensions != 1 || B.matrixDescription.nDimensions != 1) {
        [NSException raise: @"MatrixWrongFormat" format: @"We can only do one dimensional matrices at the moment."];
    }
    
    if (A.fromDimensions.count != A.toDimensions.count || B.fromDimensions.count != B.toDimensions.count ) {
        [NSException raise: @"MatrixWrongFormat" format: @"We can only do square matrices."];
    }
	
	for (NSUInteger i=0; i<A.fromDimensions.count; i++) {
		if ( ![A.fromDimensions[i] isEqualToDimension: A.toDimensions[i]] ) {
			[NSException raise: @"MatrixWrongFormat" format: @"By assumption, a linear transformation must be an endomorphism to compute the eigensystem."];
		}
	}
	
	// We need to construct a *new* eigenbasis.
	// I'm not quite sure the right definitions to use.
	NSMutableArray *eigenbasis = [NSMutableArray array];
	for (GLDimension *dim in A.fromDimensions) {
		GLDimension *newDim = [[GLDimension alloc] initDimensionWithGrid: dim.gridType nPoints: dim.nPoints domainMin:dim.domainMin length: dim.domainLength];
		if (dim.name && ![dim.name isEqualToString: @""]) {
			newDim.name = [NSString stringWithFormat: @"%@_eigen", dim.name];
		} else {
			newDim.name = @"eigen";
		}
		[eigenbasis addObject: newDim];
	}
	
	// Should check whether or not the matrix is symmetric.
	GLVariable *eigenvalues = [GLFunction functionOfComplexTypeWithDimensions: eigenbasis forEquation: A.equation];
	GLLinearTransform *eigentransform = [GLLinearTransform transformOfType: kGLSplitComplexDataFormat withFromDimensions: eigenbasis toDimensions:A.toDimensions inFormat:A.matrixFormats forEquation:A.equation matrix: nil];
    NSArray *results = @[eigenvalues, eigentransform];
    
    // http://www.nag.com/lapack-ex/node122.html
    
	NSUInteger N = [A.fromDimensions[0] nPoints];
	// first two buffers will be used to store the transpose (which will be overwritten)
	GLBuffer *buffer1 = [[GLBuffer alloc] initWithLength: A.nDataElements*sizeof(GLFloat)];
    GLBuffer *buffer2 = [[GLBuffer alloc] initWithLength: B.nDataElements*sizeof(GLFloat)];
    // third buffer will store the 'beta' part of the eigenvalues.
    GLBuffer *buffer3 = [[GLBuffer alloc] initWithLength: N*sizeof(GLFloat)];
	// fourth buffer will store the annoyingly formatted output
	GLBuffer *buffer4 = [[GLBuffer alloc] initWithLength: A.nDataElements*sizeof(GLFloat)];
	// fifth buffer is the lapack work buffer
	GLBuffer *buffer5 = [[GLBuffer alloc] initWithLength: 8*N*sizeof(GLFloat)];
	NSArray *buffers = @[buffer1, buffer2, buffer3, buffer4, buffer5];
	
	variableOperation op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
		GLFloat *A_row = (GLFloat *) [operandArray[0] bytes];
        GLFloat *B_row = (GLFloat *) [operandArray[1] bytes];
		
		GLFloat *A_col = (GLFloat *) [bufferArray[0] bytes];
        GLFloat *B_col = (GLFloat *) [bufferArray[1] bytes];
        
        GLFloat *beta = (GLFloat *) [bufferArray[2] bytes];
		GLFloat *output = (GLFloat *) [bufferArray[3] bytes];
		NSMutableData *work = bufferArray[4];
		
		GLSplitComplex v = splitComplexFromData(resultArray[0]);
		GLSplitComplex C = splitComplexFromData(resultArray[1]);
		
		// clapack takes matrices in column-major format.
		vGL_mtrans(A_row, 1, A_col, 1, N, N);
        vGL_mtrans(B_row, 1, B_col, 1, N, N);
		
		char JOBVL ='N';
		char JOBVR ='V';
		__CLPK_integer n = (__CLPK_integer) N;
		__CLPK_integer lwork = 8*n;
		__CLPK_integer info;
		
        sggev_(&JOBVL, &JOBVR, &n, A_col, &n, B_col, &n, v.realp, v.imagp, beta, NULL, &n, output, &n, work.mutableBytes, &lwork, &info);
		
		if (info != 0) {
			printf("sggev failed with error code %d\n", (int)info);
		}
		
		// Now we have to get the eigenvectors in the proper format.
		// If the j-th eigenvalue is real, then v(j) = VR(:,j), the j-th column of VR.
		// If the j-th and (j+1)-st eigenvalues form a complex conjugate pair,
		// then v(j) = VR(:,j) + i*VR(:,j+1) and v(j+1) = VR(:,j) - i*VR(:,j+1).
		// And, don't forget, we need to fix the transpose.
		
		vGL_vclr( C.imagp, 1,  N*N);
		
		NSUInteger j=0;
		for( NSUInteger i = 0; i < N; i++ ) { // i indicates which eigenvector we're copying
			if ( v.imagp[i] == (GLFloat)0.0 ) {
				for( NSUInteger k = 0; k < N; k++ ) { // k walks down the column
					C.realp[k*N+i] = output[j*N+k];
				}
				j++;
			} else {
				for( NSUInteger k = 0; k < N; k++ ) {
					C.realp[k*N+i] = output[j*N+k];
					C.imagp[k*N+i] = output[(j+1)*N+k];
					
					C.realp[k*N+(i+1)] = output[j*N+k];
					C.imagp[k*N+(i+1)] = -output[(j+1)*N+k];
				}
				j+=2;
				i++;
			}
		}
        
        vGL_vdiv(beta, 1, v.realp, 1, v.realp, 1, N);
        vGL_vdiv(beta, 1, v.imagp, 1, v.imagp, 1, N);
	};
	
	if (( self = [super initWithResult: results operand: @[A, B] buffers: buffers operation: op] )) {
        
    }
    return self;
}
@end

/************************************************/
/*		GLTensorProductOperation				*/
/************************************************/
#pragma mark -
#pragma mark GLTensorProductOperation
#pragma mark

@implementation GLTensorProductOperation

- (id) initWithLinearTransformations: (NSArray *) linearTransformations
{
	NSMutableArray *fromDimensions = [NSMutableArray array];
    NSMutableArray *toDimensions = [NSMutableArray array];
    NSMutableArray *matrixFormat = [NSMutableArray array];
    NSMutableArray *matrixBlocks = [NSMutableArray array];
    BOOL isComplex = NO;
	GLEquation *equation = [linearTransformations[0] equation];
    for (GLLinearTransform *transform in linearTransformations) {
        [fromDimensions addObject: transform.fromDimensions[0]];
		[toDimensions addObject: transform.toDimensions[0]];
        isComplex |= transform.isComplex;
        [matrixFormat addObject: transform.matrixFormats[0]];
		if (transform.matrixBlock) {
			[matrixBlocks addObject: transform.matrixBlock];
		}
    }
    GLDataFormat format = isComplex ? kGLSplitComplexDataFormat : kGLRealDataFormat;
    
	if (matrixBlocks.count == linearTransformations.count)
	{	// In this scenario, the linear transforms don't depend on any operations (they already have matrixBlocks), so we can immediately compute the tensor product.
		GLLinearTransform *tensorProduct = [GLLinearTransform transformOfType:format withFromDimensions: fromDimensions toDimensions: toDimensions inFormat:matrixFormat forEquation:equation matrix: ^( NSUInteger *row, NSUInteger *col ) {
			transformMatrix theMatrixBlock = matrixBlocks[0];
			GLFloatComplex value = theMatrixBlock(&(row[0]), &(col[0]));
			for (NSUInteger i=1;i<matrixBlocks.count;i++) {
				theMatrixBlock = matrixBlocks[i];
				value *= theMatrixBlock(&(row[i]), &(col[i]));
			}
			return value;
		}];
		
		if (( self = [super initWithResult: @[tensorProduct] operand:@[]] )) {
			
		}
		return self;
	} else {
		GLLinearTransform *tensorProduct = [GLLinearTransform transformOfType:format withFromDimensions: fromDimensions toDimensions: toDimensions inFormat:matrixFormat forEquation:equation matrix: NULL];
		
		variableOperation op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
		
		};
		
		if (( self = [super initWithResult: @[tensorProduct] operand: linearTransformations buffers: @[] operation: op] )) {
			
		}
		return self;
	}
}

@end
