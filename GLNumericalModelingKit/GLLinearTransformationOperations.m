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

NSUInteger compute_total_loops( GLMatrixDescription *matrixDescription, GLMatrixDescription *vectorDescription, NSUInteger loopIndex )
{
    NSUInteger lastNonTrivialNonLoopIndex = NSNotFound;
	NSUInteger lastTrivialIndex = NSNotFound;
	NSUInteger totalTrivialPoints = 1;
	
    for ( NSUInteger index=0; index < matrixDescription.nDimensions; index++) {
		if (matrixDescription.strides[index].matrixFormat == kGLIdentityMatrixFormat) {
            lastTrivialIndex = index;
			totalTrivialPoints *= matrixDescription.strides[index].nPoints;
        } else if (matrixDescription.strides[index].matrixFormat == kGLDiagonalMatrixFormat && index != loopIndex ) {
			lastNonTrivialNonLoopIndex = index;
		} else if (index != loopIndex) {
			[NSException raise: @"BadMatrixFormat" format:@"Cannot apply the matrix loop across a matrix of this format"];
		}
    }
	NSUInteger totalEquations = matrixDescription.nPoints / matrixDescription.strides[loopIndex].nPoints;
	NSUInteger totalOuterLoops = totalTrivialPoints;
    
    return totalEquations*totalOuterLoops;
}

void apply_matrix_loop( GLMatrixDescription *matrixDescription, GLMatrixDescription *vectorDescription, NSUInteger loopIndex, dispatch_queue_t queue, void (^block)(NSUInteger, NSUInteger, NSUInteger ))
{
	NSUInteger lastNonTrivialNonLoopIndex = NSNotFound;
	NSUInteger lastTrivialIndex = NSNotFound;
	NSUInteger totalTrivialPoints = 1;
	
    for ( NSUInteger index=0; index < matrixDescription.nDimensions; index++) {
		if (matrixDescription.strides[index].matrixFormat == kGLIdentityMatrixFormat) {
            lastTrivialIndex = index;
			totalTrivialPoints *= matrixDescription.strides[index].nPoints;
        } else if (matrixDescription.strides[index].matrixFormat == kGLDiagonalMatrixFormat && index != loopIndex ) {
			lastNonTrivialNonLoopIndex = index;
		} else if (index != loopIndex) {
			[NSException raise: @"BadMatrixFormat" format:@"Cannot apply the matrix loop across a matrix of this format"];
		}
    }
	
	// How many 'inner loop' steps we need to take depends on how many other nontrivial dimensions there are.
	// Total equations should be the product of diagonal dimensions that are not the loop index.
	NSUInteger inEquationStride = lastNonTrivialNonLoopIndex == NSNotFound ? 0 : matrixDescription.strides[lastNonTrivialNonLoopIndex].stride;
	NSUInteger totalEquations = matrixDescription.nPoints / matrixDescription.strides[loopIndex].nPoints;
	
	// Now we need the strides to match up to the inner loop
	NSUInteger outEquationStride = lastNonTrivialNonLoopIndex == NSNotFound ? 0 : vectorDescription.strides[lastNonTrivialNonLoopIndex].stride;
	
	// Finally, we need the outer loop strides and totals
	NSUInteger totalOuterLoops = totalTrivialPoints;
	NSUInteger outerLoopStride = lastTrivialIndex == NSNotFound ? 0 : vectorDescription.strides[lastTrivialIndex].stride;
	
	if (totalOuterLoops > 1 && totalEquations > 1 ) {
		// This operation has two loops.
		// The inner loop walks over non-trivial elements of the linear operator, while
		// the out loop walks over trivial elements of the linear operator, really just to new x and b elements
		dispatch_apply(totalOuterLoops, queue, ^(size_t outerIteration) {
			
			NSInteger outerOffset = outerIteration*outerLoopStride;
			
			dispatch_apply(totalEquations, queue, ^(size_t iteration) {
				
				NSUInteger inEquationPos = iteration*inEquationStride;
				NSUInteger outEquationPos = outerOffset + iteration*outEquationStride;
				
				block(outerIteration*iteration, inEquationPos, outEquationPos);
			});
		});
	} else if (totalOuterLoops == 1 && totalEquations > 1 ) {
		dispatch_apply(totalEquations, queue, ^(size_t iteration) {
			
			NSUInteger inEquationPos = iteration*inEquationStride;
			NSUInteger outEquationPos = iteration*outEquationStride;
			
			block(iteration, inEquationPos, outEquationPos);
		});
	} else if (totalOuterLoops > 1 && totalEquations == 1 ) {
		dispatch_apply(totalOuterLoops, queue, ^(size_t outerIteration) {
			
			NSInteger outerOffset = outerIteration*outerLoopStride;
			
			NSUInteger inEquationPos = 0;
			NSUInteger outEquationPos = outerOffset;
			
			block(outerIteration, inEquationPos, outEquationPos);
		});
	} else if (totalOuterLoops == 1 && totalEquations == 1 ) {
		block(0, 0, 0);
	}
}

#pragma mark -
#pragma mark Creation
#pragma mark

/************************************************/
/*		GLDiagonalTransformCreationOperation	*/
/************************************************/

/** Places the values of the function along the diagonal of matrix
 @param function The function to be placed along the diagonal.
 @returns The transformed function with fromDimensions and toDimensions matching the dimension of the function.
 */
@implementation GLDiagonalTransformCreationOperation
- (id) initWithFunction: (GLFunction *) function
{
    NSMutableArray *matrixFormats = [NSMutableArray array];
    for (GLDimension *aDim in function.dimensions) {
        [matrixFormats addObject: @(kGLDiagonalMatrixFormat)];
    }
	GLLinearTransform *result = [[GLLinearTransform alloc] initTransformOfType: function.dataFormat withFromDimensions: function.dimensions toDimensions:function.dimensions inFormat:matrixFormats forEquation:function.equation matrix:nil];
	
    if (function.name) {
        result.name = [NSString stringWithFormat: @"%@_transform", function.name];
    }
    
    NSUInteger nBytes = function.matrixDescription.nBytes;
    if ((self = [super initWithResult:@[result] operand:@[function] buffers:nil operation:^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
        NSMutableData *result = resultArray[0];
        NSMutableData *transform = operandArray[0];
        memcpy(result.mutableBytes, transform.mutableBytes, nBytes);
    }])) {
        
    }
    
    return self;
}

@end

/************************************************/
/*		GLExpandMatrixDimensionsOperation		*/
/************************************************/

@implementation GLExpandMatrixDimensionsOperation
- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform fromDimensions: (NSArray *) fromDims toDimensions: (NSArray *) toDims
{
	if (fromDims.count != toDims.count || fromDims.count < linearTransform.fromDimensions.count) {
		[NSException raise: @"BadFormat" format: @"Dimensions don't make sense."];
	}
	
	NSMutableArray *matrixFormats = [NSMutableArray array];
	NSUInteger iDim = 0;
	for (NSUInteger i=0; i<fromDims.count; i++) {
		GLDimension *fromDim = fromDims[i];
		GLDimension *toDim = toDims[i];
		if (iDim < linearTransform.fromDimensions.count && [linearTransform.fromDimensions[iDim] isEqualToDimension: fromDim] && [linearTransform.toDimensions[iDim] isEqualToDimension: toDim]) {
			[matrixFormats addObject: linearTransform.matrixFormats[iDim]];
			iDim++;
		} else {
			[matrixFormats addObject: @(kGLIdentityMatrixFormat)];
		}
	}
	
	if (iDim != linearTransform.fromDimensions.count) {
		[NSException raise: @"BadFormat" format: @"Dimensions don't make sense."];
	}
	
	GLLinearTransform *result = [[GLLinearTransform alloc] initTransformOfType: linearTransform.dataFormat withFromDimensions: fromDims toDimensions:toDims inFormat:matrixFormats forEquation:linearTransform.equation matrix:nil];
    NSUInteger nBytes = result.matrixDescription.nBytes;
	
	if (nBytes != linearTransform.matrixDescription.nBytes) {
		[NSException raise: @"BadFormat" format: @"Urg. By assumption these shoudl be the same."];
	}
	
    if ((self = [super initWithResult:@[result] operand:@[linearTransform] buffers:nil operation:^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
        NSMutableData *result = resultArray[0];
        NSMutableData *transform = operandArray[0];
        memcpy(result.mutableBytes, transform.mutableBytes, nBytes);
    }])) {
        
    }
    
    return self;
}

@end

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
    
	GLMatrixDescription *operandDescription = linearTransform.matrixDescription;

    NSUInteger triIndex = NSNotFound;
	NSUInteger numTriIndices = 0;
	for ( NSUInteger index=0; index < operandDescription.nDimensions; index++) {
        if (operandDescription.strides[index].matrixFormat == kGLTridiagonalMatrixFormat ) {
			numTriIndices++;
			triIndex = index;
		}
    }
	
    if ( numTriIndices != 1 ) {
        [NSException raise: @"TridiagonalIndexNotFound" format: @"Unable to find a tridiagonal index."];
    }
    
	BOOL isComplex = linearTransform.isComplex || function.isComplex;
	GLDataFormat format = isComplex ? kGLSplitComplexDataFormat : kGLRealDataFormat;
	GLFunction *result = [GLFunction functionOfType: format withDimensions: linearTransform.toDimensions forEquation: linearTransform.equation];
	GLMatrixDescription *resultDescription = result.matrixDescription;
	
    if (linearTransform.name && function.name) {
        result.name = [NSString stringWithFormat: @"%@_%@", function.name, linearTransform.name];
    }
    
	if (( self = [super initWithResult: @[result] operand: @[linearTransform, function]] )) {
        GLMatrixDescription *matrixDescription = linearTransform.matrixDescription;
        
        dispatch_queue_t globalQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
        
        NSUInteger inNDiagonalPoints = matrixDescription.strides[triIndex].nDiagonalPoints;
        NSUInteger inElementStride = matrixDescription.strides[triIndex].stride;
        NSUInteger inDiagonalStride = matrixDescription.strides[triIndex].diagonalStride;
        NSUInteger outElementStride = result.matrixDescription.strides[triIndex].stride;

        self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				
			apply_matrix_loop(operandDescription, resultDescription, triIndex, globalQueue, ^(NSUInteger iteration, NSUInteger inEquationPos, NSUInteger outEquationPos) {
				
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
				vGL_mmul( A, 1, B, 1, C, 1, M, N, K);
			};
		}
		else if ( linearTransform.isComplex && !function.isComplex)
		{	// (A+iB).(X) = A.X + iB.X
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLSplitComplex A = splitComplexFromData( operandArray[0] );
				GLFloat *B = (GLFloat *) [operandArray[1] bytes];
				GLSplitComplex C = splitComplexFromData( resultArray[0] );
				
				vGL_mmul( A.realp, 1, B, 1, C.realp, 1, M, N, K);
				vGL_mmul( A.imagp, 1, B, 1, C.imagp, 1, M, N, K);
			};
		}
		else if ( !linearTransform.isComplex && function.isComplex)
		{	// A.(X+iY) = A.X + iA.Y
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *A = (GLFloat *) [operandArray[0] bytes];
				GLSplitComplex B = splitComplexFromData( operandArray[1] );
				GLSplitComplex C = splitComplexFromData( resultArray[0] );
				
				vGL_mmul( A, 1, B.realp, 1, C.realp, 1, M, N, K);
				vGL_mmul( A, 1, B.imagp, 1, C.imagp, 1, M, N, K);
			};
		}
		else if ( linearTransform.isComplex && function.isComplex)
		{
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLSplitComplex A = splitComplexFromData( operandArray[0] );
				GLSplitComplex B = splitComplexFromData( operandArray[1] );
				GLSplitComplex C = splitComplexFromData( resultArray[0] );
				
				vGL_zmmul( &A, 1, &B, 1, &C, 1, M, N, K);
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
    
	GLMatrixDescription *operandDescription = linearTransform.matrixDescription;
	
    NSUInteger triIndex = NSNotFound;
	NSUInteger numTriIndices = 0;
	for ( NSUInteger index=0; index < operandDescription.nDimensions; index++) {
        if (operandDescription.strides[index].matrixFormat == kGLTridiagonalMatrixFormat ) {
			numTriIndices++;
			triIndex = index;
		}
    }
	
    if ( numTriIndices != 1 ) {
        [NSException raise: @"TridiagonalIndexNotFound" format: @"Unable to find a tridiagonal index."];
    }
	
	BOOL isComplex = linearTransform.isComplex || function.isComplex;
	GLDataFormat format = isComplex ? kGLSplitComplexDataFormat : kGLRealDataFormat;
	GLFunction *result = [GLFunction functionOfType: format withDimensions: linearTransform.fromDimensions forEquation: linearTransform.equation];
    GLMatrixDescription *resultDescription = result.matrixDescription;
	
	// This buffer is 3 times as large as it needs to be, but it makes bookkeeping easier.
	GLBuffer *buffer = [[GLBuffer alloc] initWithLength: result.nDataPoints*sizeof(GLFloat)];
	
	if (( self = [super initWithResult: @[result] operand: @[linearTransform, function] buffers: @[buffer] operation: nil] )) {
		
        dispatch_queue_t globalQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
        
        NSUInteger inNDiagonalPoints = operandDescription.strides[triIndex].nDiagonalPoints;
        NSUInteger inElementStride = operandDescription.strides[triIndex].stride;
        NSUInteger inDiagonalStride = operandDescription.strides[triIndex].diagonalStride;
        NSUInteger outElementStride = resultDescription.strides[triIndex].stride;
        
        self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
            
            apply_matrix_loop(operandDescription, resultDescription, triIndex, globalQueue, ^(NSUInteger iteration, NSUInteger inEquationPos, NSUInteger outEquationPos) {
                
                GLFloat *f = (GLFloat *) [operandArray[0] bytes];
                
                // plus the offset!!!!
                GLFloat *a = &(f[0*inDiagonalStride]);
                GLFloat *b = &(f[1*inDiagonalStride]);
                GLFloat *c = &(f[2*inDiagonalStride]);
                
                GLFloat *cprime = (GLFloat *) [bufferArray[0] mutableBytes];
                
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
    
	GLMatrixDescription *matrixDescription = linearTransform.matrixDescription;
	
    NSUInteger denseIndex = NSNotFound;
	NSUInteger numDenseIndices = 0;
	for ( NSUInteger index=0; index < matrixDescription.nDimensions; index++) {
        if (matrixDescription.strides[index].matrixFormat == kGLTridiagonalMatrixFormat ) {
			numDenseIndices++;
			denseIndex = index;
		}
    }
	
    if ( numDenseIndices != 1 ) {
        [NSException raise: @"BadInputFormat" format: @"The GLDenseMatrixSolver can only take exactly one dense dimension."];
    }
	
    
	BOOL isComplex = linearTransform.isComplex || function.isComplex;
	GLDataFormat format = isComplex ? kGLSplitComplexDataFormat : kGLRealDataFormat;
	GLFunction *result = [GLFunction functionOfType: format withDimensions: linearTransform.fromDimensions forEquation: linearTransform.equation];
	GLMatrixDescription *vectorDescription = result.matrixDescription;
    
    dispatch_queue_t globalQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
    
    NSUInteger N = matrixDescription.strides[denseIndex].nRows;
    NSUInteger inElementStride = matrixDescription.strides[denseIndex].stride;
    NSUInteger outElementStride = vectorDescription.strides[denseIndex].stride;
    
    // sgesv is only capable of handling input arrays with a stride of 1. For us, this means that if that lastNonTrivialNonDenseIndex is greater
    // than the denseIndex, we have to copy the elements to a temporary, contiguous buffer.
    if ( inElementStride != 1 ) {
        NSLog(@"Warning! The dense dimension is not the last non-trivial dimensions, which results in slower solutions times due to an extra memory copy.");
    }
    // sgesv also overwrites the b vector (in Mx=b) with the output x. So we need to copy the output as well.
    
    NSUInteger totalLoops = compute_total_loops(matrixDescription, vectorDescription, denseIndex);
    GLBuffer *buffer1 = [[GLBuffer alloc] initWithLength: totalLoops*N*sizeof(__CLPK_integer)];
    NSMutableArray *buffers = [NSMutableArray arrayWithObject: buffer1];
    if (inElementStride != 1) {
        GLBuffer *buffer2 = [[GLBuffer alloc] initWithLength: matrixDescription.nBytes];
        [buffers addObject: buffer2];
    }
    if (outElementStride != 1) {
        GLBuffer *buffer3 = [[GLBuffer alloc] initWithLength: vectorDescription.nBytes];
        [buffers addObject: buffer3];
    }
    
    variableOperation op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
        
        apply_matrix_loop(matrixDescription, vectorDescription, denseIndex, globalQueue, ^(NSUInteger iteration, NSUInteger inEquationPos, NSUInteger outEquationPos) {
            __CLPK_integer *ipivData = (__CLPK_integer *) [bufferArray[0] bytes];
            __CLPK_integer *ipiv = &(ipivData[iteration*N]);
            __CLPK_integer n = (__CLPK_integer) N;
            __CLPK_integer nrhs = 1;
            __CLPK_integer info;
            
            GLFloat *MData = (GLFloat *) [operandArray[0] bytes];
            GLFloat *bData = (GLFloat *) [operandArray[1] bytes];
            GLFloat *xData = (GLFloat *) [resultArray[0] bytes];
            
            GLFloat *M;
            GLFloat *b;
            GLFloat *x;
            
            if (inElementStride != 1) {
                // Take the strided data, copy it to the buffer unstrided.
                GLFloat *Mstrided = &(MData[inEquationPos]);
                GLFloat *MbufferData = (GLFloat *) [bufferArray[1] bytes];
                M = &(MbufferData[N*N*iteration]);
                for (NSUInteger i=0; i < matrixDescription.strides[denseIndex].nRows; i++) {
                    for (NSUInteger j=0; j < matrixDescription.strides[denseIndex].nColumns; j++) {
                        M[i+N*j] = Mstrided[i*matrixDescription.strides[denseIndex].rowStride + j*matrixDescription.strides[denseIndex].columnStride];
                    }
                }
            } else {
                M = &(MData[inEquationPos]);
            }
            
            if ( outElementStride != 1) {
                GLFloat *bstrided = &(bData[outEquationPos]);
                GLFloat *xbufferData = (GLFloat *) [[bufferArray lastObject] bytes];
                x = &(xbufferData[N*iteration]);
                for (NSUInteger i=0; i < vectorDescription.strides[denseIndex].nRows; i++) {
                    x[i] = bstrided[i*vectorDescription.strides[denseIndex].stride];
                }
            } else {
                b = &(bData[outEquationPos]);
                x = &(xData[outEquationPos]);
                memcpy( x, b, n*n*sizeof(GLFloat));
            }
            
            sgesv_( &n, &nrhs, M, &n, ipiv, x, &n, (__CLPK_integer *) &info );
            
            if (info != 0) {
                printf("sgesv failed with error code %d\n", (int)info);
            }
            
            if ( outElementStride != 1) {
                GLFloat *output = &(xData[outEquationPos]);
                for (NSUInteger i=0; i < vectorDescription.strides[denseIndex].nRows; i++) {
                    output[i*vectorDescription.strides[denseIndex].stride] = x[i];
                }
            }
            
        });
    };
    
	if (( self = [super initWithResult: @[result] operand: @[linearTransform, function] buffers: buffers operation: op] )) {
    }
    return self;
}

@end

/********************************************************/
/*		GLDiagonalMatrixMatrixMultiplicationOperation   */
/********************************************************/

#pragma mark -
#pragma mark GLDiagonalMatrixMatrixMultiplicationOperation
#pragma mark

@implementation GLDiagonalMatrixMatrixMultiplicationOperation

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
	
	A = [A copyWithDataType: format matrixFormat: A.matrixFormats ordering: kGLRowMatrixOrder];
	B = [B copyWithDataType: format matrixFormat: B.matrixFormats ordering: kGLRowMatrixOrder];
    
	GLLinearTransform *result = [GLLinearTransform transformOfType: format withFromDimensions: B.fromDimensions toDimensions:A.toDimensions inFormat:B.matrixFormats forEquation:B.equation matrix: nil];
    
	if (( self = [super initWithResult: @[result] operand: @[A, B]] )) {
                
        int M = (int) A.matrixDescription.strides[0].nRows;
        int N = (int) B.matrixDescription.strides[0].nColumns;
		int K = (int) A.matrixDescription.strides[0].nColumns;
		
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
	
	A = [A copyWithDataType: format matrixFormat: A.matrixFormats ordering: kGLRowMatrixOrder];
	B = [B copyWithDataType: format matrixFormat: B.matrixFormats ordering: kGLRowMatrixOrder];

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
	
	if (linearTransform.matrixDescription.strides[0].matrixFormat != kGLDenseMatrixFormat) {
        [NSException raise: @"MatrixWrongFormat" format: @"Can only inverse dense matrices."];
    }
	
	if (linearTransform.matrixDescription.strides[0].nRows != linearTransform.matrixDescription.strides[0].nColumns) {
        [NSException raise: @"MatrixWrongFormat" format: @"We can only do square matrices."];
    }
    
    GLLinearTransform *A = (GLLinearTransform *) linearTransform;
	
	GLDataFormat format = A.isComplex ? kGLSplitComplexDataFormat : kGLRealDataFormat;
	GLLinearTransform *result;
	NSArray *buffers;
	variableOperation op;
	
	if (A.isComplex) {
		A = [A copyWithDataType: kGLInterleavedComplexDataFormat matrixFormat: A.matrixFormats ordering: kGLColumnMatrixOrder];
		result = [[GLLinearTransform alloc] initTransformOfType: kGLInterleavedComplexDataFormat withFromDimensions: A.toDimensions toDimensions: A.fromDimensions inFormat: A.matrixFormats withOrdering: kGLColumnMatrixOrder forEquation: A.equation matrix: nil];
		
		NSUInteger N = [A.fromDimensions[0] nPoints];
		GLBuffer *buffer1 = [[GLBuffer alloc] initWithLength: N*sizeof(__CLPK_integer)];
		GLBuffer *buffer2 = [[GLBuffer alloc] initWithLength: N*N*sizeof(__CLPK_complex)];
		buffers = @[buffer1, buffer2];
		
		op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			__CLPK_complex *A = (__CLPK_complex *) [operandArray[0] bytes];
			__CLPK_complex *C = (__CLPK_complex *) [resultArray[0] bytes];
			NSMutableData *ipiv = bufferArray[0];
			NSMutableData *work = bufferArray[1];
			
			// clapack takes matrices in column-major format.
			// However, the transpose of the inverse is the inverse of the transpose---so we don't need to worry here.
			
			__CLPK_integer n = (__CLPK_integer) N;
			__CLPK_integer info;
			memcpy( C, A, n*n*sizeof(__CLPK_complex));
			cgetrf_(&n, &n, C, &n, ipiv.mutableBytes, (__CLPK_integer *)&info);
			
			if (info != 0) {
				printf("cgetrf_ failed with error code %d\n", (int)info);
			}
			
			__CLPK_integer lwork = n*n;
			cgetri_(&n, C, &n, ipiv.mutableBytes, work.mutableBytes, &lwork, (__CLPK_integer *)&info);
			
			if (info != 0) {
				printf("cgetri_ failed with error code %d\n", (int)info);
			}
		};
		
	} else {
		result = [GLLinearTransform transformOfType: format withFromDimensions: A.toDimensions toDimensions:A.fromDimensions inFormat:A.matrixFormats forEquation:A.equation matrix: nil];
		NSUInteger N = [A.fromDimensions[0] nPoints];
		GLBuffer *buffer1 = [[GLBuffer alloc] initWithLength: N*sizeof(__CLPK_integer)];
		GLBuffer *buffer2 = [[GLBuffer alloc] initWithLength: N*N*sizeof(GLFloat)];
		buffers = @[buffer1, buffer2];
		
		op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
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
	}

	if (( self = [super initWithResult: @[result] operand: @[A] buffers: buffers operation: op] )) {
        
        
    }
	
    return self;
}

@end

/************************************************/
/*		GLMatrixNormalizationOperation          */
/************************************************/

#pragma mark -
#pragma mark GLMatrixNormalizationOperation
#pragma mark

@implementation GLMatrixNormalizationOperation

// Int const*f*f = 1
- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform normalizationConstant: (GLFloat) aConst dimensionIndex: (NSUInteger) index;
{
	GLLinearTransform *A = (GLLinearTransform *) linearTransform;
	if (index >= A.fromDimensions.count) {
		[NSException raise: @"InvalidIndex" format: @"The index you provided is greater than the number of dimensions."];
	}
	
	if (aConst <= 0.0) {
		[NSException raise: @"InvalidNormalization" format: @"The normalization must be nonzero."];
	}
	
	if ( ![A.toDimensions[index] isEvenlySampled] ) {
		[NSException raise: @"NotYetImplemented" format: @"Normalization is only possible on evenly sampled dimensions."];
	}
	
	GLLinearTransform *result = [[GLLinearTransform alloc] initTransformOfType: A.dataFormat withFromDimensions: A.toDimensions toDimensions: A.fromDimensions inFormat: A.matrixFormats withOrdering: A.matrixOrder forEquation: A.equation matrix: nil];
	GLMatrixDescription *operandDescription = linearTransform.matrixDescription;
	
	NSUInteger totalVectors = operandDescription.nPoints / operandDescription.strides[index].nColumns;
	NSUInteger vectorStride = operandDescription.strides[index].columnStride;
	NSUInteger vectorLength = operandDescription.strides[index].nRows;
	NSUInteger vectorElementStride = operandDescription.strides[index].rowStride;
	NSUInteger complexStride = operandDescription.strides[index].complexStride;
	
	GLFloat deltaX = [A.toDimensions[index] sampleInterval];
	
	dispatch_queue_t queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
	variableOperation op;
	if (A.isComplex) {
		op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			dispatch_apply(totalVectors, queue, ^(size_t iteration) {
				
				NSUInteger inEquationPos = iteration*vectorStride;
				
				GLFloat *A = (GLFloat *) [operandArray[0] bytes];
				GLSplitComplex Asplit;
				Asplit.realp = &(A[inEquationPos]);
				Asplit.imagp = &(A[inEquationPos+complexStride]);
				
				GLFloat *B = (GLFloat *) [resultArray[0] bytes];
				GLSplitComplex Bsplit;
				Bsplit.realp = &(B[inEquationPos]);
				Bsplit.imagp = &(B[inEquationPos+complexStride]);
				
				// square it
				vGL_zvmags(&Asplit, vectorElementStride, Bsplit.realp, vectorElementStride, vectorLength);
								
				// sum it. sum = 1/(2*h) * (f_0 + 2*f_1 + 2*f_2 + ... + 2*f_{N-1} + 2*f_N
				GLFloat sum = 0.0;
				vGL_sve(Bsplit.realp, vectorElementStride, &sum, vectorLength);
				sum = aConst * deltaX * (sum - Bsplit.realp[0]/2.0 - Bsplit.realp[0+vectorElementStride*(vectorLength-1)]/2.0);
				
				GLFloat norm = 1/sqrt(fabs(sum));
				vGL_vsmul(Asplit.realp, vectorElementStride, &norm, Bsplit.realp, vectorElementStride, vectorLength);
				vGL_vsmul(Asplit.imagp, vectorElementStride, &norm, Bsplit.imagp, vectorElementStride, vectorLength);
			});
		};
	} else {
		op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			dispatch_apply(totalVectors, queue, ^(size_t iteration) {
				GLFloat *A = (GLFloat *) [operandArray[0] bytes];
				GLFloat *B = (GLFloat *) [resultArray[0] bytes];
				
				NSUInteger inEquationPos = iteration*vectorStride;

				// square it
				vGL_vsq(&(A[inEquationPos]), vectorElementStride, &(B[inEquationPos]), vectorElementStride, vectorLength);
				
				// sum it. sum = 1/(2*h) * (f_0 + 2*f_1 + 2*f_2 + ... + 2*f_{N-1} + 2*f_N
				GLFloat sum = 0.0;
				vGL_sve(&(B[inEquationPos]), vectorElementStride, &sum, vectorLength);
				sum = aConst * deltaX * (sum - B[inEquationPos]/2.0 - B[inEquationPos+vectorElementStride*(vectorLength-1)]/2.0);
				
				GLFloat norm = 1/sqrt(fabs(sum));
				vGL_vsmul(&(A[inEquationPos]), vectorElementStride, &norm, &(B[inEquationPos]), vectorElementStride, vectorLength);
			});
		};
	}
	
	
	if (( self = [super initWithResult: @[result] operand: @[A] buffers: @[] operation: op] )) {
        
    }
    return self;
}

// Int aFunction*f*f = 1
- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform normalizationFunction: (GLFunction *) aFunction
{
	GLLinearTransform *A = (GLLinearTransform *) linearTransform;
	NSUInteger index = [A.toDimensions indexOfObject: aFunction.dimensions[0]];
	if (index == NSNotFound) {
		[NSException raise: @"InvalidIndex" format: @"The index you provided is greater than the number of dimensions."];
	}
	
	if (aFunction.dimensions.count > 1) {
		[NSException raise: @"InvalidNormalizationFunction" format: @"The normalization function can only be one dimensional at this point."];
	}
	
	if ( ![A.toDimensions[index] isEvenlySampled] ) {
		[NSException raise: @"NotYetImplemented" format: @"Normalization is only possible on evenly sampled dimensions."];
	}
	
	if ( aFunction.isComplex ) {
		[NSException raise: @"NotYetImplemented" format: @"Normalization is only possible on real function."];
	}
	
	GLLinearTransform *result = [[GLLinearTransform alloc] initTransformOfType: A.dataFormat withFromDimensions: A.toDimensions toDimensions: A.fromDimensions inFormat: A.matrixFormats withOrdering: A.matrixOrder forEquation: A.equation matrix: nil];
	GLMatrixDescription *operandDescription = linearTransform.matrixDescription;
	
#warning The vector stride should be the columnStride for a 1D matrix, but could be 1 in a 3D matrix.
	
	NSUInteger totalVectors = operandDescription.nPoints / operandDescription.strides[index].nColumns;
	NSUInteger vectorStride = 1; //operandDescription.strides[index].columnStride;
	NSUInteger vectorLength = operandDescription.strides[index].nRows;
	NSUInteger vectorElementStride = operandDescription.strides[index].rowStride;
	NSUInteger complexStride = operandDescription.strides[index].complexStride;
	
	GLFloat deltaX = [A.toDimensions[index] sampleInterval];
	
	dispatch_queue_t queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
	variableOperation op;
	if (A.isComplex) {
		op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			dispatch_apply(totalVectors, queue, ^(size_t iteration) {
				
				NSUInteger inEquationPos = iteration*vectorStride;
				
				GLFloat *A = (GLFloat *) [operandArray[0] bytes];
				GLSplitComplex Asplit;
				Asplit.realp = &(A[inEquationPos]);
				Asplit.imagp = &(A[inEquationPos+complexStride]);
				
				GLFloat *f = (GLFloat *) [operandArray[1] bytes];
				
				GLFloat *B = (GLFloat *) [resultArray[0] bytes];
				GLSplitComplex Bsplit;
				Bsplit.realp = &(B[inEquationPos]);
				Bsplit.imagp = &(B[inEquationPos+complexStride]);
				
				// square it
				vGL_zvmags(&Asplit, vectorElementStride, Bsplit.realp, vectorElementStride, vectorLength);
				
				// Multiply by the normalization function.
				vGL_vmul( Bsplit.realp, vectorElementStride, f, 1, Bsplit.realp, vectorElementStride, vectorLength );
				
				// sum it. sum = 1/(2*h) * (f_0 + 2*f_1 + 2*f_2 + ... + 2*f_{N-1} + 2*f_N
				GLFloat sum = 0.0;
				vGL_sve(Bsplit.realp, vectorElementStride, &sum, vectorLength);
				sum = deltaX * (sum - Bsplit.realp[0]/2.0 - Bsplit.realp[0+vectorElementStride*(vectorLength-1)]/2.0);
				
				GLFloat norm = 1/sqrt(fabs(sum));
				vGL_vsmul(Asplit.realp, vectorElementStride, &norm, Bsplit.realp, vectorElementStride, vectorLength);
				vGL_vsmul(Asplit.imagp, vectorElementStride, &norm, Bsplit.imagp, vectorElementStride, vectorLength);
			});
		};
	} else {
		op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			dispatch_apply(totalVectors, queue, ^(size_t iteration) {
				GLFloat *A = (GLFloat *) [operandArray[0] bytes];
				GLFloat *f = (GLFloat *) [operandArray[1] bytes];
				GLFloat *B = (GLFloat *) [resultArray[0] bytes];
				
				NSUInteger inEquationPos = iteration*vectorStride;
				
				// square it
				vGL_vsq(&(A[inEquationPos]), vectorElementStride, &(B[inEquationPos]), vectorElementStride, vectorLength);
				
				// Multiply by the normalization function.
				vGL_vmul( &(B[inEquationPos]), vectorElementStride, f, 1, &(B[inEquationPos]), vectorElementStride, vectorLength );
				
				// sum it. sum = 1/(2*h) * (f_0 + 2*f_1 + 2*f_2 + ... + 2*f_{N-1} + 2*f_N
				GLFloat sum = 0.0;
				vGL_sve(&(B[inEquationPos]), vectorElementStride, &sum, vectorLength);
				sum = deltaX * (sum - B[inEquationPos]/2.0 - B[inEquationPos+vectorElementStride*(vectorLength-1)]/2.0);
				
				GLFloat norm = 1/sqrt(fabs(sum));
				vGL_vsmul(&(A[inEquationPos]), vectorElementStride, &norm, &(B[inEquationPos]), vectorElementStride, vectorLength);
			});
		};
	}
	
	
	if (( self = [super initWithResult: @[result] operand: @[A, aFunction] buffers: @[] operation: op] )) {
        
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
	GLLinearTransform *A = (GLLinearTransform *) linearTransform;
    for (NSUInteger i=0; i<A.fromDimensions.count; i++) {
		if ( ![A.fromDimensions[i] isEqualToDimension: A.toDimensions[i]] ) {
			[NSException raise: @"MatrixWrongFormat" format: @"By assumption, a linear transformation must be an endomorphism to compute the eigensystem."];
		}
	}
	
	GLMatrixDescription *operandDescription = linearTransform.matrixDescription;
	
    NSUInteger denseIndex = NSNotFound;
	NSUInteger numDenseIndices = 0;
	NSUInteger densifiableIndex = NSNotFound;
	NSUInteger numDensifiableIndices = 0;
	for ( NSUInteger index=0; index < operandDescription.nDimensions; index++) {
        if (operandDescription.strides[index].matrixFormat == kGLDenseMatrixFormat ) {
			numDenseIndices++;
			denseIndex = index;
			if (operandDescription.strides[index].nRows != operandDescription.strides[index].nColumns) {
				[NSException raise: @"MatrixWrongFormat" format: @"We can only do square matrices."];
			}
		} else if (operandDescription.strides[index].matrixFormat != kGLIdentityMatrixFormat && operandDescription.strides[index].matrixFormat != kGLDiagonalMatrixFormat) {
			numDensifiableIndices++;
			densifiableIndex = index;
		}
    }
	
    if ( numDenseIndices+numDensifiableIndices != 1 ) {
        [NSException raise: @"DenseIndexNotFound" format: @"Unable to find exactly one non-diagonal index."];
    }
	
	if (numDensifiableIndices) {
		NSLog(@"Warning: we are copy the matrix to a dense format. This is not efficient.");
		NSMutableArray *matrixFormats = [A.matrixFormats mutableCopy];
		matrixFormats[densifiableIndex] = @(kGLDenseMatrixFormat);
		A = [A copyWithDataType: A.dataFormat matrixFormat: matrixFormats ordering: kGLRowMatrixOrder];
		operandDescription = A.matrixDescription;
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
	
	dispatch_queue_t globalQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
	
	// Should check whether or not the matrix is symmetric.
	GLVariable *eigenvalues = [GLFunction functionOfComplexTypeWithDimensions: eigenbasis forEquation: A.equation];
	GLLinearTransform *eigentransform = [GLLinearTransform transformOfType: kGLSplitComplexDataFormat withFromDimensions: eigenbasis toDimensions:A.fromDimensions inFormat:A.matrixFormats forEquation:A.equation matrix: nil];
	GLMatrixDescription *resultMatrixDescription = eigentransform.matrixDescription;
    GLMatrixDescription *resultVectorDescription = eigenvalues.matrixDescription;
    NSArray *results = @[eigenvalues, eigentransform];
	
	// http://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/sgeev_ex.c.htm
	
	NSUInteger N = operandDescription.strides[denseIndex].nColumns;
    NSUInteger lwork_size = 8*N;
    NSUInteger totalLoops = compute_total_loops(operandDescription, resultMatrixDescription, denseIndex);
    
	// first buffer will be used to store the transpose (which will be overwritten)
	GLBuffer *buffer1 = [[GLBuffer alloc] initWithLength: operandDescription.nBytes];
	// second buffer will store the annoyingly formatted output
	GLBuffer *buffer2 = [[GLBuffer alloc] initWithLength: operandDescription.nBytes];
    // third buffer will store the unstrided eigenvalue output
	GLBuffer *buffer3 = [[GLBuffer alloc] initWithLength: resultVectorDescription.nBytes];
	// fourth buffer is the lapack work buffer
	GLBuffer *buffer4 = [[GLBuffer alloc] initWithLength: totalLoops*lwork_size*sizeof(GLFloat)];
	NSArray *buffers = @[buffer1, buffer2, buffer3, buffer4];
	
	variableOperation op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
		apply_matrix_loop(operandDescription, resultVectorDescription, denseIndex, globalQueue, ^(NSUInteger iteration, NSUInteger inEquationPos, NSUInteger outEquationPos) {
			GLFloat *A = (GLFloat *) [operandArray[0] bytes];
			
			GLFloat *B_data = (GLFloat *) [bufferArray[0] bytes];
            GLFloat *B = &(B_data[iteration*N*N]);
            
			GLFloat *output_data = (GLFloat *) [bufferArray[1] bytes];
            GLFloat *output = &(output_data[iteration*N*N]);
            
            GLFloat *output_v_data = (GLFloat *) [bufferArray[2] bytes];
            GLSplitComplex output_v;
            output_v.realp = &(output_v_data[(2*iteration)*N]);
            output_v.imagp = &(output_v_data[(2*iteration+1)*N]);
            
            GLFloat *work_data = (GLFloat *) [bufferArray[3] bytes];
			GLFloat *work = &(work_data[iteration*lwork_size]);
			
			GLSplitComplex v = splitComplexFromData(resultArray[0]);
			GLSplitComplex C = splitComplexFromData(resultArray[1]);
			
			// clapack takes matrices in column-major format.
			// To be clever, we could use the left-eigenvectors instead of the right-eigenvectors and just take the conjugate, but we need to prevent the input from being overwritten anyway.
			vGL_mtrans(&(A[inEquationPos]), operandDescription.strides[denseIndex].stride, B, 1, N, N);
			
			char JOBVL ='N';
			char JOBVR ='V';
			__CLPK_integer n = (__CLPK_integer) N;
			__CLPK_integer lwork = 8*n;
			__CLPK_integer info;
			
			sgeev_(&JOBVL, &JOBVR, &n, B, &n, output_v.realp, output_v.imagp, NULL, &n, output, &n, work, &lwork, (__CLPK_integer *)&info);
			
			if (info != 0) {
				printf("sgeev failed with error code %d\n", (int)info);
			}
			
			// Now we have to get the eigenvectors in the proper format.
			// If the j-th eigenvalue is real, then v(j) = VR(:,j), the j-th column of VR.
			// If the j-th and (j+1)-st eigenvalues form a complex conjugate pair,
			// then v(j) = VR(:,j) + i*VR(:,j+1) and v(j+1) = VR(:,j) - i*VR(:,j+1).
			// And, don't forget, we need to fix the transpose.
			vGL_vclr( &(C.imagp[inEquationPos]), resultMatrixDescription.strides[denseIndex].stride,  resultMatrixDescription.strides[denseIndex].nPoints);
			
            NSUInteger stride = resultVectorDescription.strides[denseIndex].stride;
            for( NSUInteger i = 0; i < N; i++ ) {
                v.realp[outEquationPos+i*stride] = output_v.realp[i];
                v.imagp[outEquationPos+i*stride] = output_v.imagp[i];
            }
            
            // inEquationPos should have the correct offsets because the output matrix has the same form as the input matrix.
            NSUInteger rowStride = resultMatrixDescription.strides[denseIndex].rowStride;
			NSUInteger colStride = resultMatrixDescription.strides[denseIndex].columnStride;
            NSUInteger j=0;
			for( NSUInteger i = 0; i < N; i++ ) { // i indicates which eigenvector we're copying
				if ( v.imagp[i] == (GLFloat)0.0 ) {
					for( NSUInteger k = 0; k < N; k++ ) { // k walks down the column
						C.realp[inEquationPos+k*rowStride+i*colStride] = output[j*N+k];
					}
					j++;
				} else {
					for( NSUInteger k = 0; k < N; k++ ) {
						C.realp[inEquationPos+k*rowStride+i*colStride] = output[j*N+k];
						C.imagp[inEquationPos+k*rowStride+i*colStride] = output[(j+1)*N+k];
						
						C.realp[inEquationPos+k*rowStride+(i+1)*colStride] = output[j*N+k];
						C.imagp[inEquationPos+k*rowStride+(i+1)*colStride] = -output[(j+1)*N+k];
					}
					j+=2;
					i++;
				}
			}
		});
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
	for (NSUInteger i=0; i<A.fromDimensions.count; i++) {
		if ( ![A.fromDimensions[i] isEqualToDimension: A.toDimensions[i]] ) {
			[NSException raise: @"MatrixWrongFormat" format: @"By assumption, a linear transformation must be an endomorphism to compute the eigensystem."];
		}
	}
	
	NSMutableArray *matrixFormats = [NSMutableArray array];
	NSUInteger candidateIndex = NSNotFound;
	NSUInteger numCandidateIndices = 0;
	BOOL aShouldFormatShift = NO;
	BOOL bShouldFormatShift = NO;
	for ( NSUInteger index=0; index < A.matrixDescription.nDimensions; index++) {
		GLMatrixFormat a = A.matrixDescription.strides[index].matrixFormat;
		GLMatrixFormat b = B.matrixDescription.strides[index].matrixFormat;
		BOOL aIsCandidate = (a != kGLIdentityMatrixFormat && a != kGLDiagonalMatrixFormat);
		BOOL bIsCandidate = (b != kGLIdentityMatrixFormat && b != kGLDiagonalMatrixFormat);
		
		if (aIsCandidate || bIsCandidate) {
			numCandidateIndices++;
			candidateIndex = index;
			matrixFormats[index] = @(kGLDenseMatrixFormat);
		} else { // okay, so both are either identity or diagonal
			if ( a == b ) { // if they're the same, then there's no decision to make.
				matrixFormats[index] = @(a);
			} else { // if they're different, then we need to go with the most general, e.g., diagonal.
				matrixFormats[index] = @(kGLDiagonalMatrixFormat);
			}
		}
		
		aShouldFormatShift |= (a != [matrixFormats[index] unsignedIntegerValue]);
		bShouldFormatShift |= (b != [matrixFormats[index] unsignedIntegerValue]);
	}
	
	if ( numCandidateIndices != 1 ) {
        [NSException raise: @"DenseIndexNotFound" format: @"Unable to find exactly one non-diagonal index."];
    }
	
	if (aShouldFormatShift) {
		A = [A copyWithDataType: A.dataFormat matrixFormat: matrixFormats ordering: kGLRowMatrixOrder];
	}
	
	if (bShouldFormatShift) {
		B = [B copyWithDataType: B.dataFormat matrixFormat: matrixFormats ordering: kGLRowMatrixOrder];
	}
	
	
	GLMatrixDescription *operandDescription = A.matrixDescription;
	
    NSUInteger denseIndex = candidateIndex;
	
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
	
	dispatch_queue_t globalQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
	
	// Should check whether or not the matrix is symmetric.
	GLVariable *eigenvalues = [GLFunction functionOfComplexTypeWithDimensions: eigenbasis forEquation: A.equation];
	GLLinearTransform *eigentransform = [GLLinearTransform transformOfType: kGLSplitComplexDataFormat withFromDimensions: eigenbasis toDimensions:A.toDimensions inFormat:A.matrixFormats forEquation:A.equation matrix: nil];
	GLMatrixDescription *resultMatrixDescription = eigentransform.matrixDescription;
    GLMatrixDescription *resultVectorDescription = eigenvalues.matrixDescription;
    NSArray *results = @[eigenvalues, eigentransform];
    
    // http://www.nag.com/lapack-ex/node122.html
    
	NSUInteger N = operandDescription.strides[denseIndex].nColumns;
    NSUInteger lwork_size = 8*N;
    NSUInteger totalLoops = compute_total_loops(operandDescription, resultMatrixDescription, denseIndex);
	
	// first two buffers will be used to store the transpose (which will be overwritten)
	GLBuffer *buffer1 = [[GLBuffer alloc] initWithLength: operandDescription.nBytes];
    GLBuffer *buffer2 = [[GLBuffer alloc] initWithLength: operandDescription.nBytes];
	// third buffer will store the 'alpha' part of the eigenvalues.
    GLBuffer *buffer3 = [[GLBuffer alloc] initWithLength: resultVectorDescription.nBytes];
    // third buffer will store the 'beta' part of the eigenvalues.
    GLBuffer *buffer4 = [[GLBuffer alloc] initWithLength: resultVectorDescription.nBytes];
	// fourth buffer will store the annoyingly formatted output
	GLBuffer *buffer5 = [[GLBuffer alloc] initWithLength: resultMatrixDescription.nBytes];
	// fifth buffer is the lapack work buffer
	GLBuffer *buffer6 = [[GLBuffer alloc] initWithLength: totalLoops*lwork_size*sizeof(GLFloat)];
	NSArray *buffers = @[buffer1, buffer2, buffer3, buffer4, buffer5, buffer6];
	
	variableOperation op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
		apply_matrix_loop(operandDescription, resultVectorDescription, denseIndex, globalQueue, ^(NSUInteger iteration, NSUInteger inEquationPos, NSUInteger outEquationPos) {
			GLFloat *A_row = (GLFloat *) [operandArray[0] bytes];
			GLFloat *B_row = (GLFloat *) [operandArray[1] bytes];
			
			GLFloat *A_col_data = (GLFloat *) [bufferArray[0] bytes];
			GLFloat *B_col_data = (GLFloat *) [bufferArray[1] bytes];
			GLFloat *A_col = &(A_col_data[iteration*N*N]);
			GLFloat *B_col = &(B_col_data[iteration*N*N]);
			
			GLFloat *output_v_data = (GLFloat *) [bufferArray[2] bytes];
            GLSplitComplex output_v;
            output_v.realp = &(output_v_data[(2*iteration)*N]);
            output_v.imagp = &(output_v_data[(2*iteration+1)*N]);
			
			GLFloat *beta_data = (GLFloat *) [bufferArray[3] bytes];
			GLFloat *beta = &(beta_data[iteration*N]);
			
			GLFloat *output_data = (GLFloat *) [bufferArray[4] bytes];
            GLFloat *output = &(output_data[iteration*N*N]);
			
			GLFloat *work_data = (GLFloat *) [bufferArray[5] bytes];
			GLFloat *work = &(work_data[iteration*lwork_size]);
			
			GLSplitComplex v = splitComplexFromData(resultArray[0]);
			GLSplitComplex C = splitComplexFromData(resultArray[1]);
			
			// clapack takes matrices in column-major format.
			vGL_mtrans(&(A_row[inEquationPos]), operandDescription.strides[denseIndex].stride, A_col, 1, N, N);
			vGL_mtrans(&(B_row[inEquationPos]), operandDescription.strides[denseIndex].stride, B_col, 1, N, N);
			
			char JOBVL ='N';
			char JOBVR ='V';
			__CLPK_integer n = (__CLPK_integer) N;
			__CLPK_integer lwork = 8*n;
			__CLPK_integer info;
			
			sggev_(&JOBVL, &JOBVR, &n, A_col, &n, B_col, &n, output_v.realp, output_v.imagp, beta, NULL, &n, output, &n, work, &lwork, &info);
			
			if (info != 0) {
				printf("sggev failed with error code %d\n", (int)info);
			}
			
			// Now we have to get the eigenvectors in the proper format.
			// If the j-th eigenvalue is real, then v(j) = VR(:,j), the j-th column of VR.
			// If the j-th and (j+1)-st eigenvalues form a complex conjugate pair,
			// then v(j) = VR(:,j) + i*VR(:,j+1) and v(j+1) = VR(:,j) - i*VR(:,j+1).
			// And, don't forget, we need to fix the transpose.
			vGL_vclr( &(C.imagp[inEquationPos]), resultMatrixDescription.strides[denseIndex].stride,  resultMatrixDescription.strides[denseIndex].nPoints);
			
			NSUInteger stride = resultVectorDescription.strides[denseIndex].stride;
            for( NSUInteger i = 0; i < N; i++ ) {
                v.realp[outEquationPos+i*stride] = output_v.realp[i];
                v.imagp[outEquationPos+i*stride] = output_v.imagp[i];
            }
			
			// inEquationPos should have the correct offsets because the output matrix has the same form as the input matrix.
            NSUInteger rowStride = resultMatrixDescription.strides[denseIndex].rowStride;
			NSUInteger colStride = resultMatrixDescription.strides[denseIndex].columnStride;
            NSUInteger j=0;
			for( NSUInteger i = 0; i < N; i++ ) { // i indicates which eigenvector we're copying
				if ( v.imagp[i] == (GLFloat)0.0 ) {
					for( NSUInteger k = 0; k < N; k++ ) { // k walks down the column
						C.realp[inEquationPos+k*rowStride+i*colStride] = output[j*N+k];
					}
					j++;
				} else {
					for( NSUInteger k = 0; k < N; k++ ) {
						C.realp[inEquationPos+k*rowStride+i*colStride] = output[j*N+k];
						C.imagp[inEquationPos+k*rowStride+i*colStride] = output[(j+1)*N+k];
						
						C.realp[inEquationPos+k*rowStride+(i+1)*colStride] = output[j*N+k];
						C.imagp[inEquationPos+k*rowStride+(i+1)*colStride] = -output[(j+1)*N+k];
					}
					j+=2;
					i++;
				}
			}

			
			vGL_vdiv(beta, 1, &(v.realp[outEquationPos]), stride, &(v.realp[outEquationPos]), stride, N);
			vGL_vdiv(beta, 1, &(v.imagp[outEquationPos]), stride, &(v.imagp[outEquationPos]), stride, N);
		});
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
	NSMutableArray *matrixDescriptions = [NSMutableArray array];
	NSMutableArray *numDimensions = [NSMutableArray array];
    BOOL isComplex = NO;
    for (GLLinearTransform *transform in linearTransformations) {
        [fromDimensions addObjectsFromArray: transform.fromDimensions];
		[toDimensions addObjectsFromArray: transform.toDimensions];
		[numDimensions addObject: @(fromDimensions.count)];
        isComplex |= transform.isComplex;
        [matrixFormat addObjectsFromArray: transform.matrixFormats];
		[matrixDescriptions addObject: transform.matrixDescription];
		if (transform.matrixBlock) {
			[matrixBlocks addObject: transform.matrixBlock];
		} else {
			[matrixBlocks addObject: [NSNull null]];
		}
		
    }
    GLDataFormat format = isComplex ? kGLSplitComplexDataFormat : kGLRealDataFormat;
    
	GLLinearTransform *tensorProductTransform = [GLLinearTransform transformOfType:format withFromDimensions: fromDimensions toDimensions: toDimensions inFormat:matrixFormat forEquation:[linearTransformations[0] equation] matrix: NULL];
	GLMatrixDescription *tensorProductMatrixDescription = tensorProductTransform.matrixDescription;
	
	variableOperation op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
		// create a matrix block for any transformation that doesn't have one.
		for (NSUInteger i=0;i<matrixBlocks.count;i++) {
			if (matrixBlocks[i] == [NSNull null]) {
				matrixBlocks[i] = [GLLinearTransform matrixBlockWithFormat: matrixDescriptions[i] fromData: operandArray[i]];
			}
		}
		
		// create the tensor product matrix block.
		transformMatrix tensorProduct = ^( NSUInteger *row, NSUInteger *col ) {
			NSUInteger iMatrix = 0;
			NSUInteger iDimension = 0;
			transformMatrix theMatrixBlock = matrixBlocks[iMatrix];
			GLFloatComplex value = theMatrixBlock(&(row[iDimension]), &(col[iDimension]));
			iDimension += [numDimensions[iMatrix] unsignedIntegerValue];
			for (iMatrix=1;iMatrix<matrixBlocks.count;iMatrix++) {
				theMatrixBlock = matrixBlocks[iMatrix];
				value *= theMatrixBlock(&(row[iDimension]), &(col[iDimension]));
				iDimension += [numDimensions[iMatrix] unsignedIntegerValue];
			}
			return value;
		};
		
		[GLLinearTransform writeToData: resultArray[0] withFormat: tensorProductMatrixDescription fromMatrixBlock: tensorProduct];
	};
		
	if (( self = [super initWithResult: @[tensorProductTransform] operand: linearTransformations buffers: @[] operation: op] )) {
		
	}
	
	return self;
}

@end

/************************************************/
/*		GLFormatShiftOperation					*/
/************************************************/

#pragma mark -
#pragma mark GLFormatShiftOperation
#pragma mark

@implementation GLFormatShiftOperation

- (id) initWithLinearTransformation: (GLLinearTransform *) linearTransform dataType: (GLDataFormat) dataFormat matrixFormat: (NSArray *) matrixFormats ordering: (GLMatrixOrder) ordering
{
	BOOL identityMatrix = YES;
	for (NSNumber *format in linearTransform.matrixFormats) {
		identityMatrix &= format.unsignedIntegerValue == kGLIdentityMatrixFormat;
	}
	
	GLLinearTransform *newLinearTransform = [[GLLinearTransform alloc] initTransformOfType: dataFormat withFromDimensions: linearTransform.fromDimensions toDimensions: linearTransform.toDimensions inFormat: matrixFormats withOrdering: ordering forEquation: linearTransform.equation matrix: nil];
	
	transformMatrix matrixBlock = linearTransform.matrixBlock;
	GLMatrixDescription *oldMatrixDescription = linearTransform.matrixDescription;
	GLMatrixDescription *newMatrixDescription = newLinearTransform.matrixDescription;
	
	variableOperation op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
		transformMatrix matrix = matrixBlock ? matrixBlock : [GLLinearTransform matrixBlockWithFormat: oldMatrixDescription fromData: operandArray[0]];
		[GLLinearTransform writeToData: resultArray[0] withFormat: newMatrixDescription fromMatrixBlock: matrix];
	};
	
	NSArray *operandArray = matrixBlock ? [NSArray array] : @[linearTransform];
	
	if (( self = [super initWithResult: @[newLinearTransform] operand: operandArray buffers: @[] operation: op] )) {
		
	}
	
	return self;
}

@end

