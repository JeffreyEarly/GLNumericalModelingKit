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

@implementation GLTriadiagonalOperation

// This is copy and pasted from the superclass, needs to be properly retooled.
- (id) initWithResult: (GLVariable *) resultVariable firstOperand: (GLVariable *) fOperand secondOperand: (GLVariable *) sOperand;
{
	if ( ![fOperand.class isSubclassOfClass: [GLLinearTransform class]] ) {
        [NSException raise: @"BadArgument" format: @"The first argument must be a GLLinearTransform."];
    }
    
    GLLinearTransform *linearTransform = (GLLinearTransform *) fOperand;
    
    if ( ![linearTransform.toDimensions isEqualToArray: sOperand.dimensions] ) {
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

    
	if (( self = [super init] )) {
		self.firstOperand = fOperand;
		self.secondOperand = sOperand;
		
		if (!resultVariable) {
			BOOL isComplex = fOperand.isComplex || sOperand.isComplex;
			GLDataFormat format = isComplex ? kGLSplitComplexDataFormat : kGLRealDataFormat;
			
			self.result = [GLVariable variableOfType: format withDimensions: linearTransform.fromDimensions forEquation: self.firstOperand.equation];
		} else {
			self.result = resultVariable;
		}
		
		[self setupDependencies];
        
		
        GLMatrixDescription *matrixDescription = linearTransform.matrixDescription;
        
#warning non of this deals with complex numbers yet.
		
#warning this needs the same outer loop as the transform below.
        
        // This buffer is 3 times as large as it needs to be, but it makes bookkeeping easier.
        NSMutableData *buffer = [[GLMemoryPool sharedMemoryPool] dataWithLength: self.result.nDataPoints*sizeof(GLFloat)];
		
        dispatch_queue_t globalQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
        
        NSUInteger inNDiagonalPoints = matrixDescription.strides[triIndex].nDiagonalPoints;
        NSUInteger inElementStride = matrixDescription.strides[triIndex].stride;
        NSUInteger inDiagonalStride = matrixDescription.strides[triIndex].diagonalStride;
        NSUInteger inEquationStride = lastNonTriIndex == NSNotFound ? 0 : matrixDescription.strides[lastNonTriIndex].stride;
        NSUInteger outElementStride = self.result.matrixDescription.strides[triIndex].stride;
        NSUInteger outEquationStride = lastNonTriIndex == NSNotFound ? 0 : self.result.matrixDescription.strides[lastNonTriIndex].stride;
        NSUInteger totalEquations = linearTransform.nDataPoints / matrixDescription.strides[triIndex].nPoints;
        
        self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
            
            dispatch_apply(totalEquations, globalQueue, ^(size_t iteration) {

                NSInteger inEquationPos = iteration*inEquationStride;
                NSInteger outEquationPos = iteration*outEquationStride;
                
                GLFloat *f = (GLFloat *) fOperand.bytes;
                
                // plus the offset!!!!
                GLFloat *a = &(f[0*inDiagonalStride]);
                GLFloat *b = &(f[1*inDiagonalStride]);
                GLFloat *c = &(f[2*inDiagonalStride]);
                
                GLFloat *cprime = (GLFloat *) buffer.mutableBytes;
                
                GLFloat *d = (GLFloat *) sOperand.bytes;
                GLFloat *x = (GLFloat *) result.mutableBytes;                
                
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

- (void) setupDependencies
{
	if (self.firstOperand.lastOperation && ![self.dependencies containsObject:self.firstOperand.lastOperation]) {
		[self addDependency: self.firstOperand.lastOperation];
	}
	if (self.secondOperand.lastOperation && ![self.dependencies containsObject:self.secondOperand.lastOperation]) {
		[self addDependency: self.secondOperand.lastOperation];
	}
	if (self.result.lastOperation && ![self.dependencies containsObject:self.result.lastOperation]) {
		[self addDependency: self.result.lastOperation];
	}
	
	[self.result addOperation: self];
}

@end


/************************************************/
/*		GLTriadiagonalTransformOperation		*/
/************************************************/

#pragma mark -
#pragma mark GLTriadiagonalTransformOperation
#pragma mark

@implementation GLTriadiagonalTransformOperation

// This is copy and pasted from the superclass, needs to be properly retooled.
- (id) initWithResult: (GLVariable *) resultVariable firstOperand: (GLVariable *) fOperand secondOperand: (GLVariable *) sOperand;
{
	if ( ![fOperand.class isSubclassOfClass: [GLLinearTransform class]] ) {
        [NSException raise: @"BadArgument" format: @"The first argument must be a GLLinearTransform."];
    }
    
    GLLinearTransform *linearTransform = (GLLinearTransform *) fOperand;
    
    if ( ![linearTransform.fromDimensions isEqualToArray: sOperand.dimensions] ) {
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
	
    
	if (( self = [super init] )) {
		self.firstOperand = fOperand;
		self.secondOperand = sOperand;
		
		if (!resultVariable) {
			BOOL isComplex = fOperand.isComplex || sOperand.isComplex;
			GLDataFormat format = isComplex ? kGLSplitComplexDataFormat : kGLRealDataFormat;
			
			self.result = [GLVariable variableOfType: format withDimensions: linearTransform.toDimensions forEquation: self.firstOperand.equation];
		} else {
			self.result = resultVariable;
		}
		
		[self setupDependencies];
        
		
        GLMatrixDescription *matrixDescription = linearTransform.matrixDescription;
        
#warning non of this deals with complex numbers yet.
		
        dispatch_queue_t globalQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
        
        NSUInteger inNDiagonalPoints = matrixDescription.strides[triIndex].nDiagonalPoints;
        NSUInteger inElementStride = matrixDescription.strides[triIndex].stride;
        NSUInteger inDiagonalStride = matrixDescription.strides[triIndex].diagonalStride;
		
		// How many 'inner loop' steps we need to take depends on how many other nontrivial dimensions there are.
        NSUInteger inEquationStride = lastNonTrivialNonTriIndex == NSNotFound ? 0 : matrixDescription.strides[lastNonTrivialNonTriIndex].stride;
		NSUInteger totalEquations = linearTransform.nDataPoints / matrixDescription.strides[triIndex].nPoints;
		
		// Now we need the strides to match up to the inner loop
        NSUInteger outElementStride = self.result.matrixDescription.strides[triIndex].stride;
        NSUInteger outEquationStride = lastNonTrivialNonTriIndex == NSNotFound ? 0 : self.result.matrixDescription.strides[lastNonTrivialNonTriIndex].stride;
        
		// Finally, we need the outer loop strides and totals
		NSUInteger totalOuterLoops = totalTrivialPoints;
		NSUInteger outerLoopStride = lastTrivialIndex == NSNotFound ? 0 : self.result.matrixDescription.strides[lastTrivialIndex].stride;
		
		// This operation has two loops.
		// The inner loop walks over non-trivial elements of the linear operator, while
		// the out loop walks over trivial elements of the linear operator, really just to new x and b elements
        
        self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
            
			dispatch_apply(totalOuterLoops, globalQueue, ^(size_t outerIteration) {
				
				NSInteger outerOffset = outerIteration*outerLoopStride;
				
				dispatch_apply(totalEquations, globalQueue, ^(size_t iteration) {
					
					NSInteger inEquationPos = iteration*inEquationStride;
					NSInteger outEquationPos = outerOffset + iteration*outEquationStride;
					
					GLFloat *f = (GLFloat *) fOperand.bytes;
					
					// plus the offset!!!!
					GLFloat *a = &(f[0*inDiagonalStride]);
					GLFloat *b = &(f[1*inDiagonalStride]);
					GLFloat *c = &(f[2*inDiagonalStride]);
					
					GLFloat *x = (GLFloat *) sOperand.bytes;
					GLFloat *d = (GLFloat *) result.mutableBytes;
					
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

- (void) setupDependencies
{
	if (self.firstOperand.lastOperation && ![self.dependencies containsObject:self.firstOperand.lastOperation]) {
		[self addDependency: self.firstOperand.lastOperation];
	}
	if (self.secondOperand.lastOperation && ![self.dependencies containsObject:self.secondOperand.lastOperation]) {
		[self addDependency: self.secondOperand.lastOperation];
	}
	if (self.result.lastOperation && ![self.dependencies containsObject:self.result.lastOperation]) {
		[self addDependency: self.result.lastOperation];
	}
	
	[self.result addOperation: self];
}

@end

/************************************************/
/*		GLDenseMatrixTransformOperation		*/
/************************************************/

#pragma mark -
#pragma mark GLDenseMatrixTransformOperation
#pragma mark

@implementation GLDenseMatrixTransformOperation

// This is copy and pasted from the superclass, needs to be properly retooled.
- (id) initWithResult: (GLVariable *) resultVariable firstOperand: (GLVariable *) fOperand secondOperand: (GLVariable *) sOperand;
{
	if ( ![fOperand.class isSubclassOfClass: [GLLinearTransform class]] ) {
        [NSException raise: @"BadArgument" format: @"The first argument must be a GLLinearTransform."];
    }
    
    GLLinearTransform *linearTransform = (GLLinearTransform *) fOperand;
    
    if ( ![linearTransform.fromDimensions isEqualToArray: sOperand.dimensions] ) {
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

    if (( self = [super init] )) {
		self.firstOperand = fOperand;
		self.secondOperand = sOperand;
		
		if (!resultVariable) {
			BOOL isComplex = fOperand.isComplex || sOperand.isComplex;
			GLDataFormat format = isComplex ? kGLSplitComplexDataFormat : kGLRealDataFormat;
			
			self.result = [GLVariable variableOfType: format withDimensions: linearTransform.toDimensions forEquation: self.firstOperand.equation];
		} else {
			self.result = resultVariable;
		}
		
		[self setupDependencies];
        
		
        GLMatrixDescription *matrixDescription = linearTransform.matrixDescription;
        
//        int M = (int) matrixDescription.strides[0].nRows;
//        int N = (int) matrixDescription.strides[0].nColumns;
//        self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
//            cblas_sgemv( CblasRowMajor,  CblasNoTrans, M, N, 1.0, fOperand.bytes, M, sOperand.bytes, 1.0, 1.0, result.mutableBytes, 1);
//        };
		
		int M = (int) matrixDescription.strides[0].nRows;
        int N = (int) 1;
		int K = (int) matrixDescription.strides[0].nColumns;
		if ( !self.firstOperand.isComplex && !self.secondOperand.isComplex)
		{	// C = A.X
			self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
				vDSP_mmul( (GLFloat *)fOperand.bytes, 1, (GLFloat *)sOperand.bytes, 1, result.mutableBytes, 1, M, N, K);
			};
		}
		else if ( self.firstOperand.isComplex && !self.secondOperand.isComplex)
		{	// (A+iB).(X) = A.X + iB.X
			self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
				GLSplitComplex leftComplex = splitComplexFromData( fOperand );
				GLSplitComplex destComplex = splitComplexFromData( result );
				
				vDSP_mmul( leftComplex.realp, 1, (GLFloat *)sOperand.bytes, 1, destComplex.realp, 1, M, N, K);
				vDSP_mmul( leftComplex.imagp, 1, (GLFloat *)sOperand.bytes, 1, destComplex.imagp, 1, M, N, K);
			};
		}
		else if ( !self.firstOperand.isComplex && self.secondOperand.isComplex)
		{	// A.(X+iY) = A.X + iA.Y
			self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
				GLSplitComplex rightComplex = splitComplexFromData( sOperand );
				GLSplitComplex destComplex = splitComplexFromData( result );
				
				vDSP_mmul( (GLFloat *)fOperand.bytes, 1, rightComplex.realp, 1, destComplex.realp, 1, M, N, K);
				vDSP_mmul( (GLFloat *)fOperand.bytes, 1, rightComplex.imagp, 1, destComplex.imagp, 1, M, N, K);
			};
		}
		else if ( self.firstOperand.isComplex && self.secondOperand.isComplex)
		{
			self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
				GLSplitComplex leftComplex = splitComplexFromData( fOperand );
				GLSplitComplex rightComplex = splitComplexFromData( sOperand );
				GLSplitComplex destComplex = splitComplexFromData( result );
				
				vDSP_zmmul( &leftComplex, 1, &rightComplex, 1, &destComplex, 1, M, N, K);
			};
		}
    }
    return self;
}

- (void) setupDependencies
{
	if (self.firstOperand.lastOperation && ![self.dependencies containsObject:self.firstOperand.lastOperation]) {
		[self addDependency: self.firstOperand.lastOperation];
	}
	if (self.secondOperand.lastOperation && ![self.dependencies containsObject:self.secondOperand.lastOperation]) {
		[self addDependency: self.secondOperand.lastOperation];
	}
	if (self.result.lastOperation && ![self.dependencies containsObject:self.result.lastOperation]) {
		[self addDependency: self.result.lastOperation];
	}
	
	[self.result addOperation: self];
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
- (id) initWithResult: (GLVariable *) resultVariable firstOperand: (GLVariable *) fOperand secondOperand: (GLVariable *) sOperand;
{
	if ( ![fOperand.class isSubclassOfClass: [GLLinearTransform class]] ) {
        [NSException raise: @"BadArgument" format: @"The first argument must be a GLLinearTransform."];
    }
    
    if ( ![sOperand.class isSubclassOfClass: [GLLinearTransform class]] ) {
        [NSException raise: @"BadArgument" format: @"The second argument must be a GLLinearTransform."];
    }
    
    GLLinearTransform *A = (GLLinearTransform *) fOperand;
    GLLinearTransform *B = (GLLinearTransform *) sOperand;
    
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
    
    if (( self = [super init] )) {
		self.firstOperand = fOperand;
		self.secondOperand = sOperand;
		
		if (!resultVariable) {
			BOOL isComplex = fOperand.isComplex || sOperand.isComplex;
			GLDataFormat format = isComplex ? kGLSplitComplexDataFormat : kGLRealDataFormat;
			
			self.result = [GLLinearTransform transformOfType: format withFromDimensions: B.fromDimensions toDimensions: A.toDimensions inFormat:B.matrixFormats forEquation: A.equation];
		} else {
			self.result = resultVariable;
		}
		
		[self setupDependencies];
        
		GLLinearTransform *C = (GLLinearTransform *) self.result;
		        
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
		
		if ( !self.firstOperand.isComplex && !self.secondOperand.isComplex)
		{	// C = A.X
			self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
				vDSP_mmul( (GLFloat *)fOperand.bytes, 1, (GLFloat *)sOperand.bytes, 1, result.mutableBytes, 1, M, N, K);
			};
		}
		else if ( self.firstOperand.isComplex && !self.secondOperand.isComplex)
		{	// (A+iB).(X) = A.X + iB.X
			self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
				GLSplitComplex leftComplex = splitComplexFromData( fOperand );
				GLSplitComplex destComplex = splitComplexFromData( result );
				
				vDSP_mmul( leftComplex.realp, 1, (GLFloat *)sOperand.bytes, 1, destComplex.realp, 1, M, N, K);
				vDSP_mmul( leftComplex.imagp, 1, (GLFloat *)sOperand.bytes, 1, destComplex.imagp, 1, M, N, K);
			};
		}
		else if ( !self.firstOperand.isComplex && self.secondOperand.isComplex)
		{	// A.(X+iY) = A.X + iA.Y
			self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
				GLSplitComplex rightComplex = splitComplexFromData( sOperand );
				GLSplitComplex destComplex = splitComplexFromData( result );
				
				vDSP_mmul( (GLFloat *)fOperand.bytes, 1, rightComplex.realp, 1, destComplex.realp, 1, M, N, K);
				vDSP_mmul( (GLFloat *)fOperand.bytes, 1, rightComplex.imagp, 1, destComplex.imagp, 1, M, N, K);
			};
		}
		else if ( self.firstOperand.isComplex && self.secondOperand.isComplex)
		{
			self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
				GLSplitComplex leftComplex = splitComplexFromData( fOperand );
				GLSplitComplex rightComplex = splitComplexFromData( sOperand );
				GLSplitComplex destComplex = splitComplexFromData( result );
				
				vDSP_zmmul( &leftComplex, 1, &rightComplex, 1, &destComplex, 1, M, N, K);
			};
		}
		
    }
    return self;
}

- (void) setupDependencies
{
	if (self.firstOperand.lastOperation && ![self.dependencies containsObject:self.firstOperand.lastOperation]) {
		[self addDependency: self.firstOperand.lastOperation];
	}
	if (self.secondOperand.lastOperation && ![self.dependencies containsObject:self.secondOperand.lastOperation]) {
		[self addDependency: self.secondOperand.lastOperation];
	}
	if (self.result.lastOperation && ![self.dependencies containsObject:self.result.lastOperation]) {
		[self addDependency: self.result.lastOperation];
	}
	
	[self.result addOperation: self];
}

@end

/************************************************/
/*		GLMatrixInversionOperation              */
/************************************************/

#pragma mark -
#pragma mark GLMatrixInversionOperation
#pragma mark

@implementation GLMatrixInversionOperation

- (id) initWithResult:(GLVariable *)resultVariable operand:(GLVariable *)variable
{
    if ( ![variable.class isSubclassOfClass: [GLLinearTransform class]] ) {
        [NSException raise: @"BadArgument" format: @"The argument must be a GLLinearTransform."];
    }
    
    if (variable.matrixDescription.nDimensions != 1) {
        [NSException raise: @"MatrixWrongFormat" format: @"We can only do one dimensional matrices at the moment."];
    }
    
    GLLinearTransform *A = (GLLinearTransform *) variable;
    if (A.fromDimensions.count != A.toDimensions.count ) {
        [NSException raise: @"MatrixWrongFormat" format: @"We can only do square matrices."];
    }
    
    if ((self=[super init]))
    {
        self.operand = variable;
        if (!resultVariable) {
			GLDataFormat format = self.operand.isComplex ? kGLSplitComplexDataFormat : kGLRealDataFormat;
			
			self.result = [GLLinearTransform transformOfType: format withFromDimensions: A.toDimensions toDimensions: A.fromDimensions inFormat:A.matrixFormats forEquation: A.equation];
		} else {
			self.result = resultVariable;
		}
        
        [self setupDependencies];
        
        NSUInteger N = [A.fromDimensions[0] nPoints];
        self.blockOperation = ^(NSMutableData *result, NSData *operand) {
            
            // clapack takes matrices in column-major format.
            // However, the transpose of the inverse is the inverse of the transpose---so we don't need to worry here.
            
            NSMutableData *ipiv = [[GLMemoryPool sharedMemoryPool] dataWithLength: N*sizeof(__CLPK_integer)];
            __CLPK_integer n = (__CLPK_integer) N;
            __CLPK_integer info;
            memcpy( result.mutableBytes, operand.bytes, n*n*sizeof(GLFloat));
            sgetrf_(&n, &n, result.mutableBytes, &n, ipiv.mutableBytes, (__CLPK_integer *)&info);
            
            if (info != 0) {
                printf("sgetrf failed with error code %d\n", (int)info);
            }
            
            __CLPK_integer lwork = n*n;
            NSMutableData *work = [[GLMemoryPool sharedMemoryPool] dataWithLength: lwork*sizeof(GLFloat)];
            sgetri_(&n, result.mutableBytes, &n, ipiv.mutableBytes, work.mutableBytes, &lwork, (__CLPK_integer *)&info);
            
            if (info != 0) {
                printf("sgetri failed with error code %d\n", (int)info);
            }
            
            [[GLMemoryPool sharedMemoryPool] returnData: ipiv];
            [[GLMemoryPool sharedMemoryPool] returnData: work];
        };
    }
    return self;
}

- (void) setupDependencies
{
	if (self.operand.lastOperation && ![self.dependencies containsObject:self.operand.lastOperation]) {
		[self addDependency: self.operand.lastOperation];
	}
	if (self.result.lastOperation && ![self.dependencies containsObject:self.result.lastOperation]) {
		[self addDependency: self.result.lastOperation];
	}
	
	[self.result addOperation: self];
}

@end

/************************************************/
/*		GLLinearTransformAdditionOperation		*/
/************************************************/

@implementation GLLinearTransformAdditionOperation

- (id) initWithResult: (GLVariable *) resultVariable firstOperand: (GLVariable *) fOperand secondOperand: (GLVariable *) sOperand;
{
    if ( ![fOperand.class isSubclassOfClass: [GLLinearTransform class]] ) {
        [NSException raise: @"BadArgument" format: @"The first argument must be a GLLinearTransform."];
    }
    
    if ( ![sOperand.class isSubclassOfClass: [GLLinearTransform class]] ) {
        [NSException raise: @"BadArgument" format: @"The second argument must be a GLLinearTransform."];
    }
    
    GLLinearTransform *A = (GLLinearTransform *) fOperand;
    GLLinearTransform *B = (GLLinearTransform *) sOperand;
    
    if ( ![A.fromDimensions isEqualToArray: B.fromDimensions] ) {
        [NSException raise: @"DimensionsNotEqualException" format: @"fromDimensions of A, must equal the fromDimensions of B."];
    }
    
    if ( ![A.toDimensions isEqualToArray: B.toDimensions] ) {
        [NSException raise: @"DimensionsNotEqualException" format: @"toDimensions of A, must equal the toDimensions of B."];
    }
    
    if (( self = [super init] )) {
		self.firstOperand = fOperand;
		self.secondOperand = sOperand;
		
		if (!resultVariable) {
			BOOL isComplex = fOperand.isComplex || sOperand.isComplex;
			GLDataFormat format = isComplex ? kGLSplitComplexDataFormat : kGLRealDataFormat;
			
			self.result = [GLLinearTransform transformOfType: format withFromDimensions: A.fromDimensions toDimensions: A.toDimensions inFormat:A.matrixFormats forEquation: A.equation];
		} else {
			self.result = resultVariable;
		}
		
		[self setupDependencies];
        
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
    }
    
    return self;   
}

- (BOOL) canOperateInPlace {
	return YES;
}

- (void) setupDependencies
{
	if (self.firstOperand.lastOperation && ![self.dependencies containsObject:self.firstOperand.lastOperation]) {
		[self addDependency: self.firstOperand.lastOperation];
	}
	if (self.secondOperand.lastOperation && ![self.dependencies containsObject:self.secondOperand.lastOperation]) {
		[self addDependency: self.secondOperand.lastOperation];
	}
	if (self.result.lastOperation && ![self.dependencies containsObject:self.result.lastOperation]) {
		[self addDependency: self.result.lastOperation];
	}
	
	[self.result addOperation: self];
}

@end

/************************************************/
/*		GLDenseMatrixSolver						*/
/************************************************/

#pragma mark -
#pragma mark GLDenseMatrixSolver
#pragma mark

@implementation GLDenseMatrixSolver

// This is copy and pasted from the superclass, needs to be properly retooled.
- (id) initWithResult: (GLVariable *) resultVariable firstOperand: (GLVariable *) fOperand secondOperand: (GLVariable *) sOperand;
{
	if ( ![fOperand.class isSubclassOfClass: [GLLinearTransform class]] ) {
        [NSException raise: @"BadArgument" format: @"The first argument must be a GLLinearTransform."];
    }
    
    GLLinearTransform *linearTransform = (GLLinearTransform *) fOperand;
    
    if ( ![linearTransform.toDimensions isEqualToArray: sOperand.dimensions] ) {
        [NSException raise: @"DimensionsNotEqualException" format: @"To dimensions of the linear operator must equal the operand vector."];
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
	
    
	if (( self = [super init] )) {
		self.firstOperand = fOperand;
		self.secondOperand = sOperand;
		
		if (!resultVariable) {
			BOOL isComplex = fOperand.isComplex || sOperand.isComplex;
			GLDataFormat format = isComplex ? kGLSplitComplexDataFormat : kGLRealDataFormat;
			
			self.result = [GLVariable variableOfType: format withDimensions: linearTransform.toDimensions forEquation: self.firstOperand.equation];
		} else {
			self.result = resultVariable;
		}
		
		[self setupDependencies];
        
		
        GLMatrixDescription *matrixDescription = linearTransform.matrixDescription;
        
#warning non of this deals with complex numbers yet.
		
        dispatch_queue_t globalQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
        		
		// How many 'inner loop' steps we need to take depends on how many other nontrivial dimensions there are.
        NSUInteger inEquationStride = lastNonTrivialNonDenseIndex == NSNotFound ? 0 : matrixDescription.strides[lastNonTrivialNonDenseIndex].stride;
		NSUInteger totalEquations = linearTransform.nDataPoints / matrixDescription.strides[denseIndex].nPoints;
		
		// Now we need the strides to match up to the inner loop
        NSUInteger outElementStride = self.result.matrixDescription.strides[denseIndex].stride;
        NSUInteger outEquationStride = lastNonTrivialNonDenseIndex == NSNotFound ? 0 : self.result.matrixDescription.strides[lastNonTrivialNonDenseIndex].stride;
        
		// Finally, we need the outer loop strides and totals
		// The totalTrivialPoints is the number of times we need to repeat a point
		// trivialPointStride is the distance between those trivial points
		// totalNonTrivialPoints is the number of points we need to copy
		NSUInteger totalNonTrivialPoints = self.result.nDataPoints / totalTrivialPoints;
		NSUInteger trivialPointStride = lastTrivialIndex == NSNotFound ? 0 : self.result.matrixDescription.strides[lastTrivialIndex].stride;
		
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
		
        self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
				
				dispatch_apply(totalEquations, globalQueue, ^(size_t iteration) {
					
					NSInteger inEquationPos = iteration*inEquationStride;
					NSInteger outEquationPos = iteration*outEquationStride;
					
					NSMutableData *ipiv = [[GLMemoryPool sharedMemoryPool] dataWithLength: N*sizeof(__CLPK_integer)];
					__CLPK_integer n = (__CLPK_integer) N;
					__CLPK_integer nrhs = 1;
					__CLPK_integer info;
					
					GLFloat *MData = (GLFloat *) fOperand.bytes;
					GLFloat *bData = (GLFloat *) sOperand.bytes;
					GLFloat *xData = (GLFloat *) result.bytes;
					
					GLFloat *M;
					GLFloat *b;
					GLFloat *x;
					if ( lastNonTrivialNonDenseIndex != NSNotFound && lastNonTrivialNonDenseIndex > denseIndex) {
						
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
							GLFloat *theOutput = (GLFloat *) result.mutableBytes;
							vGL_vfill( &theOutput[ outIteration*outElementStride], &theOutput[outIteration*outElementStride], trivialPointStride, totalTrivialPoints);
						});
						
					}
					
				});
        };
        
    }
    return self;
}

- (void) setupDependencies
{
	if (self.firstOperand.lastOperation && ![self.dependencies containsObject:self.firstOperand.lastOperation]) {
		[self addDependency: self.firstOperand.lastOperation];
	}
	if (self.secondOperand.lastOperation && ![self.dependencies containsObject:self.secondOperand.lastOperation]) {
		[self addDependency: self.secondOperand.lastOperation];
	}
	if (self.result.lastOperation && ![self.dependencies containsObject:self.result.lastOperation]) {
		[self addDependency: self.result.lastOperation];
	}
	
	[self.result addOperation: self];
}

@end
