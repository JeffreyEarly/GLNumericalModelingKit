//
//  GLVectorVectorOperations.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 1/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import "GLVectorVectorOperations.h"
#import "GLMatrixDescription.h"
#import "GLScalar.h"
#import "GLLinearTransform.h"

/************************************************/
/*		GLAdditionOperation						*/
/************************************************/

@implementation GLAdditionOperation

- (id) initWithFirstOperand: (GLTensor *) fOperand secondOperand: (GLTensor *) sOperand;
{
    // We order the operands so that scalars are always in the first position.
    // We can do this in this case because order doesn't matter for addition.
    GLDataFormat format = (fOperand.isPurelyReal && sOperand.isPurelyReal) ? kGLRealDataFormat : kGLSplitComplexDataFormat;
    GLTensor *op1 = (sOperand.rank < fOperand.rank) ? sOperand : fOperand;
    GLTensor *op2 = (sOperand.rank < fOperand.rank) ? fOperand : sOperand;
    GLTensor *result;
	variableOperation operation;
	NSString *graphvisDescription;
	
	// Lots of cases to parse through.
	if (op2.rank == 0)
	{	// scalar-scalar
		result = [[GLScalar alloc] initWithType: format forEquation: op1.equation];
		result.isPurelyReal = fOperand.isPurelyReal && sOperand.isPurelyReal;
		result.isPurelyImaginary = fOperand.isPurelyImaginary && sOperand.isPurelyImaginary;
		
		if ( !op1.isComplex && !op2.isComplex) {
			graphvisDescription = [NSString stringWithFormat: @"add (real scalar, real scalar)"];
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *a = (GLFloat *) [operandArray[0] bytes];
				GLFloat *b = (GLFloat *) [operandArray[1] bytes];
				GLFloat *c = (GLFloat *) [resultArray[0] bytes];
				*c = (*a) + (*b);
			};
		} else if ( !op1.isComplex && op2.isComplex) {
			graphvisDescription = [NSString stringWithFormat: @"add (real scalar, complex scalar)"];
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *a = (GLFloat *) [operandArray[0] bytes];
				GLFloatComplex *b = (GLFloatComplex *) [operandArray[1] bytes];
				GLFloatComplex *c = (GLFloatComplex *) [resultArray[0] bytes];
				*c = (*a) + (*b);
			};
		} else if ( op1.isComplex && !op2.isComplex) {
			graphvisDescription = [NSString stringWithFormat: @"add (complex scalar, real scalar)"];
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloatComplex *a = (GLFloatComplex *) [operandArray[0] bytes];
				GLFloat *b = (GLFloat *) [operandArray[1] bytes];
				GLFloatComplex *c = (GLFloatComplex *) [resultArray[0] bytes];
				*c = (*a) + (*b);
			};
		} else {
			graphvisDescription = [NSString stringWithFormat: @"add (complex scalar, complex scalar)"];
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloatComplex *a = (GLFloatComplex *) [operandArray[0] bytes];
				GLFloatComplex *b = (GLFloatComplex *) [operandArray[1] bytes];
				GLFloatComplex *c = (GLFloatComplex *) [resultArray[0] bytes];
				*c = (*a) + (*b);
			};
		}
		
	}
	else if (op2.rank == 1 || op2.rank == 2)
	{
		// Adding functions and tensors share the same implementation details, so we divide the cases up here.
		NSUInteger implementationCase = 0;
		NSUInteger numPoints;
		NSUInteger nDataElements;
		
		if (op2.rank == 1)
		{	// scalar-function or function-function
			GLVariable *func2 = (GLVariable *) op2;
			result = [GLVariable variableOfType:format withDimensions: func2.dimensions forEquation: op2.equation];
			result.isPurelyReal = fOperand.isPurelyReal && sOperand.isPurelyReal;
			result.isPurelyImaginary = fOperand.isPurelyImaginary && sOperand.isPurelyImaginary;
			
			numPoints = result.nDataPoints;
			nDataElements = result.nDataElements;
			
			if (op1.rank == 0)
			{		// scalar-function
				if ( !op1.isComplex && !op2.isComplex) {
					graphvisDescription = [NSString stringWithFormat: @"add (real scalar, real function)"];
					implementationCase = 1;
				} else if ( !op1.isComplex && op2.isComplex) {
					graphvisDescription = [NSString stringWithFormat: @"add (real scalar, complex function)"];
					implementationCase = 2;
				} else if ( op1.isComplex && !op2.isComplex) {
					graphvisDescription = [NSString stringWithFormat: @"add (complex scalar, real function)"];
					implementationCase = 3;
				} else {
					graphvisDescription = [NSString stringWithFormat: @"add (complex scalar, complex function)"];
					implementationCase = 4;
				}
			}
			else if (op1.rank == 1)
			{	// function-function
				GLVariable *func1 = (GLVariable *) op1;
				if ( ![func1.dimensions isEqualToArray: func2.dimensions] ) {
					[NSException raise: @"DimensionsNotEqualException" format: @"Cannot add two functions of different dimensions"];
				}
				
				if ( !op1.isComplex && !op2.isComplex) {
					graphvisDescription = [NSString stringWithFormat: @"add (real function, real function)"];
					implementationCase = 5;
				} else if (!op1.isComplex && op2.isComplex) {
					graphvisDescription = [NSString stringWithFormat: @"add (real function, complex function)"];
					implementationCase = 6;
				} else if (op1.isComplex && !op2.isComplex) {
					graphvisDescription = [NSString stringWithFormat: @"add (complex function, real function)"];
					implementationCase = 7;
				} else {
					graphvisDescription = [NSString stringWithFormat: @"add (complex function, complex function)"];
					implementationCase = 5;
				}
			}
			
		} else if (op2.rank == 2) {
			GLLinearTransform *B = (GLLinearTransform *) op2;
			result = [GLLinearTransform transformOfType: format withFromDimensions: B.fromDimensions toDimensions:B.toDimensions inFormat:B.matrixFormats forEquation:B.equation matrix: nil];
			result.isPurelyReal = fOperand.isPurelyReal && sOperand.isPurelyReal;
			result.isPurelyImaginary = fOperand.isPurelyImaginary && sOperand.isPurelyImaginary;
			
			numPoints = result.nDataPoints;
			nDataElements = result.nDataElements;
			
			if (op1.rank == 0 ) {
				if ( !op1.isComplex && !op2.isComplex) {
					graphvisDescription = [NSString stringWithFormat: @"add (real scalar, real matrix)"];
					implementationCase = 1;
				} else if ( !op1.isComplex && op2.isComplex) {
					graphvisDescription = [NSString stringWithFormat: @"add (real scalar, complex matrix)"];
					implementationCase = 2;
				} else if ( op1.isComplex && !op2.isComplex) {
					graphvisDescription = [NSString stringWithFormat: @"add (complex scalar, real matrix)"];
					implementationCase = 3;
				} else {
					graphvisDescription = [NSString stringWithFormat: @"add (complex scalar, complex matrix)"];
					implementationCase = 4;
				}
			}
			else if (op1.rank == 1) {
				[NSException raise: @"TensorAdditionMismatch" format: @"Cannot add a function to a linear transformation"];
			} else if (op1.rank == 2) {
				GLLinearTransform *A = (GLLinearTransform *) op1;
				if ( ![A.fromDimensions isEqualToArray: B.fromDimensions] ) {
					[NSException raise: @"DimensionsNotEqualException" format: @"When adding two matrices, the fromDimensions of A, must equal the fromDimensions of B."];
				}
				if ( ![A.toDimensions isEqualToArray: B.toDimensions] ) {
					[NSException raise: @"DimensionsNotEqualException" format: @"When adding two matrices, the toDimensions of A, must equal the toDimensions of B."];
				}
				if ( ![A.matrixDescription isEqualToMatrixDescription: B.matrixDescription] ) {
					[NSException raise: @"UnsupportedMatrixFormatException" format: @"Cannot add two matrices in different formats using this operation."];
				}
				if ( !op1.isComplex && !op2.isComplex) {
					graphvisDescription = [NSString stringWithFormat: @"add (real function, real matrix)"];
					implementationCase = 5;
				} else if (!op1.isComplex && op2.isComplex) {
					graphvisDescription = [NSString stringWithFormat: @"add (real function, complex matrix)"];
					implementationCase = 6;
				} else if (op1.isComplex && !op2.isComplex) {
					graphvisDescription = [NSString stringWithFormat: @"add (complex function, real matrix)"];
					implementationCase = 7;
				} else {
					graphvisDescription = [NSString stringWithFormat: @"add (complex function, complex matrix)"];
					implementationCase = 5;
				}
			}
		}
		
		if (implementationCase == 1) {
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				vGL_vsadd( (void *) [operandArray[1] bytes], 1, (void *) [operandArray[0] bytes], (GLFloat *) [resultArray[0] bytes], 1, numPoints);
			};
		} else if (implementationCase == 2) {
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLSplitComplex rightComplex = splitComplexFromData( operandArray[1] );
				GLSplitComplex destComplex = splitComplexFromData( resultArray[0] );
				vGL_vsadd( rightComplex.realp, 1, (void *) [operandArray[0] bytes], destComplex.realp, 1, numPoints);
				vGL_mmov( rightComplex.imagp, destComplex.imagp, numPoints, 1, numPoints, numPoints);
			};
		} else if (implementationCase == 3) {
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloatComplex *a = (GLFloatComplex *) [operandArray[0] bytes];
				GLFloat a_realp = creal(*a);
				GLFloat a_imagp = cimag(*a);
				GLFloat *b = (GLFloat *) [operandArray[1] bytes];
				GLSplitComplex c = splitComplexFromData(resultArray[0]);
				vGL_vsadd( (void *) b, 1, &a_realp, (GLFloat *) c.realp, 1, numPoints);
				vGL_vfill( &a_imagp, (GLFloat *) c.imagp, 1, numPoints);
			};
		} else if (implementationCase == 4) {
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloatComplex *a = (GLFloatComplex *) [operandArray[0] bytes];
				GLFloat a_realp = creal(*a);
				GLFloat a_imagp = cimag(*a);
				GLSplitComplex b = splitComplexFromData(operandArray[1]);
				GLSplitComplex c = splitComplexFromData(resultArray[0]);
				vGL_vsadd( b.realp, 1, &a_realp, c.realp, 1, numPoints);
				vGL_vsadd( b.imagp, 1, &a_imagp, c.imagp, 1, numPoints);
			};
		} else if (implementationCase == 5) {
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				vGL_vadd( (void *) [operandArray[0] bytes], 1, (void *) [operandArray[1] bytes], 1, (GLFloat *) [resultArray[0] bytes], 1, nDataElements);
			};
		} else if (implementationCase == 6) {
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLSplitComplex leftComplex = splitComplexFromData( operandArray[0] );
				GLSplitComplex rightComplex = splitComplexFromData( operandArray[1] );
				GLSplitComplex destComplex = splitComplexFromData( resultArray[0] );
				
				vGL_vadd( leftComplex.realp, 1, rightComplex.realp, 1, destComplex.realp, 1, numPoints);
				vGL_mmov( rightComplex.imagp, destComplex.imagp, numPoints, 1, numPoints, numPoints);
			};
		} else if (implementationCase == 7) {
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLSplitComplex leftComplex = splitComplexFromData( operandArray[0] );
				GLSplitComplex rightComplex = splitComplexFromData( operandArray[1] );
				GLSplitComplex destComplex = splitComplexFromData( resultArray[0] );
				
				vGL_vadd( leftComplex.realp, 1, rightComplex.realp, 1, destComplex.realp, 1, numPoints);
				vGL_mmov( leftComplex.imagp, destComplex.imagp, numPoints, 1, numPoints, numPoints);
			};
		} else {
			[NSException raise: @"UnableToFindImplementation" format: @"Cannot find the implementation for this operation."];
		}
	}
    
	if (( self = [super initWithResult: @[result] operand: @[op1, op2] buffers: @[] operation: operation] )) {
		self.graphvisDescription = graphvisDescription;
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

@implementation GLSubtractionOperation

- (id) initWithFirstOperand: (GLTensor *) fOperand secondOperand: (GLTensor *) sOperand;
{
    // We order the operands so that scalars are always in the first position.
    // We can do this in this case because order doesn't matter for addition.
    GLDataFormat format = (fOperand.isPurelyReal && sOperand.isPurelyReal) ? kGLRealDataFormat : kGLSplitComplexDataFormat;
    GLTensor *op1 = (sOperand.rank < fOperand.rank) ? sOperand : fOperand;
    GLTensor *op2 = (sOperand.rank < fOperand.rank) ? fOperand : sOperand;
	BOOL didSwap = sOperand.rank < fOperand.rank;
    GLTensor *result;
	variableOperation operation;
	NSString *graphvisDescription;
	
	// Lots of cases to parse through.
	if (op2.rank == 0)
	{	// scalar-scalar
		result = [[GLScalar alloc] initWithType: format forEquation: op1.equation];
		result.isPurelyReal = fOperand.isPurelyReal && sOperand.isPurelyReal;
		result.isPurelyImaginary = fOperand.isPurelyImaginary && sOperand.isPurelyImaginary;
		
		if ( !op1.isComplex && !op2.isComplex) {
			graphvisDescription = [NSString stringWithFormat: @"subtract (real scalar, real scalar)"];
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *a = (GLFloat *) [operandArray[0] bytes];
				GLFloat *b = (GLFloat *) [operandArray[1] bytes];
				GLFloat *c = (GLFloat *) [resultArray[0] bytes];
				*c = (*a) - (*b);
			};
		} else if ( !op1.isComplex && op2.isComplex) {
			graphvisDescription = [NSString stringWithFormat: @"subtract (real scalar, complex scalar)"];
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *a = (GLFloat *) [operandArray[0] bytes];
				GLFloatComplex *b = (GLFloatComplex *) [operandArray[1] bytes];
				GLFloatComplex *c = (GLFloatComplex *) [resultArray[0] bytes];
				*c = (*a) - (*b);
			};
		} else if ( op1.isComplex && !op2.isComplex) {
			graphvisDescription = [NSString stringWithFormat: @"subtract (complex scalar, real scalar)"];
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloatComplex *a = (GLFloatComplex *) [operandArray[0] bytes];
				GLFloat *b = (GLFloat *) [operandArray[1] bytes];
				GLFloatComplex *c = (GLFloatComplex *) [resultArray[0] bytes];
				*c = (*a) - (*b);
			};
		} else {
			graphvisDescription = [NSString stringWithFormat: @"subtract (complex scalar, complex scalar)"];
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloatComplex *a = (GLFloatComplex *) [operandArray[0] bytes];
				GLFloatComplex *b = (GLFloatComplex *) [operandArray[1] bytes];
				GLFloatComplex *c = (GLFloatComplex *) [resultArray[0] bytes];
				*c = (*a) - (*b);
			};
		}
		
	}
	else if (op2.rank == 1 || op2.rank == 2)
	{
		// Adding functions and tensors share the same implementation details, so we divide the cases up here.
		NSUInteger implementationCase = 0;
		NSUInteger numPoints;
		NSUInteger nDataElements;
		
		if (op2.rank == 1)
		{	// scalar-function or function-function
			GLVariable *func2 = (GLVariable *) op2;
			result = [GLVariable variableOfType:format withDimensions: func2.dimensions forEquation: op2.equation];
			result.isPurelyReal = fOperand.isPurelyReal && sOperand.isPurelyReal;
			result.isPurelyImaginary = fOperand.isPurelyImaginary && sOperand.isPurelyImaginary;
			
			numPoints = result.nDataPoints;
			nDataElements = result.nDataElements;
			
			if (op1.rank == 0)
			{		// scalar-function
				if (!didSwap) {
					if ( !op1.isComplex && !op2.isComplex) {
						// C = a_real - B_real;
						graphvisDescription = [NSString stringWithFormat: @"subtract (real scalar, real function)"];
						implementationCase = 1;
					} else if ( !op1.isComplex && op2.isComplex) {
						// C = a_real - B_complex;
						graphvisDescription = [NSString stringWithFormat: @"subtract (real scalar, complex function)"];
						implementationCase = 2;
					} else if ( op1.isComplex && !op2.isComplex) {
						// C = a_complex - B_real
						graphvisDescription = [NSString stringWithFormat: @"subtract (complex scalar, real function)"];
						implementationCase = 3;
					} else {
						// C = a_complex - B_complex
						graphvisDescription = [NSString stringWithFormat: @"subtract (complex scalar, complex function)"];
						implementationCase = 4;
					}
				} else {
					if ( !op1.isComplex && !op2.isComplex) {
						// C = - a_real + B_real;
						graphvisDescription = [NSString stringWithFormat: @"subtract (real function, real scalar)"];
						implementationCase = 8;
					} else if ( !op1.isComplex && op2.isComplex) {
						// C = - a_real + B_complex;
						graphvisDescription = [NSString stringWithFormat: @"subtract (complex function, real scalar)"];
						implementationCase = 9;
					} else if ( op1.isComplex && !op2.isComplex) {
						// C = - a_complex + B_real
						graphvisDescription = [NSString stringWithFormat: @"subtract (real function, complex scalar)"];
						implementationCase = 10;
					} else {
						// C = - a_complex + B_complex
						graphvisDescription = [NSString stringWithFormat: @"subtract (complex function, complex scalar)"];
						implementationCase = 11;
					}
				}
			}
			else if (op1.rank == 1)
			{	// function-function
				GLVariable *func1 = (GLVariable *) op1;
				if ( ![func1.dimensions isEqualToArray: func2.dimensions] ) {
					[NSException raise: @"DimensionsNotEqualException" format: @"Cannot subtract two functions of different dimensions"];
				}
				
				if ( !op1.isComplex && !op2.isComplex) {
					graphvisDescription = [NSString stringWithFormat: @"subtract (real function, real function)"];
					implementationCase = 5;
				} else if (!op1.isComplex && op2.isComplex) {
					graphvisDescription = [NSString stringWithFormat: @"subtract (real function, complex function)"];
					implementationCase = 6;
				} else if (op1.isComplex && !op2.isComplex) {
					graphvisDescription = [NSString stringWithFormat: @"subtract (complex function, real function)"];
					implementationCase = 7;
				} else {
					graphvisDescription = [NSString stringWithFormat: @"subtract (complex function, complex function)"];
					implementationCase = 5;
				}
			}
			
		} else if (op2.rank == 2) {
			GLLinearTransform *B = (GLLinearTransform *) op2;
			result = [GLLinearTransform transformOfType: format withFromDimensions: B.fromDimensions toDimensions:B.toDimensions inFormat:B.matrixFormats forEquation:B.equation matrix: nil];
			result.isPurelyReal = fOperand.isPurelyReal && sOperand.isPurelyReal;
			result.isPurelyImaginary = fOperand.isPurelyImaginary && sOperand.isPurelyImaginary;
			
			numPoints = result.nDataPoints;
			nDataElements = result.nDataElements;
			
			if (op1.rank == 0 ) {
				if (!didSwap) {
					if ( !op1.isComplex && !op2.isComplex) {
						// C = a_real - B_real;
						graphvisDescription = [NSString stringWithFormat: @"subtract (real scalar, real matrix)"];
						implementationCase = 1;
					} else if ( !op1.isComplex && op2.isComplex) {
						// C = a_real - B_complex;
						graphvisDescription = [NSString stringWithFormat: @"subtract (real scalar, complex matrix)"];
						implementationCase = 2;
					} else if ( op1.isComplex && !op2.isComplex) {
						// C = a_complex - B_real
						graphvisDescription = [NSString stringWithFormat: @"subtract (complex scalar, real matrix)"];
						implementationCase = 3;
					} else {
						// C = a_complex - B_complex
						graphvisDescription = [NSString stringWithFormat: @"subtract (complex scalar, complex matrix)"];
						implementationCase = 4;
					}
				} else {
					if ( !op1.isComplex && !op2.isComplex) {
						// C = - a_real + B_real;
						graphvisDescription = [NSString stringWithFormat: @"subtract (real matrix, real scalar)"];
						implementationCase = 8;
					} else if ( !op1.isComplex && op2.isComplex) {
						// C = - a_real + B_complex;
						graphvisDescription = [NSString stringWithFormat: @"subtract (complex matrix, real scalar)"];
						implementationCase = 9;
					} else if ( op1.isComplex && !op2.isComplex) {
						// C = - a_complex + B_real
						graphvisDescription = [NSString stringWithFormat: @"subtract (real matrix, complex scalar)"];
						implementationCase = 10;
					} else {
						// C = - a_complex + B_complex
						graphvisDescription = [NSString stringWithFormat: @"subtract (complex matrix, complex scalar)"];
						implementationCase = 11;
					}
				}
			}
			else if (op1.rank == 1) {
				[NSException raise: @"TensorAdditionMismatch" format: @"Cannot subtract a function to a linear transformation"];
			} else if (op1.rank == 2) {
				GLLinearTransform *A = (GLLinearTransform *) op1;
				if ( ![A.fromDimensions isEqualToArray: B.fromDimensions] ) {
					[NSException raise: @"DimensionsNotEqualException" format: @"When subtracting two matrices, the fromDimensions of A, must equal the fromDimensions of B."];
				}
				if ( ![A.toDimensions isEqualToArray: B.toDimensions] ) {
					[NSException raise: @"DimensionsNotEqualException" format: @"When subtracting two matrices, the toDimensions of A, must equal the toDimensions of B."];
				}
				if ( ![A.matrixDescription isEqualToMatrixDescription: B.matrixDescription] ) {
					[NSException raise: @"UnsupportedMatrixFormatException" format: @"Cannot subtract two matrices in different formats using this operation."];
				}
				if ( !op1.isComplex && !op2.isComplex) {
					// C = A_real - B_real
					graphvisDescription = [NSString stringWithFormat: @"subtract (real matrix, real matrix)"];
					implementationCase = 5;
				} else if (!op1.isComplex && op2.isComplex) {
					// C = A_real - B_complex
					graphvisDescription = [NSString stringWithFormat: @"subtract (real matrix, complex matrix)"];
					implementationCase = 6;
				} else if (op1.isComplex && !op2.isComplex) {
					// C = A_complex - B_real
					graphvisDescription = [NSString stringWithFormat: @"subtract (complex matrix, real matrix)"];
					implementationCase = 7;
				} else {
					// C = A_complex - B_complex
					graphvisDescription = [NSString stringWithFormat: @"subtract (complex matrix, complex matrix)"];
					implementationCase = 5;
				}
			}
		}
		
		if (implementationCase == 1) {
			// C = a_real - B_real;
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *a = (GLFloat *) [operandArray[0] bytes];
				GLFloat *B = (GLFloat *) [operandArray[1] bytes];
				GLFloat *C = (GLFloat *) [resultArray[0] bytes];
				vGL_vneg( B, 1, C, 1, nDataElements );
				vGL_vsadd( C, 1, a, C, 1, nDataElements);
			};
		} else if (implementationCase == 2) {
			// C = a_real - B_complex;
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *a = (GLFloat *) [operandArray[0] bytes];
				GLSplitComplex B = splitComplexFromData( operandArray[1] );
				GLSplitComplex C = splitComplexFromData( resultArray[0] );
				vGL_vneg( B.realp, 1, C.realp, 1, nDataElements ); // we negate nDataElements, as in both the real and imaginary part
				vGL_vsadd( C.realp, 1, a, C.realp, 1, numPoints); // then add the real scalar to the real part
			};
		} else if (implementationCase == 3) {
			// C = a_complex - B_real
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloatComplex *a = (GLFloatComplex *) [operandArray[0] bytes];
				GLFloat a_realp = creal(*a);
				GLFloat a_imagp = cimag(*a);
				GLFloat *B = (GLFloat *) [operandArray[1] bytes];
				GLSplitComplex C = splitComplexFromData(resultArray[0]);
				vGL_vneg( B, 1, C.realp, 1, numPoints ); // we negate just the real part
				vGL_vsadd( C.realp, 1, &a_realp, C.realp, 1, numPoints); // add the real part of the scalar
				vGL_vfill( &a_imagp, C.imagp, 1, numPoints); // then copy/fill the imaginary part of the scalar
			};
		} else if (implementationCase == 4) {
			// C = a_complex - B_complex
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloatComplex *a = (GLFloatComplex *) [operandArray[0] bytes];
				GLFloat a_realp = creal(*a);
				GLFloat a_imagp = cimag(*a);
				GLSplitComplex B = splitComplexFromData(operandArray[1]);
				GLSplitComplex C = splitComplexFromData(resultArray[0]);
				vGL_vneg( B.realp, 1, C.realp, 1, nDataElements ); // we negate nDataElements, as in both the real and imaginary part
				vGL_vsadd( C.realp, 1, &a_realp, C.realp, 1, numPoints); // then add the real scalar to the real part
				vGL_vsadd( C.imagp, 1, &a_imagp, C.imagp, 1, numPoints); // then add the imaginary part of the scalar to the imaginary part of the vector
			};
		} else if (implementationCase == 5) {
			// C = A_real - B_real
			// C = A_complex - B_complex
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *A = (GLFloat *) [operandArray[0] bytes];
				GLFloat *B = (GLFloat *) [operandArray[1] bytes];
				GLFloat *C = (GLFloat *) [resultArray[0] bytes];
				vGL_vsub( B, 1, A, 1, C, 1, nDataElements); // Note that vsub does: C = B - A
			};
		} else if (implementationCase == 6) {
			// C = A_real - B_complex
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *A = (GLFloat *) [operandArray[0] bytes];
				GLSplitComplex B = splitComplexFromData( operandArray[1] );
				GLSplitComplex C = splitComplexFromData( resultArray[0] );
				
				vGL_vsub( B.realp, 1, A, 1, C.realp, 1, numPoints); // Note that vsub does: C = B - A.
				vGL_vneg( B.imagp, 1, C.imagp, 1, numPoints ); // we negate the imaginary part
			};
		} else if (implementationCase == 7) {
			// C = A_complex - B_real
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLSplitComplex A = splitComplexFromData( operandArray[0] );
				GLFloat *B = (GLFloat *) [operandArray[1] bytes];
				GLSplitComplex C = splitComplexFromData( resultArray[0] );
				
				vGL_vsub( B, 1, A.realp, 1, C.realp, 1, numPoints); // Note that vsub does: C = B - A.
				vGL_mmov( A.imagp, C.imagp, numPoints, 1, numPoints, numPoints); // we copy the imaginary part
			};
		} else if (implementationCase == 8) {
			// C = - a_real + B_real;
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *a = (GLFloat *) [operandArray[0] bytes];
				GLFloat a_neg = -(*a);
				GLFloat *B = (GLFloat *) [operandArray[1] bytes];
				GLFloat *C = (GLFloat *) [resultArray[0] bytes];
				vGL_vsadd( B, 1, &a_neg, C, 1, nDataElements);
			};
		} else if (implementationCase == 9) {
			// C = - a_real + B_complex;
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *a = (GLFloat *) [operandArray[0] bytes];
				GLFloat a_neg = -(*a);
				GLSplitComplex B = splitComplexFromData(operandArray[1]);
				GLSplitComplex C = splitComplexFromData(resultArray[0]);
				vGL_vsadd( B.realp, 1, &a_neg, C.realp, 1, numPoints); // add the negative real scalar to the real part of B
				vGL_mmov( B.imagp, C.imagp, numPoints, 1, numPoints, numPoints); // we copy the imaginary part
			};
		} else if (implementationCase == 10) {
			// C = - a_complex + B_real
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloatComplex *a = (GLFloatComplex *) [operandArray[0] bytes];
				GLFloat a_realp = -creal(*a);
				GLFloat a_imagp = -cimag(*a);
				GLFloat *B = (GLFloat *) [operandArray[1] bytes];
				GLSplitComplex C = splitComplexFromData(resultArray[0]);
				
				vGL_vsadd( B, 1, &a_realp, C.realp, 1, numPoints); // add the negative real scalar to the to B
				vGL_vfill( &a_imagp, C.imagp, 1, numPoints); // then copy/fill the imaginary part of the scalar
			};
		} else if (implementationCase == 11) {
			// C = - a_complex + B_complex
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloatComplex *a = (GLFloatComplex *) [operandArray[0] bytes];
				GLFloat a_realp = -creal(*a);
				GLFloat a_imagp = -cimag(*a);
				GLSplitComplex B = splitComplexFromData(operandArray[1]);
				GLSplitComplex C = splitComplexFromData(resultArray[0]);
				vGL_vsadd( B.realp, 1, &a_realp, C.realp, 1, numPoints); // add the negative real scalar to the to B
				vGL_vsadd( B.imagp, 1, &a_imagp, C.imagp, 1, numPoints); // add the negative imaginary scalar to the to imaginary B
			};
		}
		else {
			[NSException raise: @"UnableToFindImplementation" format: @"Cannot find the implementation for this operation."];
		}
	}
    
	if (( self = [super initWithResult: @[result] operand: @[op1, op2] buffers: @[] operation: operation] )) {
		self.graphvisDescription = graphvisDescription;
    }
    return self;
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

- (id) initWithFirstOperand: (GLTensor *) firstOperand secondOperand: (GLTensor *) secondOperand;
{
	/********************************************************************/
	/*		vector - lower dimensional vector multiplication			*/
	/********************************************************************/
	if (firstOperand.rank == 1 && secondOperand.rank == 1)
	{
		GLVariable *fOperand = (GLVariable *) firstOperand;
		GLVariable *sOperand = (GLVariable *) firstOperand;
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
			
			if (( self = [super initWithResult: @[result] operand: @[lowerDimVariable, higherDimVariable]] )) {
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
				
				if ( !lowerDimVariable.isComplex && !higherDimVariable.isComplex )
				{
					dispatch_queue_t globalQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
					NSUInteger totalLoops = outerLoopSize*innerLoopSize;
					self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
						const GLFloat *lowerF = [operandArray[0] bytes];
						const GLFloat *higherF = [operandArray[1] bytes];
						GLFloat *higherOut = [resultArray[0] mutableBytes];
						
						dispatch_apply(totalLoops, globalQueue, ^(size_t iteration) {
							NSUInteger index = (iteration/innerLoopSize)*outerLoopStride + (iteration%innerLoopSize)*innerLoopStride;
							vGL_vmul( lowerF, 1, &(higherF[index]), multiplicationStride, &(higherOut[index]), multiplicationStride, multiplicationLength);
						});
					};
					self.canOperateInPlace = YES;
					self.graphvisDescription = [NSString stringWithFormat: @"multiplication (real %lu dim, real %lu dim)", (unsigned long)lowerDimVariable.dimensions.count, (unsigned long)higherDimVariable.dimensions.count];
				}
				else if ( lowerDimVariable.isComplex && !higherDimVariable.isComplex )
				{
					dispatch_queue_t globalQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
					NSUInteger totalLoops = outerLoopSize*innerLoopSize;
					self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
						const GLSplitComplex lowerF = splitComplexFromData(operandArray[0]);
						const GLFloat *higherF = [operandArray[1] bytes];
						const GLSplitComplex higherOut = splitComplexFromData(resultArray[0]);
						dispatch_apply(totalLoops, globalQueue, ^(size_t iteration) {
							NSUInteger index = (iteration/innerLoopSize)*outerLoopStride + (iteration%innerLoopSize)*innerLoopStride;
							// (a + i b)*(x) = a x + i b x
							vGL_vmul( lowerF.realp, 1, &(higherF[index]), multiplicationStride, &(higherOut.realp[index]), multiplicationStride, multiplicationLength);
							vGL_vmul( lowerF.imagp, 1, &(higherF[index]), multiplicationStride, &(higherOut.imagp[index]), multiplicationStride, multiplicationLength);
						});
					};
					self.canOperateInPlace = YES;
					self.graphvisDescription = [NSString stringWithFormat: @"multiplication (complex %lu dim, real %lu dim)", (unsigned long)lowerDimVariable.dimensions.count, (unsigned long)higherDimVariable.dimensions.count];
				} else {
					[NSException raise: @"MethodNotImplemented" format: @"This case has not yet been implemented."];
				}

			}
			return self;
		}
	}
	
	// We order the operands so that scalars are always in the first position.
    // We can do this in this case because order doesn't matter for addition.
    GLDataFormat format = (firstOperand.isPurelyReal && secondOperand.isPurelyReal) ? kGLRealDataFormat : kGLSplitComplexDataFormat;
    GLTensor *op1 = (secondOperand.rank < firstOperand.rank) ? secondOperand : firstOperand;
    GLTensor *op2 = (secondOperand.rank < firstOperand.rank) ? firstOperand : secondOperand;
	BOOL didSwap = secondOperand.rank < firstOperand.rank;
    GLTensor *result;
	variableOperation operation;
	NSString *graphvisDescription;
	
	// Lots of cases to parse through.
	if (op2.rank == 0)
	{	// scalar-scalar
		result = [[GLScalar alloc] initWithType: format forEquation: op1.equation];
		result.isPurelyReal = (op1.isPurelyReal && op2.isPurelyReal) || (op1.isPurelyImaginary && op2.isPurelyImaginary);
		result.isPurelyImaginary = (op1.isPurelyReal && op2.isPurelyImaginary) || (op1.isPurelyImaginary && op2.isPurelyReal);
		
		if ( !op1.isComplex && !op2.isComplex) {
			graphvisDescription = [NSString stringWithFormat: @"multiply (real scalar, real scalar)"];
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *a = (GLFloat *) [operandArray[0] bytes];
				GLFloat *b = (GLFloat *) [operandArray[1] bytes];
				GLFloat *c = (GLFloat *) [resultArray[0] bytes];
				*c = (*a) * (*b);
			};
		} else if ( !op1.isComplex && op2.isComplex) {
			graphvisDescription = [NSString stringWithFormat: @"multiply (real scalar, complex scalar)"];
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *a = (GLFloat *) [operandArray[0] bytes];
				GLFloatComplex *b = (GLFloatComplex *) [operandArray[1] bytes];
				GLFloatComplex *c = (GLFloatComplex *) [resultArray[0] bytes];
				*c = (*a) * (*b);
			};
		} else if ( op1.isComplex && !op2.isComplex) {
			graphvisDescription = [NSString stringWithFormat: @"multiply (complex scalar, real scalar)"];
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloatComplex *a = (GLFloatComplex *) [operandArray[0] bytes];
				GLFloat *b = (GLFloat *) [operandArray[1] bytes];
				GLFloatComplex *c = (GLFloatComplex *) [resultArray[0] bytes];
				*c = (*a) * (*b);
			};
		} else {
			graphvisDescription = [NSString stringWithFormat: @"multiply (complex scalar, complex scalar)"];
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloatComplex *a = (GLFloatComplex *) [operandArray[0] bytes];
				GLFloatComplex *b = (GLFloatComplex *) [operandArray[1] bytes];
				GLFloatComplex *c = (GLFloatComplex *) [resultArray[0] bytes];
				*c = (*a) * (*b);
			};
		}
		
	}
	else if (op2.rank == 1 || op2.rank == 2)
	{
		// Adding functions and tensors share the same implementation details, so we divide the cases up here.
		NSUInteger implementationCase = 0;
		
		if (op2.rank == 1)
		{	// scalar-function or function-function
			GLVariable *func2 = (GLVariable *) op2;
			result = [GLVariable variableOfType:format withDimensions: func2.dimensions forEquation: op2.equation];
			result.isPurelyReal = (op1.isPurelyReal && op2.isPurelyReal) || (op1.isPurelyImaginary && op2.isPurelyImaginary);
			result.isPurelyImaginary = (op1.isPurelyReal && op2.isPurelyImaginary) || (op1.isPurelyImaginary && op2.isPurelyReal);
			
			if (op1.rank == 0)
			{		// scalar-function
				if ( !op1.isComplex && !op2.isComplex) {
					// C = a_real * B_real;
					graphvisDescription = [NSString stringWithFormat: @"multiply (real scalar, real function)"];
					implementationCase = 1;
				} else if ( !op1.isComplex && op2.isComplex) {
					// C = a_real * B_complex;
					graphvisDescription = [NSString stringWithFormat: @"multiply (real scalar, complex function)"];
					implementationCase = 1;
				} else if ( op1.isComplex && !op2.isComplex) {
					// C = a_complex * B_real
					graphvisDescription = [NSString stringWithFormat: @"multiply (complex scalar, real function)"];
					implementationCase = 3;
				} else {
					// C = a_complex * B_complex
					graphvisDescription = [NSString stringWithFormat: @"multiply (complex scalar, complex function)"];
					implementationCase = 4;
				}
			}
			else if (op1.rank == 1)
			{	// function-function
				GLVariable *func1 = (GLVariable *) op1;
				if ( ![func1.dimensions isEqualToArray: func2.dimensions] ) {
					[NSException raise: @"DimensionsNotEqualException" format: @"Cannot multiply two functions of different dimensions"];
				}
				
				if ( !op1.isComplex && !op2.isComplex) {
					// C = A_real * B_real
					graphvisDescription = [NSString stringWithFormat: @"multiply (real function, real function)"];
					implementationCase = 5;
				} else if (!op1.isComplex && op2.isComplex) {
					// C = A_real * B_complex
					graphvisDescription = [NSString stringWithFormat: @"multiply (real function, complex function)"];
					implementationCase = 6;
				} else if (op1.isComplex && !op2.isComplex) {
					// C = A_complex * B_real
					graphvisDescription = [NSString stringWithFormat: @"multiply (complex function, real function)"];
					implementationCase = 7;
				} else {
					if (op1.isPurelyReal && op2.isPurelyReal) {
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex real function, complex real function)"];
						implementationCase = 1;
					} else if (op1.isPurelyReal && op2.isPurelyImaginary) {
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex real function, complex imaginary function)"];
						implementationCase = 9;
					} else if (op1.isPurelyReal) {
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex real function, complex function)"];
						implementationCase = 10;
					} else if (op1.isPurelyImaginary && op2.isPurelyReal) {
						GLTensor *tmp = op1;
						op1 = op2;
						op2 = tmp;
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex imaginary function, complex real function)"];
						implementationCase = 9;
					} else if (op1.isPurelyImaginary && op2.isPurelyImaginary) {
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex imaginary function, complex imaginary function)"];
						implementationCase = 11;
					} else if (op1.isPurelyImaginary) {
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex imaginary function, complex function)"];
						implementationCase = 12;
					} else if (op2.isPurelyReal) {
						GLTensor *tmp = op1;
						op1 = op2;
						op2 = tmp;
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex real function, complex function)"];
						implementationCase = 10;
					} else if (op2.isPurelyImaginary) {
						GLTensor *tmp = op1;
						op1 = op2;
						op2 = tmp;
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex imaginary function, complex function)"];
						implementationCase = 11;
					} else {
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex function, complex function)"];
						implementationCase = 13;
					}
				}
			}
			
		} else if (op2.rank == 2) {
			GLLinearTransform *B = (GLLinearTransform *) op2;
			if (op1.rank == 0 ) {
				// C^i_j = a * B^i_j
				result = [GLLinearTransform transformOfType: format withFromDimensions: B.fromDimensions toDimensions:B.toDimensions inFormat:B.matrixFormats forEquation:B.equation matrix: nil];
				result.isPurelyReal = (op1.isPurelyReal && op2.isPurelyReal) || (op1.isPurelyImaginary && op2.isPurelyImaginary);
				result.isPurelyImaginary = (op1.isPurelyReal && op2.isPurelyImaginary) || (op1.isPurelyImaginary && op2.isPurelyReal);
	
				if ( !op1.isComplex && !op2.isComplex) {
					// C = a_real * B_real;
					graphvisDescription = [NSString stringWithFormat: @"multiply (real scalar, real matrix)"];
					implementationCase = 1;
				} else if ( !op1.isComplex && op2.isComplex) {
					// C = a_real * B_complex;
					graphvisDescription = [NSString stringWithFormat: @"multiply (real scalar, complex matrix)"];
					implementationCase = 1;
				} else if ( op1.isComplex && !op2.isComplex) {
					// C = a_complex * B_real
					graphvisDescription = [NSString stringWithFormat: @"multiply (complex scalar, real matrix)"];
					implementationCase = 3;
				} else {
					// C = a_complex * B_complex
					graphvisDescription = [NSString stringWithFormat: @"multiply (complex scalar, complex matrix)"];
					implementationCase = 4;
				}
			}
			else if (op1.rank == 1) {
				GLVariable *func1 = (GLVariable *) op1;
				if (!didSwap) {
					[NSException raise: @"TensorMultiplicationMismatch" format: @"Cannot left-multiply a function with a linear transformation"];
				}
				if ( ![B.fromDimensions isEqualToArray: func1.dimensions] ) {
					[NSException raise: @"DimensionsNotEqualException" format: @"When multiplying transforming a vector, the fromDimensions of A, must equal the dimensions of b."];
				}
				GLMatrixDescription *matrix = B.matrixDescription;
				for (NSUInteger i=0; i<matrix.nDimensions; i++) {
					if (matrix.strides[i].format != kGLDiagonalMatrixFormat) {
						[NSException raise: @"BadFormat" format: @"This operation type can only transform with matrices in a diagonal format."];
					}
				}
				
				result = [GLVariable variableOfType:format withDimensions: B.toDimensions forEquation: op2.equation];
				result.isPurelyReal = (op1.isPurelyReal && op2.isPurelyReal) || (op1.isPurelyImaginary && op2.isPurelyImaginary);
				result.isPurelyImaginary = (op1.isPurelyReal && op2.isPurelyImaginary) || (op1.isPurelyImaginary && op2.isPurelyReal);
				
				if ( !op1.isComplex && !op2.isComplex) {
					// C = A_real * B_real
					graphvisDescription = [NSString stringWithFormat: @"multiply (real matrix, real function)"];
					implementationCase = 5;
				} else if (!op1.isComplex && op2.isComplex) {
					// C = A_real * B_complex
					graphvisDescription = [NSString stringWithFormat: @"multiply (real matrix, complex function)"];
					implementationCase = 6;
				} else if (op1.isComplex && !op2.isComplex) {
					// C = A_complex * B_real
					graphvisDescription = [NSString stringWithFormat: @"multiply (complex matrix, real function)"];
					implementationCase = 7;
				} else {
					if (op1.isPurelyReal && op2.isPurelyReal) {
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex real matrix, complex real function)"];
						implementationCase = 1;
					} else if (op1.isPurelyReal && op2.isPurelyImaginary) {
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex real matrix, complex imaginary function)"];
						implementationCase = 9;
					} else if (op1.isPurelyReal) {
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex real matrix, complex function)"];
						implementationCase = 10;
					} else if (op1.isPurelyImaginary && op2.isPurelyReal) {
						GLTensor *tmp = op1;
						op1 = op2;
						op2 = tmp;
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex imaginary matrix, complex real function)"];
						implementationCase = 9;
					} else if (op1.isPurelyImaginary && op2.isPurelyImaginary) {
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex imaginary matrix, complex imaginary function)"];
						implementationCase = 11;
					} else if (op1.isPurelyImaginary) {
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex imaginary matrix, complex function)"];
						implementationCase = 12;
					} else if (op2.isPurelyReal) {
						GLTensor *tmp = op1;
						op1 = op2;
						op2 = tmp;
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex real matrix, complex function)"];
						implementationCase = 10;
					} else if (op2.isPurelyImaginary) {
						GLTensor *tmp = op1;
						op1 = op2;
						op2 = tmp;
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex imaginary matrix, complex function)"];
						implementationCase = 11;
					} else {
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex matrix, complex function)"];
						implementationCase = 13;
					}
				}
				
			} else if (op1.rank == 2) {
				// C^i_j = A^i_k * B^k_j
				GLLinearTransform *A = (GLLinearTransform *) op1;
				if ( ![A.fromDimensions isEqualToArray: B.toDimensions] ) {
					[NSException raise: @"DimensionsNotEqualException" format: @"When multiplying two matrices, the fromDimensions of A, must equal the toDimensions of B."];
				}
				if ( ![A.matrixDescription isEqualToMatrixDescription: B.matrixDescription] ) {
					[NSException raise: @"UnsupportedMatrixFormatException" format: @"Cannot multiply two matrices in different formats using this operation."];
				}
				GLMatrixDescription *matrix = A.matrixDescription;
				for (NSUInteger i=0; i<matrix.nDimensions; i++) {
					if (matrix.strides[i].format != kGLDiagonalMatrixFormat) {
						[NSException raise: @"BadFormat" format: @"This operation type can only transform with matrices in a diagonal format."];
					}
				}
				
				result = [GLLinearTransform transformOfType: format withFromDimensions: B.fromDimensions toDimensions:A.toDimensions inFormat:B.matrixFormats forEquation:B.equation matrix: nil];
				result.isPurelyReal = (op1.isPurelyReal && op2.isPurelyReal) || (op1.isPurelyImaginary && op2.isPurelyImaginary);
				result.isPurelyImaginary = (op1.isPurelyReal && op2.isPurelyImaginary) || (op1.isPurelyImaginary && op2.isPurelyReal);
				
				if ( !op1.isComplex && !op2.isComplex) {
					// C = A_real * B_real
					graphvisDescription = [NSString stringWithFormat: @"multiply (real matrix, real matrix)"];
					implementationCase = 5;
				} else if (!op1.isComplex && op2.isComplex) {
					// C = A_real * B_complex
					graphvisDescription = [NSString stringWithFormat: @"multiply (real matrix, complex matrix)"];
					implementationCase = 6;
				} else if (op1.isComplex && !op2.isComplex) {
					// C = A_complex * B_real
					graphvisDescription = [NSString stringWithFormat: @"multiply (complex matrix, real matrix)"];
					implementationCase = 7;
				} else {
					if (op1.isPurelyReal && op2.isPurelyReal) {
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex real matrix, complex real matrix)"];
						implementationCase = 1;
					} else if (op1.isPurelyReal && op2.isPurelyImaginary) {
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex real matrix, complex imaginary matrix)"];
						implementationCase = 9;
					} else if (op1.isPurelyReal) {
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex real matrix, complex matrix)"];
						implementationCase = 10;
					} else if (op1.isPurelyImaginary && op2.isPurelyReal) {
						GLTensor *tmp = op1;
						op1 = op2;
						op2 = tmp;
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex imaginary matrix, complex real matrix)"];
						implementationCase = 9;
					} else if (op1.isPurelyImaginary && op2.isPurelyImaginary) {
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex imaginary matrix, complex imaginary matrix)"];
						implementationCase = 11;
					} else if (op1.isPurelyImaginary) {
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex imaginary matrix, complex matrix)"];
						implementationCase = 12;
					} else if (op2.isPurelyReal) {
						GLTensor *tmp = op1;
						op1 = op2;
						op2 = tmp;
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex real matrix, complex matrix)"];
						implementationCase = 10;
					} else if (op2.isPurelyImaginary) {
						GLTensor *tmp = op1;
						op1 = op2;
						op2 = tmp;
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex imaginary matrix, complex matrix)"];
						implementationCase = 11;
					} else {
						graphvisDescription = [NSString stringWithFormat: @"multiply (complex matrix, complex matrix)"];
						implementationCase = 13;
					}
				}
			}
		}
		
		NSUInteger numPoints = result.nDataPoints;
		NSUInteger nDataElements = result.nDataElements;
		
		if (implementationCase == 1) {
			// C = a_real * B_real;		(a)*(x) = a x
			// C = a_real * B_complex;
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *a = (GLFloat *) [operandArray[0] bytes];
				GLFloat *B = (GLFloat *) [operandArray[1] bytes];
				GLFloat *C = (GLFloat *) [resultArray[0] bytes];
				vGL_vsmul( B, 1, a, C, 1, nDataElements);
			};
		} else if (implementationCase == 3) {
			// C = a_complex * B_real	(a + i b)*(x) = a x + i b x
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloatComplex *a = (GLFloatComplex *) [operandArray[0] bytes];
				GLFloat a_realp = creal(*a);
				GLFloat a_imagp = cimag(*a);
				GLFloat *B = (GLFloat *) [operandArray[1] bytes];
				GLSplitComplex C = splitComplexFromData(resultArray[0]);
				vGL_vsmul( B, 1, &a_realp, C.realp, 1, numPoints);
				vGL_vsmul( B, 1, &a_imagp, C.imagp, 1, numPoints);
			};
		} else if (implementationCase == 4) {
			// C = a_complex * B_complex	(a + i b)*(x + i y) = (a x - b y) + i (b x + ay)
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLSplitComplex a = splitComplexFromData(operandArray[0]); // This is a complex float, the implementation should be the same.
				GLSplitComplex B = splitComplexFromData(operandArray[1]);
				GLSplitComplex C = splitComplexFromData(resultArray[0]);
				vGL_zvzsml( &B, 1, &a, &C, 1, numPoints);
			};
		} else if (implementationCase == 5) {
			// C = A_real * B_real
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *A = (GLFloat *) [operandArray[0] bytes];
				GLFloat *B = (GLFloat *) [operandArray[1] bytes];
				GLFloat *C = (GLFloat *) [resultArray[0] bytes];
				vGL_vmul( A, 1, B, 1, C, 1, numPoints);
			};
		} else if (implementationCase == 6) {
			// C = A_real * B_complex
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *A = (GLFloat *) [operandArray[0] bytes];
				GLSplitComplex B = splitComplexFromData( operandArray[1] );
				GLSplitComplex C = splitComplexFromData( resultArray[0] );
				vGL_vmul( A, 1, B.realp, 1, C.realp, 1, numPoints);
				vGL_vmul( A, 1, B.imagp, 1, C.imagp, 1, numPoints);
			};
		} else if (implementationCase == 7) {
			// C = A_complex * B_real
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLSplitComplex A = splitComplexFromData( operandArray[0] );
				GLFloat *B = (GLFloat *) [operandArray[1] bytes];
				GLSplitComplex C = splitComplexFromData( resultArray[0] );
				vGL_vmul( A.realp, 1, B, 1, C.realp, 1, numPoints);
				vGL_vmul( A.imagp, 1, B, 1, C.imagp, 1, numPoints);
			};
		} else if (implementationCase == 9) {
			// (a)*(i y) = i a y
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLSplitComplex A = splitComplexFromData(operandArray[0]);
				GLSplitComplex B = splitComplexFromData(operandArray[1]);
				GLSplitComplex C = splitComplexFromData(resultArray[0]);
				vGL_vclr( C.realp, 1, numPoints);
				vGL_vmul( A.realp, 1, B.imagp, 1, C.imagp, 1, numPoints);
			};
		} else if (implementationCase == 10) {
			// a *(x + i y) = a x + i a y
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLSplitComplex A = splitComplexFromData(operandArray[0]);
				GLSplitComplex B = splitComplexFromData(operandArray[1]);
				GLSplitComplex C = splitComplexFromData(resultArray[0]);
				vGL_vmul( A.realp, 1, B.imagp, 1, C.imagp, 1, numPoints);
				vGL_vmul( A.realp, 1, B.imagp, 1, C.imagp, 1, numPoints);
			};
		} else if (implementationCase == 11) {
			// (i b)*(i y) = - b y
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLSplitComplex A = splitComplexFromData(operandArray[0]);
				GLSplitComplex B = splitComplexFromData(operandArray[1]);
				GLSplitComplex C = splitComplexFromData(resultArray[0]);
				vGL_vmul( A.imagp, 1, B.imagp, 1, C.realp, 1, numPoints);
				vGL_vneg( C.realp, 1, C.realp, 1, numPoints );
				vGL_vclr( C.imagp, 1, numPoints);
			};
		} else if (implementationCase == 12) {
			// i b *(x + i y) = - b y + i b x
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLSplitComplex A = splitComplexFromData(operandArray[0]);
				GLSplitComplex B = splitComplexFromData(operandArray[1]);
				GLSplitComplex C = splitComplexFromData(resultArray[0]);
				vGL_vmul( A.imagp, 1, B.imagp, 1, C.realp, 1, numPoints);
				vGL_vneg( C.realp, 1, C.realp, 1, numPoints );
				vGL_vmul( A.imagp, 1, B.realp, 1, C.imagp, 1, numPoints);
			};
		} else if (implementationCase == 13) {
			// (a + i b)*(x + i y) = (a x - b y) + i (b x + ay)
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLSplitComplex A = splitComplexFromData(operandArray[0]);
				GLSplitComplex B = splitComplexFromData(operandArray[1]);
				GLSplitComplex C = splitComplexFromData(resultArray[0]);
				vGL_zvmul(&A, 1, &B, 1, &C, 1, numPoints, 1);
			};
		}
		else {
			[NSException raise: @"UnableToFindImplementation" format: @"Cannot find the implementation for this operation."];
		}
	}
    
	result.name = [NSString stringWithFormat: @"%@_%@", secondOperand.name, firstOperand.name];
	if (( self = [super initWithResult: @[result] operand: @[op1, op2] buffers: @[] operation: operation] )) {
		self.graphvisDescription = graphvisDescription;
    }
    return self;
}

@end

/************************************************/
/*		GLDivisionOperation						*/
/************************************************/
// variable = leftVariable / rightVariable
@implementation GLDivisionOperation

- (id) initWithFirstOperand: (GLTensor *) A secondOperand: (GLTensor *) B {
	return [self initWithFirstOperand: A secondOperand: B shouldUseComplexArithmetic: YES];
}

- (id) initWithFirstOperand: (GLTensor *) op1 secondOperand: (GLTensor *) op2 shouldUseComplexArithmetic: (BOOL) useComplexArithmetic
{
	GLDataFormat format = (op1.isPurelyReal && op2.isPurelyReal) ? kGLRealDataFormat : kGLSplitComplexDataFormat;
	GLTensor *result;
	variableOperation operation;
	NSString *graphvisDescription;
	
	if (op1.rank == 0 && op2.rank == 0)
	{	// c = a / b
		result = [[GLScalar alloc] initWithType: format forEquation: op1.equation];
		result.isPurelyReal = (op1.isPurelyReal && op2.isPurelyReal) || (op1.isPurelyImaginary && op2.isPurelyImaginary);
		result.isPurelyImaginary = (op1.isPurelyReal && op2.isPurelyImaginary) || (op1.isPurelyImaginary && op2.isPurelyReal);
		
		if ( !op1.isComplex && !op2.isComplex) {
			graphvisDescription = [NSString stringWithFormat: @"division (real scalar, real scalar)"];
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *a = (GLFloat *) [operandArray[0] bytes];
				GLFloat *b = (GLFloat *) [operandArray[1] bytes];
				GLFloat *c = (GLFloat *) [resultArray[0] bytes];
				*c = (*a) / (*b);
			};
		} else if ( !op1.isComplex && op2.isComplex) {
			graphvisDescription = [NSString stringWithFormat: @"division (real scalar, complex scalar)"];
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *a = (GLFloat *) [operandArray[0] bytes];
				GLFloatComplex *b = (GLFloatComplex *) [operandArray[1] bytes];
				GLFloatComplex *c = (GLFloatComplex *) [resultArray[0] bytes];
				*c = (*a) / (*b);
			};
		} else if ( op1.isComplex && !op2.isComplex) {
			graphvisDescription = [NSString stringWithFormat: @"division (complex scalar, real scalar)"];
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloatComplex *a = (GLFloatComplex *) [operandArray[0] bytes];
				GLFloat *b = (GLFloat *) [operandArray[1] bytes];
				GLFloatComplex *c = (GLFloatComplex *) [resultArray[0] bytes];
				*c = (*a) / (*b);
			};
		} else {
			graphvisDescription = [NSString stringWithFormat: @"division (complex scalar, complex scalar)"];
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloatComplex *a = (GLFloatComplex *) [operandArray[0] bytes];
				GLFloatComplex *b = (GLFloatComplex *) [operandArray[1] bytes];
				GLFloatComplex *c = (GLFloatComplex *) [resultArray[0] bytes];
				*c = (*a) / (*b);
			};
		}
	}
	else if (op1.rank == 0 && op2.rank == 1)
	{	// C^i = a / B^i
		GLVariable *func2 = (GLVariable *) op2;
		result = [GLVariable variableOfType:format withDimensions: func2.dimensions forEquation: op2.equation];
		result.isPurelyReal = (op1.isPurelyReal && op2.isPurelyReal) || (op1.isPurelyImaginary && op2.isPurelyImaginary);
		result.isPurelyImaginary = (op1.isPurelyReal && op2.isPurelyImaginary) || (op1.isPurelyImaginary && op2.isPurelyReal);
		
		if ( !op1.isComplex && !op2.isComplex) {
			// C = a_real / B_real;
			graphvisDescription = [NSString stringWithFormat: @"division (real scalar, real function)"];
		} else if ( !op1.isComplex && op2.isComplex) {
			// C = a_real / B_complex;
			graphvisDescription = [NSString stringWithFormat: @"division (real scalar, complex function)"];
		} else if ( op1.isComplex && !op2.isComplex) {
			// C = a_complex / B_real
			graphvisDescription = [NSString stringWithFormat: @"division (complex scalar, real function)"];
		} else {
			// C = a_complex / B_complex
			graphvisDescription = [NSString stringWithFormat: @"division (complex scalar, complex function)"];
		}
	}
	else if (op1.rank == 0 && op2.rank == 1)
	{	// C^i = A^i / b
		// This is just multiplication by a scalar 1/b.
		
	}
	else if (op1.rank == 1 && op2.rank == 1)
	{	// C^i = A^i / B^i
		GLVariable *func1 = (GLVariable *) op1;
		GLVariable *func2 = (GLVariable *) op2;
		
		if ( ![func1.dimensions isEqualToArray: func2.dimensions] ) {
			[NSException raise: @"DimensionsNotEqualException" format: @"Cannot add two functions of different dimensions"];
		}
		NSUInteger nDataPoints = op1.nDataPoints;
		if ( !op1.isComplex && !op2.isComplex) {
			// C = A_real / B_real;
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *A = (GLFloat *) [operandArray[0] bytes];
				GLFloat *B = (GLFloat *) [operandArray[1] bytes];
				GLFloat *C = (GLFloat *) [resultArray[0] bytes];
				vGL_vdiv( B, 1, A, 1, C, 1, nDataPoints); // Note that vdiv does: C = B / A
			};
			graphvisDescription = [NSString stringWithFormat: @"division (real function, real function)"];
		} else if ( !op1.isComplex && op2.isComplex) {
			// C = A_real / B_complex;
			graphvisDescription = [NSString stringWithFormat: @"division (real function, complex function)"];
		} else if ( op1.isComplex && !op2.isComplex) {
			// C = A_complex / B_real
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLSplitComplex A = splitComplexFromData([operandArray[0] bytes]);
				GLFloat *B = (GLFloat *) [operandArray[1] bytes];
				GLSplitComplex C = splitComplexFromData([resultArray[0] bytes]);
				vGL_vdiv( B, 1, A.realp, 1, C.realp, 1, nDataPoints); // Note that vdiv does: C = B / A
				vGL_vdiv( B, 1, A.imagp, 1, C.imagp, 1, nDataPoints); // Note that vdiv does: C = B / A
			};
			graphvisDescription = [NSString stringWithFormat: @"division (complex function, real function)"];
		} else {
			// C = A_complex / B_complex
			if (useComplexArithmetic) {
				operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
					GLSplitComplex A = splitComplexFromData([operandArray[0] bytes]);
					GLSplitComplex B = splitComplexFromData([operandArray[1] bytes]);
					GLSplitComplex C = splitComplexFromData([resultArray[0] bytes]);
					vGL_zvdiv( &B, 1, &A, 1, &C, 1, nDataPoints); // Note that vdiv does: C = B / A
				};
				graphvisDescription = [NSString stringWithFormat: @"division (complex function, complex function)"];
			} else {
				operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
					GLSplitComplex A = splitComplexFromData([operandArray[0] bytes]);
					GLSplitComplex B = splitComplexFromData([operandArray[1] bytes]);
					GLSplitComplex C = splitComplexFromData([resultArray[0] bytes]);
					vGL_vdiv( B.realp, 1, A.realp, 1, C.realp, 1, nDataPoints);
					vGL_vdiv( B.imagp, 1, A.imagp, 1, C.imagp, 1, nDataPoints);
				};
				graphvisDescription = [NSString stringWithFormat: @"element-wise division (complex function, complex function)"];
			}
		}
	}
	
	if (!operation) {
		[NSException raise: @"UnableToFindImplementation" format: @"Cannot find the implementation for this operation."];
	}
	
	if (( self = [super initWithResult: @[result] operand: @[op1, op2] buffers: @[] operation: operation] )) {
		self.graphvisDescription = graphvisDescription;
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
/*		GLAbsoluteLargestOperation				*/
/************************************************/
// variable = max( abs(leftVariable), abs(rightVariable ) element-wise

@implementation GLAbsoluteLargestOperation

- (id) initWithFirstOperand: (GLTensor *) fOperand secondOperand: (GLTensor *) sOperand;
{
	if (sOperand.rank != fOperand.rank) {
		[NSException raise: @"RankMismatch" format: @"Both tensors must be of the same rank."];
	}
	
	if (fOperand.isComplex != sOperand.isComplex) {
		[NSException raise: @"FormatMismatch" format: @"Both tensors must be in the same format."];
	}
	
    // We order the operands so that scalars are always in the first position.
    // We can do this in this case because order doesn't matter for addition.
    GLTensor *op1 = (sOperand.rank < fOperand.rank) ? sOperand : fOperand;
    GLTensor *op2 = (sOperand.rank < fOperand.rank) ? fOperand : sOperand;
    GLTensor *result;
	NSUInteger nDataElements = op1.nDataElements;
	NSString *graphvisDescription;
	
	if (fOperand.rank == 0)
	{
		result = [[GLScalar alloc] initWithType: fOperand.dataFormat forEquation: op1.equation];
		graphvisDescription = @"element-wise max (scalar)";
	}
	else if (fOperand.rank == 1)
	{
		GLVariable *func1 = (GLVariable *) fOperand;
		GLVariable *func2 = (GLVariable *) sOperand;
		if ( ![func1.dimensions isEqualToArray: func2.dimensions] ) {
			[NSException raise: @"DimensionsNotEqualException" format: @"Cannot compare two functions of different dimensions"];
		}
		result = [[GLVariable alloc] initVariableOfType: func1.dataFormat withDimensions: func1.dimensions forEquation:func1.equation];
		graphvisDescription = @"element-wise max (function)";
	}
	else if (fOperand.rank == 2)
	{
		GLLinearTransform *A = (GLLinearTransform *) fOperand;
		GLLinearTransform *B = (GLLinearTransform *) sOperand;
		if ( ![A.fromDimensions isEqualToArray: B.fromDimensions] ) {
			[NSException raise: @"DimensionsNotEqualException" format: @"When comparing two matrices, the fromDimensions of A, must equal the fromDimensions of B."];
		}
		if ( ![A.toDimensions isEqualToArray: B.toDimensions] ) {
			[NSException raise: @"DimensionsNotEqualException" format: @"When comparing two matrices, the toDimensions of A, must equal the toDimensions of B."];
		}
		if ( ![A.matrixDescription isEqualToMatrixDescription: B.matrixDescription] ) {
			[NSException raise: @"UnsupportedMatrixFormatException" format: @"Cannot compare two matrices in different formats using this operation."];
		}
		result = [GLLinearTransform transformOfType: B.dataFormat withFromDimensions: B.fromDimensions toDimensions:B.toDimensions inFormat:B.matrixFormats forEquation:B.equation matrix: nil];
		graphvisDescription = @"element-wise max (matrix)";
	}
	result.isPurelyReal = fOperand.isPurelyReal && sOperand.isPurelyReal;
	result.isPurelyImaginary = fOperand.isPurelyImaginary && sOperand.isPurelyImaginary;
	
	variableOperation operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
		GLFloat *A = (GLFloat *) [operandArray[0] bytes];
		GLFloat *B = (GLFloat *) [operandArray[1] bytes];
		GLFloat *C = (GLFloat *) [resultArray[0] bytes];
		vGL_vmaxmg( A, 1, B, 1, C, 1, nDataElements);
	};
	
	if (( self = [super initWithResult: @[result] operand: @[op1, op2] buffers: @[] operation: operation] )) {
		self.graphvisDescription = graphvisDescription;
    }
    return self;
}

- (BOOL) canOperateInPlace {
	return YES;
}

@end

/************************************************/
/*		GLDotProductOperation					*/
/************************************************/

@implementation GLDotProductOperation

- (id) initWithFirstOperand: (GLVariable *) fOperand secondOperand: (GLVariable *) sOperand {
	
	if ( ![fOperand.dimensions isEqualToArray: sOperand.dimensions] ) {
		[NSException raise: @"DimensionsNotEqualException" format: @"Cannot dot two functions of different dimensions"];
	}
	
	GLDataFormat format = fOperand.isComplex || sOperand.isComplex ? kGLSplitComplexDataFormat : kGLRealDataFormat;
	GLScalar *result = [[GLScalar alloc] initWithType: format forEquation: fOperand.equation];
	
	GLVariable *op1 = (GLVariable *) fOperand;
	GLVariable *op2 = (GLVariable *) sOperand;
	variableOperation operation;
	NSString *graphvisDescription;
	NSUInteger nDataPoints = op1.nDataPoints;
	
	if ( !op1.isComplex && !op2.isComplex) {
		// C = A_real • B_real;
		operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			GLFloat *A = (GLFloat *) [operandArray[0] bytes];
			GLFloat *B = (GLFloat *) [operandArray[1] bytes];
			GLFloat *C = (GLFloat *) [resultArray[0] bytes];
			vGL_dotpr( A, 1, B, 1, C, nDataPoints);
		};
		graphvisDescription = [NSString stringWithFormat: @"dot (real function, real function)"];
	} else if ( !op1.isComplex && op2.isComplex) {
		// C = A_real • B_complex;
		operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			GLFloat *A = (GLFloat *) [operandArray[0] bytes];
			GLSplitComplex B = splitComplexFromData([operandArray[1] bytes]);
			GLFloatComplex *c = (GLFloatComplex *) [resultArray[0] bytes];
			GLFloat *c_realp = (GLFloat *) c;
			vGL_dotpr( A, 1, B.realp, 1, c_realp, nDataPoints);
		};
		graphvisDescription = [NSString stringWithFormat: @"dot (real function, complex function)"];
	} else if ( op1.isComplex && !op2.isComplex) {
		// C = A_complex / B_real
		operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			GLSplitComplex A = splitComplexFromData([operandArray[0] bytes]);
			GLFloat *B = (GLFloat *) [operandArray[1] bytes];
			GLFloatComplex *c = (GLFloatComplex *) [resultArray[0] bytes];
			GLFloat *c_realp = (GLFloat *) c;
			vGL_dotpr( A.realp, 1, B, 1, c_realp, nDataPoints);
		};
		graphvisDescription = [NSString stringWithFormat: @"dot (complex function, real function)"];
	} else {
		operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			GLSplitComplex A = splitComplexFromData([operandArray[0] bytes]);
			GLSplitComplex B = splitComplexFromData([operandArray[1] bytes]);
			GLSplitComplex C = splitComplexFromData([resultArray[0] bytes]);
			vGL_zdotpr( &A, 1, &B, 1, &C, nDataPoints );
		};
		graphvisDescription = [NSString stringWithFormat: @"dot (complex function, complex function)"];
	}
	
	if (( self = [super initWithResult: @[result] operand: @[op1, op2] buffers: @[] operation: operation] )) {
		self.graphvisDescription = graphvisDescription;
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
	GLVariable *result= [GLVariable variableOfType: format withDimensions: variable.dimensions forEquation: variable.equation];
	
	if (( self = [super initWithResult: @[result] operand: @[variable, aScalarVariable]]))
	{
        self.indexString = indexString;
        
        NSUInteger numBytes = result.nDataElements*sizeof(GLFloat);
        
        result.name = variable.name;
		result.isPurelyReal = variable.isPurelyReal;
		result.isPurelyImaginary = variable.isPurelyImaginary;
        
        if ( variable.dimensions.count == 1 )
		{
            NSRange fastRange = [ranges[0] rangeValue];
            
            NSUInteger startIndex = fastRange.location;
            NSUInteger fastIndexLength = fastRange.length;
            
            self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *A = (GLFloat *) [operandArray[0] bytes];
				GLFloat *b = (GLFloat *) [operandArray[1] bytes];
                GLFloat *toData = (GLFloat *) [resultArray[0] bytes];
                // first copy the data
                memcpy( toData, A,  numBytes );
                // then replace the value at the desired indices
                vGL_vfill( b, &toData[startIndex], 1, fastIndexLength);
            };
            self.graphvisDescription = [NSString stringWithFormat: @"set leftVar=rightVar (1 dim)"];
		}
		else if ( variable.dimensions.count == 2 )
		{
			NSRange fastRange = [ranges[1] rangeValue];
            NSUInteger fastDimLength = [variable.dimensions[1] nPoints];
            
            NSRange slowRange = [ranges[0] rangeValue];
			
            self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *A = (GLFloat *) [operandArray[0] bytes];
				GLFloat *b = (GLFloat *) [operandArray[1] bytes];
                GLFloat *toData = (GLFloat *) [resultArray[0] bytes];
                // first copy the data
                memcpy( toData, A,  numBytes );
				
                dispatch_apply(slowRange.length, dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0), ^(size_t iteration) {
                    // then replace the value at the desired indices
                    vGL_vfill( b, &(toData[(slowRange.location + iteration)*fastDimLength + fastRange.location]), 1, fastRange.length);
                    
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
            
            self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *A = (GLFloat *) [operandArray[0] bytes];
				GLFloat *aScalar = (GLFloat *) [operandArray[1] bytes];
                GLFloat *toData = (GLFloat *) [resultArray[0] bytes];
				// first copy the data
                memcpy( toData, A,  numBytes );
				
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


