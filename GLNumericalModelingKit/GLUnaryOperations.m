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
#import "GLLinearTransform.h"

/************************************************/
/*		GLNegationOperation                     */
/************************************************/

@implementation GLNegationOperation

- (GLNegationOperation *) initWithVariable: (GLVariable *) variable
{
	if (( self = [super initWithOperand: @[variable] ]))
	{
		GLFunction *resultVariable = self.result[0];
		GLFunction *operandVariable = self.operand[0];
		
		NSUInteger numElements = resultVariable.nDataElements;
		resultVariable.isPurelyReal = operandVariable.isPurelyReal;
		resultVariable.isPurelyImaginary = operandVariable.isPurelyImaginary;
		self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			NSMutableData *result = resultArray[0];
			NSMutableData *operand = operandArray[0];
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


- (GLAbsoluteValueOperation *) initWithVariable: (GLVariable *) variable
{
	GLVariable *resultVar;
	
    // We have to override the initializer in this case because the result *will* be real.
    if (variable.rank == 0) {
        GLScalar *scalar = (GLScalar *) variable;
        resultVar=[[GLScalar alloc] initWithType: kGLRealDataFormat forEquation:scalar.equation];
    } else if (variable.rank == 1) {
        GLFunction *function = (GLFunction *) variable;
        resultVar=[[function class] functionOfType: kGLRealDataFormat withDimensions: function.dimensions forEquation: function.equation];
    }  else if (variable.rank == 2) {
        GLLinearTransform *matrix = (GLLinearTransform *) variable;
        resultVar=[GLLinearTransform transformOfType: kGLRealDataFormat withFromDimensions: matrix.fromDimensions toDimensions: matrix.toDimensions inFormat: matrix.matrixFormats forEquation:matrix.equation matrix:nil];
    } else {
		return nil;
	}
    
	if (( self = [super initWithResult: @[resultVar] operand: @[variable] ] ))
	{
		GLFunction *resultVariable = self.result[0];
		GLFunction *operandVariable = self.operand[0];
		
		resultVariable.isPurelyReal = YES;
		
		if (operandVariable.isComplex) {
			NSUInteger numPoints = resultVariable.nDataPoints;
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				NSMutableData *result = resultArray[0];
				GLSplitComplex fromSplit = splitComplexFromData(operandArray[0]);
				vGL_zvabs( &fromSplit, 1, result.mutableBytes, 1, numPoints );
			};
            self.graphvisDescription = @"abs (complex)";
		} else {
			NSUInteger numElements = resultVariable.nDataElements;
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				NSMutableData *result = resultArray[0];
				NSMutableData *operand = operandArray[0];
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

@end

/************************************************/
/*		GLRealPartOperation						*/
/************************************************/

@implementation GLRealPartOperation

#warning this still needs to test for the obvious case (that the variable is already real).
- (GLRealPartOperation *) initWithVariable: (GLVariable *) variable
{
	GLVariable *resultVar;
	
	// We have to override the initializer in this case because the result *will* be real.
	if (variable.rank == 0) {
		GLScalar *scalar = (GLScalar *) variable;
		resultVar=[[GLScalar alloc] initWithType: kGLRealDataFormat forEquation:scalar.equation];
	} else if (variable.rank == 1) {
		GLFunction *function = (GLFunction *) variable;
		resultVar=[[function class] functionOfType: kGLRealDataFormat withDimensions: function.dimensions forEquation: function.equation];
	}  else if (variable.rank == 2) {
		GLLinearTransform *matrix = (GLLinearTransform *) variable;
		resultVar=[GLLinearTransform transformOfType: kGLRealDataFormat withFromDimensions: matrix.fromDimensions toDimensions: matrix.toDimensions inFormat: matrix.matrixFormats forEquation:matrix.equation matrix:nil];
	} else {
		return nil;
	}
	
	if (( self = [super initWithResult: @[resultVar] operand: @[variable] ] ))
	{
		GLFunction *resultVariable = self.result[0];
		GLFunction *operandVariable = self.operand[0];
		
		resultVariable.isPurelyReal = YES;
		
		if (operandVariable.isComplex) {
			NSUInteger numPoints = resultVariable.nDataPoints;
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				NSMutableData *result = resultArray[0];
				GLSplitComplex fromSplit = splitComplexFromData(operandArray[0]);
				memcpy( result.mutableBytes, fromSplit.realp,  numPoints*sizeof(GLFloat));
			};
			self.graphvisDescription = @"real (complex)";
		} else {
			NSUInteger numElements = resultVariable.nDataElements;
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				NSMutableData *result = resultArray[0];
				NSMutableData *operand = operandArray[0];
				memcpy( result.mutableBytes, operand.mutableBytes,  numElements*sizeof(GLFloat));
			};
			self.graphvisDescription = @"real (real)";
		}
	}
	
	return self;
}
@end

/************************************************/
/*		GLExponentialOperation                  */
/************************************************/

@implementation GLExponentialOperation

- (GLExponentialOperation *) initWithVariable: (GLVariable *) variable
{
	// If the operand is purely real, we don't need a complex number
	GLDataFormat format = variable.isPurelyReal ? kGLRealDataFormat : kGLSplitComplexDataFormat;
	GLVariable *resultVar;
    
    if (variable.rank == 0) {
        GLScalar *scalar = (GLScalar *) variable;
        resultVar=[[GLScalar alloc] initWithType: format forEquation:scalar.equation];
    } else if (variable.rank == 1) {
        GLFunction *function = (GLFunction *) variable;
        resultVar=[[function class] functionOfType: format withDimensions: function.dimensions forEquation: function.equation];
    }  else if (variable.rank == 2) {
        [NSException raise: @"NotYetImplemented" format: @"Exponentiation is not yet implemented for linear transformations. This should be trivial for diagonal matrices"];
    }
    
	if (( self = [super initWithResult: @[resultVar] operand: @[variable] ] ))
	{
		GLFunction *resultVariable = self.result[0];
		GLFunction *operandVariable = self.operand[0];
		
		const int numElements = (int) resultVariable.nDataElements;
		const int numPoints = (int) resultVariable.nDataPoints;
		resultVariable.isPurelyReal = operandVariable.isPurelyReal;
				
		if (operandVariable.isComplex) {
			if (operandVariable.isPurelyReal) {
				self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
					NSMutableData *result = resultArray[0];
					NSMutableData *operand = operandArray[0];
					vGL_vvexp( result.mutableBytes, operand.bytes, &numElements );
				};
                self.graphvisDescription = @"exp (complex, purely real)";
			} else if (operandVariable.isPurelyImaginary) {
				self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
					GLSplitComplex fromSplit = splitComplexFromData(operandArray[0]);
					GLSplitComplex toSplit = splitComplexFromData(resultArray[0]);
					
					vGL_vvsincos( toSplit.imagp, toSplit.realp, fromSplit.imagp, &numPoints);
				};
                self.graphvisDescription = @"exp (complex, purely imaginary)";
			} else {
				self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
					GLSplitComplex fromSplit = splitComplexFromData(operandArray[0]);
					GLSplitComplex toSplit = splitComplexFromData(resultArray[0]);
                    GLFloat *buffer = (GLFloat *) [bufferArray[0] mutableBytes];
					
					vGL_vvsincos( toSplit.imagp, toSplit.realp, fromSplit.imagp, &numPoints);
					vGL_vvexp( buffer, fromSplit.realp, &numPoints );
					
					vGL_vmul( toSplit.realp, 1, buffer, 1, toSplit.realp, 1, numPoints);
					vGL_vmul( toSplit.imagp, 1, buffer, 1, toSplit.imagp, 1, numPoints);
				};
                self.buffer = @[ [[GLBuffer alloc] initWithLength: numPoints*sizeof(GLFloat)] ];
                self.graphvisDescription = @"exp (complex)";
			}
		} else {
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				NSMutableData *result = resultArray[0];
				NSMutableData *operand = operandArray[0];
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
/*		GLLogarithmOperation					*/
/************************************************/

@implementation GLLogarithmOperation

- (GLLogarithmOperation *) initWithVariable: (GLVariable *) variable
{
    if (!variable.isPurelyReal) {
        [NSException raise: @"NotYetImplemented" format: @"You can only take the sine of real numbers."];
    }
    
	GLVariable *resultVar;
    
    if (variable.rank == 0) {
        GLScalar *scalar = (GLScalar *) variable;
        resultVar=[[GLScalar alloc] initWithType: kGLRealDataFormat forEquation:scalar.equation];
    } else if (variable.rank == 1) {
        GLFunction *function = (GLFunction *) variable;
        resultVar=[[function class] functionOfType: kGLRealDataFormat withDimensions: function.dimensions forEquation: function.equation];
    }  else if (variable.rank == 2) {
        [NSException raise: @"BadMath" format: @"I don't know how to take the sine of a linear transformation."];
    }
    
	if (( self = [super initWithResult: @[resultVar] operand: @[variable] ]))
	{
		GLFunction *resultVariable = self.result[0];
		GLFunction *operandVariable = self.operand[0];
		
		const int numElements = (int) resultVariable.nDataElements;
		resultVariable.isPurelyReal = operandVariable.isPurelyReal;
		resultVariable.isPurelyImaginary = operandVariable.isPurelyImaginary;
		self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			NSMutableData *result = resultArray[0];
			NSMutableData *operand = operandArray[0];
			vGL_vvlog( result.mutableBytes, operand.bytes, &numElements );
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
/*		GLSineOperation							*/
/************************************************/

@implementation GLSineOperation

- (GLSineOperation *) initWithVariable: (GLVariable *) variable
{
    if (!variable.isPurelyReal) {
        [NSException raise: @"NotYetImplemented" format: @"You can only take the sine of real numbers."];
    }

	GLVariable *resultVar;
    
    if (variable.rank == 0) {
        GLScalar *scalar = (GLScalar *) variable;
        resultVar=[[GLScalar alloc] initWithType: kGLRealDataFormat forEquation:scalar.equation];
    } else if (variable.rank == 1) {
        GLFunction *function = (GLFunction *) variable;
        resultVar=[[function class] functionOfType: kGLRealDataFormat withDimensions: function.dimensions forEquation: function.equation];
    }  else if (variable.rank == 2) {
        [NSException raise: @"BadMath" format: @"I don't know how to take the sine of a linear transformation."];
    }
    
	if (( self = [super initWithResult: @[resultVar] operand: @[variable] ]))
	{
		GLFunction *resultVariable = self.result[0];
		GLFunction *operandVariable = self.operand[0];
		
		const int numElements = (int) resultVariable.nDataElements;
		resultVariable.isPurelyReal = operandVariable.isPurelyReal;
		resultVariable.isPurelyImaginary = operandVariable.isPurelyImaginary;
		self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			NSMutableData *result = resultArray[0];
			NSMutableData *operand = operandArray[0];
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

- (GLCosineOperation *) initWithVariable: (GLVariable *) variable
{
    if (!variable.isPurelyReal) {
        [NSException raise: @"NotYetImplemented" format: @"You can only take the cosine of real numbers."];
    }
    
	GLVariable *resultVar;
    
    if (variable.rank == 0) {
        GLScalar *scalar = (GLScalar *) variable;
        resultVar=[[GLScalar alloc] initWithType: kGLRealDataFormat forEquation:scalar.equation];
    } else if (variable.rank == 1) {
        GLFunction *function = (GLFunction *) variable;
        resultVar=[[function class] functionOfType: kGLRealDataFormat withDimensions: function.dimensions forEquation: function.equation];
    }  else if (variable.rank == 2) {
        [NSException raise: @"BadMath" format: @"I don't know how to take the cosine of a linear transformation."];
    }
    
	if (( self = [super initWithResult: @[resultVar] operand: @[variable] ]))
	{
		GLFunction *resultVariable = self.result[0];
		GLFunction *operandVariable = self.operand[0];
		
		const int numElements = (int) resultVariable.nDataElements;
		resultVariable.isPurelyReal = operandVariable.isPurelyReal;
		resultVariable.isPurelyImaginary = operandVariable.isPurelyImaginary;
		self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			NSMutableData *result = resultArray[0];
			NSMutableData *operand = operandArray[0];
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
/*		GLHyperbolicSineOperation				*/
/************************************************/

@implementation GLHyperbolicSineOperation

- (GLHyperbolicSineOperation *) initWithVariable: (GLVariable *) variable
{
    if (!variable.isPurelyReal) {
        [NSException raise: @"NotYetImplemented" format: @"You can only take the hyperbolic sine of real numbers."];
    }
    
	GLVariable *resultVar;
    
    if (variable.rank == 0) {
        GLScalar *scalar = (GLScalar *) variable;
        resultVar=[[GLScalar alloc] initWithType: kGLRealDataFormat forEquation:scalar.equation];
    } else if (variable.rank == 1) {
        GLFunction *function = (GLFunction *) variable;
        resultVar=[[function class] functionOfType: kGLRealDataFormat withDimensions: function.dimensions forEquation: function.equation];
    }  else if (variable.rank == 2) {
        [NSException raise: @"BadMath" format: @"I don't know how to take the tanh of a linear transformation."];
    }
    
	if (( self = [super initWithResult: @[resultVar] operand: @[variable] ]))
	{
		GLFunction *resultVariable = self.result[0];
		GLFunction *operandVariable = self.operand[0];
		
		const int numElements = (int) resultVariable.nDataElements;
		resultVariable.isPurelyReal = operandVariable.isPurelyReal;
		resultVariable.isPurelyImaginary = operandVariable.isPurelyImaginary;
		self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			NSMutableData *result = resultArray[0];
			NSMutableData *operand = operandArray[0];
			vGL_vvsinh( result.mutableBytes, operand.bytes, &numElements );
		};
        self.graphvisDescription = @"sinh";
	}
	
    return self;
}

- (BOOL) canOperateInPlace {
	return YES;
}

@end

/************************************************/
/*		GLInverseHyperbolicSineOperation		*/
/************************************************/

@implementation GLInverseHyperbolicSineOperation

- (GLInverseHyperbolicSineOperation *) initWithVariable: (GLVariable *) variable
{
    if (!variable.isPurelyReal) {
        [NSException raise: @"NotYetImplemented" format: @"You can only take the hyperbolic sine of real numbers."];
    }
    
	GLVariable *resultVar;
    
    if (variable.rank == 0) {
        GLScalar *scalar = (GLScalar *) variable;
        resultVar=[[GLScalar alloc] initWithType: kGLRealDataFormat forEquation:scalar.equation];
    } else if (variable.rank == 1) {
        GLFunction *function = (GLFunction *) variable;
        resultVar=[[function class] functionOfType: kGLRealDataFormat withDimensions: function.dimensions forEquation: function.equation];
    }  else if (variable.rank == 2) {
        [NSException raise: @"BadMath" format: @"I don't know how to take the tanh of a linear transformation."];
    }
    
	if (( self = [super initWithResult: @[resultVar] operand: @[variable] ]))
	{
		GLFunction *resultVariable = self.result[0];
		GLFunction *operandVariable = self.operand[0];
		
		const int numElements = (int) resultVariable.nDataElements;
		resultVariable.isPurelyReal = operandVariable.isPurelyReal;
		resultVariable.isPurelyImaginary = operandVariable.isPurelyImaginary;
		self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			NSMutableData *result = resultArray[0];
			NSMutableData *operand = operandArray[0];
			vGL_vvasinh( result.mutableBytes, operand.bytes, &numElements );
		};
        self.graphvisDescription = @"asinh";
	}
	
    return self;
}

- (BOOL) canOperateInPlace {
	return YES;
}

@end

/************************************************/
/*		GLHyperbolicCosineOperation				*/
/************************************************/

@implementation GLHyperbolicCosineOperation

- (GLHyperbolicCosineOperation *) initWithVariable: (GLVariable *) variable
{
    if (!variable.isPurelyReal) {
        [NSException raise: @"NotYetImplemented" format: @"You can only take the hyperbolic cosine of real numbers."];
    }
    
	GLVariable *resultVar;
    
    if (variable.rank == 0) {
        GLScalar *scalar = (GLScalar *) variable;
        resultVar=[[GLScalar alloc] initWithType: kGLRealDataFormat forEquation:scalar.equation];
    } else if (variable.rank == 1) {
        GLFunction *function = (GLFunction *) variable;
        resultVar=[[function class] functionOfType: kGLRealDataFormat withDimensions: function.dimensions forEquation: function.equation];
    }  else if (variable.rank == 2) {
        [NSException raise: @"BadMath" format: @"I don't know how to take the cosh of a linear transformation."];
    }
    
	if (( self = [super initWithResult: @[resultVar] operand: @[variable] ]))
	{
		GLFunction *resultVariable = self.result[0];
		GLFunction *operandVariable = self.operand[0];
		
		const int numElements = (int) resultVariable.nDataElements;
		resultVariable.isPurelyReal = operandVariable.isPurelyReal;
		resultVariable.isPurelyImaginary = operandVariable.isPurelyImaginary;
		self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			NSMutableData *result = resultArray[0];
			NSMutableData *operand = operandArray[0];
			vGL_vvcosh( result.mutableBytes, operand.bytes, &numElements );
		};
        self.graphvisDescription = @"cosh";
	}
	
    return self;
}

- (BOOL) canOperateInPlace {
	return YES;
}

@end

/************************************************/
/*		GLInverseHyperbolicCosineOperation		*/
/************************************************/

@implementation GLInverseHyperbolicCosineOperation

- (GLInverseHyperbolicCosineOperation *) initWithVariable: (GLVariable *) variable
{
    if (!variable.isPurelyReal) {
        [NSException raise: @"NotYetImplemented" format: @"You can only take the inverse hyperbolic cosine of real numbers."];
    }
    
	GLVariable *resultVar;
    
    if (variable.rank == 0) {
        GLScalar *scalar = (GLScalar *) variable;
        resultVar=[[GLScalar alloc] initWithType: kGLRealDataFormat forEquation:scalar.equation];
    } else if (variable.rank == 1) {
        GLFunction *function = (GLFunction *) variable;
        resultVar=[[function class] functionOfType: kGLRealDataFormat withDimensions: function.dimensions forEquation: function.equation];
    }  else if (variable.rank == 2) {
        [NSException raise: @"BadMath" format: @"I don't know how to take the acosh of a linear transformation."];
    }
    
	if (( self = [super initWithResult: @[resultVar] operand: @[variable] ]))
	{
		GLFunction *resultVariable = self.result[0];
		GLFunction *operandVariable = self.operand[0];
		
		const int numElements = (int) resultVariable.nDataElements;
		resultVariable.isPurelyReal = operandVariable.isPurelyReal;
		resultVariable.isPurelyImaginary = operandVariable.isPurelyImaginary;
		self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			NSMutableData *result = resultArray[0];
			NSMutableData *operand = operandArray[0];
			vGL_vvacosh( result.mutableBytes, operand.bytes, &numElements );
		};
        self.graphvisDescription = @"acosh";
	}
	
    return self;
}

- (BOOL) canOperateInPlace {
	return YES;
}

@end

/************************************************/
/*		GLHyperbolicTangentOperation            */
/************************************************/

@implementation GLHyperbolicTangentOperation

- (GLHyperbolicTangentOperation *) initWithVariable: (GLVariable *) variable
{
    if (!variable.isPurelyReal) {
        [NSException raise: @"NotYetImplemented" format: @"You can only take the hyperbolic tangent of real numbers."];
    }
    
	GLVariable *resultVar;
    
    if (variable.rank == 0) {
        GLScalar *scalar = (GLScalar *) variable;
        resultVar=[[GLScalar alloc] initWithType: kGLRealDataFormat forEquation:scalar.equation];
    } else if (variable.rank == 1) {
        GLFunction *function = (GLFunction *) variable;
        resultVar=[[function class] functionOfType: kGLRealDataFormat withDimensions: function.dimensions forEquation: function.equation];
    }  else if (variable.rank == 2) {
        [NSException raise: @"BadMath" format: @"I don't know how to take the tanh of a linear transformation."];
    }
    
	if (( self = [super initWithResult: @[resultVar] operand: @[variable] ]))
	{
		GLFunction *resultVariable = self.result[0];
		GLFunction *operandVariable = self.operand[0];
		
		const int numElements = (int) resultVariable.nDataElements;
		resultVariable.isPurelyReal = operandVariable.isPurelyReal;
		resultVariable.isPurelyImaginary = operandVariable.isPurelyImaginary;
		self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			NSMutableData *result = resultArray[0];
			NSMutableData *operand = operandArray[0];
			vGL_vvtanh( result.mutableBytes, operand.bytes, &numElements );
		};
        self.graphvisDescription = @"tanh";
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

- (GLInverseTangentOperation *) initWithVariable: (GLVariable *) variable
{
    if (!variable.isPurelyReal) {
        [NSException raise: @"NotYetImplemented" format: @"You can only take the inverse tangent of real numbers."];
    }
    
	GLVariable *resultVar;
    
    if (variable.rank == 0) {
        GLScalar *scalar = (GLScalar *) variable;
        resultVar=[[GLScalar alloc] initWithType: kGLRealDataFormat forEquation:scalar.equation];
    } else if (variable.rank == 1) {
        GLFunction *function = (GLFunction *) variable;
        resultVar=[[function class] functionOfType: kGLRealDataFormat withDimensions: function.dimensions forEquation: function.equation];
    }  else if (variable.rank == 2) {
        [NSException raise: @"BadMath" format: @"I don't know how to take the inverse tangent of a linear transformation."];
    }
    
	if (( self = [super initWithResult: @[resultVar] operand: @[variable] ]))
	{
		GLFunction *resultVariable = self.result[0];
		GLFunction *operandVariable = self.operand[0];
		
		const int numElements = (int) resultVariable.nDataElements;
		resultVariable.isPurelyReal = operandVariable.isPurelyReal;
		resultVariable.isPurelyImaginary = operandVariable.isPurelyImaginary;
		self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			NSMutableData *result = resultArray[0];
			NSMutableData *operand = operandArray[0];
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
- (GLSquareRootOperation *) initWithVariable: (GLVariable *) variable
{
    if (!variable.isPurelyReal) {
        [NSException raise: @"NotYetImplemented" format: @"You can only take the square root of real numbers, for now."];
    }
    
	GLVariable *resultVar;
    
    if (variable.rank == 0) {
        GLScalar *scalar = (GLScalar *) variable;
        resultVar=[[GLScalar alloc] initWithType: kGLRealDataFormat forEquation:scalar.equation];
    } else if (variable.rank == 1) {
        GLFunction *function = (GLFunction *) variable;
        resultVar=[[function class] functionOfType: kGLRealDataFormat withDimensions: function.dimensions forEquation: function.equation];
    }  else if (variable.rank == 2) {
        [NSException raise: @"BadMath" format: @"I don't know how to take the square root of a linear transformation."];
    }
    
	if (( self = [super initWithResult: @[resultVar] operand: @[variable] ]))
	{
		GLFunction *resultVariable = self.result[0];
		GLFunction *operandVariable = self.operand[0];
		
		const int numElements = (int) resultVariable.nDataElements;
		resultVariable.isPurelyReal = operandVariable.isPurelyReal;
		resultVariable.isPurelyImaginary = operandVariable.isPurelyImaginary;
		self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			NSMutableData *result = resultArray[0];
			NSMutableData *operand = operandArray[0];
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

- (GLFourierTransformOperation *) initWithFunction: (GLFunction *) variable
{
	NSMutableArray *transformedDimensions = [NSMutableArray array];
	for (GLDimension *aDim in variable.dimensions) {
		[transformedDimensions addObject: [aDim fourierTransformedDimension]];
	}
	
	// TODO - for the moment I'm assuming that anything NOT in the frequency domain is real.
	GLDataFormat format = !variable.isFrequencyDomain ? kGLSplitComplexDataFormat : kGLRealDataFormat;
	GLFunction *resultVariable = [[variable class] functionOfType: format withDimensions: transformedDimensions forEquation: variable.equation];
	
	if (( self = [super initWithResult: @[resultVariable] operand: @[variable]] ))
	{
		NSNumber *key = [[GLFourierTransformPool sharedFourierTransformPool] keyForDimensions: resultVariable.dimensions];
		if (resultVariable.isFrequencyDomain) {
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFourierTransform *transform = [[GLFourierTransformPool sharedFourierTransformPool] transformForKey: key];
				GLSplitComplex fbar = splitComplexFromData(resultArray[0]);
				GLFloat *f = (void *) ((NSMutableData *)operandArray[0]).bytes;
				
				[transform transform: f forward: &fbar];
				
				[[GLFourierTransformPool sharedFourierTransformPool] returnTransform: transform];
			};
		} else  {
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFourierTransform *transform = [[GLFourierTransformPool sharedFourierTransformPool] transformForKey: key];
				GLSplitComplex fbar = splitComplexFromData(operandArray[0]);
				GLFloat *f = (void *) ((NSMutableData *)resultArray[0]).bytes;
				
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

- (GLSwapComplexOperation *) initWithFunction: (GLFunction *) variable
{
	GLFunction *resultVariable = [[variable class] functionOfComplexTypeWithDimensions: variable.dimensions forEquation: variable.equation];
	resultVariable.isPurelyImaginary = variable.isPurelyReal;
	resultVariable.isPurelyReal = variable.isPurelyImaginary;
	
	if (( self = [super initWithResult: @[resultVariable] operand: @[variable]] ))
	{
		GLFunction *operandVariable = self.operand[0];
		
		NSUInteger numBytes = resultVariable.nDataPoints*sizeof(GLFloat);
		if (operandVariable.isComplex) {
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLSplitComplex fromSplit = splitComplexFromData(operandArray[0]);
				GLSplitComplex toSplit = splitComplexFromData(resultArray[0]);
				
				memcpy( toSplit.realp, fromSplit.imagp, numBytes );
				memcpy( toSplit.imagp, fromSplit.realp, numBytes );	
			};
            self.graphvisDescription = @"swap complex (complex)";
		} else {
			NSUInteger numPoints = resultVariable.nDataPoints;
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLSplitComplex toSplit = splitComplexFromData(resultArray[0]);
				NSMutableData *operand = operandArray[0];
				
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

- (GLCopyVariableOperation *) initWithFunction: (GLFunction *) variable
{
	if (( self = [super initWithOperand: @[variable]] ))
	{
		GLFunction *resultVariable = self.result[0];
		GLFunction *operandVariable = self.operand[0];
		
        resultVariable.name = variable.name;
		NSUInteger numBytes = resultVariable.nDataElements*sizeof(GLFloat);
		resultVariable.isPurelyReal = operandVariable.isPurelyReal;
		resultVariable.isPurelyImaginary = operandVariable.isPurelyImaginary;
		self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			NSMutableData *result = resultArray[0];
			NSMutableData *operand = operandArray[0];
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

- (GLMaxOperation *) initWithFunction: (GLFunction *) variable
{
	GLScalar *resultVariable = [GLScalar scalarWithType: kGLRealDataFormat forEquation: variable.equation];
	
	if (( self = [super initWithResult: @[resultVariable] operand: @[variable]] ))
	{		
		NSUInteger nDataElements = ((GLFunction *)self.operand[0]).nDataElements;
		self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			NSMutableData *result = resultArray[0];
			NSMutableData *operand = operandArray[0];
			vGL_maxv( (GLFloat *) operand.bytes, 1, result.mutableBytes, nDataElements);
		};
        self.graphvisDescription = @"max";
	}
	
    return self;
}

@end

/************************************************/
/*		GLMinOperation							*/
/************************************************/

@implementation GLMinOperation

- (GLMinOperation *) initWithFunction: (GLFunction *) variable
{
	GLScalar *resultVariable = [GLScalar scalarWithType: kGLRealDataFormat forEquation: variable.equation];
	
	if (( self = [super initWithResult: @[resultVariable] operand: @[variable]] ))
	{
		NSUInteger nDataElements = ((GLFunction *)self.operand[0]).nDataElements;
		self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			NSMutableData *result = resultArray[0];
			NSMutableData *operand = operandArray[0];
			vGL_minv( (GLFloat *) operand.bytes, 1, result.mutableBytes, nDataElements);
		};
        self.graphvisDescription = @"min";
	}
	
    return self;
}

@end

/************************************************/
/*		GLAverageOperation						*/
/************************************************/

@implementation GLAverageOperation

- (GLAverageOperation *) initWithFunction: (GLFunction *) variable dimensionIndex: (NSUInteger) index
{
    GLDimension *dim = variable.dimensions[index];
    return [self initWithFunction: variable dimensionIndex: index range: NSMakeRange(0, dim.nPoints)];
}

- (GLAverageOperation *) initWithFunction: (GLFunction *) variable dimensionIndex: (NSUInteger) index range: (NSRange) aRange;
{
	NSMutableArray *newDimensions = [variable.dimensions mutableCopy];
	[newDimensions removeObjectAtIndex: index];
	
    
	GLVariable *resultVariable;
    
    if (newDimensions.count) {
        resultVariable = [GLFunction functionOfType: variable.dataFormat withDimensions: newDimensions forEquation: variable.equation];
    } else {
        resultVariable = [GLScalar scalarWithType: variable.dataFormat forEquation: variable.equation];
    }
	
    // index = i*ny*nz + j*nz + k
    // For example: for a given (i,k), we want to sum over all j's. The stride between j elements, is nz and there are clearly ny of them.
    
    // The spacing between the elements that we want to sum over.
    // If we're collapsing the last index
    NSUInteger summingStride=variable.matrixDescription.strides[index].stride;
    
    // The number of elements that we want to sum over
    NSUInteger summingPoints = variable.matrixDescription.strides[index].nPoints;
    
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
    
    GLDimension *dim = variable.dimensions[index];
    
    variableOperation op;
    NSString *graphvisDescription;
    NSArray *buffers = @[];
    if (dim.isEvenlySampled)
    {
        if (resultVariable.dataFormat == kGLRealDataFormat)
        {
            op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
                NSMutableData *result = resultArray[0];
                NSMutableData *operand = operandArray[0];
                dispatch_apply( nDataPoints, dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0), ^(size_t iteration) {
                    NSUInteger index = (iteration/divisor)*a + (iteration%divisor)*b + aRange.location * summingStride;
                    GLFloat *op = (GLFloat *) operand.bytes;
                    GLFloat *res = (GLFloat *) result.mutableBytes;
                    vGL_meanv(&op[index], summingStride, &res[iteration], aRange.length);
                });
                
            };
            graphvisDescription = [NSString stringWithFormat: @"average dimension %lu (real formatting)", index] ;
        }
        else if (resultVariable.dataFormat == kGLSplitComplexDataFormat)
        {
            op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
                
                dispatch_apply( nDataPoints, dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0), ^(size_t iteration) {
                    NSUInteger index = (iteration/divisor)*a + (iteration%divisor)*b+ aRange.location * summingStride;
                    GLSplitComplex toSplit = splitComplexFromData(resultArray[0]);
                    GLSplitComplex fromSplit = splitComplexFromData(operandArray[0]);
                    vGL_meanv(fromSplit.realp + index*sizeof(GLFloat), summingStride, toSplit.realp + index*sizeof(GLFloat), aRange.length);
                    vGL_meanv(fromSplit.imagp + index*sizeof(GLFloat), summingStride, toSplit.imagp + index*sizeof(GLFloat), aRange.length);
                    
                });
                
            };
            graphvisDescription = [NSString stringWithFormat: @"average dimension %lu (split complex formatting)", index] ;
        }
    } else {
        NSData *dimensionData = dim.data;
        NSData *dimensionDiffData = [[GLMemoryPool sharedMemoryPool] dataWithLength: dimensionData.length];
        GLFloat *dimPoints = (GLFloat *) dimensionData.bytes;
        GLFloat *dimDiffPoints = (GLFloat *) dimensionDiffData.bytes;
        vGL_vsub(&(dimPoints[0]), 1, &(dimPoints[1]), 1, dimDiffPoints, 1, dim.nPoints-1); // vsub does C = B - A
        GLFloat length = [dim valueAtIndex: aRange.location+aRange.length-1] - [dim valueAtIndex: aRange.location];
        GLBuffer *aBuffer = [[GLBuffer alloc] initWithLength: variable.dataBytes];
        buffers = @[aBuffer];
#warning check this!!!
        if (resultVariable.dataFormat == kGLRealDataFormat)
        {
            op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
                NSMutableData *result = resultArray[0];
                NSMutableData *operand = operandArray[0];
                NSMutableData *buffer = bufferArray[0];
                dispatch_apply( nDataPoints, dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0), ^(size_t iteration) {
                    NSUInteger index = (iteration/divisor)*a + (iteration%divisor)*b + aRange.location * summingStride;
                    GLFloat *op = (GLFloat *) operand.bytes;
                    GLFloat *res = (GLFloat *) result.mutableBytes;
                    GLFloat *buf = (GLFloat *) buffer.mutableBytes;
            
                    // now integrate: 1/2 sum_{n=0}^{n-2} (x_{n+1} - x_n) * ( f(n+1) + f(n) )
                    vGL_vadd( &op[index], summingStride, &op[index+summingStride], summingStride, &buf[index], summingStride, aRange.length-1 );
                    vGL_vmul(&buf[index], summingStride, &(dimDiffPoints[aRange.location]), 1, &buf[index], summingStride, aRange.length-1);
                    GLFloat sum = 0.0;
                    vGL_sve(&(buf[index]), summingStride, &sum, aRange.length-1);
                    res[iteration] = sum/(2*length);
                });
                
            };
            graphvisDescription = [NSString stringWithFormat: @"average dimension %lu (real formatting)", index] ;
        }
        else if (resultVariable.dataFormat == kGLSplitComplexDataFormat)
        {
            [NSException raise:@"NotYetImplemented" format:@"NotYetImplemented"];
        }
    }
    
	if (( self = [super initWithResult: @[resultVariable] operand: @[variable] buffers: buffers operation: op] ))
	{
        self.dimIndex = index;
        self.graphvisDescription = graphvisDescription;
	}
	
    return self;
}

- (GLAverageOperation *) initWithFunction: (GLFunction *) variable
{
	GLVariable *resultVariable = [GLScalar scalarWithType: variable.dataFormat forEquation: variable.equation];
	
	if (( self = [super initWithResult: @[resultVariable] operand: @[variable]] ))
	{
        self.dimIndex = NSNotFound;
		NSUInteger nPoints = variable.matrixDescription.nPoints;
		if (resultVariable.dataFormat == kGLRealDataFormat)
		{
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *a = [operandArray[0] mutableBytes];
				GLFloat *b = [resultArray[0] mutableBytes];
				vGL_meanv(a, 1, b, nPoints);
				
			};
            self.graphvisDescription = [NSString stringWithFormat: @"average all dimensions (real formatting)"] ;
		}
        else if (resultVariable.dataFormat == kGLSplitComplexDataFormat)
		{
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLSplitComplex toSplit = splitComplexFromData(resultArray[0]);
				GLSplitComplex fromSplit = splitComplexFromData(operandArray[0]);
				vGL_meanv(fromSplit.realp,1, toSplit.realp, nPoints);
				vGL_meanv(fromSplit.imagp,1, toSplit.imagp, nPoints);
				
			};
            self.graphvisDescription = [NSString stringWithFormat: @"average all dimensions (split complex formatting)"] ;
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

/************************************************/
/*		GLSummationOperation					*/
/************************************************/

@implementation GLSummationOperation

- (GLSummationOperation *) initWithFunction: (GLFunction *) variable dimensionIndex: (NSUInteger) index
{
	NSMutableArray *newDimensions = [variable.dimensions mutableCopy];
	[newDimensions removeObjectAtIndex: index];
	
    
	GLVariable *resultVariable;
    
    if (newDimensions.count) {
        resultVariable = [GLFunction functionOfType: variable.dataFormat withDimensions: newDimensions forEquation: variable.equation];
    } else {
        resultVariable = [GLScalar scalarWithType: variable.dataFormat forEquation: variable.equation];
    }
	
	if (( self = [super initWithResult: @[resultVariable] operand: @[variable]] ))
	{
        self.dimIndex = index;
        
		// index = i*ny*nz + j*nz + k
		// For example: for a given (i,k), we want to sum over all j's. The stride between j elements, is nz and there are clearly ny of them.
		
		// The spacing between the elements that we want to sum over.
		// If we're collapsing the last index
		NSUInteger summingStride=variable.matrixDescription.strides[index].stride;
        
		// The number of elements that we want to sum over
		NSUInteger summingPoints = variable.matrixDescription.strides[index].nPoints;
		
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
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				NSMutableData *result = resultArray[0];
				NSMutableData *operand = operandArray[0];
				dispatch_apply( nDataPoints, dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0), ^(size_t iteration) {
                    NSUInteger index = (iteration/divisor)*a + (iteration%divisor)*b;
					GLFloat *op = (GLFloat *) operand.bytes;
					GLFloat *res = (GLFloat *) result.mutableBytes;
					vGL_sve(&op[index], summingStride, &res[iteration], summingPoints);
				});
				
			};
            self.graphvisDescription = [NSString stringWithFormat: @"average dimension %lu (real formatting)", index] ;
		}
        else if (resultVariable.dataFormat == kGLSplitComplexDataFormat)
		{
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				
				dispatch_apply( nDataPoints, dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0), ^(size_t iteration) {
                    NSUInteger index = (iteration/divisor)*a + (iteration%divisor)*b;
                    GLSplitComplex toSplit = splitComplexFromData(resultArray[0]);
                    GLSplitComplex fromSplit = splitComplexFromData(operandArray[0]);
					vGL_sve(fromSplit.realp + index*sizeof(GLFloat), summingStride, toSplit.realp + index*sizeof(GLFloat), summingPoints);
                    vGL_sve(fromSplit.imagp + index*sizeof(GLFloat), summingStride, toSplit.imagp + index*sizeof(GLFloat), summingPoints);
                    
				});
				
			};
            self.graphvisDescription = [NSString stringWithFormat: @"average dimension %lu (split complex formatting)", index] ;
		}
        
	}
	
    return self;
}

- (GLSummationOperation *) initWithFunction: (GLFunction *) variable
{
	GLVariable *resultVariable = [GLScalar scalarWithType: variable.dataFormat forEquation: variable.equation];
	
	if (( self = [super initWithResult: @[resultVariable] operand: @[variable]] ))
	{
        self.dimIndex = NSNotFound;
		NSUInteger nPoints = variable.matrixDescription.nPoints;
		if (resultVariable.dataFormat == kGLRealDataFormat)
		{
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLFloat *a = [operandArray[0] mutableBytes];
				GLFloat *b = [resultArray[0] mutableBytes];
				vGL_sve(a, 1, b, nPoints);
				
			};
            self.graphvisDescription = [NSString stringWithFormat: @"average all dimensions (real formatting)"] ;
		}
        else if (resultVariable.dataFormat == kGLSplitComplexDataFormat)
		{
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLSplitComplex toSplit = splitComplexFromData(resultArray[0]);
				GLSplitComplex fromSplit = splitComplexFromData(operandArray[0]);
				vGL_sve(fromSplit.realp,1, toSplit.realp, nPoints);
				vGL_sve(fromSplit.imagp,1, toSplit.imagp, nPoints);
				
			};
            self.graphvisDescription = [NSString stringWithFormat: @"average all dimensions (split complex formatting)"] ;
		}
		
	}
	
    return self;
}

- (BOOL) isEqualToOperation: (id) otherOperation {
    if ( ![super isEqualToOperation: otherOperation] )  {
        return NO;
    }
    
    GLSummationOperation * op = otherOperation;
    if (self.dimIndex != op.dimIndex) {
        return NO;
    }
    
    return YES;
}

@end

/************************************************/
/*		GLIntegrationToLimitsOperation			*/
/************************************************/

@implementation GLIntegrationToLimitsOperation

- (GLIntegrationToLimitsOperation *) initWithFunction: (GLFunction *) variable
{
	if (variable.dataFormat != kGLRealDataFormat || variable.dimensions.count != 1 || ![variable.dimensions[0] isEvenlySampled]) {
        [NSException raise: @"BadFormat" format: @"This operation can only take evenly sampled one-dimensional real variables."];
    }

	GLScalar *result = [[GLScalar alloc] initWithType: kGLRealDataFormat forEquation: variable.equation];
	NSUInteger N = variable.nDataElements;
	GLFloat deltaX = [variable.dimensions[0] sampleInterval];
	variableOperation op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
		GLFloat *A = (GLFloat *) [operandArray[0] bytes];
		GLFloat *B = (GLFloat *) [resultArray[0] bytes];
	
		// sum it. sum = 1/(2*h) * (f_0 + 2*f_1 + 2*f_2 + ... + 2*f_{N-1} + 2*f_N
		vGL_sve(A, 1, B, N);
		*B = deltaX*(*B - A[0]/2.0 - A[N-1]/2.0);
	};
	
	
	if (( self = [super initWithResult: @[result] operand: @[variable] buffers: @[] operation: op] )) {
        
    }
    return self;
}

@end

/************************************************/
/*		GLIntegrationOperation		*/
/************************************************/

@implementation GLIntegrationOperation

- (GLIntegrationOperation *) initWithFunction: (GLFunction *) variable
{
    if (variable.dataFormat != kGLRealDataFormat || variable.dimensions.count != 1 ) {
        [NSException raise: @"BadFormat" format: @"This operation can only take one-dimensional real variables."];
    }
    
    GLFunction *result = (GLFunction *) [GLVariable variableWithPrototype: variable];
    GLDimension *dim = variable.dimensions[0];
    
    variableOperation op;
    if ([variable.dimensions[0] isEvenlySampled]) {
        NSUInteger N = variable.nDataElements;
        GLFloat deltaX = [variable.dimensions[0] sampleInterval];
        op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
            GLFloat *A = (GLFloat *) [operandArray[0] bytes];
            GLFloat *B = (GLFloat *) [resultArray[0] bytes];
            vGL_vtrapz( A, 1, &deltaX, B, 1, N);
        };
    } else {
        NSUInteger nPointsMinusOne = dim.nPoints-1;
        NSUInteger nPoints = dim.nPoints;
        NSData *dimensionData = dim.data;
        NSData *dimensionDiffData = [[GLMemoryPool sharedMemoryPool] dataWithLength: dimensionData.length];
        GLFloat *dimPoints = (GLFloat *) dimensionData.bytes;
        GLFloat *dimDiffPoints = (GLFloat *) dimensionDiffData.bytes;
        vGL_vsub(&(dimPoints[0]), 1, &(dimPoints[1]), 1, dimDiffPoints, 1, nPointsMinusOne); // vsub does C = B - A
        
        op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
            GLFloat *A = (GLFloat *) [operandArray[0] bytes];
            GLFloat *B = (GLFloat *) [resultArray[0] bytes];
            GLFloat *dimDiff = (GLFloat *) dimensionDiffData.bytes;
            GLFloat half = 0.5;
            
            vGL_vadd(A, 1, &(A[1]), 1, &(B[1]), 1, nPointsMinusOne);
            vGL_vmul(&(B[1]), 1, dimDiff, 1, &(B[1]), 1, nPointsMinusOne);
            vGL_vrsum(&(B[0]), 1, &half, &(B[0]), 1, nPoints);
        };
    }
    
    if (( self = [super initWithResult: @[result] operand: @[variable] buffers: @[] operation: op] )) {
        
    }
    return self;
}

@end

/************************************************/
/*		GLInterleavedToSplitComplexOperation    */
/************************************************/

@implementation GLInterleavedToSplitComplexOperation

- (GLInterleavedToSplitComplexOperation *) initWithVariable: (GLVariable *) variable
{
    if (variable.dataFormat != kGLInterleavedComplexDataFormat) {
        [NSException raise: @"BadFormat" format: @"This operation can only take variables in interleaved complex format"];
    }
    
	GLVariable *result;
    if (variable.rank == 0) {
        GLScalar *scalar = (GLScalar *) variable;
        result = [[GLScalar alloc] initWithType: kGLSplitComplexDataFormat forEquation:scalar.equation];
    } else if (variable.rank == 1) {
        GLFunction *function = (GLFunction *) variable;
        result = [[function class] functionOfType: kGLSplitComplexDataFormat withDimensions: function.dimensions forEquation: function.equation];
    }  else if (variable.rank == 2) {
        GLLinearTransform *matrix = (GLLinearTransform *) variable;
        result = [GLLinearTransform transformOfType: kGLSplitComplexDataFormat withFromDimensions: matrix.fromDimensions toDimensions: matrix.toDimensions inFormat: matrix.matrixFormats forEquation:matrix.equation matrix:nil];
    }
    
	if (( self = [super initWithResult: @[result] operand: @[variable] ] ))
	{
        NSUInteger numPoints = variable.nDataPoints;
        self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
            GLFloatComplex *A = [operandArray[0] mutableBytes];
            GLSplitComplex B = splitComplexFromData(resultArray[0]);
            vGL_ctoz( (DSPComplex *)A, 2, &B, 1, numPoints );
        };
        self.graphvisDescription = @"interleaved->split";
	}
	
    return self;
}

- (BOOL) canOperateInPlace {
	return YES;
}

@end

/************************************************/
/*		GLSplitToInterleavedComplexOperation    */
/************************************************/

@implementation GLSplitToInterleavedComplexOperation

- (GLSplitToInterleavedComplexOperation *) initWithVariable: (GLVariable *) variable
{
    if (variable.dataFormat != kGLSplitComplexDataFormat) {
        [NSException raise: @"BadFormat" format: @"This operation can only take variables in interleaved complex format"];
    }
    
	GLVariable *result;
    if (variable.rank == 0) {
        GLScalar *scalar = (GLScalar *) variable;
        result = [[GLScalar alloc] initWithType: kGLInterleavedComplexDataFormat forEquation:scalar.equation];
    } else if (variable.rank == 1) {
        GLFunction *function = (GLFunction *) variable;
        result = [[function class] functionOfType: kGLInterleavedComplexDataFormat withDimensions: function.dimensions forEquation: function.equation];
    }  else if (variable.rank == 2) {
        GLLinearTransform *matrix = (GLLinearTransform *) variable;
        result = [GLLinearTransform transformOfType: kGLInterleavedComplexDataFormat withFromDimensions: matrix.fromDimensions toDimensions: matrix.toDimensions inFormat: matrix.matrixFormats forEquation:matrix.equation matrix:nil];
    }
    
	if (( self = [super initWithResult: @[result] operand: @[variable] ] ))
	{
        NSUInteger numPoints = variable.nDataPoints;
        self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
            GLFloatComplex *A = [resultArray[0] mutableBytes];
            GLSplitComplex B = splitComplexFromData(operandArray[0]);
            vGL_ztoc( &B, 1, (DSPComplex *)A, 2, numPoints );
        };
        self.graphvisDescription = @"split->interleaved";
	}
	
    return self;
}

- (BOOL) canOperateInPlace {
	return YES;
}

@end

/************************************************/
/*		GLDataTransposeOperation				*/
/************************************************/

@implementation GLDataTransposeOperation

- (GLDataTransposeOperation *) initWithLinearTransform: (GLLinearTransform *) transform
{
	NSUInteger denseIndex = NSNotFound;
	NSUInteger lastNonTrivialNonDenseIndex = NSNotFound;
	NSUInteger lastTrivialIndex = NSNotFound;
	NSUInteger totalTrivialPoints = 1;
	NSUInteger numDenseIndices = 0;
    for ( NSUInteger index=0; index < transform.matrixFormats.count; index++) {
        GLMatrixFormat format = [transform.matrixFormats[index] unsignedIntegerValue];
        if (format == kGLDenseMatrixFormat) {
            denseIndex = index;
			numDenseIndices++;
        } else if (format == kGLIdentityMatrixFormat) {
            lastTrivialIndex = index;
			totalTrivialPoints *= [transform.fromDimensions[index] nPoints];
        }
        if (format != kGLIdentityMatrixFormat && format != kGLDenseMatrixFormat ) {
			lastNonTrivialNonDenseIndex = index;
		}
    }

	if (numDenseIndices != 1) {
		[NSException raise: @"BadInputFormat" format: @"The GLDataTransposeOperation can only take exactly one dense dimension."];
	}
	
	GLLinearTransform *result = [GLLinearTransform transformOfType: transform.dataFormat withFromDimensions: transform.fromDimensions toDimensions: transform.toDimensions inFormat: transform.matrixFormats forEquation:transform.equation matrix:nil];
    result.matrixOrder = transform.matrixOrder == kGLRowMatrixOrder ? kGLColumnMatrixOrder : kGLRowMatrixOrder;
	result.matrixDescription = [[GLMatrixDescription alloc] initWithLinearTransform: result];
	
	// How many 'inner loop' steps we need to take depends on how many other nontrivial dimensions there are.
	NSUInteger loopStride = lastNonTrivialNonDenseIndex == NSNotFound ? 0 : transform.matrixDescription.strides[lastNonTrivialNonDenseIndex].stride;
	NSUInteger totalLoops = transform.nDataPoints / transform.matrixDescription.strides[denseIndex].nPoints;
	NSUInteger nColumns = transform.matrixDescription.strides[denseIndex].nColumns;
	NSUInteger nRows = transform.matrixDescription.strides[denseIndex].nRows;
	NSUInteger denseStride = transform.matrixDescription.strides[denseIndex].stride;
	dispatch_queue_t globalQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
	
	variableOperation op;
	NSString *graphvisDescription;
	if (result.dataFormat == kGLRealDataFormat) {
		op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			dispatch_apply(totalLoops, globalQueue, ^(size_t iteration) {
				NSInteger startPosition = iteration*loopStride;
				GLFloat *A = [operandArray[0] mutableBytes];
				GLFloat *M = &(A[startPosition]);
				
				GLFloat *B = [resultArray[0] mutableBytes];
				GLFloat *N = &(B[startPosition]);
				
				vGL_mtrans( M, denseStride, N, denseStride, nColumns, nRows);
			});
        };
		graphvisDescription = @"data transpose (real)";
	} else if (result.dataFormat == kGLInterleavedComplexDataFormat) {
		op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			dispatch_apply(totalLoops, globalQueue, ^(size_t iteration) {
				NSInteger startPosition = iteration*loopStride;
				GLFloat *A = [operandArray[0] mutableBytes];
				GLFloat *M = &(A[startPosition]);
				GLFloat *iM = &(A[startPosition+1]);
				
				GLFloat *B = [resultArray[0] mutableBytes];
				GLFloat *N = &(B[startPosition]);
				GLFloat *iN = &(B[startPosition+1]);
				
				vGL_mtrans( M, denseStride, N, denseStride, nColumns, nRows);
				vGL_mtrans( iM, denseStride, iN, denseStride, nColumns, nRows);
			});
        };
		graphvisDescription = @"data transpose (interleaved complex)";
	} else if (result.dataFormat == kGLSplitComplexDataFormat) {
		op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			dispatch_apply(totalLoops, globalQueue, ^(size_t iteration) {
				NSInteger startPosition = iteration*loopStride;
				GLSplitComplex A = splitComplexFromData(operandArray[0]);
				GLFloat *M = &(A.realp[startPosition]);
				GLFloat *iM = &(A.imagp[startPosition]);
				
				GLSplitComplex B = splitComplexFromData(resultArray[0]);
				GLFloat *N = &(B.realp[startPosition]);
				GLFloat *iN = &(B.imagp[startPosition+1]);
				
				vGL_mtrans( M, denseStride, N, denseStride, nColumns, nRows);
				vGL_mtrans( iM, denseStride, iN, denseStride, nColumns, nRows);
			});
        };
		graphvisDescription = @"data transpose (split complex)";
    } else {
        [NSException raise: @"InvalidFormat" format: @"NotImplemented"];
    }
	
	if (( self = [super initWithResult: @[result] operand: @[transform] ] ))
	{
		self.operation = op;
        self.graphvisDescription = graphvisDescription;
	}
	
    return self;
}

- (BOOL) canOperateInPlace {
	return YES;
}

@end

/************************************************/
/*		GLDenseMatrixOperation					*/
/************************************************/

@implementation GLDenseMatrixOperation

- (GLDenseMatrixOperation *) initWithLinearTransform: (GLLinearTransform *) transform
{
	if (transform.matrixFormats.count != 1) {
		[NSException raise: @"BadInputFormat" format: @"The GLDenseMatrixOperation can only take exactly one dimension at this time."];
	}
	
	GLLinearTransform *result = [GLLinearTransform transformOfType: transform.dataFormat withFromDimensions: transform.fromDimensions toDimensions: transform.toDimensions inFormat: @[@(kGLDenseMatrixFormat)] forEquation:transform.equation matrix:nil];
	
	NSUInteger nDiagonalPoints = transform.matrixDescription.strides[0].nDiagonalPoints;
	NSUInteger nDiagonals = transform.matrixDescription.strides[0].nDiagonals;
	NSUInteger diagonalStride = transform.matrixDescription.strides[0].diagonalStride;
	
	NSUInteger nRows = result.matrixDescription.strides[0].nRows;
	NSUInteger nColumns = result.matrixDescription.strides[0].nColumns;
	NSUInteger rowStride = result.matrixDescription.strides[0].rowStride;
	NSUInteger colStride = result.matrixDescription.strides[0].columnStride;
	
	variableOperation op;
	NSString *graphvisDescription;
	if (result.dataFormat == kGLRealDataFormat) {
		op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			GLFloat *A = [operandArray[0] mutableBytes];
			GLFloat *B = [resultArray[0] mutableBytes];
			
			vGL_vclr(B,1,nRows*nColumns);
			
			NSUInteger bandwidth = (nDiagonals-1)/2;
			for (NSUInteger iDiagonal=0; iDiagonal<nDiagonals; iDiagonal++) {
				NSInteger iStart = (NSInteger)bandwidth - (NSInteger) iDiagonal > 0 ? (NSInteger)bandwidth - (NSInteger) iDiagonal : 0;
				NSInteger iEnd = (NSInteger)iDiagonal - (NSInteger) bandwidth > 0 ? (NSInteger)iDiagonal - (NSInteger) bandwidth : 0;
				for (NSUInteger i=iStart; i<nDiagonalPoints-iEnd; i++) {
						NSUInteger row = i; // The i-th point in the diagonal always corresponds to the i-th row.
						NSUInteger col = i + (iDiagonal-bandwidth);
						B[row*rowStride+col*colStride] = A[iDiagonal*diagonalStride+i];
				}
			}
        };
		graphvisDescription = @"matrix densification (real)";
	}
	
	if (( self = [super initWithResult: @[result] operand: @[transform] ] ))
	{
		self.operation = op;
        self.graphvisDescription = graphvisDescription;
	}
	
    return self;
}

- (BOOL) canOperateInPlace {
	return YES;
}

@end


/************************************************/
/*		GLMakeHermitianOperation                */
/************************************************/

@implementation GLMakeHermitianOperation

- (GLMakeHermitianOperation *) initWithFunction: (GLFunction *) variable
{
	if (!variable.isComplex) {
		[NSException raise: @"StupidImplementationException" format: @"Can't deal with real variables"];
	}
    GLFunction *resultVariable = (GLFunction *) [GLVariable variableWithPrototype: variable];
    
    NSUInteger nDataElements = resultVariable.nDataElements;
    BOOL halfComplex2D = resultVariable.dimensions.count == 2 && [resultVariable.dimensions[0] basisFunction] == kGLExponentialBasis && [resultVariable.dimensions[1] basisFunction] == kGLExponentialBasis && [resultVariable.dimensions[1] isStrictlyPositive];
    BOOL halfComplex3D = resultVariable.dimensions.count == 3 && [resultVariable.dimensions[0] basisFunction] == kGLDeltaBasis && [resultVariable.dimensions[1] basisFunction] == kGLExponentialBasis && [resultVariable.dimensions[2] basisFunction] == kGLExponentialBasis && [resultVariable.dimensions[2] isStrictlyPositive];
    if (nDataElements % 2 == 1) [NSException raise: @"StupidImplementationException" format: @"Can't deal with odd numbers"];
    
    variableOperation op;
    NSString *graphvisDescription;
    if (halfComplex2D)
    {
        NSUInteger kDimNPoints = [resultVariable.dimensions[0] nPoints];
        NSUInteger lDimNPoints = [resultVariable.dimensions[1] nPoints];
        op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
            
            GLSplitComplex C = splitComplexFromData( resultArray[0] );
            
            // Hermitian conjugates
            for (NSUInteger i=1; i<kDimNPoints/2; i++) {
                C.realp[(kDimNPoints-i)*lDimNPoints] = C.realp[i*lDimNPoints];
                C.imagp[(kDimNPoints-i)*lDimNPoints] = -C.imagp[i*lDimNPoints];
                
                C.realp[(kDimNPoints-i)*lDimNPoints+(lDimNPoints-1)] = C.realp[i*lDimNPoints+(lDimNPoints-1)];
                C.imagp[(kDimNPoints-i)*lDimNPoints+(lDimNPoints-1)] = -C.imagp[i*lDimNPoints+(lDimNPoints-1)];
            }
            
            // For the four self-conjugate components, that means that there can be no imaginary component
            C.imagp[0] = 0;
            C.imagp[(lDimNPoints-1)] = 0;
            C.imagp[(kDimNPoints/2)*lDimNPoints+0] = 0;
            C.imagp[(kDimNPoints/2)*lDimNPoints+(lDimNPoints-1)] = 0;
            
            // But that their real components should be doubled, to make all else equal.
            C.realp[0] = 2*C.realp[0];
            C.realp[(lDimNPoints-1)] = 2*C.realp[(lDimNPoints-1)];
            C.realp[(kDimNPoints/2)*lDimNPoints+0] = 2*C.realp[(kDimNPoints/2)*lDimNPoints+0];
            C.realp[(kDimNPoints/2)*lDimNPoints+(lDimNPoints-1)] = 2*C.realp[(kDimNPoints/2)*lDimNPoints+(lDimNPoints-1)];
        };
        graphvisDescription = [NSString stringWithFormat: @"make hermitian (2D)"];
    }
    else if (halfComplex3D)
    {
        NSUInteger zDimNPoints = [resultVariable.dimensions[0] nPoints];
        NSUInteger kDimNPoints = [resultVariable.dimensions[1] nPoints];
        NSUInteger lDimNPoints = [resultVariable.dimensions[2] nPoints];
        op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
            
            memcpy( [resultArray[0] mutableBytes], [operandArray[0] bytes], 2*zDimNPoints*kDimNPoints*lDimNPoints*sizeof(GLFloat));
            GLSplitComplex C = splitComplexFromData( resultArray[0] );
            
            // index=i*ny*nz+j*nz+k
            // index=(i*ny+j)*nz+k
            // my notation, (z*kDimNPoints+k)*lDimNPoints+l
            for (NSUInteger z=0; z<zDimNPoints; z++) {
                // Hermitian conjugates
                for (NSUInteger i=1; i<kDimNPoints/2; i++) {
                    C.realp[(z*kDimNPoints+(kDimNPoints-i))*lDimNPoints] = C.realp[(z*kDimNPoints+i)*lDimNPoints];
                    C.imagp[(z*kDimNPoints+(kDimNPoints-i))*lDimNPoints] = -C.imagp[(z*kDimNPoints+i)*lDimNPoints];
                    
                    C.realp[(z*kDimNPoints+(kDimNPoints-i))*lDimNPoints+(lDimNPoints-1)] = C.realp[(z*kDimNPoints+i)*lDimNPoints+(lDimNPoints-1)];
                    C.imagp[(z*kDimNPoints+(kDimNPoints-i))*lDimNPoints+(lDimNPoints-1)] = -C.imagp[(z*kDimNPoints+i)*lDimNPoints+(lDimNPoints-1)];
                }
                
                // For the four self-conjugate components, that means that there can be no imaginary component
                C.imagp[z*kDimNPoints*lDimNPoints+0] = 0;
                C.imagp[z*kDimNPoints*lDimNPoints+(lDimNPoints-1)] = 0;
                C.imagp[(z*kDimNPoints+(kDimNPoints/2))*lDimNPoints+0] = 0;
                C.imagp[(z*kDimNPoints+(kDimNPoints/2))*lDimNPoints+(lDimNPoints-1)] = 0;
                
                // But that their real components should be doubled, to make all else equal.
//                C.realp[z*kDimNPoints*lDimNPoints+0] *= 2;
//                C.realp[z*kDimNPoints*lDimNPoints+(lDimNPoints-1)] *= 2;
//                C.realp[(z*kDimNPoints+(kDimNPoints/2))*lDimNPoints+0] *= 2;
//                C.realp[(z*kDimNPoints+(kDimNPoints/2))*lDimNPoints+(lDimNPoints-1)] *= 2;
            }
        };
        graphvisDescription = [NSString stringWithFormat: @"make hermitian (3D)"];
    } else {
        [NSException raise: @"CannotMakeHermitian" format:@"Invalid case"];
    }
    
	if (( self = [super initWithResult: @[resultVariable] operand: @[variable] buffers: @[] operation: op ]))
	{
        self.graphvisDescription = graphvisDescription;
	}
	
    return self;
}

- (BOOL) canOperateInPlace {
	return YES;
}

@end
