//
//  GLVectorScalarOperations.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 1/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import "GLVectorScalarOperations.h"
#import "GLLinearTransform.h"

/************************************************/
/*		GLVectorScalarOperation					*/
/************************************************/

@implementation GLVectorScalarOperation

- (id) initWithVectorOperand: (GLTensor *) vOperand scalarOperand: (GLFloat) sOperand
{
	if (( self = [super initWithOperand: @[vOperand]] ))
	{
		self.scalarOperand = sOperand;
	}
	return self;
}

@synthesize scalarOperand;

- (BOOL) isEqualToOperation: (id) otherOperation {
    if ( ![super isEqualToOperation: otherOperation] )  {
        return NO;
    }
    
    GLVectorScalarOperation * op = otherOperation;
    if (self.scalarOperand != op.scalarOperand) {
        return NO;
    }
    
    return YES;
}

@end


/************************************************/
/*		GLScalarMultiplyOperation				*/
/************************************************/

@implementation GLScalarMultiplyOperation

- (id) initWithVectorOperand: (GLTensor *) vOperand scalarOperand: (GLFloat) sOperand
{
	if (( self = [super initWithVectorOperand: vOperand scalarOperand: sOperand] ))
	{
		GLTensor *resultVariable = self.result[0];
		GLTensor *operandVariable = self.operand[0];
		
		resultVariable.isPurelyReal = operandVariable.isPurelyReal;
		resultVariable.isPurelyImaginary = operandVariable.isPurelyImaginary;
		
		NSUInteger numElements = resultVariable.nDataElements;
		GLFloat localScalar = self.scalarOperand;
		self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			NSMutableData *result = resultArray[0];
			NSMutableData *operand = operandArray[0];
			vGL_vsmul( operand.bytes, 1, &localScalar, result.mutableBytes, 1, numElements);
		};
        self.graphvisDescription = [NSString stringWithFormat: @"scalar multiply %+1.1g", sOperand];
	}
	
    return self;
}

- (BOOL) canOperateInPlace {
	return YES;
}

@end

/************************************************/
/*		GLScalarAddOperation					*/
/************************************************/

@implementation GLScalarAddOperation

- (id) initWithVectorOperand: (GLTensor *) vOperand scalarOperand: (GLFloat) sOperand
{
	if (( self = [super initWithVectorOperand: vOperand scalarOperand: sOperand] ))
	{
		GLTensor *resultVariable = self.result[0];
		GLTensor *operandVariable = self.operand[0];
		
		resultVariable.isPurelyReal = operandVariable.isPurelyReal;
		resultVariable.isPurelyImaginary = NO; // because we just added a (hopefully nonzero) real number
		
		NSUInteger numPoints = resultVariable.nDataPoints;
		GLFloat localScalar = self.scalarOperand;
		self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			NSMutableData *result = resultArray[0];
			NSMutableData *operand = operandArray[0];
			vGL_vsadd( (void *) operand.bytes, 1, (void *) &localScalar, result.mutableBytes, 1, numPoints);
		};
        self.graphvisDescription = [NSString stringWithFormat: @"scalar add %+1.1g", sOperand];
	}
	
    return self;
}

- (BOOL) canOperateInPlace {
	return YES;
}

@end

/************************************************/
/*		GLScalarDivideOperation					*/
/************************************************/

@interface GLScalarDivideOperation ()
@property(readwrite) BOOL canOperateInPlace;
@end

@implementation GLScalarDivideOperation
@synthesize canOperateInPlace;
- (id) initWithVectorOperand: (GLTensor *) vOperand scalarOperand: (GLFloat) sOperand
{
    GLTensor *result;
    if (vOperand.rank == 0) {
        GLScalar *scalar = (GLScalar *) vOperand;
        result = [[GLScalar alloc] initWithType: scalar.dataFormat forEquation:scalar.equation];
    } else if (vOperand.rank == 1) {
        GLVariable *function = (GLVariable *) vOperand;
        result = [[function class] variableOfType: function.dataFormat withDimensions: function.dimensions forEquation: function.equation];
    }  else if (vOperand.rank == 2) {
        GLLinearTransform *aTransform = (GLLinearTransform *) vOperand;
        GLMatrixDescription *matrix = aTransform.matrixDescription;
        for (NSUInteger i=0; i<matrix.nDimensions; i++) {
            if (matrix.strides[i].format != kGLDiagonalMatrixFormat) {
                [NSException raise: @"BadFormat" format: @"This operation type can only transform with matrices in a diagonal format."];
            }
        }
        
        result = [GLLinearTransform transformOfType: aTransform.dataFormat withFromDimensions: aTransform.toDimensions toDimensions:aTransform.fromDimensions inFormat:aTransform.matrixFormats forEquation:aTransform.equation matrix: nil];
    }
	
	if (( self = [super initWithResult: @[result] operand: @[vOperand]] ))
	{
        self.scalarOperand = sOperand;
        
		GLVariable *resultVariable = self.result[0];
		GLVariable *operandVariable = self.operand[0];
        
		resultVariable.isPurelyReal = operandVariable.isPurelyReal;
		resultVariable.isPurelyImaginary = operandVariable.isPurelyImaginary;
				
		if (operandVariable.isComplex) {
			NSUInteger nPoints = operandVariable.nDataPoints;
			if (operandVariable.isPurelyReal) {
				self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
					GLSplitComplex resultComplex = splitComplexFromData(resultArray[0]);
					GLSplitComplex operandComplex = splitComplexFromData(operandArray[0]);
					GLFloat localScalar = sOperand;
					
					vGL_svdiv( (void *) &localScalar, operandComplex.realp, 1, resultComplex.realp, 1, nPoints);
					vGL_vclr( resultComplex.imagp, 1, nPoints);
				};
				self.canOperateInPlace = YES;
                self.graphvisDescription = [NSString stringWithFormat: @"scalar divide %+1.1g (complex, purely real)", sOperand];
			}  else if (operandVariable.isPurelyImaginary) {
				self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
					GLSplitComplex resultComplex = splitComplexFromData(resultArray[0]);
					GLSplitComplex operandComplex = splitComplexFromData(operandArray[0]);
					GLFloat localScalar = -sOperand;
					
					vGL_svdiv( &localScalar, operandComplex.imagp, 1, resultComplex.imagp, 1, nPoints);
					vGL_vclr( resultComplex.realp, 1, nPoints);
				};
				self.canOperateInPlace = YES;
                self.graphvisDescription = [NSString stringWithFormat: @"scalar divide %+1.1g (complex, purely imaginary)", sOperand];
			} else {
				self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
					GLSplitComplex resultComplex = splitComplexFromData(resultArray[0]);
					GLSplitComplex operandComplex = splitComplexFromData(operandArray[0]);
					GLFloat localScalar = sOperand;
					
					vGL_zvmags( &operandComplex, 1, resultComplex.imagp, 1, nPoints );
					// resultComplex.imagp contains (a^2+b^2)
					
					localScalar = 1/localScalar;
					vGL_vsmul( resultComplex.imagp, 1, &localScalar, resultComplex.imagp, 1, nPoints);
					// resultComplex.imagp now contains (a^2+b^2)/k
					
					// Takes C = B/A
					vGL_vdiv( resultComplex.imagp, 1, operandComplex.realp, 1, resultComplex.realp, 1, nPoints );
					// resultComplex.realp = k*a/(a^2+b^2)
					
					vGL_vdiv( resultComplex.imagp, 1, operandComplex.imagp, 1, resultComplex.imagp, 1, nPoints );
					// resultComplex.imagp = k*b/(a^2+b^2)
					
					vGL_vneg( resultComplex.imagp, 1, resultComplex.imagp, 1, nPoints);
					// resultComplex.imagp = -k*b/(a^2+b^2)
				};
				self.canOperateInPlace = NO;
                self.graphvisDescription = [NSString stringWithFormat: @"scalar divide %+1.1g (complex)", sOperand];
			}
		} else {
			NSUInteger nDataElements = resultVariable.nDataElements;
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				NSMutableData *result = resultArray[0];
				NSMutableData *operand = operandArray[0];
				GLFloat localScalar = sOperand;
				vGL_svdiv( &localScalar, (void *) operand.bytes, 1, result.mutableBytes, 1, nDataElements);
			};
			self.canOperateInPlace = YES;
            self.graphvisDescription = [NSString stringWithFormat: @"scalar divide %+1.1g (real)", sOperand];
		}
		
		
	}
	
    return self;
}

@end

/************************************************/
/*		GLPowerOperation						*/
/************************************************/

@implementation GLPowerOperation

- (id) initWithVectorOperand: (GLTensor *) vOperand scalarOperand: (GLFloat) sOperand
{
	if ( [[vOperand class] isSubclassOfClass: [GLLinearTransform class]] ) {
		[NSException raise: @"InvalidOperation" format: @"Cannot take the power of a matrix, yet."];
	}
	
	if (( self = [super initWithVectorOperand: vOperand scalarOperand: sOperand] ))
	{
		GLVariable *resultVariable = self.result[0];
		NSUInteger nDataElements = resultVariable.nDataElements;
		
		GLFloat localScalar = self.scalarOperand;
		self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			NSMutableData *result = resultArray[0];
			NSMutableData *operand = operandArray[0];
            GLFloat *inVar = (GLFloat *) operand.bytes;
            GLFloat *outVar = result.mutableBytes;
            for (NSUInteger i=0; i<nDataElements; i++) {
                outVar[i] = pow( inVar[i], localScalar );
            }
		};
		self.graphvisDescription = [NSString stringWithFormat: @"power %+1.1g", sOperand];
		
	}
	
	return self;
}

- (BOOL) canOperateInPlace {
	return YES;
}

@end

/************************************************/
/*		GLScalarThresholdOperation				*/
/************************************************/
// variable = min( operand, scalar )
@implementation GLScalarThresholdOperation

- (id) initWithVectorOperand: (GLTensor *) vOperand scalarOperand: (GLFloat) sOperand
{
	if (( self = [super initWithVectorOperand: vOperand scalarOperand: sOperand] ))
	{
		GLVariable *resultVariable = self.result[0];
		NSUInteger nDataElements = resultVariable.nDataElements;
		GLFloat localScalar = self.scalarOperand;
		self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			NSMutableData *result = resultArray[0];
			NSMutableData *operand = operandArray[0];
			vGL_vthr( (void *) operand.bytes, 1, (void *) &localScalar, result.mutableBytes, 1, nDataElements);
		};
        self.graphvisDescription = [NSString stringWithFormat: @"min(operand, %+1.1g)", sOperand];
	}

	return self;
}

- (BOOL) canOperateInPlace {
	return YES;
}

@end

/************************************************/
/*		GLZeroThresholdOperation				*/
/************************************************/
// variable = operand > scalar ? operand : 0.0
@implementation GLZeroThresholdOperation

- (id) initWithVectorOperand: (GLTensor *) vOperand scalarOperand: (GLFloat) sOperand
{
	if (( self = [super initWithVectorOperand: vOperand scalarOperand: sOperand] ))
	{
		GLVariable *resultVariable = self.result[0];
		NSUInteger nDataElements = resultVariable.nDataElements;
		GLFloat localScalar = self.scalarOperand;
		self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			NSMutableData *result = resultArray[0];
			NSMutableData *operand = operandArray[0];
			vGL_vthres( (void *) operand.bytes, 1, (void *) &localScalar, result.mutableBytes, 1, nDataElements);
		};
        self.graphvisDescription = [NSString stringWithFormat: @"threshold: operand > %+1.1g ? operand : 0.0", sOperand];
	}
	
	return self;
}

- (BOOL) canOperateInPlace {
	return YES;
}

@end

/************************************************/
/*		GLClipOperation				*/
/************************************************/
// variable = max( operand, scalar )
@implementation GLClipOperation 

- (id) initWithVectorOperand: (GLTensor *) vOperand firstScalarOperand: (GLFloat) fsOperand secondScalarOperand: (GLFloat) ssOperand
{
	if (( self = [super initWithOperand: @[vOperand]] ))
	{
        self.firstScalarOperand = fsOperand;
        self.secondScalarOperand = ssOperand;
		GLVariable *resultVariable = self.result[0];
		NSUInteger nDataElements = resultVariable.nDataElements;
		GLFloat min = fsOperand;
		GLFloat max = ssOperand;
		self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			NSMutableData *result = resultArray[0];
			NSMutableData *operand = operandArray[0];
			vGL_vclip( (void *) operand.bytes, 1, (void *) &min, (void *) &max, result.mutableBytes, 1, nDataElements);
		};
        self.graphvisDescription = [NSString stringWithFormat: @"clip (%+1.1g, %+1.1g)", fsOperand, ssOperand];
	}
	return self;
}

- (BOOL) isEqualToOperation: (id) otherOperation {
    if ( ![super isEqualToOperation: otherOperation] )  {
        return NO;
    }
    
    GLClipOperation * op = otherOperation;
    if (self.firstScalarOperand != op.firstScalarOperand) {
        return NO;
    }
    if (self.secondScalarOperand != op.secondScalarOperand) {
        return NO;
    }
    
    return YES;
}

@end

/************************************************/
/*		GLSetValueOperation						*/
/************************************************/

@implementation GLSetValueOperation

- (id) initWithVectorOperand: (GLTensor *) aTensor scalarOperand: (GLFloat) aScalar indexString: (NSString *) indexString
{
	if ( ![[aTensor class] isSubclassOfClass: [GLVariable class]] ) {
		[NSException raise: @"InvalidOperation" format: @"Can only set the value of a function."];
	}
	
	GLVariable *variable = (GLVariable *) aTensor;
    NSArray *ranges = [GLDimension rangesFromIndexString: indexString usingDimensions: variable.dimensions];
    
    if (( self = [super initWithVectorOperand: variable scalarOperand: aScalar] ))
	{
		GLVariable *resultVariable = self.result[0];
		GLVariable *operandVariable = self.operand[0];

        self.indexString = indexString;
        
        NSUInteger numBytes = resultVariable.nDataElements*sizeof(GLFloat);
        
        resultVariable.name = variable.name;
		resultVariable.isPurelyReal = operandVariable.isPurelyReal;
		resultVariable.isPurelyImaginary = operandVariable.isPurelyImaginary;
        
		self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			NSMutableData *result = resultArray[0];
			NSMutableData *operand = operandArray[0];
			memcpy( result.mutableBytes, operand.bytes,  numBytes );
		};
        
        if ( variable.dimensions.count == 1 )
		{
            NSRange fastRange = [ranges[0] rangeValue];
            
            NSUInteger startIndex = fastRange.location;
            NSUInteger fastIndexLength = fastRange.length;
            
            GLFloat scalarValue = aScalar;
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				NSMutableData *result = resultArray[0];
				NSMutableData *operand = operandArray[0];
                GLFloat *toData = (GLFloat *) result.mutableBytes;
                // first copy the data
                memcpy( result.mutableBytes, operand.bytes,  numBytes );
                // then replace the value at the desired indices
                vGL_vfill( (GLFloat *) &scalarValue, &toData[startIndex], 1, fastIndexLength);
                
            };
            self.graphvisDescription = [NSString stringWithFormat: @"set=%+1.1g (1 dim)", aScalar];
		}
		else if ( variable.dimensions.count == 2 )
		{
			NSRange fastRange = [ranges[1] rangeValue];
            NSUInteger fastDimLength = [variable.dimensions[1] nPoints];
            
            NSRange slowRange = [ranges[0] rangeValue];
            			
            GLFloat scalarValue = aScalar;
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				NSMutableData *result = resultArray[0];
				NSMutableData *operand = operandArray[0];
				// first copy the data
				memcpy( result.mutableBytes, operand.bytes,  numBytes );
				
                GLFloat *toData = (GLFloat *) result.mutableBytes;
                dispatch_apply(slowRange.length, dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0), ^(size_t iteration) {
                    // then replace the value at the desired indices
                    vGL_vfill( (GLFloat *)&scalarValue, &(toData[(slowRange.location + iteration)*fastDimLength + fastRange.location]), 1, fastRange.length);
                    
                });
            };
            self.graphvisDescription = [NSString stringWithFormat: @"set=%+1.1g (2 dim)", aScalar];
        }
        else if ( variable.dimensions.count == 3 )
        {            
            NSUInteger ny = [variable.dimensions[1] nPoints];
            NSUInteger nz = [variable.dimensions[2] nPoints];
            
            NSRange xrange = [ranges[0] rangeValue];
            NSRange yrange = [ranges[1] rangeValue];
            NSRange zrange = [ranges[2] rangeValue];
            
			self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				NSMutableData *result = resultArray[0];
				NSMutableData *operand = operandArray[0];
				// first copy the data
				memcpy( result.mutableBytes, operand.bytes,  numBytes );
				
                GLFloat *toData = (GLFloat *) result.mutableBytes;
                for (NSUInteger i=xrange.location; i<xrange.location+xrange.length; i++) {
                    for (NSUInteger j=yrange.location; j<yrange.location+yrange.length; j++) {
                        for (NSUInteger k=zrange.location; k<zrange.location+zrange.length; k++) {
                            toData[(i*ny+j)*nz+k] = aScalar;
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
    
    GLSetValueOperation * op = otherOperation;
    if ( [self.indexString isEqualToString: op.indexString]) {
        return NO;
    }
    
    return YES;
}

@end

/************************************************/
/*		GLScaleOperation						*/
/************************************************/

@implementation GLScaleOperation

-(id) initWithOperand: (GLTensor *) aTensor variableScale: (GLFloat) sOperand units: (NSString *) varUnits dimensionScale: (GLFloat) dimScale translation: (GLFloat) delta units: (NSString *) dimUnits
{
	if ( ![[aTensor class] isSubclassOfClass: [GLVariable class]] ) {
		[NSException raise: @"InvalidOperation" format: @"Can only set the value of a function."];
	}
	
	GLVariable *vOperand = (GLVariable *) aTensor;
	NSMutableArray *newDims = [NSMutableArray arrayWithCapacity: vOperand.dimensions.count];
	for (GLDimension *dim in vOperand.dimensions) {
		[newDims addObject: [dim scaledBy: dimScale translatedBy: delta withUnits: dimUnits]];
	}
	GLVariable *result = [[vOperand class] variableOfType: vOperand.dataFormat withDimensions: newDims forEquation: vOperand.equation];
	result.units = varUnits;
	if (( self = [super initWithResult: @[result] operand: @[vOperand]] ))
	{
		GLVariable *resultVariable = self.result[0];
		GLVariable *operandVariable = self.operand[0];
		
		self.scalarOperand = sOperand;
		self.dimScale = dimScale;
		self.dimTranslation = delta;
		resultVariable.isPurelyReal = operandVariable.isPurelyReal;
		resultVariable.isPurelyImaginary = operandVariable.isPurelyImaginary;
		
		NSUInteger numElements = resultVariable.nDataElements;
		GLFloat localScalar = self.scalarOperand;
		self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			NSMutableData *result = resultArray[0];
			NSMutableData *operand = operandArray[0];
			vGL_vsmul( operand.bytes, 1, &localScalar, result.mutableBytes, 1, numElements);
		};
        self.graphvisDescription = [NSString stringWithFormat: @"dimensionalize by %+1.1g", sOperand];
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
    
    GLScaleOperation * op = otherOperation;
    if ( self.scalarOperand != op.scalarOperand) {
        return NO;
    } else if ( self.dimScale != op.dimScale) {
        return NO;
    } else if ( self.dimTranslation != op.dimTranslation) {
        return NO;
    }
    
    return YES;
}

@end