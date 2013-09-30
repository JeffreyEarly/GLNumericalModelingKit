//
//  GLVectorScalarOperations.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 1/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import "GLVectorScalarOperations.h"

/************************************************/
/*		GLVectorScalarOperation					*/
/************************************************/

@implementation GLVectorScalarOperation

- (id) initWithVectorOperand: (GLVariable *) vOperand scalarOperand: (GLFloat) sOperand
{
	if (( self = [super initWithOperand: vOperand] ))
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

- (id) initWithVectorOperand: (GLVariable *) vOperand scalarOperand: (GLFloat) sOperand
{
	if (( self = [super initWithVectorOperand: vOperand scalarOperand: sOperand] ))
	{
		self.result.isPurelyReal = self.operand.isPurelyReal;
		self.result.isPurelyImaginary = self.operand.isPurelyImaginary;
		
		NSUInteger numElements = self.result.nDataElements;
		GLFloat localScalar = self.scalarOperand;
		self.blockOperation = ^(NSMutableData *result, NSData *operand) {
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

- (id) initWithVectorOperand: (GLVariable *) vOperand scalarOperand: (GLFloat) sOperand
{
	if (( self = [super initWithVectorOperand: vOperand scalarOperand: sOperand] ))
	{
		self.result.isPurelyReal = self.operand.isPurelyReal;
		self.result.isPurelyImaginary = NO; // because we just added a (hopefully nonzero) real number
		
		NSUInteger numPoints = self.result.nDataPoints;
		GLFloat localScalar = self.scalarOperand;
		self.blockOperation = ^(NSMutableData *result, NSData *operand) {
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
- (id) initWithVectorOperand: (GLVariable *) vOperand scalarOperand: (GLFloat) sOperand
{
	if (( self = [super initWithVectorOperand: vOperand scalarOperand: sOperand] ))
	{
		self.result.isPurelyReal = self.operand.isPurelyReal;
		self.result.isPurelyImaginary = self.operand.isPurelyImaginary;
				
		if (self.operand.isComplex) {
			NSUInteger nPoints = self.operand.nDataPoints;
			if (self.operand.isPurelyReal) {
				self.blockOperation = ^(NSMutableData *result, NSData *operand) {
					GLSplitComplex resultComplex = splitComplexFromData(result);
					GLSplitComplex operandComplex = splitComplexFromData(operand);
					GLFloat localScalar = sOperand;
					
					vGL_svdiv( (void *) &localScalar, operandComplex.realp, 1, resultComplex.realp, 1, nPoints);
					vGL_vclr( resultComplex.imagp, 1, nPoints);
				};
				self.canOperateInPlace = YES;
                self.graphvisDescription = [NSString stringWithFormat: @"scalar divide %+1.1g (complex, purely real)", sOperand];
			}  else if (self.operand.isPurelyImaginary) {
				self.blockOperation = ^(NSMutableData *result, NSData *operand) {
					GLSplitComplex resultComplex = splitComplexFromData(result);
					GLSplitComplex operandComplex = splitComplexFromData(operand);
					GLFloat localScalar = -sOperand;
					
					vGL_svdiv( &localScalar, operandComplex.imagp, 1, resultComplex.imagp, 1, nPoints);
					vGL_vclr( resultComplex.realp, 1, nPoints);
				};
				self.canOperateInPlace = YES;
                self.graphvisDescription = [NSString stringWithFormat: @"scalar divide %+1.1g (complex, purely imaginary)", sOperand];
			} else {
				self.blockOperation = ^(NSMutableData *result, NSData *operand) {
					// This seems fairly unoptimized.
					GLSplitComplex resultComplex = splitComplexFromData(result);
					GLSplitComplex operandComplex = splitComplexFromData(operand);
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
			NSUInteger nDataElements = self.result.nDataElements;
			self.blockOperation = ^(NSMutableData *result, NSData *operand) {
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

#warning GLPowerOperation remains unoptimized
- (id) initWithVectorOperand: (GLVariable *) vOperand scalarOperand: (GLFloat) sOperand
{
	if (( self = [super initWithVectorOperand: vOperand scalarOperand: sOperand] ))
	{
		const int nDataElements = (int) self.result.nDataElements;
		GLFloat localScalar = self.scalarOperand;
		self.blockOperation = ^(NSMutableData *result, NSData *operand) {
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

- (id) initWithVectorOperand: (GLVariable *) vOperand scalarOperand: (GLFloat) sOperand
{
	if (( self = [super initWithVectorOperand: vOperand scalarOperand: sOperand] ))
	{
		NSUInteger nDataElements = self.result.nDataElements;
		GLFloat localScalar = self.scalarOperand;
		self.blockOperation = ^(NSMutableData *result, NSData *operand) {
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

- (id) initWithVectorOperand: (GLVariable *) vOperand scalarOperand: (GLFloat) sOperand
{
	if (( self = [super initWithVectorOperand: vOperand scalarOperand: sOperand] ))
	{
		NSUInteger nDataElements = self.result.nDataElements;
		GLFloat localScalar = self.scalarOperand;
		self.blockOperation = ^(NSMutableData *result, NSData *operand) {
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

- (id) initWithVectorOperand: (GLVariable *) vOperand firstScalarOperand: (GLFloat) fsOperand secondScalarOperand: (GLFloat) ssOperand
{
	if (( self = [super initWithOperand: vOperand] ))
	{
        self.firstScalarOperand = fsOperand;
        self.secondScalarOperand = ssOperand;
		NSUInteger nDataElements = self.result.nDataElements;
		GLFloat min = fsOperand;
		GLFloat max = ssOperand;
		self.blockOperation = ^(NSMutableData *result, NSData *operand) {
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

- (id) initWithVectorOperand: (GLVariable *) variable scalarOperand: (GLFloat) aScalar indexString: (NSString *) indexString
{
    NSArray *ranges = [GLDimension rangesFromIndexString: indexString usingDimensions: variable.dimensions];
    
    if (( self = [super initWithVectorOperand: variable scalarOperand: aScalar] ))
	{
        self.indexString = indexString;
        
        NSUInteger numBytes = self.result.nDataElements*sizeof(GLFloat);
        
        self.result.name = variable.name;
		self.result.isPurelyReal = self.operand.isPurelyReal;
		self.result.isPurelyImaginary = self.operand.isPurelyImaginary;
        
		self.blockOperation = ^(NSMutableData *result, NSData *operand) {
			memcpy( result.mutableBytes, operand.bytes,  numBytes );
            
		};
        
        if ( variable.dimensions.count == 1 )
		{
            NSRange fastRange = [ranges[0] rangeValue];
            
            NSUInteger startIndex = fastRange.location;
            NSUInteger fastIndexLength = fastRange.length;
            
            GLFloat scalarValue = aScalar;
            self.blockOperation = ^(NSMutableData *result, NSData *operand) {
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
            self.blockOperation = ^(NSMutableData *result, NSData *operand) {
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
            
            self.blockOperation = ^(NSMutableData *result, NSData *operand) {
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

-(id) initWithOperand: (GLVariable *) vOperand variableScale: (GLFloat) sOperand units: (NSString *) varUnits dimensionScale: (GLFloat) dimScale translation: (GLFloat) delta units: (NSString *) dimUnits
{
	NSMutableArray *newDims = [NSMutableArray arrayWithCapacity: vOperand.dimensions.count];
	for (GLDimension *dim in vOperand.dimensions) {
		[newDims addObject: [dim scaledBy: dimScale translatedBy: delta withUnits: dimUnits]];
	}
	GLVariable *result = [[vOperand class] variableOfType: vOperand.dataFormat withDimensions: newDims forEquation: vOperand.equation];
	result.units = varUnits;
	if (( self = [super initWithResult: result operand: vOperand] ))
	{
		self.scalarOperand = sOperand;
		self.dimScale = dimScale;
		self.dimTranslation = delta;
		self.result.isPurelyReal = self.operand.isPurelyReal;
		self.result.isPurelyImaginary = self.operand.isPurelyImaginary;
		
		NSUInteger numElements = self.result.nDataElements;
		GLFloat localScalar = self.scalarOperand;
		self.blockOperation = ^(NSMutableData *result, NSData *operand) {
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