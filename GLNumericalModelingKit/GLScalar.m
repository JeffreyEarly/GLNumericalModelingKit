//
//  GLScalar.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 10/11/13.
//  Copyright (c) 2013 Jeffrey J. Early. All rights reserved.
//

#import "GLScalar.h"
#import "GLVectorVectorOperations.h"
#import "GLUnaryOperations.h"

@interface GLScalar ()
@property(readwrite, assign, nonatomic) NSUInteger nDataPoints;
@property(readwrite, assign, nonatomic) NSUInteger nDataElements;
@property(readwrite, assign, nonatomic) NSUInteger dataBytes;
@end

@implementation GLScalar

@synthesize nDataPoints = _nDataPoints;
@synthesize nDataElements = _nDataElements;
@synthesize dataBytes = _dataBytes;

+ (GLScalar *) scalarWithValue: (GLFloatComplex) aValue forEquation: (GLEquation *) anEquation;
{
	return [[GLScalar alloc] initWithValue: aValue forEquation: anEquation];
}

+ (GLScalar *) scalarWithType: (GLDataFormat) format forEquation: (GLEquation *) anEquation;
{
	return [[GLScalar alloc] initWithType: format forEquation: anEquation];
}

- (GLScalar *) initWithType: (GLDataFormat) format forEquation: (GLEquation *) anEquation
{
	if ((self = [super initWithType: format withEquation:anEquation])) {
		_nDataPoints = 1;
        if (self.dataFormat == kGLSplitComplexDataFormat || self.dataFormat == kGLInterleavedComplexDataFormat) {
            _nDataElements = 2;
        } else {
			_nDataElements = 1;
		}
		_dataBytes = _nDataElements*sizeof(GLFloat);
	}
	
	return self;
}

- (GLScalar *) initWithValue: (GLFloatComplex) aValue forEquation: (GLEquation *) anEquation
{
	GLDataFormat format = cimag(aValue) == 0.0 ? kGLRealDataFormat : kGLSplitComplexDataFormat;
	
	if ((self = [self initWithType: format forEquation: anEquation])) {
		if (format == kGLSplitComplexDataFormat) {
			self.splitComplex.realp[0] = creal(aValue);
			self.splitComplex.imagp[0] = cimag(aValue);
		} else {
			self.pointerValue[0] = creal(aValue);
		}
	}
	
	return self;
}

- (NSUInteger) rank {
	return 0;
}

- (id) dividedBy: (id) otherVariable
{
	GLDivisionOperation *operation = [[GLDivisionOperation alloc] initWithFirstOperand: self secondOperand: otherVariable];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result[0];
}

- (NSString *) graphvisDescription
{
    NSMutableString *extra = [NSMutableString stringWithFormat: @""];
    if (self.name) [extra appendFormat: @"%@:", self.name];
	[extra appendString: self.isComplex ? @"complex" : @"real"];
	if (self.isRealPartZero) [extra appendString:@", zero real part"];
    if (self.isComplex && self.isImaginaryPartZero) [extra appendString:@", zero imaginary part"];
	[extra appendString: @"scalar"];
    
    return extra;
}

- (GLScalar *) exponentiate
{
	GLExponentialOperation *operation = [[GLExponentialOperation alloc] initWithVariable: self];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result[0];
}

- (GLScalar *) log
{
	GLLogarithmOperation *operation = [[GLLogarithmOperation alloc] initWithVariable: self];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result[0];
}

- (GLScalar *) sin
{
	GLSineOperation *operation = [[GLSineOperation alloc] initWithVariable: self];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result[0];
}

- (GLScalar *) cos
{
	GLCosineOperation *operation = [[GLCosineOperation alloc] initWithVariable: self];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result[0];
}

- (GLScalar *) atan
{
	GLInverseTangentOperation *operation = [[GLInverseTangentOperation alloc] initWithVariable: self];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result[0];
}

- (GLScalar *) sqrt
{
	GLSquareRootOperation *operation = [[GLSquareRootOperation alloc] initWithVariable: self];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result[0];
}

- (GLScalar *) negate
{
	GLNegationOperation *operation = [[GLNegationOperation alloc] initWithVariable: self];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result[0];
}

- (GLScalar *) abs
{
	GLAbsoluteValueOperation *operation = [[GLAbsoluteValueOperation alloc] initWithVariable: self];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result[0];
}

@end
