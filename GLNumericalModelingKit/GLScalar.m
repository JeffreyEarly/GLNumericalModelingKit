//
//  GLScalar.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 10/11/13.
//  Copyright (c) 2013 Jeffrey J. Early. All rights reserved.
//

#import "GLScalar.h"

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
		self.splitComplex.realp[0] = creal(aValue);
		self.splitComplex.imagp[0] = cimag(aValue);
	}
	
	return self;
}

- (NSUInteger) rank {
	return 0;
}

@end
