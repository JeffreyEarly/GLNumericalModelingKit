//
//  GLDSPFourierTransform.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 10/11/11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/Precision.h>
#import "GLFourierTransform.h"

@interface GLDSPFourierTransform : GLFourierTransform
{
	NSArray *dimensions;
	NSMutableData *zerosBuffer;
	GLFFTSetup setup;
	NSMutableArray *logDimensions;
	NSUInteger dataPoints;
}

- (GLFourierTransform *) initWithDimensions: (NSArray *) dimensions;

- (void) transform: (GLFloat *) f forward: (GLSplitComplex *) fbar;
- (void) transform: (GLSplitComplex *) fbar inverse: (GLFloat *) f;

@property(strong) NSArray *dimensions;
@property(strong) NSMutableData *zerosBuffer;
@property GLFFTSetup setup;
@property(strong) NSMutableArray *logDimensions;
@property NSUInteger dataPoints;

@end
