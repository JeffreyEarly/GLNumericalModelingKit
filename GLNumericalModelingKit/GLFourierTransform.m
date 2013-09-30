//
//  GLFourierTransform.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 4/18/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import "GLFourierTransform.h"

//#import "GLDSPFourierTransform.h"
//#import "GLMatrixFourierTransform.h"
//#import "GLOpenCLFourierTransform.h"
#import "GLFFTWFourierTransform.h"

@implementation GLFourierTransform

+ (GLFourierTransform *) fourierTransformWithDimensions: (NSArray *) dimensions
{
//	return [[GLOpenCLFourierTransform alloc] initWithDimensions: dimensions];
//	return [[GLMatrixFourierTransform alloc] initWithDimensions: dimensions];
	return [[GLFFTWFourierTransform alloc] initWithDimensions: dimensions];
//	return [[GLDSPFourierTransform alloc] initWithDimensions: dimensions];
}

- (void) transform: (GLFloat *) f forward: (GLSplitComplex *) fbar
{

}

- (void) transform: (GLSplitComplex *) fbar inverse: (GLFloat *) f
{

}

- (NSArray *) dimensions
{
	return [[NSArray alloc] init];
}

@end
