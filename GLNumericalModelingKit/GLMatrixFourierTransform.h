//
//  GLMatrixFourierTransform.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 10/11/11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#import "GLFourierTransform.h"

@interface GLMatrixFourierTransform : GLFourierTransform
{
	NSArray *dimensions;
	NSMutableData *zerosBuffer;
	NSUInteger dataPoints;
}

- (GLFourierTransform *) initWithDimensions: (NSArray *) dimensions;

- (void) transform: (GLFloat *) f forward: (GLSplitComplex *) fbar;
- (void) transform: (GLSplitComplex *) fbar inverse: (GLFloat *) f;

@property(readonly, strong, nonatomic) NSArray *dimensions;
@property(readonly, strong, nonatomic) NSMutableData *zerosBuffer;
@property(readonly, assign, nonatomic) NSUInteger dataPoints;


@end
