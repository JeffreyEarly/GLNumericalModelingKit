//
//  GLFFTWFourierTransform.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 5/16/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import "GLFourierTransform.h"

@interface GLFFTWFourierTransform : GLFourierTransform

- (GLFourierTransform *) initWithDimensions: (NSArray *) dimensions;

- (void) transform: (GLFloat *) f forward: (GLSplitComplex *) fbar;
- (void) transform: (GLSplitComplex *) fbar inverse: (GLFloat *) f;

@property(readwrite, strong, nonatomic) NSArray *dimensions;
@property(readwrite, strong, nonatomic) NSMutableData *zerosBuffer;
@property(readwrite, assign, nonatomic) NSUInteger dataPoints;


@end
