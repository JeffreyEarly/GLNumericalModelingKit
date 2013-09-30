//
//  GLFourierTransform.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 4/18/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/Precision.h>

@interface GLFourierTransform : NSObject

+ (GLFourierTransform *) fourierTransformWithDimensions: (NSArray *) dimensions;

- (void) transform: (GLFloat *) f forward: (GLSplitComplex *) fbar;
- (void) transform: (GLSplitComplex *) fbar inverse: (GLFloat *) f;

- (NSArray *) dimensions;

@end
