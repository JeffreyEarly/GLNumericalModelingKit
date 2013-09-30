//
//  GLSpectralDifferentialOperatorPool.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 4/26/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "GLDifferentialOperatorPool.h"
#import "GLSpectralDifferentialOperator.h"

// This class overrides -setDimensions to make sure they're in the correct domain.

// This is a bad idea. These differential operators truly can operate on two different sets of dimensions
// and this should be made clear.

@interface GLSpectralDifferentialOperatorPool : GLDifferentialOperatorPool

- (GLSpectralDifferentialOperator *) operatorWithDerivatives: (NSArray *) derivatives;

- (GLSpectralDifferentialOperator *) harmonicOperator;

- (GLSpectralDifferentialOperator *) harmonicOperatorOfOrder: (NSUInteger) order;

// The cutoff should be between 0 and 1.
- (GLSpectralDifferentialOperator *) spectralVanishingViscosityFilter;

// These are thrown in here to surpress compiler warnings these (and higher order derivatives)
// will be resolved dynam/Users/jearly/Dropbox/Documents/Manuscriptsically at runtime if they exist (given the dimensions).
@property(readonly) GLSpectralDifferentialOperator* x;
@property(readonly) GLSpectralDifferentialOperator* y;
@property(readonly) GLSpectralDifferentialOperator* xx;
@property(readonly) GLSpectralDifferentialOperator* xy;
@property(readonly) GLSpectralDifferentialOperator* yy;
@property(readonly) GLSpectralDifferentialOperator* xxx;
@property(readonly) GLSpectralDifferentialOperator* xxy;
@property(readonly) GLSpectralDifferentialOperator* xyy;
@property(readonly) GLSpectralDifferentialOperator* yyy;

@end


@interface GLVariable (SpectralDifferentiationExtensions)

@property(readonly) GLVariable* x;
@property(readonly) GLVariable* y;
@property(readonly) GLVariable* xx;
@property(readonly) GLVariable* xy;
@property(readonly) GLVariable* yy;
@property(readonly) GLVariable* xxx;
@property(readonly) GLVariable* xxy;
@property(readonly) GLVariable* xyy;
@property(readonly) GLVariable* yyy;

@end