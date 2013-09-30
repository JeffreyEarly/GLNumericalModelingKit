//
//  GLSpectralDifferentialOperator.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 4/22/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLDifferentialOperator.h>
#import <GLNumericalModelingKit/GLVariable.h>

// A spectral differentiation operator is just complex variable that lives in fourier space.
// To do its differentiation, it just multiplies!
@interface GLSpectralDifferentialOperator : GLDifferentialOperator {

}

@property(readonly, nonatomic) BOOL basisTransformationRequired;

@end
