//
//  GLInterpolationOperations.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 5/2/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import <GLNumericalModelingKit/GLVariableOperations.h>

/************************************************/
/*		GLInterpolationOperation				*/
/************************************************/

// variable = leftVariable interpolated at points rightVariable
// The result variable and the right variable will be of the same dimension.
// The rightVariable vector must be the same size as the dimensions of the left variable.
// Typical usage:
//
// Find the value of the two variables u and v at the positions defined by the vector (xPos, yPos)
// GLInterpolationOperation *interp = [[GLInterpolationOperation alloc] initWithFirstOperand: @[u,v] secondOperand: @[xPos, yPos]];
//
// You can interpolate as many functions as desired at this points, so
// GLInterpolationOperation *interp = [[GLInterpolationOperation alloc] initWithFirstOperand: @[u,v,rv,pv] secondOperand: @[xPos, yPos]];
// is just as valid.
//
// The guiding assumption in the above is that all the functions u,v,rv,pv have two dimensions.
//
// This also works fine for one-dimensional functions,
// GLInterpolationOperation *interp = [[GLInterpolationOperation alloc] initWithFirstOperand: @[f1,f2] secondOperand: @[pos]];
//
@interface GLInterpolationOperation : GLBinaryVectorOperation
@end