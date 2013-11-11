//
//  GLNetCDFConcatenationOperation.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 3/13/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import <GLNumericalModelingKit/GLVariableOperations.h>

@class GLMutableNetCDFVariable;

// This operation mutates the firstOperand but concatenating the data from the second operand.
// The mutable dimension grows in accordance.
@interface GLNetCDFConcatenationOperation :  GLVariableOperation

// The first operand and the second operand are both n-dimensional.
// The dimensions at indices other than dimIndex must have the same number of points.
// The dimensions at dimIndex do not need to match. If the first operand's dimension is evenly spaced
// then it will be extended, otherwise the values from the second operand's dimension will be added.
- (id) initWithFirstOperand: (GLMutableNetCDFVariable *) fOperand secondOperand: (GLVariable *) sOperand alongDimensionAtIndex: (NSUInteger) mutableDimensionIndex;

// The first operand is n-dimensional and the second operand is (n-1)-dimensional.
// The (n-1) dimensions of the two variables must have the same number of points.
// If the mutableDimension is evenly spaced, then it will be extended to length pointIndex+1, if necessary.
// If the mutableDimension is not evenly spaced, then it must already have the correct value.
- (id) initWithFirstOperand: (GLMutableNetCDFVariable *) fOperand lowerDimensionalSecondOperand: (GLVariable *) sOperand alongDimensionAtIndex: (NSUInteger) mutableDimensionIndex index: (NSUInteger) pointIndex;

@property(strong) NSMutableArray *indexRanges;

@end
