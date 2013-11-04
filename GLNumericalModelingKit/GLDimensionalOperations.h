//
//  GLDimensionalOperations.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 1/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import <GLNumericalModelingKit/GLVectorVectorOperations.h>

/************************************************/
/*		GLAddDimensionOperation                 */
/************************************************/

// Returns nil if the variable already has that dimension.
@interface GLAddDimensionOperation : GLVariableOperation

- (id) initWithOperand: (GLVariable *) variable dimension: (GLDimension *) dim;

// If the dimension is mutable, we need to log the number of points at the time of creation.
@property(readwrite, strong, nonatomic) GLDimension *theDimension;
@property NSUInteger nPoints;

@end



/************************************************/
/*		GLSubdomainOperation					*/
/************************************************/

// Returns a variable with only the elements indicated by by the array of ranges.
// The size of the ranges array must match the number of dimensions.
// If the -shouldFlatten flag is set, this will cause any dimensions of length 1 to be eliminated.
@interface GLSubdomainOperation : GLVariableOperation

- (id) initWithOperand: (GLVariable *) variable indexRange: (NSArray *) ranges flatten: (BOOL) aFlag;

@property(readwrite, strong, nonatomic) NSArray *theRanges;
@property(readwrite, assign, nonatomic) BOOL shouldFlatten;


@end

/****************************************************/
/*		GLExistingDimensionConcatenationOperation	*/
/****************************************************/

// This operation is minimilly restrictive.
// You can concatenate an n dimensional variable with an n or n-1 dimensional variable.
// The dimensionIndex refers to the dimensions of the firstOperand.
// The dimension of concatenation must be evenly sampled.
// The other dimensions must have the same number of points.
// 
// The result variable will have dimensions based on the recievers dimensions.
@interface GLExistingDimensionConcatenationOperation : GLVariableOperation

- (id) initWithFirstOperand: (GLVariable *) fOperand secondOperand: (GLVariable *) sOperand dimensionIndex: (NSUInteger) dimIndex;

@property(readwrite, strong, nonatomic) GLDimension *theDimension;

@end

/************************************************/
/*		GLNewDimensionConcatenationOperation    */
/************************************************/

// The new dimension should have nPoints=2, otherwise a new dimension based on the same properties will be created.
// The dimensions in the other two variables must be the same.
@interface GLNewDimensionConcatenationOperation : GLVariableOperation

- (id) initWithFirstOperand: (GLVariable *) fOperand secondOperand: (GLVariable *) sOperand dimension: (GLDimension *) dim;

@end

/************************************************/
/*		GLHalfToSplitComplexOperation			*/
/************************************************/
@interface GLHalfToSplitComplexOperation : GLVariableOperation

@end

/************************************************/
/*		GLZeroPadOperation                      */
/************************************************/

@interface GLZeroPadOperation : GLVariableOperation

- (id) initWithOperand: (GLVariable *) variable newDimensions: (NSArray *) newDimensions basis: (NSArray *) basis;
@property(readwrite, strong, nonatomic) NSArray *projectedDimensions;
@end
