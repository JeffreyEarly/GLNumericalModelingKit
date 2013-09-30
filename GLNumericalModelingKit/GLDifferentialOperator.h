//
//  GLDifferentialOperator.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 4/27/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLVariable.h>
#import <GLNumericalModelingKit/GLVariableOperations.h>

/************************************************/
/*		GLDifferentialOperator					*/
/************************************************/

#pragma mark -
#pragma mark GLDifferentialOperator
#pragma mark

// A differential operator can be implemented in any number of ways, but at it's most basic level
// it just takes a variable and returns a differentiated version.

@interface GLDifferentialOperator : GLVariable

@property(readwrite, strong) NSArray *fromDimensions;
@property(readwrite, strong) NSArray *toDimensions;

// Returns a fully initialized differentiation operation that will differentiate the operand.
- (GLBinaryOperation *) differentiationOperationFromVariable: (GLVariable *) operand;

// Returns the result variable of the differentiation operation.
- (GLVariable *) differentiateVariable: (GLVariable *) operand;

@end
