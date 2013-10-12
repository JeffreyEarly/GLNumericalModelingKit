//
//  GLNullaryOperations.h
//  GLNumericalModelingKitCopy
//
//  Created by Jeffrey J. Early on 11/27/12.
//
//

#import <GLNumericalModelingKit/GLVariableOperations.h>

/************************************************/
/*		GLRandomNumberOperation					*/
/************************************************/
// variable is populated with random numbers between the two values.
@interface GLRandomNumberOperation : GLVariableOperation

- (id) initWithResult: (GLVariable *) result firstScalarOperand: (GLFloat) fsOperand secondScalarOperand: (GLFloat) ssOperand;

@property GLFloat firstScalarOperand;
@property GLFloat secondScalarOperand;

@end
