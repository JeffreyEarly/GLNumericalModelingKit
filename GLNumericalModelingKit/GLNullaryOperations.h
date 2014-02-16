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

- (id) initWithResult: (GLFunction *) result firstScalarOperand: (GLFloat) fsOperand secondScalarOperand: (GLFloat) ssOperand;

@property GLFloat firstScalarOperand;
@property GLFloat secondScalarOperand;

@end

/************************************************/
/*		GLNormalDistributionOperation			*/
/************************************************/

/** C = randn. This randomly generates numbers with a standard normal distribution
 @discussion If the function is in half complex format, this operation will make it hermitian, including a doubling of the four purely real components.
 @discussion mean + expectation*randn will scale this correctly.
 @returns A GLNormalDistributionOperation object.
 */
@interface GLNormalDistributionOperation : GLVariableOperation
- (id) initWithResult: (GLFunction *) result;
- (id) initWithResult: (GLFunction *) result seed: (GLScalar *) seed;
@end