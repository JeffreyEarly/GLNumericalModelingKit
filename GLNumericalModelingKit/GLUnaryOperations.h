//
//  GLUnaryOperations.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 1/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import <GLNumericalModelingKit/GLVariableOperations.h>

/************************************************/
/*		GLNegationOperation                     */
/************************************************/

// Returns the negative of the variable.

@interface GLNegationOperation : GLVariableOperation
- (GLNegationOperation *) initWithFunction: (GLVariable *) variable;
@end

/************************************************/
/*		GLAbsoluteValueOperation               */
/************************************************/

// Simply returns the negative of the variable.

@interface GLAbsoluteValueOperation : GLVariableOperation

- (GLAbsoluteValueOperation *) initWithFunction: (GLVariable *) variable;
// You can optionally set useComplexDivision to NO, and it will simply do an element-wise divide, ignoring complex math.
- (id) initWithOperand: (GLVariable *) operand shouldUseComplexArithmetic: (BOOL) useComplexArithmetic;

@property BOOL useComplexArithmetic;

@end

/************************************************/
/*		GLExponentialOperation					*/
/************************************************/
// variable = exponential( variable )
@interface GLExponentialOperation : GLVariableOperation
- (GLExponentialOperation *) initWithFunction: (GLVariable *) variable;
@end

/************************************************/
/*		GLSineOperation							*/
/************************************************/
// variable = sin( variable )
@interface GLSineOperation : GLVariableOperation
- (GLSineOperation *) initWithFunction: (GLVariable *) variable;
@end

/************************************************/
/*		GLCosineOperation						*/
/************************************************/
// variable = cos( variable )
@interface GLCosineOperation : GLVariableOperation
- (GLCosineOperation *) initWithFunction: (GLVariable *) variable;
@end

/************************************************/
/*		GLInverseTangentOperation				*/
/************************************************/
// variable = atan( variable )
@interface GLInverseTangentOperation : GLVariableOperation
- (GLInverseTangentOperation *) initWithFunction: (GLVariable *) variable;
@end

/************************************************/
/*		GLSquareRootOperation						*/
/************************************************/
// variable = sqrt( variable )
@interface GLSquareRootOperation : GLVariableOperation
- (GLSquareRootOperation *) initWithFunction: (GLVariable *) variable;
@end

/************************************************/
/*		GLFourierTransformOperation             */
/************************************************/

@interface GLFourierTransformOperation : GLVariableOperation
- (GLFourierTransformOperation *) initWithFunction: (GLVariable *) variable;
@end

/************************************************/
/*		GLSwapComplexOperation                  */
/************************************************/

@interface GLSwapComplexOperation : GLVariableOperation
- (GLSwapComplexOperation *) initWithFunction: (GLVariable *) variable;
@end

/************************************************/
/*		GLCopyVariableOperation                 */
/************************************************/

@interface GLCopyVariableOperation : GLVariableOperation
- (GLCopyVariableOperation *) initWithFunction: (GLVariable *) variable;
@end

/************************************************/
/*		GLMaxOperation							*/
/************************************************/

@interface GLMaxOperation : GLVariableOperation
- (GLMaxOperation *) initWithFunction: (GLVariable *) variable;
@end

/************************************************/
/*		GLAverageOperation						*/
/************************************************/

@interface GLAverageOperation : GLVariableOperation

- (GLAverageOperation *) initWithFunction: (GLVariable *) variable dimensionIndex: (NSUInteger) index;
@property NSUInteger dimIndex;

@end
