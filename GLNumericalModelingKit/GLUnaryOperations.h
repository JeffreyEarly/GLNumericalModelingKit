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

@end

/************************************************/
/*		GLAbsoluteValueOperation               */
/************************************************/

// Simply returns the negative of the variable.

@interface GLAbsoluteValueOperation : GLVariableOperation

// You can optionally set useComplexDivision to NO, and it will simply do an element-wise divide, ignoring complex math.
- (id) initWithOperand: (GLVariable *) operand shouldUseComplexArithmetic: (BOOL) useComplexArithmetic;

@property BOOL useComplexArithmetic;

@end

/************************************************/
/*		GLExponentialOperation					*/
/************************************************/
// variable = exponential( variable )
@interface GLExponentialOperation : GLVariableOperation

@end

/************************************************/
/*		GLSineOperation							*/
/************************************************/
// variable = sin( variable )
@interface GLSineOperation : GLVariableOperation

@end

/************************************************/
/*		GLCosineOperation						*/
/************************************************/
// variable = cos( variable )
@interface GLCosineOperation : GLVariableOperation

@end

/************************************************/
/*		GLInverseTangentOperation				*/
/************************************************/
// variable = atan( variable )
@interface GLInverseTangentOperation : GLVariableOperation

@end

/************************************************/
/*		GLSquareRootOperation						*/
/************************************************/
// variable = sqrt( variable )
@interface GLSquareRootOperation : GLVariableOperation

@end

/************************************************/
/*		GLFourierTransformOperation             */
/************************************************/

@interface GLFourierTransformOperation : GLVariableOperation

@end

/************************************************/
/*		GLSwapComplexOperation                  */
/************************************************/

@interface GLSwapComplexOperation : GLVariableOperation

@end

/************************************************/
/*		GLCopyVariableOperation                 */
/************************************************/

@interface GLCopyVariableOperation : GLVariableOperation

@end

/************************************************/
/*		GLMaxOperation							*/
/************************************************/

@interface GLMaxOperation : GLVariableOperation

@end

/************************************************/
/*		GLAverageOperation						*/
/************************************************/

@interface GLAverageOperation : GLVariableOperation

- (id) initWithOperand: (GLVariable *) variable dimensionIndex: (NSUInteger) index;
@property NSUInteger dimIndex;

@end
