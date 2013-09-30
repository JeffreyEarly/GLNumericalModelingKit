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

@interface GLNegationOperation : GLUnaryOperation

@end

/************************************************/
/*		GLAbsoluteValueOperation               */
/************************************************/

// Simply returns the negative of the variable.

@interface GLAbsoluteValueOperation : GLUnaryOperation

// You can optionally set useComplexDivision to NO, and it will simply do an element-wise divide, ignoring complex math.
- (id) initWithOperand: (GLVariable *) operand shouldUseComplexArithmetic: (BOOL) useComplexArithmetic;

@property BOOL useComplexArithmetic;

@end

/************************************************/
/*		GLExponentialOperation					*/
/************************************************/
// variable = exponential( variable )
@interface GLExponentialOperation : GLUnaryOperation

@end

/************************************************/
/*		GLSineOperation							*/
/************************************************/
// variable = sin( variable )
@interface GLSineOperation : GLUnaryOperation

@end

/************************************************/
/*		GLCosineOperation						*/
/************************************************/
// variable = cos( variable )
@interface GLCosineOperation : GLUnaryOperation

@end

/************************************************/
/*		GLInverseTangentOperation				*/
/************************************************/
// variable = atan( variable )
@interface GLInverseTangentOperation : GLUnaryOperation

@end

/************************************************/
/*		GLSquareRootOperation						*/
/************************************************/
// variable = sqrt( variable )
@interface GLSquareRootOperation : GLUnaryOperation

@end

/************************************************/
/*		GLFourierTransformOperation             */
/************************************************/

@interface GLFourierTransformOperation : GLUnaryOperation

@end

/************************************************/
/*		GLSwapComplexOperation                  */
/************************************************/

@interface GLSwapComplexOperation : GLUnaryOperation

@end

/************************************************/
/*		GLCopyVariableOperation                 */
/************************************************/

@interface GLCopyVariableOperation : GLUnaryOperation

@end

/************************************************/
/*		GLMaxOperation							*/
/************************************************/

@interface GLMaxOperation : GLUnaryOperation

@end

/************************************************/
/*		GLAverageOperation						*/
/************************************************/

@interface GLAverageOperation : GLUnaryOperation

- (id) initWithOperand: (GLVariable *) variable dimensionIndex: (NSUInteger) index;
@property NSUInteger dimIndex;

@end
