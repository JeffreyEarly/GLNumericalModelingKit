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
- (GLNegationOperation *) initWithFunction: (GLFunction *) variable;
@end

/************************************************/
/*		GLAbsoluteValueOperation               */
/************************************************/

// Simply returns the negative of the variable.

@interface GLAbsoluteValueOperation : GLVariableOperation

- (GLAbsoluteValueOperation *) initWithFunction: (GLFunction *) variable;
// You can optionally set useComplexDivision to NO, and it will simply do an element-wise divide, ignoring complex math.
- (id) initWithOperand: (GLFunction *) operand shouldUseComplexArithmetic: (BOOL) useComplexArithmetic;

@property BOOL useComplexArithmetic;

@end

/************************************************/
/*		GLExponentialOperation					*/
/************************************************/
// variable = exponential( variable )
@interface GLExponentialOperation : GLVariableOperation
- (GLExponentialOperation *) initWithFunction: (GLFunction *) variable;
@end

/************************************************/
/*		GLSineOperation							*/
/************************************************/
// variable = sin( variable )
@interface GLSineOperation : GLVariableOperation
- (GLSineOperation *) initWithFunction: (GLFunction *) variable;
@end

/************************************************/
/*		GLCosineOperation						*/
/************************************************/
// variable = cos( variable )
@interface GLCosineOperation : GLVariableOperation
- (GLCosineOperation *) initWithFunction: (GLFunction *) variable;
@end

/************************************************/
/*		GLInverseTangentOperation				*/
/************************************************/
// variable = atan( variable )
@interface GLInverseTangentOperation : GLVariableOperation
- (GLInverseTangentOperation *) initWithFunction: (GLFunction *) variable;
@end

/************************************************/
/*		GLSquareRootOperation						*/
/************************************************/
// variable = sqrt( variable )
@interface GLSquareRootOperation : GLVariableOperation
- (GLSquareRootOperation *) initWithFunction: (GLFunction *) variable;
@end

/************************************************/
/*		GLFourierTransformOperation             */
/************************************************/

@interface GLFourierTransformOperation : GLVariableOperation
- (GLFourierTransformOperation *) initWithFunction: (GLFunction *) variable;
@end

/************************************************/
/*		GLSwapComplexOperation                  */
/************************************************/

@interface GLSwapComplexOperation : GLVariableOperation
- (GLSwapComplexOperation *) initWithFunction: (GLFunction *) variable;
@end

/************************************************/
/*		GLCopyVariableOperation                 */
/************************************************/

@interface GLCopyVariableOperation : GLVariableOperation
- (GLCopyVariableOperation *) initWithFunction: (GLFunction *) variable;
@end

/************************************************/
/*		GLMaxOperation							*/
/************************************************/

@interface GLMaxOperation : GLVariableOperation
- (GLMaxOperation *) initWithFunction: (GLFunction *) variable;
@end

/************************************************/
/*		GLAverageOperation						*/
/************************************************/

@interface GLAverageOperation : GLVariableOperation

- (GLAverageOperation *) initWithFunction: (GLFunction *) variable dimensionIndex: (NSUInteger) index;
@property NSUInteger dimIndex;

@end
