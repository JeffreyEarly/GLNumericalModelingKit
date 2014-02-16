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
- (GLNegationOperation *) initWithVariable: (GLVariable *) variable;
@end

/************************************************/
/*		GLAbsoluteValueOperation               */
/************************************************/

// Returns the negative of the variable, which will always be real.

@interface GLAbsoluteValueOperation : GLVariableOperation
- (GLAbsoluteValueOperation *) initWithVariable: (GLVariable *) variable;
@end

/************************************************/
/*		GLExponentialOperation					*/
/************************************************/
// variable = exponential( variable )
@interface GLExponentialOperation : GLVariableOperation
- (GLExponentialOperation *) initWithVariable: (GLVariable *) variable;
@end

/************************************************/
/*		GLLogarithmOperation					*/
/************************************************/
// variable = log( variable )
@interface GLLogarithmOperation : GLVariableOperation
- (GLLogarithmOperation *) initWithVariable: (GLVariable *) variable;
@end

/************************************************/
/*		GLSineOperation							*/
/************************************************/
// variable = sin( variable )
@interface GLSineOperation : GLVariableOperation
- (GLSineOperation *) initWithVariable: (GLVariable *) variable;
@end

/************************************************/
/*		GLCosineOperation						*/
/************************************************/
// variable = cos( variable )
@interface GLCosineOperation : GLVariableOperation
- (GLCosineOperation *) initWithVariable: (GLVariable *) variable;
@end

/************************************************/
/*		GLHyperbolicTangentOperation			*/
/************************************************/
// variable = sin( variable )
@interface GLHyperbolicTangentOperation : GLVariableOperation
- (GLHyperbolicTangentOperation *) initWithVariable: (GLVariable *) variable;
@end

/************************************************/
/*		GLInverseTangentOperation				*/
/************************************************/
// variable = atan( variable )
@interface GLInverseTangentOperation : GLVariableOperation
- (GLInverseTangentOperation *) initWithVariable: (GLVariable *) variable;
@end

/************************************************/
/*		GLSquareRootOperation						*/
/************************************************/
// variable = sqrt( variable )
@interface GLSquareRootOperation : GLVariableOperation
- (GLSquareRootOperation *) initWithVariable: (GLVariable *) variable;
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
/*		GLMinOperation							*/
/************************************************/

@interface GLMinOperation : GLVariableOperation
- (GLMinOperation *) initWithFunction: (GLFunction *) variable;
@end

/************************************************/
/*		GLAverageOperation						*/
/************************************************/

@interface GLAverageOperation : GLVariableOperation

- (GLAverageOperation *) initWithFunction: (GLFunction *) variable dimensionIndex: (NSUInteger) index;
- (GLAverageOperation *) initWithFunction: (GLFunction *) variable;
@property NSUInteger dimIndex;
@end

/************************************************/
/*		GLSummationOperation					*/
/************************************************/

@interface GLSummationOperation : GLVariableOperation

- (GLSummationOperation *) initWithFunction: (GLFunction *) variable dimensionIndex: (NSUInteger) index;
- (GLSummationOperation *) initWithFunction: (GLFunction *) variable;
@property NSUInteger dimIndex;
@end

/************************************************/
/*		GLIntegrationOperation					*/
/************************************************/

@interface GLIntegrationOperation : GLVariableOperation
- (GLIntegrationOperation *) initWithFunction: (GLFunction *) variable;
@end

/************************************************/
/*		GLInterleavedToSplitComplexOperation    */
/************************************************/

@interface GLInterleavedToSplitComplexOperation : GLVariableOperation
- (GLInterleavedToSplitComplexOperation *) initWithVariable: (GLVariable *) variable;
@end

/************************************************/
/*		GLSplitToInterleavedComplexOperation    */
/************************************************/

@interface GLSplitToInterleavedComplexOperation : GLVariableOperation
- (GLSplitToInterleavedComplexOperation *) initWithVariable: (GLVariable *) variable;
@end

/************************************************/
/*		GLDataTransposeOperation				*/
/************************************************/

@interface GLDataTransposeOperation : GLVariableOperation
- (GLDataTransposeOperation *) initWithLinearTransform: (GLLinearTransform *) transform;
@end
/************************************************/
/*		GLDenseMatrixOperation					*/
/************************************************/

// Takes the existing format, and converts to a dense matrix.
@interface GLDenseMatrixOperation : GLVariableOperation
- (GLDenseMatrixOperation *) initWithLinearTransform: (GLLinearTransform *) transform;
@end

/************************************************/
/*		GLMakeHermitianOperation                     */
/************************************************/

// Fixes the components to make hermitian.

@interface GLMakeHermitianOperation : GLVariableOperation
- (GLMakeHermitianOperation *) initWithFunction: (GLVariable *) variable;
@end