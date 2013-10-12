//
//  GLBasisTransformationOperations.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 6/14/12.
//
//

#import <GLNumericalModelingKit/GLVariableOperations.h>

/************************************************/
/*		GLBasisTransformOperation               */
/************************************************/

// Return nil if no transformation is necessary. Raises an exception if it can't do it.
@interface GLBasisTransformOperation : GLVariableOperation
+ (id) basisTransformationWithOperand: (GLVariable *) variable destinationBasis: (NSArray *) toBasis;
+ (void) readWisdom;
+ (void) saveWisdom;

// If set to yes, this will apply a two-thirds filter before an inverse-transformation.
+ (void) setShouldAntialias: (BOOL) flag;
+ (BOOL) shouldAntialias;
+ (GLVariable *) antialiasingFilterFromDimensions: (NSArray *) dimensions forEquation: (GLEquation *) equation;

@property(readwrite, strong, nonatomic) NSArray *fromBasis;
@property(readwrite, strong, nonatomic) NSArray *toBasis;

@end

/************************************************/
/*		GLRealToRealTransformOperation          */
/************************************************/

// Transform real variables between delta function basis, and sine/cosine basis.
@interface GLRealToRealTransformOperation : GLBasisTransformOperation
- (id) initWithOperand: (GLVariable *) variable destinationBasis: (NSArray *) toBasis;
@end

/************************************************/
/*		GLRealToComplexTransformOperation       */
/************************************************/

// Transform real variables between delta function basis, and exponential basis.
// Also does the inverse for Hermitian variables.
@interface GLRealToComplexTransformOperation : GLBasisTransformOperation
- (id) initWithOperand: (GLVariable *) variable destinationBasis: (NSArray *) toBasis;
@end

/************************************************/
/*		GLComplexToComplexTransformOperation    */
/************************************************/

// Transform complex variables between delta function basis, and exponential basis.
@interface GLComplexToComplexTransformOperation : GLBasisTransformOperation
- (id) initWithOperand: (GLVariable *) variable destinationBasis: (NSArray *) toBasis;
@end


/************************************************/
/*		GLMatrixFFTTransformOperation           */
/************************************************/

// Return nil if no transformation is necessary. Raises an exception if it can't do it.
@interface GLMatrixFFTTransformOperation : GLVariableOperation
+ (id) basisTransformationWithOperand: (GLVariable *) variable destinationBasis: (NSArray *) toBasis;
@end