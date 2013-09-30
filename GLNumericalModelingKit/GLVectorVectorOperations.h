//
//  GLVectorVectorOperations.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 1/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import <GLNumericalModelingKit/GLVariableOperations.h>

/************************************************/
/*		GLAdditionOperation						*/
/************************************************/
// variable = leftVariable + rightVariable
@interface GLAdditionOperation : GLBinaryOperation
@end


/************************************************/
/*		GLSubtractionOperation					*/
/************************************************/
// variable = leftVariable - rightVariable
@interface GLSubtractionOperation : GLBinaryOperation
@end


/************************************************/
/*		GLMultiplicationOperation				*/
/************************************************/
// variable = leftVariable * rightVariable
@interface GLMultiplicationOperation : GLBinaryOperation
@end

/************************************************/
/*		GLAbsoluteLargestOperation				*/
/************************************************/
// variable = max( abs(leftVariable), abs(rightVariable ) element-wise
@interface GLAbsoluteLargestOperation : GLBinaryOperation
@end

/************************************************/
/*		GLDivisionOperation						*/
/************************************************/
// variable = leftVariable / rightVariable
@interface GLDivisionOperation : GLBinaryOperation

// You can optionally set useComplexDivision to NO, and it will simply do an element-wise divide, ignoring complex math.
- (id) initWithFirstOperand: (GLVariable *) fOperand secondOperand: (GLVariable *) sOperand shouldUseComplexArithmetic: (BOOL) useComplexArithmetic;
@property BOOL useComplexArithmetic;

@end

/************************************************/
/*		GLDotProductOperation					*/
/************************************************/

// variable = leftVariable dot rightVariable
@interface GLDotProductOperation : GLBinaryOperation
@end

/************************************************/
/*		GLSetVariableValueOperation				*/
/************************************************/

@interface GLSetVariableValueOperation : GLBinaryOperation
- (id) initWithVectorOperand: (GLVariable *) variable scalarVariableOperand: (GLVariable *) aScalarVariable indexString: (NSString *) indexString;
@property(copy) NSString * indexString;
@end