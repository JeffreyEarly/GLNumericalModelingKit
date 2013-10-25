//
//  GLVectorScalarOperations.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 1/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import <GLNumericalModelingKit/GLVariableOperations.h>

/************************************************/
/*		GLVectorScalarOperation					*/
/************************************************/

// This subclass of GLVariableOperation takes a single variable operand and has one, constant, scalar operand.
// Two vector-scalar operations are equal only if their scalar operands are equal.

@interface GLVectorScalarOperation : GLVariableOperation

// Init with only the operand to automatically create the target variable.
- (id) initWithVectorOperand: (GLTensor *) vOperand scalarOperand: (GLFloat) sOperand;

@property(readwrite, assign, nonatomic) GLFloat scalarOperand;

@end


/************************************************/
/*		GLScalarMultiplyOperation				*/
/************************************************/

@interface GLScalarMultiplyOperation : GLVectorScalarOperation

@end

/************************************************/
/*		GLScalarAddOperation					*/
/************************************************/

@interface GLScalarAddOperation : GLVectorScalarOperation

@end

/************************************************/
/*		GLScalarDivideOperation					*/
/************************************************/
// variable = scalar / variable
@interface GLScalarDivideOperation : GLVectorScalarOperation

@end

/************************************************/
/*		GLPowerOperation						*/
/************************************************/

@interface GLPowerOperation : GLVectorScalarOperation

@end

/************************************************/
/*		GLScalarThresholdOperation				*/
/************************************************/
// variable = max( operand, scalar )
@interface GLScalarThresholdOperation : GLVectorScalarOperation

@end

/************************************************/
/*		GLZeroThresholdOperation				*/
/************************************************/
// variable = operand > scalar ? operand : 0.0
@interface GLZeroThresholdOperation : GLVectorScalarOperation

@end

/************************************************/
/*		GLClipOperation				*/
/************************************************/
// variable = clip(operand, min, max)
@interface GLClipOperation : GLVariableOperation

- (id) initWithVectorOperand: (GLVariable *) vOperand firstScalarOperand: (GLFloat) fsOperand secondScalarOperand: (GLFloat) ssOperand;

@property GLFloat firstScalarOperand;
@property GLFloat secondScalarOperand;

@end

/************************************************/
/*		GLSetValueOperation						*/
/************************************************/

@interface GLSetValueOperation : GLVectorScalarOperation

- (id) initWithVectorOperand: (GLVariable *) vOperand scalarOperand: (GLFloat) sOperand indexString: (NSString *) indexString;
@property(copy) NSString * indexString;

@end

/************************************************/
/*		GLScaleOperation				*/
/************************************************/

@interface GLScaleOperation : GLVectorScalarOperation
-(id) initWithOperand: (GLVariable *) vOperand variableScale: (GLFloat) sOperand units: (NSString *) varUnits dimensionScale: (GLFloat) dimScale translation: (GLFloat) delta units: (NSString *) dimUnits;
@property(readwrite, assign, nonatomic) GLFloat dimScale;
@property(readwrite, assign, nonatomic) GLFloat dimTranslation;
@end