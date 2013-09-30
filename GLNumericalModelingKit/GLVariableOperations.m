//
//  GLVariableOperations.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 4/21/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import "GLVariableOperations.h"
#import "GLDimension.h"
#import "GLVariable.h"
#import "GLSpectralDifferentialOperator.h"
#import "GLDifferentialOperator.h"

/************************************************/
/*		GLVariableOperation						*/
/************************************************/

@implementation GLVariableOperation : NSOperation

- (id) init
{
    if ((self=[super init])) {
        self.graphvisDescription = @"unnamed";
    }
    return self;
}

@synthesize operationType;
@synthesize isEnqueued;

- (BOOL) canOperateInPlace {
	return NO;
}

- (BOOL) isEqualToOperation: (id) otherOperation {
    if ( ![[self class] isSubclassOfClass: [otherOperation class]] )  {
        return NO;
    }
    return YES;
}

@end

/************************************************/
/*		GLNullaryOperation						*/
/************************************************/

@interface GLNullaryOperation ()

// To be called after the result and operands are set.
- (void) setupDependencies;

// To be called after the operation is complete.
- (void) tearDownDependencies;

@end

@implementation GLNullaryOperation

- (id) initWithResult: (GLVariable *) resultVariable
{
	if (( self = [super init] ))
	{
		self.result = resultVariable;
		[self setupDependencies];
	}
	
    return self;
}

- (void) main
{
	self.blockOperation( self.result.data);
	
	[self tearDownDependencies];
}

@synthesize result;
@synthesize blockOperation;

- (void) setupDependencies
{
	if (self.result.lastOperation && ![self.dependencies containsObject:self.result.lastOperation]) {
		[self addDependency: self.result.lastOperation];
	}
	
	[self.result addOperation: self];
}

- (void) tearDownDependencies
{
	[self.result removeOperation: self];
	self.result = nil;
}

- (GLOperationType) operationType {
	return kGLNullaryOperation;
}

- (BOOL) isEqualToOperation: (id) otherOperation {
    if ( ![super isEqualToOperation: otherOperation] )  {
        return NO;
    }
    
    return YES;
}

@end

/************************************************/
/*		GLUnaryOperation						*/
/************************************************/

@interface GLUnaryOperation ()

// To be called after the result and operands are set.
- (void) setupDependencies;

// To be called after the operation is complete.
- (void) tearDownDependencies;

@end

@implementation GLUnaryOperation

- (id) initWithOperand: (GLVariable *) variable
{
	if ((self=[self initWithResult: nil operand: variable])) {
		
	}
	return self;
}

- (id) initWithResult: (GLVariable *) resultVariable operand: (GLVariable *) variable
{
	if (( self = [super init] ))
	{
		if (!resultVariable) {
			self.result = [[variable class] variableOfType: variable.dataFormat withDimensions: variable.dimensions forEquation: variable.equation];
		} else {
			self.result = resultVariable;
		}
		
		if ([[self.result class] isSubclassOfClass: [GLSpectralDifferentialOperator class]]) {
			[(GLDifferentialOperator *) self.result setToDimensions: [(GLDifferentialOperator *) variable toDimensions]];
			[(GLDifferentialOperator *) self.result setFromDimensions: [(GLDifferentialOperator *) variable fromDimensions]];
		}
		
		self.operand = variable;
		
		[self setupDependencies];
	}
	
    return self;
}

- (void) main
{	
	self.blockOperation( self.result.data, self.operand.data );
	
	[self tearDownDependencies];
}

@synthesize result;
@synthesize operand;
@synthesize blockOperation;

- (void) setupDependencies
{
	if (self.operand.lastOperation && ![self.dependencies containsObject:self.operand.lastOperation]) {
		[self addDependency: self.operand.lastOperation];
	}
	if (self.result.lastOperation && ![self.dependencies containsObject:self.result.lastOperation]) {
		[self addDependency: self.result.lastOperation];
	}
	
	[self.result addOperation: self];
}

- (void) tearDownDependencies
{	
	[self.result removeOperation: self];
	self.result = nil;
	self.operand = nil;
}

- (GLOperationType) operationType {
	return kGLUnaryOperation;
}

- (BOOL) isEqualToOperation: (id) otherOperation {
    if ( ![super isEqualToOperation: otherOperation] )  {
        return NO;
    }
    
    GLUnaryOperation * op = otherOperation;
    if ( self.operand !=  op.operand) {
        return NO;
    }
    return YES;
}

@end

/************************************************/
/*		GLBinaryOperation						*/
/************************************************/

@implementation GLBinaryOperation

- (id) initWithFirstOperand: (GLVariable *) fOperand secondOperand: (GLVariable *) sOperand
{
	if ( (self = [self initWithResult: nil firstOperand: fOperand secondOperand: sOperand])) {
		
	}
	return self;
}

- (id) initWithResult: (GLVariable *) resultVariable firstOperand: (GLVariable *) fOperand secondOperand: (GLVariable *) sOperand;
{
	// Generally, the dimensions of one variable, must contain the dimensions of another variable
	
	BOOL firstOpIsScalar = !fOperand.dimensions.count;
	BOOL secondOpIsScalar = !sOperand.dimensions.count;
	BOOL scalarOp = firstOpIsScalar || secondOpIsScalar;

	if ( !scalarOp && ![fOperand.dimensions isEqualToArray: sOperand.dimensions] ) {
        [NSException raise: @"DimensionsNotEqualException" format: @"Cannot binary operate two variables of different dimensions, unless one is a scalar"];
    }

	#warning this is a total mess.
	if ([[fOperand class] isSubclassOfClass: [GLSpectralDifferentialOperator class]] && [[sOperand class] isSubclassOfClass: [GLSpectralDifferentialOperator class]] && [[resultVariable class] isSubclassOfClass: [GLSpectralDifferentialOperator class]]) {
		NSArray *fToDims = [(GLSpectralDifferentialOperator *) fOperand toDimensions];
		NSArray *sToDims = [(GLSpectralDifferentialOperator *) sOperand toDimensions];
		
		if ( fToDims && sToDims && ![fToDims isEqualToArray: sToDims] ) {
			[NSException raise: @"DimensionsNotEqualException" format: @"Differential operators must have the same destination dimension."];
		}
	}
	

#warning this is actuallly still true in some cases.
//	if ( fOperand.isComplex != sOperand.isComplex ) {
//		[NSException raise: @"MethodNotYetImplemented" format: @"Cannot perform a binary operation on two variables that are not either both complex or both real."];
//	}
	
	
    
	if (( self = [super init] )) {
		self.firstOperand = fOperand;
		self.secondOperand = sOperand;
		
		if (!resultVariable) {
			NSArray *dimensions;
			Class aClass;
			if (firstOpIsScalar) {
				dimensions = [NSArray arrayWithArray: sOperand.dimensions];
				aClass = [sOperand class];
			} else if (secondOpIsScalar) {
				dimensions = [NSArray arrayWithArray: fOperand.dimensions];
				aClass = [fOperand class];
			} else { // Prevent extra differential operators from being stamped out.
				if ( [[fOperand class] isSubclassOfClass: [sOperand class]] && ![[sOperand class] isSubclassOfClass: [fOperand class]] ) {
					dimensions = [NSArray arrayWithArray: sOperand.dimensions];
					aClass = [sOperand class];
				} else if ( ![[fOperand class] isSubclassOfClass: [sOperand class]] && [[sOperand class] isSubclassOfClass: [fOperand class]] ) {
					dimensions = [NSArray arrayWithArray: fOperand.dimensions];
					aClass = [fOperand class];
				} else if ( [[fOperand class] isSubclassOfClass: [sOperand class]] && [[sOperand class] isSubclassOfClass: [fOperand class]] ) {
					dimensions = [NSArray arrayWithArray: fOperand.dimensions];
					aClass = [fOperand class];
				}
			}
			
			BOOL isComplex = fOperand.isComplex || sOperand.isComplex;
			GLDataFormat format = isComplex ? kGLSplitComplexDataFormat : kGLRealDataFormat;
			
			self.result = [aClass variableOfType: format withDimensions: dimensions forEquation: self.firstOperand.equation];
		} else {
			self.result = resultVariable;
		}
		
		if ([[self.result class] isSubclassOfClass: [GLSpectralDifferentialOperator class]]) {
			[(GLDifferentialOperator *) self.result setToDimensions: [(GLDifferentialOperator *) fOperand toDimensions]];
			[(GLDifferentialOperator *) self.result setFromDimensions: [(GLDifferentialOperator *) fOperand fromDimensions]];
		}
		
		[self setupDependencies];
    }
    return self;
}

- (void) main
{	
	self.blockOperation( self.result.data, self.firstOperand.data, self.secondOperand.data );
	
	[self tearDownDependencies];
}

@synthesize result=_result;
@synthesize firstOperand=_firstOperand;
@synthesize secondOperand=_secondOperand;
@synthesize blockOperation;

- (void) setupDependencies
{
	if (self.firstOperand.lastOperation && ![self.dependencies containsObject:self.firstOperand.lastOperation]) {
		[self addDependency: self.firstOperand.lastOperation];
	}
	if (self.secondOperand.lastOperation && ![self.dependencies containsObject:self.secondOperand.lastOperation]) {
		[self addDependency: self.secondOperand.lastOperation];
	}
	if (self.result.lastOperation && ![self.dependencies containsObject:self.result.lastOperation]) {
		[self addDependency: self.result.lastOperation];
	}
	
	[self.result addOperation: self];
}

- (void) tearDownDependencies
{	
	[self.result removeOperation: self];
	self.result = nil;
	self.firstOperand = nil;
	self.secondOperand = nil;
}

- (GLOperationType) operationType {
	return kGLBinaryOperation;
}

- (BOOL) isEqualToOperation: (id) otherOperation {
    if ( ![super isEqualToOperation: otherOperation] )  {
        return NO;
    }
    
    GLBinaryOperation * op = otherOperation;
    if ( self.firstOperand !=  op.firstOperand) {
        return NO;
    }
    if ( self.secondOperand !=  op.secondOperand) {
        return NO;
    }
    return YES;
}

@end


/************************************************/
/*		GLUnaryVectorOperation						*/
/************************************************/

@interface GLUnaryVectorOperation ()

// To be called after the result and operands are set.
- (void) setupDependencies;

// To be called after the operation is complete.
- (void) tearDownDependencies;

@end

@implementation GLUnaryVectorOperation

- (id) initWithOperand: (NSArray *) variable
{
	return [self initWithResult: nil operand: variable];
}

- (id) initWithResult: (NSArray *) resultVariable operand: (NSArray *) operandVariable
{
	if (( self = [super init] ))
	{
		if (!resultVariable) {
			NSMutableArray *array = [[NSMutableArray alloc] initWithCapacity: operandVariable.count];
			for (GLVariable *variable in operandVariable) {
				[array addObject: [[variable class] variableOfType: variable.dataFormat withDimensions: variable.dimensions forEquation: variable.equation]];
			}
			self.result = array;
		} else {
			self.result = resultVariable;
		}
		
		self.operand = operandVariable;
		
		[self setupDependencies];
	}
	
    return self;
}

- (void) main
{
	NSMutableArray *resultBuffer = [[NSMutableArray alloc] initWithCapacity: self.result.count];
	NSMutableArray *operandBuffer = [[NSMutableArray alloc] initWithCapacity: self.operand.count];
	for (GLVariable *variable in self.result) {
		[resultBuffer addObject: variable.data];
	}
	for (GLVariable *variable in self.operand) {
		[operandBuffer addObject: variable.data];
	}
	self.blockOperation( resultBuffer, operandBuffer );
	
	[self tearDownDependencies];
}

@synthesize result;
@synthesize operand;
@synthesize blockOperation;

- (void) setupDependencies
{
	for (GLVariable *variable in self.operand) {
		if (variable.lastOperation && ![self.dependencies containsObject:variable.lastOperation]) {
			[self addDependency: variable.lastOperation];
		}
	}
	for (GLVariable *variable in self.result) {
		if (variable.lastOperation && ![self.dependencies containsObject:variable.lastOperation]) {
			[self addDependency: variable.lastOperation];
		}
	}
	for (GLVariable *variable in self.result) {
		[variable addOperation: self];
	}	
}

- (void) tearDownDependencies
{	
	for (GLVariable *variable in self.result) {
		[variable removeOperation: self];
	}
	self.result = nil;
	self.operand = nil;
}

- (GLOperationType) operationType {
	return kGLUnaryVectorOperation;
}

- (BOOL) isEqualToOperation: (id) otherOperation {
    if ( ![super isEqualToOperation: otherOperation] )  {
        return NO;
    }
    
    GLUnaryVectorOperation * op = otherOperation;
    if (self.operand.count != op.operand.count) {
        return NO;
    }
    for ( NSUInteger i=0; i<self.operand.count; i++) {
        if ( self.operand[i] != op.operand[i]) {
            return NO;
        }
    }

    return YES;
}

@end


/************************************************/
/*		GLBinaryVectorOperation						*/
/************************************************/

@interface GLBinaryVectorOperation ()

// To be called after the result and operands are set.
- (void) setupDependencies;

// To be called after the operation is complete.
- (void) tearDownDependencies;

@end

@implementation GLBinaryVectorOperation

- (id) initWithFirstOperand: (NSArray *) fOperand secondOperand: (NSArray *) sOperand
{
	return [self initWithResult: nil firstOperand: fOperand secondOperand: sOperand];
}

- (id) initWithResult: (NSArray *) resultVariable firstOperand: (NSArray *) fOperand secondOperand: (NSArray *) sOperand
{
	if (( self = [super init] ))
	{
		if (!resultVariable) {
			NSMutableArray *array = [[NSMutableArray alloc] initWithCapacity: fOperand.count];
			for (GLVariable *variable in fOperand) {
				[array addObject: [[variable class] variableOfType: variable.dataFormat withDimensions: variable.dimensions forEquation: variable.equation]];
			}
			self.result = array;
		} else {
			self.result = resultVariable;
		}
		
		self.firstOperand = fOperand;
		self.secondOperand = sOperand;
		
		[self setupDependencies];
	}
	
    return self;
}

- (void) main
{
	NSMutableArray *resultBuffer = [[NSMutableArray alloc] initWithCapacity: self.result.count];
	NSMutableArray *firstOperandBuffer = [[NSMutableArray alloc] initWithCapacity: self.firstOperand.count];
	NSMutableArray *secondOperandBuffer = [[NSMutableArray alloc] initWithCapacity: self.secondOperand.count];
	for (GLVariable *variable in self.result) {
		[resultBuffer addObject: variable.data];
	}
	for (GLVariable *variable in self.firstOperand) {
		[firstOperandBuffer addObject: variable.data];
	}
	for (GLVariable *variable in self.secondOperand) {
		[secondOperandBuffer addObject: variable.data];
	}
	self.blockOperation( resultBuffer, firstOperandBuffer, secondOperandBuffer );
	
	[self tearDownDependencies];
}

@synthesize result;
@synthesize firstOperand;
@synthesize secondOperand;
@synthesize blockOperation;

- (void) setupDependencies
{
	for (GLVariable *variable in self.firstOperand) {
		if (variable.lastOperation && ![self.dependencies containsObject:variable.lastOperation]) {
			[self addDependency: variable.lastOperation];
		}
	}
	for (GLVariable *variable in self.secondOperand) {
		if (variable.lastOperation && ![self.dependencies containsObject:variable.lastOperation]) {
			[self addDependency: variable.lastOperation];
		}
	}
	for (GLVariable *variable in self.result) {
		if (variable.lastOperation && ![self.dependencies containsObject:variable.lastOperation]) {
			[self addDependency: variable.lastOperation];
		}
	}
	for (GLVariable *variable in self.result) {
		[variable addOperation: self];
	}	
}

- (void) tearDownDependencies
{	
	for (GLVariable *variable in self.result) {
		[variable removeOperation: self];
	}
	self.result = nil;
	self.firstOperand = nil;
	self.secondOperand = nil;
}

- (GLOperationType) operationType {
	return kGLBinaryVectorOperation;
}

- (BOOL) isEqualToOperation: (id) otherOperation {
    if ( ![super isEqualToOperation: otherOperation] )  {
        return NO;
    }
    
    GLBinaryVectorOperation * op = otherOperation;
    if (self.firstOperand.count != op.firstOperand.count) {
        return NO;
    }
    if (self.secondOperand.count != op.secondOperand.count) {
        return NO;
    }
    for ( NSUInteger i=0; i<self.firstOperand.count; i++) {
        if ( self.firstOperand[i] != op.firstOperand[i]) {
            return NO;
        }
    }
    for ( NSUInteger i=0; i<self.secondOperand.count; i++) {
        if ( self.secondOperand[i] != op.secondOperand[i]) {
            return NO;
        }
    }
    
    return YES;
}

@end












