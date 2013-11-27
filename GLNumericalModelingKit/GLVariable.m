//
//  GLTensor.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 10/11/13.
//  Copyright (c) 2013 Jeffrey J. Early. All rights reserved.
//

#import "GLVariable.h"

#import "GLMemoryPool.h"
#import "GLNetCDFFile.h"
#import "GLEquation.h"
#import "GLLinearTransform.h"
#import "GLVectorVectorOperations.h"
#import "GLVectorScalarOperations.h"
#include <mach/mach_time.h>

@interface GLVariable ()
@property(readwrite, strong, nonatomic) NSMutableDictionary *metadata;
@property(readwrite, assign, nonatomic) NSUInteger uniqueID;
@property(readwrite, strong, nonatomic) NSMutableData *data;
@property(readwrite, strong, nonatomic) NSMutableArray *pendingOperations;
@property(readwrite, strong) NSMutableArray *existingOperations;
@property(readwrite, assign, nonatomic) GLDataFormat dataFormat;
@property(readwrite, assign, nonatomic) BOOL isComplex;
@end

GLSplitComplex splitComplexFromData( NSData *data )
{
	GLSplitComplex fbar;
	fbar.realp = (void *) data.bytes;
	fbar.imagp = (void *) (data.bytes + data.length/2);
	return fbar;
}

@implementation GLVariable

+ (GLVariable *) variableWithPrototype: (GLVariable *) variable
{
    if (variable.rank == 0) {
        GLScalar *scalar = (GLScalar *) variable;
        return [[GLScalar alloc] initWithType: scalar.dataFormat forEquation:scalar.equation];
    } else if (variable.rank == 1) {
        GLFunction *function = (GLFunction *) variable;
        return [[function class] functionOfType: function.dataFormat withDimensions: function.dimensions forEquation: function.equation];
    }  else if (variable.rank == 2) {
        GLLinearTransform *matrix = (GLLinearTransform *) variable;
        return [GLLinearTransform transformOfType: matrix.dataFormat withFromDimensions: matrix.fromDimensions toDimensions: matrix.toDimensions inFormat: matrix.matrixFormats forEquation:matrix.equation matrix:nil];
    }
    return nil;
}

- (id) initWithType: (GLDataFormat) dataFormat withEquation: (GLEquation *) theEquation
{
	if (!theEquation) {
		NSLog(@"Attempted to initialize GLTensor without an equation!!!");
		return nil;
	}
	
	if ((self = [super init])) {
		_pendingOperations = [[NSMutableArray alloc] init];
		_existingOperations = [[NSMutableArray alloc] init];
		_equation = theEquation;
		_uniqueID = mach_absolute_time();
		_dataFormat = dataFormat;
		_isComplex = dataFormat != kGLRealDataFormat;
		_isImaginaryPartZero = dataFormat == kGLRealDataFormat;
	}
	
	return self;
}


- (NSMutableDictionary *) metadata {
	if (!_metadata) {
		_metadata = [NSMutableDictionary dictionary];
	}
	return _metadata;
}

- (void) setUnits:(NSString *)theUnits
{
	_units = theUnits;
	[self.metadata setValue: theUnits forKey: @"units"];
}

@dynamic isPurelyReal;
@dynamic isPurelyImaginary;

- (BOOL) isPurelyReal {
    return self.isImaginaryPartZero;
}

- (void) setIsPurelyReal:(BOOL)isPurelyReal {
	self.isImaginaryPartZero = isPurelyReal;
}

- (BOOL) isPurelyImaginary {
    return self.isRealPartZero;
}

- (void) setIsPurelyImaginary:(BOOL)isPurelyImaginary {
	self.isRealPartZero = isPurelyImaginary;
}

/************************************************/
/*		Data									*/
/************************************************/

#pragma mark -
#pragma mark Data
#pragma mark

// for some reason the compiler refuses to create this ivar, with this name, automatically.
@synthesize data=_data;

- (NSMutableData *) data
{
    if (!_data && self.dataBytes) {
		_data =[[GLMemoryPool sharedMemoryPool] dataWithLength: self.dataBytes];
	}
    return _data;
}

- (void) setData: (NSData *) newData
{
	if (!_data) {
		_data =[[GLMemoryPool sharedMemoryPool] dataWithLength: self.dataBytes];
	}
	
	memcpy( _data.mutableBytes, newData.bytes, newData.length);
}

- (BOOL) hasData {
	return _data != nil ? YES : NO;
}

- (GLFloat *) pointerValue
{
	[self solve];
    GLFloat *f = self.data.mutableBytes;
    return f;
}

- (GLSplitComplex) splitComplex
{
	[self solve];
	GLSplitComplex fbar;
	if (self.isComplex) {
		fbar.realp = self.data.mutableBytes;
		fbar.imagp = self.data.mutableBytes + self.data.length/2;
	} else {
		fbar.realp = NULL;
		fbar.imagp = NULL;
		NSLog(@"Error! Requesting splitComplex from a real variable!!!");
	}
	return fbar;
}

- (void) solve
{
	[self.equation solveForVariable: self];
}

- (void) zero
{
	vGL_vclr( self.pointerValue, 1, self.nDataElements);
}

/************************************************/
/*		Superclass Overrides					*/
/************************************************/

#pragma mark -
#pragma mark Superclass Overrides
#pragma mark

- (BOOL) isEqual: (id) otherObject
{
	return ([[self class] isSubclassOfClass: [otherObject class]] && self.uniqueID == [(GLVariable *)otherObject uniqueID]);
}

- (NSUInteger)hash {
    return _uniqueID;
}

/************************************************/
/*		Operations								*/
/************************************************/

#pragma mark -
#pragma mark Operations
#pragma mark

- (GLVariable *) plus: (id) otherVariable
{
	GLVariableOperation *operation;
	if ([[otherVariable class] isSubclassOfClass: [NSNumber class]]) {
		operation  = [[GLScalarAddOperation alloc] initWithVectorOperand: self scalarOperand: [(NSNumber *) otherVariable doubleValue]];
	} else {
		operation = [[GLAdditionOperation alloc] initWithFirstOperand: self secondOperand: otherVariable];
	}
    operation = [self replaceWithExistingOperation: operation];
	return operation.result[0];
}

- (GLVariable *) minus: (id) otherVariable
{
	GLVariableOperation *operation;
	if ([[otherVariable class] isSubclassOfClass: [NSNumber class]]) {
		operation  = [[GLScalarAddOperation alloc] initWithVectorOperand: self scalarOperand: -[(NSNumber *) otherVariable doubleValue]];
	} else {
		operation = [[GLSubtractionOperation alloc] initWithFirstOperand: self secondOperand: otherVariable];
	}
    operation = [self replaceWithExistingOperation: operation];
	return operation.result[0];
}

- (GLVariable *) multiply: (id) otherVariable
{
	GLVariableOperation *operation;
	if ([[otherVariable class] isSubclassOfClass: [NSNumber class]]) {
		operation  = [[GLScalarMultiplyOperation alloc] initWithVectorOperand: self scalarOperand: [(NSNumber *) otherVariable doubleValue]];
	} else {
		if ([[self class] isSubclassOfClass: [GLScalar class]] || [[otherVariable class] isSubclassOfClass: [GLScalar class]]) {
			operation = [[GLMultiplicationOperation alloc] initWithFirstOperand: self secondOperand: otherVariable];
		} else if ([[self class] isSubclassOfClass: [GLFunction class]] && [[otherVariable class] isSubclassOfClass: [GLFunction class]]) {
			operation = [[GLMultiplicationOperation alloc] initWithFirstOperand: self secondOperand: otherVariable];
		} else if ([[self class] isSubclassOfClass: [GLFunction class]] && [[otherVariable class] isSubclassOfClass: [GLLinearTransform class]]) {
			[NSException raise: @"InvalidOperation" format: @"You cannot left-multipy a function by a linear transformation"];
		} else if ([[self class] isSubclassOfClass: [GLLinearTransform class]] && [[otherVariable class] isSubclassOfClass: [GLFunction class]]) {
			return [(GLLinearTransform *) self transform: otherVariable];
		} else if ([[self class] isSubclassOfClass: [GLLinearTransform class]] && [[otherVariable class] isSubclassOfClass: [GLLinearTransform class]]) {
			return [(GLLinearTransform *) self matrixMultiply: (GLLinearTransform *)otherVariable];
		}
	}
    
    if (operation) {
        operation = [self replaceWithExistingOperation: operation];
        return operation.result[0];
    }
	return nil;
}

- (GLVariable *) times: (id) otherVariable
{
	GLVariableOperation *operation;
	if ([[otherVariable class] isSubclassOfClass: [NSNumber class]]) {
		operation  = [[GLScalarMultiplyOperation alloc] initWithVectorOperand: self scalarOperand: [(NSNumber *) otherVariable doubleValue]];
	} else {
		if ([[self class] isSubclassOfClass: [GLScalar class]] || [[otherVariable class] isSubclassOfClass: [GLScalar class]]) {
			operation = [[GLMultiplicationOperation alloc] initWithFirstOperand: self secondOperand: otherVariable];
		} else if ([[self class] isSubclassOfClass: [GLFunction class]] && [[otherVariable class] isSubclassOfClass: [GLFunction class]]) {
			operation = [[GLMultiplicationOperation alloc] initWithFirstOperand: [(GLFunction*)self spatialDomain] secondOperand: [(GLFunction*)otherVariable spatialDomain]];
		} else if ([[self class] isSubclassOfClass: [GLFunction class]] && [[otherVariable class] isSubclassOfClass: [GLLinearTransform class]]) {
			[NSException raise: @"InvalidOperation" format: @"You cannot left-multipy a function by a linear transformation"];
		} else if ([[self class] isSubclassOfClass: [GLLinearTransform class]] && [[otherVariable class] isSubclassOfClass: [GLFunction class]]) {
			return [(GLLinearTransform *) self transform: otherVariable];
		} else if ([[self class] isSubclassOfClass: [GLLinearTransform class]] && [[otherVariable class] isSubclassOfClass: [GLLinearTransform class]]) {
			return [(GLLinearTransform *) self matrixMultiply: (GLLinearTransform *)otherVariable];
		}
	}
    if (operation) {
        operation = [self replaceWithExistingOperation: operation];
        return operation.result[0];
    }
	return nil;
}

/************************************************/
/*		Reading & Writing						*/
/************************************************/

#pragma mark -
#pragma mark Reading & Writing
#pragma mark

- (void) dumpToConsole
{
	return NSLog(@"%@", [NSString stringWithFormat: @"Uninitialized rank %lu tensor", self.rank]);
}

- (NSString *) graphvisDescription
{
    return [NSString stringWithFormat: @"Uninitialized rank %lu tensor", self.rank];
}


/************************************************/
/*		Private									*/
/************************************************/

#pragma mark -
#pragma mark Private
#pragma mark

- (NSArray *) pendingOperations
{
	NSArray *returnArray;
	@synchronized (self) {
		returnArray = [NSArray arrayWithArray: _pendingOperations];
	}
	return returnArray;
}

- (void) addOperation: (id) operation
{
	@synchronized (self) {
		if ( ![_pendingOperations containsObject: operation] ) {
			[_pendingOperations addObject: operation];
		}
	}
}

- (void) removeOperation: (id) operation
{
	@synchronized (self) {
		[_pendingOperations removeObject: operation];
	}
}

- (GLVariableOperation *) lastOperation
{
	GLVariableOperation *lastOp;
	@synchronized (self) {
		lastOp = [_pendingOperations lastObject];
	}
	return lastOp;
}

- (void) dealloc
{
	if (_data) {
		[[GLMemoryPool sharedMemoryPool] returnData: _data];
	}
}

- (id) replaceWithExistingOperation: (GLVariableOperation *) newOperation
{
    GLVariableOperation *existingOperation;
    for (GLVariableOperation *op in self.existingOperations) {
        if ( [op isEqualToOperation: newOperation] ) {
            existingOperation = op;
        }
    }
    
    if (existingOperation) {
        return existingOperation;
    } else {
        [self.existingOperations addObject: newOperation];
    }
    return newOperation;
}

@end
