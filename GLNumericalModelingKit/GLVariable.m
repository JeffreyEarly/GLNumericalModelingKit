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
#import "GLUnaryOperations.h"
#include <mach/mach_time.h>

@interface GLVariable ()
@property(readwrite, strong, nonatomic) NSMutableDictionary *metadata;
@property(readwrite, assign, nonatomic) NSUInteger uniqueID;
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

static NSString *GLVariableEquationKey = @"GLVariableEquationKey";
static NSString *GLVariableNameKey = @"GLVariableNameKey";
static NSString *GLVariableUnitsKey = @"GLVariableUnitsKey";
static NSString *GLVariableMetadataKey = @"GLVariableMetadataKey";
static NSString *GLVariableUniqueIDKey = @"GLVariableUniqueIDKey";
static NSString *GLVariableIsImaginaryPartZeroKey = @"GLVariableIsImaginaryPartZeroKey";
static NSString *GLVariableIsRealPartZeroKey = @"GLVariableIsRealPartZeroKey";
static NSString *GLVariableNDataPointsKey = @"GLVariableNDataPointsKey";
static NSString *GLVariableNDataElementsKey = @"GLVariableNDataElementsKey";
static NSString *GLVariableDataFormatKey = @"GLVariableDataFormatKey";
static NSString *GLVariableDataBytesKey = @"GLVariableDataBytesKey";
static NSString *GLVariableDataKey = @"GLVariableDataKey";
static NSString *GLVariableMatrixDescriptionKey = @"GLVariableMatrixDescriptionKey";

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

+ (GLVariable *) variableOfRealTypeWithPrototype: (GLVariable *) variable
{
    if (variable.rank == 0) {
        GLScalar *scalar = (GLScalar *) variable;
        return [[GLScalar alloc] initWithType: kGLRealDataFormat forEquation:scalar.equation];
    } else if (variable.rank == 1) {
        GLFunction *function = (GLFunction *) variable;
        return [[function class] functionOfType: kGLRealDataFormat withDimensions: function.dimensions forEquation: function.equation];
    }  else if (variable.rank == 2) {
        GLLinearTransform *matrix = (GLLinearTransform *) variable;
        return [GLLinearTransform transformOfType: kGLRealDataFormat withFromDimensions: matrix.fromDimensions toDimensions: matrix.toDimensions inFormat: matrix.matrixFormats forEquation:matrix.equation matrix:nil];
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

- (void)encodeWithCoder:(NSCoder *)coder
{
    [coder encodeObject: self.equation forKey:GLVariableEquationKey];
    
    if (self.name) [coder encodeObject:self.name forKey:GLVariableNameKey];
    if (self.units) [coder encodeObject:self.units forKey:GLVariableUnitsKey];
    if (_metadata) [coder encodeObject:_metadata forKey:GLVariableMetadataKey];
    
    [coder encodeObject: @(self.uniqueID) forKey:GLVariableUniqueIDKey];
    [coder encodeObject: @(self.isImaginaryPartZero) forKey:GLVariableIsImaginaryPartZeroKey];
    [coder encodeObject: @(self.isRealPartZero) forKey:GLVariableIsRealPartZeroKey];
    
    [coder encodeObject: @(self.nDataPoints) forKey:GLVariableNDataPointsKey];
    [coder encodeObject: @(self.nDataElements) forKey:GLVariableNDataElementsKey];
    [coder encodeObject: @(self.dataFormat) forKey:GLVariableDataFormatKey];
    
    [coder encodeObject: @(self.dataBytes) forKey:GLVariableDataBytesKey];
    [coder encodeObject: self.data forKey:GLVariableDataKey];
    
    [coder encodeObject: self.matrixDescription forKey:GLVariableMatrixDescriptionKey];
}

- (id)initWithCoder:(NSCoder *)decoder
{
    if ((self=[super init])) {
        _equation = [decoder decodeObjectForKey:GLVariableEquationKey];
        _name = [decoder decodeObjectForKey: GLVariableNameKey];
        _units = [decoder decodeObjectForKey: GLVariableUnitsKey];
        _metadata = [decoder decodeObjectForKey: GLVariableMetadataKey];
        
        _uniqueID = [[decoder decodeObjectForKey: GLVariableUniqueIDKey] unsignedIntegerValue];
        _isImaginaryPartZero = [[decoder decodeObjectForKey: GLVariableIsImaginaryPartZeroKey] boolValue];
        _isRealPartZero = [[decoder decodeObjectForKey: GLVariableIsRealPartZeroKey] boolValue];
        
        _nDataPoints = [[decoder decodeObjectForKey: GLVariableNDataPointsKey] unsignedIntegerValue];
        _nDataElements = [[decoder decodeObjectForKey: GLVariableNDataElementsKey] unsignedIntegerValue];
        _dataFormat = [[decoder decodeObjectForKey: GLVariableDataFormatKey] unsignedIntegerValue];
        
        _dataBytes = [[decoder decodeObjectForKey: GLVariableDataBytesKey] unsignedIntegerValue];
        _data = [decoder decodeObjectForKey: GLVariableDataKey];
        
        _matrixDescription = [decoder decodeObjectForKey: GLVariableMatrixDescriptionKey];
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

- (void) setData: (NSMutableData *) newData
{
	if (newData.length < self.dataBytes) {
		[NSException raise:@"InsufficientBufferSize" format:@"You are setting a buffer that is too small"];
	}
	_data=newData;
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

- (GLFloatComplex *) floatComplex
{
	[self solve];
    GLFloatComplex *f = self.data.mutableBytes;
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

- (instancetype) interleavedFormat
{
    if (!self.isComplex) {
        [NSException raise: @"NotYetImplemented" format: @"Not yet able to force interleaved format from real format"];
        return nil;
    }
    
    if (self.dataFormat == kGLInterleavedComplexDataFormat) {
        return self;
    } else {
        GLVariableOperation *operation = [[GLSplitToInterleavedComplexOperation alloc] initWithVariable: self];
        operation = [self replaceWithExistingOperation: operation];
        return operation.result[0];
    }
}

- (instancetype) splitFormat
{
    if (!self.isComplex) {
        [NSException raise: @"NotYetImplemented" format: @"Not yet able to force split format from real format"];
        return nil;
    }
    
    if (self.dataFormat == kGLSplitComplexDataFormat) {
        return self;
    } else {
        GLVariableOperation *operation = [[GLInterleavedToSplitComplexOperation alloc] initWithVariable: self];
        operation = [self replaceWithExistingOperation: operation];
        return operation.result[0];
    }
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

- (instancetype) makeRealIfPossible
{
    if (!self.isComplex) {
        return self;
    }
    
    BOOL isReal = YES;
    GLFloat *value = self.pointerValue;
    NSUInteger es = self.matrixDescription.elementStride;
    NSUInteger cs = self.matrixDescription.complexStride;
    GLFloat prec = sizeof(GLFloat) == 4 ? 1e-6 : 1e-14;
    for (NSUInteger i=0; i<self.matrixDescription.nPoints; i++ ) {
        GLFloat imagp = fabs(value[i*es+cs]);
        isReal &= (imagp < prec || imagp/cabs(value[i*es] + value[i*es+cs]*I) < prec);
        if (!isReal) {
            break;
        }
    }
    
    if (isReal) {
        GLVariable *result = [GLVariable variableOfRealTypeWithPrototype: self];
        GLFloat *b = result.pointerValue;
        for (NSUInteger i=0; i<self.matrixDescription.nPoints; i++ ) {
            b[i] = value[i*es];
        }
        return result;
    } else {
        return self;
    }
}

/************************************************/
/*		Reading & Writing						*/
/************************************************/

#pragma mark -
#pragma mark Reading & Writing
#pragma mark

void createMatrixDescriptionString( GLMatrixDescription *matrixDescription, NSMutableString *descrip, GLFloat *a, NSUInteger rank)
{
	NSUInteger lastNonTrivialDimension = NSNotFound;
	for (NSUInteger i=0; i<matrixDescription.nDimensions; i++) {
		if (matrixDescription.strides[i].matrixFormat != kGLIdentityMatrixFormat) {
			lastNonTrivialDimension = i;
		}
	}
	
	GLFloat max, min;
	vGL_maxv( a, 1, &max, matrixDescription.nElements);
	vGL_minv( a, 1, &min, matrixDescription.nElements);
	
	if ( fabs(min) > max) {
		max = fabs(min);
	}
	
	GLFloat divisor = pow(10, floor(log10(max)));
	if ( divisor == 0.0) divisor = 1;
	
	[descrip appendFormat: @"%g * \n", divisor];
	createMatrixDescriptionStringRecursion(matrixDescription, descrip, @"(", @"(", 0, a, divisor, lastNonTrivialDimension, rank);
}

void createMatrixDescriptionStringRecursion( GLMatrixDescription *matrixDescription, NSMutableString *descrip, NSString *leftPos, NSString *rightPos, NSUInteger iDim, GLFloat *a, GLFloat divisor, NSUInteger lastNonTrivialDimension, NSUInteger rank)
{
	NSString *index;
	if (iDim==0) {
		index=@"i";
	} else if (iDim==1) {
		index=@"j";
	} else if (iDim==2) {
		index=@"k";
	} else {
		index=@"m";
	}
	NSString *newLeftPos = leftPos;
	NSString *newRightPos = rightPos;
	
    if (matrixDescription.strides[iDim].matrixFormat == kGLIdentityMatrixFormat) {
        newLeftPos = [leftPos stringByAppendingString:index];
		newRightPos = [rightPos stringByAppendingString:index];
		
		if (iDim+1<matrixDescription.nDimensions) {
			newLeftPos = [newLeftPos stringByAppendingString:@","];
			newRightPos = [newRightPos stringByAppendingString:@","];
		} else if (iDim+1==matrixDescription.nDimensions) {
			newLeftPos = [newLeftPos stringByAppendingString:@")"];
			newRightPos = [newRightPos stringByAppendingString:@")"];
		}
		
		if (iDim+1 < matrixDescription.nDimensions) {
			createMatrixDescriptionStringRecursion(matrixDescription, descrip, newLeftPos, newRightPos, iDim+1, a, divisor, lastNonTrivialDimension, rank);
		} else {
			[descrip appendString: newLeftPos];
			[descrip appendString: newRightPos];
		}
    } else if (matrixDescription.strides[iDim].matrixFormat == kGLDenseMatrixFormat) {
        NSUInteger rs = matrixDescription.strides[iDim].rowStride;
        NSUInteger cs = matrixDescription.strides[iDim].columnStride;
        NSUInteger is = matrixDescription.complexStride;
        
		NSString *columnEndCharacter = matrixDescription.strides[iDim].nColumns > 1 ? @"\n" : @"";
        
        for (NSUInteger i=0; i<matrixDescription.strides[iDim].nRows; i++) {
			
			if (rank == 1 && iDim+1<matrixDescription.nDimensions) {
				if (iDim+2<matrixDescription.nDimensions) {
					newLeftPos = [leftPos stringByAppendingFormat:@"%lu,", i];
				} else if (iDim+2==matrixDescription.nDimensions) {
					newLeftPos = [leftPos stringByAppendingFormat:@"%lu,0:%lu)\t", i, matrixDescription.strides[iDim+1].nRows-1];
				}
			}
			
            for (NSUInteger j=0; j<matrixDescription.strides[iDim].nColumns; j++) {
                
                if (iDim == lastNonTrivialDimension) {
					if (iDim+1 < matrixDescription.nDimensions) {
						createMatrixDescriptionStringRecursion(matrixDescription, descrip, newLeftPos, newRightPos, iDim+1, &(a[i*rs+j*cs]), divisor, lastNonTrivialDimension, rank);
					}
                    if (!is) [descrip appendFormat: @"%6.6f\t", a[i*rs+j*cs]/divisor];
                    else [descrip appendFormat: @"%6.2f + %6.2fi\t", a[i*rs+j*cs]/divisor, a[i*rs+j*cs + is]/divisor];
                } else {
                    [descrip appendFormat: @"\n"];
					if (rank == 1 && iDim+2==matrixDescription.nDimensions) [descrip appendString: newLeftPos];
					newLeftPos = [leftPos stringByAppendingFormat:@"%lu,",i];
					newRightPos = [rightPos stringByAppendingFormat:@"%lu,", j];
                    createMatrixDescriptionStringRecursion(matrixDescription, descrip, newLeftPos, newRightPos, iDim+1, &(a[i*rs+j*cs]), divisor, lastNonTrivialDimension, rank);
                }
                
            }
            [descrip appendString:columnEndCharacter];
        }
        
    } else {
        
		NSString *diagonalEndCharacter = matrixDescription.strides[iDim].nDiagonals > 1 ? @"\n" : @"";
		
        NSUInteger ds = matrixDescription.strides[iDim].diagonalStride;
        NSUInteger es = matrixDescription.strides[iDim].stride;
        NSUInteger is = matrixDescription.complexStride;
        
		if (iDim==lastNonTrivialDimension) {
			if (matrixDescription.strides[iDim].matrixFormat == kGLDiagonalMatrixFormat) {
				newLeftPos = [leftPos stringByAppendingFormat:@"%@=0:%lu",index,matrixDescription.strides[iDim].nDiagonalPoints-1];
				newRightPos = [rightPos stringByAppendingFormat:@"%@", index];
			} else if (matrixDescription.strides[iDim].matrixFormat == kGLSubdiagonalMatrixFormat) {
				newLeftPos = [leftPos stringByAppendingFormat:@"%@=0:%lu",index,matrixDescription.strides[iDim].nDiagonalPoints-1];
				newRightPos = [rightPos stringByAppendingFormat:@"%@-1", index];
			} else if (matrixDescription.strides[iDim].matrixFormat == kGLSuperdiagonalMatrixFormat) {
				newLeftPos = [leftPos stringByAppendingFormat:@"%@=0:%lu",index,matrixDescription.strides[iDim].nDiagonalPoints-1];
				newRightPos = [rightPos stringByAppendingFormat:@"%@+1", index];
			}
			if (iDim+1<matrixDescription.nDimensions) {
				newLeftPos = [newLeftPos stringByAppendingString:@","];
				newRightPos = [newRightPos stringByAppendingString:@","];
				createMatrixDescriptionStringRecursion(matrixDescription, descrip, newLeftPos, newRightPos, iDim+1, a, divisor, lastNonTrivialDimension, rank);
			} else if (iDim+1==matrixDescription.nDimensions) {
				newLeftPos = [newLeftPos stringByAppendingString:@")"];
				newRightPos = [newRightPos stringByAppendingString:@")"];
				[descrip appendString: newLeftPos];
				[descrip appendString: newRightPos];
			}
		}
		
        for (NSUInteger i=0; i<matrixDescription.strides[iDim].nDiagonalPoints; i++) {
			
			if (iDim==lastNonTrivialDimension && matrixDescription.strides[iDim].matrixFormat == kGLTridiagonalMatrixFormat) {
				if (iDim==lastNonTrivialDimension) {
					newLeftPos = [leftPos stringByAppendingFormat:@"%lu",i];
					newRightPos = [rightPos stringByAppendingFormat:@"%ld:%ld", (NSInteger)i-1, i+1];
				}
				if (iDim+1<matrixDescription.nDimensions) {
					newLeftPos = [newLeftPos stringByAppendingString:@","];
					newRightPos = [newRightPos stringByAppendingString:@","];
					createMatrixDescriptionStringRecursion(matrixDescription, descrip, newLeftPos, newRightPos, iDim+1, a, divisor, lastNonTrivialDimension, rank);
				} else if (iDim+1==matrixDescription.nDimensions) {
					newLeftPos = [newLeftPos stringByAppendingString:@")"];
					newRightPos = [newRightPos stringByAppendingString:@")"];
					[descrip appendString: newLeftPos];
					[descrip appendString: newRightPos];
				}
			}
			
            for (NSUInteger iDiagonal=0; iDiagonal<matrixDescription.strides[iDim].nDiagonals; iDiagonal++) {
                if (iDim == lastNonTrivialDimension) {
                    if (!is) [descrip appendFormat: @"%6.2f\t", a[iDiagonal*ds+i*es]/divisor];
                    else [descrip appendFormat: @"%6.2f + %6.2fi\t", a[iDiagonal*ds+i*es]/divisor, a[iDiagonal*ds+i*es+is]/divisor];
                } else {
                    [descrip appendFormat: @"\n"];
					if (matrixDescription.strides[iDim].matrixFormat == kGLDiagonalMatrixFormat) {
						newLeftPos = [leftPos stringByAppendingFormat:@"%lu,",i];
						newRightPos = [rightPos stringByAppendingFormat:@"%lu,", i];
					} else if (matrixDescription.strides[iDim].matrixFormat == kGLTridiagonalMatrixFormat) {
						newLeftPos = [leftPos stringByAppendingFormat:@"%lu,",i];
						newRightPos = [rightPos stringByAppendingFormat:@"%ld,", i+iDiagonal-1];
					} else if (matrixDescription.strides[iDim].matrixFormat == kGLSubdiagonalMatrixFormat) {
						newLeftPos = [leftPos stringByAppendingFormat:@"%lu,",i];
						newRightPos = [rightPos stringByAppendingFormat:@"%ld,", i-1];
					} else if (matrixDescription.strides[iDim].matrixFormat == kGLSuperdiagonalMatrixFormat) {
						newLeftPos = [leftPos stringByAppendingFormat:@"%lu,",i];
						newRightPos = [rightPos stringByAppendingFormat:@"%ld,", i+1];
					}
                    createMatrixDescriptionStringRecursion(matrixDescription, descrip, newLeftPos, newRightPos, iDim+1, &(a[iDiagonal*ds+i*es]), divisor, lastNonTrivialDimension, rank);
                }
            }
            [descrip appendString:diagonalEndCharacter];
        }
    }
}

- (NSString *) matrixDescriptionString
{
	NSMutableString *descrip = [NSMutableString string];
	createMatrixDescriptionString(self.matrixDescription, descrip, self.pointerValue, self.rank);
	return descrip;
}

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
