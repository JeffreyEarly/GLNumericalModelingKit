//
//  GLVariable.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 4/21/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import "GLVariable.h"
#import "GLEquation.h"
#import "GLDimension.h"

#import "GLVariableOperations.h"
#import "GLNullaryOperations.h"
#import "GLUnaryOperations.h"
#import "GLVectorScalarOperations.h"
#import "GLVectorVectorOperations.h"
#import "GLBasisTransformationOperations.h"
#import "GLInterpolationOperations.h"
#import "GLSimpleInterpolationOperations.h"
#import "GLDimensionalOperations.h"

#import "GLSpectralDifferentialOperatorPool.h"
#import "GLMemoryPool.h"
#import "GLNetCDFFile.h"

#import "GLMatrixDescription.h"

#include <mach/mach_time.h>

#define MatrixFFT 0

/************************************************/
/*		GLVariable								*/
/************************************************/

#pragma mark -
#pragma mark GLVariable
#pragma mark

@interface GLFunction ()
@property(readwrite, assign, nonatomic) NSUInteger nDataPoints;
@property(readwrite, assign, nonatomic) NSUInteger nDataElements;
@property(readwrite, assign, nonatomic) NSUInteger dataBytes;
@end

static BOOL prefersSpatialMultiplication = YES;

@implementation GLVariable

/************************************************/
/*		Class Methods							*/
/************************************************/

#pragma mark -
#pragma mark Class Methods
#pragma mark

+ (void) setPrefersSpatialMultiplication: (BOOL) aFlag
{
	prefersSpatialMultiplication = aFlag;
}

+ (BOOL) prefersSpatialMultiplication
{
	return prefersSpatialMultiplication;
}

// To convert from an element to an index: vec[(i*ny+j)*nz+k] = matrix[i][j][k]
+ (NSArray *) indicesFromElement: (NSUInteger) theIndex forDimensions: (NSArray *) theDimensions
{
	NSMutableArray *indices = [NSMutableArray arrayWithCapacity: theDimensions.count];
	
	NSUInteger theRemainder = 0;
	for (GLDimension *aDim in theDimensions.reverseObjectEnumerator)
	{
		theRemainder = theIndex % aDim.nPoints;
		[indices addObject: [NSNumber numberWithUnsignedInteger: theRemainder]];
		theIndex = (theIndex - theRemainder)/aDim.nPoints;
	}
	
	return indices.reverseObjectEnumerator.allObjects;
}

// Returns a variable built from the given dimension
+ (id) variableOfRealTypeFromDimension: (GLDimension *) aDimension withDimensions: (NSArray *) theDimensions forEquation: (GLEquation *) equation
{
	GLVariable *aVariable = [self variableOfRealTypeWithDimensions: theDimensions forEquation: equation];
	aVariable.name = aDimension.name;
	aVariable.units = aDimension.units;
	GLFloat *f = aVariable.pointerValue;
	
	NSUInteger indexOfKeyDimension = [aVariable.dimensions indexOfObject: aDimension];
	NSUInteger modulus = aDimension.nPoints;
	NSUInteger divisor = 1;
	for (NSUInteger i=indexOfKeyDimension+1; i<theDimensions.count; i++)
	{
		GLDimension *aDim = [aVariable.dimensions objectAtIndex: i];
		divisor *= aDim.nPoints;
	}
	
	GLFloat *dimValue = aDimension.data.mutableBytes;
	for (NSUInteger i=0; i<aVariable.nDataPoints; i++)
	{
		NSUInteger dimIndex = (i/divisor) % modulus;
		f[i] = dimValue[dimIndex];
	}	
	
	return aVariable;
}

+ (id) variableOfComplexTypeFromDimension: (GLDimension *) aDimension withDimensions: (NSArray *) theDimensions forEquation: (GLEquation *) equation
{
	GLVariable *aVariable = [self variableOfComplexTypeWithDimensions: theDimensions forEquation: equation];
	aVariable.name = aDimension.name;
	aVariable.units = aDimension.units;
	GLSplitComplex f= aVariable.splitComplex;
	
	NSUInteger indexOfKeyDimension = [aVariable.dimensions indexOfObject: aDimension];
	NSUInteger modulus = aDimension.nPoints;
	NSUInteger divisor = 1;
	for (NSUInteger i=indexOfKeyDimension+1; i<theDimensions.count; i++)
	{
		GLDimension *aDim = [aVariable.dimensions objectAtIndex: i];
		divisor *= aDim.nPoints;
	}
	
	GLFloat *dimValue = aDimension.data.mutableBytes;
	for (NSUInteger i=0; i<aVariable.nDataPoints; i++)
	{
		NSUInteger dimIndex = (i/divisor) % modulus;
		f.realp[i] = dimValue[dimIndex];
        f.imagp[i] = 0.0;
	}	
	
	return aVariable;
}

//+ (id) variableOfHalfComplexTypeFromDimension: (GLDimension *) aDimension withDimensions: (NSArray *) theDimensions forEquation: (GLEquation *) equation
//{
//	GLVariable *aVariable = [self variableOfType: kGLHalfComplexDataFormat withDimensions: theDimensions forEquation: equation];
//	aVariable.name = aDimension.name;
//	aVariable.units = aDimension.units;
//	GLFloat *f = aVariable.pointerValue;
//	
//	NSUInteger indexOfKeyDimension = [aVariable.dimensions indexOfObject: aDimension];
//	NSUInteger modulus = 2*(aDimension.nPoints-1);
//    NSUInteger n = aDimension.nPoints-1;
//	NSUInteger divisor = 1;
//	for (NSUInteger i=indexOfKeyDimension+1; i<theDimensions.count; i++)
//	{
//		GLDimension *aDim = [aVariable.dimensions objectAtIndex: i];
//		divisor *= 2*(aDim.nPoints-1);
//	}
//	
//	GLFloat *dimValue = aDimension.data.mutableBytes;
//	for (NSUInteger i=0; i<aVariable.nDataElements; i++)
//	{
//		NSUInteger dimIndex = (i/divisor) % modulus;
//        if (dimIndex <= n) {
//            f[i] = dimValue[dimIndex];
//        } else {
//            f[i] = dimValue[2*n-dimIndex];
//        }
//	}	
//	
//	return aVariable;
//}

+ (id) variableOfType: (GLDataFormat) format fromDimension: (GLDimension *) aDimension withDimensions: (NSArray *) theDimensions forEquation: (GLEquation *) equation
{
    if (format == kGLRealDataFormat) {
        return [self variableOfRealTypeFromDimension: aDimension withDimensions: theDimensions forEquation: equation];
    } else if (format == kGLSplitComplexDataFormat) {
        return [self variableOfComplexTypeFromDimension: aDimension withDimensions: theDimensions forEquation: equation];
    }
    
	GLVariable *aVariable = [self variableOfType: format withDimensions: theDimensions forEquation: equation];
	aVariable.name = aDimension.name;
	aVariable.units = aDimension.units;
	GLFloat *f = aVariable.pointerValue;
	
	NSUInteger indexOfKeyDimension = [aVariable.dimensions indexOfObject: aDimension];
	NSUInteger modulus = 2*(aDimension.nPoints-1);
    NSUInteger n = aDimension.nPoints-1;
	NSUInteger divisor = 1;
	for (NSUInteger i=indexOfKeyDimension+1; i<theDimensions.count; i++)
	{
		GLDimension *aDim = [aVariable.dimensions objectAtIndex: i];
		divisor *= 2*(aDim.nPoints-1);
	}
	
	GLFloat *dimValue = aDimension.data.mutableBytes;
	for (NSUInteger i=0; i<aVariable.nDataElements; i++)
	{
		NSUInteger dimIndex = (i/divisor) % modulus;
        if (dimIndex <= n) {
            f[i] = dimValue[dimIndex];
        } else {
            f[i] = dimValue[2*n-dimIndex];
        }
	}	
	
	return aVariable;
}

+ (id) variableWithRandomValuesBetween: (GLFloat) min and: (GLFloat) max withDimensions: (NSArray *) theDimensions forEquation: (GLEquation *) equation
{
	GLVariable *var = [self variableOfRealTypeWithDimensions: theDimensions forEquation: equation];
	GLRandomNumberOperation *rand = [[GLRandomNumberOperation alloc] initWithResult: var firstScalarOperand: min secondScalarOperand: max];
	return rand.result;
}

+ (id) variableOfRealTypeWithDimensions: (NSArray *) theDimensions forEquation: (GLEquation *) equation
{
	return [self variableOfType: kGLRealDataFormat withDimensions: theDimensions forEquation: equation];
}

+ (id) variableOfComplexTypeWithDimensions: (NSArray *) theDimensions forEquation: (GLEquation *) equation
{
	return [self variableOfType: kGLSplitComplexDataFormat withDimensions: theDimensions forEquation: equation];
}

+ (id) variableOfType: (GLDataFormat) dataFormat withDimensions: (NSArray *) theDimensions forEquation: (GLEquation *) equation;
{	
	BOOL isMutable = NO;
	for (GLDimension *dimension in theDimensions) {
		isMutable |= dimension.isMutable;
	}
	
	Class GLVariableClass = isMutable ? [GLMutableVariable class] : self;
	
	return [[GLVariableClass alloc] initVariableOfType: dataFormat withDimensions: theDimensions forEquation: equation];
}

// Copies the data from the other variable (now!) not a delayed operation.
+ (id) variableFromVariable: (GLVariable *) otherVariable
{
    [otherVariable solve];
    
	GLVariable *newVar = [[[self class] alloc] initVariableOfType: otherVariable.dataFormat withDimensions:otherVariable.dimensions forEquation:otherVariable.equation];
	newVar.data = otherVariable.data;
    for ( NSUInteger i=0; i<otherVariable.dimensions.count; i++) {
        newVar.realSymmetry[i] = otherVariable.realSymmetry[i];
        newVar.imaginarySymmetry[i] = otherVariable.imaginarySymmetry[i];
    }
	
	return newVar;
}

- (NSUInteger) rank {
	return 1;
}

/************************************************/
/*		Initialization							*/
/************************************************/

#pragma mark -
#pragma mark Initialization
#pragma mark

//- (id) init
//{
//	[NSException raise: @"BadInitialization" format: @"Cannot initialize GLVariable with -init method. Use -initComplexVariable:withDimensions:forEquation instead."];
//	
//	return self;
//}

@synthesize nDataPoints = _nDataPoints;
@synthesize nDataElements = _nDataElements;
@synthesize dataBytes = _dataBytes;

- (id) initVariableOfType: (GLDataFormat) dataFormat withDimensions: (NSArray *) theDimensions forEquation: (GLEquation *) theEquation
{
	if (!theDimensions) {
		NSLog(@"Attempted to initialize GLVariable without an equation or dimensions!!!");
		return nil;
	}
	
	if ((self = [super initWithType:dataFormat withEquation:theEquation])) {
		_dimensions = theDimensions;
		
		_isFrequencyDomain = 0;
		
		// We loop through the dimensions and allocate enough memory for the variable
		// defined on each dimension.
		_nDataPoints = 1;
		_nDataElements = 1;
		_realSymmetry = [NSMutableArray array];
		_imaginarySymmetry = [NSMutableArray array];
		for ( GLDimension *aDim in theDimensions ) {
            NSUInteger idx = [theDimensions indexOfObject: aDim];
            
            _nDataElements *= aDim.nPoints;
            _nDataPoints *= aDim.nPoints;
            
            if (aDim.basisFunction == kGLDeltaBasis) {
                self.realSymmetry[idx] = @(kGLNoSymmetry);
                self.imaginarySymmetry[idx] = (dataFormat == kGLRealDataFormat ? @(kGLZeroSymmetry) : @(kGLNoSymmetry));
            } else if (aDim.basisFunction == kGLCosineBasis || aDim.basisFunction == kGLCosineHalfShiftBasis) {
                self.realSymmetry[idx] = @(kGLEvenSymmetry);
                self.imaginarySymmetry[idx] = @(kGLZeroSymmetry);
            } else if (aDim.basisFunction == kGLSineBasis || aDim.basisFunction == kGLSineHalfShiftBasis) {
                self.realSymmetry[idx] = @(kGLOddSymmetry);
                self.imaginarySymmetry[idx] = @(kGLZeroSymmetry);
            }
            
			_isFrequencyDomain |= aDim.isFrequencyDomain;
		}
        
        if (dataFormat == kGLSplitComplexDataFormat) {
            _nDataElements *= 2;
        }
		
		_dataBytes = _nDataElements*sizeof(GLFloat);
        
        self.matrixDescription = [[GLMatrixDescription alloc] initWithVariable: self];
	}	
	
	return self;
}

/************************************************/
/*		Dimensionality							*/
/************************************************/

#pragma mark -
#pragma mark Dimensionality
#pragma mark

- (BOOL) isMutable {
	return NO;
}

- (BOOL) isHermitian {
	if (self.hermitianDimension) {
		NSUInteger i = [self.dimensions indexOfObject: self.hermitianDimension];
		return ( ([self.realSymmetry[i] unsignedIntegerValue] == kGLEvenSymmetry || self.isRealPartZero) && ([self.imaginarySymmetry[i] unsignedIntegerValue] == kGLOddSymmetry || self.isImaginaryPartZero) );
	}
	return NO;	
}

- (NSArray *) basis
{
	NSMutableArray *basis = [NSMutableArray arrayWithCapacity: self.dimensions.count];
	for (GLDimension *dim in self.dimensions) {
		[basis addObject: @(dim.basisFunction)];
	}
	return basis;
}

- (NSArray *) dimensionsTransformedToBasis: (NSArray *) basis
{
	if ( basis.count == 1 && self.dimensions.count > 1) {
        NSMutableArray *array = [NSMutableArray array];
        for (NSUInteger i=0; i < self.dimensions.count; i++) {
            [array addObject: basis.lastObject];
        }
        basis = array;
    }
	
	
	NSMutableArray *transformedDimensions = [NSMutableArray array];
	NSUInteger lastExponentialIndex = NSNotFound;
	
#if MatrixFFT
    
#else
	if (self.isPurelyReal) {
		for (NSUInteger i=0; i < self.dimensions.count; i++) {
			if ([basis[i] unsignedIntegerValue] == kGLExponentialBasis) {
				lastExponentialIndex = i;
			}
		}		
	}
#endif
    
	for (NSUInteger i=0; i < self.dimensions.count; i++) {
		BOOL strictlyPositive = ( lastExponentialIndex == i || [basis[i] unsignedIntegerValue] > 1) ? YES : NO;
		GLDimension *dim = [[GLDimension alloc] initAsDimension: self.dimensions[i] transformedToBasis: [basis[i] unsignedIntegerValue] strictlyPositive: strictlyPositive];
		[transformedDimensions addObject: dim];
	}
	
	return transformedDimensions;
}

/************************************************/
/*		Data									*/
/************************************************/

#pragma mark -
#pragma mark Data
#pragma mark

- (void) zero
{
	vGL_vclr( self.pointerValue, 1, self.nDataElements);
    
    for (NSUInteger i=0; i<self.dimensions.count; i++) {
        self.realSymmetry[i] = @(kGLZeroSymmetry);
        self.imaginarySymmetry[i] = @(kGLZeroSymmetry);
    }
}

- (void) rand: (GLFloat) amp
{
	for (NSUInteger i=0; i<self.nDataElements; i++) {
		self.pointerValue[i] = 2*amp*((double) rand())/( (double) RAND_MAX ) - amp;
	}
}

- (BOOL) isFinite
{
	for (NSUInteger i=0; i<self.nDataElements; i++) {
		if ( !isfinite(self.pointerValue[i])) return NO;
	}
	return YES;
}

- (GLFloat) maxNow
{
	GLFloat max = 0.0;
	vGL_maxv( (float *) self.data.bytes, 1, &max, self.nDataElements);
	return max;
}



/************************************************/
/*		Operations								*/
/************************************************/

#pragma mark -
#pragma mark Operations
#pragma mark

- (GLVariable *) plus: (GLVariable *) otherVariable
{
	GLAdditionOperation *operation = [[GLAdditionOperation alloc] initWithFirstOperand: self secondOperand: otherVariable];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result;
}

- (GLVariable *) minus: (GLVariable *) otherVariable
{
	GLSubtractionOperation *operation = [[GLSubtractionOperation alloc] initWithFirstOperand: self secondOperand: otherVariable];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result;
}

- (GLVariable *) times: (GLVariable *) otherVariable
{
	GLMultiplicationOperation *operation;
	if (prefersSpatialMultiplication) {
		operation = [[GLMultiplicationOperation alloc] initWithFirstOperand: [self spatialDomain] secondOperand: [otherVariable spatialDomain]];
	} else {
		operation = [[GLMultiplicationOperation alloc] initWithFirstOperand: self secondOperand: otherVariable];
	}
    operation = [self replaceWithExistingOperation: operation];
	return operation.result;	
}

- (GLVariable *) multiply: (GLVariable *) otherVariable
{
	GLMultiplicationOperation *operation = [[GLMultiplicationOperation alloc] initWithFirstOperand: self secondOperand: otherVariable];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result;	
}

- (GLVariable *) dividedBy: (GLVariable *) otherVariable
{
	GLDivisionOperation *operation = [[GLDivisionOperation alloc] initWithFirstOperand: self secondOperand: otherVariable];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result;
}

- (GLVariable *) dot: (GLVariable *) otherVariable
{
	GLDotProductOperation *operation = [[GLDotProductOperation alloc] initWithFirstOperand: self secondOperand: otherVariable];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result;	
}

- (GLVariable *) absMax: (GLVariable *) otherVariable
{
	GLAbsoluteLargestOperation *operation = [[GLAbsoluteLargestOperation alloc] initWithFirstOperand: self secondOperand: otherVariable];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result;	
}

- (GLVariable *) exponentiate
{	
	GLExponentialOperation *operation = [[GLExponentialOperation alloc] initWithOperand: self];
    operation = [self replaceWithExistingOperation: operation];
    return operation.result;
}

- (GLVariable *) sin
{	
	GLSineOperation *operation = [[GLSineOperation alloc] initWithOperand: self];
    operation = [self replaceWithExistingOperation: operation];
    return operation.result;
}

- (GLVariable *) cos
{	
	GLCosineOperation *operation = [[GLCosineOperation alloc] initWithOperand: self];
    operation = [self replaceWithExistingOperation: operation];
    return operation.result;
}

- (GLVariable *) atan
{
	GLInverseTangentOperation *operation = [[GLInverseTangentOperation alloc] initWithOperand: self];
    operation = [self replaceWithExistingOperation: operation];
    return operation.result;
}

- (GLVariable *) sqrt
{	
	GLSquareRootOperation *operation = [[GLSquareRootOperation alloc] initWithOperand: self];
    operation = [self replaceWithExistingOperation: operation];
    return operation.result;
}

- (GLVariable *) negate
{
	GLNegationOperation *operation = [[GLNegationOperation alloc] initWithOperand: self];
    operation = [self replaceWithExistingOperation: operation];
    return operation.result;
}

- (GLVariable *) abs
{
	GLAbsoluteValueOperation *operation = [[GLAbsoluteValueOperation alloc] initWithOperand: self];
    operation = [self replaceWithExistingOperation: operation];
    return operation.result;
}

- (GLVariable *) scalarAdd: (GLFloat) aScalar
{
	GLScalarAddOperation *operation = [[GLScalarAddOperation alloc] initWithVectorOperand: self scalarOperand: aScalar];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result;
}

- (GLVariable *) scalarMultiply: (GLFloat) aScalar
{	
	GLScalarMultiplyOperation *operation = [[GLScalarMultiplyOperation alloc] initWithVectorOperand: self scalarOperand: aScalar];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result;
}

- (GLVariable *) scalarDivide: (GLFloat) aScalar
{
	GLScalarDivideOperation *operation = [[GLScalarDivideOperation alloc] initWithVectorOperand: self scalarOperand: aScalar];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result;	
}

- (GLVariable *) scalarMax: (GLFloat) aScalar
{
	GLScalarThresholdOperation *operation = [[GLScalarThresholdOperation alloc] initWithVectorOperand: self scalarOperand: aScalar];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result;	
}

- (GLVariable *) zeroThreshold: (GLFloat) aScalar
{
	GLZeroThresholdOperation *operation = [[GLZeroThresholdOperation alloc] initWithVectorOperand: self scalarOperand: aScalar];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result;
}

- (id) clipToMin: (GLFloat) min max: (GLFloat) max
{
	GLClipOperation *operation = [[GLClipOperation alloc] initWithVectorOperand: self firstScalarOperand:min secondScalarOperand:max];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result;
}

//- (id) randomMin: (GLFloat) min max: (GLFloat) max
//{
//	GLRandomNumberOperation *operation = [[GLRandomNumberOperation alloc] initWithVectorOperand: self firstScalarOperand:min secondScalarOperand:max];
//    operation = [self replaceWithExistingOperation: operation];
//	return operation.result;
//}

- (GLVariable *) pow: (GLFloat) aScalar
{
	GLPowerOperation *operation = [[GLPowerOperation alloc] initWithVectorOperand: self scalarOperand: aScalar];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result;
}

- (GLVariable *) max
{
	GLMaxOperation *operation = [[GLMaxOperation alloc] initWithOperand: self];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result;
}

- (GLVariable *) mean: (NSUInteger) index
{
	GLAverageOperation *operation = [[GLAverageOperation alloc] initWithOperand: self dimensionIndex: index];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result;
}

- (id) setValue: (GLFloat) aScalar atIndices: (NSString *) string
{
    GLSetValueOperation *operation = [[GLSetValueOperation alloc] initWithVectorOperand: self scalarOperand:aScalar indexString: string];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result;
}

- (id) setVariableValue: (GLVariable *) aScalarVariable atIndices: (NSString *) string
{
    GLSetVariableValueOperation *operation = [[GLSetVariableValueOperation alloc] initWithVectorOperand: self scalarVariableOperand: aScalarVariable indexString: string];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result;
}

- (GLVariable *) transformToBasis: (NSArray *) orderedBasis
{
#if MatrixFFT
    GLBasisTransformOperation *operation = [GLMatrixFFTTransformOperation basisTransformationWithOperand: self destinationBasis: orderedBasis];
#else
	GLBasisTransformOperation *operation = [GLBasisTransformOperation basisTransformationWithOperand: self destinationBasis: orderedBasis];
#endif
	if (!operation) {
		return self;
	}
//    operation = [self replaceWithExistingOperation: operation];
	return operation.result;
}

- (GLVariable *) fourierTransform
{
    GLVariable *result;
    if (self.isFrequencyDomain) {
        result = [self transformToBasis: @[@(kGLDeltaBasis)]];
    } else {
        NSArray *array = self.differentiationBasis ?: [self.equation defaultDifferentiationBasisForOrder: self.dimensions.count];
        result = [self transformToBasis: array];
    }		
	return result;
}

- (GLVariable *) frequencyDomain
{
	if (self.isFrequencyDomain) {
		return self;
	} else {
		return [self fourierTransform];
	}
}

- (GLVariable *) spatialDomain
{
	if (!self.isFrequencyDomain) {
		return self;
	} else {
		return [self fourierTransform];
	}
}

- (GLVariable *) swapComplex
{
	GLSwapComplexOperation *operation = [[GLSwapComplexOperation alloc] initWithOperand: self];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result;
}

- (GLVariable *) duplicate
{
	GLCopyVariableOperation *operation = [[GLCopyVariableOperation alloc] initWithOperand: self];
    return operation.result;
}

- (GLVariable *) interpolateAtPoints: (NSArray *) otherVariables
{
	GLSimpleInterpolationOperation *operation = [[GLSimpleInterpolationOperation alloc] initWithFirstOperand: [NSArray arrayWithObject: self] secondOperand: otherVariables];
    operation = [self replaceWithExistingOperation: operation];
	return operation.result[0];
}

- (id) scaleVariableBy: (GLFloat) varScale withUnits: (NSString *) varUnits dimensionsBy: (GLFloat) dimScale units: (NSString *) dimUnits
{
	GLScaleOperation *operation = [[GLScaleOperation alloc] initWithOperand: self variableScale: varScale units: varUnits dimensionScale: dimScale translation: 0.0 units: dimUnits];
	operation = [self replaceWithExistingOperation: operation];
	return operation.result;
}

- (id) projectOntoDimensions: (NSArray *) dims usingSpectralBasis: (NSArray *) basis
{
    GLZeroPadOperation *operation = [[GLZeroPadOperation alloc] initWithOperand: self newDimensions: dims basis: basis];
    return operation.result;
}

/************************************************/
/*		Dimension Gymnastics					*/
/************************************************/

#pragma mark -
#pragma mark Dimension Gymnastics
#pragma mark

// This prepends the newDimension to the dimension list. The dimension should have 1 point, if it has more, the values will simply be copied.
- (GLVariable *) variableByAddingDimension: (GLDimension *) newDimension
{
	if ([self.dimensions containsObject: newDimension]) {
		return self;
	} else {
		GLAddDimensionOperation *operation = [[GLAddDimensionOperation alloc] initWithOperand: self dimension: newDimension];
		return operation.result;
	}
}

- (GLVariable *) convertToSplitComplex
{
	GLHalfToSplitComplexOperation *operation = [[GLHalfToSplitComplexOperation alloc] initWithOperand: self];
	return operation.result;
}

- (GLVariable *) variableFromIndexRangeString: (NSString *) indexString
{
    NSArray *ranges = [GLDimension rangesFromIndexString: indexString usingDimensions: self.dimensions];
	GLSubdomainOperation *operation = [[GLSubdomainOperation alloc] initWithOperand: self indexRange: ranges flatten: YES];
	return operation.result;
}

- (GLVariable *) variableFromIndexRange: (NSArray *) ranges
{
	GLSubdomainOperation *operation = [[GLSubdomainOperation alloc] initWithOperand: self indexRange: ranges flatten: YES];
	return operation.result;
}

- (GLVariable *) variableByConcatenatingWithVariable: (GLVariable *) otherVariable alongExistingDimension: (GLDimension *) aDim
{
	NSUInteger index = [self.dimensions indexOfObject: aDim];
	GLExistingDimensionConcatenationOperation *operation = [[GLExistingDimensionConcatenationOperation alloc] initWithFirstOperand: self secondOperand: otherVariable dimensionIndex:index];
	return operation.result;
}

- (GLVariable *) variableByConcatenatingWithVariable: (GLVariable *) otherVariable alongNewDimension: (GLDimension *) aDim
{
	GLNewDimensionConcatenationOperation *operation = [[GLNewDimensionConcatenationOperation alloc] initWithFirstOperand: self secondOperand: otherVariable dimension: aDim];
	return operation.result;
}

- (GLVariable *) variableByConcatenatingWithVariable: (GLVariable *) otherVariable alongDimension: (GLDimension *) aDim
{
	if ([self.dimensions containsObject: aDim]) {
		return [self variableByConcatenatingWithVariable: otherVariable alongExistingDimension:aDim];
	} else {
		return [self variableByConcatenatingWithVariable: otherVariable alongNewDimension:aDim];
	}
}

/************************************************/
/*		Differential Operations	- Primitive		*/
/************************************************/

#pragma mark -
#pragma mark Differential Operations - Primitive
#pragma mark

- (GLVariable *) differentiateWithOperator: (GLDifferentialOperator *) diffOperator
{
    GLVariable *transformedVariable = [self transformToBasis: diffOperator.basis];
    
    GLBinaryOperation *operation = [diffOperator differentiationOperationFromVariable: transformedVariable];
    operation = [self replaceWithExistingOperation: operation];

	return operation.result;
}

/************************************************/
/*		Differential Operations					*/
/************************************************/

#pragma mark -
#pragma mark Differential Operations
#pragma mark

@synthesize differentialOperatorPool=_differentialOperatorPool;

- (id) differentialOperatorPool
{
	if(!_differentialOperatorPool) {
		return [self.equation defaultDifferentialOperatorPoolForVariable:self];
	} else {
		return _differentialOperatorPool;
	}
}

- (GLVariable *) diff: (NSString *) operatorName
{
	GLDifferentialOperator *diffOperator = [self.differentialOperatorPool differentialOperatorWithName: operatorName];
	if (!diffOperator) {
		[NSException raise: @"DifferentialOperatorDoesNotExist" format: @"The differential operator %@ does not exist.", operatorName];
	}
	
	return [self differentiateWithOperator: diffOperator];
}

- (GLVariable *) differentiate: (NSString *) operatorName byTransformingToBasis: (NSArray *) orderedBasis
{
    GLVariable *transformedVariable = [self transformToBasis: orderedBasis];
    transformedVariable.differentiationBasis = orderedBasis;
    return [transformedVariable diff: operatorName];
}

- (NSMethodSignature *)methodSignatureForSelector:(SEL)selector
{
	NSMethodSignature *signature = [super methodSignatureForSelector:selector];
	if (!signature)
	{
		NSString *possibleDiff = NSStringFromSelector(selector);
		GLDifferentialOperator * diffOperator = [self.differentialOperatorPool differentialOperatorWithName: possibleDiff];
		
		if (diffOperator)
		{
			signature = [(NSObject *) self methodSignatureForSelector:@selector(differentiateWithOperator:)];
		}
	}
	
	return signature;
}

- (void)forwardInvocation:(NSInvocation *)anInvocation
{
	NSString *possibleDiff = NSStringFromSelector([anInvocation selector]);
	GLDifferentialOperator * diffOperator = [self.differentialOperatorPool differentialOperatorWithName: possibleDiff];
	
	if (diffOperator)
	{		
		[anInvocation setSelector: @selector(differentiateWithOperator:)];
		[anInvocation setArgument: &diffOperator atIndex:2];
		[anInvocation invokeWithTarget: self];
	}
	else
	{
		[super forwardInvocation:anInvocation];
	}
}

/************************************************/
/*		Reading & Writing						*/
/************************************************/

#pragma mark -
#pragma mark Reading & Writing
#pragma mark

// These methods create a new file and write out the variable.
- (BOOL) writeToNetCDFFile: (NSURL *) anURL
{
	GLNetCDFFile *netcdf = [[GLNetCDFFile alloc] initWithURL: anURL forEquation: self.equation overwriteExisting:YES];
	
	[self solve];
	[netcdf addVariable: self];
	[netcdf waitUntilAllOperationsAreFinished];
	[netcdf close];
	
	return YES;
}

- (NSString *) matrixDescriptionString
{
	NSUInteger n = [[self.dimensions lastObject] nPoints];
	n = n==0?1:n;
	NSMutableString *descrip = [NSMutableString string];
	
	GLFloat max, min;
	vGL_maxv( self.data.mutableBytes, 1, &max, self.nDataElements);
	vGL_minv( self.data.mutableBytes, 1, &min, self.nDataElements);
	
	if ( fabs(min) > max) {
		max = fabs(min);
	}
	
	GLFloat divisor = pow(10, floor(log10(max)));
	if ( divisor == 0.0) divisor = 1;
	
	if (0 && self.dataFormat == kGLSplitComplexDataFormat)
	{
		GLSplitComplex splitComplex = self.splitComplex;
		[descrip appendFormat: @"%f * ", divisor];
//		for (NSUInteger i=0; i<self.nDataPoints; i++)
//		{
//			if ( i % n == 0 ) {
//				[descrip appendFormat: @"\n"];
//			}
//			[descrip appendFormat: @"%1.1f ", sqrt(fabs(splitComplex.realp[i] * splitComplex.realp[i] - splitComplex.imagp[i] * splitComplex.imagp[i]))/divisor];
//		}
		
		for (NSUInteger i=0; i<self.nDataPoints; i++)
		{
			if ( i % n == 0 ) {
				[descrip appendFormat: @"\n"];
			}
			[descrip appendFormat: @"%+1.1f ", splitComplex.realp[i]/divisor];
		}
		
		[descrip appendFormat: @" imagp \n"];
		for (NSUInteger i=0; i<self.nDataPoints; i++)
		{
			if ( i % n == 0 ) {
				[descrip appendFormat: @"\n"];
			}
			[descrip appendFormat: @"%+1.1f ", splitComplex.imagp[i]/divisor];
		}
	}
    if ( self.dimensions.count == 3)
    {
        NSUInteger m = [self.dimensions[2] nPoints] * [self.dimensions[1] nPoints];
        GLFloat *f = self.pointerValue;
		[descrip appendFormat: @"%g * ", divisor];
		for (NSUInteger i=0; i<self.nDataElements; i++)
		{
			if ( i % m == 0 ) {
				[descrip appendFormat: @"\n"];
			}
            if ( i % n == 0 ) {
				[descrip appendFormat: @"\n"];
			}
            
			[descrip appendFormat: @"%+1.1f ", f[i]/divisor];
		}
    }
	else
	{
		GLFloat *f = self.pointerValue;
		[descrip appendFormat: @"%g * ", divisor];
		for (NSUInteger i=0; i<self.nDataElements; i++)
		{
			if ( i % n == 0 ) {
				[descrip appendFormat: @"\n"];
			}
			[descrip appendFormat: @"%+1.1f ", f[i]/divisor];
		}
	}
	
	return descrip;
}

- (void) dumpToConsole
{
	NSLog(@"%@", [self matrixDescriptionString]);
}

- (NSString *) graphvisDescription
{
    NSMutableString *extra = [NSMutableString stringWithFormat: @""];
    if (self.name) [extra appendFormat: @"%@:", self.name];
	[extra appendString: self.isComplex ? @"complex" : @"real"];
	if (self.isRealPartZero) [extra appendString:@", zero real part"];
    if (self.isComplex && self.isImaginaryPartZero) [extra appendString:@", zero imaginary part"];
	if (self.isHermitian) [extra appendString: @"hermitian"];
    for (GLDimension *dim in self.dimensions) {
        [extra appendFormat: @"\\n%@", dim.graphvisDescription];
    }
    
    return extra;
}

- (NSString *) description
{
//	return [NSString stringWithFormat: @"%@ <0x%lx> (%@: %lu points)", NSStringFromClass([self class]), (NSUInteger) self, self.name, self.nDataPoints];
	NSMutableString *extra = [NSMutableString stringWithFormat: @""];
	[extra appendString: self.isComplex ? @"complex variable with" : @"real variable with"];
	[extra appendString: self.isRealPartZero ? @" zero real part" : @" nonzero real part"];
	[extra appendString: self.isImaginaryPartZero ? @", zero imaginary part" : @", nonzero imaginary part"];
	[extra appendString: self.isHermitian ? @" and has hermitian symmetry." : @"."];
	
    //return [NSString stringWithFormat: @"%@ <0x%lx> (%@: %lu points) %@\n", NSStringFromClass([self class]), (NSUInteger) self, self.name, self.nDataPoints, extra];
    
	return [NSString stringWithFormat: @"%@ <0x%lx> (%@: %lu points) %@\n%@", NSStringFromClass([self class]), (NSUInteger) self, self.name, self.nDataPoints, extra, [self matrixDescriptionString]];
}

@end

@implementation GLMutableVariable

- (BOOL) isMutable {
	return YES;
}

- (void) concatenateWithVariable: (GLVariable *) otherVariable alongDimensionAtIndex: (NSUInteger) mutableDimensionIndex
{
	NSLog(@"Method not yet implemented: - (void) concatenateWithVariable: (GLVariable *) otherVariable alongDimensionAtIndex: (NSUInteger) mutableDimensionIndex");
}

- (void) concatenateWithLowerDimensionalVariable: (GLVariable *) otherVariable alongDimensionAtIndex: (NSUInteger) mutableDimensionIndex toIndex: (NSUInteger) pointIndex
{
	NSLog(@"Method not yet implemented:- (void) concatenateWithLowerDimensionalVariable: (GLVariable *) otherVariable alongDimensionAtIndex: (NSUInteger) mutableDimensionIndex toIndex: (NSUInteger) pointIndex;");
}

@end

