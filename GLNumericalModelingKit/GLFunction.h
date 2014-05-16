//
//  GLVariable.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 4/21/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLVariable.h>
#import <GLNumericalModelingKit/GLDimension.h>
#import <GLNumericalModelingKit/GLScalar.h>

// If a symmetry exists, then function is assumed to be defined for negative
// values with the given symmetry. For example, if the dimension has even symmetry,
// then any function defined on this dimension is assumed to have values f(-x), x<0.
typedef NS_ENUM(NSUInteger, GLVariableSymmetry) {
	kGLNoSymmetry = 0,          // No known symmetry
    kGLZeroSymmetry = 1,        // The value is zero (and therefore both even and odd symmetric)
	kGLEvenSymmetry = 2,        // Even symmetry, H(-f)=H(f)
	kGLOddSymmetry = 3          // Odd symmetry, H(-f)=-H(f)
};

@class GLEquation, GLVariableOperation, GLDifferentialOperator;
@interface GLFunction : GLVariable

/************************************************/
/*		Initialization							*/
/************************************************/

#pragma mark -
#pragma mark Initialization
#pragma mark

// Returns a variable built from the given dimension, either a GLVariable or a GLMutableVariable.
+ (id) functionOfRealTypeFromDimension: (GLDimension *) aDim withDimensions: (NSArray *) theDimensions forEquation: (GLEquation *) equation;
+ (id) functionOfComplexTypeFromDimension: (GLDimension *) aDimension withDimensions: (NSArray *) theDimensions forEquation: (GLEquation *) equation;

+ (id) functionWithRandomValuesBetween: (GLFloat) min and: (GLFloat) max withDimensions: (NSArray *) theDimensions forEquation: (GLEquation *) equation;
+ (id) functionWithNormallyDistributedValueWithDimensions: (NSArray *) theDimensions forEquation: (GLEquation *) equation;

// Returns an empty variable (no value) from the given dimensions.
+ (id) functionOfRealTypeWithDimensions: (NSArray *) theDimensions forEquation: (GLEquation *) equation;
+ (id) functionOfComplexTypeWithDimensions: (NSArray *) theDimensions forEquation: (GLEquation *) equation;
+ (id) functionOfType: (GLDataFormat) dataFormat withDimensions: (NSArray *) theDimensions forEquation: (GLEquation *) equation;

// Copies the data from the other variable (now!) not a delayed operation.
+ (id) functionFromFunction: (GLFunction *) otherVariable;


/************************************************/
/*		Dimensionality							*/
/************************************************/

#pragma mark -
#pragma mark Dimensionality
#pragma mark

// The most fundamental requirement of a variable is that is has some dimensions.
// This is an array of GLDimensionObjects.
@property(readonly, strong, nonatomic) NSArray *dimensions;

// Derived property that returns YES if any one of the dimensions is in the frequency domain.
@property(readonly, assign, nonatomic) BOOL isFrequencyDomain;

// This symmetry stuff isn't used and might be crap.

// The symmetry (none, even, or odd) of each dimension.
@property(readwrite, strong, nonatomic) NSMutableArray *realSymmetry;
@property(readwrite, strong, nonatomic) NSMutableArray *imaginarySymmetry;

// Returns YES if the variable is Hermitian,  H(-f)=H†(f).
@property(readonly, assign, nonatomic) BOOL isHermitian;

// Return the dimension restricted to positive values that can be used to recover
// the negative values. Returns nil if the variable is not hermitian.
@property(readonly, strong, nonatomic) GLDimension *hermitianDimension;

// Given an array basis functions, returns the associated dimensions.
- (NSArray *) dimensionsTransformedToBasis: (NSArray *) basis;

// Ordered array of basis functions
@property(readonly, strong, nonatomic) NSArray *basis;

/************************************************/
/*		Data									*/
/************************************************/

#pragma mark -
#pragma mark Data
#pragma mark

// Sets the value at each point to a different random number between [-amp, amp]
- (void) rand: (GLFloat) amp;

// Unoptimized check to see if all elements are finite.
- (BOOL) isFinite;

// Returns the max value, NOW.
- (GLFloat) maxNow;
- (GLFloat) minNow;

/************************************************/
/*		Operations								*/
/************************************************/

#pragma mark -
#pragma mark Operations
#pragma mark

/** Divides the function by the second variable: result = receiving variable / otherFunctionOrScalar.
 @discussion This operation is treated differently depending on the otherFunctionOrScalar class.
 @discussion The receiving function can be divided by a scalar (given as a GLScalar or NSNumber) or a function.
 @discussion Note that a different algorithm is used depending on whether the scalar is given as a subclass of NSNumber or GLVariable. In the former case, the value is assumed constant and can't be altered in subsequent uses, while in the latter case it can vary. It is generally most computationally efficient to use the constant value.
 @param otherFunctionOrScalar An input scalar or function.
 @returns A GLFunction object.
 */
- (GLFunction *) dividedBy: (id) otherFunctionOrScalar;

// C = A * B
- (id) dot: (GLFunction *) otherVariable;

// C = max( abs(A), abs(B) ) element-wise
- (id) absMax: (GLFunction *) otherVariable;

/// C = -A
- (GLFunction *) negate;

/// C = abs(A)
- (GLFunction *) abs;

/// C = exp(A)
- (GLFunction *) exponentiate;

/// C = log(A)
- (GLFunction *) log;

/// C = sin(A)
- (GLFunction *) sin;

/// C = cos(A)
- (GLFunction *) cos;

/// C = atan(A)
- (GLFunction *) atan;

/// C = atan2(A, B)
- (GLFunction *) atan2: (GLFunction *) x;

/// C = sinh(A)
- (GLFunction *) sinh;

/// C = asinh(A)
- (GLFunction *) asinh;

/// C = tanh(A)
- (GLFunction *) tanh;

/// C = sqrt(A)
- (GLFunction *) sqrt;

// C = A + k
- (GLFunction *) scalarAdd: (GLFloat) aScalar;

// C = k*A
- (GLFunction *) scalarMultiply: (GLFloat) aScalar;

// C = k/A
- (GLFunction *) scalarDivide: (GLFloat) aScalar;

// C = max( A, k )
- (id) scalarMax: (GLFloat) aScalar;

// C = A > k ? A : 0.0;
- (id) zeroThreshold: (GLFloat) aScalar;

- (id) clipToMin: (GLFloat) min max: (GLFloat) max;

// Set each element to a random number between min and max.
//- (id) randomMin: (GLFloat) min max: (GLFloat) max;

- (id) pow: (GLFloat) aScalar;

// Returns a scalar GLVariable.
- (GLScalar *) max;

- (GLScalar *) min;


// Computes the mean along the specified index, effectively collapsing the variable along that index.
- (id) mean: (NSUInteger) index range: (NSRange) range;
- (id) mean: (NSUInteger) index;
- (GLScalar *) mean;

- (id) sum: (NSUInteger) index;
- (GLScalar *) sum;

- (GLScalar *) integrate;

- (GLFunction *) makeHermitian;

// Matlab style. String should be @"start:end,start:end,..."
- (id) setValue: (GLFloat) aScalar atIndices: (NSString *) string;
- (id) setVariableValue: (GLScalar *) aScalarVariable atIndices: (NSString *) string;

- (id) transformToBasis: (NSArray *) orderedBasis;

- (id) fourierTransform;


// Performs a fourier transform if necessary, otherwise returns self.
- (id) frequencyDomain;
- (id) spatialDomain;

// Returns a variable with the real and imaginary part swapped.
- (GLFunction *) swapComplex;

- (id) duplicate;

- (GLFunction *) interpolateAtPoints: (NSArray *) otherVariables;

// Multiplies the variables and the dimensions by this scale factor.
- (id) scaleVariableBy: (GLFloat) varScale withUnits: (NSString *) varUnits dimensionsBy: (GLFloat) dimScale units: (NSString *) dimUnits;

- (id) projectOntoDimensions: (NSArray *) dims usingSpectralBasis: (NSArray *) basis;

/************************************************/
/*		Dimension Gymnastics					*/
/************************************************/

#pragma mark -
#pragma mark Dimension Gymnastics
#pragma mark

// This prepends the newDimension to the dimension list.
// The dimension should have 1 point, if it has more, the values will simply be copied.
// The returned variable may be mutable if the dimensions dictate.
- (id) variableByAddingDimension: (GLDimension *) newDimension;

// Returns a variable with only the elements indicated by by the array of ranges.
// The size of the ranges array must match the number of dimensions.
- (GLFunction *) variableFromIndexRange: (NSArray *) ranges;
- (GLFunction *) variableFromIndexRangeString: (NSString *) indexString;

//- (GLVariable *) convertToSplitComplex;

// This operation is minimilly restrictive.
// You can concatenate an n dimensional variable with an n or n-1 dimensional variable.
// The dimensionIndex refers to the dimensions of the firstOperand.
// The dimension of concatenation must be evenly sampled.
// The other dimensions must have the same number of points.
// 
// The result variable will have dimensions based on the receiver's dimensions.
- (GLFunction *) variableByConcatenatingWithVariable: (GLFunction *) otherVariable alongExistingDimension: (GLDimension *) aDim;

// The new dimension should have nPoints=2, otherwise a new dimension based on the same properties will be created.
// The dimensions in the other two variables must be the same.
- (GLFunction *) variableByConcatenatingWithVariable: (GLFunction *) otherVariable alongNewDimension: (GLDimension *) aDim;

// Short cut to one of the two operations above, depending on whether or not the dimension exists.
- (GLFunction *) variableByConcatenatingWithVariable: (GLFunction *) otherVariable alongDimension: (GLDimension *) aDim;

/************************************************/
/*		Reading & Writing						*/
/************************************************/

#pragma mark -
#pragma mark Reading & Writing
#pragma mark

// These methods create a new file and write out the variable.
- (BOOL) writeToNetCDFFile: (NSURL *) anURL;

/************************************************/
/*		Differential Operations	- Primitive		*/
/************************************************/

#pragma mark -
#pragma mark Differential Operations - Primitive
#pragma mark

// Apply a differential operation. Transforms to the correct basis, then uses the diffOperator to transform the variable.
- (GLFunction *) differentiateWithOperator: (GLLinearTransform *) diffOperator;

/************************************************/
/*		Differential Operations					*/
/************************************************/

#pragma mark -
#pragma mark Differential Operations
#pragma mark

// All operators are assumed to use the differentiationBasis. So this means that calling "xxx" or any other string will cause,
// 1. This variable to transform to the differentiationBasis, if needed
// 2. A search for a GLLinearTransform with that name and -fromDimensions
// 3. The operator to be applied to this function, and a new function to be returned.

// Occasionally more convenient is to set a "Differential Operator Pool" to be associated with this variable.
// The differential operator pool is a collection of *named* differential operators that can conviently
// be referenced. For example, if you set the SpectralDifferentiationPool, then calling
//	1. [var diff: @"xx"]; or,
//	2. [var xx];
// extracts the differential operator named @"xx" from the pool, then calls -differentiate.

// If it's set to nil (or an empty), it will use the dimension default.
@property(readwrite, strong, nonatomic) NSArray *differentiationBasis;

// Differentiate with one of the built-in (or saved) operators from the pool, e.g. "xxy".
- (GLFunction *) diff: (NSString *) operatorName;

- (GLFunction *) differentiate: (NSString *) operatorName byTransformingToBasis: (NSArray *) orderedBasis;


// The one and only initializer for a variable.
- (id) initVariableOfType: (GLDataFormat) dataFormat withDimensions: (NSArray *) theDimensions forEquation: (GLEquation *) theEquation;

@end


@interface GLMutableVariable : GLFunction

/************************************************/
/*		Dimension Gymnastics					*/
/************************************************/

#pragma mark -
#pragma mark Dimension Gymnastics
#pragma mark

// These operations can only be applied along an existing mutable dimension.
// The result of these operation is that the mutable dimension will increase in length,
// and the variable will add new data.

// Both variables must be n-dimensional.
// The receiver's dimensions must include aDim (at mutableDimensionIndex).
// The dimensions at indices other than mutableDimensionIndex must have the same number of points.
// The dimensions at mutableDimensionIndex do not need to match. If the first operand's dimension is evenly spaced
// then it will be extended, otherwise the values from the second operand's dimension will be added.
- (void) concatenateWithVariable: (GLFunction *) otherVariable alongDimensionAtIndex: (NSUInteger) mutableDimensionIndex;

// The receiver must be n-dimensional and the other variable must be (n-1)-dimensional.
// The (n-1) dimensions of the two variables must have the same number of points.
// If the mutableDimension is evenly spaced, then it will be extended to length pointIndex+1, if necessary.
// If the mutableDimension is not evenly spaced, then it must already have the correct value.
- (void) concatenateWithLowerDimensionalVariable: (GLFunction *) otherVariable alongDimensionAtIndex: (NSUInteger) mutableDimensionIndex toIndex: (NSUInteger) pointIndex;
@end

// These are thrown in here to surpress compiler warnings these (and higher order derivatives)
// will be resolved dynamically at runtime if they exist (given the dimensions).
@interface GLFunction (DifferentiationExtensions)

@property(readonly) GLFunction* x;
@property(readonly) GLFunction* y;
@property(readonly) GLFunction* xx;
@property(readonly) GLFunction* xy;
@property(readonly) GLFunction* yy;
@property(readonly) GLFunction* xxx;
@property(readonly) GLFunction* xxy;
@property(readonly) GLFunction* xyy;
@property(readonly) GLFunction* yyy;

@end

