//
//  GLVariable.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 4/21/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLTensor.h>
#import <GLNumericalModelingKit/GLDimension.h>
#import <GLNumericalModelingKit/GLScalar.h>

// If a symmetry exists, then function is assumed to be defined for negative
// values with the given symmetry. For example, if the dimension has even symmetry,
// then any function defined on this dimension is assumed to have values f(-x), x<0.
enum {
	kGLNoSymmetry = 0,          // No known symmetry
    kGLZeroSymmetry = 1,        // The value is zero (and therefore both even and odd symmetric)
	kGLEvenSymmetry = 2,        // Even symmetry, H(-f)=H(f)
	kGLOddSymmetry = 3          // Odd symmetry, H(-f)=-H(f)
};
typedef NSUInteger GLVariableSymmetry;

@class GLEquation, GLVariableOperation, GLDifferentialOperator;
@interface GLVariable : GLTensor

// Variables do not get computed immediately and should only be computed when absolutely needed or a choke point has been reached.
// This minimizes the amount of memory and computation required. If we computed a variable's value immediately (like [psi x]),
// then we would have to create a strong reference to that variable so that we don't end up recomputing it. As it stands, we create
// weak references to variables.
//
// May need to remove dependencies after computation if it's true that they are strongly referenced.
//
// However, we do *strongly* reference the transformed variable. Does this create a bad situation? Maybe we need operations to only
// weakly reference the instance variable dependencies. Hmm, make sure we set those values to nil.

/************************************************/
/*		Behavior								*/
/************************************************/

#pragma mark -
#pragma mark Behavior
#pragma mark

+ (void) setPrefersSpatialMultiplication: (BOOL) aFlag;
+ (BOOL) prefersSpatialMultiplication;

/************************************************/
/*		Initialization							*/
/************************************************/

#pragma mark -
#pragma mark Initialization
#pragma mark

// Returns a variable built from the given dimension, either a GLVariable or a GLMutableVariable.
+ (id) variableOfRealTypeFromDimension: (GLDimension *) aDim withDimensions: (NSArray *) theDimensions forEquation: (GLEquation *) equation;
+ (id) variableOfComplexTypeFromDimension: (GLDimension *) aDimension withDimensions: (NSArray *) theDimensions forEquation: (GLEquation *) equation;

+ (id) variableWithRandomValuesBetween: (GLFloat) min and: (GLFloat) max withDimensions: (NSArray *) theDimensions forEquation: (GLEquation *) equation;

// Returns an empty variable (no value) from the given dimensions.
+ (id) variableOfRealTypeWithDimensions: (NSArray *) theDimensions forEquation: (GLEquation *) equation;
+ (id) variableOfComplexTypeWithDimensions: (NSArray *) theDimensions forEquation: (GLEquation *) equation;
+ (id) variableOfType: (GLDataFormat) dataFormat withDimensions: (NSArray *) theDimensions forEquation: (GLEquation *) equation;

// Copies the data from the other variable (now!) not a delayed operation.
+ (id) variableFromVariable: (GLVariable *) otherVariable;


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

// The symmetry (none, even, or odd) of each dimension.
@property(readwrite, assign, nonatomic) NSMutableArray *realSymmetry;
@property(readwrite, assign, nonatomic) NSMutableArray *imaginarySymmetry;

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

/************************************************/
/*		Operations								*/
/************************************************/

#pragma mark -
#pragma mark Operations
#pragma mark

// C = A * B
- (id) multiply: (GLVariable *) otherVariable;

- (id) dividedBy: (GLVariable *) otherVariable;

// C = A * B
- (id) dot: (GLVariable *) otherVariable;

// C = max( abs(A), abs(B) ) element-wise
- (id) absMax: (GLVariable *) otherVariable;

// C = -A
- (id) negate;

// C = abs(A)
- (id) abs;

// C = exp(A)
- (id) exponentiate;

// C = sin(A)
- (id) sin;

// C = cos(A)
- (id) cos;

// C = atan(A)
- (id) atan;

// C = sqrt(A)
- (id) sqrt;

// C = A + k
- (id) scalarAdd: (GLFloat) aScalar;

// C = k*A
- (id) scalarMultiply: (GLFloat) aScalar;

// C = k/A
- (id) scalarDivide: (GLFloat) aScalar;

// C = max( A, k )
- (id) scalarMax: (GLFloat) aScalar;

// C = A > k ? A : 0.0;
- (id) zeroThreshold: (GLFloat) aScalar;

- (id) clipToMin: (GLFloat) min max: (GLFloat) max;

// Set each element to a random number between min and max.
//- (id) randomMin: (GLFloat) min max: (GLFloat) max;

- (id) pow: (GLFloat) aScalar;

// Returns a scalar GLVariable.
- (id) max;

// Computes the mean along the specified index, effectively collapsing the variable along that index.
- (id) mean: (NSUInteger) index;

// Matlab style. String should be @"start:end,start:end,..."
- (id) setValue: (GLFloat) aScalar atIndices: (NSString *) string;
- (id) setVariableValue: (GLScalar *) aScalarVariable atIndices: (NSString *) string;

- (id) transformToBasis: (NSArray *) orderedBasis;

- (id) fourierTransform;


// Performs a fourier transform if necessary, otherwise returns self.
- (id) frequencyDomain;
- (id) spatialDomain;

// Returns a variable with the real and imaginary part swapped.
- (id) swapComplex;

- (id) duplicate;

- (GLVariable *) interpolateAtPoints: (NSArray *) otherVariables;

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
- (GLVariable *) variableFromIndexRange: (NSArray *) ranges;
- (GLVariable *) variableFromIndexRangeString: (NSString *) indexString;

//- (GLVariable *) convertToSplitComplex;

// This operation is minimilly restrictive.
// You can concatenate an n dimensional variable with an n or n-1 dimensional variable.
// The dimensionIndex refers to the dimensions of the firstOperand.
// The dimension of concatenation must be evenly sampled.
// The other dimensions must have the same number of points.
// 
// The result variable will have dimensions based on the receiver's dimensions.
- (GLVariable *) variableByConcatenatingWithVariable: (GLVariable *) otherVariable alongExistingDimension: (GLDimension *) aDim;

// The new dimension should have nPoints=2, otherwise a new dimension based on the same properties will be created.
// The dimensions in the other two variables must be the same.
- (GLVariable *) variableByConcatenatingWithVariable: (GLVariable *) otherVariable alongNewDimension: (GLDimension *) aDim;

// Short cut to one of the two operations above, depending on whether or not the dimension exists.
- (GLVariable *) variableByConcatenatingWithVariable: (GLVariable *) otherVariable alongDimension: (GLDimension *) aDim;

/************************************************/
/*		Differential Operations	- Primitive		*/
/************************************************/

#pragma mark -
#pragma mark Differential Operations - Primitive
#pragma mark

// Apply a differential operation. Transforms to the correct basis, then uses the diffOperator to transform the variable.
- (GLVariable *) differentiateWithOperator: (GLLinearTransform *) diffOperator;

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
- (GLVariable *) diff: (NSString *) operatorName;

- (GLVariable *) differentiate: (NSString *) operatorName byTransformingToBasis: (NSArray *) orderedBasis;


// The one and only initializer for a variable.
- (id) initVariableOfType: (GLDataFormat) dataFormat withDimensions: (NSArray *) theDimensions forEquation: (GLEquation *) theEquation;

@end


@interface GLMutableVariable : GLVariable

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
- (void) concatenateWithVariable: (GLVariable *) otherVariable alongDimensionAtIndex: (NSUInteger) mutableDimensionIndex;

// The receiver must be n-dimensional and the other variable must be (n-1)-dimensional.
// The (n-1) dimensions of the two variables must have the same number of points.
// If the mutableDimension is evenly spaced, then it will be extended to length pointIndex+1, if necessary.
// If the mutableDimension is not evenly spaced, then it must already have the correct value.
- (void) concatenateWithLowerDimensionalVariable: (GLVariable *) otherVariable alongDimensionAtIndex: (NSUInteger) mutableDimensionIndex toIndex: (NSUInteger) pointIndex;
@end

// These are thrown in here to surpress compiler warnings these (and higher order derivatives)
// will be resolved dynamically at runtime if they exist (given the dimensions).
@interface GLVariable (DifferentiationExtensions)

@property(readonly) GLVariable* x;
@property(readonly) GLVariable* y;
@property(readonly) GLVariable* xx;
@property(readonly) GLVariable* xy;
@property(readonly) GLVariable* yy;
@property(readonly) GLVariable* xxx;
@property(readonly) GLVariable* xxy;
@property(readonly) GLVariable* xyy;
@property(readonly) GLVariable* yyy;

@end

