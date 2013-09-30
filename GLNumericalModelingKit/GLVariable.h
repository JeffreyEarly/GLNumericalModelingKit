//
//  GLVariable.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 4/21/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>

#import <GLNumericalModelingKit/Precision.h>
#import <GLNumericalModelingKit/GLDimension.h>

GLSplitComplex splitComplexFromData( NSData *data );

// This specifies how the data is organized in the memory buffer.
// kGLRealDataFormat means that there is no memory allocated for the imaginary part.
// kGLSplitComplexDataFormat means that the imaginary nPoints follow the real nPoints in the buffer.
// The half complex format means that the reversed imaginary part follows real part *in that dimension*.
// Unlike the split complex format then, you must specify which dimension is responsible for imaginary parts.
enum {
	kGLRealDataFormat = 0,
    kGLSplitComplexDataFormat = 1,
    kGLInterleavedComplexDataFormat = 2,
	kGLHalfComplex1DDataFormat = 10,
    kGLHalfComplex2DDataFormat = 11,
    kGLHalfComplex3DDataFormat = 12
};
typedef NSUInteger GLDataFormat;

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

@class GLEquation, GLVariableOperation, GLDifferentialOperator, GLMatrixDescription;
@interface GLVariable : NSObject {
    NSArray *_dimensions;
	NSUInteger _nDataPoints;
	NSUInteger _nDataElements;
	BOOL _isFrequencyDomain;
	BOOL _isComplex;
	BOOL _isImaginaryPartZero;
	NSMutableDictionary *_metadata;
	NSUInteger _uniqueID;
	
	NSMutableData *_data;
	NSUInteger _dataBytes;
	GLDataFormat _dataFormat;
	
    __weak GLEquation *_equation;
	NSMutableArray *_existingOperations;
	NSMutableArray *_pendingOperations;
	GLVariable *_transformedVariable;
	id _differentialOperatorPool;
	NSMutableSet *_variableDifferentialOperationMaps;
}

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

@property(readwrite, copy, nonatomic) NSString *name;
@property(readwrite, copy, nonatomic) NSString *units;

// Any metadata that should follow around the variable. Units property is automicatically added to this.
@property(readonly, strong, nonatomic) NSMutableDictionary *metadata;

// An attempt to make a fairly unique variable id. Copies of this variable have the same id.
@property(readonly, assign, nonatomic) NSUInteger uniqueID;

@property(readonly) NSString *graphvisDescription;
@property(readwrite, strong) GLMatrixDescription *matrixDescription;

/************************************************/
/*		Dimensionality							*/
/************************************************/

#pragma mark -
#pragma mark Dimensionality
#pragma mark

// The most fundamental requirement of a variable is that is has some dimensions.
// This is an array of GLDimensionObjects.
@property(readonly, strong, nonatomic) NSArray *dimensions;

// Derived by computing the product of the number of data points in each dimension.
@property(readonly, assign, nonatomic) NSUInteger nDataPoints;

// For a real number, the number of data points is equal to the number of elements.
// For a split complex number, there are twice as many elements as points.
@property(readonly, assign, nonatomic) NSUInteger nDataElements;

// Derived property that returns YES if any one of the dimensions is in the frequency domain.
@property(readonly, assign, nonatomic) BOOL isFrequencyDomain;

// Determines whether the data is holding a real or complex number.
// Variables in the frequency domain are always assumed to be complex.
@property(readonly, assign, nonatomic) BOOL isComplex;

// Returns NO.
@property(readonly, assign, nonatomic) BOOL isMutable;

// The symmetry (none, even, or odd) of each dimension.
@property(readwrite, assign, nonatomic) NSMutableArray *realSymmetry;
@property(readwrite, assign, nonatomic) NSMutableArray *imaginarySymmetry;

@property(readwrite, assign, nonatomic) BOOL isRealPartZero;
@property(readwrite, assign, nonatomic) BOOL isImaginaryPartZero;

@property(readwrite, assign, nonatomic) BOOL isPurelyReal;
@property(readwrite, assign, nonatomic) BOOL isPurelyImaginary;

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

// Access to the raw computed data. If this variable is dependent on others, you should call
// have the equation -solveForVariable first, otherwise an empty (non-zeroed!) chunk of data will be returned.
@property(readonly, strong, nonatomic) NSMutableData *data;
@property(readonly, assign, nonatomic) NSUInteger dataBytes;
@property(readonly, assign, nonatomic) BOOL hasData;

// The data format for each dimension corresponds 1-1 with the dataFormats array.
// If the data formats are homogenous, -dataFormat will return the value, otherwise
// it will return kGLMixedDataFormat.
@property(readonly, assign, nonatomic) GLDataFormat dataFormat;

// This will return a (GLSplitComplex *) pointing to the data the variable is complex,
// or it will return a (GLFloat *) pointing to the data otherwise. Request the right one!
@property(readonly, assign, nonatomic) GLFloat *pointerValue;
@property(readonly, assign, nonatomic) GLSplitComplex splitComplex;

- (void) solve;

// Set the value to zero everywhere.
- (void) zero;

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

// This method takes an operation and checks the GLVariable object's internal cache
// to see if the equivalent operation has already been computed.
- (id) replaceWithExistingOperation: (GLVariableOperation *) newOperation;

// C = A + B
- (id) plus: (GLVariable *) otherVariable NS_RETURNS_NOT_RETAINED;

// C = A - B
- (id) minus: (GLVariable *) otherVariable NS_RETURNS_NOT_RETAINED;

// C = A * B (will switch to spatial domain, if default)
- (id) times: (GLVariable *) otherVariable NS_RETURNS_NOT_RETAINED;

// C = A * B
- (id) multiply: (GLVariable *) otherVariable NS_RETURNS_NOT_RETAINED;

- (id) dividedBy: (GLVariable *) otherVariable;

// C = A * B
- (id) dot: (GLVariable *) otherVariable NS_RETURNS_NOT_RETAINED;

// C = max( abs(A), abs(B) ) element-wise
- (id) absMax: (GLVariable *) otherVariable NS_RETURNS_NOT_RETAINED;

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
- (id) scalarAdd: (GLFloat) aScalar NS_RETURNS_NOT_RETAINED;

// C = k*A
- (id) scalarMultiply: (GLFloat) aScalar NS_RETURNS_NOT_RETAINED;

// C = k/A
- (id) scalarDivide: (GLFloat) aScalar NS_RETURNS_NOT_RETAINED;

// C = max( A, k )
- (id) scalarMax: (GLFloat) aScalar NS_RETURNS_NOT_RETAINED;

// C = A > k ? A : 0.0;
- (id) zeroThreshold: (GLFloat) aScalar NS_RETURNS_NOT_RETAINED;

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
- (id) setVariableValue: (GLVariable *) aScalarVariable atIndices: (NSString *) string;

- (id) transformToBasis: (NSArray *) orderedBasis;

- (id) fourierTransform;


// Performs a fourier transform if necessary, otherwise returns self.
- (id) frequencyDomain;
- (id) spatialDomain;

// Returns a variable with the real and imaginary part swapped.
- (id) swapComplex;

- (id) duplicate;

- (GLVariable *) interpolateAtPoints: (NSArray *) otherVariables NS_RETURNS_NOT_RETAINED;

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

- (GLVariable *) convertToSplitComplex;

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

// Apply a differential operation.
- (GLVariable *) differentiateWithOperator: (GLDifferentialOperator *) diffOperator;

/************************************************/
/*		Differential Operations					*/
/************************************************/

#pragma mark -
#pragma mark Differential Operations
#pragma mark

// Occasionally more convenient is to set a "Differential Operator Pool" to be associated with this variable.
// The differential operator pool is a collection of *named* differential operators that can conviently
// be referenced. For example, if you set the SpectralDifferentiationPool, then calling
//	1. [var diff: @"xx"]; or,
//	2. [var xx];
// extracts the differential operator named @"xx" from the pool, then calls -differentiate.

@property(readwrite, strong, nonatomic) id differentialOperatorPool;

// If it's set to nil (or an empty), it will use the equation default.
@property(readwrite, strong, nonatomic) NSArray *differentiationBasis;

// Differentiate with one of the built-in (or saved) operators from the pool, e.g. "xxy".
- (GLVariable *) diff: (NSString *) operatorName;

- (GLVariable *) differentiate: (NSString *) operatorName byTransformingToBasis: (NSArray *) orderedBasis;

/************************************************/
/*		Reading & Writing						*/
/************************************************/

#pragma mark -
#pragma mark Reading & Writing
#pragma mark

// These methods create a new file and write out the variable.
- (BOOL) writeToNetCDFFile: (NSURL *) anURL;
- (void) dumpToConsole;

/************************************************/
/*		Private									*/
/************************************************/

#pragma mark -
#pragma mark Private
#pragma mark

// If I try to use an operation a second time, I'm going to find it won't work.

@property(readonly, weak, nonatomic) GLEquation *equation;
@property(readonly, strong, nonatomic) NSMutableArray *pendingOperations;

// Operations upon which this variable depends.
- (void) addOperation: (id) operation;
- (void) removeOperation: (id) operation;

// The last operation upon which this variable is dependent.
- (GLVariableOperation *) lastOperation;

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

