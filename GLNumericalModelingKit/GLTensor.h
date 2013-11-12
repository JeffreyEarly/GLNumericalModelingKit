//
//  GLTensor.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 10/11/13.
//  Copyright (c) 2013 Jeffrey J. Early. All rights reserved.
//

#import <GLNumericalModelingKit/Precision.h>
#import <GLNumericalModelingKit/GLMatrixDescription.h>

// This specifies how the data is organized in the memory buffer.
// kGLRealDataFormat means that there is no memory allocated for the imaginary part.
// kGLSplitComplexDataFormat means that the imaginary nPoints follow the real nPoints in the buffer.
// kGLInterleavedComplexFormat means that the imaginary part of a point immediately follows the real part.
typedef NS_ENUM(NSUInteger, GLDataFormat) {
	kGLRealDataFormat = 0,
    kGLSplitComplexDataFormat = 1,
    kGLInterleavedComplexDataFormat = 2
};

GLSplitComplex splitComplexFromData( NSData *data );

@class GLVariableOperation, GLEquation;
@interface GLTensor : NSObject

+ (GLTensor *) variableWithPrototype: (GLTensor *) anotherVariable;

// Not to be called directly. Only to be used by the subclasses.
- (id) initWithType: (GLDataFormat) dataFormat withEquation: (GLEquation *) theEquation;

/************************************************/
/*		Properties								*/
/************************************************/

#pragma mark -
#pragma mark Properties
#pragma mark

// Naming a tensor or giving it units is useful when writing output to files, and may even be mandatory, depending on the format.
@property(readwrite, copy, nonatomic) NSString *name;
@property(readwrite, copy, nonatomic) NSString *units;

// Any metadata that should follow around the variable. Units property is automicatically added to this.
@property(readonly, strong, nonatomic) NSMutableDictionary *metadata;

// An attempt to make a fairly unique variable id. Copies of this variable have the same id.
@property(readonly, assign, nonatomic) NSUInteger uniqueID;

// Rank of the tensor. 0-scalar, 1-vector, 2-linear transformation
@property(readonly, assign, nonatomic) NSUInteger rank;

// Determines whether the data is holding a real or complex number.
// Variables in the frequency domain are always assumed to be complex.
@property(readonly, assign, nonatomic) BOOL isComplex;

@property(readwrite, assign, nonatomic) BOOL isRealPartZero;
@property(readwrite, assign, nonatomic) BOOL isImaginaryPartZero;

@property(readwrite, assign, nonatomic) BOOL isPurelyReal;
@property(readwrite, assign, nonatomic) BOOL isPurelyImaginary;

// Returns NO.
@property(readonly, assign, nonatomic) BOOL isMutable;

/************************************************/
/*		Data									*/
/************************************************/

#pragma mark -
#pragma mark Data
#pragma mark

/// The total number of components saved to the data. Note that some formats may not save zeros, for example, while others may save extra points than required.
@property(readonly, assign, nonatomic) NSUInteger nDataPoints;

/// For a real number, the number of data points is equal to the number of elements. For a split complex number, there are twice as many elements as points.
@property(readonly, assign, nonatomic) NSUInteger nDataElements;

/// Access to the raw computed data. If this variable is dependent on others, you should call -solve first, otherwise an empty (non-zeroed!) chunk of data will be returned.
@property(readonly, strong, nonatomic) NSMutableData *data;

/// Number of bytes allocated in the raw data, essentially the nDataElements*sizeof(GLFloat)
@property(readonly, assign, nonatomic) NSUInteger dataBytes;

/// Indicates whether or not an NSMutableData object has been assigned yet.
@property(readonly, assign, nonatomic) BOOL hasData;

// The data format for each dimension corresponds 1-1 with the dataFormats array.
// If the data formats are homogenous, -dataFormat will return the value, otherwise
// it will return kGLMixedDataFormat.
@property(readonly, assign, nonatomic) GLDataFormat dataFormat;

@property(readwrite, strong) GLMatrixDescription *matrixDescription;

// This will return a (GLSplitComplex *) pointing to the data the variable is complex,
// or it will return a (GLFloat *) pointing to the data otherwise. Request the right one!
@property(readonly, assign, nonatomic) GLFloat *pointerValue;
@property(readonly, assign, nonatomic) GLSplitComplex splitComplex;

- (void) solve;

// Set the value to zero everywhere.
- (void) zero;

/************************************************/
/*		Operations								*/
/************************************************/

#pragma mark -
#pragma mark Operations
#pragma mark

/** Adds two variables together: result = receiving variable + otherVariableOrScalar.
 @discussion Adding two artibrary variables together doesn't always make sense and may cause an exception to be thrown. For example, adding a function to a matrix is undefined, and adding two matrices or two functions with different dimensions is also undefined. Adding a scalar always makes sense.
 @discussion Note that a different algorithm is used depending on whether the scalar is given as a subclass of NSNumber or GLVariable. In the former case, the value is assumed constant and can't be altered in subsequent uses, while in the latter case it can vary. It is generally most computationally efficient to use the constant value.
 @param otherVariableOrScalar An input scalar, vector or matrix.
 @returns A subclass of GLVariable with the same dimensions as the highest rank input.
 */
- (instancetype) plus: (id) otherVariableOrScalar;

/** Subtract two variables: result = receiving variable - otherVariableOrScalar.
 @discussion Subtracting two artibrary variables together doesn't always make sense and may cause an exception to be thrown. For example, subtracting a function to a matrix is undefined, and subtracting two matrices or two functions with different dimensions is also undefined. Subtracting a scalar always makes sense.
 @discussion Note that a different algorithm is used depending on whether the scalar is given as a subclass of NSNumber or GLVariable. In the former case, the value is assumed constant and can't be altered in subsequent uses, while in the latter case it can vary. It is generally most computationally efficient to use the constant value.
 @param otherVariableOrScalar An input scalar, vector or matrix.
 @returns A subclass of GLVariable with the same dimensions as the highest rank input.
 */
- (instancetype) minus: (id) otherVariableOrScalar;

/** Multiplies two variables together: result = receiving variable * otherVariableOrScalar.
 @discussion This operation is treated differently depending on the receiver and otherVariable class.
 @discussion A scalar can multiply another scalar, a function, or a linear transformation.
 @discussion Two functions can be multiplied together if they have the same dimension, e.g., h(x)=f(x)*g(x), or if one dimension is missing, e.g. h(x,y) = f(x)*g(x,y).
 @discussion A linear transformation 'times' a function is treated as a transformation of that function, and is only valid if the fromDimensions of the linear transformation match the dimensions of the function. Left-sided transformations are not supported.
 @discussion Two linear transformations are treated as a matrix multiplication and the fromDimensions of the receiver must match the the toDimensions of the otherVariable linear transformation.
 @discussion Note that a different algorithm is used depending on whether the scalar is given as a subclass of NSNumber or GLVariable. In the former case, the value is assumed constant and can't be altered in subsequent uses, while in the latter case it can vary. It is generally most computationally efficient to use the constant value.
 @param otherVariableOrScalar An input scalar, vector or matrix.
 @returns A subclass of GLVariable.
 */
- (id) multiply: (id) otherVariableOrScalar;

/** Multiplies two variables together: result = receiving variable * otherVariableOrScalar, after transforming back into the spatial domain.
 @discussion This is the same as multiply, but transforms variables into the spatial domain (if necessary) before multiplying.
  @param otherVariableOrScalar An input scalar, vector or matrix.
 @returns A subclass of GLVariable.
 */
- (id) times: (id) otherVariableOrScalar;

/************************************************/
/*		Reading & Writing						*/
/************************************************/

#pragma mark -
#pragma mark Reading & Writing
#pragma mark

// These methods create a new file and write out the variable.
- (BOOL) writeToNetCDFFile: (NSURL *) anURL;

// Description of the tensor appropriate for output in an object graph.
@property(readonly) NSString *graphvisDescription;

// Essentially just print the long-form description of this variable to the console.
- (void) dumpToConsole;

/************************************************/
/*		Private									*/
/************************************************/

#pragma mark -
#pragma mark Private
#pragma mark

// If I try to use an operation a second time, I'm going to find it won't work.

@property(readonly, weak, nonatomic) GLEquation *equation;

// Operations upon which this variable depends.
@property(readonly, strong, nonatomic) NSMutableArray *pendingOperations;

// Operations upon which this variable depends.
- (void) addOperation: (id) operation;
- (void) removeOperation: (id) operation;

// The last operation upon which this variable is dependent.
- (GLVariableOperation *) lastOperation;

// Operations for which this variable is an operand
@property(readonly, strong) NSMutableArray *existingOperations;

// This method takes an operation and checks the GLVariable object's internal cache
// to see if the equivalent operation has already been computed.
- (id) replaceWithExistingOperation: (GLVariableOperation *) newOperation;

@end
