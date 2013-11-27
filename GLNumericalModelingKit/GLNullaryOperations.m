//
//  GLNullaryOperations.m
//  GLNumericalModelingKitCopy
//
//  Created by Jeffrey J. Early on 11/27/12.
//
//

#import "GLNullaryOperations.h"

/************************************************/
/*		GLRandomNumberOperation					*/
/************************************************/
// variable = clip(operand, min, max)
@implementation GLRandomNumberOperation

- (id) initWithResult: (GLFunction *) resultVariable firstScalarOperand: (GLFloat) fsOperand secondScalarOperand: (GLFloat) ssOperand
{
	if (( self = [super initWithResult: @[resultVariable] operand: @[]] ))
	{
        self.firstScalarOperand = fsOperand;
        self.secondScalarOperand = ssOperand;
		NSUInteger nDataElements = resultVariable.nDataElements;
		
		GLFloat amp = (ssOperand - fsOperand);
		GLFloat offset = fsOperand;
		
        BOOL halfComplex = resultVariable.dimensions.count == 2 && [resultVariable.dimensions[0] basisFunction] == kGLExponentialBasis && [resultVariable.dimensions[1] basisFunction] == kGLExponentialBasis && [resultVariable.dimensions[1] isStrictlyPositive];
        
        if (halfComplex)
        {
            NSUInteger kDimNPoints = [resultVariable.dimensions[0] nPoints];
            NSUInteger lDimNPoints = [resultVariable.dimensions[1] nPoints];
            self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				NSMutableData *result = resultArray[0];
                GLFloat *pointerValue = (GLFloat *) result.bytes;
                for (NSUInteger i=0; i<nDataElements; i++) {
                    pointerValue[i] = amp*((double) rand())/( (double) RAND_MAX ) + offset;
                }
                for (NSUInteger i=1; i<kDimNPoints/2; i++) {
                    pointerValue[(kDimNPoints-i)*lDimNPoints] = -pointerValue[i*lDimNPoints];
                }
            };
            self.graphvisDescription = [NSString stringWithFormat: @"hermitian rand (%f, %f)", fsOperand, ssOperand];
        }
        else
        {
            self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				NSMutableData *result = resultArray[0];
                GLFloat *pointerValue = (GLFloat *) result.bytes;
                for (NSUInteger i=0; i<nDataElements; i++) {
                    pointerValue[i] = amp*((double) rand())/( (double) RAND_MAX ) + offset;
                }
            };
            self.graphvisDescription = [NSString stringWithFormat: @"rand (%f, %f)", fsOperand, ssOperand];
        }
	}
	return self;
}

- (BOOL) isEqualToOperation: (id) otherOperation {
    if ( ![super isEqualToOperation: otherOperation] )  {
        return NO;
    }
    
    GLRandomNumberOperation * op = otherOperation;
    if (self.firstScalarOperand != op.firstScalarOperand) {
        return NO;
    }
    if (self.secondScalarOperand != op.secondScalarOperand) {
        return NO;
    }
    
    return YES;
}

@end


/************************************************/
/*		GLNormalDistributionOperation			*/
/************************************************/
// variable = clip(operand, min, max)
@implementation GLNormalDistributionOperation

- (id) initWithResult: (GLFunction *) resultVariable
{
	if (( self = [super initWithResult: @[resultVariable] operand: @[]] ))
	{

		NSUInteger nDataElements = resultVariable.nDataElements;
		self.buffer = @[ [[GLBuffer alloc] initWithLength: nDataElements*sizeof(GLFloat)] ];
        BOOL halfComplex = resultVariable.dimensions.count == 2 && [resultVariable.dimensions[0] basisFunction] == kGLExponentialBasis && [resultVariable.dimensions[1] basisFunction] == kGLExponentialBasis && [resultVariable.dimensions[1] isStrictlyPositive];
        if (nDataElements % 2 == 1) [NSException raise: @"StupidImplementationException" format: @"Can't deal with odd numbers"];
		
		variableOperation op = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
			const int halfTheElements = (int) nDataElements/2;
			
			GLFloat *A = (GLFloat *) [resultArray[0] bytes];
			GLFloat *A2 = &(A[halfTheElements]);
			GLFloat *B = (GLFloat *) [bufferArray[0] bytes];
			GLFloat *B2 = &(B[halfTheElements]);
			
			// Put random numbers (0,1) into the buffer
			for (NSUInteger i=0; i<nDataElements; i++) {
				B[i] = ((double) rand())/( (double) RAND_MAX );
			}
			
			// A should now be fully populated with the sin and cos of the first half of B
			GLFloat twoPi = 2*M_PI;
			vGL_vsmul( B, 1, &twoPi, B, 1, halfTheElements);
			vGL_vvsincos( A, A2, B, &halfTheElements);
			
			// B2 will contain sqrt(-2*log(B2))
			vGL_vvlog( B2, B2, &halfTheElements );
			vGL_vneg( B2, 1, B2, 1, halfTheElements );
			GLFloat two = 2.;
			vGL_vsmul( B2, 1, &two, B2, 1, halfTheElements);
			vGL_vvsqrt( B2, B2, &halfTheElements );
			
			// Now we multiply
			vGL_vmul( B2, 1, A, 1, A, 1, halfTheElements);
			vGL_vmul( B2, 1, A2, 1, A2, 1, halfTheElements);
		};
		
        if (halfComplex)
        {
            NSUInteger kDimNPoints = [resultVariable.dimensions[0] nPoints];
            NSUInteger lDimNPoints = [resultVariable.dimensions[1] nPoints];
            self.operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				
				op(resultArray, operandArray, bufferArray);
				
				GLSplitComplex C = splitComplexFromData( resultArray[0] );
				
				// Hermitian conjugates
                for (NSUInteger i=1; i<kDimNPoints/2; i++) {
                    C.realp[(kDimNPoints-i)*lDimNPoints] = C.realp[i*lDimNPoints];
					C.imagp[(kDimNPoints-i)*lDimNPoints] = -C.imagp[i*lDimNPoints];
					
					C.realp[(kDimNPoints-i)*lDimNPoints+(lDimNPoints-1)] = C.realp[i*lDimNPoints+(lDimNPoints-1)];
					C.imagp[(kDimNPoints-i)*lDimNPoints+(lDimNPoints-1)] = -C.imagp[i*lDimNPoints+(lDimNPoints-1)];
                }
				
				// For the four self-conjugate components, that means that there can be no imaginary component
				C.imagp[0] = 0;
				C.imagp[(lDimNPoints-1)] = 0;
				C.imagp[(kDimNPoints/2)*lDimNPoints+0] = 0;
				C.imagp[(kDimNPoints/2)*lDimNPoints+(lDimNPoints-1)] = 0;
				
				// But that their real components should be doubled, to make all else equal.
				C.realp[0] = 2*C.realp[0];
				C.realp[(lDimNPoints-1)] = 2*C.realp[(lDimNPoints-1)];
				C.realp[(kDimNPoints/2)*lDimNPoints+0] = 2*C.realp[(kDimNPoints/2)*lDimNPoints+0];
				C.realp[(kDimNPoints/2)*lDimNPoints+(lDimNPoints-1)] = 2*C.realp[(kDimNPoints/2)*lDimNPoints+(lDimNPoints-1)];
            };
            self.graphvisDescription = [NSString stringWithFormat: @"hermitian randn"];
        }
        else
        {
            self.operation = op;
            self.graphvisDescription = [NSString stringWithFormat: @"randn"];
        }
	}
	return self;
}

@end