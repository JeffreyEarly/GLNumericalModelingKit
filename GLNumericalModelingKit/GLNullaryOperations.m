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

- (id) initWithResult: (GLVariable *) resultVariable firstScalarOperand: (GLFloat) fsOperand secondScalarOperand: (GLFloat) ssOperand
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
