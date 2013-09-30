//
//  GLSpectralDifferentiationVariable.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 4/22/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import "GLSpectralDifferentialOperator.h"
#import "GLVectorVectorOperations.h"
#import "GLUnaryOperations.h"

/************************************************/
/*		GLSpectralDifferentiationOperation		*/
/************************************************/
// variable = leftVariable * rightVariable
// This differs from the multiplication operation in that it knows how to transform to an alternate basis, if necessary
@interface GLSpectralDifferentiationOperation : GLBinaryOperation
- (id) initWithSpectralOperator: (GLSpectralDifferentialOperator *) fOperand operand: (GLVariable *) sOperand;
@end

/************************************************/
/*		GLSpectralDifferentialOperator			*/
/************************************************/


@implementation GLSpectralDifferentialOperator

- (BOOL) basisTransformationRequired
{
	BOOL match = YES;
	for (NSUInteger i=0; i<self.fromDimensions.count; i++) {
		match &= [self.fromDimensions[i] isEqualToDimension: self.toDimensions[i]];
	}
	return !match;
}

- (GLVariableOperation *) differentiationOperationFromVariable: (GLVariable *) operand
{
    GLBinaryOperation *operation;
	if (self.basisTransformationRequired) {
		operation = [[GLSpectralDifferentiationOperation alloc] initWithSpectralOperator: self operand: operand];
	} else {
		operation = [[GLMultiplicationOperation alloc] initWithFirstOperand: self secondOperand: operand];
	}
    
    return operation;
}

@end

/************************************************/
/*		GLSpectralDifferentiationOperation		*/
/************************************************/

@interface GLSpectralDifferentiationOperation ()
@property(readwrite) BOOL canOperateInPlace;
@end

@implementation GLSpectralDifferentiationOperation

@synthesize canOperateInPlace;

- (id) initWithSpectralOperator: (GLSpectralDifferentialOperator *) fOperand operand: (GLVariable *) sOperand
{
	BOOL isComplex = fOperand.isComplex || sOperand.isComplex;
	GLDataFormat format = isComplex ? kGLSplitComplexDataFormat : kGLRealDataFormat;
	GLVariable *result = [GLVariable variableOfType: format withDimensions: fOperand.toDimensions forEquation: fOperand.equation];
	
    if (sOperand.name && fOperand.name) {
        result.name = [NSString stringWithFormat: @"%@_%@", sOperand.name, fOperand.name];
    }
    
	if (( self = [super initWithResult: result firstOperand:fOperand secondOperand:sOperand] )) {
		self.result.isPurelyReal = (self.firstOperand.isPurelyReal && self.secondOperand.isPurelyReal) || (self.firstOperand.isPurelyImaginary && self.secondOperand.isPurelyImaginary);
		self.result.isPurelyImaginary= (self.firstOperand.isPurelyReal && self.secondOperand.isPurelyImaginary) || (self.firstOperand.isPurelyImaginary && self.secondOperand.isPurelyReal);
		
		// First we determine how many elements over we need to shift the result.
		// This is necessary because shifting between a sine and cosine basis results
		// in a different set of wavenumbers available. A positive shift indicates that the
		// result array shifts forward, a negative shift indicates that the operator shifts forward.
		NSInteger totalShift = 0;
		NSInteger stride = 1;
		for (NSInteger i=fOperand.fromDimensions.count-1; i>=0; i--) {
			GLBasisFunction fromBasis = [fOperand.fromDimensions[i] basisFunction];
			GLBasisFunction toBasis = [fOperand.toDimensions[i] basisFunction];
			if ( (fromBasis == kGLCosineBasis || fromBasis == kGLCosineHalfShiftBasis) && (toBasis == kGLSineBasis || toBasis == kGLSineHalfShiftBasis))
			{
				totalShift -= stride;
			}
			else if ( (fromBasis == kGLSineBasis || fromBasis == kGLSineHalfShiftBasis) && (toBasis == kGLCosineBasis || toBasis == kGLCosineHalfShiftBasis))
			{
				totalShift += stride;
			}
			
			stride = stride * [fOperand.fromDimensions[i] nPoints];
		}
		
		BOOL shiftResult = totalShift > 0 ? YES : NO; // Yes if we shift the result, no if we shift the operand
		totalShift = labs(totalShift); // total number of elements to shift by
		NSInteger totalShiftBytes = totalShift*sizeof(GLFloat); // total number of bytes to shift by
		NSInteger nDataElements = self.result.nDataElements - totalShift;
		NSInteger nDataPoints = self.result.nDataPoints - totalShift;
		NSInteger totalEndShiftBytes = nDataPoints*sizeof(GLFloat);
		
		if ( !self.firstOperand.isComplex && !self.secondOperand.isComplex )
		{
			if (shiftResult) {
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					vGL_vmul( fOperand.bytes, 1, sOperand.bytes, 1, result.mutableBytes+totalShiftBytes, 1, nDataElements);
					vGL_vclr( result.mutableBytes, 1, totalShift );
				};
			} else {
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					vGL_vmul( fOperand.bytes+totalShiftBytes, 1, sOperand.bytes+totalShiftBytes, 1, result.mutableBytes, 1, nDataElements);
					vGL_vclr( result.mutableBytes+totalEndShiftBytes, 1, totalShift );
				};
			}
			self.canOperateInPlace = YES;
		}
		else if ( !self.firstOperand.isComplex && self.secondOperand.isComplex )
		{
			if (shiftResult) {
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					GLSplitComplex rightComplex = splitComplexFromData( sOperand );
					GLSplitComplex destComplex = splitComplexFromData( result );
					vGL_vmul( rightComplex.realp, 1, fOperand.bytes, 1, ((void*)destComplex.realp)+totalShiftBytes, 1, nDataPoints);
					vGL_vclr( destComplex.realp, 1, totalShift );
					
					vGL_vmul( rightComplex.imagp, 1, fOperand.bytes, 1, ((void*)destComplex.imagp)+totalShiftBytes, 1, nDataPoints);
					vGL_vclr( destComplex.imagp, 1, totalShift );
				};
			} else {
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					
					GLSplitComplex rightComplex = splitComplexFromData( sOperand );
					GLSplitComplex destComplex = splitComplexFromData( result );
										
					vGL_vmul( ((void*)rightComplex.realp)+totalShiftBytes, 1, fOperand.bytes+totalShiftBytes, 1, destComplex.realp, 1, nDataPoints);
					vGL_vclr( ((void*)destComplex.realp)+totalEndShiftBytes, 1, totalShift );
					
					vGL_vmul( ((void*)rightComplex.imagp)+totalShiftBytes, 1, fOperand.bytes+totalShiftBytes, 1, destComplex.imagp, 1, nDataPoints);
					vGL_vclr( ((void*)destComplex.imagp)+totalEndShiftBytes, 1, totalShift );
				};
			}
			
			self.canOperateInPlace = YES;
		} else {
			[NSException raise: @"MethodNotImplemented" format: @"MethodNotImplemented"];
		}
		
		return self;
		
		if ( self.firstOperand.isComplex && !self.secondOperand.isComplex ) {
			NSUInteger nDataPoints = self.result.nDataPoints;
			self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
				GLSplitComplex leftComplex = splitComplexFromData( fOperand );
				GLSplitComplex destComplex = splitComplexFromData( result );
				// (a + i b)*(x) = a x + i b x
				vGL_vmul( leftComplex.realp, 1, sOperand.bytes, 1, destComplex.realp, 1, nDataPoints);
				vGL_vmul( leftComplex.imagp, 1, sOperand.bytes, 1, destComplex.imagp, 1, nDataPoints);
			};
			
			self.canOperateInPlace = YES;
		} else {
			NSUInteger nDataPoints = self.result.nDataPoints;
			
			if ((!self.firstOperand.isPurelyReal && !self.firstOperand.isPurelyImaginary) && self.secondOperand.isPurelyReal) {
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					
					GLSplitComplex leftComplex = splitComplexFromData( fOperand );
					GLSplitComplex rightComplex = splitComplexFromData( sOperand );
					GLSplitComplex destComplex = splitComplexFromData( result );
					
					// (a + i b)*(x) = a x + i b x
					vGL_vmul( leftComplex.realp, 1, rightComplex.realp, 1, destComplex.realp, 1, nDataPoints);
					vGL_vmul( leftComplex.imagp, 1, rightComplex.realp, 1, destComplex.imagp, 1, nDataPoints);					
				};
				self.canOperateInPlace = NO;
			} else if ((!self.firstOperand.isPurelyReal && !self.firstOperand.isPurelyImaginary) && self.secondOperand.isPurelyImaginary) {
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					
					GLSplitComplex leftComplex = splitComplexFromData( fOperand );
					GLSplitComplex rightComplex = splitComplexFromData( sOperand );
					GLSplitComplex destComplex = splitComplexFromData( result );
					
					// (a + i b)*(i y) = - b y + i a y
					vGL_vmul( leftComplex.imagp, 1, rightComplex.imagp, 1, destComplex.realp, 1, nDataPoints);
					vGL_vneg( destComplex.realp, 1, destComplex.realp, 1, nDataPoints );
					vGL_vmul( leftComplex.realp, 1, rightComplex.imagp, 1, destComplex.imagp, 1, nDataPoints);					
				};
				self.canOperateInPlace = NO;
			} else if (self.firstOperand.isPurelyReal && (!self.secondOperand.isPurelyReal && !self.secondOperand.isPurelyImaginary)) {
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					
					GLSplitComplex leftComplex = splitComplexFromData( fOperand );
					GLSplitComplex rightComplex = splitComplexFromData( sOperand );
					GLSplitComplex destComplex = splitComplexFromData( result );
					
					// a *(x + i y) = a x + i a y
					vGL_vmul( leftComplex.realp, 1, rightComplex.realp, 1, destComplex.realp, 1, nDataPoints);
					vGL_vmul( leftComplex.realp, 1, rightComplex.imagp, 1, destComplex.imagp, 1, nDataPoints);					
				};
				self.canOperateInPlace = NO;
			} else if (self.firstOperand.isPurelyImaginary && (!self.secondOperand.isPurelyReal && !self.secondOperand.isPurelyImaginary)) {
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					
					GLSplitComplex leftComplex = splitComplexFromData( fOperand );
					GLSplitComplex rightComplex = splitComplexFromData( sOperand );
					GLSplitComplex destComplex = splitComplexFromData( result );
					
					// i b *(x + i y) = - b y + i b x
					vGL_vmul( leftComplex.imagp, 1, rightComplex.imagp, 1, destComplex.realp, 1, nDataPoints);
					vGL_vneg( destComplex.realp, 1, destComplex.realp, 1, nDataPoints );
					vGL_vmul( leftComplex.imagp, 1, rightComplex.realp, 1, destComplex.imagp, 1, nDataPoints);					
				};
				self.canOperateInPlace = NO;
			} else {
				self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
					
					GLSplitComplex leftComplex = splitComplexFromData( fOperand );
					GLSplitComplex rightComplex = splitComplexFromData( sOperand );
					GLSplitComplex destComplex = splitComplexFromData( result );
					
					// (a + i b)*(x + i y) = (a x - b y) + i (b x + ay)
					vGL_vmmsb( leftComplex.realp, 1, rightComplex.realp, 1, leftComplex.imagp, 1, rightComplex.imagp, 1, destComplex.realp, 1, nDataPoints);
					vGL_vmma( leftComplex.imagp, 1, rightComplex.realp, 1, leftComplex.realp, 1, rightComplex.imagp, 1, destComplex.imagp, 1, nDataPoints);		
				};
				self.canOperateInPlace = NO;
			}
			//#warning overriding this for testing purposes.
			//			self.canOperateInPlace = NO;
		}
	}
    return self;
}

@end



