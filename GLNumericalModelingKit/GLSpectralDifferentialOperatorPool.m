//
//  GLSpectralDifferentialOperatorPool.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 4/26/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import "GLSpectralDifferentialOperatorPool.h"
#import "GLDimension.h"
#import "GLMemoryPool.h"
#import "GLVariable.h"
#import "GLSpectralMatrices.h"
#import "GLEquation.h"
#import "GLVectorVectorOperations.h"
#import "GLBasisTransformationOperations.h"

@implementation GLSpectralDifferentialOperatorPool

- (GLDifferentialOperator *) differentialOperatorWithName: (NSString *) opName
{
	GLDifferentialOperator *diffOp = [super differentialOperatorWithName: opName];
	
	if (!diffOp)
	{
        if ([opName isEqualToString: @"harmonicOperator"]) {
            diffOp = [self harmonicOperator];
        } else {
            NSMutableArray *derivatives = [NSMutableArray arrayWithCapacity: self.differentiationDimensions.count];
            NSString *reducedString = [opName copy];
            for (GLDimension *dim in self.differentiationDimensions)
            {
                NSString *newReducedString = [reducedString stringByReplacingOccurrencesOfString: dim.name withString: @""];
                [derivatives addObject: [NSNumber numberWithUnsignedInteger: reducedString.length - newReducedString.length]];
                reducedString = newReducedString;
            }
            // We check to make sure we can account for the entire string
            if (reducedString.length == 0)
            {
                diffOp = [self operatorWithDerivatives: derivatives];
            }
        }
		// We check to make sure we can account for the entire string
		if (diffOp)
		{
            diffOp.name = opName;
			[self setDifferentialOperator: diffOp forName: opName];
		}
		else
		{
			[NSException raise: @"BadDifferentiationRequest" format: @"Cannot find the differential operator: %@", opName];
		}
	}
	
	return diffOp;
}

- (GLSpectralDifferentialOperator *) operatorWithDerivatives: (NSArray *) derivatives
{
	GLSpectralDifferentialOperator *diffOp;
	NSMutableArray *destinationBasis = [NSMutableArray array];
	NSUInteger iFactors = 0; // total number of time sqrt(-1) gets included.
	NSUInteger negFactors = 0; // total number of times a -1 gets included.
	
	for (NSUInteger i=0; i<derivatives.count; i++)
	{
		NSUInteger numDeriv = [derivatives[i] unsignedIntegerValue];
		GLDimension *dim = self.transformedDimensions[i];
		
		if (!numDeriv) {
			[destinationBasis addObject: dim];
		} else {
			if (dim.basisFunction == kGLExponentialBasis) {
				iFactors += numDeriv;
				[destinationBasis addObject: dim];
			} else if (dim.basisFunction == kGLSineBasis || dim.basisFunction == kGLSineHalfShiftBasis) {
				negFactors += (numDeriv % 4 > 1);
				if (numDeriv % 2 == 1) {
					GLBasisFunction newBasis = dim.basisFunction == kGLSineBasis ? kGLCosineBasis : kGLCosineHalfShiftBasis;
					GLDimension *newDim = [[GLDimension alloc] initAsDimension: dim transformedToBasis: newBasis strictlyPositive: YES];
					[destinationBasis addObject: newDim];
				} else {
					[destinationBasis addObject: dim];
				}
			}  else if (dim.basisFunction == kGLCosineBasis || dim.basisFunction == kGLCosineHalfShiftBasis) {
				negFactors += (numDeriv % 4 == 1 || numDeriv % 4 == 2) ? 1 : 0;
				if (numDeriv % 2 == 1) {
					GLBasisFunction newBasis = dim.basisFunction == kGLCosineBasis ? kGLSineBasis : kGLSineHalfShiftBasis;
					GLDimension *newDim = [[GLDimension alloc] initAsDimension: dim transformedToBasis: newBasis strictlyPositive: YES];
					[destinationBasis addObject: newDim];
				} else {
					[destinationBasis addObject: dim];
				}
			}
		}
	}
	
	
	for (NSUInteger i=0; i<derivatives.count; i++)
	{
		NSUInteger numDeriv = [derivatives[i] unsignedIntegerValue];
		GLDimension *dim = self.transformedDimensions[i];
		
		if (numDeriv) {
            GLSpectralDifferentialOperator *var;
            if (numDeriv % 2 == 1 && (dim.basisFunction == kGLExponentialBasis || dim.basisFunction == kGLSineBasis || dim.basisFunction == kGLSineHalfShiftBasis))
            { // The goal here is to stamp out (set to zero) the n/2 component
				if (dim.basisFunction == kGLExponentialBasis)
				{
					GLFloat *f = dim.data.mutableBytes;
					NSUInteger index = dim.isStrictlyPositive ? dim.nPoints-1 : dim.nPoints/2;
					GLFloat oldValue = f[index];
					f[index] = 0.0;
					var = [GLSpectralDifferentialOperator variableOfRealTypeFromDimension: dim withDimensions: self.transformedDimensions forEquation: self.equation];
					f[index] = oldValue;
				}
				else if (dim.basisFunction == kGLSineBasis || dim.basisFunction == kGLSineHalfShiftBasis)
				{	// Stamp the nyquist frequency
					GLFloat *f = dim.data.mutableBytes;
					NSUInteger index = dim.nPoints-1;
					GLFloat oldValue = f[index];
					f[index] = 0.0;
					var = [GLSpectralDifferentialOperator variableOfRealTypeFromDimension: dim withDimensions: self.transformedDimensions forEquation: self.equation];
					f[index] = oldValue;
				}
            }
            else {
                var = [GLSpectralDifferentialOperator variableOfRealTypeFromDimension: dim withDimensions: self.transformedDimensions forEquation: self.equation];
            }
			
			var.fromDimensions = self.transformedDimensions;
			var.toDimensions = destinationBasis;
			
            var = [[var scalarMultiply: 2*M_PI] pow: numDeriv];
            
			if (!diffOp) {
				diffOp = (GLSpectralDifferentialOperator *) var;
			} else {
				diffOp = (GLSpectralDifferentialOperator *) [diffOp multiply: var];
			}
		}
	}
	
	if ((iFactors % 2) == 1) {
		diffOp = [diffOp swapComplex];
	}
    
    if (iFactors % 4 > 1) {
        negFactors += 1;
    }
	
	if (negFactors % 2 ) {
		diffOp = [diffOp negate];
	}
	
	diffOp.isPurelyReal = (iFactors % 2) == 0;
	diffOp.isPurelyImaginary = (iFactors % 2) == 1;
	diffOp.fromDimensions = self.transformedDimensions;
	diffOp.toDimensions = destinationBasis;
	
	return diffOp;
}

//- (GLSpectralDifferentialOperator *) operatorWithDerivatives: (NSArray *) derivatives
//{
//	if (self.differentiationDimensions.count != 2 || derivatives.count != 2 ) {
//		[NSException raise: @"BadDiffMatrixDimensionsException" format: @"Cannot transform anything other than two dimensions!"];
//		return nil;
//	}
//	
//	GLSpectralDifferentialOperator *diffOp = [[GLSpectralDifferentialOperator alloc] initVariableOfType: kGLSplitComplexDataFormat withDimensions: self.transformedDimensions forEquation: self.equation];
//	
//	GLSplitComplex split = diffOp.splitComplex;
//	int diffX = [[derivatives objectAtIndex: 0] intValue];
//	int diffY = [[derivatives objectAtIndex: 1] intValue];
//	int nx = (int) [[self.differentiationDimensions objectAtIndex: 0] nPoints];
//	int ny = (int) [[self.differentiationDimensions objectAtIndex: 1] nPoints];
//	GLFloat xdomain = [[self.differentiationDimensions objectAtIndex: 0] domainLength];
//	GLFloat ydomain = [[self.differentiationDimensions objectAtIndex: 1] domainLength];
//	
//	spectralDifferentiationMatrixComplex( &split, diffX, diffY, nx, ny, xdomain, ydomain);
//	
//	diffOp.isPurelyReal = ((diffX+diffY) % 2) == 0;
//	diffOp.isPurelyImaginary = ((diffX+diffY) % 2) == 1;
//	
//	return diffOp;
//}


- (GLSpectralDifferentialOperator *) harmonicOperator
{
	return [self harmonicOperatorOfOrder: 1];
}

- (GLSpectralDifferentialOperator *) harmonicOperatorOfOrder: (NSUInteger) order
{	
	NSMutableArray *zeros = [NSMutableArray array];
	for (NSUInteger i=0; i<self.differentiationDimensions.count; i++) {
		[zeros addObject: @0];
	}
	
	// Build the operators 
	GLVariable *harmonicOperator=nil;;
	for (NSUInteger i=0; i<self.differentiationDimensions.count; i++) {
		NSMutableArray *deriv = [zeros mutableCopy];
		deriv[i] = @2;
		if (harmonicOperator) {
			harmonicOperator = [harmonicOperator plus: [self operatorWithDerivatives: deriv]];
		} else {
			harmonicOperator = [self operatorWithDerivatives: deriv];
		}
	}
	
	for (NSUInteger i=1; i<order; i++) {
		GLMultiplicationOperation *operation = [[GLMultiplicationOperation alloc] initWithFirstOperand: harmonicOperator secondOperand: harmonicOperator];
		harmonicOperator = operation.result;
	}
	
	[self.equation solveForVariable: harmonicOperator waitUntilFinished: YES];
	
//	GLSpectralDifferentialOperator *diffOp = [GLSpectralDifferentialOperator variableFromVariable: harmonicOperator];
	
    harmonicOperator.name = [NSString stringWithFormat: @"nabla^%lu", order];
	return harmonicOperator;
}

//- (GLSpectralDifferentialOperator *) spectralVanishingViscosityOperatorWithCutoff: (GLFloat) sigma
//{
//#warning Only valid for two dimensions at the moment
//	GLSpectralDifferentialOperator *diffOp = [[GLSpectralDifferentialOperator alloc] initVariableOfType: kGLRealDataFormat withDimensions: self.transformedDimensions forEquation: self.equation];
//	
//	int nx = (int) [[self.transformedDimensions objectAtIndex: 0] nPoints];
//	int ny = (int) [[self.transformedDimensions objectAtIndex: 1] nPoints];
//	
//	spectralVanishingViscoityMatrix(diffOp.pointerValue, sigma, nx, ny);
//	
//	return diffOp;
//}

- (GLSpectralDifferentialOperator *) spectralVanishingViscosityFilter
{	
	GLFloat minNyquist;
	GLFloat minSampleInterval;
	
	GLSpectralDifferentialOperator *bigK;
	for (GLDimension *dim in self.transformedDimensions) {
		GLSpectralDifferentialOperator *k = [GLSpectralDifferentialOperator variableOfRealTypeFromDimension: dim withDimensions: self.transformedDimensions forEquation: self.equation];
		GLFloat nyquist = dim.domainMin + dim.domainLength;
		if (!bigK) {
			bigK = [k multiply: k];
			minNyquist = nyquist;
			minSampleInterval = dim.sampleInterval;
		} else {
			bigK = [bigK plus: [k multiply: k]];
			if (nyquist < minNyquist) {
				minNyquist = nyquist;
				minSampleInterval = dim.sampleInterval;
			}
		}
	}
	bigK = [bigK sqrt];
	
	if ([GLBasisTransformOperation shouldAntialias]) {
		minNyquist = 2*minNyquist/3;
	}
	
    // Note that the 'minSampleInterval' is deltaK (not deltaX)
	GLFloat wavenumberCutoff = minSampleInterval * pow(minNyquist/minSampleInterval, 0.75);
	
	
	// expf( -alpha * powf( (k-max)/(k-cutoff), p) );
	
	GLSpectralDifferentialOperator *filter = [[bigK scalarAdd: -minNyquist] multiply: [[bigK scalarAdd: -wavenumberCutoff] scalarDivide: 1.0]];
	filter = [[[filter pow: 2] negate] exponentiate];
	
	[filter solve];
	GLFloat *kk = [bigK pointerValue];
	GLFloat *f = [filter pointerValue];
	for (NSUInteger i=0; i<filter.nDataPoints; i++) {
		if ( kk[i] < wavenumberCutoff ) {
			f[i] = 0;
		} else if ( kk[i] > minNyquist ) {
			f[i] = 1.0;
		}
	}
	
	return filter;
}

/************************************************/
/*		Object Messaging Voodoo					*/
/************************************************/

#pragma mark -
#pragma mark Object Messaging Voodoo
#pragma mark

- (NSMethodSignature *)methodSignatureForSelector:(SEL)selector
{
	NSMethodSignature *signature = [super methodSignatureForSelector:selector];
	if (!signature)
	{
		if ([self differentialOperatorWithName:  NSStringFromSelector(selector)])
		{
			signature = [self methodSignatureForSelector:@selector(differentialOperatorWithName:)];
		}
	}
	
	return signature;
}

- (void)forwardInvocation:(NSInvocation *)anInvocation
{
	// We're forwarding ALL unknown invocations as if they're attempts at differentiating.
	NSString *argument = NSStringFromSelector([anInvocation selector]);
	[anInvocation setSelector: @selector(differentialOperatorWithName:)];
	[anInvocation setArgument: &argument atIndex:2];
	[anInvocation invokeWithTarget: self];
}

@dynamic x, y, xx, xy, yy, xxx, xxy, xyy, yyy;


@end

@implementation GLVariable (SpectralDifferentiationExtensions)

@dynamic x, y, xx, xy, yy, xxx, xxy, xyy, yyy;

@end
