//
//  GLNetCDFConcatenationOperation.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 3/13/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import "GLNetCDFConcatenationOperation.h"
#import "GLNetCDFVariable.h"
#import "GLLowLevelNetCDF.h"

@implementation GLNetCDFConcatenationOperation

@synthesize indexRanges;

- (id) initWithFirstOperand: (GLMutableNetCDFVariable *) fOperand secondOperand: (GLVariable *) sOperand alongDimensionAtIndex: (NSUInteger) mutableDimensionIndex;
{
	// At this point we don't know how to handle variables that aren't of the same type.
	if ( fOperand.isComplex != sOperand.isComplex ) {
		[NSException raise: @"MethodNotYetImplemented" format: @"Cannot perform a binary operation on two variables that are not either both complex or both real."];
		return nil;
	}
	
	// By assumption these variables have the same number of dimensions.
	if ( fOperand.dimensions.count != sOperand.dimensions.count) {
		[NSException raise: @"DimensionsException" format: @"Incompatible number of dimensions for concatenation."];
		return nil;
	}
	
	// The mutable dimension had better be a dimension of the first operand.
	if (mutableDimensionIndex >= fOperand.dimensions.count) {
		[NSException raise: @"DimensionsException" format: @"Mutable dimension not found."];
		return nil;
	}
	
	GLMutableDimension *mutableDimension = [fOperand.dimensions objectAtIndex: mutableDimensionIndex];
	if (!mutableDimension.isMutable) {
		[NSException raise: @"DimensionsException" format: @"Mutable dimension not found."];
		return nil;
	}
	
	// We need to determine the number of points by which we will be extending the mutable dimension,
	NSUInteger nPointsAfterMutation=0;
	NSArray *newPoints;
	
	// and we need to determine the index ranges where the data in the secondOperand will be written.
	self.indexRanges = [NSMutableArray array];
	for (NSUInteger i=0; i<fOperand.dimensions.count; i++)
	{
		GLDimension *leftDim = [fOperand.dimensions objectAtIndex: i];
		GLDimension *rightDim = [sOperand.dimensions objectAtIndex: i];
		
		if (i==mutableDimensionIndex)
		{
			if (leftDim.isEvenlySampled) {
				nPointsAfterMutation = leftDim.nPoints+rightDim.nPoints;
			} else {
				newPoints = rightDim.points;
			}
			
			[self.indexRanges addObject: [NSValue valueWithRange: NSMakeRange(leftDim.nPoints, rightDim.nPoints)]];
		}
		else
		{	// We only require the number of points match.
			if (leftDim.nPoints != rightDim.nPoints) {
				[NSException raise: @"DimensionsException" format: @"The dimensional lengths do not match."];
				return nil;
			}
			[self.indexRanges addObject: [NSValue valueWithRange: NSMakeRange(0, rightDim.nPoints)]];
		}
	}	

	if (( self = [super initWithResult: @[fOperand] operand: @[fOperand,sOperand]] ))
	{
		// Now we mutate the dimension
		if (newPoints) {
			[mutableDimension addPointsFromArray: newPoints];
		} else {
			mutableDimension.nPoints = nPointsAfterMutation;
		}
	}
    
    return self;
}

- (id) initWithFirstOperand: (GLMutableNetCDFVariable *) fOperand lowerDimensionalSecondOperand: (GLVariable *) sOperand alongDimensionAtIndex: (NSUInteger) mutableDimensionIndex index: (NSUInteger) pointIndex
{
	// The mutable dimension had better be a dimension of the first operand.
	if (mutableDimensionIndex >= fOperand.dimensions.count) {
		[NSException raise: @"DimensionsException" format: @"Mutable dimension not found."];
		return nil;
	}
	
	GLMutableDimension *mutableDimension = [fOperand.dimensions objectAtIndex: mutableDimensionIndex];
	if (!mutableDimension.isMutable) {
		[NSException raise: @"DimensionsException" format: @"Mutable dimension not found."];
		return nil;
	}
	
	// At this point we don't know how to handle variables that aren't of the same type.
	if ( fOperand.isComplex != sOperand.isComplex ) {
		[NSException raise: @"MethodNotYetImplemented" format: @"Cannot perform a binary operation on two variables that are not either both complex or both real."];
		return nil;
	}
	
	// By assumption the second operand has one less dimension
	if ( fOperand.dimensions.count != sOperand.dimensions.count+1) {
		[NSException raise: @"DimensionsException" format: @"Incompatible number of dimensions for concatenation."];
		return nil;
	}
	
	// We *may* need to increase the number of points, but then again, maybe not.
	NSUInteger nPointsAfterMutation = pointIndex+1>mutableDimension.nPoints ? pointIndex+1:mutableDimension.nPoints;
	
	// This only handles evenly sampled dimensions, unless they're already extended.
	if (!mutableDimension.isEvenlySampled && nPointsAfterMutation != mutableDimension.nPoints) {
		[NSException raise: @"DimensionsException" format: @"If the mutable dimension is not evenly sampled, it must already be long enough to handle to new data."];
		return nil;
	}
	
	NSMutableArray *reducedDimensions = [NSMutableArray arrayWithArray: fOperand.dimensions];
	[reducedDimensions removeObjectAtIndex: mutableDimensionIndex];
	
	self.indexRanges = [NSMutableArray array];
	for (NSUInteger i=0; i<reducedDimensions.count; i++)
	{
		GLDimension *leftDim = [reducedDimensions objectAtIndex: i];
		GLDimension *rightDim = [sOperand.dimensions objectAtIndex: i];
		
		// We only require the number of points match.
		if (leftDim.nPoints != rightDim.nPoints) {
			[NSException raise: @"DimensionsException" format: @"The dimensional lengths do not match."];
			return nil;
		}
		[self.indexRanges addObject: [NSValue valueWithRange: NSMakeRange(0, rightDim.nPoints)]];
	}
	
	[self.indexRanges insertObject:[NSValue valueWithRange: NSMakeRange(pointIndex, 1)] atIndex:mutableDimensionIndex];
	
	if (( self = [super initWithResult: @[fOperand] operand: @[fOperand,sOperand]] ))
	{
		mutableDimension.nPoints = nPointsAfterMutation;
	}
    
    return self;
}

@synthesize result;
@synthesize firstOperand;
@synthesize secondOperand;

- (void) setupDependencies
{
	for ( NSOperation *op in self.firstOperand.pendingOperations ) {
		[self addDependency: op];
	}
	for ( NSOperation *op in self.secondOperand.pendingOperations ) {
		[self addDependency: op];
	}

	[self.result[0] addOperation: self];
}

- (void) tearDownDependencies
{	
	[self.result[0] removeOperation: self];
	self.result = nil;
	self.firstOperand = nil;
	self.secondOperand = nil;
}

-(void) main
{
	if (self.firstOperand.isComplex)
	{
		GLSplitComplex split = self.secondOperand.splitComplex;
		NSData *imagp = [NSData dataWithBytes:split.imagp length:self.secondOperand.nDataPoints*sizeof(GLFloat)];
		if (sizeof(GLFloat)==sizeof(float)) {
			[self.firstOperand.file writeFloatData:self.secondOperand.data toVariableWithID:self.firstOperand.variableID atIndexRange:self.indexRanges];
			[self.firstOperand.file writeFloatData:imagp toVariableWithID:self.firstOperand.imagpVariableID atIndexRange:self.indexRanges];
		} else {
			[self.firstOperand.file writeDoubleData:self.secondOperand.data toVariableWithID:self.firstOperand.variableID atIndexRange:self.indexRanges];
			[self.firstOperand.file writeDoubleData:imagp toVariableWithID:self.firstOperand.imagpVariableID atIndexRange:self.indexRanges];
		}
	}
	else
	{
		if (sizeof(GLFloat)==sizeof(float)) {
			[self.firstOperand.file writeFloatData:self.secondOperand.data toVariableWithID:self.firstOperand.variableID atIndexRange:self.indexRanges];
		} else {
			[self.firstOperand.file writeDoubleData:self.secondOperand.data toVariableWithID:self.firstOperand.variableID atIndexRange:self.indexRanges];
		}
	}
	[self tearDownDependencies];
}

@end
