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

- (id) initWithFirstOperand: (GLMutableNetCDFVariable *) fOperand secondOperand: (GLFunction *) sOperand alongDimensionAtIndex: (NSUInteger) mutableDimensionIndex;
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
[NSException raise: @"NotYetImplemented" format: @"This operation needs to mimic the operation below"];
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
    if ([[sOperand class] isSubclassOfClass: [GLFunction class]]) {
        GLFunction *sFunction = (GLFunction *) sOperand;
        if ( fOperand.dimensions.count != sFunction.dimensions.count+1) {
            [NSException raise: @"DimensionsException" format: @"Incompatible number of dimensions for concatenation."];
            return nil;
        }
    } else if ([[sOperand class] isSubclassOfClass: [GLScalar class]]) {
        if ( fOperand.dimensions.count != 1) {
            [NSException raise: @"DimensionsException" format: @"Incompatible number of dimensions for concatenation."];
            return nil;
        }
    } else {
        [NSException raise: @"DimensionsException" format: @"Cannot concatenate a linear transform."];
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
    if ([[sOperand class] isSubclassOfClass: [GLFunction class]]) {
        GLFunction *sFunction = (GLFunction *) sOperand;
        for (NSUInteger i=0; i<reducedDimensions.count; i++)
        {
            GLDimension *leftDim = [reducedDimensions objectAtIndex: i];
            GLDimension *rightDim = [sFunction.dimensions objectAtIndex: i];
            
            // We only require the number of points match.
            if (leftDim.nPoints != rightDim.nPoints) {
                [NSException raise: @"DimensionsException" format: @"The dimensional lengths do not match."];
                return nil;
            }
            [self.indexRanges addObject: [NSValue valueWithRange: NSMakeRange(0, rightDim.nPoints)]];
        }
    }
	
	[self.indexRanges insertObject:[NSValue valueWithRange: NSMakeRange(pointIndex, 1)] atIndex:mutableDimensionIndex];
	
	variableOperation operation;
	GLLowLevelNetCDF *file = fOperand.file;
	NSArray *indexRanges = self.indexRanges;
	if (fOperand.isComplex)
	{
		int variableID = fOperand.variableID;
		int imagpVariableID = fOperand.imagpVariableID;
		NSUInteger nDataPoints = sOperand.nDataPoints;
		
		if (sizeof(GLFloat)==sizeof(float)) {
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLSplitComplex split = splitComplexFromData(operandArray[0]);
				NSData *realp = [NSData dataWithBytes:split.realp length: nDataPoints*sizeof(GLFloat)];
				NSData *imagp = [NSData dataWithBytes:split.imagp length: nDataPoints*sizeof(GLFloat)];
				[file writeFloatData: realp toVariableWithID:variableID atIndexRange:indexRanges];
				[file writeFloatData: imagp toVariableWithID:imagpVariableID atIndexRange: indexRanges];
				[fOperand removeOperation: self];
			};
		} else {
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				GLSplitComplex split = splitComplexFromData(operandArray[0]);
				NSData *realp = [NSData dataWithBytes:split.realp length: nDataPoints*sizeof(GLFloat)];
				NSData *imagp = [NSData dataWithBytes:split.imagp length: nDataPoints*sizeof(GLFloat)];
				[file writeDoubleData: realp toVariableWithID:variableID atIndexRange:indexRanges];
				[file writeDoubleData: imagp toVariableWithID:imagpVariableID atIndexRange: indexRanges];
				[fOperand removeOperation: self];
			};
		}
	}
	else
	{
		int variableID = fOperand.variableID;
		if (sizeof(GLFloat)==sizeof(float)) {
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				[file writeFloatData: operandArray[0] toVariableWithID:variableID atIndexRange:indexRanges];
				[fOperand removeOperation: self];
			};
		} else {
			operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
				[file writeDoubleData: operandArray[0] toVariableWithID:variableID atIndexRange:indexRanges];
				[fOperand removeOperation: self];
			};
		}
	}
	
	if (( self = [super initWithResult: @[] operand: @[sOperand] buffers: @[] operation: operation] ))
	{
		mutableDimension.nPoints = nPointsAfterMutation;
		[fOperand addOperation: self];
	}
    
    return self;
}

@end
