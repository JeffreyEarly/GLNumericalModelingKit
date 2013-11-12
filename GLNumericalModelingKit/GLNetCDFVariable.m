//
//  GLNetCDFVariable.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 10/21/11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#import "GLNetCDFVariable.h"

#import "GLDimension.h"
#import "GLVariableOperations.h"
#import "GLMemoryPool.h"
#import "GLNetCDFFile.h"
#import "GLNetCDFFetchDataOperation.h"
#import "GLNetCDFConcatenationOperation.h"
#import "GLEquation.h"

@implementation GLNetCDFVariable
{
	GLLowLevelNetCDF *file;
	int variableID;
}

/************************************************/
/*		Class Methods							*/
/************************************************/

#pragma mark -
#pragma mark Class Methods
#pragma mark

+ (id) variableOfType: (GLDataFormat) dataFormat withDimensions: (NSArray *) theDimensions forEquation: (GLEquation *) equation;
{
	BOOL isMutable = NO;
	for (GLDimension *dimension in theDimensions) {
		isMutable |= dimension.isMutable;
	}
	
	Class GLVariableClass = isMutable ? [GLMutableNetCDFVariable class] : [GLNetCDFVariable class];
	
	return [[GLVariableClass alloc] initVariableOfType: dataFormat withDimensions: theDimensions forEquation: equation];
}

+ (id) variableWithVariable: (GLVariable *) existing
{
	GLNetCDFVariable *variable = [self variableOfType: existing.isComplex withDimensions: existing.dimensions forEquation: existing.equation];
	variable.name = existing.name;
	variable.units = existing.units;
	[variable.metadata addEntriesFromDictionary: existing.metadata];
	variable.uniqueID = existing.uniqueID;

	return variable;
}

/************************************************/
/*		Properties								*/
/************************************************/

#pragma mark -
#pragma mark Properties
#pragma mark

@synthesize file;
@synthesize variableID;
@synthesize imagpVariableID;
@synthesize uniqueID;

- (void) setupDependency
{
	NSMutableArray *ranges = [NSMutableArray array];
	for (GLDimension *dimension in self.dimensions) {
		[ranges addObject: [NSValue valueWithRange: NSMakeRange(0, dimension.nPoints)]];
	}
	GLNetCDFFetchDataOperation *operation = [[GLNetCDFFetchDataOperation alloc] initWithResult: self netCDFVariable: self indexRange: ranges flatten: NO];
	[self addOperation: operation];
}

/************************************************/
/*		NSCopying								*/
/************************************************/

#pragma mark -
#pragma mark NSCopying
#pragma mark

- (id)copyWithZone:(NSZone *)zone
{
	GLNetCDFVariable *copy = [GLNetCDFVariable variableWithVariable: self];
	copy.file = self.file;
	copy.variableID = self.variableID;
	copy.imagpVariableID = self.imagpVariableID;
	copy.uniqueID = self.uniqueID;
	[copy setupDependency];
	
	return copy;
}

/************************************************/
/*		GLVariable Overrides					*/
/************************************************/

#pragma mark -
#pragma mark GLVariable Overrides
#pragma mark

- (GLVariable *) variableFromIndexRange: (NSArray *) ranges
{
	GLNetCDFFetchDataOperation *operation = [[GLNetCDFFetchDataOperation alloc] initWithNetCDFVariable: self indexRange: ranges flatten: YES];
	return operation.result[0];
}

- (NSString *) matrixDescription
{
	NSUInteger n = [[self.dimensions lastObject] nPoints];
	NSMutableString *descrip = [NSMutableString string];
	
	GLFloat max, min;
	vGL_maxv( self.data.mutableBytes, 1, &max, self.nDataElements);
	vGL_minv( self.data.mutableBytes, 1, &min, self.nDataElements);
	
	if ( fabs(min) > max) {
		max = fabs(min);
	}
	
	GLFloat divisor = pow(10, floor(log10(max)));
	if ( divisor == 0.0) divisor = 1;
	
	if (0 && self.dataFormat == kGLSplitComplexDataFormat)
	{
		GLSplitComplex splitComplex = self.splitComplex;
		[descrip appendFormat: @"%f * ", divisor];
		//		for (NSUInteger i=0; i<self.nDataPoints; i++)
		//		{
		//			if ( i % n == 0 ) {
		//				[descrip appendFormat: @"\n"];
		//			}
		//			[descrip appendFormat: @"%1.1f ", sqrt(fabs(splitComplex.realp[i] * splitComplex.realp[i] - splitComplex.imagp[i] * splitComplex.imagp[i]))/divisor];
		//		}
		
		for (NSUInteger i=0; i<self.nDataPoints; i++)
		{
			if ( i % n == 0 ) {
				[descrip appendFormat: @"\n"];
			}
			[descrip appendFormat: @"%+1.1f ", splitComplex.realp[i]/divisor];
		}
		
		[descrip appendFormat: @" imagp \n"];
		for (NSUInteger i=0; i<self.nDataPoints; i++)
		{
			if ( i % n == 0 ) {
				[descrip appendFormat: @"\n"];
			}
			[descrip appendFormat: @"%+1.1f ", splitComplex.imagp[i]/divisor];
		}
	}
	else
	{
		GLFloat *f = self.pointerValue;
		[descrip appendFormat: @"%g * ", divisor];
		for (NSUInteger i=0; i<self.nDataElements; i++)
		{
			if ( i % n == 0 ) {
				[descrip appendFormat: @"\n"];
			}
			[descrip appendFormat: @"%+1.1f ", f[i]/divisor];
		}
	}
	
	return descrip;
}

- (void) dumpToConsole
{
	NSLog(@"%@", [self matrixDescription]);
}



- (NSString *) description
{
	//	return [NSString stringWithFormat: @"%@ <0x%lx> (%@: %lu points)", NSStringFromClass([self class]), (NSUInteger) self, self.name, self.nDataPoints];
	NSMutableString *extra = [NSMutableString stringWithFormat: @""];
	[extra appendString: self.isComplex ? @"complex variable with" : @"real variable with"];
	[extra appendString: self.isRealPartZero ? @" zero real part" : @" nonzero real part"];
	[extra appendString: self.isImaginaryPartZero ? @", zero imaginary part" : @", nonzero imaginary part"];
	[extra appendString: self.isHermitian ? @" and has hermitian symmetry." : @"."];
	
	return [NSString stringWithFormat: @"%@ <0x%lx> (%@: %lu points) %@\n%@", NSStringFromClass([self class]), (NSUInteger) self, self.name, self.nDataPoints, extra, [self matrixDescription]];
}

@end

@implementation GLMutableNetCDFVariable

- (BOOL) isMutable {
	return YES;
}

- (void) concatenateWithVariable: (GLVariable *) otherVariable alongDimensionAtIndex: (NSUInteger) mutableDimensionIndex
{
	GLVariableOperation *operation = [[GLNetCDFConcatenationOperation alloc] initWithFirstOperand: self secondOperand: otherVariable alongDimensionAtIndex: mutableDimensionIndex];
	[self addOperation: operation];
	[self.equation solveForVariable: self waitUntilFinished: YES];
}

- (void) concatenateWithLowerDimensionalVariable: (GLVariable *) otherVariable alongDimensionAtIndex: (NSUInteger) mutableDimensionIndex toIndex: (NSUInteger) pointIndex;
{
	GLVariableOperation *operation = [[GLNetCDFConcatenationOperation alloc] initWithFirstOperand: self lowerDimensionalSecondOperand: otherVariable alongDimensionAtIndex:mutableDimensionIndex index:pointIndex];
	[self addOperation: operation];
	[self.equation solveForVariable: self waitUntilFinished: YES];
}

/************************************************/
/*		NSCopying								*/
/************************************************/

#pragma mark -
#pragma mark NSCopying
#pragma mark

- (id)copyWithZone:(NSZone *)zone
{
	GLNetCDFVariable *copy = [GLNetCDFVariable variableWithVariable: self];
	copy.file = self.file;
	copy.variableID = self.variableID;
	copy.imagpVariableID = self.imagpVariableID;
	copy.uniqueID = self.uniqueID;
	[copy setupDependency];
	
	return copy;
}

@end
