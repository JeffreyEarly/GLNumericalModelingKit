//
//  GLNetCDFFetchDataOperation.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 1/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import "GLNetCDFFetchDataOperation.h"

@implementation GLNetCDFFetchDataOperation

- (id) initWithNetCDFVariable: (GLNetCDFVariable *) variable indexRange: (NSArray *) ranges flatten: (BOOL) aFlag
{
	return [self initWithResult: nil netCDFVariable:variable indexRange:ranges flatten:aFlag];
}


- (id) initWithResult: (GLVariable *) aResult netCDFVariable: (GLNetCDFVariable *) variable indexRange: (NSArray *) ranges flatten: (BOOL) aFlag
{
	self.shouldFlatten = aFlag;
	
	// Sanity check first
	if ( variable.dimensions.count != ranges.count) return nil;
	
	// What do the reduced dimensions look like?
	NSMutableArray *newDimensions = [[NSMutableArray alloc] init];
	BOOL invalidRange = NO;
	for ( GLDimension *dim in variable.dimensions ) {
		NSRange range = [[ranges objectAtIndex: [variable.dimensions indexOfObject: dim]] rangeValue];
		GLDimension *aNewDimension = [dim subdimensionWithRange: range];
		if (aNewDimension) {
			if ( !(self.shouldFlatten && aNewDimension.nPoints == 1) ) {
				[newDimensions addObject: aNewDimension];
			}
		} else {
			invalidRange = YES;
			NSLog(@"GLNetCDFFetchDataOperation reports index out of bounds.");
		}
	}
	if (invalidRange) return nil;
	
	if (( self = [super init] ))
	{
		if (aResult) {
			self.result = aResult;
		} else {
			self.result = [[GLVariable alloc] initVariableOfType: variable.dataFormat withDimensions: newDimensions forEquation: variable.equation];
		}
		self.file = variable.file;
		self.isComplex = variable.isComplex;
		self.variableID = variable.variableID;
		self.imagpVariableID = variable.imagpVariableID;
		self.theRanges = ranges;
		
		[self setupDependencies];
	}
	
    return self;
}

@synthesize result;
@synthesize file;
@synthesize isComplex;
@synthesize variableID;
@synthesize imagpVariableID;

- (void) setupDependencies
{
	[self.result addOperation: self];
}

- (void) tearDownDependencies
{	
	[self.result removeOperation: self];
	self.file = nil;
}

@synthesize theRanges;
@synthesize shouldFlatten;

- (void) main
{
	if (self.isComplex) {
		GLSplitComplex split = self.result.splitComplex;
		if ( sizeof(GLFloat) == sizeof(double) ) {
			[self.file readDoubleVariableWithID: self.variableID intoBuffer: split.realp indexRange: self.theRanges];
			[self.file readDoubleVariableWithID: self.imagpVariableID intoBuffer: split.imagp indexRange: self.theRanges];
		} else {
			[self.file readFloatVariableWithID: self.variableID intoBuffer: split.realp indexRange: self.theRanges];
			[self.file readFloatVariableWithID: self.imagpVariableID intoBuffer: split.imagp indexRange: self.theRanges];
		}
	} else {
		if ( sizeof(GLFloat) == sizeof(double) ) {
			[self.file readDoubleVariableWithID: self.variableID intoData: self.result.data indexRange: self.theRanges];
		} else {
			[self.file readFloatVariableWithID: self.variableID intoData: self.result.data indexRange: self.theRanges];
		}
	}
	
	[self tearDownDependencies];
}


@end
