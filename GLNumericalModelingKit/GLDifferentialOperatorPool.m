//
//  GLDifferentialOperatorPool.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 4/25/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import "GLEquation.h"
#import "GLDifferentialOperatorPool.h"

@interface GLDifferentialOperatorPool ()
@property(strong) NSMutableDictionary *operatorDictionary;
@property(strong, readwrite) GLEquation *equation;
@end

@implementation GLDifferentialOperatorPool

- (id) initWithDifferentiationDimensions: (NSArray *) dDims transformDimensions: (NSArray *) tDims forEquation: (GLEquation *) equation
{
    if ((self=[super init]))
    {
        _operatorDictionary = [NSMutableDictionary dictionary];
		_differentiationDimensions = dDims;
		_transformedDimensions = tDims;
		_equation = equation;
		
		NSMutableArray *array = [NSMutableArray arrayWithCapacity: tDims.count];
		for (GLDimension *dim in tDims) {
			[array addObject: @(dim.basisFunction)];
		}
		
		_transformedBasis = array;
		
		array = [NSMutableArray arrayWithCapacity: dDims.count];
		for (GLDimension *dim in dDims) {
			[array addObject: @(dim.basisFunction)];
		}
		
		_differentiationBasis = array;
    }
    return self;
}

- (BOOL) canOperateOnVariable: (GLVariable *) aVariable
{
	NSArray *transformedDimensions = [aVariable dimensionsTransformedToBasis: self.transformedBasis];
	
	return [self.transformedDimensions isEqualToArray: transformedDimensions];
}

// Request an empty (but not necessarily zeroed) data object of a particular size
- (GLDifferentialOperator *) differentialOperatorWithName: (NSString *) opName;
{
    return [self.operatorDictionary objectForKey: opName];
}

- (void) setDifferentialOperator: (GLDifferentialOperator *) diffOp forName: (NSString *) name
{
	[(GLVariable *) diffOp solve];
	
	BOOL containsNan = NO;
	for (NSUInteger i=0; i<diffOp.nDataElements; i++) {
		if ( !isfinite(diffOp.pointerValue[i])) containsNan = YES;
	}
	
	if ( containsNan ) {
		NSLog(@"Warning! Your differential operator \"%@\" contains at least one non-finite value. This is likely not what you intended.", name);
	}
	
    [(GLVariable *) diffOp setName: name];
	
	// We don't check to see if the dimensions are appropriate, that's the responsibility of the user.
	[self.operatorDictionary setObject: diffOp forKey: name];
}

- (void) destroy
{
    self.operatorDictionary = [NSMutableDictionary dictionary];
}


@end
