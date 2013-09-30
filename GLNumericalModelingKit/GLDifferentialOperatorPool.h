//
//  GLDifferentialOperatorPool.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 4/25/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLDifferentialOperator.h>

/************************************************/
/*		GLDifferentialOperatorPool				*/
/************************************************/

#pragma mark -
#pragma mark GLDifferentialOperatorPool
#pragma mark

// A differential operator pool is a collection of differential operators appropriate for some dimensions.
// Subclasses of GLDifferentialOperatorPool implement standard operators (like d/dx) for different
// differentiation methods (such as spectral and finite difference).
//
// Note that the pool may be valid for more than one type of dimensions.

@class GLEquation;
@interface GLDifferentialOperatorPool : NSObject

- (id) initWithDifferentiationDimensions: (NSArray *) dDims transformDimensions: (NSArray *) tDims forEquation: (GLEquation *) equation;

// The dimensions that are being differentiated.
@property(strong, readonly) NSArray *differentiationDimensions;
@property(strong, readonly) NSArray *differentiationBasis;

// The dimensions that variables must be transformed to in order to differentiate.
@property(strong, readonly) NSArray *transformedDimensions;
@property(strong, readonly) NSArray *transformedBasis;

@property(strong, readonly) GLEquation *equation;

- (BOOL) canOperateOnVariable: (GLVariable *) aVariable;

// Returns the differential operator of the given name.
- (GLDifferentialOperator *) differentialOperatorWithName: (NSString *) opName;

// Add an operator to the pool so that it can be used later.
- (void) setDifferentialOperator: (GLDifferentialOperator *) diffOp forName: (NSString *) name;


- (GLVariable *) differentiationMatrixWithName: (NSString *) opName;
- (void) setDifferentiationMatrix: (GLVariable *) diffVar withFinalDimensions: (NSArray *) finalDims forName: (NSString *) name;

// Destroy the entire pool.
- (void) destroy;

@end