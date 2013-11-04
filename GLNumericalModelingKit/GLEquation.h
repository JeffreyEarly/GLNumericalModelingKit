//
//  GLEquation.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 4/15/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>

#import <GLNumericalModelingKit/Precision.h>
#import <GLNumericalModelingKit/GLVariableOperations.h>

@class GLDimension, GLVariable, GLTensor;
@interface GLEquation : NSObject

// An equation can be thought of as the "controller" of the variables, in the MVC paradigm.

// We can define 'dynamical' variables associated with the equation. These are the things that we keep track of.
// The equation then needs a way

// Key issue:	1. Often we operate with variables that are nondimensionalized, but we want to dimensionalize them before we write them out.
//				2. Often we have some derived variables, like relative vorticity, that we might want to write out.
// In theory these are both the same, they are variables dependent on the ones being computed, but don't need to be computed each time step.
//
// Okay, so, let's call these "Derived Variables".
//
// Why do we want the controller to know which variables are available? Because this allows simple inspection and inquiry about what can be displayed or recorded.
// Furthermore, the inspector could also possibly trace back in time, digging into the netcdf file


/************************************************/
/*		Dimensions								*/
/************************************************/

#pragma mark -
#pragma mark Dimensions
#pragma mark

/************************************************/
/*		Variables								*/
/************************************************/

#pragma mark -
#pragma mark Variables
#pragma mark

// Perform the computation. You don't need (or want) to do this after every operation, but instead
// call this when you reach a choke point. A choke point would be a place where you actually finally
// need the data from the variable, or a place where you *know* you won't benefit from parallization.
- (void) solveForVariable: (GLTensor *) aVariable;

- (void) solveForVariable: (GLTensor *) aVariable waitUntilFinished: (BOOL) shouldWait;

- (void) solveForOperation: (GLVariableOperation *) anOperation waitUntilFinished: (BOOL) shouldWait;

- (void) waitUntilAllOperationsAreFinished;

/************************************************/
/*		Differentiation							*/
/************************************************/

#pragma mark -
#pragma mark Differentiation
#pragma mark

/**  Returns the linear transformation of the given name, if it can be found, nil otherwise.
 @discussion This interface simply serves as storage for the GLLinearTransform class methods, to prevent duplicate transform generation.
 @param name The name of the linear transformation.
 @param dimensions The dimensions that we are transforming from.
 @returns A GLLinearTransform with fromDimensions that match dimensions, and toDimensions that may be different.
 */
- (GLLinearTransform *) linearTransformWithName: (NSString *) name forDimensions: (NSArray *) dimensions;

/**  Add an operator to the pool so that it can be used later.
 @discussion This interface simply serves as storage for the GLLinearTransform class methods, to prevent duplicate transform generation.
 @param diffOp The linear transformation to be saved.
 @param name The name of the linear transformation.
 */
- (void) setLinearTransform: (GLLinearTransform *) diffOp withName: (NSString *) name;

@end



