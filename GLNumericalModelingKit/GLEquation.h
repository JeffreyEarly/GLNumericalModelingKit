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
#import <GLNumericalModelingKit/GLIntegrationOperations.h>

@class GLDimension, GLVariable;
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
// Furthermore, the inspector could also possibly trace back in time, digging into the netcdf file.



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
- (void) solveForVariable: (GLVariable *) aVariable;

- (void) solveForVariable: (GLVariable *) aVariable waitUntilFinished: (BOOL) shouldWait;

- (void) solveForOperation: (GLVariableOperation *) anOperation waitUntilFinished: (BOOL) shouldWait;

- (void) waitUntilAllOperationsAreFinished;

/************************************************/
/*		Differentiation							*/
/************************************************/

#pragma mark -
#pragma mark Differentiation
#pragma mark

// Returns the default set of differential operators that will be used by the variables.
// The value set here will be the default, but can be overridden on a case by case basis
// with individual variables.

// Depends on whether or not the variable is real or complex, and the variables dimensions.
- (id) defaultDifferentialOperatorPoolForVariable: (GLVariable *) aVariable;

// Set the default differentiation basis for different orders (the size of the array should match the order).
// If none is specified for a given order, one will be created for the basis of order 1, or if that doesn't exist
// it will use the kGLExponentialBasis.
- (void) setDefaultDifferentiationBasis:(NSArray *) aBasis forOrder: (NSUInteger) order;
- (NSArray *) defaultDifferentiationBasisForOrder: (NSUInteger) order;

/************************************************/
/*		Integration / Time Stepping				*/
/************************************************/

#pragma mark -
#pragma mark Integration / Time Stepping
#pragma mark

// Model equation: dy/dt=f

// In order to time step forward, you need to be able to compute the value of f, given y.
// The FfromY block type takes y as an argument and expects you to return f.
//typedef GLVariable * (^FfromY)(GLVariable *);
//typedef NSArray * (^FfromYVector)(NSArray *);

// Given dy/dt=f at some t(n), this computes y at t(n+1)
// Inputs of y and f (at t=n), a time step delta, and a function to compute f given y.
// Outputs a variable at t=n+delta.
- (GLVariable *) rungeKuttaAdvanceY: (GLVariable *) y withF: (GLVariable *) f stepSize: (GLFloat) deltaT fFromY: (FfromY) aBlock;

// Test op.
- (GLVariable *) rungeKuttaAdvanceY: (GLVariable *) y stepSize: (GLFloat) deltaT fFromY: (FfromY) aBlock;
- (NSArray *) rungeKuttaAdvanceYVector: (NSArray *) y stepSize: (GLFloat) deltaT fFromY: (FfromYVector) aBlock;


@end



