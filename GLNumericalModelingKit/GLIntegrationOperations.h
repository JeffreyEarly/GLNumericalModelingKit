//
//  GLIntegrationOperations.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 4/19/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import <GLNumericalModelingKit/GLVariableOperations.h>

// In order to time step forward, you need to be able to compute the value of f, given y.
// The FfromY block type takes y as an argument and expects you to return f.
typedef GLVariable * (^FfromY)(GLVariable *);
typedef NSArray * (^FfromYVector)(NSArray *);

// Here the first argument is *scalar* variable t (time).
typedef GLVariable * (^FfromTY)(GLVariable *, GLVariable *);
typedef NSArray * (^FfromTYVector)(GLVariable *, NSArray *);

// Return an appropriate time step given y. The variable returned must be a scalar variable.
typedef GLVariable * (^tFromY)(GLVariable *);
typedef GLVariable * (^tFromYVector)(NSArray *);

/************************************************/
/*		GLIntegrationOperation					*/
/************************************************/

#pragma mark -
#pragma mark GLIntegrationOperation
#pragma mark

@interface GLIntegrationOperation : GLUnaryOperation

// Classical 4th order Runge-Kutta time stepping with fixed time step.
+ (id) rungeKutta4AdvanceY: (GLVariable *) y stepSize: (double) deltaT fFromY: (FfromY) aBlock;

// Cash-Karp 5th order Runge-Kutta time stepping with fixed time step.
+ (id) rungeKutta5AdvanceY: (GLVariable *) y stepSize: (double) deltaT fFromY: (FfromY) aBlock;

// Classical 4th order Runge-Kutta time stepping with fixed time step; also take t as an argument.
+ (id) rungeKutta4AdvanceY: (GLVariable *) y stepSize: (double) deltaT fFromTY: (FfromTY) fFromTY;

// The step size is fixed upon initialization.
@property(readonly) GLFloat stepSize;

// Size of the previous time step
@property(readonly) GLFloat lastStepSize;

// The current time is automatically updated each time the operation is stepped forward.
@property(readonly) double currentTime;

// The initial time defaults to zero.
@property(nonatomic) double initialTime;

// If yes, the simulation will stop iterating if the max value is NaN.
@property(nonatomic) BOOL exitOnBlowUp;

// Steps forward one time-step.
- (GLVariable *) stepForward: (GLVariable *) y;

// Steps foward until time is reached or exceeded.
// No attempt is made to shorten the time step to hit the time exactly.
- (GLVariable *) stepForward: (GLVariable *) y toTime: (double) time;

@end

/************************************************/
/*		GLAdaptiveIntegrationOperation			*/
/************************************************/

#pragma mark -
#pragma mark GLAdaptiveIntegrationOperation
#pragma mark

@interface GLAdaptiveIntegrationOperation : GLIntegrationOperation

+ (id) rungeKutta23AdvanceY: (GLVariable *) y stepSize: (GLFloat) deltaT fFromY: (FfromY) fFromY;

// Cash-Karp 4th order Runge-Kutta time stepping with *adaptive* time step.
+ (id) rungeKutta45AdvanceY: (GLVariable *) y stepSize: (GLFloat) deltaT fFromY: (FfromY) aBlock;

// The time step gets computed at the start of the time step operation.
+ (id) rungeKutta4AdvanceY:(GLVariable *)y dynamicStepSize: (tFromY) deltaTfromY fFromY:(FfromY)aBlock;

// The integrator will attempt to satisfy,
// |y_err|/max( relTolerance*|y|, absoluteTolerance ) <= 1
// Relative tolerance reflects the relative accuracy of each of the output value.
// Absolute tolerance essentially defines what you consider zero. For single precision floating
// point, this should be about 1e-7 times a normal value.
@property GLFloat relativeTolerance;
@property GLFloat absoluteTolerance;

@end

/************************************************/
/*		GLVectorIntegrationOperation			*/
/************************************************/

#pragma mark -
#pragma mark GLVectorIntegrationOperation
#pragma mark

@interface GLVectorIntegrationOperation : GLUnaryVectorOperation

// Classical 4th order Runge-Kutta time stepping with fixed time step.
+ (id) rungeKutta4AdvanceY: (NSArray *) y stepSize: (GLFloat) deltaT fFromY: (FfromYVector) aBlock;

// Cash-Karp 5th order Runge-Kutta time stepping with fixed time step.
+ (id) rungeKutta5AdvanceY: (NSArray *) y stepSize: (GLFloat) deltaT fFromY: (FfromYVector) fFromY;

// Classical 4th order Runge-Kutta time stepping with fixed time step.
+ (id) rungeKutta4AdvanceY: (NSArray *) y stepSize: (GLFloat) deltaT fFromTY: (FfromTYVector) fFromY;

// The step size is fixed upon initialization.
@property(readonly) GLFloat stepSize;

// Size of the previous time step
@property(readonly) GLFloat lastStepSize;

// The current time is automatically updated each time the operation is stepped forward.
@property(readonly) GLFloat currentTime;

// The initial time defaults to zero.
@property(nonatomic) GLFloat initialTime;

// If yes, the simulation will stop iterating if the max value is NaN.
@property(nonatomic) BOOL exitOnBlowUp;

// Steps forward one time-step.
- (NSArray *) stepForward: (NSArray *) y;

// Steps foward until time is reached or exceeded.
// No attempt is made to shorten the time step to hit the time exactly.
- (NSArray *) stepForward: (NSArray *) y toTime: (GLFloat) time;

@end

/************************************************/
/*		GLAdaptiveVectorIntegrationOperation	*/
/************************************************/

#pragma mark -
#pragma mark GLAdaptiveVectorIntegrationOperation
#pragma mark

@interface GLAdaptiveVectorIntegrationOperation : GLVectorIntegrationOperation

// Cash-Karp 4th order Runge-Kutta time stepping with *adaptive* time step.
+ (id) rungeKutta45AdvanceY: (NSArray *) y stepSize: (GLFloat) deltaT fFromY: (FfromYVector) aBlock;

// The integrator will attempt to satisfy,
// |y_err|/max( relTolerance*|y|, absoluteTolerance ) <= 1
// Relative tolerance reflects the relative accuracy of each of the output value.
// Absolute tolerance essentially defines what you consider zero. For single precision floating
// point, this should be about 1e-7 times a normal value.
// These are arrays of NSNumbers
@property NSArray *relativeTolerance;
@property NSArray *absoluteTolerance;

@end
