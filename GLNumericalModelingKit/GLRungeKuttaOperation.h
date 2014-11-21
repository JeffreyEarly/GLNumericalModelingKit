//
//  GLRungeKuttaOperation.h
//  GLNumericalModelingKitCopy
//
//  Created by Jeffrey J. Early on 10/26/12.
//
//

#import <GLNumericalModelingKit/GLVariableOperations.h>
#import <GLNumericalModelingKit/GLScalar.h>
#import <GLNumericalModelingKit/GLNetCDFFile.h>

// You are computing dy/dt=f(t,y). This block is how you specify f, given the current time (t)
// and the current point (y).
// Here the first argument is *scalar* variable t (time).
typedef NSArray * (^FfromTYVector)(GLScalar *, NSArray *);

@class GLRungeKuttaOperation;
typedef NSDictionary * (^OutputFromInputVector)(GLScalar *, NSArray *, GLRungeKuttaOperation* rkint);

/************************************************/
/*		GLRungeKuttaOperation					*/
/************************************************/

#pragma mark -
#pragma mark GLRungeKuttaOperation
#pragma mark

@interface GLRungeKuttaOperation : GLVariableOperation

// Classical 4th order Runge-Kutta time stepping with fixed time step.
+ (instancetype) rungeKutta4AdvanceY: (NSArray *) y stepSize: (GLFloat) deltaT fFromTY: (FfromTYVector) aBlock;

// Cash-Karp 5th order Runge-Kutta time stepping with fixed time step.
+ (instancetype) rungeKutta5AdvanceY: (NSArray *) y stepSize: (GLFloat) deltaT fFromTY: (FfromTYVector) aBlock;


// Steps forward one time-step and returns the new value of y at that point.
- (NSArray *) stepForward;

// Steps forward until time is reached or exceeded and returns the new value of y at that point.
// No attempt is made to shorten the time step to hit the time exactly.
- (NSArray *) stepForwardToTime: (GLFloat) time;

// Creates an operation dependent on time, that returns the integrated value.
// What happens when we request something too far back in time? This will happen if this is embedded in an another RK iterator.
//- (NSArray *) valueAtTime: (GLScalar *) time;

// Integrates from the initial point to the end point of the dimension.
// Note that this function forces variables to be solved. Very different than the usual paradigm.
- (NSArray *) integrateAlongDimension: (GLDimension *) tDim;

// Returns the value at the last time point after it's finished integrating and writing ot file.
- (void) integrateAlongDimension: (GLDimension *) tDim0 toFile: (GLNetCDFFile *) file withTimeScale: (GLFloat) timeScale variables: (OutputFromInputVector) aBlock;

// The size of the input vector y.
@property(readonly)	NSUInteger nInputs;

// The step size is fixed upon initialization.
@property(readonly) GLFloat stepSize;

// Size of the previous time step
@property(readonly) GLFloat lastStepSize;

// The current time is automatically updated each time the operation is stepped forward.
@property(readonly) GLFloat currentTime;

// The initial time defaults to zero.
@property(nonatomic) GLFloat initialTime;

// Total number of time steps taken.
@property(nonatomic) NSUInteger totalIterations;

// If yes, the simulation will stop iterating if the max value is NaN.
@property(nonatomic) BOOL exitOnBlowUp;

// If you set this value only if you've modified y in some way.
@property(strong, nonatomic) NSArray *currentY;

// Sets new current Y values, and also resets the current and initial time. This was already done at initialization,
// so this only needs to be done if you've modified y in some way.
//- (void) setY: (NSArray *) y atTime: (GLFloat) time;

@end


/************************************************/
/*		GLAdaptiveRungeKuttaOperation			*/
/************************************************/

#pragma mark -
#pragma mark GLAdaptiveRungeKuttaOperation
#pragma mark

@interface GLAdaptiveRungeKuttaOperation : GLRungeKuttaOperation


+ (instancetype) rungeKutta23AdvanceY: (NSArray *) y stepSize: (GLFloat) deltaT fFromTY: (FfromTYVector) fFromY;

// Cash-Karp 4th order Runge-Kutta time stepping with *adaptive* time step.
+ (instancetype) rungeKutta45AdvanceY: (NSArray *) y stepSize: (GLFloat) deltaT fFromTY: (FfromTYVector) fFromY;

// The integrator will attempt to satisfy,
// |y_err|/max( relTolerance*|y|, absoluteTolerance ) <= 1
// Relative tolerance reflects the relative accuracy of each of the output value.
// Absolute tolerance essentially defines what you consider zero. For single precision floating
// point, this should be about 1e-7 times a normal value.
// These are arrays of NSNumbers
@property(strong, readwrite, nonatomic) NSArray *relativeTolerance;
@property(strong, readwrite, nonatomic) NSArray *absoluteTolerance;

@end
