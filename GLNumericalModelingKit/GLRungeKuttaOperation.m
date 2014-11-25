//
//  GLRungeKuttaOperation.m
//  GLNumericalModelingKitCopy
//
//  Created by Jeffrey J. Early on 10/26/12.
//
//

#import "GLRungeKuttaOperation.h"

#import "GLMemoryPool.h"
#import "GLUnaryOperations.h"
#import "GLVectorVectorOperations.h"
#import "GLOperationVisualizer.h"
#import "GLOperationOptimizer.h"

@interface GLElementErrorOperation : GLVariableOperation
- (id) initWithFirstOperand: (GLVariable *) yerrVar secondOperand: (GLVariable *) yVar relativeError: (GLScalar *) relErr absoluteError: (GLScalar *) absErr;
@end

@implementation GLElementErrorOperation

// *element-wise*
// |y_err|/max( relTolerance*|y|, absoluteTolerance ) <= 1
- (id) initWithFirstOperand: (GLVariable *) yerrVar secondOperand: (GLVariable *) yVar relativeError: (GLScalar *) relErr absoluteError: (GLScalar *) absErr
{
    GLScalar *errorOut = [[GLScalar alloc] initWithType: kGLRealDataFormat forEquation: yerrVar.equation];
    GLBuffer *aBuffer = [[GLBuffer alloc] initWithLength: yerrVar.dataBytes];
    NSUInteger nDataElements = yerrVar.nDataElements;
    variableOperation operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
        GLFloat *yerr = (GLFloat *) [operandArray[0] bytes];
        GLFloat *y = (GLFloat *) [operandArray[1] bytes];
        GLFloat *relativeError = (GLFloat *) [operandArray[2] bytes];
        GLFloat *absoluteError = (GLFloat *) [operandArray[3] bytes];
        GLFloat *max = (GLFloat *) [resultArray[0] bytes];
        GLFloat *buf = (GLFloat *) [bufferArray[0] bytes];
        
        
        vGL_vsmul( y, 1, relativeError, buf, 1, nDataElements); // buf = relTolerance*y
        vGL_vabs( buf, 1, buf, 1, nDataElements ); // buf = relTolerance*|y|
        vGL_vthr( buf, 1, absoluteError, buf, 1, nDataElements); // buf = max( relTolerance*|y|, absoluteTolerance )
        vGL_vdiv( buf, 1, yerr, 1, buf, 1, nDataElements); // buf = y_err/max( relTolerance*|y|, absoluteTolerance )
        vGL_vabs( buf, 1, buf, 1, nDataElements ); // |y_err|/max( relTolerance*|y|, absoluteTolerance )
        vGL_maxv( buf, 1, max, nDataElements);
    };
    
    if (( self = [super initWithResult: @[errorOut] operand: @[yerrVar, yVar, relErr, absErr] buffers: @[aBuffer] operation: operation] )) {
		self.graphvisDescription = @"error";
    }
    
    return self;
}

@end

/************************************************/
/*		GLRungeKuttaOperation					*/
/************************************************/

#pragma mark -
#pragma mark GLRungeKuttaOperation
#pragma mark

// This type is designed to take an array "y", and then write "f" to the appropriate data chunks that are fixed.
typedef void (^stagePrepOperation)(NSArray *, NSArray *);

BOOL isFinite( NSData *data, NSUInteger nDataElements )
{
	GLFloat *pointerValue = (GLFloat *) data.bytes;
	for (NSUInteger i=0; i<nDataElements; i++) {
		if ( !isfinite(pointerValue[i])) return NO;
	}
	return YES;
}

BOOL isZero( NSNumber *a )
{
	return fabs( a.doubleValue ) < 1e-7;
}

BOOL isOne( NSNumber *a )
{
	return fabs( a.doubleValue - 1.0) < 1e-7;
}

@interface GLRungeKuttaOperation ()
- (id) initMethodWithCoefficientsA: (NSArray *) a b: (NSArray *) b c: (NSArray *) c y: (NSArray *) y stepSize: (GLFloat) deltaT fFromTY: (FfromTYVector) aBlock;

@property(strong) NSMutableArray *dataBuffers;
@property GLFloat stepSize;
@property NSUInteger nInputs;

@property(nonatomic) GLFloat requestedTime;

@property(strong) GLScalar *stepSizeVariable;
@property(strong) GLScalar *lastStepSizeVariable;
@property(strong) GLScalar *currentTimeVariable;
@property(strong) GLScalar *requestedTimeVariable;

@property(strong) NSMutableData *stepSizeData;
@property(strong) NSMutableData *lastStepSizeData;
@property(strong) NSMutableData *currentTimeData;
@property(strong) NSMutableData *requestedTimeData;

// Whether or not to hold onto the previous value of y.
@property BOOL shouldRetainPreviousY;
@property(strong) NSMutableArray *previousY;
@property(strong) NSMutableArray *previousYData;

// Indicates whether or not the first stage is the same as the last.
@property BOOL isFSAL;
@property(strong) NSMutableArray *firstStageVariables;
@property(strong) NSMutableArray *lastStageVariables;
@property(strong) NSMutableArray *firstStageData;
@property(strong) NSMutableArray *lastStageData;
// The operation computes f and stores the results in lastStageData.
// This will need to be called anytime that y is reset.
@property BOOL stageIsPrepped;
@property(copy) stagePrepOperation prepStageOperation;
@property(strong) NSMutableArray *prepStageDataBuffers;

@property(strong) NSArray *initialY;

@end

@implementation GLRungeKuttaOperation

+ (id) rungeKutta4AdvanceY: (NSArray *) y stepSize: (GLFloat) deltaT fFromTY: (FfromTYVector) aBlock
{
	NSUInteger numStages = 4;
	
	NSMutableArray *a = [NSMutableArray arrayWithCapacity:numStages];
	a[0] = @(0.0);
	a[1] = @(1./2.);
	a[2] = @(1./2.);
	a[3] = @(1.);

	
	NSMutableArray *b = [NSMutableArray arrayWithCapacity:numStages];
	for (NSUInteger i=0; i<numStages; i++) {
		b[i] = [NSMutableArray arrayWithCapacity: i+1];
	}
	b[0][0] = @(0.0);
	b[1][0] = @(1./2.); b[1][1] = @(0.0);
	b[2][0] = @(0.0); b[2][1] = @(1./2.); b[2][2] = @(0.0);
	b[3][0] = @(0.0); b[3][1] = @(0.0); b[3][2] = @(1.0); b[3][3] = @(0.0);

	
	NSMutableArray *c = [NSMutableArray arrayWithCapacity:numStages];
	c[0] = @(1./6.);
	c[1] = @(1./3.);
	c[2] = @(1./3.);
	c[3] = @(1./6.);
	
	GLRungeKuttaOperation *rk = [[GLRungeKuttaOperation alloc] initMethodWithCoefficientsA: a b: b c: c y: y stepSize: deltaT fFromTY: aBlock];
	
	return rk;
}

+ (id) rungeKutta5AdvanceY: (NSArray *) y stepSize: (GLFloat) deltaT fFromTY: (FfromTYVector) aBlock
{
	NSUInteger numStages = 6;
	
	NSMutableArray *a = [NSMutableArray arrayWithCapacity:numStages];
	a[0] = @(0.0);
	a[1] = @(1./5.);
	a[2] = @(3./10.);
	a[3] = @(3./5.);
	a[4] = @(1.);
	a[5] = @(7./8.);
	
	NSMutableArray *b = [NSMutableArray arrayWithCapacity:numStages];
	for (NSUInteger i=0; i<numStages; i++) {
		b[i] = [NSMutableArray arrayWithCapacity: i+1];
	}
	b[0][0] = @(0.0);
	b[1][0] = @(1./5.); b[1][1] = @(0.0);
	b[2][0] = @(3./40.); b[2][1] = @(9./40.); b[2][2] = @(0.0);
	b[3][0] = @(3./10.); b[3][1] = @(-9./10); b[3][2] = @(6./5.); b[3][3] = @(0.0);
	b[4][0] = @(-11./54.); b[4][1] = @(5./2.); b[4][2] = @(-70./27.); b[4][3] = @(35./27.); b[4][4] = @(0.0);
	b[5][0] = @(1631./55296.); b[5][1] = @(175./512.); b[5][2] = @(575./13824.); b[5][3] = @(44275./110592.); b[5][4] = @(253./4096.); b[5][5] = @(0.0);
	
	NSMutableArray *c = [NSMutableArray arrayWithCapacity:numStages];
	c[0] = @(37./378.);
	c[1] = @(0.0);
	c[2] = @(250./621.);
	c[3] = @(125./594.);
	c[4] = @(0.0);
	c[5] = @(512./1771.);
	
	GLRungeKuttaOperation *rk = [[GLRungeKuttaOperation alloc] initMethodWithCoefficientsA: a b: b c: c y: y stepSize: deltaT fFromTY: aBlock];
	
	return rk;
}

- (id) initWithResult:(NSArray *)result operand:(NSArray *)operand
{
	if ((self=[super initWithResult: result operand: operand])) {
		GLEquation *equation =[operand[0] equation];
		
		self.currentTimeVariable = [[GLScalar alloc] initWithType: kGLRealDataFormat forEquation: equation];
		self.stepSizeVariable = [[GLScalar alloc] initWithType: kGLRealDataFormat forEquation: equation];
		self.lastStepSizeVariable = [[GLScalar alloc] initWithType: kGLRealDataFormat forEquation: equation];
		self.requestedTimeVariable = [[GLScalar alloc] initWithType: kGLRealDataFormat forEquation: equation];
		
		self.currentTimeData = self.currentTimeVariable.data;
		self.stepSizeData = self.stepSizeVariable.data;
		self.lastStepSizeData = self.lastStepSizeVariable.data;
		self.requestedTimeData = self.requestedTimeVariable.data;
		
		self.nInputs = operand.count;
		
		self.exitOnBlowUp = YES;
		
		self.previousY = [NSMutableArray array];
		self.previousYData = [NSMutableArray array];
		for (GLVariable *variable in operand ) {
			[self.previousY addObject: [GLVariable variableWithPrototype: variable]];
			[self.previousYData addObject: [[self.previousY lastObject] data]];
		}
	}
	
	return self;
}

// This assumes *fixed* step size, and no error correction.
- (id) initMethodWithCoefficientsA: (NSArray *) a b: (NSArray *) b c: (NSArray *) c y: (NSArray *) y stepSize: (GLFloat) deltaT fFromTY: (FfromTYVector) fFromY
{
	NSUInteger numStages = a.count;
	
	GLScalar *time = [[GLScalar alloc] initWithType: kGLRealDataFormat forEquation: [y[0] equation]];
	NSMutableArray *t = [NSMutableArray arrayWithCapacity: numStages];
	for (NSUInteger i=0; i<a.count; i++) {
		t[i] = [time plus: @([a[i] doubleValue]*deltaT)];
	}
	
	// Store the value at each stage point
	NSMutableArray *yp = [NSMutableArray arrayWithCapacity: numStages];
	yp[0] = y;
	
	// Store the slope at each stage point
	NSMutableArray *f = [NSMutableArray arrayWithCapacity: numStages];
	f[0] = fFromY(t[0], yp[0]);
	
	// Walk through all the stages (i)
	for (NSUInteger i=1; i<numStages; i++) {
		// Compute y at this stage location
		yp[i] = [NSMutableArray arrayWithArray: y];
		for (NSUInteger n=0; n<y.count; n++) {
			for (NSUInteger j=0; j<i; j++) {
				if (!isZero(b[i][j])) {
					yp[i][n] = [yp[i][n] plus:[f[j][n] scalarMultiply: deltaT*[b[i][j] doubleValue]]];
				}
			}
		}
		// Compute the slope at this stage location
		f[i] = fFromY(t[i], yp[i]);
	}
	
	// Now compute the value of y at the very end
	NSMutableArray *yout = [NSMutableArray arrayWithArray: y];
	for (NSUInteger n=0; n<y.count; n++) {
		for (NSUInteger i=0; i<numStages; i++) {
			if (!isZero(c[i])) {
				yout[n] = [yout[n] plus:[f[i][n] scalarMultiply: deltaT*[c[i] doubleValue]]];
			}
		}
	}
	
	NSMutableArray *yfull = [NSMutableArray arrayWithArray: y];
	[yfull addObject: time];
    
	GLOperationOptimizer *optimizer = [[GLOperationOptimizer alloc] initWithTopVariables: yfull bottomVariables: yout];
    GLOperationVisualizer *vizualizer = [[GLOperationVisualizer alloc] initWithTopVariables: yfull bottomVariables: yout];
	
	if ((self=[self initWithResult: optimizer.resultVariablePrototypes operand: y])) {
		
		variableOperation vectorBlock = optimizer.operationBlock;
		NSMutableArray *operandArray = [NSMutableArray arrayWithCapacity: yfull.count];
		NSMutableData *tData = self.currentTimeData;
		self.operation = ^(NSArray *result, NSArray *operand, NSArray *bufferArray) {
			[operandArray addObjectsFromArray: operand];
			[operandArray addObject: tData];
			
			vectorBlock( result, operandArray, bufferArray);
			
			[operandArray removeAllObjects];
		};
		
        // Go ahead and allocate the memory for those buffers.
        self.dataBuffers = [NSMutableArray array];
		for (GLBuffer *aBuffer in optimizer.internalDataBuffers) {
			[self.dataBuffers addObject: [[GLMemoryPool sharedMemoryPool] dataWithLength: aBuffer.numBytes]];
		}
        
		self.initialY = y;
		self.currentY = y;
		self.stepSize = deltaT;
		self.lastStepSize = deltaT;
		self.exitOnBlowUp = YES;
        self.graphvisDescription = vizualizer.graphvisDescription;
	}
	
	return self;
}

- (void) updateCurrentTime
{
	self.currentTime = self.stepSize * ( (GLFloat) self.totalIterations) + self.initialTime;
}

- (NSArray *) stepForward
{
	return [self stepForwardToTime: self.stepSize * (( (GLFloat) self.totalIterations ) + 0.5) + self.initialTime];
}

//- (void) addRequestedTime: (NSData *) newTimeData
//{
//    
//}
//
//- (NSArray *) valueAtTime: (GLScalar *) time
//{
//    // integrator objects keeps track of the last X (4?) requests for values.
//    // The integrator stores the y values that span all 4 of those points. That's a minimum of 2, maybe lots more
//    // This needs to be an operation that depends on the time scalar, and returns the result array.
//    // This operation doesn't depend on the RK operation (in that we don't add it as a parent), but it does hold onto it.
//    NSMutableArray *yout = [[NSMutableArray alloc] initWithCapacity: self.result.count];
//    for (GLVariable *variable in self.result) {
//        [yout addObject: [GLVariable variableWithPrototype: variable]];
//    }
//    
//    GLVariableOperation *op = [[GLVariableOperation alloc] initWithResult: yout operand: @[time]];
//    op.operation = ^(NSArray *result, NSArray *operand, NSArray *bufferArray) {
//        [self addRequestedTime: operand[0]];
//        NSArray *precedingResult = [self storedResultPrecedingTime: operand[0]];
//        NSArray *followingResult = [self storedResultFollowingTime: operand[0]];
//        
//    };
//    return op.result;
//}

- (NSArray *) stepForwardToTime: (GLFloat) time
{
	self.requestedTime = time;
	
	// If we've already exceeded the time, just return the same variable.
	if ( self.currentTime >= time) {
		return self.currentY;
	}
		
	// If we haven't, then take one time step forward...
	[self.currentY makeObjectsPerformSelector: @selector(solve)];
	NSMutableArray *yout = [[NSMutableArray alloc] initWithCapacity: self.result.count];
	NSMutableArray *resultBuffer = [[NSMutableArray alloc] initWithCapacity: self.result.count];
	NSMutableArray *operandBuffer = [[NSMutableArray alloc] initWithCapacity: self.operand.count];
	for (GLVariable *variable in self.result) {
		GLVariable *newVar = [GLVariable variableWithPrototype: variable];
		[yout addObject: newVar];
		[resultBuffer addObject: newVar.data];
	}
	for (GLFunction *variable in self.currentY) {
		[operandBuffer addObject: variable.data];
	}
	
	// If the algorithm isFSAL, then we plan on using the last stage from the previous time step.
	// However, if there was no previous time step, then the stage won't actually be there and we need to compute it!
	if (self.isFSAL && !self.stageIsPrepped) {
		self.prepStageOperation( operandBuffer, self.prepStageDataBuffers );
		self.stageIsPrepped = YES;
	}
	
	// Take one step forward in time
	self.operation( resultBuffer, operandBuffer, self.dataBuffers );
	self.totalIterations = self.totalIterations + 1;
	[self updateCurrentTime];
	
    if (self.currentTime < time)
    {
        // We can't feed the same databuffer to the input and output, because we'll overwrite
        // the input. This could be bad if happened to take too large of a step.
        NSMutableArray *youtOut = [[NSMutableArray alloc] initWithCapacity: self.result.count];
        NSMutableArray *resultBufferOut= [[NSMutableArray alloc] initWithCapacity: self.result.count];
        for (GLVariable *variable in self.result) {
            GLVariable *newVar = [GLVariable variableWithPrototype: variable];
            [youtOut addObject: newVar];
            [resultBufferOut addObject: newVar.data];
        }
        
        while ( self.currentTime < time )
        {
            self.operation( resultBufferOut, resultBuffer, self.dataBuffers );
            self.totalIterations = self.totalIterations + 1;
            [self updateCurrentTime];
            
            NSMutableArray *tmp = resultBuffer;
            resultBuffer = resultBufferOut;
            resultBufferOut = tmp;
            tmp = yout;
            yout = youtOut;
            youtOut = tmp;
        }
    }
	
    self.currentY = yout;
    return yout;

}

- (NSArray *) integrateAlongDimension: (GLDimension *) tDim
{
	self.totalIterations = 0;
	self.currentTime = 0;
	self.currentY = self.initialY;
	
	NSMutableArray *yout = [[NSMutableArray alloc] initWithCapacity: self.result.count];
	for (GLVariable *variable in self.result) {
		GLFunction *newFunction;
		if (variable.rank == 0) {
			GLScalar *scalar = (GLScalar *) variable;
			newFunction = [GLFunction functionOfType: scalar.dataFormat withDimensions: @[tDim] forEquation: scalar.equation];
			[yout addObject: newFunction];
			GLScalar *initialY = self.currentY[[self.result indexOfObject: variable]];
			[initialY solve];
			memcpy(newFunction.pointerValue, initialY.pointerValue, initialY.dataBytes);
		} else if (variable.rank == 1) {
			GLFunction *function = (GLFunction *) variable;
			NSMutableArray *newDimensions = [NSMutableArray arrayWithObject: tDim];
			[newDimensions addObjectsFromArray: function.dimensions];
			newFunction = [[function class] functionOfType: function.dataFormat withDimensions: newDimensions forEquation: function.equation];
			[yout addObject: newFunction];
			GLFunction *initialY = self.currentY[[self.result indexOfObject: variable]];
			[initialY solve];
			memcpy(newFunction.pointerValue, initialY.pointerValue, initialY.dataBytes);
		}  else if (variable.rank == 2) {
			[NSException raise:@"BadFormat" format: @"Cannot add a time dimension to a linear transformation"];
		}
	}
	
	for (NSUInteger iPoint=1; iPoint<tDim.nPoints; iPoint++) {
		@autoreleasepool {
			NSArray *y = [self stepForwardToTime: [tDim valueAtIndex: iPoint]];
			[y makeObjectsPerformSelector: @selector(solve)];
			[[y[0] equation] solveForVariable: y[0] waitUntilFinished: YES];
			for (GLVariable *variable in y) {
				if (variable.rank == 0) {
					GLScalar *scalar = (GLScalar *) variable;
					GLFunction *newFunction = yout[[y indexOfObject: variable]];
					memcpy((void *)newFunction.pointerValue + iPoint*scalar.dataBytes, scalar.pointerValue, scalar.dataBytes);
				} else if (variable.rank == 1) {
					GLFunction *function = (GLFunction *) variable;
					GLFunction *newFunction = yout[[y indexOfObject: variable]];
					memcpy((void *)newFunction.pointerValue + iPoint*function.dataBytes, function.pointerValue, function.dataBytes);
				}
			}
		}
	}
	
	return yout;
}

- (void) integrateAlongDimension: (GLDimension *) tDim0 withTimeScale: (GLFloat) timeScale file: (GLNetCDFFile *) file output: (OutputFromInputVector) aBlock;
{
	GLEquation *equation = [self.currentY.lastObject equation];
	
    if (!tDim0.name) {
        [NSException raise:@"BadFormat" format: @"tDim must have a name because it will be saved to file."];
    }
	GLMutableDimension *tDim = [[GLMutableDimension alloc] initWithPoints: @[@([tDim0 valueAtIndex: 0])]];
    tDim.name = tDim0.name;
	
	// Write the initial values to file
	GLScalar *t0 = [GLScalar scalarWithValue: [tDim valueAtIndex: 0] forEquation: equation];
	NSDictionary *initialY_scaled = aBlock( t0, self.currentY );
	NSMutableDictionary *variableHistories = [NSMutableDictionary dictionary];
	for (NSString *key in initialY_scaled) {
		GLVariable *variable = initialY_scaled[key];
		GLMutableVariable *variableHistory;
		if (variable.rank == 0) {
			GLScalar *scalar = (GLScalar *) variable;
			GLFunction *newFunction = [GLFunction functionOfType: scalar.dataFormat withDimensions: @[] forEquation: scalar.equation];
//            newFunction.name = key;
//            newFunction.units = variable.units;
			variableHistory = [newFunction variableByAddingDimension: tDim];
		} else if (variable.rank == 1) {
			GLFunction *function = (GLFunction *) variable;
			variableHistory = [function variableByAddingDimension: tDim];
		} else if (variable.rank == 2) {
			[NSException raise:@"BadFormat" format: @"Cannot add a time dimension to a linear transformation"];
		}
		variableHistory.name = key;
		variableHistory.units = variable.units;
		variableHistory = [file addVariable: variableHistory];
		variableHistories[key] = variableHistory;
	}
    
    if (self.shouldDisplayProgress) {
        NSLog(@"Logging at time: %f, step size: %f.", timeScale*self.currentTime, self.lastStepSize*timeScale);
    }
	
	for (NSUInteger iPoint=1; iPoint<tDim0.nPoints; iPoint++) {
		@autoreleasepool {
			GLFloat time = [tDim0 valueAtIndex: iPoint];
			NSArray *y = [self stepForwardToTime: time];
            
            if (self.shouldDisplayProgress) {
                NSLog(@"Logging at time: %f, step size: %f.", timeScale*self.currentTime, self.lastStepSize*timeScale);
            }
            
			GLScalar *tn = [GLScalar scalarWithValue: time forEquation: equation];
			NSDictionary *y_scaled =aBlock( tn, y );
			
			[tDim addPoint: @(time*timeScale)];
			
			for (NSString *key in y_scaled) {
				GLVariable *variable = y_scaled[key];
				GLMutableVariable *variableHistory = variableHistories[key];
				[variableHistory concatenateWithLowerDimensionalVariable: variable alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
			}
			
			[file waitUntilAllOperationsAreFinished];
		}
	}	
}

- (void) setCurrentY:(NSArray *)currentY {
	_currentY = currentY;
	self.stageIsPrepped = NO;
}

- (void) setStepSize: (GLFloat)size
{
	GLFloat *a = (GLFloat *) self.stepSizeData.mutableBytes;
	*a = size;
}

- (GLFloat) stepSize
{
	GLFloat *a = (GLFloat *) self.stepSizeData.mutableBytes;
	return *a;
}

- (void) setLastStepSize: (GLFloat)size
{
	GLFloat *a = (GLFloat *) self.lastStepSizeData.mutableBytes;
	*a = size;
}

- (GLFloat) lastStepSize
{
	GLFloat *a = (GLFloat *) self.lastStepSizeData.mutableBytes;
	return *a;
}

- (void) setCurrentTime: (GLFloat)time
{
	GLFloat *a = (GLFloat *) self.currentTimeData.mutableBytes;
	*a = time;
}

- (GLFloat) currentTime
{
	GLFloat *a = (GLFloat *) self.currentTimeData.mutableBytes;
	return *a;
}

- (void) setRequestedTime: (GLFloat)time
{
	GLFloat *a = (GLFloat *) self.requestedTimeData.mutableBytes;
	*a = time;
}

- (GLFloat) requestedTime
{
	GLFloat *a = (GLFloat *) self.requestedTimeData.mutableBytes;
	return *a;
}

@end

/************************************************/
/*		GLAdaptiveRungeKuttaOperation			*/
/************************************************/

#pragma mark -
#pragma mark GLAdaptiveRungeKuttaOperation
#pragma mark

@interface GLAdaptiveRungeKuttaOperation ()

- (id) initMethodWithCoefficientsA: (NSArray *) a b: (NSArray *) b c: (NSArray *) c d: (NSArray *) d y: (NSArray *) y stepSize: (GLFloat) deltaT fsal: (BOOL) isFSAL retainPreviousY: (BOOL) shouldRetainPreviousY  order: (NSUInteger) order fFromTY: (FfromTYVector) fFromY;
- (NSArray *) interpolateAtTime: (GLScalar *) t currentTime: (GLScalar *) tNow lastStepSize: (GLScalar *) h yNow: (NSArray *) yNow  fNow: (NSArray *) fNow p: (NSArray *) polyCoeffs;



// Error variables---one for each input
@property(strong) NSMutableArray *relativeToleranceVariables;
@property(strong) NSMutableArray *absoluteToleranceVariables;
@property(strong) NSMutableArray *errorVariables;

@property(strong) NSMutableArray *relativeToleranceData;
@property(strong) NSMutableArray *absoluteToleranceData;
@property(strong) NSMutableArray *errorData;

// Hold onto all the stages, in case we want to interpolate
@property(strong) NSArray *stageVariables;

@property(copy) variableOperation interpolationOperation;
@property(strong) NSMutableArray *interpolationBuffers;

// Re-declared

@property(nonatomic) GLFloat requestedTime;

@property(strong) GLScalar *stepSizeVariable;
@property(strong) GLScalar *lastStepSizeVariable;
@property(strong) GLScalar *currentTimeVariable;
@property(strong) GLScalar *requestedTimeVariable;

@property(strong) NSMutableData *stepSizeData;
@property(strong) NSMutableData *lastStepSizeData;
@property(strong) NSMutableData *currentTimeData;
@property(strong) NSMutableData *requestedTimeData;

// Whether or not to hold onto the previous value of y.
@property BOOL shouldRetainPreviousY;
@property(strong) NSMutableArray *previousY;
@property(strong) NSMutableArray *previousYData;

// Indicates whether or not the first stage is the same as the last.
// These variables should only EVER be used internally---never handed back to the user.
// The only reason we hold onto the variables
@property BOOL isFSAL;
@property(strong) NSMutableArray *firstStageVariables;
@property(strong) NSMutableArray *lastStageVariables;
@property(strong) NSMutableArray *firstStageData;
@property(strong) NSMutableArray *lastStageData;
// The operation computes f and stores the results in lastStageData.
// This will need to be called anytime that y is reset.
@property BOOL stageIsPrepped;
@property(copy) stagePrepOperation prepStageOperation;

@end

@implementation GLAdaptiveRungeKuttaOperation

+ (id) rungeKutta45AdvanceY: (NSArray *) y stepSize: (GLFloat) deltaT fFromTY: (FfromTYVector) aBlock
{
	NSUInteger numStages = 6;
	NSUInteger order = 4;
	
	NSMutableArray *a = [NSMutableArray arrayWithCapacity:numStages];
	a[0] = @(0.0);
	a[1] = @(1./5.);
	a[2] = @(3./10.);
	a[3] = @(3./5.);
	a[4] = @(1.);
	a[5] = @(7./8.);
	
	NSMutableArray *b = [NSMutableArray arrayWithCapacity:numStages];
	for (NSUInteger i=0; i<numStages; i++) {
		b[i] = [NSMutableArray arrayWithCapacity: i+1];
	}
	b[0][0] = @(0.0);
	b[1][0] = @(1./5.); b[1][1] = @(0.0);
	b[2][0] = @(3./40.); b[2][1] = @(9./40.); b[2][2] = @(0.0);
	b[3][0] = @(3./10.); b[3][1] = @(-9./10); b[3][2] = @(6./5.); b[3][3] = @(0.0);
	b[4][0] = @(-11./54.); b[4][1] = @(5./2.); b[4][2] = @(-70./27.); b[4][3] = @(35./27.); b[4][4] = @(0.0);
	b[5][0] = @(1631./55296.); b[5][1] = @(175./512.); b[5][2] = @(575./13824.); b[5][3] = @(44275./110592.); b[5][4] = @(253./4096.); b[5][5] = @(0.0);
	
	NSMutableArray *c = [NSMutableArray arrayWithCapacity:numStages];
	c[0] = @(37./378.);
	c[1] = @(0.0);
	c[2] = @(250./621.);
	c[3] = @(125./594.);
	c[4] = @(0.0);
	c[5] = @(512./1771.);
	
	NSMutableArray *d = [NSMutableArray arrayWithCapacity:numStages];
	d[0] = @(2825./27648.);
	d[1] = @(0.0);
	d[2] = @(18575./48384.);
	d[3] = @(13525./55296.);
	d[4] = @(277./14336.);
	d[5] = @(1./4.);
	
	GLAdaptiveRungeKuttaOperation *rk = [[GLAdaptiveRungeKuttaOperation alloc] initMethodWithCoefficientsA: a b: b c: c d: d y: y stepSize: deltaT fsal: NO retainPreviousY: NO order: order fFromTY: aBlock];
	for (GLFunction *rt in rk.relativeToleranceVariables) {
		*(rt.pointerValue) = 1e-5;
	}
	for (GLFunction *at in rk.absoluteToleranceVariables) {
		*(at.pointerValue) = 1e-6;
	}
	
	return rk;
}

+ (id) rungeKutta23AdvanceY: (NSArray *) y stepSize: (GLFloat) deltaT fFromTY: (FfromTYVector) aBlock
{
	NSUInteger numStages = 4;
	NSUInteger order = 2;
	
	NSMutableArray *a = [NSMutableArray arrayWithCapacity:numStages];
	a[0] = @(0.0);
	a[1] = @(1./2.);
	a[2] = @(3./4.);
	a[3] = @(1.);
	
	NSMutableArray *b = [NSMutableArray arrayWithCapacity:numStages];
	for (NSUInteger i=0; i<numStages; i++) {
		b[i] = [NSMutableArray arrayWithCapacity: i+1];
	}
	b[0][0] = @(0.0);
	b[1][0] = @(1./2.); b[1][1] = @(0.0);
	b[2][0] = @(0.); b[2][1] = @(3./4.); b[2][2] = @(0.0);
	b[3][0] = @(2./9.); b[3][1] = @(1./3.); b[3][2] = @(4./9.); b[3][3] = @(0.0);
	
	NSMutableArray *c = [NSMutableArray arrayWithCapacity:numStages];
	c[0] = @(2./9.);
	c[1] = @(1./3.);
	c[2] = @(4./9.);
	c[3] = @(0.0);
	
	NSMutableArray *d = [NSMutableArray arrayWithCapacity:numStages];
	d[0] = @(7./24.);
	d[1] = @(1./4.);
	d[2] = @(1./3.);
	d[3] = @(1./8.);
	
	GLAdaptiveRungeKuttaOperation *rk = [[GLAdaptiveRungeKuttaOperation alloc] initMethodWithCoefficientsA: a b: b c: c d: d y: y stepSize: deltaT fsal: YES retainPreviousY: YES order: order fFromTY: aBlock];
	for (GLFunction *rt in rk.relativeToleranceVariables) {
		*(rt.pointerValue) = 1e-3;
	}
	for (GLFunction *at in rk.absoluteToleranceVariables) {
		*(at.pointerValue) = 1e-6;
	}
	
	// Note: this is an absolute mess at the moment
	// We may want a public interface, or internal interface which includes some combination of the following,
	// The internal yn and yn-1
	// -lastInternalTime
	// -lastInternalY
	// -internalTime
	// -internalY
	// The last user requested:
	// -time
	// -y
	
	NSMutableArray *p = [NSMutableArray array];
	p[0] = [NSMutableArray array];
	p[1] = [NSMutableArray array];
	for (NSUInteger i=0; i<rk.nInputs; i++) {
		p[0][i] = [[[rk.previousY[i] minus: rk.currentY[i]] times: @3.0] plus: [[rk.firstStageVariables[i] plus: [rk.lastStageVariables[i] times: @2.0]] times: rk.lastStepSizeVariable]];
		p[1][i] = [[[rk.previousY[i] minus: rk.currentY[i]] times: @2.0] plus: [[rk.firstStageVariables[i] plus: rk.lastStageVariables[i]] times: rk.lastStepSizeVariable]];
	}
	
	NSArray *yout = [rk interpolateAtTime: rk.requestedTimeVariable currentTime: rk.currentTimeVariable lastStepSize: rk.lastStepSizeVariable yNow: rk.currentY fNow: rk.lastStageVariables p: p];
	
	NSMutableArray *topVariables = [[NSMutableArray alloc] init];
	[topVariables addObjectsFromArray: rk.currentY];
	[topVariables addObjectsFromArray: rk.previousY];
	[topVariables addObjectsFromArray: rk.lastStageVariables];
	[topVariables addObjectsFromArray: rk.firstStageVariables];
	[topVariables addObject: rk.requestedTimeVariable];
	[topVariables addObject: rk.currentTimeVariable];
	[topVariables addObject: rk.lastStepSizeVariable];

	NSMutableArray *bottomVariables = [[NSMutableArray alloc] init];
	[bottomVariables addObjectsFromArray: yout];
	
	GLOperationOptimizer *optimizer = [[GLOperationOptimizer alloc] initWithTopVariables: topVariables bottomVariables: bottomVariables];
	variableOperation vectorBlock = optimizer.operationBlock;
	
//	GLOperationVisualizer *visualizer = [[GLOperationVisualizer alloc] initWithTopVariables: topVariables bottomVariables: bottomVariables];
//	NSLog(@"%@", visualizer.graphvisDescription);
	
	NSMutableArray *operandBuffer = [[NSMutableArray alloc] init];
	[operandBuffer addObjectsFromArray: rk.previousYData];
	[operandBuffer addObjectsFromArray: rk.lastStageData];
	[operandBuffer addObjectsFromArray: rk.firstStageData];
	[operandBuffer addObject: rk.requestedTimeData];
	[operandBuffer addObject: rk.currentTimeData];
	[operandBuffer addObject: rk.lastStepSizeData];
	
	NSMutableArray *fullOperandBuffer = [[NSMutableArray alloc] init];
	variableOperation interpBlock = ^(NSArray *result, NSArray *operand, NSArray *bufferArray) {
		[fullOperandBuffer addObjectsFromArray: operand];
		[fullOperandBuffer addObjectsFromArray: operandBuffer];
		
		vectorBlock( result, fullOperandBuffer, bufferArray );
		[fullOperandBuffer removeAllObjects];
	};
	
	rk.interpolationOperation = interpBlock;
    rk.interpolationBuffers = [NSMutableArray array];
    for (GLBuffer *aBuffer in optimizer.internalDataBuffers) {
        [rk.interpolationBuffers addObject: [[GLMemoryPool sharedMemoryPool] dataWithLength: aBuffer.numBytes]];
    }
	
	[yout makeObjectsPerformSelector:@selector(solve)];
	
	return rk;
}

- (id) initWithResult:(NSArray *)result operand:(NSArray *)operand
{
	if ((self=[super initWithResult: result operand: operand])) {
		
		self.relativeToleranceData = [[NSMutableArray alloc] initWithCapacity: operand.count];
		self.absoluteToleranceData = [[NSMutableArray alloc] initWithCapacity: operand.count];
		self.errorVariables = [[NSMutableArray alloc] initWithCapacity: operand.count];
		self.relativeToleranceVariables = [[NSMutableArray alloc] initWithCapacity: operand.count];
		self.absoluteToleranceVariables = [[NSMutableArray alloc] initWithCapacity: operand.count];
		self.errorData = [[NSMutableArray alloc] initWithCapacity: operand.count];
		for (GLFunction *variable in operand ) {
			[self.relativeToleranceVariables addObject: [[GLScalar alloc] initWithType: kGLRealDataFormat forEquation: [operand[0] equation]]];
			[self.absoluteToleranceVariables addObject: [[GLScalar alloc] initWithType: kGLRealDataFormat forEquation: [operand[0] equation]]];
			[self.errorVariables addObject: [GLVariable variableWithPrototype: variable]];
			[self.relativeToleranceData addObject: [[self.relativeToleranceVariables lastObject] data]];
			[self.absoluteToleranceData addObject: [[self.absoluteToleranceVariables lastObject] data]];
			[self.errorData addObject: [[self.errorVariables lastObject] data]];
		}
	}
	
	return self;
}

// coefficients a -- indicate where in the time step
// coefficients b -- scale f at each stage
// coefficients c -- combine all stages to get the highest order solution
// coefficients d -- combine all stages to get the highest order-1 solution
- (instancetype) initMethodWithCoefficientsA: (NSArray *) a b: (NSArray *) b c: (NSArray *) c d: (NSArray *) d y: (NSArray *) y stepSize: (GLFloat) deltaT fsal: (BOOL) isFSAL retainPreviousY: (BOOL) shouldRetainPreviousY order: (NSUInteger) order fFromTY: (FfromTYVector) fFromY
{
	NSUInteger numStages = a.count;
	
	GLScalar *timeStep = [[GLScalar alloc] initWithType: kGLRealDataFormat forEquation: [y[0] equation]]; timeStep.name = @"time_step";
	GLScalar *time = [[GLScalar alloc] initWithType: kGLRealDataFormat forEquation: [y[0] equation]]; time.name = @"time";
	
	NSMutableArray *t = [NSMutableArray arrayWithCapacity: numStages];
	for (NSUInteger i=0; i<a.count; i++) {
		t[i] = [time plus: [timeStep times: a[i]]];
	}
	
	// Store the value at each stage point
	NSMutableArray *yp = [NSMutableArray arrayWithCapacity: numStages];
	yp[0] = y;
	
	NSMutableArray *yold;
	if (shouldRetainPreviousY) {
        NSUInteger i=0;
		yold = [NSMutableArray array];
		for (GLFunction *var in y ) {
            GLFunction *vardup = [var duplicate];
			[yold addObject: vardup];
            var.name = [NSString stringWithFormat:@"yin_%lu",i];
            vardup.name = [NSString stringWithFormat:@"yin_%lu_old",i];
            i++;
		}
	}
	
	// Store the slope at each stage point
	NSMutableArray *f = [NSMutableArray arrayWithCapacity: numStages];
	f[0] = fFromY(t[0], yp[0]);
	
	// Walk through all the stages (i)
	for (NSUInteger i=1; i<numStages; i++) {
		// Compute y at this stage location
		yp[i] = [NSMutableArray arrayWithArray: y];
		for (NSUInteger n=0; n<y.count; n++) {
			GLFunction *ftotal;
			for (NSUInteger j=0; j<i; j++) {
				if (!isZero(b[i][j])) {
					if ( isOne(b[i][j])) ftotal = !ftotal ? f[j][n] : [ftotal plus: f[j][n]];
					else ftotal = !ftotal ? [f[j][n] scalarMultiply: [b[i][j] doubleValue]] : [ftotal plus:[f[j][n] scalarMultiply: [b[i][j] doubleValue]]];
				}
			}
			yp[i][n] = [yp[i][n] plus: [ftotal multiply: timeStep]];
		}

		// Compute the slope at this stage location
		f[i] = fFromY(t[i], yp[i]);
	}
	
	// Now compute the value of y at the very end
	NSMutableArray *yout = [NSMutableArray arrayWithArray: y];
	for (NSUInteger n=0; n<y.count; n++) {
		GLFunction *ftotal;
		for (NSUInteger i=0; i<numStages; i++) {
			if (!isZero(c[i])) {
				if ( isOne( c[i]) ) ftotal = !ftotal ? f[i][n] : [ftotal plus: f[i][n]];
				else ftotal = !ftotal ? [f[i][n] scalarMultiply: [c[i] doubleValue]] : [ftotal plus:[f[i][n] scalarMultiply: [c[i] doubleValue]]];
			}
		}
		yout[n] = [yout[n] plus: [ftotal multiply: timeStep]];
	}
	
	// compute the error
	NSMutableArray *yerr = [NSMutableArray arrayWithCapacity: y.count];
	for (NSUInteger n=0; n<y.count; n++) {
		GLFunction *ftotal;
		for (NSUInteger i=0; i<numStages; i++) {
			NSNumber *dc = @([c[i] doubleValue] - [d[i] doubleValue]);
			if (!isZero(dc)) {
				if ( isOne( dc) ) ftotal = !ftotal ? f[i][n] : [ftotal plus: f[i][n]];
				else ftotal = !ftotal ? [f[i][n] scalarMultiply: [dc doubleValue]] : [ftotal plus:[f[i][n] scalarMultiply: [dc doubleValue]]];
			}
		}
		yerr[n] = [ftotal multiply: timeStep];
	}
	
	// Define the tolerance variables
	NSMutableArray *relativeToleranceArray = [[NSMutableArray alloc] initWithCapacity: y.count];
	NSMutableArray *absoluteToleranceArray = [[NSMutableArray alloc] initWithCapacity: y.count];
	for (NSUInteger i=0; i<y.count; i++) {
		[relativeToleranceArray addObject: [[GLScalar alloc] initWithType: kGLRealDataFormat forEquation: [y[0] equation]]];
		[absoluteToleranceArray addObject: [[GLScalar alloc] initWithType: kGLRealDataFormat forEquation: [y[0] equation]]];
        [(GLFunction *) relativeToleranceArray[i] setName: [NSString stringWithFormat: @"rel_tolerance_yin[%ld]",(unsigned long)i]];
        [(GLFunction *) absoluteToleranceArray[i] setName: [NSString stringWithFormat: @"abs_tolerance_yin[%ld]",(unsigned long)i]];
	}
	
	// *element-wise*
	// |y_err|/max( relTolerance*|y|, absoluteTolerance ) <= 1
	NSMutableArray *error = [[NSMutableArray alloc] initWithCapacity: y.count];
	NSMutableArray *errorVectorData = [[NSMutableArray alloc] initWithCapacity: y.count];
	for (NSUInteger i=0; i < y.count; i++) {
        GLElementErrorOperation *theError = [[GLElementErrorOperation alloc] initWithFirstOperand:yerr[i] secondOperand:yout[i] relativeError:relativeToleranceArray[i] absoluteError:absoluteToleranceArray[i]];
		[error addObject: theError.result[0]];
		[errorVectorData addObject: [theError.result[0] data]];
	}
	
	// (*) The optimized operation graph takes a whole bunch of inputs
	NSMutableArray *topVariables = [[NSMutableArray alloc] init];
	[topVariables addObjectsFromArray: y];
	[topVariables addObject: time];
	[topVariables addObject: timeStep];
	[topVariables addObjectsFromArray: relativeToleranceArray];
	[topVariables addObjectsFromArray: absoluteToleranceArray];
	if (isFSAL) {
		[topVariables addObjectsFromArray: f[0]];
	}
    for (GLVariable *var in f[0]) {
        var.name = @"f[0]";
    }
	
	// (%) but only has a few outputs
	NSMutableArray *bottomVariables = [[NSMutableArray alloc] init];
	[bottomVariables addObjectsFromArray: yout];
	[bottomVariables addObjectsFromArray: error];
	if (isFSAL) {
		[bottomVariables addObjectsFromArray: f[numStages-1]];
	}
	if (shouldRetainPreviousY) {
		[bottomVariables addObjectsFromArray: yold];
	}
    for (GLVariable *var in yout) {
        var.name = @"y_out";
    }
    for (GLVariable *var in error) {
        var.name = @"error_out";
    }
    for (GLVariable *var in f[numStages-1]) {
        var.name = @"f[n]";
    }
    
    // We may have identified the same variable as being at both the top and bottom of the tree.
    // This is a bit of a hack, and I think we can do better than this.
    for (NSUInteger i=0; i<bottomVariables.count; i++) {
        if ( [topVariables containsObject: bottomVariables[i]]) {
            bottomVariables[i] = [bottomVariables[i] duplicate];
        }
    }
	
	GLOperationOptimizer *optimizer = [[GLOperationOptimizer alloc] initWithTopVariables: topVariables bottomVariables: bottomVariables];
	GLOperationVisualizer *vizualizer = [[GLOperationVisualizer alloc] initWithTopVariables: topVariables bottomVariables: bottomVariables];
    
    NSMutableArray *resultPrototypes = [NSMutableArray array];
    for (GLVariable *var in yout) {
        [resultPrototypes addObject: [GLVariable variableWithPrototype: var]];
    }
    
	if ((self=[self initWithResult: resultPrototypes operand: y])) {
		self.graphvisDescription = vizualizer.graphvisDescription;
		variableOperation vectorBlock = optimizer.operationBlock;
		// Go ahead and allocate the memory for those buffers.
        self.dataBuffers = [NSMutableArray array];
		for (GLBuffer *aBuffer in optimizer.internalDataBuffers) {
			[self.dataBuffers addObject: [[GLMemoryPool sharedMemoryPool] dataWithLength: aBuffer.numBytes]];
		}
		self.currentY = y;
		self.stepSize = deltaT;
		self.isFSAL = isFSAL;
		self.shouldRetainPreviousY = shouldRetainPreviousY;
		
		if (self.isFSAL) {
			self.firstStageVariables = [NSMutableArray array];
			self.firstStageData = [NSMutableArray array];
			self.lastStageVariables = [NSMutableArray array];
			self.lastStageData = [NSMutableArray array];
			for (GLFunction *variable in f[0]) {
				[self.firstStageVariables addObject: [GLVariable variableWithPrototype: variable]];
				[self.firstStageData addObject: [[self.firstStageVariables lastObject] data]];
			}
			for (GLFunction *variable in f[numStages-1]) {
				[self.lastStageVariables addObject: [GLVariable variableWithPrototype: variable]];
				[self.lastStageData addObject: [[self.lastStageVariables lastObject] data]];
			}
			
			// If the method isFSAL, we need a way to quickly compute the derivative (fFromTY) at a point.
			NSMutableArray *stagePrepTopVariables = [[NSMutableArray alloc] init];
			[stagePrepTopVariables addObjectsFromArray: y];
			[stagePrepTopVariables addObject: time];
			[stagePrepTopVariables addObject: timeStep];
			
			NSMutableArray *stagePrepBottomVariables = [[NSMutableArray alloc] init];
			[stagePrepBottomVariables addObjectsFromArray: f[0]];
			
			GLOperationOptimizer *stagePrepOptimizer = [[GLOperationOptimizer alloc] initWithTopVariables: stagePrepTopVariables bottomVariables: stagePrepBottomVariables];
			variableOperation prepOperation = stagePrepOptimizer.operationBlock;
			
			NSMutableArray *operandBuffer = [[NSMutableArray alloc] init];
			NSMutableData *timeData = self.currentTimeData;
			NSMutableData *timeStepData = self.stepSizeData;
			NSArray *lData = [NSArray arrayWithArray: self.lastStageData];
			self.prepStageOperation = ^(NSArray *operand, NSArray *theBuffers) {
				[operandBuffer addObjectsFromArray: operand];
				[operandBuffer addObject: timeData];
				[operandBuffer addObject: timeStepData];
				
				prepOperation( lData, operandBuffer, theBuffers );
				[operandBuffer removeAllObjects];
			};
			
            self.prepStageDataBuffers = [NSMutableArray array];
            for (GLBuffer *aBuffer in stagePrepOptimizer.internalDataBuffers) {
                [self.prepStageDataBuffers addObject: [[GLMemoryPool sharedMemoryPool] dataWithLength: aBuffer.numBytes]];
            }
		}
		
		
		// We will need to keep track of how 
		NSMutableArray *nDataElementsArray = [[NSMutableArray alloc] initWithCapacity: error.count];
		for (NSUInteger i=0; i<error.count; i++) {
			[nDataElementsArray addObject: [NSNumber numberWithUnsignedInteger: [[error objectAtIndex: i] nDataElements]]];
		}
		
		NSMutableArray *nBytesArray = [[NSMutableArray alloc] initWithCapacity: self.firstStageVariables.count];
		for (NSUInteger i=0; i<self.firstStageVariables.count; i++) {
			[nBytesArray addObject: [NSNumber numberWithUnsignedInteger: [[self.firstStageVariables objectAtIndex: i] nDataElements]*sizeof(GLFloat)]];
		}
		
		// We are very careful not to copy self or other unnecessary objects with the block we create,
		// hence the extra boilerplate here.
		NSUInteger num = y.count;
		GLFloat safety = 0.5;
		GLFloat grow = 3.5;
		GLFloat errcon = pow(grow/safety, -(order+1.));
		NSMutableArray *operandBuffer = [[NSMutableArray alloc] init];
		NSMutableArray *resultBuffer = [[NSMutableArray alloc] init];
		NSMutableData *timeData = self.currentTimeData;
		NSMutableData *timeStepData = self.stepSizeData;
		NSMutableData *lastTimeStepData = self.lastStepSizeData;
		NSArray *rtData = [NSArray arrayWithArray: self.relativeToleranceData];
		NSArray *atData = [NSArray arrayWithArray: self.absoluteToleranceData];
		NSArray *errData = [NSArray arrayWithArray: self.errorData];
		NSMutableArray *fsData = [NSMutableArray arrayWithArray: self.firstStageData];
		NSMutableArray *lsData = [NSMutableArray arrayWithArray: self.lastStageData];
		NSArray *previousYData = [NSArray arrayWithArray: self.previousYData];
		
//        __block NSUInteger overshoots = 0;
		self.operation = ^(NSArray *result, NSArray *operand, NSArray *bufferArray) {
			// The operand buffer must match the order of (*) above
			[operandBuffer addObjectsFromArray: operand];
			[operandBuffer addObject: timeData];
			[operandBuffer addObject: timeStepData];
			[operandBuffer addObjectsFromArray: rtData];
			[operandBuffer addObjectsFromArray: atData];
			if (isFSAL) {
				[operandBuffer addObjectsFromArray: fsData];
				// Kind of stupid. But, we can't simply point the input to the output because we can't guarantee the order of operations.
				for (NSUInteger i=0; i<fsData.count; i++) {
					memcpy( [fsData[i] mutableBytes], [lsData[i] bytes], [nBytesArray[i] unsignedIntegerValue]);
				}
			}
			
			
			// The result buffer must match the order of (%) above
			[resultBuffer addObjectsFromArray: result];
			[resultBuffer addObjectsFromArray: errData];
			if (isFSAL) {
				[resultBuffer addObjectsFromArray: lsData];
			}
			if (shouldRetainPreviousY) {
				[resultBuffer addObjectsFromArray: previousYData];
			}
			
			GLFloat error=0.0;
			GLFloat *step = (GLFloat *) timeStepData.mutableBytes;
			GLFloat *lastStep = (GLFloat *) lastTimeStepData.mutableBytes;
			
			// This loop is pretty much straight out of Price et. al, Chapter 16.2
			while (1) {
				// Time step
				vectorBlock( resultBuffer, operandBuffer, bufferArray );
				
				// Find the max error
				for (NSUInteger i=0; i<num; i++) {
					GLFloat *localError = (GLFloat *) [errData[i] bytes];
					if (i==0 || *localError > error) {
						error=*localError;
					}
				}
				
				if (error <= 1.0) {
//                    NSLog(@"Succeeded with error: %f. Step size was %f", error, *step);
					break;
				} else {
//                    GLFloat aaa = *step;
					GLFloat htemp = safety* (*step) * pow(error, -1./order);
					(*step) = 0.1 * (*step) >= htemp ? 0.1 * (*step) : htemp;
//                    NSLog(@"Overshot with error: %f. Step size was %f, trying %f", error, aaa, *step);
//                    overshoots++;
//                    NSLog(@"Overshots: %lu", overshoots);
				}
			}
			
			(*lastStep) = (*step);
			
			if ( error > errcon) {
				(*step) = safety*(*step)*pow(error, -1./(order+1));
			} else {
				(*step) = grow * (*step);
			}
			
			[operandBuffer removeAllObjects];
			[resultBuffer removeAllObjects];
		};

	}
	
	return self;
}

- (void) updateCurrentTime
{
	self.currentTime = self.currentTime + self.lastStepSize;
}

- (NSArray *) stepForwardToTime: (GLFloat)time
{
	[super stepForwardToTime: time];
	
	if (self.interpolationOperation) {
		NSMutableArray *yout = [[NSMutableArray alloc] initWithCapacity: self.result.count];
		NSMutableArray *resultBuffer = [[NSMutableArray alloc] initWithCapacity: self.result.count];
		NSMutableArray *operandBuffer = [[NSMutableArray alloc] initWithCapacity: self.operand.count];
		for (GLVariable *variable in self.result) {
			GLVariable *newVar = [GLVariable variableWithPrototype: variable];
			[yout addObject: newVar];
			[resultBuffer addObject: newVar.data];
		}
		for (GLFunction *variable in self.currentY) {
			[operandBuffer addObject: variable.data];
		}
		
//		GLVariable *ydiff = [[self.currentY[0] minus: self.previousY[0]] spatialDomain];
//		GLVariable *fdiff = [[self.lastStageVariables[0] minus: self.firstStageVariables[0]] spatialDomain];
//		[ydiff solve];
//		[fdiff solve];
		
		self.interpolationOperation( resultBuffer, operandBuffer, self.interpolationBuffers);
		
//		ydiff = [[self.previousY[0] minus: yout[0]] spatialDomain];
//		[ydiff solve];
		
		return yout;
	} else {
		return self.currentY;
	}

}

// This roughly matches "evali" in RKSuite.
- (NSArray *) interpolateAtTime: (GLScalar *) t currentTime: (GLScalar *) tNow lastStepSize: (GLScalar *) h yNow: (NSArray *) yNow  fNow: (NSArray *) fNow p: (NSArray *) polyCoeffs
{
	GLScalar *sigma = [[t minus: tNow] dividedBy: h];
	
	NSUInteger nCoeffs = polyCoeffs.count;
	NSUInteger nInputs = yNow.count;
	
	if (nCoeffs < 1) {
		NSLog(@"Your polynomial must have at least one coefficient");
		return nil;
	}
	
	// We start by multiply the highest coefficient by sigma (say, p[n]*sigma)
	NSMutableArray *yout = [NSMutableArray arrayWithCapacity: nInputs];
	for (NSUInteger i=0; i<nInputs; i++) {
		yout[i] = [polyCoeffs[nCoeffs-1][i] times: sigma];
	}
	
	// Now we multiply lower and lower coefficients, polys = (((p[n]*sigma + p[n-1])*sigma +p[n-2])*sigma + p[n-3]*sigma)...
	for (NSInteger n=nCoeffs-2; n>=0; n--) {
		for (NSUInteger i=0; i<nInputs; i++) {
			yout[i] = [[yout[i] plus: polyCoeffs[n][i]] times: sigma];
		}
	}
	
	// Finally, we add to those  (polys + h*f)*sigma + y
	for (NSUInteger i=0; i<nInputs; i++) {
		yout[i] = [[[yout[i] plus: [fNow[i] times: h]] times: sigma] plus: yNow[i]];
	}
	
	return yout;
}

@dynamic requestedTime;

@dynamic stepSizeVariable;
@dynamic lastStepSizeVariable;
@dynamic currentTimeVariable;
@dynamic requestedTimeVariable;

@dynamic stepSizeData;
@dynamic lastStepSizeData;
@dynamic currentTimeData;
@dynamic requestedTimeData;

@dynamic shouldRetainPreviousY;
@dynamic previousY;
@dynamic previousYData;

@dynamic isFSAL;
@dynamic firstStageVariables;
@dynamic lastStageVariables;
@dynamic firstStageData;
@dynamic lastStageData;
@dynamic stageIsPrepped;
@dynamic prepStageOperation;



- (void) setRelativeTolerance: (NSArray *) relativeTolerance
{
	for (NSUInteger i=0; i<relativeTolerance.count; i++) {
		NSMutableData *data = [self.relativeToleranceData objectAtIndex: i];
		GLFloat *a = (GLFloat *) data.mutableBytes;
		*a = [[relativeTolerance objectAtIndex: i] doubleValue];
	}
}

- (NSArray *) relativeTolerance
{
	NSMutableArray *array = [[NSMutableArray alloc] initWithCapacity: self.relativeToleranceData.count];
	for (NSUInteger i=0; i< self.relativeToleranceData.count; i++) {
		NSMutableData *data = [self.relativeToleranceData objectAtIndex: i];
		GLFloat *a = (GLFloat *) data.mutableBytes;
		[array addObject: [NSNumber numberWithDouble: *a]];
	}
	
	return array;
}

- (void) setAbsoluteTolerance: (NSArray *) absoluteTolerance
{
	for (NSUInteger i=0; i<absoluteTolerance.count; i++) {
		NSMutableData *data = [self.absoluteToleranceData objectAtIndex: i];
		GLFloat *a = (GLFloat *) data.mutableBytes;
		*a = [[absoluteTolerance objectAtIndex: i] doubleValue];
	}
}

- (NSArray *) absoluteTolerance
{
	NSMutableArray *array = [[NSMutableArray alloc] initWithCapacity: self.absoluteToleranceData.count];
	for (NSUInteger i=0; i< self.absoluteToleranceData.count; i++) {
		NSMutableData *data = [self.absoluteToleranceData objectAtIndex: i];
		GLFloat *a = (GLFloat *) data.mutableBytes;
		[array addObject: [NSNumber numberWithDouble: *a]];
	}
	
	return array;
}

@end
