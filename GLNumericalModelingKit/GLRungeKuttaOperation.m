//
//  GLRungeKuttaOperation.m
//  GLNumericalModelingKitCopy
//
//  Created by Jeffrey J. Early on 10/26/12.
//
//

#import "GLRungeKuttaOperation.h"

#import "GLOperationOptimizer.h"
#import "GLMemoryPool.h"
#import "GLUnaryOperations.h"
#import "GLVectorVectorOperations.h"
#import "GLOperationVisualizer.h"

// This type is designed to take an array "y", and then write "f" to the appropriate data chunks that are fixed.
typedef void (^stagePrepOperation)(NSArray *);


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

@property(strong) NSArray *dataBuffers;
@property GLFloat stepSize;
@property NSUInteger nInputs;

@property(nonatomic) GLFloat requestedTime;

@property(strong) GLVariable *stepSizeVariable;
@property(strong) GLVariable *lastStepSizeVariable;
@property(strong) GLVariable *currentTimeVariable;
@property(strong) GLVariable *requestedTimeVariable;

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
		self.currentTimeVariable = [GLVariable variableOfRealTypeWithDimensions: @[] forEquation: [operand[0] equation]];
		self.stepSizeVariable = [GLVariable variableOfRealTypeWithDimensions: @[] forEquation: [operand[0] equation]];
		self.lastStepSizeVariable = [GLVariable variableOfRealTypeWithDimensions: @[] forEquation: [operand[0] equation]];
		self.requestedTimeVariable = [GLVariable variableOfRealTypeWithDimensions: @[] forEquation: [operand[0] equation]];
		
		self.currentTimeData = self.currentTimeVariable.data;
		self.stepSizeData = self.stepSizeVariable.data;
		self.lastStepSizeData = self.lastStepSizeVariable.data;
		self.requestedTimeData = self.requestedTimeVariable.data;
		
		self.nInputs = operand.count;
		
		self.exitOnBlowUp = YES;
		
		self.previousY = [NSMutableArray array];
		self.previousYData = [NSMutableArray array];
		for (GLVariable *variable in operand ) {
			[self.previousY addObject: [GLVariable variableOfType: variable.dataFormat withDimensions: variable.dimensions forEquation: variable.equation]];
			[self.previousYData addObject: [[self.previousY lastObject] data]];
		}
	}
	
	return self;
}

// This assumes *fixed* step size, and no error correction.
- (id) initMethodWithCoefficientsA: (NSArray *) a b: (NSArray *) b c: (NSArray *) c y: (NSArray *) y stepSize: (GLFloat) deltaT fFromTY: (FfromTYVector) fFromY
{
	NSUInteger numStages = a.count;
	
	GLVariable *time = [GLVariable variableOfRealTypeWithDimensions: @[] forEquation: [y[0] equation]];
	NSMutableArray *t = [NSMutableArray arrayWithCapacity: numStages];
	for (NSUInteger i=0; i<a.count; i++) {
		t[i] = [time scalarAdd: [a[0] doubleValue]*deltaT];
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
	
	// Re-create similar variables, but without the dependencies.
	NSMutableArray *result = [[NSMutableArray alloc] initWithCapacity: yout.count];
	for (GLVariable *variable in yout) {
		[result addObject: [GLVariable variableOfType: variable.dataFormat withDimensions: variable.dimensions forEquation: variable.equation]];
	}
	
	if ((self=[self initWithResult: result operand: y])) {
		
		unaryVectorOperation vectorBlock = optimizer.unaryVectorOperationBlock;
		NSMutableArray *operandArray = [NSMutableArray arrayWithCapacity: yfull.count];
		NSMutableData *tData = self.currentTimeData;
		self.blockOperation = ^(NSArray *result, NSArray *operand) {
			[operandArray addObjectsFromArray: operand];
			[operandArray addObject: tData];
			
			vectorBlock( result, operandArray);
			
			[operandArray removeAllObjects];
		};
		
		self.dataBuffers = optimizer.internalDataBuffers;
		self.currentY = y;
		self.stepSize = deltaT;
		self.lastStepSize = deltaT;
		self.exitOnBlowUp = YES;
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
		GLVariable *newVar = [GLVariable variableOfType: variable.dataFormat withDimensions: variable.dimensions forEquation: variable.equation];
		[yout addObject: newVar];
		[resultBuffer addObject: newVar.data];
	}
	for (GLVariable *variable in self.currentY) {
		[operandBuffer addObject: variable.data];
	}
	
	// If the algorithm isFSAL, then we plan on using the last stage from the previous time step.
	// However, if there was no previous time step, then the stage won't actually be there and we need to compute it!
	if (self.isFSAL && !self.stageIsPrepped) {
		self.prepStageOperation( operandBuffer );
		self.stageIsPrepped = YES;
	}
	
	// Take one step forward in time
	self.blockOperation( resultBuffer, operandBuffer );
	self.totalIterations = self.totalIterations + 1;
	[self updateCurrentTime];
	
    if (self.currentTime < time)
    {
        // We can't feed the same databuffer to the input and output, because we'll overwrite
        // the input. This could be bad if happened to take too large of a step.
        NSMutableArray *youtOut = [[NSMutableArray alloc] initWithCapacity: self.result.count];
        NSMutableArray *resultBufferOut= [[NSMutableArray alloc] initWithCapacity: self.result.count];
        for (GLVariable *variable in self.result) {
            GLVariable *newVar = [GLVariable variableOfType: variable.dataFormat withDimensions: variable.dimensions forEquation: variable.equation];
            [youtOut addObject: newVar];
            [resultBufferOut addObject: newVar.data];
        }
        
        while ( self.currentTime < time )
        {
            self.blockOperation( resultBufferOut, resultBuffer );
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

@interface GLAdaptiveRungeKuttaOperation ()

- (id) initMethodWithCoefficientsA: (NSArray *) a b: (NSArray *) b c: (NSArray *) c d: (NSArray *) d y: (NSArray *) y stepSize: (GLFloat) deltaT fsal: (BOOL) isFSAL retainPreviousY: (BOOL) shouldRetainPreviousY  order: (NSUInteger) order fFromTY: (FfromTYVector) fFromY;
- (NSArray *) interpolateAtTime: (GLVariable *) t currentTime: (GLVariable *) tNow lastStepSize: (GLVariable *) h yNow: (NSArray *) yNow  fNow: (NSArray *) fNow p: (NSArray *) polyCoeffs;



// Error variables---one for each input
@property(strong) NSMutableArray *relativeToleranceVariables;
@property(strong) NSMutableArray *absoluteToleranceVariables;
@property(strong) NSMutableArray *errorVariables;

@property(strong) NSMutableArray *relativeToleranceData;
@property(strong) NSMutableArray *absoluteToleranceData;
@property(strong) NSMutableArray *errorData;

// Hold onto all the stages, in case we want to interpolate
@property(strong) NSArray *stageVariables;

@property(copy) unaryVectorOperation interpolationOperation;

// Re-declared

@property(nonatomic) GLFloat requestedTime;

@property(strong) GLVariable *stepSizeVariable;
@property(strong) GLVariable *lastStepSizeVariable;
@property(strong) GLVariable *currentTimeVariable;
@property(strong) GLVariable *requestedTimeVariable;

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
	for (GLVariable *rt in rk.relativeToleranceVariables) {
		*(rt.pointerValue) = 1e-5;
	}
	for (GLVariable *at in rk.absoluteToleranceVariables) {
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
	for (GLVariable *rt in rk.relativeToleranceVariables) {
		*(rt.pointerValue) = 1e-3;
	}
	for (GLVariable *at in rk.absoluteToleranceVariables) {
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
		p[0][i] = [[[rk.previousY[i] minus: rk.currentY[i]] scalarMultiply: 3.0] plus: [[rk.firstStageVariables[i] plus: [rk.lastStageVariables[i] scalarMultiply: 2]] multiply: rk.lastStepSizeVariable]];
		p[1][i] = [[[rk.previousY[i] minus: rk.currentY[i]] scalarMultiply: 2.0] plus: [[rk.firstStageVariables[i] plus: rk.lastStageVariables[i]] multiply: rk.lastStepSizeVariable]];
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
	unaryVectorOperation vectorBlock = optimizer.unaryVectorOperationBlock;
	
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
	unaryVectorOperation interpBlock = ^(NSArray *result, NSArray *operand) {
		[fullOperandBuffer addObjectsFromArray: operand];
		[fullOperandBuffer addObjectsFromArray: operandBuffer];
		
		vectorBlock( result, fullOperandBuffer );
		[fullOperandBuffer removeAllObjects];
	};
	
	rk.interpolationOperation = interpBlock;
	NSMutableArray *internalBuffers = [NSMutableArray arrayWithArray: optimizer.internalDataBuffers];
	[internalBuffers addObjectsFromArray: rk.dataBuffers];
	rk.dataBuffers = internalBuffers;
	
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
		for (GLVariable *variable in operand ) {
			[self.relativeToleranceVariables addObject: [GLVariable variableOfRealTypeWithDimensions: @[] forEquation: [operand[0] equation]]];
			[self.absoluteToleranceVariables addObject: [GLVariable variableOfRealTypeWithDimensions: @[] forEquation: [operand[0] equation]]];
			[self.errorVariables addObject: [GLVariable variableOfType: variable.dataFormat withDimensions: variable.dimensions forEquation: variable.equation]];
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
- (id) initMethodWithCoefficientsA: (NSArray *) a b: (NSArray *) b c: (NSArray *) c d: (NSArray *) d y: (NSArray *) y stepSize: (GLFloat) deltaT fsal: (BOOL) isFSAL retainPreviousY: (BOOL) shouldRetainPreviousY order: (NSUInteger) order fFromTY: (FfromTYVector) fFromY
{
	NSUInteger numStages = a.count;
	
	GLVariable *timeStep = [GLVariable variableOfRealTypeWithDimensions: @[] forEquation: [y[0] equation]]; timeStep.name = @"time_step";
	GLVariable *time = [GLVariable variableOfRealTypeWithDimensions: @[] forEquation: [y[0] equation]]; time.name = @"time";
	
	NSMutableArray *t = [NSMutableArray arrayWithCapacity: numStages];
	for (NSUInteger i=0; i<a.count; i++) {
		t[i] = [time plus: [timeStep scalarMultiply: [a[0] doubleValue]]];
	}
	
	// Store the value at each stage point
	NSMutableArray *yp = [NSMutableArray arrayWithCapacity: numStages];
	yp[0] = y;
	
	NSMutableArray *yold;
	if (shouldRetainPreviousY) {
		yold = [NSMutableArray array];
		for (GLVariable *var in y ) {
			[yold addObject: [var duplicate]];
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
			GLVariable *ftotal;
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
		GLVariable *ftotal;
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
		GLVariable *ftotal;
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
		[relativeToleranceArray addObject: [GLVariable variableOfRealTypeWithDimensions: @[] forEquation: [y[0] equation]]];
		[absoluteToleranceArray addObject: [GLVariable variableOfRealTypeWithDimensions: @[] forEquation: [y[0] equation]]];
        [(GLVariable *) relativeToleranceArray[i] setName: [NSString stringWithFormat: @"rel_tolerance_yin[%ld]",(unsigned long)i]];
        [(GLVariable *) absoluteToleranceArray[i] setName: [NSString stringWithFormat: @"abs_tolerance_yin[%ld]",(unsigned long)i]];
	}
	
	// *element-wise*
	// |y_err|/max( relTolerance*|y|, absoluteTolerance ) <= 1
	NSMutableArray *error = [[NSMutableArray alloc] initWithCapacity: y.count];
	NSMutableArray *errorVectorData = [[NSMutableArray alloc] initWithCapacity: y.count];
	for (NSUInteger i=0; i < y.count; i++) {
		GLAbsoluteValueOperation *absRelErr = [[GLAbsoluteValueOperation alloc] initWithOperand: [yout[i] multiply: relativeToleranceArray[i]] shouldUseComplexArithmetic: NO];
		GLAbsoluteValueOperation *absYErr = [[GLAbsoluteValueOperation alloc] initWithOperand: yerr[i] shouldUseComplexArithmetic: NO];
		GLDivisionOperation * op = [[GLDivisionOperation alloc] initWithFirstOperand: absYErr.result secondOperand: [absRelErr.result absMax: absoluteToleranceArray[i]] shouldUseComplexArithmetic: NO];
		[error addObject: op.result];
		[errorVectorData addObject: op.result.data];
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
    
    // We may have identified the same variable as being at both the top and bottom of the tree.
#warning This is a bit of a hack, and I think we can do better than this.
    for (NSUInteger i=0; i<bottomVariables.count; i++) {
        if ( [topVariables containsObject: bottomVariables[i]]) {
            bottomVariables[i] = [bottomVariables[i] duplicate];
        }
    }
	
	GLOperationOptimizer *optimizer = [[GLOperationOptimizer alloc] initWithTopVariables: topVariables bottomVariables: bottomVariables];
	
	// Re-create result variables, but without the dependencies.
	NSMutableArray *resultVars = [[NSMutableArray alloc] initWithCapacity: yout.count];
	for (GLVariable *variable in yout) {
		[resultVars addObject: [GLVariable variableOfType: variable.dataFormat withDimensions: variable.dimensions forEquation: variable.equation]];
	}
    
//        GLOperationVisualizer *visualizer = [[GLOperationVisualizer alloc] initWithTopVariables: topVariables bottomVariables:bottomVariables];
//        NSLog(@"%@", visualizer.graphvisDescription);
	
	if ((self=[self initWithResult: resultVars operand: y])) {
		
		unaryVectorOperation vectorBlock = optimizer.unaryVectorOperationBlock;
		self.dataBuffers = optimizer.internalDataBuffers;
		self.currentY = y;
		self.stepSize = deltaT;
		self.isFSAL = isFSAL;
		self.shouldRetainPreviousY = shouldRetainPreviousY;
		
		if (self.isFSAL) {
			self.firstStageVariables = [NSMutableArray array];
			self.firstStageData = [NSMutableArray array];
			self.lastStageVariables = [NSMutableArray array];
			self.lastStageData = [NSMutableArray array];
			for (GLVariable *variable in f[0]) {
				[self.firstStageVariables addObject: [GLVariable variableOfType: variable.dataFormat withDimensions: variable.dimensions forEquation: variable.equation]];
				[self.firstStageData addObject: [[self.firstStageVariables lastObject] data]];
			}
			for (GLVariable *variable in f[numStages-1]) {
				[self.lastStageVariables addObject: [GLVariable variableOfType: variable.dataFormat withDimensions: variable.dimensions forEquation: variable.equation]];
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
			unaryVectorOperation prepOperation = stagePrepOptimizer.unaryVectorOperationBlock;
			
			NSMutableArray *operandBuffer = [[NSMutableArray alloc] init];
			NSMutableData *timeData = self.currentTimeData;
			NSMutableData *timeStepData = self.stepSizeData;
			NSArray *lData = [NSArray arrayWithArray: self.lastStageData];
			self.prepStageOperation = ^(NSArray *operand) {
				[operandBuffer addObjectsFromArray: operand];
				[operandBuffer addObject: timeData];
				[operandBuffer addObject: timeStepData];
				
				prepOperation( lData, operandBuffer );
				[operandBuffer removeAllObjects];
			};
			
			NSMutableArray *internalBuffers = [NSMutableArray arrayWithArray: stagePrepOptimizer.internalDataBuffers];
			[internalBuffers addObjectsFromArray: self.dataBuffers];
			self.dataBuffers = internalBuffers;
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
		GLFloat safety = 0.8;
		GLFloat errcon = pow(5./safety, -(order+1.));
		NSMutableArray *operandBuffer = [[NSMutableArray alloc] init];
		NSMutableArray *resultBuffer = [[NSMutableArray alloc] init];
		NSMutableData *timeData = self.currentTimeData;
		NSMutableData *timeStepData = self.stepSizeData;
		NSMutableData *lastTimeStepData = self.lastStepSizeData;
		NSArray *rtData = [NSArray arrayWithArray: self.relativeToleranceData];
		NSArray *atData = [NSArray arrayWithArray: self.absoluteToleranceData];
		NSArray *errData = [NSArray arrayWithArray: self.errorData];
		NSMutableArray *fsData = [NSArray arrayWithArray: self.firstStageData];
		NSMutableArray *lsData = [NSArray arrayWithArray: self.lastStageData];
		NSArray *previousYData = [NSArray arrayWithArray: self.previousYData];
		
        __block NSUInteger overshoots = 0;
		self.blockOperation = ^(NSArray *result, NSArray *operand) {
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
			
			GLFloat error;
			GLFloat *step = (GLFloat *) timeStepData.mutableBytes;
			GLFloat *lastStep = (GLFloat *) lastTimeStepData.mutableBytes;
			
			// This loop is pretty much straight out of Price et. al, Chapter 16.2
			while (1) {
				// Time step
				vectorBlock( resultBuffer, operandBuffer );
				
				// Find the max error
				for (NSUInteger i=0; i<num; i++) {
					NSData *errorData = errData[i];
					NSUInteger nDataElements = [[nDataElementsArray objectAtIndex: i] unsignedIntegerValue];
					GLFloat localError;
					vGL_maxv( (float *) errorData.bytes, 1, &localError, nDataElements);
					if (i==0 || localError > error) {
						error=localError;
					}
				}
				
				if (error <= 1.0) {
					break;
				} else {
//                    GLFloat aaa = *step;
					GLFloat htemp = safety* (*step) * pow(error, -1./order);
					(*step) = 0.1 * (*step) >= htemp ? 0.1 * (*step) : htemp;
//                    NSLog(@"damn, overshot. was %f, trying %f", aaa, *step);
//                    overshoots++;
//                    NSLog(@"Overshots: %lu", overshoots);
				}
			}
			
			(*lastStep) = (*step);
			
			if ( error > errcon) {
				(*step) = safety*(*step)*pow(error, -1./(order+1));
			} else {
				(*step) = 5.0 * (*step);
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
			GLVariable *newVar = [GLVariable variableOfType: variable.dataFormat withDimensions: variable.dimensions forEquation: variable.equation];
			[yout addObject: newVar];
			[resultBuffer addObject: newVar.data];
		}
		for (GLVariable *variable in self.currentY) {
			[operandBuffer addObject: variable.data];
		}
		
//		GLVariable *ydiff = [[self.currentY[0] minus: self.previousY[0]] spatialDomain];
//		GLVariable *fdiff = [[self.lastStageVariables[0] minus: self.firstStageVariables[0]] spatialDomain];
//		[ydiff solve];
//		[fdiff solve];
		
		self.interpolationOperation( resultBuffer, operandBuffer);
		
//		ydiff = [[self.previousY[0] minus: yout[0]] spatialDomain];
//		[ydiff solve];
		
		return yout;
	} else {
		return self.currentY;
	}

}

// This roughly matches "evali" in RKSuite.
- (NSArray *) interpolateAtTime: (GLVariable *) t currentTime: (GLVariable *) tNow lastStepSize: (GLVariable *) h yNow: (NSArray *) yNow  fNow: (NSArray *) fNow p: (NSArray *) polyCoeffs
{
	GLVariable *sigma = [[t minus: tNow] dividedBy: h];
	
	NSUInteger nCoeffs = polyCoeffs.count;
	NSUInteger nInputs = yNow.count;
	
	if (nCoeffs < 1) {
		NSLog(@"Your polynomial must have at least one coefficient");
		return nil;
	}
	
	// We start by multiply the highest coefficient by sigma (say, p[n]*sigma)
	NSMutableArray *yout = [NSMutableArray arrayWithCapacity: nInputs];
	for (NSUInteger i=0; i<nInputs; i++) {
		yout[i] = [polyCoeffs[nCoeffs-1][i] multiply: sigma];
	}
	
	// Now we multiply lower and lower coefficients, polys = (((p[n]*sigma + p[n-1])*sigma +p[n-2])*sigma + p[n-3]*sigma)...
	for (NSInteger n=nCoeffs-2; n>=0; n--) {
		for (NSUInteger i=0; i<nInputs; i++) {
			yout[i] = [[yout[i] plus: polyCoeffs[n][i]] multiply: sigma];
		}
	}
	
	// Finally, we add to those  (polys + h*f)*sigma + y
	for (NSUInteger i=0; i<nInputs; i++) {
		yout[i] = [[[yout[i] plus: [fNow[i] multiply: h]] multiply: sigma] plus: yNow[i]];
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
