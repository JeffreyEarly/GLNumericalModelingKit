//
//  GLIntegrationOperations.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 4/19/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import "GLIntegrationOperations.h"
#import "GLOperationOptimizer.h"
#import "GLMemoryPool.h"
#import "GLUnaryOperations.h"
#import "GLVectorVectorOperations.h"

BOOL isFinite( NSData *data, NSUInteger nDataElements )
{
	GLFloat *pointerValue = (GLFloat *) data.bytes;
	for (NSUInteger i=0; i<nDataElements; i++) {
		if ( !isfinite(pointerValue[i])) return NO;
	}
	return YES;
}

@interface GLIntegrationOperation ()
@property(strong) NSArray *dataBuffers;
@end

/************************************************/
/*		GLIntegrationOperation					*/
/************************************************/

#pragma mark -
#pragma mark GLIntegrationOperation
#pragma mark

@interface GLIntegrationOperation ()
@property GLFloat stepSize;
@property double currentTime;
@property GLFloat lastStepSize;
@property NSMutableData *currentTimeData;
@end

@implementation GLIntegrationOperation
{
	double _initialTime;
}
@synthesize lastStepSize;
@synthesize stepSize;
@synthesize currentTime;
@synthesize initialTime=_initialTime;
@synthesize exitOnBlowUp;
@synthesize currentTimeData;

+ (id) rungeKutta4AdvanceY: (GLVariable *) y stepSize: (double) deltaT fFromY: (FfromY) fFromY
{
	GLVariable *f = fFromY(y);
	// Half a step forward and find the new slope (f) at this point
	GLVariable *f2 = fFromY([[f scalarMultiply: 0.5*deltaT] plus: y]);
	// Find the new slope at this 2nd point
	GLVariable *f3 = fFromY([[f2 scalarMultiply: 0.5*deltaT] plus: y]);	
	// Go a full step forward with the 2nd new slope
	GLVariable *f4 = fFromY([[f3 scalarMultiply: deltaT] plus: y]);
	// yout = y + (Delta/6)*(
	GLVariable *yout = [y plus: [[[f plus: f4] plus: [[f3 plus: f2] scalarMultiply: 2.0]] scalarMultiply: deltaT/6.0]];
	
	GLOperationOptimizer *optimizer = [[GLOperationOptimizer alloc] initWithTopVariables: @[y] bottomVariables: @[yout]];
	
	// We can't use yout as the result variable, because it has a long string of irrelevant dependencies.
	GLVariable *result = [GLVariable variableOfType: yout.dataFormat withDimensions: yout.dimensions forEquation: yout.equation];
	GLIntegrationOperation *integrationOperation = [[GLIntegrationOperation alloc] initWithResult: result operand: y];
	integrationOperation.blockOperation = optimizer.unaryOperationBlock;
	integrationOperation.dataBuffers = optimizer.internalDataBuffers;
	integrationOperation.stepSize = deltaT;
	integrationOperation.currentTime = 0;
	integrationOperation.initialTime = 0;
	integrationOperation.exitOnBlowUp = YES;
	
	if (!integrationOperation.blockOperation) {
		return nil;
	}
	
	return integrationOperation;
}

+ (id) rungeKutta4AdvanceY: (GLVariable *) y stepSize: (double) deltaT fFromTY: (FfromTY) fFromTY
{
	GLVariable *t1 = [GLVariable variableOfRealTypeWithDimensions: @[] forEquation: y.equation];
	GLVariable *t2 = [t1 scalarAdd: 0.5*deltaT];
	GLVariable *t4 = [t1 scalarAdd: deltaT];
	
	GLVariable *f = fFromTY(t1, y);
	// Half a step forward and find the new slope (f) at this point
	GLVariable *f2 = fFromTY(t2, [[f scalarMultiply: 0.5*deltaT] plus: y]);
	// Find the new slope at this 2nd point
	GLVariable *f3 = fFromTY(t2, [[f2 scalarMultiply: 0.5*deltaT] plus: y]);
	// Go a full step forward with the 2nd new slope
	GLVariable *f4 = fFromTY(t4, [[f3 scalarMultiply: deltaT] plus: y]);
	// yout = y + (Delta/6)*(
	GLVariable *yout = [y plus: [[[f plus: f4] plus: [[f3 plus: f2] scalarMultiply: 2.0]] scalarMultiply: deltaT/6.0]];
	
	GLOperationOptimizer *optimizer = [[GLOperationOptimizer alloc] initWithTopVariables: @[t1, y] bottomVariables: @[yout]];
	
	// We can't use yout as the result variable, because it has a long string of irrelevant dependencies.
	GLVariable *result = [GLVariable variableOfType: yout.dataFormat withDimensions: yout.dimensions forEquation: yout.equation];
	GLIntegrationOperation *integrationOperation = [[GLIntegrationOperation alloc] initWithResult: result operand: y];
	integrationOperation.blockOperation = optimizer.unaryOperationBlock;
	integrationOperation.dataBuffers = optimizer.internalDataBuffers;
	integrationOperation.stepSize = deltaT;
	integrationOperation.currentTime = 0;
	integrationOperation.initialTime = 0;
	integrationOperation.exitOnBlowUp = YES;
	integrationOperation.currentTimeData = t1.data;
	
	unaryVectorOperation vectorBlock = optimizer.unaryVectorOperationBlock;
	NSMutableData *tData = t1.data;
	integrationOperation.blockOperation = ^(NSMutableData *result, NSData *operand) {
		vectorBlock( @[result], @[tData, operand]);
		
		GLFloat *a = (GLFloat *) tData.mutableBytes;
		*a += deltaT;
	};
	
	if (!integrationOperation.blockOperation) {
		return nil;
	}
	
	return integrationOperation;
}

+ (id) rungeKutta5AdvanceY: (GLVariable *) y stepSize: (double) deltaT fFromY: (FfromY) fFromY
{
	GLFloat b21=1./5.;
	GLFloat b31=3./40., b32=9./40.;
	GLFloat b41=3./10., b42=-9./10., b43=6./5.;
	GLFloat b51=-11./54., b52=5./2., b53=-70./27., b54=35./27.;
	GLFloat b61=1631./55296., b62=175./512., b63=575./13824., b64=44275./110592., b65=253./4096.;
	GLFloat c1=37./378., c3=250./621., c4=125./594., c6=512./1771.;
	
	GLVariable *f = fFromY(y);
	GLVariable *f2 = fFromY([y plus: [f scalarMultiply: b21*deltaT]]);
	GLVariable *f3 = fFromY([y plus: [[[f scalarMultiply: b31] plus: [f2 scalarMultiply: b32]] scalarMultiply: deltaT]]);	
	GLVariable *f4 = fFromY([y plus: [[[[f scalarMultiply: b41] plus: [f2 scalarMultiply: b42]] plus: [f3 scalarMultiply: b43]] scalarMultiply: deltaT]]);
	GLVariable *f5 = fFromY([y plus: [[[[f scalarMultiply: b51] plus: [f2 scalarMultiply: b52]] plus: [[f3 scalarMultiply: b53] plus: [f4 scalarMultiply: b54]]] scalarMultiply: deltaT]]);
	GLVariable *f6 = fFromY([y plus: [[[[[f scalarMultiply: b61] plus: [f2 scalarMultiply: b62]] plus: [[f3 scalarMultiply: b63] plus: [f4 scalarMultiply: b64]]] plus: [f5 scalarMultiply: b65]] scalarMultiply: deltaT]]);
	GLVariable *yout = [y plus: [[[[f scalarMultiply: c1] plus: [f3 scalarMultiply: c3]] plus: [[f4 scalarMultiply: c4] plus: [f6 scalarMultiply: c6]]] scalarMultiply: deltaT]];
	
	GLOperationOptimizer *optimizer = [[GLOperationOptimizer alloc] initWithTopVariables: @[y] bottomVariables: @[yout]];
	
	// We can't use yout as the result variable, because it has a long string of irrelevant dependencies.
	GLVariable *result = [GLVariable variableOfType: yout.dataFormat withDimensions: yout.dimensions forEquation: yout.equation];
	GLIntegrationOperation *integrationOperation = [[GLIntegrationOperation alloc] initWithResult: result operand: y];
	integrationOperation.blockOperation = optimizer.unaryOperationBlock;
	integrationOperation.dataBuffers = optimizer.internalDataBuffers;
	integrationOperation.stepSize = deltaT;
	integrationOperation.currentTime = 0;
	integrationOperation.initialTime = 0;
	integrationOperation.lastStepSize = 0;
	integrationOperation.exitOnBlowUp = YES;
	
	if (!integrationOperation.blockOperation) {
		return nil;
	}
	
	return integrationOperation;
}

- (void) setInitialTime:(double)initialTime
{
	self.currentTime = initialTime + (self.currentTime-_initialTime);
	_initialTime=initialTime;
	if (self.currentTimeData) {
		double *a = (double *) self.currentTimeData.mutableBytes;
		*a = self.currentTime;
	}
}

- (GLVariable *) stepForward: (GLVariable *) y
{	
	[y solve];
	
	self.lastStepSize = self.stepSize;
	GLVariable *yout = [GLVariable variableOfType: self.result.dataFormat withDimensions: self.result.dimensions forEquation: self.result.equation];
	self.blockOperation( yout.data, y.data );
	self.currentTime = self.currentTime + self.stepSize;
	
	if (self.exitOnBlowUp) {
		if (!isFinite(yout.data, yout.nDataElements)) {
			NSLog(@"Blow up detected! Exiting...");
			return nil;
		}
	}
	
	return yout;
}

- (GLVariable *) stepForward: (GLVariable *) y toTime: (double) time
{
	// If we've already exceeded the time, just return the same variable.
	if (self.currentTime >= time) {
		return y;
	}
	
	// If we haven't, then take one time step forward...
	[y solve];
	
	self.lastStepSize = self.stepSize;
	GLVariable *yout = [GLVariable variableOfType: self.result.dataFormat withDimensions: self.result.dimensions forEquation: self.result.equation];
	self.blockOperation( yout.data, y.data );
	self.currentTime = self.currentTime + self.stepSize;
	
	if (self.exitOnBlowUp) {
		if (!isFinite(yout.data, yout.nDataElements)) {
			NSLog(@"Blow up detected! Exiting...");
			return nil;
		}
	}
	
	// ...and as many as necessary
	while (self.currentTime < time) {
        
        // The last output is the new input.
		self.blockOperation( yout.data, yout.data );
		self.currentTime = self.currentTime + self.stepSize;
		
		if (self.exitOnBlowUp) {
			if (!isFinite(yout.data, yout.nDataElements)) {
				NSLog(@"Blow up detected! Exiting...");
				return nil;
			}
		}
	}
	
	return yout;
}

@synthesize dataBuffers;

- (void) dealloc
{
	if (self.dataBuffers) {
		for (NSMutableData *data in self.dataBuffers) {
			[[GLMemoryPool sharedMemoryPool] returnData: data];
		}
	}
}

@end

/************************************************/
/*		GLAdaptiveIntegrationOperation			*/
/************************************************/

#pragma mark -
#pragma mark GLAdaptiveIntegrationOperation
#pragma mark

@interface GLAdaptiveIntegrationOperation ( )

@property(strong) NSArray *dataBuffers;
@property NSMutableData *stepSizeData;
@property NSMutableData *lastStepSizeData;
@property NSMutableData *relativeToleranceData;
@property NSMutableData *absoluteToleranceData;
@property(copy) unaryVectorOperation vectorBlockOperation;
@property GLVariable *errorVector;
@property double lastStepSize;

@end

@implementation GLAdaptiveIntegrationOperation
{
	NSArray *_dataBuffers;
	NSMutableData *_stepSizeData;
	NSMutableData *_lastStepSizeData;
	NSMutableData *_relativeToleranceData;
	NSMutableData *_absoluteToleranceData;
	
	double _initialTime;
	double _currentTime;
}

@synthesize dataBuffers=_dataBuffers;
@synthesize stepSizeData=_stepSizeData;
@synthesize lastStepSizeData=_lastStepSizeData;
@synthesize relativeToleranceData=_relativeToleranceData;
@synthesize absoluteToleranceData=_absoluteToleranceData;
@synthesize vectorBlockOperation;
@synthesize errorVector;

@synthesize initialTime=_initialTime;
@synthesize currentTime=_currentTime;

+ (id) rungeKutta4AdvanceY:(GLVariable *)y dynamicStepSize: (tFromY) deltaTfromY fFromY:(FfromY)fFromY
{
	GLVariable *deltaT = deltaTfromY(y);
	
	GLVariable *f = fFromY(y);
	// Half a step forward and find the new slope (f) at this point
	GLVariable *f2 = fFromY([[f multiply: [deltaT scalarMultiply: 0.5]] plus: y]);
	// Find the new slope at this 2nd point
	GLVariable *f3 = fFromY([[f2 multiply: [deltaT scalarMultiply: 0.5]] plus: y]);
	// Go a full step forward with the 2nd new slope
	GLVariable *f4 = fFromY([[f3 multiply: deltaT] plus: y]);
	// yout = y + (Delta/6)*(
	GLVariable *yout = [y plus: [[[[f plus: f4] plus: [[f3 plus: f2] scalarMultiply: 2.0]] scalarMultiply: 1/6.0] multiply: deltaT]];
	
	GLOperationOptimizer *optimizer = [[GLOperationOptimizer alloc] initWithTopVariables: @[y] bottomVariables: @[yout, deltaT]];
	
	// We can't use yout as the result variable, because it has a long string of irrelevant dependencies.
	GLVariable *result = [GLVariable variableOfType: yout.dataFormat withDimensions: yout.dimensions forEquation: yout.equation];
	GLAdaptiveIntegrationOperation *integrationOperation = [[GLAdaptiveIntegrationOperation alloc] initWithResult: result operand: y];
	integrationOperation.vectorBlockOperation = optimizer.unaryVectorOperationBlock;
	integrationOperation.dataBuffers = optimizer.internalDataBuffers;
	integrationOperation.stepSizeData = [[GLMemoryPool sharedMemoryPool] dataWithLength: deltaT.data.length];
	integrationOperation.lastStepSizeData = integrationOperation.stepSizeData;
	integrationOperation.exitOnBlowUp = YES;
	
	[deltaT solve];
	integrationOperation.stepSize = *(deltaT.pointerValue);
	
	unaryVectorOperation vectorBlock = integrationOperation.vectorBlockOperation;
	NSMutableData *deltaTData = integrationOperation.stepSizeData;
	integrationOperation.blockOperation = ^(NSMutableData *result, NSData *operand) {
		vectorBlock( @[result, deltaTData], @[operand]);
	};
	
	if (!integrationOperation.blockOperation) {
		return nil;
	}
	
	return integrationOperation;
}

+ (id) rungeKutta23AdvanceY: (GLVariable *) y stepSize: (GLFloat) deltaT fFromY: (FfromY) fFromY
{
	//	GLFloat a2=1./5., a3=3./10., a4=3./5., a5=1., a6=7./8.;
	GLFloat b21=1./2.;
	GLFloat b32=3./4.;
	GLFloat b41=2./9., b42=1./3., b43=4./9.;
	GLFloat c1=2./9., c2=1./3., c3=4./9., c4=0.;
	GLFloat dc1=c1-7./24., dc2 = c2 - 1./4., dc3=c3-1./3., dc4=c4-1./8.;
	
	GLVariable *timeStep = [GLVariable variableOfRealTypeWithDimensions: @[] forEquation: y.equation];
	*(timeStep.pointerValue) = deltaT;
	GLVariable *lastTimeStep = [GLVariable variableOfRealTypeWithDimensions: @[] forEquation: y.equation];
	*(lastTimeStep.pointerValue) = 0;
	GLVariable *relativeTolerance = [GLVariable variableOfRealTypeWithDimensions: @[] forEquation: y.equation];
	*(relativeTolerance.pointerValue) = 1e-4;
	GLVariable *absoluteTolerance = [GLVariable variableOfRealTypeWithDimensions: @[] forEquation: y.equation];
	*(absoluteTolerance.pointerValue) = 1e-7;
	
	GLVariable *f = fFromY(y);
	GLVariable *f2 = fFromY([y plus: [[f scalarMultiply: b21] multiply: timeStep]]);
	GLVariable *f3 = fFromY([y plus: [[f2 scalarMultiply: b32] multiply: timeStep]]);
	GLVariable *f4 = fFromY([y plus: [[[[f scalarMultiply: b41] plus: [f2 scalarMultiply: b42]] plus: [f3 scalarMultiply: b43]] multiply: timeStep]]);

	GLVariable *yout = [y plus: [[[[f scalarMultiply: c1] plus: [f2 scalarMultiply: c2]] plus: [f3 scalarMultiply: c3]] multiply: timeStep]];
	GLVariable *yerr = [[[[f scalarMultiply: dc1] plus: [f2 scalarMultiply: dc2]] plus: [[f3 scalarMultiply: dc3] plus: [f4 scalarMultiply: dc4]]] multiply: timeStep];
	
	// *element-wise*
	// |y_err|/max( relTolerance*|y|, absoluteTolerance ) <= 1
	GLAbsoluteValueOperation *absRelErr = [[GLAbsoluteValueOperation alloc] initWithOperand: [yout multiply: relativeTolerance] shouldUseComplexArithmetic: NO];
	GLAbsoluteValueOperation *absYErr = [[GLAbsoluteValueOperation alloc] initWithOperand: yerr shouldUseComplexArithmetic: NO];
	GLDivisionOperation * op = [[GLDivisionOperation alloc] initWithFirstOperand: absYErr.result secondOperand: [absRelErr.result absMax: absoluteTolerance] shouldUseComplexArithmetic: NO];
	GLVariable *error = op.result;
	
	GLOperationOptimizer *optimizer = [[GLOperationOptimizer alloc] initWithTopVariables: @[y, timeStep, relativeTolerance, absoluteTolerance] bottomVariables: @[yout, error]];
	
	// We can't use yout as the result variable, because it has a long string of irrelevant dependencies.
	GLVariable *result = [GLVariable variableOfType: yout.dataFormat withDimensions: yout.dimensions forEquation: yout.equation];
	GLAdaptiveIntegrationOperation *integrationOperation = [[GLAdaptiveIntegrationOperation alloc] initWithResult: result operand: y];
	
	integrationOperation.vectorBlockOperation = optimizer.unaryVectorOperationBlock;
	integrationOperation.dataBuffers = optimizer.internalDataBuffers;
	integrationOperation.errorVector = error;
	integrationOperation.stepSizeData = timeStep.data;
	integrationOperation.lastStepSizeData = lastTimeStep.data;
	integrationOperation.relativeToleranceData = relativeTolerance.data;
	integrationOperation.absoluteToleranceData = absoluteTolerance.data;
	integrationOperation.currentTime = 0;
	integrationOperation.initialTime = 0;
	integrationOperation.exitOnBlowUp = YES;
	
	NSMutableData *stepSizeData = integrationOperation.stepSizeData;
	NSMutableData *lastStepSizeData = integrationOperation.lastStepSizeData;
	NSData *errorVectorData = integrationOperation.errorVector.data;
	unaryVectorOperation vectorBlock = integrationOperation.vectorBlockOperation;
	NSData *relativeToleranceData = integrationOperation.relativeToleranceData;
	NSData *absoluteToleranceData = integrationOperation.absoluteToleranceData;
	NSUInteger nDataElements = integrationOperation.errorVector.nDataElements;
	
	integrationOperation.blockOperation = ^(NSMutableData *result, NSData *operand){
		NSMutableArray *operandBuffer = [[NSMutableArray alloc] initWithCapacity: 4];
		[operandBuffer addObject: operand];
		[operandBuffer addObject: stepSizeData];
		[operandBuffer addObject: relativeToleranceData];
		[operandBuffer addObject: absoluteToleranceData];
		
		NSMutableArray *resultBuffer = [[NSMutableArray alloc] initWithCapacity: 2];
		[resultBuffer addObject: result];
		[resultBuffer addObject: errorVectorData];
		
		GLFloat error;
		GLFloat *step = (GLFloat  *) stepSizeData.mutableBytes;
		GLFloat *lastStep = (GLFloat *) lastStepSizeData.mutableBytes;
		
		GLFloat order = 2;
		GLFloat safety = 0.9;
		GLFloat errcon = pow(5./safety, -(order+1));
		while (1) {
			// Time step
			vectorBlock( resultBuffer, operandBuffer );
			
			// Find the max error
			vGL_maxv( (float *) errorVectorData.bytes, 1, &error, nDataElements);
			
			if (error <= 1.0) {
				break;
			} else {
				GLFloat htemp = safety* (*step) * pow(error, -1./order);
				(*step) = 0.1 * (*step) >= htemp ? 0.1 * (*step) : htemp;
			}
		}
		
		(*lastStep) = (*step);
		
		if ( error > errcon) {
			(*step) = safety*(*step)*pow(error, -1./(order+1));
		} else {
			(*step) = 5.0 * (*step);
		}
	};
	
	return integrationOperation;
}

+ (id) rungeKutta45AdvanceY: (GLVariable *) y stepSize: (GLFloat) deltaT fFromY: (FfromY) fFromY
{
	//	GLFloat a2=1./5., a3=3./10., a4=3./5., a5=1., a6=7./8.;
	GLFloat b21=1./5.;
	GLFloat b31=3./40., b32=9./40.;
	GLFloat b41=3./10., b42=-9./10., b43=6./5.;
	GLFloat b51=-11./54., b52=5./2., b53=-70./27., b54=35./27.;
	GLFloat b61=1631./55296., b62=175./512., b63=575./13824., b64=44275./110592., b65=253./4096.;
	GLFloat c1=37./378., c3=250./621., c4=125./594., c6=512./1771.;
	GLFloat dc1=c1-2825./27648., dc3=c3-18575./48384., dc4=c4-13525./55296., dc5=-277./14336., dc6=c6-1./4.;
	
	GLVariable *timeStep = [GLVariable variableOfRealTypeWithDimensions: [NSArray array] forEquation: y.equation];
	*(timeStep.pointerValue) = (double) deltaT;
	GLVariable *lastTimeStep = [GLVariable variableOfRealTypeWithDimensions: [NSArray array] forEquation: y.equation];
	*(lastTimeStep.pointerValue) = 0;
	GLVariable *relativeTolerance = [GLVariable variableOfRealTypeWithDimensions: [NSArray array] forEquation: y.equation];
	*(relativeTolerance.pointerValue) = 1e-4;
	GLVariable *absoluteTolerance = [GLVariable variableOfRealTypeWithDimensions: [NSArray array] forEquation: y.equation];
	*(absoluteTolerance.pointerValue) = 1e-7;
	
	GLVariable *f = fFromY(y);
	GLVariable *f2 = fFromY([y plus: [[f scalarMultiply: b21] multiply: timeStep]]);
	GLVariable *f3 = fFromY([y plus: [[[f scalarMultiply: b31] plus: [f2 scalarMultiply: b32]] multiply: timeStep]]);	
	GLVariable *f4 = fFromY([y plus: [[[[f scalarMultiply: b41] plus: [f2 scalarMultiply: b42]] plus: [f3 scalarMultiply: b43]] multiply: timeStep]]);
	GLVariable *f5 = fFromY([y plus: [[[[f scalarMultiply: b51] plus: [f2 scalarMultiply: b52]] plus: [[f3 scalarMultiply: b53] plus: [f4 scalarMultiply: b54]]] multiply: timeStep]]);
	GLVariable *f6 = fFromY([y plus: [[[[[f scalarMultiply: b61] plus: [f2 scalarMultiply: b62]] plus: [[f3 scalarMultiply: b63] plus: [f4 scalarMultiply: b64]]] plus: [f5 scalarMultiply: b65]] multiply: timeStep]]);	
	GLVariable *yout = [y plus: [[[[f scalarMultiply: c1] plus: [f3 scalarMultiply: c3]] plus: [[f4 scalarMultiply: c4] plus: [f6 scalarMultiply: c6]]] multiply: timeStep]];
	GLVariable *yerr = [[[[[f scalarMultiply: dc1] plus: [f3 scalarMultiply: dc3]] plus: [[f4 scalarMultiply: dc4] plus: [f6 scalarMultiply: dc6]]] plus: [f5 scalarMultiply: dc5]] multiply: timeStep];
	
	// *element-wise*
	// |y_err|/max( relTolerance*|y|, absoluteTolerance ) <= 1
	GLAbsoluteValueOperation *absRelErr = [[GLAbsoluteValueOperation alloc] initWithOperand: [yout multiply: relativeTolerance] shouldUseComplexArithmetic: NO];
	GLAbsoluteValueOperation *absYErr = [[GLAbsoluteValueOperation alloc] initWithOperand: yerr shouldUseComplexArithmetic: NO];
	GLDivisionOperation * op = [[GLDivisionOperation alloc] initWithFirstOperand: absYErr.result secondOperand: [absRelErr.result absMax: absoluteTolerance] shouldUseComplexArithmetic: NO];
	GLVariable *error = op.result;
	
	GLOperationOptimizer *optimizer = [[GLOperationOptimizer alloc] initWithTopVariables: [NSArray arrayWithObjects: y, timeStep, relativeTolerance, absoluteTolerance,  nil] bottomVariables: [NSArray arrayWithObjects: yout, error, nil]];
	
	// We can't use yout as the result variable, because it has a long string of irrelevant dependencies.
	GLVariable *result = [GLVariable variableOfType: yout.dataFormat withDimensions: yout.dimensions forEquation: yout.equation];	
	GLAdaptiveIntegrationOperation *integrationOperation = [[GLAdaptiveIntegrationOperation alloc] initWithResult: result operand: y];
	
	integrationOperation.vectorBlockOperation = optimizer.unaryVectorOperationBlock;
	integrationOperation.dataBuffers = optimizer.internalDataBuffers;
	integrationOperation.errorVector = error;
	integrationOperation.stepSizeData = timeStep.data;
	integrationOperation.lastStepSizeData = lastTimeStep.data;
	integrationOperation.relativeToleranceData = relativeTolerance.data;
	integrationOperation.absoluteToleranceData = absoluteTolerance.data;
	integrationOperation.currentTime = 0;
	integrationOperation.initialTime = 0;
	integrationOperation.exitOnBlowUp = YES;
	
	NSMutableData *stepSizeData = integrationOperation.stepSizeData;
	NSMutableData *lastStepSizeData = integrationOperation.lastStepSizeData;
	NSData *errorVectorData = integrationOperation.errorVector.data;
	unaryVectorOperation vectorBlock = integrationOperation.vectorBlockOperation;
	NSData *relativeToleranceData = integrationOperation.relativeToleranceData;
	NSData *absoluteToleranceData = integrationOperation.absoluteToleranceData;
	NSUInteger nDataElements = integrationOperation.errorVector.nDataElements;
	
	integrationOperation.blockOperation = ^(NSMutableData *result, NSData *operand){ 
		NSMutableArray *operandBuffer = [[NSMutableArray alloc] initWithCapacity: 4];
		[operandBuffer addObject: operand];
		[operandBuffer addObject: stepSizeData];
		[operandBuffer addObject: relativeToleranceData];
		[operandBuffer addObject: absoluteToleranceData];
		
		NSMutableArray *resultBuffer = [[NSMutableArray alloc] initWithCapacity: 2];
		[resultBuffer addObject: result];
		[resultBuffer addObject: errorVectorData];
		
		GLFloat error;
		GLFloat *step = (GLFloat *) stepSizeData.mutableBytes;
		GLFloat *lastStep = (GLFloat *) lastStepSizeData.mutableBytes;
		
		while (1) {
			// Time step
			vectorBlock( resultBuffer, operandBuffer );
			
			// Find the max error
			vGL_maxv( (float *) errorVectorData.bytes, 1, &error, nDataElements);
			
			if (error <= 1.0) {
				break;
			} else {
				GLFloat htemp = 0.9* (*step) * pow(error, -0.25);
				(*step) = 0.1 * (*step) >= htemp ? 0.1 * (*step) : htemp;
			}
		}
		
		(*lastStep) = (*step);
		
		if ( error > 1.89e-4) {
			(*step) = 0.9*(*step)*pow(error, -0.2);
		} else {
			(*step) = 5.0 * (*step);
		}
	};
	
	return integrationOperation;
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

- (double) lastStepSize
{
	double *a = (double *) self.lastStepSizeData.mutableBytes;
	return *a;
}

- (void) setRelativeTolerance: (GLFloat)size
{
	GLFloat *a = (GLFloat *) self.relativeToleranceData.mutableBytes;
	*a = size;
}

- (GLFloat) relativeTolerance
{
	GLFloat *a = (GLFloat *) self.relativeToleranceData.mutableBytes;
	return *a;
}

- (void) setAbsoluteTolerance: (GLFloat)size
{
	GLFloat *a = (GLFloat *) self.absoluteToleranceData.mutableBytes;
	*a = size;
}

- (GLFloat) absoluteTolerance
{
	GLFloat *a = (GLFloat *) self.absoluteToleranceData.mutableBytes;
	return *a;
}

@end

/************************************************/
/*		GLVectorIntegrationOperation			*/
/************************************************/

#pragma mark -
#pragma mark GLVectorIntegrationOperation
#pragma mark

@interface GLVectorIntegrationOperation ()
@property(strong) NSArray *dataBuffers;
@property GLFloat stepSize;
@property GLFloat currentTime;
@property GLFloat lastStepSize;
@property NSMutableData *currentTimeData;
@end

@implementation GLVectorIntegrationOperation
{
	GLFloat _initialTime;
}
@synthesize stepSize;
@synthesize lastStepSize;
@synthesize currentTime;
@synthesize initialTime=_initialTime;
@synthesize exitOnBlowUp;

+ (id) rungeKutta4AdvanceY: (NSArray *) y stepSize: (GLFloat) deltaT fFromY: (FfromYVector) fFromY
{
	NSUInteger num = y.count;
	
	NSArray *f = fFromY(y);
	
	// Half a step forward and find the new slope (f) at this point
	NSMutableArray *y2 = [[NSMutableArray alloc] initWithCapacity: num];
	for (NSUInteger i=0; i < num; i++) {
		[y2 addObject: [[[f objectAtIndex: i] scalarMultiply: 0.5*deltaT] plus: [y objectAtIndex: i]]];
	}
	NSArray *f2 = fFromY( y2 );
	
	// Find the new slope at this 2nd point
	NSMutableArray *y3 = [[NSMutableArray alloc] initWithCapacity: num];
	for (NSUInteger i=0; i < num; i++) {
		[y3 addObject: [[[f2 objectAtIndex: i] scalarMultiply: 0.5*deltaT] plus: [y objectAtIndex: i]]];
	}
	NSArray *f3 = fFromY( y3 );
	
	// Go a full step forward with the 2nd new slope
	NSMutableArray *y4 = [[NSMutableArray alloc] initWithCapacity: num];
	for (NSUInteger i=0; i < num; i++) {
		[y4 addObject: [[[f3 objectAtIndex: i] scalarMultiply: deltaT] plus: [y objectAtIndex: i]]];
	}
	NSArray *f4 = fFromY( y4 );
	
	NSMutableArray *yout = [[NSMutableArray alloc] initWithCapacity: num];
	for (NSUInteger i=0; i < num; i++) {
		[yout addObject: [[y objectAtIndex: i] plus: [[[[f objectAtIndex: i] plus: [f4 objectAtIndex: i]] plus: [[[f3 objectAtIndex: i] plus: [f2 objectAtIndex: i]] scalarMultiply: 2.0]] scalarMultiply: deltaT/6.0]]];
	}
	
	GLOperationOptimizer *optimizer = [[GLOperationOptimizer alloc] initWithTopVariables: y bottomVariables: yout];
	
	NSMutableArray *result = [[NSMutableArray alloc] initWithCapacity: yout.count];
	for (GLVariable *variable in yout) {
		[result addObject: [GLVariable variableOfType: variable.dataFormat withDimensions: variable.dimensions forEquation: variable.equation]];
	}
	
	GLVectorIntegrationOperation *integrationOperation = [[GLVectorIntegrationOperation alloc] initWithResult: result operand: y];
	integrationOperation.blockOperation = optimizer.unaryVectorOperationBlock;
	integrationOperation.dataBuffers = optimizer.internalDataBuffers;
	integrationOperation.stepSize = deltaT;
	integrationOperation.currentTime = 0;
	integrationOperation.initialTime = 0;
	integrationOperation.lastStepSize = 0;
	integrationOperation.exitOnBlowUp = YES;
	
	if (!integrationOperation.blockOperation) {
		return nil;
	}
	
	return integrationOperation;
}

+ (id) rungeKutta4AdvanceY: (NSArray *) y stepSize: (GLFloat) deltaT fFromTY: (FfromTYVector) fFromY
{
	NSUInteger num = y.count;
	GLVariable *t1 = [GLVariable variableOfRealTypeWithDimensions: @[] forEquation: [y[0] equation]];
	GLVariable *t2 = [t1 scalarAdd: 0.5*deltaT];
	GLVariable *t4 = [t1 scalarAdd: deltaT];
	
	NSArray *f = fFromY( t1, y);
	
	// Half a step forward and find the new slope (f) at this point
	NSMutableArray *y2 = [[NSMutableArray alloc] initWithCapacity: num];
	for (NSUInteger i=0; i < num; i++) {
		[y2 addObject: [[[f objectAtIndex: i] scalarMultiply: 0.5*deltaT] plus: [y objectAtIndex: i]]];
	}
	NSArray *f2 = fFromY( t2, y2 );
	
	// Find the new slope at this 2nd point
	NSMutableArray *y3 = [[NSMutableArray alloc] initWithCapacity: num];
	for (NSUInteger i=0; i < num; i++) {
		[y3 addObject: [[[f2 objectAtIndex: i] scalarMultiply: 0.5*deltaT] plus: [y objectAtIndex: i]]];
	}
	NSArray *f3 = fFromY( t2, y3 );
	
	// Go a full step forward with the 2nd new slope
	NSMutableArray *y4 = [[NSMutableArray alloc] initWithCapacity: num];
	for (NSUInteger i=0; i < num; i++) {
		[y4 addObject: [[[f3 objectAtIndex: i] scalarMultiply: deltaT] plus: [y objectAtIndex: i]]];
	}
	NSArray *f4 = fFromY( t4, y4 );
	
	NSMutableArray *yout = [[NSMutableArray alloc] initWithCapacity: num];
	for (NSUInteger i=0; i < num; i++) {
		[yout addObject: [[y objectAtIndex: i] plus: [[[[f objectAtIndex: i] plus: [f4 objectAtIndex: i]] plus: [[[f3 objectAtIndex: i] plus: [f2 objectAtIndex: i]] scalarMultiply: 2.0]] scalarMultiply: deltaT/6.0]]];
	}
	
	NSMutableArray *yfull = [NSMutableArray arrayWithArray: y];
	[yfull addObject: t1];
	GLOperationOptimizer *optimizer = [[GLOperationOptimizer alloc] initWithTopVariables: yfull bottomVariables: yout];
	
	NSMutableArray *result = [[NSMutableArray alloc] initWithCapacity: yout.count];
	for (GLVariable *variable in yout) {
		[result addObject: [GLVariable variableOfType: variable.dataFormat withDimensions: variable.dimensions forEquation: variable.equation]];
	}
	
	GLVectorIntegrationOperation *integrationOperation = [[GLVectorIntegrationOperation alloc] initWithResult: result operand: y];
	integrationOperation.dataBuffers = optimizer.internalDataBuffers;
	integrationOperation.stepSize = deltaT;
	integrationOperation.currentTime = 0;
	integrationOperation.initialTime = 0;
	integrationOperation.lastStepSize = 0;
	integrationOperation.exitOnBlowUp = YES;
	integrationOperation.currentTimeData = t1.data;
	
	unaryVectorOperation vectorBlock = optimizer.unaryVectorOperationBlock;
	NSMutableData *tData = t1.data;
	NSMutableArray *operandArray = [NSMutableArray arrayWithCapacity: yfull.count];
	integrationOperation.blockOperation = ^(NSArray *result, NSArray *operand) {
		[operandArray addObjectsFromArray: operand];
		[operandArray addObject: tData];
		
		vectorBlock( result, operandArray);
		
		GLFloat *a = (GLFloat *) tData.mutableBytes;
		*a += deltaT;
		
		[operandArray removeAllObjects];
	};
	
	if (!integrationOperation.blockOperation) {
		return nil;
	}
	
	return integrationOperation;
}


+ (id) rungeKutta5AdvanceY: (NSArray *) y stepSize: (GLFloat) deltaT fFromY: (FfromYVector) fFromY
{
	GLFloat b21=1./5.;
	GLFloat b31=3./40., b32=9./40.;
	GLFloat b41=3./10., b42=-9./10., b43=6./5.;
	GLFloat b51=-11./54., b52=5./2., b53=-70./27., b54=35./27.;
	GLFloat b61=1631./55296., b62=175./512., b63=575./13824., b64=44275./110592., b65=253./4096.;
	GLFloat c1=37./378., c3=250./621., c4=125./594., c6=512./1771.;
	
	NSUInteger num = y.count;
	
	NSArray *f = fFromY(y);
	
	NSMutableArray *y2 = [[NSMutableArray alloc] initWithCapacity: num];
	for (NSUInteger i=0; i < num; i++) {
		[y2 addObject: [y[i] plus: [f[i] scalarMultiply: b21*deltaT]]];
	}
	NSArray *f2 = fFromY( y2 );
	
	NSMutableArray *y3 = [[NSMutableArray alloc] initWithCapacity: num];
	for (NSUInteger i=0; i < num; i++) {
		[y3 addObject: [y[i] plus: [[[f[i] scalarMultiply: b31] plus: [f2[i] scalarMultiply: b32]] scalarMultiply: deltaT]]];
	}
	NSArray *f3 = fFromY( y3 );
	
	NSMutableArray *y4 = [[NSMutableArray alloc] initWithCapacity: num];
	for (NSUInteger i=0; i < num; i++) {
		[y4 addObject: [y[i] plus: [[[[f[i] scalarMultiply: b41] plus: [f2[i] scalarMultiply: b42]] plus: [f3[i] scalarMultiply: b43]] scalarMultiply: deltaT]]];
	}
	NSArray *f4 = fFromY( y4 );
	
	NSMutableArray *y5 = [[NSMutableArray alloc] initWithCapacity: num];
	for (NSUInteger i=0; i < num; i++) {
		[y5 addObject: [y[i] plus: [[[[f[i] scalarMultiply: b51] plus: [f2[i] scalarMultiply: b52]] plus: [[f3[i] scalarMultiply: b53] plus: [f4[i] scalarMultiply: b54]]] scalarMultiply: deltaT]]];
	}
	NSArray *f5 = fFromY( y5 );
	
	NSMutableArray *y6 = [[NSMutableArray alloc] initWithCapacity: num];
	for (NSUInteger i=0; i < num; i++) {
		[y6 addObject: [y[i] plus: [[[[[f[i] scalarMultiply: b61] plus: [f2[i] scalarMultiply: b62]] plus: [[f3[i] scalarMultiply: b63] plus: [f4[i] scalarMultiply: b64]]] plus: [f5[i] scalarMultiply: b65]] scalarMultiply: deltaT]]];
	}
	NSArray *f6 = fFromY( y6 );
	
	NSMutableArray *yout = [[NSMutableArray alloc] initWithCapacity: num];
	for (NSUInteger i=0; i < num; i++) {
		[yout addObject: [y[i] plus: [[[[f[i] scalarMultiply: c1] plus: [f3[i] scalarMultiply: c3]] plus: [[f4[i] scalarMultiply: c4] plus: [f6[i] scalarMultiply: c6]]] scalarMultiply: deltaT]]];
	}
	
	GLOperationOptimizer *optimizer = [[GLOperationOptimizer alloc] initWithTopVariables: y bottomVariables: yout];
	
	NSMutableArray *result = [[NSMutableArray alloc] initWithCapacity: yout.count];
	for (GLVariable *variable in yout) {
		[result addObject: [GLVariable variableOfType: variable.dataFormat withDimensions: variable.dimensions forEquation: variable.equation]];
	}
	
	GLVectorIntegrationOperation *integrationOperation = [[GLVectorIntegrationOperation alloc] initWithResult: result operand: y];
	integrationOperation.blockOperation = optimizer.unaryVectorOperationBlock;
	integrationOperation.dataBuffers = optimizer.internalDataBuffers;
	integrationOperation.stepSize = deltaT;
	integrationOperation.currentTime = 0;
	integrationOperation.initialTime = 0;
	integrationOperation.lastStepSize = 0;
	integrationOperation.exitOnBlowUp = YES;
	
	if (!integrationOperation.blockOperation) {
		return nil;
	}
	
	return integrationOperation;
}

- (void) setInitialTime:(GLFloat)initialTime
{
	self.currentTime = initialTime + (self.currentTime-_initialTime);
	_initialTime=initialTime;
}

- (NSArray *) stepForward: (NSArray *) y
{
	[y makeObjectsPerformSelector: @selector(solve)];
	self.lastStepSize = self.stepSize;
	NSMutableArray *yout = [[NSMutableArray alloc] initWithCapacity: self.result.count];
	NSMutableArray *resultBuffer = [[NSMutableArray alloc] initWithCapacity: self.result.count];
	NSMutableArray *operandBuffer = [[NSMutableArray alloc] initWithCapacity: self.operand.count];
	
	for (GLVariable *variable in self.result) {
		GLVariable *newVar = [GLVariable variableOfType: variable.dataFormat withDimensions: variable.dimensions forEquation: variable.equation];
		[yout addObject: newVar];
		[resultBuffer addObject: newVar.data];
	}
	for (GLVariable *variable in y) {
		[operandBuffer addObject: variable.data];
	}
	self.blockOperation( resultBuffer, operandBuffer );
	self.currentTime = self.currentTime + self.stepSize;
	
//	if (self.exitOnBlowUp) {
//		for (GLVariable *aVariable in yout) {
//			if (![aVariable isFinite]) {
//				NSLog(@"Blow up detected! Exiting...");
//				return nil;
//			}
//		}
//	}
	
	return yout;
}

- (NSArray *) stepForward: (NSArray *) y toTime: (GLFloat) time
{
	self.lastStepSize = self.stepSize;
	// If we've already exceeded the time, just return the same variable.
	if (self.currentTime >= time) {
		return y;
	}
	
	// If we haven't, then take one time step forward...
	[y makeObjectsPerformSelector: @selector(solve)];
	NSMutableArray *yout = [[NSMutableArray alloc] initWithCapacity: self.result.count];
	NSMutableArray *resultBuffer = [[NSMutableArray alloc] initWithCapacity: self.result.count];
	NSMutableArray *operandBuffer = [[NSMutableArray alloc] initWithCapacity: self.operand.count];
	for (GLVariable *variable in self.result) {
		GLVariable *newVar = [GLVariable variableOfType: variable.dataFormat withDimensions: variable.dimensions forEquation: variable.equation];
		[yout addObject: newVar];
		[resultBuffer addObject: newVar.data];
	}
	for (GLVariable *variable in y) {
		[operandBuffer addObject: variable.data];
	}
	self.blockOperation( resultBuffer, operandBuffer );
	self.currentTime = self.currentTime + self.stepSize;
    
//	if (self.exitOnBlowUp) {
//		for (GLVariable *aVariable in yout) {
//			if (![aVariable isFinite]) {
//				NSLog(@"Blow up detected! Exiting...");
//				return nil;
//			}
//		}
//	}
	
	// ...and as many as necessary
	while ( self.currentTime < time - self.stepSize/2) {
        
        // Last iterations output is the new input
		self.blockOperation( resultBuffer, resultBuffer );
		self.currentTime = self.currentTime + self.stepSize;
		
//		if (self.exitOnBlowUp) {
//			for (GLVariable *aVariable in yout) {
//				if (![aVariable isFinite]) {
//					NSLog(@"Blow up detected! Exiting...");
//					return nil;
//				}
//			}
//		}
	}
	
	return yout;
}

@synthesize dataBuffers;

- (void) dealloc
{
	if (self.dataBuffers) {
		for (NSMutableData *data in self.dataBuffers) {
			[[GLMemoryPool sharedMemoryPool] returnData: data];
		}
	}
}

@end


/************************************************/
/*		GLAdaptiveVectorIntegrationOperation	*/
/************************************************/

#pragma mark -
#pragma mark GLAdaptiveVectorIntegrationOperation
#pragma mark

@interface GLAdaptiveVectorIntegrationOperation ()
@property NSMutableData *stepSizeData;
@property NSMutableData *lastStepSizeData;
@property NSArray *relativeToleranceData;
@property NSArray *absoluteToleranceData;
@property NSArray *errorVector;
@end

@implementation GLAdaptiveVectorIntegrationOperation
{
	NSMutableData *_stepSizeData;
	NSMutableData *_lastStepSizeData;
	NSArray *_relativeToleranceData;
	NSArray *_absoluteToleranceData;
	NSArray *_errorVector;
}

@synthesize stepSizeData=_stepSizeData;
@synthesize lastStepSizeData=_lastStepSizeData;
@synthesize relativeToleranceData=_relativeToleranceData;
@synthesize absoluteToleranceData=_absoluteToleranceData;
@synthesize errorVector=_errorVector;

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

+ (id) rungeKutta45AdvanceY: (NSArray *) y stepSize: (GLFloat) deltaT fFromY: (FfromYVector) fFromY
{
	//	GLFloat a2=1./5., a3=3./10., a4=3./5., a5=1., a6=7./8.;
	GLFloat b21=1./5.;
	GLFloat b31=3./40., b32=9./40.;
	GLFloat b41=3./10., b42=-9./10., b43=6./5.;
	GLFloat b51=-11./54., b52=5./2., b53=-70./27., b54=35./27.;
	GLFloat b61=1631./55296., b62=175./512., b63=575./13824., b64=44275./110592., b65=253./4096.;
	GLFloat c1=37./378., c3=250./621., c4=125./594., c6=512./1771.;
	GLFloat dc1=c1-2825./27648., dc3=c3-18575./48384., dc4=c4-13525./55296., dc5=-277./14336., dc6=c6-1./4.;
	
	GLEquation *equation = [[y lastObject] equation];
	
	GLVariable *timeStep = [GLVariable variableOfRealTypeWithDimensions: [NSArray array] forEquation: equation];
	*(timeStep.pointerValue) = deltaT;
	GLVariable *lastTimeStep = [GLVariable variableOfRealTypeWithDimensions: [NSArray array] forEquation: equation];
	*(lastTimeStep.pointerValue) = 0;
	
	NSMutableArray *relativeToleranceArray = [[NSMutableArray alloc] initWithCapacity: y.count];
	NSMutableArray *absoluteToleranceArray = [[NSMutableArray alloc] initWithCapacity: y.count];
	NSMutableArray *relativeToleranceDataArray = [[NSMutableArray alloc] initWithCapacity: y.count];
	NSMutableArray *absoluteToleranceDataArray = [[NSMutableArray alloc] initWithCapacity: y.count];
	for (NSUInteger i=0; i<y.count; i++) {
		GLVariable *relativeTolerance = [GLVariable variableOfRealTypeWithDimensions: [NSArray array] forEquation: equation];
		*(relativeTolerance.pointerValue) = 1e-4;
		GLVariable *absoluteTolerance = [GLVariable variableOfRealTypeWithDimensions: [NSArray array] forEquation: equation];
		*(absoluteTolerance.pointerValue) = 1e-7;
		[relativeToleranceArray addObject: relativeTolerance];
		[absoluteToleranceArray addObject: absoluteTolerance];
		[relativeToleranceDataArray addObject: relativeTolerance.data];
		[absoluteToleranceDataArray addObject: absoluteTolerance.data];
	}
	
	NSUInteger num = y.count;
	
	NSArray *f = fFromY(y);
	
	NSMutableArray *y2 = [[NSMutableArray alloc] initWithCapacity: num];
	for (NSUInteger i=0; i < num; i++) {
		[y2 addObject: [[y objectAtIndex: i] plus: [[[f objectAtIndex: i] scalarMultiply: b21] multiply: timeStep]]];
	}
	NSArray *f2 = fFromY( y2 );
	
	NSMutableArray *y3 = [[NSMutableArray alloc] initWithCapacity: num];
	for (NSUInteger i=0; i < num; i++) {
		[y3 addObject: [[y objectAtIndex: i] plus: [[[[f objectAtIndex: i] scalarMultiply: b31] plus: [[f2 objectAtIndex: i] scalarMultiply: b32]] multiply: timeStep]]];
	}
	NSArray *f3 = fFromY( y3 );
	
	NSMutableArray *y4 = [[NSMutableArray alloc] initWithCapacity: num];
	for (NSUInteger i=0; i < num; i++) {
		[y4 addObject: [[y objectAtIndex: i] plus: [[[[[f objectAtIndex: i] scalarMultiply: b41] plus: [[f2 objectAtIndex: i] scalarMultiply: b42]] plus: [[f3 objectAtIndex: i] scalarMultiply: b43]] multiply: timeStep]]];
	}
	NSArray *f4 = fFromY( y4 );
	
	NSMutableArray *y5 = [[NSMutableArray alloc] initWithCapacity: num];
	for (NSUInteger i=0; i < num; i++) {
		[y5 addObject: [[y objectAtIndex: i] plus: [[[[[f objectAtIndex: i] scalarMultiply: b51] plus: [[f2 objectAtIndex: i] scalarMultiply: b52]] plus: [[[f3 objectAtIndex: i] scalarMultiply: b53] plus: [[f4 objectAtIndex: i] scalarMultiply: b54]]] multiply: timeStep]]];
	}
	NSArray *f5 = fFromY( y5 );
	
	NSMutableArray *y6 = [[NSMutableArray alloc] initWithCapacity: num];
	for (NSUInteger i=0; i < num; i++) {
		[y6 addObject: [[y objectAtIndex: i] plus: [[[[[[f objectAtIndex: i] scalarMultiply: b61] plus: [[f2 objectAtIndex: i] scalarMultiply: b62]] plus: [[[f3 objectAtIndex: i] scalarMultiply: b63] plus: [[f4 objectAtIndex: i] scalarMultiply: b64]]] plus: [[f5 objectAtIndex: i] scalarMultiply: b65]] multiply: timeStep]]];
	}
	NSArray *f6 = fFromY( y6 );
	
	NSMutableArray *yout = [[NSMutableArray alloc] initWithCapacity: num];
	for (NSUInteger i=0; i < num; i++) {
		[yout addObject: [[y objectAtIndex: i] plus: [[[[[f objectAtIndex: i] scalarMultiply: c1] plus: [[f3 objectAtIndex: i] scalarMultiply: c3]] plus: [[[f4 objectAtIndex: i] scalarMultiply: c4] plus: [[f6 objectAtIndex: i] scalarMultiply: c6]]] multiply: timeStep]]];
	}
	
	NSMutableArray *yerr = [[NSMutableArray alloc] initWithCapacity: num];
	for (NSUInteger i=0; i < num; i++) {
		[yerr addObject: [[[[[[f objectAtIndex: i] scalarMultiply: dc1] plus: [[f3 objectAtIndex: i] scalarMultiply: dc3]] plus: [[[f4 objectAtIndex: i] scalarMultiply: dc4] plus: [[f6 objectAtIndex: i] scalarMultiply: dc6]]] plus: [[f5 objectAtIndex: i] scalarMultiply: dc5]] multiply: timeStep]];
	}
	
	// *element-wise*
	// |y_err|/max( relTolerance*|y|, absoluteTolerance ) <= 1
	NSMutableArray *error = [[NSMutableArray alloc] initWithCapacity: num];
	NSMutableArray *errorVectorData = [[NSMutableArray alloc] initWithCapacity: num];
	for (NSUInteger i=0; i < num; i++) {
		GLAbsoluteValueOperation *absRelErr = [[GLAbsoluteValueOperation alloc] initWithOperand: [[yout objectAtIndex: i] multiply: [relativeToleranceArray objectAtIndex: i]] shouldUseComplexArithmetic: NO];
		GLAbsoluteValueOperation *absYErr = [[GLAbsoluteValueOperation alloc] initWithOperand: [yerr objectAtIndex: i] shouldUseComplexArithmetic: NO];
		GLDivisionOperation * op = [[GLDivisionOperation alloc] initWithFirstOperand: absYErr.result secondOperand: [absRelErr.result absMax: [absoluteToleranceArray objectAtIndex: i]] shouldUseComplexArithmetic: NO];
		[error addObject: op.result];
		[errorVectorData addObject: op.result.data];
	}
	
	NSMutableArray *topVariables = [[NSMutableArray alloc] init];
	[topVariables addObjectsFromArray: y];
	[topVariables addObject: timeStep];
	[topVariables addObjectsFromArray: relativeToleranceArray];
	[topVariables addObjectsFromArray: absoluteToleranceArray];
	
	NSMutableArray *bottomVariables = [[NSMutableArray alloc] init];
	[bottomVariables addObjectsFromArray: yout];
	[bottomVariables addObjectsFromArray: error];
	
	GLOperationOptimizer *optimizer = [[GLOperationOptimizer alloc] initWithTopVariables: topVariables bottomVariables: bottomVariables];
	
	// We can't use yout as the result variable, because it has a long string of irrelevant dependencies.
	NSMutableArray *result = [[NSMutableArray alloc] initWithCapacity: yout.count];
	for (GLVariable *variable in yout) {
		[result addObject: [GLVariable variableOfType: variable.dataFormat withDimensions: variable.dimensions forEquation: variable.equation]];
	}
	
	GLAdaptiveVectorIntegrationOperation *integrationOperation = [[GLAdaptiveVectorIntegrationOperation alloc] initWithResult: result operand: y];
	
	integrationOperation.dataBuffers = optimizer.internalDataBuffers;
	integrationOperation.errorVector = error;
	integrationOperation.stepSizeData = timeStep.data;
	integrationOperation.lastStepSizeData = lastTimeStep.data;
	integrationOperation.relativeToleranceData = relativeToleranceDataArray;
	integrationOperation.absoluteToleranceData = absoluteToleranceDataArray;
	integrationOperation.currentTime = 0;
	integrationOperation.initialTime = 0;
	integrationOperation.exitOnBlowUp = YES;
	
	NSMutableData *stepSizeData = integrationOperation.stepSizeData;
	NSMutableData *lastStepSizeData = integrationOperation.lastStepSizeData;
	unaryVectorOperation vectorBlock = optimizer.unaryVectorOperationBlock;
	NSMutableArray *nDataElementsArray = [[NSMutableArray alloc] initWithCapacity: error.count];
	for (NSUInteger i=0; i<error.count; i++) {
		[nDataElementsArray addObject: [NSNumber numberWithUnsignedInteger: [[error objectAtIndex: i] nDataElements]]];
	}
	
	integrationOperation.blockOperation = ^(NSArray *result, NSArray *operand){ 
		NSMutableArray *operandBuffer = [[NSMutableArray alloc] initWithCapacity: 3*num + 1];
		[operandBuffer addObjectsFromArray: operand];
		[operandBuffer addObject: stepSizeData];
		[operandBuffer addObjectsFromArray: relativeToleranceDataArray];
		[operandBuffer addObjectsFromArray: absoluteToleranceDataArray];
		
		NSMutableArray *resultBuffer = [[NSMutableArray alloc] initWithCapacity: 2*num];
		[resultBuffer addObjectsFromArray: result];
		[resultBuffer addObjectsFromArray: errorVectorData];
		
		GLFloat error;
		GLFloat *step = (GLFloat *) stepSizeData.mutableBytes;
		GLFloat *lastStep = (GLFloat *) lastStepSizeData.mutableBytes;
		
		while (1) {
			// Time step
			vectorBlock( resultBuffer, operandBuffer );
			
			// Find the max error
			for (NSUInteger i=0; i<num; i++) {
				NSData *errorData = [errorVectorData objectAtIndex: i];
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
				GLFloat htemp = 0.9* (*step) * pow(error, -0.25);
				(*step) = 0.1 * (*step) >= htemp ? 0.1 * (*step) : htemp;
			}
		}
		
		(*lastStep) = (*step);
		
		if ( error > 1.89e-4) {
			(*step) = 0.9*(*step)*pow(error, -0.2);
		} else {
			(*step) = 5.0 * (*step);
		}
	};
	
	return integrationOperation;
}


@end
















