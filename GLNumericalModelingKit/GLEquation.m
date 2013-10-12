//
//  GLEquation.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 4/15/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import "GLEquation.h"
#import "GLVariable.h"
#import "GLDimension.h"
#import "GLVariableOperations.h"
#import "GLSpectralDifferentialOperatorPool.h"
#import <fftw3.h>

/************************************************/
/*		Private Methods 						*/
/************************************************/

#pragma mark -
#pragma mark Private Methods
#pragma mark

@interface GLEquation ()

// All variables get computed in this queue
@property(strong, readwrite) NSOperationQueue *operationQueue;
@property(strong, readwrite) NSMapTable *dimensionOperatorPoolMapTable;
@property(strong, readwrite) NSMutableDictionary *basisDictionary;
@property dispatch_queue_t serialDispatchQueue;

@end

/************************************************/
/*		Init								*/
/************************************************/

#pragma mark -
#pragma mark Init
#pragma mark

@implementation GLEquation
{
	NSOperationQueue *operationQueue;
	dispatch_queue_t serialDispatchQueue;
}

- (id) init
{
	if ((self = [super init])) {
        fftwf_init_threads();
		operationQueue = [[NSOperationQueue alloc] init];
		serialDispatchQueue = dispatch_queue_create( "com.EarlyInnovations.SerialEquationQueue", NULL );
		self.dimensionOperatorPoolMapTable = [NSMapTable mapTableWithKeyOptions: (NSMapTableCopyIn | NSMapTableStrongMemory) valueOptions:(NSMapTableObjectPointerPersonality | NSMapTableStrongMemory)];
		self.basisDictionary = [NSMutableDictionary dictionary];
	}
	return self;
}

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

- (void) removeDependenciesFromOperation: (NSOperation *) operation
{
	for (NSOperation *depOp in operation.dependencies) {
		[operation removeDependency: depOp];
		[self removeDependenciesFromOperation: depOp];
	}
}

// Used to pull out all the dependencies for a particular operation.
// This allows us to only add those variables to the queue that really need to be solved for.
// If a dependency is recursively defined, we're absolutely screwed.
- (void) addDependenciesForOperation: (NSOperation *) anOperation toSet: (NSMutableSet *) aSet
{
	[aSet addObjectsFromArray: anOperation.dependencies];
	for (NSOperation *aDependency in anOperation.dependencies) {
		[self addDependenciesForOperation: aDependency toSet: aSet];
	}
}

//- (void) addUnenqueuedDependenciesForOperation: (NSOperation *) anOperation toSet: (NSMutableSet *) aSet basedOnCurrentQueue: (NSArray *) queue
//{
//	for (NSOperation *aDependency in anOperation.dependencies) {
//		if ( !(aDependency.isFinished || aDependency.isExecuting  || [queue containsObject:aDependency]) ) {
//			[aSet addObject: aDependency];
//		}
//		[self addUnenqueuedDependenciesForOperation: aDependency toSet: aSet basedOnCurrentQueue: queue];
//	}
//}

- (void) addUnenqueuedDependenciesForOperation: (GLVariableOperation *) anOperation toSet: (NSMutableSet *) aSet basedOnCurrentQueue: (NSArray *) queue
{
	for (GLVariableOperation *aDependency in anOperation.dependencies) {
		if ( !aDependency.isEnqueued ) {
			[aSet addObject: aDependency];
			aDependency.isEnqueued=YES;
			[self addUnenqueuedDependenciesForOperation: aDependency toSet: aSet basedOnCurrentQueue: queue];
		}
	}
}

// Operations automatically add their dependencies to other operations.
// When they complete, they remove themselves as 

- (void) solveForVariable: (GLTensor *) aVariable
{
	[self solveForVariable: aVariable waitUntilFinished:YES];
}

- (void) solveForVariable: (GLTensor *) aVariable waitUntilFinished: (BOOL) shouldWait
{
	if (!aVariable.pendingOperations.count) return;
	
	__block GLVariableOperation *lastOperation;
	
	dispatch_sync(self.serialDispatchQueue, ^{
		NSMutableSet *unenquedOperations = [NSMutableSet set];
		NSArray *queue = self.operationQueue.operations;
		for ( GLVariableOperation *operation in aVariable.pendingOperations ) {
			lastOperation = operation;
			if ( !operation.isEnqueued ) {
				[unenquedOperations addObject: operation];
				operation.isEnqueued=YES;
				[self addUnenqueuedDependenciesForOperation:operation toSet:unenquedOperations basedOnCurrentQueue: queue];
			}
		 }
		if (unenquedOperations.count) [self.operationQueue addOperations: [unenquedOperations allObjects] waitUntilFinished:NO];
	});
	
	if (shouldWait) [lastOperation waitUntilFinished];
}

- (void) solveForOperation: (GLVariableOperation *) anOperation waitUntilFinished: (BOOL) shouldWait
{	
	dispatch_sync(self.serialDispatchQueue, ^{
		NSMutableSet *unenquedOperations = [NSMutableSet set];
		NSArray *queue = self.operationQueue.operations;
		if ( !anOperation.isEnqueued ) {
            [unenquedOperations addObject: anOperation];
            anOperation.isEnqueued=YES;
            [self addUnenqueuedDependenciesForOperation:anOperation toSet:unenquedOperations basedOnCurrentQueue: queue];
        }
		if (unenquedOperations.count) [self.operationQueue addOperations: [unenquedOperations allObjects] waitUntilFinished:NO];
	});
    
    if (shouldWait) [anOperation waitUntilFinished];
}

- (void) waitUntilAllOperationsAreFinished
{
	[self.operationQueue waitUntilAllOperationsAreFinished];
}

/************************************************/
/*		Differentiation							*/
/************************************************/

#pragma mark -
#pragma mark Differentiation
#pragma mark

@synthesize dimensionOperatorPoolMapTable;

// Depends on whether or not the variable is real or complex, and the variables dimensions.
- (id) defaultDifferentialOperatorPoolForVariable: (GLVariable *) aVariable
{
	NSArray *differentiationBasis = aVariable.differentiationBasis;
	if (!differentiationBasis) {
		differentiationBasis = [self defaultDifferentiationBasisForOrder: aVariable.dimensions.count];
	}
	NSArray *transformedDimensions = [aVariable dimensionsTransformedToBasis: differentiationBasis];
	
	// This should do an individual isEquals on each dimension.
	GLDifferentialOperatorPool *pool = [self.dimensionOperatorPoolMapTable objectForKey: transformedDimensions];
	
	if (!pool) {
		// Presumably the user is trying to differentiate the spatial dimensions
		NSArray *diffDimensions = [aVariable dimensionsTransformedToBasis:  @[@(kGLDeltaBasis)]];
		
		GLSpectralDifferentialOperatorPool *spectralPool = [[GLSpectralDifferentialOperatorPool alloc] initWithDifferentiationDimensions: diffDimensions transformDimensions: transformedDimensions forEquation: self];
		pool = spectralPool;
		[self.dimensionOperatorPoolMapTable setObject: spectralPool forKey: transformedDimensions];
	}
	
	return pool;
	
}

- (void) setDefaultDifferentiationBasis:(NSArray *) aBasis forOrder: (NSUInteger) order
{
	[self.basisDictionary setObject: aBasis forKey: @(order)];
}

- (NSArray *) defaultDifferentiationBasisForOrder: (NSUInteger) order
{
	NSMutableArray *possibleMatch = [self.basisDictionary objectForKey: @(order)];
	if (!possibleMatch) {
		NSArray *orderOne = [self.basisDictionary objectForKey: @(1)];
		NSNumber *aBasis = orderOne.count ? orderOne.lastObject : @(kGLExponentialBasis);
		possibleMatch = [NSMutableArray arrayWithCapacity: order];
		for (NSUInteger i=0; i<order; i++) {
			[possibleMatch addObject: aBasis];
		}
	}
	return possibleMatch;
}

/************************************************/
/*		Time Stepping							*/
/************************************************/

#pragma mark -
#pragma mark Time Stepping
#pragma mark

- (GLVariable *) rungeKuttaAdvanceY: (GLVariable *) y withF: (GLVariable *) f stepSize: (GLFloat) deltaT fFromY: (FfromY) fFromY;
{
	if (!f) {
		f = fFromY(y);
	}
	
	// Half a step forward and find the new slope (f) at this point
	GLVariable *f2 = fFromY([[f scalarMultiply: 0.5*deltaT] plus: y]);
	[self solveForVariable: f2];
	
	// Find the new slope at this 2nd point
	GLVariable *f3 = fFromY([[f2 scalarMultiply: 0.5*deltaT] plus: y]);
	[self solveForVariable: f3];
	
	// Go a full step forward with the 2nd new slope
	GLVariable *f4 = fFromY([[f3 scalarMultiply: deltaT] plus: y]);
	
	// yout = y + (Delta/6)*(
	GLVariable *yout = [y plus: [[[f plus: f4] plus: [[f3 plus: f2] scalarMultiply: 2.0]] scalarMultiply: deltaT/6.0]];
	[self solveForVariable: yout];
    
	return yout;
}


- (GLVariable *) rungeKuttaAdvanceY: (GLVariable *) y stepSize: (GLFloat) deltaT fFromY: (FfromY) fFromY
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
	
	return yout;
}

- (NSArray *) rungeKuttaAdvanceYVector: (NSArray *) y stepSize: (GLFloat) deltaT fFromY: (FfromYVector) fFromY
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
	
	return yout;
}

/************************************************/
/*		Private Methods 						*/
/************************************************/

#pragma mark -
#pragma mark Private Methods
#pragma mark

@synthesize operationQueue;
@synthesize serialDispatchQueue;

@end
