//
//  GLEquation.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 4/15/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import "GLEquation.h"
#import "GLFunction.h"
#import "GLDimension.h"
#import "GLVariableOperations.h"
#import "GLLinearTransform.h"
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

static NSString *GLEquationDimensionOperatorPoolMapTableKey = @"GLEquationDimensionOperatorPoolMapTableKey";
- (void)encodeWithCoder:(NSCoder *)coder
{
    [coder encodeObject: self.dimensionOperatorPoolMapTable forKey:GLEquationDimensionOperatorPoolMapTableKey];
}

- (id)initWithCoder:(NSCoder *)decoder
{
    if ((self=[super init])) {
        fftwf_init_threads();
		operationQueue = [[NSOperationQueue alloc] init];
		serialDispatchQueue = dispatch_queue_create( "com.EarlyInnovations.SerialEquationQueue", NULL );
        self.dimensionOperatorPoolMapTable = [decoder decodeObjectForKey: GLEquationDimensionOperatorPoolMapTableKey];
        self.basisDictionary = [NSMutableDictionary dictionary];
    }
    return self;
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

- (void) solveForVariable: (GLVariable *) aVariable
{
	[self solveForVariable: aVariable waitUntilFinished:YES];
}

- (void) solveForVariable: (GLVariable *) aVariable waitUntilFinished: (BOOL) shouldWait
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

- (GLLinearTransform *) linearTransformWithName: (NSString *) name forDimensions: (NSArray *) dimensions
{
	NSMutableDictionary *pool = [self.dimensionOperatorPoolMapTable objectForKey: dimensions];
	if (pool) {
		return [pool objectForKey: name];
	} else {
		return nil;
	}
}

- (void) setLinearTransform: (GLLinearTransform *) diffOp withName: (NSString *) name
{
	[diffOp solve];
	
	BOOL containsNan = NO;
	for (NSUInteger i=0; i<diffOp.nDataElements; i++) {
		if ( !isfinite(diffOp.pointerValue[i])) containsNan = YES;
	}
	
	if ( containsNan ) {
		NSLog(@"Warning! Your differential operator \"%@\" contains at least one non-finite value. This is likely not what you intended.", name);
	}
	
    [diffOp setName: name];
	
	// This should do an individual isEquals on each dimension.
	NSMutableDictionary *pool = [self.dimensionOperatorPoolMapTable objectForKey: diffOp.fromDimensions];
	
	if (!pool) {
		pool = [NSMutableDictionary dictionary];
		[self.dimensionOperatorPoolMapTable setObject: pool forKey: diffOp.fromDimensions];
	}
	
	[pool setObject: diffOp forKey: name];
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
