//
//  GLSimpleOperationOptimizer.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 4/18/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import "GLOperationOptimizer.h"
#import "GLMemoryPool.h"

#import "GLBasisTransformationOperations.h"

#define CONSERVATIVE_MEMORY_BUFFERS NO

@interface GLOperationOptimizer ()
@property(strong) NSArray *topVariables;
@property(strong) NSArray *bottomVariables;
@property BOOL alreadyInitialized;
@property(strong) dispatch_queue_t queue;
@end

@implementation GLOperationOptimizer
{
	NSArray *_topVariables;
	NSArray *_bottomVariables;
	NSMutableArray *_internalDataBuffers;
	NSMutableArray *_internalVariablesAndBuffers;
	NSMutableArray *_allVariablesAndBuffers;
    NSHashTable *_hasDataBuffer;
}

@synthesize topVariables = _topVariables;
@synthesize bottomVariables = _bottomVariables;
@synthesize internalDataBuffers=_internalDataBuffers;
@synthesize internalVariablesAndBuffers=_internalVariablesAndBuffers;
@synthesize allVariablesAndBuffers=_allVariablesAndBuffers;
@synthesize alreadyInitialized;

/************************************************/
/*		Convenience Methods						*/
/************************************************/

#pragma mark -
#pragma mark Convenience Methods
#pragma mark

- (void) createOptimizedOperation
{
    if (! [self createExecutionPlan]) {
		return;
	}
	if (! [self assignMemoryBuffers]) {
		return;
	}
	if (! [self createExecutionBlocks]) {
		return;
	}
	
	//NSLog(@"Total memory buffers allocated: %lu", self.totalMemoryBuffersAllocated);
	
//    self.operationSerialBlockCountMap = nil;
//    self.operationParallelBlockCountMap = nil;
//    self.operationParallelGroupCountMap = nil;
    
	if ( self.bottomVariables.count && self.topVariables.count )
	{
		NSMutableArray *groups = [self groupResponsibilitiesForOperation: (GLVariableOperation *) self];
		NSMapTable *parallelExecutionBlocks = [self parallelResponsibilitiesForOperation: (GLVariableOperation *) self];
		if (groups.count != 0) {
			NSLog(@"The top variable should never have groups to manage. Something went wrong.");
		}
		if (parallelExecutionBlocks.count != 0 ) {
			NSLog(@"The top variable should never have groups to manage. Something went wrong.");
		}
		
		// All five of these objects will be copied/referenced by the block.
		dispatch_queue_t globalQueue = self.queue;
		NSMutableArray *serialExecutionBlocks = [self serialResponsibilitiesForOperation: (GLVariableOperation *) self];
		dispatch_group_t bottomGroup = self.bottomVariableGroup;
		NSMutableArray *allGroups = [self allGroups];
		NSArray *precomputedBuffers = self.precomputedVariablesData;
        
		variableOperation aBlock = ^(NSArray *bottomBuffers, NSArray *topBuffers, NSArray *internalBuffers) {
			
			// A unary operation only has one bottom variable, so we only need to entire it once.
			// The group variable is in the group array, so it will automatically get incremented once.
			// 1. Increment the group count once.
			for (dispatch_group_t group in allGroups) {
				dispatch_group_enter(group);
			}
            
			// Follow our strict buffer order.
			NSMutableArray *dataBuffers = [[NSMutableArray alloc] initWithArray: bottomBuffers];
			[dataBuffers addObjectsFromArray: topBuffers];
			[dataBuffers addObjectsFromArray: precomputedBuffers];
			[dataBuffers addObjectsFromArray: internalBuffers];
            
			for ( executionBlock anExecutionBlock in serialExecutionBlocks ) {
				dispatch_async( globalQueue, ^{
					anExecutionBlock( dataBuffers );
				});
			}
			
			dispatch_group_wait(bottomGroup, DISPATCH_TIME_FOREVER);
		};
		
        // 1. set the operation block
		self.operationBlock = aBlock;
        
        // 2. the internal data buffer array should already be created during the -createExecutionBlocks method.
        
        // 3. create the operand variable prototypes.
        self.operandVariablePrototypes = [NSMutableArray array];
        for (GLVariable *variable in self.topVariables) {
            [self.operandVariablePrototypes addObject: [GLVariable variableWithPrototype: variable]];
        }
        
        // 4. create the result variable prototypes.
        self.resultVariablePrototypes = [NSMutableArray array];
        for (GLVariable *variable in self.bottomVariables) {
            [self.resultVariablePrototypes addObject: [GLVariable variableWithPrototype: variable]];
        }
        
	}
	else
	{
		NSLog(@"Incorrect number of top variables (%lu) or bottom variables (%lu) to create a unaryOperation", self.topVariables.count, self.bottomVariables.count);
	}
}

- (GLOperationOptimizer *) initWithTopVariables: (NSArray *) topVariables bottomVariables: (NSArray *) bottomVariables
{
	if ((self=[super init])) {
		
		self.topVariables = topVariables;
		self.bottomVariables = bottomVariables;
		self.precomputedVariableDataMap = [NSMapTable mapTableWithKeyOptions: NSPointerFunctionsObjectPointerPersonality valueOptions:NSPointerFunctionsObjectPointerPersonality];
		self.internalVariableBufferMap = [NSMapTable mapTableWithKeyOptions: NSPointerFunctionsObjectPointerPersonality valueOptions:NSPointerFunctionsObjectPointerPersonality];
        self.internalBufferBufferMap = [NSMapTable mapTableWithKeyOptions: NSPointerFunctionsObjectPointerPersonality valueOptions:NSPointerFunctionsObjectPointerPersonality];
		self.internalBufferArray = [NSMutableArray array];
		
		self.operationSerialBlockCountMap = [NSMapTable mapTableWithKeyOptions: NSPointerFunctionsObjectPointerPersonality valueOptions:NSPointerFunctionsObjectPointerPersonality];
		self.operationParallelBlockCountMap = [NSMapTable mapTableWithKeyOptions: NSPointerFunctionsObjectPointerPersonality valueOptions:NSPointerFunctionsObjectPointerPersonality];
		self.operationParallelGroupCountMap = [NSMapTable mapTableWithKeyOptions: NSPointerFunctionsObjectPointerPersonality valueOptions:NSPointerFunctionsObjectPointerPersonality];
		self.finishedMappingOperations = [NSHashTable hashTableWithOptions: NSHashTableObjectPointerPersonality];
		
		self.hasDataBuffer = [NSHashTable hashTableWithOptions: NSHashTableObjectPointerPersonality];
		
		self.finishedOperations = [NSHashTable hashTableWithOptions: NSHashTableObjectPointerPersonality];
		self.operationSerialDependencyBlockArrayMap = [NSMapTable mapTableWithKeyOptions: NSPointerFunctionsObjectPointerPersonality valueOptions:NSPointerFunctionsObjectPointerPersonality];
		self.operationParallelDependencyBlockArrayMap = [NSMapTable mapTableWithKeyOptions: NSPointerFunctionsObjectPointerPersonality valueOptions:NSPointerFunctionsObjectPointerPersonality];
		self.operationGroupArrayMap = [NSMapTable mapTableWithKeyOptions: NSPointerFunctionsObjectPointerPersonality valueOptions:NSPointerFunctionsObjectPointerPersonality];
		
		self.queue = dispatch_queue_create("com.earlyinnovations.OperationOptimizer", DISPATCH_QUEUE_CONCURRENT);
		
        [self createOptimizedOperation];
	}
	return self;
}

/************************************************/
/*		Create Execution Plan					*/
/************************************************/

#pragma mark -
#pragma mark Create Execution Plan
#pragma mark

- (BOOL) createExecutionPlan
{
	if (self.alreadyInitialized) {
		return NO;
	}
	self.alreadyInitialized = YES;
	
	NSMutableSet *bottomOperations = [[NSMutableSet alloc] init];
	for (GLFunction *bottomVariable in self.bottomVariables) {
		if (bottomVariable.lastOperation) [bottomOperations addObject: bottomVariable.lastOperation];
	}
	for (GLVariableOperation *operation in bottomOperations) {
		[self incrementParallelGroupCountForOperation: operation];
	}
	
	
	BOOL success = NO;
    while ( !success ) {
        self.finishedMappingOperations = [NSHashTable hashTableWithOptions: NSHashTableObjectPointerPersonality];
        success = YES;
        for (GLVariableOperation *operation in bottomOperations) {
            success &= [self precomputeVariablesWithOperation: operation forTopVariables: self.topVariables];
        }
    }
    
    self.finishedMappingOperations = [NSHashTable hashTableWithOptions: NSHashTableObjectPointerPersonality];
	for (GLVariableOperation *operation in bottomOperations) {
		success &= [self mapPreliminariesWithOperation: operation forTopVariables: self.topVariables];
	}
	
	return success;
}

- (BOOL) mapPreliminariesWithOperation: (GLVariableOperation *) operation forTopVariables: (NSArray *) topVariables
{
	if ( [self.finishedMappingOperations containsObject: operation] ) {
		return YES;
	} else {
		[self.finishedMappingOperations addObject: operation];
	}
	
	// We now determine the responsibilities of each parent.
	// See the documentation for an explanation of the algorithm.
	BOOL success = NO;
	
	NSMutableArray *topVariableOperands = [[NSMutableArray alloc] init];
	NSMutableArray *precomputedVariableOperands = [[NSMutableArray alloc] init];
	NSMutableArray *otherVariableOperandOperations = [[NSMutableArray alloc] init];
	
	for (GLFunction *aVariable in operation.operand) {
		if ( [topVariables containsObject: aVariable] ) {
			[topVariableOperands addObject: aVariable];
		} else if ( aVariable.pendingOperations.count == 0 ) {
			[precomputedVariableOperands addObject: aVariable];
		} else {
			if ( ![otherVariableOperandOperations containsObject: aVariable.lastOperation] ) {
				[otherVariableOperandOperations addObject: aVariable.lastOperation];
			}
		}
	}
	
    if (!operation.operation) {
        [NSException raise:@"OperationWithoutImplementation" format:@"This operation does not have an implementation!"];
    }
    
    //	static NSUInteger basisTransformCount = 0;
    //	if ( [[operation class] isSubclassOfClass: [GLBasisTransformOperation class]]) {
    //		basisTransformCount++;
    //		NSLog(@"%lu -- %@", basisTransformCount, [operation description]);
    //	}
	
	if ( precomputedVariableOperands.count && !topVariableOperands.count && !otherVariableOperandOperations.count ) {
		NSLog(@"This operation depends only on precomputed variables. This should not happen.");
	} else if ( operation.operand.count == 0 ) {
		[self incrementSerialBlockCountForOperation: (GLVariableOperation *) self];
	} else if ( topVariableOperands.count && !otherVariableOperandOperations.count ) {
		[self incrementSerialBlockCountForOperation: (GLVariableOperation *) self];
	} else if ( otherVariableOperandOperations.count == 1 ) {
		[self incrementSerialBlockCountForOperation: [otherVariableOperandOperations objectAtIndex: 0]];
	} else if ( otherVariableOperandOperations.count > 1 ){
		[self incrementParallelBlockCountForOperation: [otherVariableOperandOperations objectAtIndex: 0]];
		for (NSUInteger i=1; i<otherVariableOperandOperations.count; i++) {
			[self incrementParallelGroupCountForOperation: [otherVariableOperandOperations objectAtIndex: i]];
		}
	} else {
		NSLog(@"Ack!!! This case should never happen if my logic is correct");
	}
	
	success = YES;
	for (GLVariableOperation *anOperation in otherVariableOperandOperations) {
		success &= [self mapPreliminariesWithOperation: anOperation forTopVariables: topVariables];
	}
	
	return success;
}

- (BOOL) precomputeVariablesWithOperation: (GLVariableOperation *) operation forTopVariables: (NSArray *) topVariables
{
	if ( [self.finishedMappingOperations containsObject: operation] ) {
		return YES;
	} else {
		[self.finishedMappingOperations addObject: operation];
	}
	
	// We now determine the responsibilities of each parent.
	// See the documentation for an explanation of the algorithm.
	NSMutableArray *topVariableOperands = [[NSMutableArray alloc] init];
	NSMutableArray *precomputedVariableOperands = [[NSMutableArray alloc] init];
	NSMutableArray *otherVariableOperandOperations = [[NSMutableArray alloc] init];
	
	for (GLFunction *aVariable in operation.operand) {
		if ( [topVariables containsObject: aVariable] ) {
			[topVariableOperands addObject: aVariable];
		} else if ( aVariable.pendingOperations.count == 0 ) {
			[precomputedVariableOperands addObject: aVariable];
		} else {
			if ( ![otherVariableOperandOperations containsObject: aVariable.lastOperation] ) {
				[otherVariableOperandOperations addObject: aVariable.lastOperation];
			}
		}
	}
	
    BOOL success = YES;
	if ( precomputedVariableOperands.count && !topVariableOperands.count && !otherVariableOperandOperations.count ) {
		//NSLog(@"This operation depends only on precomputed variables. It will be computed now.");
        GLEquation *equation = [topVariables.lastObject equation];
        [equation solveForOperation: operation waitUntilFinished: YES];
        success = NO;
	}
	
	for (GLVariableOperation *anOperation in otherVariableOperandOperations) {
		success &= [self precomputeVariablesWithOperation: anOperation forTopVariables: topVariables];
	}
	
	return success;
}

@synthesize internalVariableBufferMap;
@synthesize operationSerialBlockCountMap;
@synthesize operationParallelBlockCountMap;
@synthesize operationParallelGroupCountMap;
@synthesize finishedMappingOperations;

- (void) incrementParallelGroupCountForOperation: (GLVariableOperation *) operation
{	
	NSNumber *number = [self.operationParallelGroupCountMap objectForKey: operation];
	if (!number) {
		number = [NSNumber numberWithUnsignedInteger: 1];
	} else {
		number = [NSNumber numberWithUnsignedInteger: number.unsignedIntegerValue+1];
	}
	
	[self.operationParallelGroupCountMap setObject: number forKey: operation];
}

- (NSUInteger) parallelGroupCountForOperation: (GLVariableOperation *) operation
{
	NSNumber *number = [self.operationParallelGroupCountMap objectForKey: operation];
	if (!number) {
		return 0;
	} else {
		return number.unsignedIntegerValue;
	}
}

- (void) incrementSerialBlockCountForOperation: (GLVariableOperation *) operation
{
	NSNumber *number = [self.operationSerialBlockCountMap objectForKey: operation];
	if (!number) {
		number = [NSNumber numberWithUnsignedInteger: 1];
	} else {
		number = [NSNumber numberWithUnsignedInteger: number.unsignedIntegerValue+1];
	}
	
	[self.operationSerialBlockCountMap setObject: number forKey: operation];
}

- (NSUInteger) serialBlockCountForOperation: (GLVariableOperation *) operation
{
	NSNumber *number = [self.operationSerialBlockCountMap objectForKey: operation];
	if (!number) {
		return 0;
	} else {
		return number.unsignedIntegerValue;
	}
}

- (void) incrementParallelBlockCountForOperation: (GLVariableOperation *) operation
{
	NSNumber *number = [self.operationParallelBlockCountMap objectForKey: operation];
	if (!number) {
		number = [NSNumber numberWithUnsignedInteger: 1];
	} else {
		number = [NSNumber numberWithUnsignedInteger: number.unsignedIntegerValue+1];
	}
	
	[self.operationParallelBlockCountMap setObject: number forKey: operation];
}

- (NSUInteger) parallelBlockCountForOperation: (GLVariableOperation *) operation
{
	NSNumber *number = [self.operationParallelBlockCountMap objectForKey: operation];
	if (!number) {
		return 0;
	} else {
		return number.unsignedIntegerValue;
	}
}

- (NSUInteger) totalDependencies
{
	NSUInteger total = 0;
	for ( id key in self.operationParallelBlockCountMap ) {
		NSNumber *number = [self.operationParallelBlockCountMap objectForKey: key];
		total += number.unsignedIntegerValue;
	}
	for ( id key in self.operationSerialBlockCountMap ) {
		NSNumber *number = [self.operationSerialBlockCountMap objectForKey: key];
		total += number.unsignedIntegerValue;
	}
	
	return total;
}


/************************************************/
/*		Assign Memory Buffers					*/
/************************************************/

#pragma mark -
#pragma mark Assign Memory Buffers
#pragma mark

- (BOOL) assignMemoryBuffers
{
	BOOL success = YES;
	for (GLFunction *variable in self.bottomVariables) {
		success &= [self assignUnoptimizedMemoryBufferToVariable: variable forTopVariables: self.topVariables bottomVariables:self.bottomVariables];
	}
	return success;
}

//- (BOOL) assignMemoryBuffers
//{
//	GLMemoryOptimizer *memoryOptimizer = [[GLMemoryOptimizer alloc] initWithTopVariables: self.topVariables bottomVariables: self.bottomVariables];
//
//	for (GLVariable *variable in self.bottomVariables) {
//		[memoryOptimizer mapAllGraphCyclesStartingWithVariable: variable];
//	}
//
//	for (GLVariable *variable in self.bottomVariables) {
//		[memoryOptimizer mapChildrenStartingWithVariable: variable];
//	}
//
//	for (GLVariable *variable in self.bottomVariables) {
//		[memoryOptimizer mapAllReachableVariablesStartingWithVariable: variable];
//	}
//
//	for (GLVariable *variable in self.bottomVariables) {
//		[memoryOptimizer assignMemoryBufferToParentsAndVariable: variable];
//	}
//
//	self.variableDataMap = memoryOptimizer.variableDataMap;
//
//	return YES;
//}

@synthesize hasDataBuffer=_hasDataBuffer;
@synthesize totalMemoryBuffersAllocated;

- (BOOL) assignUnoptimizedMemoryBufferToVariable: (GLFunction *) variable forTopVariables: (NSArray *) topVariables bottomVariables: (NSArray *) bottomVariables
{
	if ( [self.hasDataBuffer containsObject: variable] )
	{	// Already assigned a buffer to this variable, and therefore its parents as well.
		return YES;
	}
	else if ([topVariables containsObject: variable])
	{	// Top variables don't need buffers and we can exit the recusion.
		return YES;
	}
	else if (variable.pendingOperations.count == 0)
	{
		// Precomputed variables are a dead end. It's values are fixed and already populated.
		// So we simply need to grab a copy of its data and exit the recusion.
        [self.precomputedVariableDataMap setObject: [variable.data copy] forKey: variable];
        [self.hasDataBuffer addObject: variable];
		return YES;
	}
	else if (variable.pendingOperations.count == 1)
	{
		if (![bottomVariables containsObject: variable])
		{	// It's not a top variable or a bottom variable and thus it needs a buffer.
			
			self.totalMemoryBuffersAllocated = self.totalMemoryBuffersAllocated + 1;
			
			// Here we fetch a data object, of appropriate size, for the variable.
			[self.internalVariableBufferMap setObject: [[GLBuffer alloc] initWithLength: variable.dataBytes] forKey: variable];
			[self.hasDataBuffer addObject: variable];
		}
		
		// Now we work further up the tree and deal with the parents.
		GLVariableOperation *operation = variable.lastOperation;
        
        // But first, let's handle the buffers for this operation
        if ( ![self.hasDataBuffer containsObject: operation])
        {
            for (GLBuffer *buffer in operation.buffer) {
                [self.internalBufferBufferMap setObject: [[GLBuffer alloc] initWithLength: buffer.numBytes] forKey: buffer];
                self.totalMemoryBuffersAllocated = self.totalMemoryBuffersAllocated + 1;
            }
            [self.hasDataBuffer addObject: operation];
        }

		BOOL success = YES;
		for (GLFunction *aVariable in operation.operand) {
			success &= [self assignUnoptimizedMemoryBufferToVariable: aVariable forTopVariables: topVariables bottomVariables: bottomVariables];
		}
		
		// If one of the result variables of this operation does not get used, it will not appear in this graph. However, it still needs a data buffer as the operation needs to write somewhere.
		for (GLFunction *aVariable in operation.result) {
			if ( ![self.hasDataBuffer containsObject: aVariable] )
			{
				self.totalMemoryBuffersAllocated = self.totalMemoryBuffersAllocated + 1;
				
				// Here we fetch a data object, of appropriate size, for the variable.
				[self.internalVariableBufferMap setObject: [[GLBuffer alloc] initWithLength: aVariable.dataBytes] forKey: aVariable];
				[self.hasDataBuffer addObject: aVariable];
			}
		}
		
		return success;

		
		return NO;
	} else {
		NSLog(@"Trying to assign memory buffers, but pending operations for this variable are greater than 1. This violates our current assumptions.");
		return NO;
	}
}

/************************************************/
/*		Execution Plan Creation					*/
/************************************************/

#pragma mark -
#pragma mark Execution Plan Creation
#pragma mark

- (BOOL) createExecutionBlocks
{
	// Now we move the internal map to two arrays that are in 1-1 correspondence.
	self.internalVariablesAndBuffers = [NSMutableArray array];
	self.internalDataBuffers = [NSMutableArray array];
	for ( GLFunction *aVariable in self.internalVariableBufferMap ) {
		[self.internalVariablesAndBuffers addObject: aVariable];
		[self.internalDataBuffers addObject: [self.internalVariableBufferMap objectForKey: aVariable]];
	}
    for ( GLFunction *aVariable in self.internalBufferBufferMap ) {
		[self.internalVariablesAndBuffers addObject: aVariable];
		[self.internalDataBuffers addObject: [self.internalBufferBufferMap objectForKey: aVariable]];
	}
	
    // Also move this internal map to two arrays that are in 1-1 correspondence.
    self.precomputedVariables = [NSMutableArray array];
	self.precomputedVariablesData = [NSMutableArray array];
    for ( GLFunction *aVariable in self.precomputedVariableDataMap ) {
		[self.precomputedVariables addObject: aVariable];
		[self.precomputedVariablesData addObject: [self.precomputedVariableDataMap objectForKey: aVariable]];
	}
    
	// We build a carefully ordered array of all the variables to be referenced.
	self.allVariablesAndBuffers = [[NSMutableArray alloc] initWithArray: self.bottomVariables];
	[self.allVariablesAndBuffers addObjectsFromArray: self.topVariables];
    [self.allVariablesAndBuffers addObjectsFromArray: self.precomputedVariables];
	[self.allVariablesAndBuffers addObjectsFromArray: self.internalVariablesAndBuffers];
	
	// All bottom variables are responsible for exiting a group that can be monitored to determine when we're done executing.
	self.bottomVariableGroup = dispatch_group_create();
	NSMutableSet *bottomOperations = [[NSMutableSet alloc] init];
	for (GLFunction *bottomVariable in self.bottomVariables) {
		if (bottomVariable.lastOperation) [bottomOperations addObject: bottomVariable.lastOperation];
	}
	for (GLVariableOperation *operation in bottomOperations) {
		[self addGroupResponsibility: self.bottomVariableGroup forOperation: operation];
	}
	
	BOOL success = NO;
	NSUInteger rationalCutoff = [self totalDependencies]+1;
	NSUInteger totalLoops = 0;
	while (!success && totalLoops < rationalCutoff) {
		success = YES;
		for (GLVariableOperation *operation in bottomOperations) {
			success &= [self createExecutionBlockFromOperation: operation forTopVariables: self.topVariables bottomVariables: self.bottomVariables];
		}
		totalLoops++;
	}
	
	if (!success) {
		NSLog(@"Apparently the algorithm failed to converge.");
	}
	
	return success;
}

@synthesize finishedOperations;
@synthesize operationGroupArrayMap;
@synthesize operationSerialDependencyBlockArrayMap;
@synthesize operationParallelDependencyBlockArrayMap;
@synthesize bottomVariableGroup;
@synthesize totalTopVariablesCreated;

- (void) addGroupResponsibility: (dispatch_group_t) aGroup forOperation: (GLVariableOperation *) operation
{
	NSMutableArray *array = [self.operationGroupArrayMap objectForKey: operation];
	if (!array) {
		array = [NSMutableArray array];
		[self.operationGroupArrayMap setObject: array forKey: operation];
	}
	
	[array addObject: aGroup];
}

- (NSMutableArray *) groupResponsibilitiesForOperation: (GLVariableOperation *) operation
{
	NSMutableArray *array = [self.operationGroupArrayMap objectForKey: operation];
	if (!array) {
        array = [NSMutableArray array];
        [self.operationGroupArrayMap setObject: array forKey: operation];
	}
	
	return array;
}

- (NSMutableArray *) allGroups
{
	NSMutableArray *allGroups = [NSMutableArray array];
	for ( GLFunction *variable in self.operationGroupArrayMap ) {
		[allGroups addObjectsFromArray: [self.operationGroupArrayMap objectForKey: variable]];
	}
		
	return allGroups;
}

- (void) addSerialResponsibility: (executionBlock) aBlock forOperation: (GLVariableOperation *) operation
{
	NSMutableArray *array = [self.operationSerialDependencyBlockArrayMap objectForKey: operation];
	if (!array) {
		array = [NSMutableArray array];
		[self.operationSerialDependencyBlockArrayMap setObject: array forKey: operation];
	}
	
	[array addObject: aBlock];
}

- (NSMutableArray *) serialResponsibilitiesForOperation: (GLVariableOperation *) operation
{
	NSMutableArray *array = [self.operationSerialDependencyBlockArrayMap objectForKey: operation];
	if (!array) {
        array = [NSMutableArray array];
        [self.operationSerialDependencyBlockArrayMap setObject: array forKey: operation];
	}
	
	return array;
}

- (void) addParallelResponsibility: (executionBlock) aBlock withGroup: (dispatch_group_t) group forOperation: (GLVariableOperation *) operation;
{
	NSMapTable *mapTable = [self.operationParallelDependencyBlockArrayMap objectForKey: operation];
	if (!mapTable) {
		mapTable = [NSMapTable mapTableWithKeyOptions: NSPointerFunctionsObjectPointerPersonality valueOptions:NSPointerFunctionsObjectPointerPersonality];
		[self.operationParallelDependencyBlockArrayMap setObject: mapTable forKey: operation];
	}
	
	[mapTable setObject: group forKey: aBlock];
}

- (NSMapTable *) parallelResponsibilitiesForOperation: (GLVariableOperation *) operation;
{
	NSMapTable *mapTable = [self.operationParallelDependencyBlockArrayMap objectForKey: operation];
	if (!mapTable) {
		mapTable = [NSMapTable mapTableWithKeyOptions: NSPointerFunctionsObjectPointerPersonality valueOptions:NSPointerFunctionsObjectPointerPersonality];
        [self.operationParallelDependencyBlockArrayMap setObject: mapTable forKey: operation];
	}
	
	return mapTable;
}

- (void) sanityCheck
{
    NSMutableSet *executionBlockSet = [[NSMutableSet alloc] init];
    for (GLVariableOperation *operation in self.operationSerialDependencyBlockArrayMap) {
        NSArray *serialBlocks = [self.operationSerialDependencyBlockArrayMap objectForKey: operation];
        for (executionBlock aBlock in serialBlocks) {
            if ([executionBlockSet containsObject: aBlock]) {
                [NSException raise:@"NotGood" format:@"VeryBad"];
            } else {
                [executionBlockSet addObject: aBlock];
            }
        }
    }
    
    for (GLVariableOperation *operation  in self.operationParallelDependencyBlockArrayMap) {
        NSMapTable *mapTable = [self.operationParallelDependencyBlockArrayMap objectForKey: operation];
        for (executionBlock aBlock in mapTable) {
            if ([executionBlockSet containsObject: aBlock]) {
                [NSException raise:@"NotGood" format:@"VeryBad"];
            } else {
                [executionBlockSet addObject: aBlock];
            }
        }
    }
    
    NSLog(@"total dependencies: %lu, total blocks: %lu, finished operations: %lu", self.totalDependencies, executionBlockSet.count, self.finishedOperations.count);
    
    
}

- (BOOL) isTopOperationReady {
	GLVariableOperation * operation = (GLVariableOperation *) self;
	
	// If this test passes, then it means we're ready to construct the actual planBlock that needs to be run
	NSMutableArray *groups = [self groupResponsibilitiesForOperation: operation];
	NSMutableArray *serialExecutionBlocks = [self serialResponsibilitiesForOperation: operation];
	NSMapTable *parallelExecutionBlocks = [self parallelResponsibilitiesForOperation: operation];
	
	// The executionBlock for this operation will only be created when it has the expected number
	// of dependent children blocks and the expected number of groups.
	BOOL canCreate = ([self serialBlockCountForOperation: operation] == serialExecutionBlocks.count &&
					  [self parallelGroupCountForOperation: operation] == groups.count &&
					  [self parallelBlockCountForOperation: operation] == parallelExecutionBlocks.count);
	
	return canCreate;
}

// Returns yes its executionBlock and all the parents are successfully created, no otherwise.
- (BOOL) createExecutionBlockFromOperation: (GLVariableOperation *) operation forTopVariables: (NSArray *) topVariables bottomVariables: (NSArray *) bottomVariables
{
	if (self.isTopOperationReady) {
		return YES;
	}
	
	// Sort the operands into categories
	NSMutableArray *topVariableOperands = [[NSMutableArray alloc] init];
	NSMutableArray *precomputedVariableOperands = [[NSMutableArray alloc] init];
	NSMutableArray *otherVariableOperandOperations = [[NSMutableArray alloc] init];
	
	for (GLFunction *aVariable in operation.operand) {
		if ( [topVariables containsObject: aVariable] ) {
			[topVariableOperands addObject: aVariable];
		} else if ( aVariable.pendingOperations.count == 0 ) {
			[precomputedVariableOperands addObject: aVariable];
		} else {
			if ( ![otherVariableOperandOperations containsObject: aVariable.lastOperation] ) {
				[otherVariableOperandOperations addObject: aVariable.lastOperation];
			}
		}
	}
	
	// Any variable already in the variables set is already part of the execution plan.
	BOOL executionBlockAlreadyCreated = [self.finishedOperations containsObject: operation];
	
	// If it's not already created, we need to check and see if we can create it.
	if (!executionBlockAlreadyCreated)
	{
		// If this test passes, then it means we're ready to construct the actual planBlock that needs to be run
		NSMutableArray *groups = [self groupResponsibilitiesForOperation: operation];
		NSMutableArray *serialExecutionBlocks = [self serialResponsibilitiesForOperation: operation];
		NSMapTable *parallelExecutionBlocks = [self parallelResponsibilitiesForOperation: operation];
		
		// The executionBlock for this operation will only be created when it has the expected number
		// of dependent children blocks and the expected number of groups.
		BOOL canCreate = ([self serialBlockCountForOperation: operation] == serialExecutionBlocks.count &&
						  [self parallelGroupCountForOperation: operation] == groups.count &&
						  [self parallelBlockCountForOperation: operation] == parallelExecutionBlocks.count);
		
		if (canCreate)
		{
			[self.finishedOperations addObject: operation];
			dispatch_queue_t globalQueue = self.queue;
            
			// The childrenBlock is responsible for sending off the children to perform their duties.
			executionBlock childrenBlock = ^( NSArray *dataBuffers ) {
				// a. serial executions first
				for ( executionBlock anExecutionBlock in serialExecutionBlocks ) {
					dispatch_async( globalQueue, ^{
						anExecutionBlock( dataBuffers );
					});
				}
				
				// b. parallel executions second
                for ( executionBlock anExecutionBlock in parallelExecutionBlocks ) {
                    dispatch_group_t group = [parallelExecutionBlocks objectForKey: anExecutionBlock];
                    dispatch_group_notify( group, globalQueue, ^{
                        anExecutionBlock( dataBuffers );
                    });
                }
                
                // c. Allow any thread joins/pending binary operations to get going.
                for (dispatch_group_t group in groups) {
                    dispatch_group_leave(group);
                }
			};
			
			// Figure out the indices
			NSMutableArray *resultIndices = [[NSMutableArray alloc] init];
			NSMutableArray *operandIndices = [[NSMutableArray alloc] init];
			NSMutableArray *bufferIndices = [[NSMutableArray alloc] init];
			for (GLFunction *aVariable in operation.result) {
				NSUInteger anIndex = [self.allVariablesAndBuffers indexOfObject: aVariable];
				if (anIndex == NSNotFound) {
					[NSException raise: @"Invalid index." format: @"The operation is malformed."];
				} else {
					[resultIndices addObject: @(anIndex)];
				}
			}
			for (GLFunction *aVariable in operation.operand) {
				NSUInteger anIndex = [self.allVariablesAndBuffers indexOfObject: aVariable];
				if (anIndex == NSNotFound) {
					[NSException raise: @"Invalid index." format: @"The operation is malformed."];
				} else {
					[operandIndices addObject: @(anIndex)];
				}
			}
			for (GLBuffer *aBuffer in operation.buffer) {
				NSUInteger anIndex = [self.allVariablesAndBuffers indexOfObject: aBuffer];
				if (anIndex == NSNotFound) {
					[NSException raise: @"Invalid index." format: @"The operation is malformed."];
				} else {
					[bufferIndices addObject: @(anIndex)];
				}
			}
			
            NSMutableArray *result = [[NSMutableArray alloc] init];
            NSMutableArray *operand = [[NSMutableArray alloc] init];
            NSMutableArray *buffer = [[NSMutableArray alloc] init];
			
			variableOperation operationBlock = operation.operation;            
			executionBlock theExecutionBlock = ^( NSArray *dataBuffers ) {
                // Initializing these arrays and adding objects at each time is VERY inefficient. It's the largest source of delay in the code currently.
//                NSMutableArray *result = [[NSMutableArray alloc] init];
//                NSMutableArray *operand = [[NSMutableArray alloc] init];
//                NSMutableArray *buffer = [[NSMutableArray alloc] init];
                
				// Note that we do not call -objectsAtIndexes because order is important.
				for (NSNumber *anIndex in resultIndices) {
                    NSUInteger intValue = anIndex.unsignedIntegerValue;
                    NSData * dataObject = [dataBuffers objectAtIndex: intValue];
					[result addObject: dataObject];
				}
				
				for (NSNumber *anIndex in operandIndices) {
                    NSUInteger intValue = anIndex.unsignedIntegerValue;
					[operand addObject: [dataBuffers objectAtIndex: intValue]];
				}
				
				for (NSNumber *anIndex in bufferIndices) {
                    NSUInteger intValue = anIndex.unsignedIntegerValue;
					[buffer addObject: [dataBuffers objectAtIndex: intValue]];
				}
				
				// 1. Compute our own operation
				operationBlock( result, operand, buffer );
				
                [result removeAllObjects];
                [operand removeAllObjects];
                [buffer removeAllObjects];
                
				// 2. Send off the children
				childrenBlock( dataBuffers );
			};
			
			if ( precomputedVariableOperands.count && !topVariableOperands.count && !otherVariableOperandOperations.count ) {
				NSLog(@"This operation depends only on precomputed variables. We have not yet implemented the appropriate optimization to deal with this.");
			} else if (operation.operand.count == 0) {
				[self addSerialResponsibility: theExecutionBlock forOperation: (GLVariableOperation *) self];
			} else if ( topVariableOperands.count && !otherVariableOperandOperations.count ) {
				[self addSerialResponsibility: theExecutionBlock forOperation: (GLVariableOperation *) self];
//				NSLog(@"<0x%x> %@ serial from top variable",(unsigned int) operation, NSStringFromClass([operation class]) );
			} else if ( otherVariableOperandOperations.count == 1 ) {
				[self addSerialResponsibility: theExecutionBlock forOperation: [otherVariableOperandOperations objectAtIndex: 0]];
//				NSLog(@"<0x%x> %@ serial from <0x%x> %@",(unsigned int) operation, NSStringFromClass([operation class]), (unsigned int) [otherVariableOperandOperations objectAtIndex: 0], NSStringFromClass( [[otherVariableOperandOperations objectAtIndex: 0] class]) );
			} else if ( otherVariableOperandOperations.count > 1 ){
				dispatch_group_t group = dispatch_group_create();
				[self addParallelResponsibility: theExecutionBlock withGroup: group forOperation: [otherVariableOperandOperations objectAtIndex: 0]];
//				NSLog(@"<0x%x> %@ async from <0x%x> %@",(unsigned int) operation, NSStringFromClass([operation class]), (unsigned int) [otherVariableOperandOperations objectAtIndex: 0], NSStringFromClass( [[otherVariableOperandOperations objectAtIndex: 0] class]) );
				for (NSUInteger i=1; i<otherVariableOperandOperations.count; i++) {
					[self addGroupResponsibility: group forOperation: [otherVariableOperandOperations objectAtIndex: i]];
//					NSLog(@"<0x%x> %@ group left from <0x%x> %@",(unsigned int) operation, NSStringFromClass([operation class]), (unsigned int) [otherVariableOperandOperations objectAtIndex: i], NSStringFromClass( [[otherVariableOperandOperations objectAtIndex: i] class]) );
				}
			} else {
				NSLog(@"Aaaaack!!!");
			}	
		} else {
			// if (canCreate)
			//NSLog(@"<0x%x> (expected, found)--Serial: (%lu, %lu), Parallel: (%lu, %lu), Parallel groups: (%lu, %lu), %@",(unsigned int) operation, [self serialBlockCountForOperation: operation], serialExecutionBlocks.count, [self parallelBlockCountForOperation: operation], parallelExecutionBlocks.count, [self parallelGroupCountForOperation: operation], groups.count, NSStringFromClass([operation class]) );
			return NO;
		}
	}
	
	// If we've gotten this far, the execution block for this variable was either already created, or we just created it.
	// So now we need to check how the parents are doing!	
	BOOL success = YES;
	for (GLVariableOperation *anOperation in otherVariableOperandOperations) {
		success &= [self createExecutionBlockFromOperation: anOperation forTopVariables: topVariables bottomVariables: bottomVariables];
	}
	
	return success;
}


@end
