//
//  GLSimpleOperationOptimizer.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 4/18/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import "GLOperationOptimizer.h"
#import "GLMemoryPool.h"
#import "GLMemoryOptimizer.h"

#import "GLBasisTransformationOperations.h"

#define CONSERVATIVE_MEMORY_BUFFERS NO

@interface GLOperationOptimizer ()
@property(strong) NSArray *topVariables;
@property(strong) NSArray *bottomVariables;
@property BOOL alreadyInitialized;
@end

@implementation GLOperationOptimizer
{
	NSArray *_topVariables;
	NSArray *_bottomVariables;
	NSMutableArray *_internalDataBuffers;
	NSMutableArray *_internalVariables;
	NSMutableArray *_allVariables;
    NSHashTable *_hasDataBuffer;
    dispatch_queue_t _childrenQueue;
}

@synthesize topVariables = _topVariables;
@synthesize bottomVariables = _bottomVariables;
@synthesize internalDataBuffers=_internalDataBuffers;
@synthesize internalVariables=_internalVariables;
@synthesize allVariables=_allVariables;
@synthesize alreadyInitialized;
@synthesize childrenQueue=_childrenQueue;

/************************************************/
/*		Convenience Methods						*/
/************************************************/

#pragma mark -
#pragma mark Convenience Methods
#pragma mark


+ (void) addVariable: (GLVariable *) variable toTopVariables: (NSMutableSet *) set
{
	GLVariableOperation *operation = variable.lastOperation;
	
	if (!operation) {
		return;
	} else if ( operation.operationType == kGLUnaryOperation) {
		GLUnaryOperation *aUnaryOperation = (GLUnaryOperation *) operation;
		if (aUnaryOperation.operand.pendingOperations.count) {
			[self addVariable: aUnaryOperation.operand toTopVariables: set];
		} else {
			[set addObject: variable];
		}
	} else if ( operation.operationType == kGLBinaryOperation ) {
		GLBinaryOperation *aBinaryOperation = (GLBinaryOperation *) operation;
		if (aBinaryOperation.firstOperand.pendingOperations.count) {
			[self addVariable: aBinaryOperation.firstOperand toTopVariables: set];
		} else if (aBinaryOperation.secondOperand.pendingOperations.count) {
			[self addVariable: aBinaryOperation.secondOperand toTopVariables: set];
		} else {
			[set addObject: variable];
		}
	}
}

+ (NSSet *) topVariablesFromVariable: (GLVariable *) variable
{
	NSMutableSet *set = [[NSMutableSet alloc] init];
	
	[self addVariable: variable toTopVariables: set];
	
	return set;
}

- (GLOperationOptimizer *) initWithTopVariables: (NSArray *) topVariables bottomVariables: (NSArray *) bottomVariables
{
	if ((self=[super init])) {
		self.variableDataMap = [NSMapTable mapTableWithKeyOptions: NSMapTableObjectPointerPersonality valueOptions:NSMapTableObjectPointerPersonality];
		self.operationSerialBlockCountMap = [NSMapTable mapTableWithKeyOptions: NSMapTableObjectPointerPersonality valueOptions:NSMapTableObjectPointerPersonality];
		self.operationParallelBlockCountMap = [NSMapTable mapTableWithKeyOptions: NSMapTableObjectPointerPersonality valueOptions:NSMapTableObjectPointerPersonality];
		self.operationParallelGroupCountMap = [NSMapTable mapTableWithKeyOptions: NSMapTableObjectPointerPersonality valueOptions:NSMapTableObjectPointerPersonality];
		self.finishedMappingOperations = [NSHashTable hashTableWithOptions: NSHashTableObjectPointerPersonality];
		
		self.hasDataBuffer = [NSHashTable hashTableWithOptions: NSHashTableObjectPointerPersonality];
		
		self.finishedOperations = [NSHashTable hashTableWithOptions: NSHashTableObjectPointerPersonality];
		self.operationSerialDependencyBlockArrayMap = [NSMapTable mapTableWithKeyOptions: NSMapTableObjectPointerPersonality valueOptions:NSMapTableObjectPointerPersonality];
		self.operationParallelDependencyBlockArrayMap = [NSMapTable mapTableWithKeyOptions: NSMapTableObjectPointerPersonality valueOptions:NSMapTableObjectPointerPersonality];
		self.operationGroupArrayMap = [NSMapTable mapTableWithKeyOptions: NSMapTableObjectPointerPersonality valueOptions:NSMapTableObjectPointerPersonality];
		
        self.childrenQueue = dispatch_queue_create("com.earlyinnovations.operationOptimizer.childrenQueue", 0);
        
		self.topVariables = topVariables;
		self.bottomVariables = bottomVariables;		
	}
	return self;
}

- (BOOL) createExecutionPlan
{
	if (self.alreadyInitialized) {
		return NO;
	}
	self.alreadyInitialized = YES;
	
	NSMutableSet *bottomOperations = [[NSMutableSet alloc] init];
	for (GLVariable *bottomVariable in self.bottomVariables) {
#warning This if statement might screw things up.
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

- (BOOL) assignMemoryBuffers
{	
	BOOL success = YES;
	for (GLVariable *variable in self.bottomVariables) {
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

- (BOOL) createExecutionBlocks
{
	// Now we move the internal map to two arrays that are in 1-1 correspondence.
	self.internalVariables = [[NSMutableArray alloc] initWithCapacity: self.variableDataMap.count];
	self.internalDataBuffers = [[NSMutableArray alloc] initWithCapacity: self.variableDataMap.count];
	for ( GLVariable *aVariable in self.variableDataMap ) {
		[self.internalVariables addObject: aVariable];
		[self.internalDataBuffers addObject: [self.variableDataMap objectForKey: aVariable]];
	}
	
	// We build a carefully ordered array of all the variables to be referenced.
	self.allVariables = [[NSMutableArray alloc] initWithArray: self.bottomVariables];
	[self.allVariables addObjectsFromArray: self.topVariables];
	[self.allVariables addObjectsFromArray: self.internalVariables];
	
	// All bottom variables are responsible for exiting a group that can be monitored to determine when we're done executing.
	self.bottomVariableGroup = dispatch_group_create();
	NSMutableSet *bottomOperations = [[NSMutableSet alloc] init];
	for (GLVariable *bottomVariable in self.bottomVariables) {
#warning This if statement might screw things up.
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

- (unaryOperation) unaryOperationBlock
{
	if (! [self createExecutionPlan]) {
		return nil;
	}
	if (! [self assignMemoryBuffers]) {
		return nil;
	}
	if (! [self createExecutionBlocks]) {
		return nil;
	}
	
	if ( self.bottomVariables.count == 1 && self.topVariables.count == 1)
	{		
		NSHashTable *groups = [self groupResponsibilitiesForOperation: (GLVariableOperation *) self];
		NSMapTable *parallelExecutionBlocks = [self parallelResponsibilitiesForOperation: (GLVariableOperation *) self];
		if (groups.count != 0) {
			NSLog(@"The top variable should never have groups to manage. Something went wrong.");
		}
		if (parallelExecutionBlocks.count != 0 ) {
			NSLog(@"The top variable should never have groups to manage. Something went wrong.");
		}
		
		// All five of these objects will be copied/referenced by the block.
		dispatch_queue_t globalQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
		NSHashTable *serialExecutionBlocks = [self serialResponsibilitiesForOperation: (GLVariableOperation *) self];
		dispatch_group_t bottomGroup = self.bottomVariableGroup;
		NSHashTable *allGroups = [self allGroups];
		NSArray *internalDataBuffers = self.internalDataBuffers;
		
		unaryOperation aBlock = ^(NSMutableData *bottomResult, NSData *topOperand) {			
			// A unary operation only has one bottom variable, so we only need to entire it once.
			// The group variable is in the group array, so it will automatically get incremented once.
			// 1. Increment the group count once.
			for (dispatch_group_t group in allGroups) {
				dispatch_group_enter(group);
			}
			
			NSMutableArray *dataBuffers = [[NSMutableArray alloc] initWithObjects: bottomResult, topOperand, nil];
			[dataBuffers addObjectsFromArray: internalDataBuffers];
			
			for ( executionBlock anExecutionBlock in serialExecutionBlocks ) {
				dispatch_async( globalQueue, ^{
					anExecutionBlock( dataBuffers );
				});
			}
			
			dispatch_group_wait(bottomGroup, DISPATCH_TIME_FOREVER);
		};
		
		return aBlock;
	}
	else
	{
		NSLog(@"Incorrect number of top variables (%lu) or bottom variables (%lu) to create a unaryOperation", self.topVariables.count, self.bottomVariables.count);
		return nil;
	}
}

- (binaryOperation) binaryOperationBlock
{
	if (! [self createExecutionPlan]) {
		return nil;
	}
	if (! [self assignMemoryBuffers]) {
		return nil;
	}
	if (! [self createExecutionBlocks]) {
		return nil;
	}
	
	if ( self.bottomVariables.count == 1 && self.topVariables.count == 2)
	{		
		NSHashTable *groups = [self groupResponsibilitiesForOperation: (GLVariableOperation *) self];
		NSMapTable *parallelExecutionBlocks = [self parallelResponsibilitiesForOperation: (GLVariableOperation *) self];
		if (groups.count != 0) {
			NSLog(@"The top variable should never have groups to manage. Something went wrong.");
		}
		if (parallelExecutionBlocks.count != 0 ) {
			NSLog(@"The top variable should never have groups to manage. Something went wrong.");
		}
		
		// All five of these objects will be copied/referenced by the block.
		dispatch_queue_t globalQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
		NSHashTable *serialExecutionBlocks = [self serialResponsibilitiesForOperation: (GLVariableOperation *) self];
		dispatch_group_t bottomGroup = self.bottomVariableGroup;
		NSHashTable *allGroups = [self allGroups];
		NSArray *internalDataBuffers = self.internalDataBuffers;
		
		binaryOperation aBlock = ^(NSMutableData *bottomResult, NSData *firstOperand, NSData *secondOperand) {			
			// A unary operation only has one bottom variable, so we only need to entire it once.
			// The group variable is in the group array, so it will automatically get incremented once.
			// 1. Increment the group count once.
			for (dispatch_group_t group in allGroups) {
				dispatch_group_enter(group);
			}
			
			NSMutableArray *dataBuffers = [[NSMutableArray alloc] initWithObjects: bottomResult, firstOperand, secondOperand, nil];
			[dataBuffers addObjectsFromArray: internalDataBuffers];
			
			for ( executionBlock anExecutionBlock in serialExecutionBlocks ) {
				dispatch_async( globalQueue, ^{
					anExecutionBlock( dataBuffers );
				});
			}
			
			dispatch_group_wait(bottomGroup, DISPATCH_TIME_FOREVER);
		};
		
		return aBlock;
	}
	else
	{
		NSLog(@"Incorrect number of top variables (%lu) or bottom variables (%lu) to create a unaryOperation", self.topVariables.count, self.bottomVariables.count);
		return nil;
	}
}

- (unaryVectorOperation) unaryVectorOperationBlock
{
	if (! [self createExecutionPlan]) {
		return nil;
	}
	if (! [self assignMemoryBuffers]) {
		return nil;
	}
	if (! [self createExecutionBlocks]) {
		return nil;
	}
	
	NSLog(@"Total memory buffers allocated: %lu", self.totalMemoryBuffersAllocated);
	
	if ( self.bottomVariables.count && self.topVariables.count )
	{		
		NSHashTable *groups = [self groupResponsibilitiesForOperation: (GLVariableOperation *) self];
		NSMapTable *parallelExecutionBlocks = [self parallelResponsibilitiesForOperation: (GLVariableOperation *) self];
		if (groups.count != 0) {
			NSLog(@"The top variable should never have groups to manage. Something went wrong.");
		}
		if (parallelExecutionBlocks.count != 0 ) {
			NSLog(@"The top variable should never have groups to manage. Something went wrong.");
		}
		
		// All five of these objects will be copied/referenced by the block.
		dispatch_queue_t globalQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
		NSHashTable *serialExecutionBlocks = [self serialResponsibilitiesForOperation: (GLVariableOperation *) self];
		dispatch_group_t bottomGroup = self.bottomVariableGroup;
		NSHashTable *allGroups = [self allGroups];
		NSArray *internalDataBuffers = self.internalDataBuffers;
				
		unaryVectorOperation aBlock = ^(NSArray *bottomBuffers, NSArray *topBuffers) {
			
			// A unary operation only has one bottom variable, so we only need to entire it once.
			// The group variable is in the group array, so it will automatically get incremented once.
			// 1. Increment the group count once.
			for (dispatch_group_t group in allGroups) {
				dispatch_group_enter(group);
			}

			
			NSMutableArray *dataBuffers = [[NSMutableArray alloc] initWithArray: bottomBuffers];
			[dataBuffers addObjectsFromArray: topBuffers];
			[dataBuffers addObjectsFromArray: internalDataBuffers];
            
			for ( executionBlock anExecutionBlock in serialExecutionBlocks ) {
				dispatch_async( globalQueue, ^{
					anExecutionBlock( dataBuffers );
				});
			}
			
			dispatch_group_wait(bottomGroup, DISPATCH_TIME_FOREVER);
		};
		
		return aBlock;
	}
	else
	{
		NSLog(@"Incorrect number of top variables (%lu) or bottom variables (%lu) to create a unaryOperation", self.topVariables.count, self.bottomVariables.count);
		return nil;
	}
}

/************************************************/
/*		Preliminary Mapping						*/
/************************************************/

#pragma mark -
#pragma mark Preliminary Mapping
#pragma mark

@synthesize variableDataMap;
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
	NSMutableArray *operands = [[NSMutableArray alloc] init];
	
	if (operation.operationType == kGLUnaryOperation) {
		GLUnaryOperation *aUnaryOperation = (GLUnaryOperation *) operation;
		[operands addObject: aUnaryOperation.operand];
	} else if (operation.operationType == kGLBinaryOperation) {
		GLBinaryOperation *aBinaryOperation = (GLBinaryOperation *) operation;
		[operands addObject: aBinaryOperation.firstOperand];
		[operands addObject: aBinaryOperation.secondOperand];
	} else if (operation.operationType == kGLUnaryVectorOperation) {
		GLUnaryVectorOperation *aUnaryVectorOperation = (GLUnaryVectorOperation *) operation;
		[operands addObjectsFromArray: aUnaryVectorOperation.operand];
	} else if (operation.operationType == kGLBinaryVectorOperation) {
		GLBinaryVectorOperation *aBinarVectorOperation = (GLBinaryVectorOperation *) operation;
		[operands addObjectsFromArray: aBinarVectorOperation.firstOperand];
		[operands addObjectsFromArray: aBinarVectorOperation.secondOperand];
	}
	
	NSMutableArray *topVariableOperands = [[NSMutableArray alloc] init];
	NSMutableArray *precomputedVariableOperands = [[NSMutableArray alloc] init];
	NSMutableArray *otherVariableOperandOperations = [[NSMutableArray alloc] init];
	
	for (GLVariable *aVariable in operands) {
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
	
//	static NSUInteger basisTransformCount = 0;
//	if ( [[operation class] isSubclassOfClass: [GLBasisTransformOperation class]]) {
//		basisTransformCount++;
//		NSLog(@"%lu -- %@", basisTransformCount, [operation description]);
//	}
	
	if ( precomputedVariableOperands.count && !topVariableOperands.count && !otherVariableOperandOperations.count ) {
		NSLog(@"This operation depends only on precomputed variables. We have not yet implemented the appropriate optimization to deal with this.");
	} else if ( operation.operationType == kGLNullaryOperation ) {
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
	NSMutableArray *operands = [[NSMutableArray alloc] init];
	
	if (operation.operationType == kGLUnaryOperation) {
		GLUnaryOperation *aUnaryOperation = (GLUnaryOperation *) operation;
		[operands addObject: aUnaryOperation.operand];
	} else if (operation.operationType == kGLBinaryOperation) {
		GLBinaryOperation *aBinaryOperation = (GLBinaryOperation *) operation;
		[operands addObject: aBinaryOperation.firstOperand];
		[operands addObject: aBinaryOperation.secondOperand];
	} else if (operation.operationType == kGLUnaryVectorOperation) {
		GLUnaryVectorOperation *aUnaryVectorOperation = (GLUnaryVectorOperation *) operation;
		[operands addObjectsFromArray: aUnaryVectorOperation.operand];
	} else if (operation.operationType == kGLBinaryVectorOperation) {
		GLBinaryVectorOperation *aBinarVectorOperation = (GLBinaryVectorOperation *) operation;
		[operands addObjectsFromArray: aBinarVectorOperation.firstOperand];
		[operands addObjectsFromArray: aBinarVectorOperation.secondOperand];
	}
	
	NSMutableArray *topVariableOperands = [[NSMutableArray alloc] init];
	NSMutableArray *precomputedVariableOperands = [[NSMutableArray alloc] init];
	NSMutableArray *otherVariableOperandOperations = [[NSMutableArray alloc] init];
	
	for (GLVariable *aVariable in operands) {
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
		NSLog(@"This operation depends only on precomputed variables. It will be computed now.");
        GLEquation *equation = [topVariables.lastObject equation];
        [equation solveForOperation: operation waitUntilFinished: YES];
        success = NO;
	}
	
	for (GLVariableOperation *anOperation in otherVariableOperandOperations) {
		success &= [self precomputeVariablesWithOperation: anOperation forTopVariables: topVariables];
	}	
	
	return success;
}



/************************************************/
/*		Internal Data Buffers					*/
/************************************************/

#pragma mark -
#pragma mark Internal Data Buffers
#pragma mark

@synthesize hasDataBuffer=_hasDataBuffer;
@synthesize totalMemoryBuffersAllocated;

- (BOOL) assignUnoptimizedMemoryBufferToVariable: (GLVariable *) variable forTopVariables: (NSArray *) topVariables bottomVariables: (NSArray *) bottomVariables
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
		[self.variableDataMap setObject: [variable.data copy] forKey: variable];
		return YES;
	}
	else if (variable.pendingOperations.count == 1)
	{
		if (![bottomVariables containsObject: variable])
		{	// It's not a top variable or a bottom variable and thus it needs a buffer.
			
			self.totalMemoryBuffersAllocated = self.totalMemoryBuffersAllocated + 1;
			
			// Here we fetch a data object, of appropriate size, for the variable.
			[self.variableDataMap setObject: [[GLMemoryPool sharedMemoryPool] dataWithLength: variable.dataBytes] forKey: variable];
			[self.hasDataBuffer addObject: variable];
		}
		
		// Now we work further up the tree and deal with the parents.
		GLVariableOperation *operation = variable.lastOperation;
		if (operation.operationType == kGLNullaryOperation)
		{
			return YES;
		}
		else if (operation.operationType == kGLUnaryOperation)
		{
			GLUnaryOperation *aUnaryOperation = (GLUnaryOperation *) operation;
			return [self assignUnoptimizedMemoryBufferToVariable: aUnaryOperation.operand forTopVariables: topVariables bottomVariables: bottomVariables];
		}
		else if (operation.operationType == kGLBinaryOperation)
		{
			GLBinaryOperation *aBinaryOperation = (GLBinaryOperation *) operation;
			
			// Do NOT put these methods into one if (a && b) call, otherwise the second one may not get executed. They
			// need to execute every time!
			BOOL a = [self assignUnoptimizedMemoryBufferToVariable: aBinaryOperation.firstOperand forTopVariables: topVariables bottomVariables: bottomVariables];
			BOOL b = [self assignUnoptimizedMemoryBufferToVariable: aBinaryOperation.secondOperand forTopVariables: topVariables bottomVariables: bottomVariables];
			
			return a && b;
		}
		else if (operation.operationType == kGLUnaryVectorOperation) {
			GLUnaryVectorOperation *aUnaryVectorOperation = (GLUnaryVectorOperation *) operation;
			BOOL success = YES;
			for (GLVariable *aVariable in aUnaryVectorOperation.operand) {
				success &= [self assignUnoptimizedMemoryBufferToVariable: aVariable forTopVariables: topVariables bottomVariables: bottomVariables];
			}
            
            // If one of the result variables of an operation does not get used, it will not appear in this graph. However, it still needs a data buffer as the operation needs to write somewhere.
            for (GLVariable *aVariable in aUnaryVectorOperation.result) {
				success &= [self assignUnoptimizedMemoryBufferToVariable: aVariable forTopVariables: topVariables bottomVariables: bottomVariables];
			}
            
			return success;
		} else if (operation.operationType == kGLBinaryVectorOperation) {
			GLBinaryVectorOperation *aBinaryVectorOperation = (GLBinaryVectorOperation *) operation;
			BOOL success = YES;
			for (GLVariable *aVariable in aBinaryVectorOperation.firstOperand) {
				success &= [self assignUnoptimizedMemoryBufferToVariable: aVariable forTopVariables: topVariables bottomVariables: bottomVariables];
			}
			for (GLVariable *aVariable in aBinaryVectorOperation.secondOperand) {
				success &= [self assignUnoptimizedMemoryBufferToVariable: aVariable forTopVariables: topVariables bottomVariables: bottomVariables];
			}
			return success;
		}
		
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

@synthesize finishedOperations;
@synthesize operationGroupArrayMap;
@synthesize operationSerialDependencyBlockArrayMap;
@synthesize operationParallelDependencyBlockArrayMap;
@synthesize bottomVariableGroup;
@synthesize totalTopVariablesCreated;

- (void) addGroupResponsibility: (dispatch_group_t) aGroup forOperation: (GLVariableOperation *) operation
{
	NSHashTable *array = [self.operationGroupArrayMap objectForKey: operation];
	if (!array) {
		array = [NSHashTable hashTableWithOptions: NSHashTableObjectPointerPersonality];
		[self.operationGroupArrayMap setObject: array forKey: operation];
	}
	
	[array addObject: aGroup];
}

- (NSHashTable *) groupResponsibilitiesForOperation: (GLVariableOperation *) operation
{
	NSHashTable *array = [self.operationGroupArrayMap objectForKey: operation];
	if (!array) {
		return [NSArray array];
	}
	
	return array;
}

- (NSHashTable *) allGroups
{
	NSHashTable *allGroups = [NSHashTable hashTableWithOptions: NSHashTableObjectPointerPersonality];
	for ( GLVariable *variable in self.operationGroupArrayMap ) {
		[allGroups unionHashTable: [self.operationGroupArrayMap objectForKey: variable]];
	}
		
	return allGroups;
}

- (void) addSerialResponsibility: (executionBlock) aBlock forOperation: (GLVariableOperation *) operation
{
	NSHashTable *array = [self.operationSerialDependencyBlockArrayMap objectForKey: operation];
	if (!array) {
		array = [NSHashTable hashTableWithOptions: NSHashTableObjectPointerPersonality];
		[self.operationSerialDependencyBlockArrayMap setObject: array forKey: operation];
	}
	
	[array addObject: [aBlock copy]];
}

- (NSHashTable *) serialResponsibilitiesForOperation: (GLVariableOperation *) operation
{
	NSHashTable *array = [self.operationSerialDependencyBlockArrayMap objectForKey: operation];
	if (!array) {
		return [NSHashTable hashTableWithOptions: NSHashTableObjectPointerPersonality];
	}
	
	return array;
}

- (void) addParallelResponsibility: (executionBlock) aBlock withGroup: (dispatch_group_t) group forOperation: (GLVariableOperation *) operation;
{
	NSMapTable *mapTable = [self.operationParallelDependencyBlockArrayMap objectForKey: operation];
	if (!mapTable) {
		mapTable = [NSMapTable mapTableWithKeyOptions: NSMapTableObjectPointerPersonality valueOptions:NSMapTableObjectPointerPersonality];
		[self.operationParallelDependencyBlockArrayMap setObject: mapTable forKey: operation];
	}
	
	[mapTable setObject: group forKey: [aBlock copy]];
}

- (NSMapTable *) parallelResponsibilitiesForOperation: (GLVariableOperation *) operation;
{
	NSMapTable *mapTable = [self.operationParallelDependencyBlockArrayMap objectForKey: operation];
	if (!mapTable) {
		mapTable = [NSMapTable mapTableWithKeyOptions: NSMapTableObjectPointerPersonality valueOptions:NSMapTableObjectPointerPersonality];
	}
	
	return mapTable;
}

- (BOOL) isTopOperationReady {
	GLVariableOperation * operation = (GLVariableOperation *) self;
	
	// If this test passes, then it means we're ready to construct the actual planBlock that needs to be run
	NSHashTable *groups = [self groupResponsibilitiesForOperation: operation];
	NSHashTable *serialExecutionBlocks = [self serialResponsibilitiesForOperation: operation];
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
	
	NSMutableArray *operands = [[NSMutableArray alloc] init];
	
	if (operation.operationType == kGLUnaryOperation) {
		GLUnaryOperation *aUnaryOperation = (GLUnaryOperation *) operation;
		[operands addObject: aUnaryOperation.operand];
	} else if (operation.operationType == kGLBinaryOperation) {
		GLBinaryOperation *aBinaryOperation = (GLBinaryOperation *) operation;
		[operands addObject: aBinaryOperation.firstOperand];
		[operands addObject: aBinaryOperation.secondOperand];
	} else if (operation.operationType == kGLUnaryVectorOperation) {
		GLUnaryVectorOperation *aUnaryVectorOperation = (GLUnaryVectorOperation *) operation;
		[operands addObjectsFromArray: aUnaryVectorOperation.operand];
	} else if (operation.operationType == kGLBinaryVectorOperation) {
		GLBinaryVectorOperation *aBinarVectorOperation = (GLBinaryVectorOperation *) operation;
		[operands addObjectsFromArray: aBinarVectorOperation.firstOperand];
		[operands addObjectsFromArray: aBinarVectorOperation.secondOperand];
	}
	
	// Sort the operands into categories
	NSMutableArray *topVariableOperands = [[NSMutableArray alloc] init];
	NSMutableArray *precomputedVariableOperands = [[NSMutableArray alloc] init];
	NSMutableArray *otherVariableOperandOperations = [[NSMutableArray alloc] init];
	
	for (GLVariable *aVariable in operands) {
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
		NSHashTable *groups = [self groupResponsibilitiesForOperation: operation];
		NSHashTable *serialExecutionBlocks = [self serialResponsibilitiesForOperation: operation];
		NSMapTable *parallelExecutionBlocks = [self parallelResponsibilitiesForOperation: operation];
		
		// The executionBlock for this operation will only be created when it has the expected number
		// of dependent children blocks and the expected number of groups.
		BOOL canCreate = ([self serialBlockCountForOperation: operation] == serialExecutionBlocks.count &&
						  [self parallelGroupCountForOperation: operation] == groups.count &&
						  [self parallelBlockCountForOperation: operation] == parallelExecutionBlocks.count);
		
		if (canCreate)
		{
			[self.finishedOperations addObject: operation];
			dispatch_queue_t globalQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
			dispatch_queue_t childrenQueue = self.childrenQueue;
            
			executionBlock childrenBlock = ^( NSArray *dataBuffers ) {
				// a. serial executions first
				for ( executionBlock anExecutionBlock in serialExecutionBlocks ) {
					dispatch_async( globalQueue, ^{
						anExecutionBlock( dataBuffers );
					});
				}
				
				// b. parallel executions second
                dispatch_async(childrenQueue, ^{
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
                });
			};
			
			executionBlock theExecutionBlock;
			
			if (operation.operationType == kGLNullaryOperation)
			{
				GLNullaryOperation *aNullaryOperation = (GLNullaryOperation *) operation;
				nullaryOperation operationBlock = [aNullaryOperation blockOperation];
				
				NSUInteger resultIndex = [self.allVariables indexOfObject: aNullaryOperation.result];
				
				theExecutionBlock = ^( NSArray *dataBuffers ) {
					NSMutableData *result = [dataBuffers objectAtIndex: resultIndex];
					
					// 1. Compute our own operation
					operationBlock( result );
					
					// 2. Send off the children
					childrenBlock( dataBuffers );
				};
			}
			else if (operation.operationType == kGLUnaryOperation)
			{
				GLUnaryOperation *aUnaryOperation = (GLUnaryOperation *) operation;
				unaryOperation operationBlock = [aUnaryOperation blockOperation];
				
				NSUInteger resultIndex = [self.allVariables indexOfObject: aUnaryOperation.result];
				NSUInteger operandIndex = [self.allVariables indexOfObject: aUnaryOperation.operand];
				
				theExecutionBlock = ^( NSArray *dataBuffers ) {					
					NSMutableData *result = [dataBuffers objectAtIndex: resultIndex];
					NSData *operand = [dataBuffers objectAtIndex: operandIndex];
					
					// 1. Compute our own operation
					operationBlock( result, operand );
					
					// 2. Send off the children
					childrenBlock( dataBuffers );
				};					
			}
			else if (operation.operationType == kGLBinaryOperation)
			{
				GLBinaryOperation *aBinaryOperation = (GLBinaryOperation *) operation;
				binaryOperation operationBlock = [aBinaryOperation blockOperation];
				
				NSUInteger resultIndex = [self.allVariables indexOfObject: aBinaryOperation.result];
				NSUInteger firstOperandIndex = [self.allVariables indexOfObject: aBinaryOperation.firstOperand];
				NSUInteger secondOperandIndex = [self.allVariables indexOfObject: aBinaryOperation.secondOperand];
				
				if (resultIndex == NSNotFound || firstOperandIndex == NSNotFound || secondOperandIndex == NSNotFound) {
					[NSException raise: @"Invalid index." format: @"The binary operation is malformed."];
				}
				
				theExecutionBlock = ^( NSArray *dataBuffers ) {
					NSMutableData *result = [dataBuffers objectAtIndex: resultIndex];
					NSData *firstOperand = [dataBuffers objectAtIndex: firstOperandIndex];
					NSData *secondOperand = [dataBuffers objectAtIndex: secondOperandIndex];
					
					// 1. Compute our own operation
					operationBlock( result, firstOperand, secondOperand );
					
					// 2. Send off the children
					childrenBlock( dataBuffers );
				};
			}
			else if (operation.operationType == kGLUnaryVectorOperation)
			{
				GLUnaryVectorOperation *aUnaryVectorOperation = (GLUnaryVectorOperation *) operation;
				unaryVectorOperation operationBlock = [aUnaryVectorOperation blockOperation];
				
				NSMutableArray *resultIndices = [[NSMutableArray alloc] init];
				NSMutableArray *operandIndices = [[NSMutableArray alloc] init];
				for (GLVariable *aVariable in aUnaryVectorOperation.result) {
                    NSUInteger anIndex = [self.allVariables indexOfObject: aVariable];
                    if (anIndex == NSNotFound) {
                        [NSException raise: @"Invalid index." format: @"The operation is malformed."];
                    } else {
                        [resultIndices addObject: @(anIndex)];
                    }     
				}
				for (GLVariable *aVariable in aUnaryVectorOperation.operand) {
                    NSUInteger anIndex = [self.allVariables indexOfObject: aVariable];
                    if (anIndex == NSNotFound) {
                        [NSException raise: @"Invalid index." format: @"The operation is malformed."];
                    } else {
                        [operandIndices addObject: @(anIndex)];
                    }
				}
				
				NSMutableArray *result = [[NSMutableArray alloc] init];
				NSMutableArray *operand = [[NSMutableArray alloc] init];
				
				theExecutionBlock = ^( NSArray *dataBuffers ) {
					for (NSNumber *anIndex in resultIndices) {
						[result addObject: [dataBuffers objectAtIndex: anIndex.unsignedIntegerValue]];
					}
					
					for (NSNumber *anIndex in operandIndices) {
						[operand addObject: [dataBuffers objectAtIndex: anIndex.unsignedIntegerValue]];
					}
                    
					// 1. Compute our own operation
					operationBlock( result, operand );
					
					// 2. Send off the children
					childrenBlock( dataBuffers );
                    
                    [result removeAllObjects];
                    [operand removeAllObjects];
				};
			}
			else if (operation.operationType == kGLBinaryVectorOperation)
			{
				GLBinaryVectorOperation *aBinaryVectorOperation = (GLBinaryVectorOperation *) operation;
				binaryVectorOperation operationBlock = [aBinaryVectorOperation blockOperation];
				
				NSMutableArray *resultIndices = [[NSMutableArray alloc] init];
				NSMutableArray *firstOperandIndices = [[NSMutableArray alloc] init];
				NSMutableArray *secondOperandIndices = [[NSMutableArray alloc] init];
				for (GLVariable *aVariable in aBinaryVectorOperation.result) {
					[resultIndices addObject: [NSNumber numberWithUnsignedInteger: [self.allVariables indexOfObject: aVariable]]];
				}
				for (GLVariable *aVariable in aBinaryVectorOperation.firstOperand) {
					[firstOperandIndices addObject: [NSNumber numberWithUnsignedInteger: [self.allVariables indexOfObject: aVariable]]];
				}
				for (GLVariable *aVariable in aBinaryVectorOperation.secondOperand) {
					[secondOperandIndices addObject: [NSNumber numberWithUnsignedInteger: [self.allVariables indexOfObject: aVariable]]];
				}
				
				NSMutableArray *result = [[NSMutableArray alloc] init];
				NSMutableArray *firstOperand = [[NSMutableArray alloc] init];
				NSMutableArray *secondOperand = [[NSMutableArray alloc] init];
				
				theExecutionBlock = ^( NSArray *dataBuffers ) {
					for (NSNumber *anIndex in resultIndices) {
						[result addObject: [dataBuffers objectAtIndex: anIndex.unsignedIntegerValue]];
					}
					
					for (NSNumber *anIndex in firstOperandIndices) {
						[firstOperand addObject: [dataBuffers objectAtIndex: anIndex.unsignedIntegerValue]];
					}
					
					for (NSNumber *anIndex in secondOperandIndices) {
						[secondOperand addObject: [dataBuffers objectAtIndex: anIndex.unsignedIntegerValue]];
					}
					
					// 1. Compute our own operation
					operationBlock( result, firstOperand, secondOperand );
					// 2. Send off the children
					childrenBlock( dataBuffers );
                    
                    [result removeAllObjects];
                    [firstOperand removeAllObjects];
                    [secondOperand removeAllObjects];
				};
			}
			
			if ( precomputedVariableOperands.count && !topVariableOperands.count && !otherVariableOperandOperations.count ) {
				NSLog(@"This operation depends only on precomputed variables. We have not yet implemented the appropriate optimization to deal with this.");
			} else if (operation.operationType == kGLNullaryOperation) {
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