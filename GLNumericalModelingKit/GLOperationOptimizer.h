//
//  GLSimpleOperationOptimizer.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 4/18/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLEquation.h>
#import <GLNumericalModelingKit/GLVariable.h>
#import <GLNumericalModelingKit/GLVariableOperations.h>

typedef void (^executionBlock)(NSArray *);

// This class has extensive external documention describing its implementation.

@interface GLOperationOptimizer : NSObject

/************************************************/
/*		Convenience Methods						*/
/************************************************/

#pragma mark -
#pragma mark Convenience Methods
#pragma mark

// Returns the top most variables
//+ (NSSet *) topVariablesFromVariable: (GLVariable *) variable;

//+ (unaryOperation) unaryOperationBlockFromTopVariables: (NSArray *) topVariables bottomVariables: (NSArray *) bottomVariables;
//
//+ (binaryOperation) binaryOperationBlockFromTopVariables: (NSArray *) topVariables bottomVariables: (NSArray *) bottomVariables;
//
//+ (unaryVectorOperation) unaryVectorOperationBlockFromTopVariables: (NSArray *) topVariables bottomVariables: (NSArray *) bottomVariables;

/************************************************/
/*		Primary Interface						*/
/************************************************/

#pragma mark -
#pragma mark Primary Interface
#pragma mark

// 1. Initialize with the variables at the top of the tree, and the bottom of the tree.
- (GLOperationOptimizer *) initWithTopVariables: (NSArray *) topVariables bottomVariables: (NSArray *) bottomVariables;

// 5c. Returns an optimized unary vector block build from the execution blocks.
// This will return nil if the number of top variables isn't the same as the number
// of bottom variables.
- (variableOperation) operationBlock;

// 6. Optional.
// The NSMutableData objects required by the block operation, ordered exactly
// the same as the internalVariables array. The block operation has an identical
// array, so they won't disappear until the block does, however if you want to return them
// to the sharedMemoryPool, you'll need to hold onto the objects in the array.
@property(strong) NSMutableArray *internalDataBuffers;

/************************************************/
/*		Internal Interface						*/
/************************************************/

#pragma mark -
#pragma mark Internal Interface
#pragma mark

// 2. Returns YES if it was able to create an execution plan. If NO, then don't proceed.
- (BOOL) createExecutionPlan;

// 3. Returns YES if it was able to properly assign internal memory buffers. If NO, then don't proceed.
- (BOOL) assignMemoryBuffers;

// 4. Returns YES if it was able to create the execution blocks. If NO, then don't proceed.
- (BOOL) createExecutionBlocks;

// An ordered list of all the internal variables.
@property(strong) NSMutableArray *internalVariables;

// bottomVariables, topVariables, internalVariables;
@property(strong) NSMutableArray *allVariables;

@property dispatch_queue_t childrenQueue;

/************************************************/
/*		Preliminary Mapping						*/
/************************************************/

#pragma mark -
#pragma mark Preliminary Mapping
#pragma mark

// Hunts down variables which need to be precomputed.
- (BOOL) precomputeVariablesWithOperation: (GLVariableOperation *) operation forTopVariables: (NSArray *) topVariables;

// This function does a first run through of the tree starting from the bottom node and
// working to the top. It adds up the number of responsibilities that each *operation* will
// take care of. 1) executing serial blocks, 2) executing parallel blocks, 3) exiting parallel groups
// The one exception is that the top-most responsibilities will be assigned to self.
- (BOOL) mapPreliminariesWithOperation: (GLVariableOperation *) operation forTopVariables: (NSArray *) topVariables;

@property NSMapTable *variableDataMap;
@property NSHashTable *finishedMappingOperations;
@property NSMapTable *operationSerialBlockCountMap;
@property NSMapTable *operationParallelBlockCountMap;
@property NSMapTable *operationParallelGroupCountMap;

- (void) incrementSerialBlockCountForOperation: (GLVariableOperation *) operation;
- (NSUInteger) serialBlockCountForOperation: (GLVariableOperation *) operation;

- (void) incrementParallelBlockCountForOperation: (GLVariableOperation *) operation;
- (NSUInteger) parallelBlockCountForOperation: (GLVariableOperation *) operation;

- (void) incrementParallelGroupCountForOperation: (GLVariableOperation *) operation;
- (NSUInteger) parallelGroupCountForOperation: (GLVariableOperation *) operation;

- (NSUInteger) totalDependencies;

/************************************************/
/*		Internal Data Buffers					*/
/************************************************/

#pragma mark -
#pragma mark Internal Data Buffers
#pragma mark

@property NSHashTable *hasDataBuffer;

@property NSUInteger totalMemoryBuffersAllocated;

// It also creates an appropriately sized data object that each operation can use.
// At the moment the data allocation is completely unoptimized.
- (BOOL) assignUnoptimizedMemoryBufferToVariable: (GLVariable *) variable forTopVariables: (NSArray *) topVariables bottomVariables: (NSArray *) bottomVariables;

/************************************************/
/*		Execution Plan Creation					*/
/************************************************/

#pragma mark -
#pragma mark Execution Plan Creation
#pragma mark

@property NSUInteger totalTopVariablesCreated;

// Returns yes its executionBlock and all the parents are successfully created, no otherwise.
- (BOOL) createExecutionBlockFromOperation: (GLVariableOperation *) operation forTopVariables: (NSArray *) topVariables bottomVariables: (NSArray *) bottomVariables;

- (void) addSerialResponsibility: (executionBlock) aBlock forOperation: (GLVariableOperation *) operation;
- (void) addParallelResponsibility: (executionBlock) aBlock withGroup: (dispatch_group_t) aGroup forOperation: (GLVariableOperation *) operation;
- (void) addGroupResponsibility: (dispatch_group_t) aGroup forOperation: (GLVariableOperation *) operation;

- (NSMutableArray *) serialResponsibilitiesForOperation: (GLVariableOperation *) operation;
- (NSMapTable *) parallelResponsibilitiesForOperation: (GLVariableOperation *) operation;
- (NSMutableArray *) groupResponsibilitiesForOperation: (GLVariableOperation *) operation;

- (NSMutableArray *) allGroups;

@property NSHashTable *finishedOperations;
@property NSMapTable *operationGroupArrayMap;
@property NSMapTable *operationSerialDependencyBlockArrayMap;
@property NSMapTable *operationParallelDependencyBlockArrayMap;
@property dispatch_group_t bottomVariableGroup;

@end
