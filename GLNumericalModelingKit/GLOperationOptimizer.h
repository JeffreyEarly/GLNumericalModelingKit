//
//  GLSimpleOperationOptimizer.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 4/18/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLEquation.h>
#import <GLNumericalModelingKit/GLFunction.h>
#import <GLNumericalModelingKit/GLVariableOperations.h>

typedef void (^executionBlock)(NSArray *);

// This class has extensive external documention describing its implementation.

@interface GLOperationOptimizer : NSObject

/************************************************/
/*		Primary Interface						*/
/************************************************/

#pragma mark -
#pragma mark Primary Interface
#pragma mark

// Initialize with the variables at the top of the tree, and the bottom of the tree.
- (GLOperationOptimizer *) initWithTopVariables: (NSArray *) topVariables bottomVariables: (NSArray *) bottomVariables;

// 5c. Returns an optimized unary vector block build from the execution blocks.
// This will return nil if the number of top variables isn't the same as the number
// of bottom variables.
@property(copy) variableOperation operationBlock;

// Ordered list of GLBuffer objects.
@property(strong) NSMutableArray *internalDataBuffers;

@property(strong) NSMutableArray *operandVariablePrototypes;

@property(strong) NSMutableArray *resultVariablePrototypes;

/************************************************/
/*		Create Execution Plan					*/
/************************************************/

#pragma mark -
#pragma mark Create Execution Plan
#pragma mark

// Step 1---returns YES if it was able to create an execution plan. If NO, then don't proceed.
- (BOOL) createExecutionPlan;

// Hunts down variables which need to be precomputed.
- (BOOL) precomputeVariablesWithOperation: (GLVariableOperation *) operation forTopVariables: (NSArray *) topVariables;

// This function does a first run through of the tree starting from the bottom node and
// working to the top. It adds up the number of responsibilities that each *operation* will
// take care of. 1) executing serial blocks, 2) executing parallel blocks, 3) exiting parallel groups
// The one exception is that the top-most responsibilities will be assigned to self.
- (BOOL) mapPreliminariesWithOperation: (GLVariableOperation *) operation forTopVariables: (NSArray *) topVariables;

/// A hash table containing operations which have had their connections fully mapped.
@property NSHashTable *finishedMappingOperations;

// The following store all the information about those dependencies.
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
/*		Assign Memory Buffers					*/
/************************************************/

#pragma mark -
#pragma mark Assign Memory Buffers
#pragma mark

// Step 2---returns YES if it was able to properly assign internal memory buffers. If NO, then don't proceed.
- (BOOL) assignMemoryBuffers;

// Creates a memory buffer for each internal variable and copies the data from pre-computed variables.
// It also creates an appropriately sized data object that each operation can use.
// This is the most inefficient memory assignment possible.
- (BOOL) assignUnoptimizedMemoryBufferToVariable: (GLFunction *) variable forTopVariables: (NSArray *) topVariables bottomVariables: (NSArray *) bottomVariables;

/// A hash table containing variables and operations which have a data buffer assigned.
@property NSHashTable *hasDataBuffer;

/// Map from precomputed variables to their NSData chunk.
@property NSMapTable *precomputedVariableDataMap;

/// Map from *internal* variables to a GLBuffer. Two variables may point to the same buffer, depending on the sophistication of the memory optimizer.
@property NSMapTable *internalVariableBufferMap;

/// Map from the *internal* buffers to a new GLBuffer.
@property NSMapTable *internalBufferBufferMap;

/// Total number of GLBuffer objects created for this operation.
@property NSUInteger totalMemoryBuffersAllocated;

/************************************************/
/*		Create Execution Blocks					*/
/************************************************/

#pragma mark -
#pragma mark Execution Plan Creation
#pragma mark

// Step 3---returns YES if it was able to create the execution blocks. If NO, then don't proceed.
- (BOOL) createExecutionBlocks;

// An ordered list of all the internal variables, followed by the internal buffers.
// This list maps one-to-one to the -internalDataBuffer array, which contains corresponding buffers.
@property(strong) NSMutableArray *internalVariablesAndBuffers;

// These take the -precomputedVariableDataMap and move it to one-to-one arrays.
@property(strong) NSMutableArray *precomputedVariables;
@property(strong) NSMutableArray *precomputedVariablesData;

// bottomVariables, topVariables, precomputedVariables, internalVariables and buffers;
@property(strong) NSMutableArray *allVariablesAndBuffers;

@property dispatch_queue_t childrenQueue;

@property NSUInteger totalTopVariablesCreated;

/// Array of all internal buffers the operations depend on.
@property NSMutableArray *internalBufferArray;

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
