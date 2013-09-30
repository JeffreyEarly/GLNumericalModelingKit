//
//  GLMemoryOptimizer.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 4/22/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLVariable.h>
#import <GLNumericalModelingKit/GLVariableOperations.h>

@interface GLMemoryOptimizer : NSObject

// 1. Initialize with the variables at the top of the tree, and the bottom of the tree.
- (GLMemoryOptimizer *) initWithTopVariables: (NSArray *) topVariables bottomVariables: (NSArray *) bottomVariables;

@property(strong) NSArray *topVariables;
@property(strong) NSArray *bottomVariables;

/************************************************/
/*		Finding Children						*/
/************************************************/

#pragma mark -
#pragma mark Finding Children
#pragma mark

@property NSMapTable *variableChildrenMap;

- (void) addChild: (GLVariable *) child toVariable: (GLVariable *) parent;
- (NSMutableSet *) childrenOfVariable: (GLVariable *) parent;

- (void) mapChildrenStartingWithVariable: (GLVariable *) variable;

/************************************************/
/*		Finding Thread Joins/Splits             */
/************************************************/

#pragma mark -
#pragma mark Finding Thread Joins/Splits
#pragma mark

@property NSMapTable *urChildUrParentMap; // variable key, returns variable object
@property NSMapTable *urParentUrChildMap; // variable key, returns variable object
@property NSMapTable *urParentIntermediateVariableMap; // variable key, returns mutable set
@property NSMapTable *urChildIntermediateVariableMap; // variable key, returns mutable set

// If the varibles last operation is,
// 1) a top variable, precomputed variable or nil, it returns nil
// 2) a unary operation, it returns the operand
// 3) a binary operation, it returns the first variable reachable by both operands.
- (GLVariable *) firstCommonAnscestorOfVariable: (GLVariable *) variable;

// Adds all reachable internal anscestor variables to the set.
- (void) addParentsAndVariable: (GLVariable *) variable toSet: (NSMutableSet *) reachableVariables;

// Adds all reachable internal descendent variables to the set.
- (void) addVariableAndDecendents: (GLVariable *) variable toSet: (NSMutableSet *) reachableVariables;

// Removes all reachable internal variables from the set.
- (void) removeParentsAndVariable: (GLVariable *) variable fromSet: (NSMutableSet *) reachableVariables;

// Given a set of variables, this returns the first anscestor of the variable in that set.
// It does this by removing all variables from the set that it was able to reach. Any other path performing better
// will then return its variables.
- (GLVariable *) firstAnscestorOfVariable: (GLVariable *) variable byClearingSet: (NSMutableSet *) reachableVariables;

//	Creates the maps of thread joins and thread splits.
- (void) mapAllGraphCyclesStartingWithVariable: (GLVariable *) variable;


/************************************************************/
/*		Categorizing Reachable Variables					*/
/************************************************************/

#pragma mark -
#pragma mark Categorizing Reachable Variables
#pragma mark

// For each variable this returns the set of variables which are allowed to use the same data buffer.
@property NSMapTable *variableSharableVariablesMap;

// Creates the list of reachable-shareable and reachable-forbidden variables.
- (void) mapAllReachableVariablesStartingWithVariable: (GLVariable *) variable;

/************************************************/
/*		Diagnostic Functions                    */
/************************************************/

#pragma mark -
#pragma mark Diagnostic Functions
#pragma mark

// It's assumed that we're working up the tree, so it
// doesn't check whether the variable is a bottom variable (that would actually mess us up).
- (BOOL) isInternalVariable: (GLVariable *) variable;

// Returns yes if the variable is a member of this data's sharable variables OR
// if this variable's sharable variables is a superset of all the variables using this data.
- (BOOL) canVariable: (GLVariable *) variable useDataBuffer: (NSMutableData *) dataBuffer;

/************************************************/
/*		Assigning Memory Buffers                */
/************************************************/

#pragma mark -
#pragma mark Assigning Memory Buffers
#pragma mark

@property NSMapTable *variableDataMap; // variable key, returns variable object
@property NSMapTable *dataVariableMap; // variable key, return mutable set

// Starting from the bottom of the tree, this function is the last sweep---and makes certain
// all variables have a data buffer. When it creates a memory buffer, it tries to offer usage
// of that buffer to variables higher up the tree.
- (BOOL) assignMemoryBufferToParentsAndVariable: (GLVariable *) variable;

// For each data buffer this returns a set of variables that are allowed to use the data.
@property NSMapTable *dataSharableVariablesMap;

/************************************************/
/*		Brute Force Optimization                */
/************************************************/

#pragma mark -
#pragma mark Brute Force Optimization
#pragma mark

- (void) computeAllPossibleDataThreads;

@end
