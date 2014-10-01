//
//  GLMemoryOptimizer.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 4/22/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import "GLMemoryOptimizer.h"
#import "GLMemoryPool.h"
#import "GLVariableOperations.h"

@implementation GLMemoryOptimizer
{
	NSArray *_topVariables;
	NSArray *_bottomVariables;
}

@synthesize topVariables = _topVariables;
@synthesize bottomVariables = _bottomVariables;
@synthesize urChildUrParentMap;
@synthesize urParentUrChildMap;
@synthesize urParentIntermediateVariableMap;
@synthesize urChildIntermediateVariableMap;

- (GLMemoryOptimizer *) initWithTopVariables: (NSArray *) topVariables bottomVariables: (NSArray *) bottomVariables
{
    if ((self=[super init])) {
		self.variableChildrenMap = [NSMapTable mapTableWithKeyOptions: NSMapTableObjectPointerPersonality valueOptions:NSMapTableObjectPointerPersonality];
		
        self.urChildUrParentMap = [NSMapTable mapTableWithKeyOptions: NSMapTableObjectPointerPersonality valueOptions:NSMapTableObjectPointerPersonality];
		self.urParentUrChildMap = [NSMapTable mapTableWithKeyOptions: NSMapTableObjectPointerPersonality valueOptions:NSMapTableObjectPointerPersonality];
		self.urParentIntermediateVariableMap = [NSMapTable mapTableWithKeyOptions: NSMapTableObjectPointerPersonality valueOptions:NSMapTableObjectPointerPersonality];
		self.urChildIntermediateVariableMap = [NSMapTable mapTableWithKeyOptions: NSMapTableObjectPointerPersonality valueOptions:NSMapTableObjectPointerPersonality];	

        self.variableDataMap = [NSMapTable mapTableWithKeyOptions: NSMapTableObjectPointerPersonality valueOptions:NSMapTableObjectPointerPersonality];
		self.dataVariableMap = [NSMapTable mapTableWithKeyOptions: NSMapTableObjectPointerPersonality valueOptions:NSMapTableObjectPointerPersonality];
		
		self.variableSharableVariablesMap = [NSMapTable mapTableWithKeyOptions: NSMapTableObjectPointerPersonality valueOptions:NSMapTableObjectPointerPersonality];
		self.dataSharableVariablesMap = [NSMapTable mapTableWithKeyOptions: NSMapTableObjectPointerPersonality valueOptions:NSMapTableObjectPointerPersonality];
		
		self.topVariables = topVariables;
		self.bottomVariables = bottomVariables;		
	}
	return self;
}

/************************************************/
/*		Finding Children						*/
/************************************************/

#pragma mark -
#pragma mark Finding Children
#pragma mark

@synthesize variableChildrenMap;

- (void) addChild: (GLVariable *) child toVariable: (GLVariable *) parent
{
	NSMutableSet *set = [self.variableChildrenMap objectForKey: parent];
	if (!set) {
		set = [NSMutableSet set];
		[self.variableChildrenMap setObject: set forKey: parent];
	}
	[set addObject: child];
}

- (NSMutableSet *) childrenOfVariable: (GLVariable *) parent
{
	NSMutableSet *set = [self.variableChildrenMap objectForKey: parent];
	if (!set) {
		return [NSMutableSet set];
	}
	
	return set;
}

- (void) mapChildrenStartingWithVariable: (GLVariable *) variable
{
	if ( ![self isInternalVariable: variable] )
	{
		return;
	}
	else
	{				
		GLVariableOperation *operation = variable.lastOperation;
		if (operation.operationType == kGLUnaryOperation)
		{
			GLUnaryOperation *aUnaryOperation = (GLUnaryOperation *) operation;
			[self addChild: variable toVariable: aUnaryOperation.operand];
			[self mapChildrenStartingWithVariable: aUnaryOperation.operand];
		}
		else if (operation.operationType == kGLBinaryOperation)
		{
			GLBinaryOperation *aBinaryOperation = (GLBinaryOperation *) operation;
			[self addChild: variable toVariable: aBinaryOperation.firstOperand];
			[self addChild: variable toVariable: aBinaryOperation.secondOperand];
			[self mapChildrenStartingWithVariable: aBinaryOperation.firstOperand];
			[self mapChildrenStartingWithVariable: aBinaryOperation.secondOperand];
		}
	}
}

/************************************************/
/*		Finding Thread Joins/Splits             */
/************************************************/

#pragma mark -
#pragma mark Finding Thread Joins/Splits
#pragma mark

- (void) addVariableAndDecendents: (GLVariable *) variable toSet: (NSMutableSet *) reachableVariables
{
	[reachableVariables addObject: variable];
	
	NSSet *children = [self childrenOfVariable: variable];
	for ( GLVariable *child in children) {
		[self addVariableAndDecendents: child toSet: reachableVariables];
	}
}

- (void) addParentsAndVariable: (GLVariable *) variable toSet: (NSMutableSet *) reachableVariables
{
	if ( ![self isInternalVariable: variable] ) {
		return;
	}
	
	// This is a valid, reachable variable
	[reachableVariables addObject: variable];
	
	// So now let's traverse a bit higher up
	GLVariableOperation *operation = variable.lastOperation;
	if (operation.operationType == kGLUnaryOperation)
	{
		GLUnaryOperation *aUnaryOperation = (GLUnaryOperation *) operation;
		[self addParentsAndVariable: aUnaryOperation.operand toSet: reachableVariables];
	}
	else if (operation.operationType == kGLBinaryOperation)
	{
		GLBinaryOperation *aBinaryOperation = (GLBinaryOperation *) operation;
		[self addParentsAndVariable: aBinaryOperation.firstOperand toSet: reachableVariables];
		[self addParentsAndVariable: aBinaryOperation.secondOperand toSet: reachableVariables];
	}
	else if (operation.operationType == kGLUnaryVectorOperation)
	{
		NSLog(@"We cannot handle unary vector operations at this time. This will fail.");
	}
}

- (void) removeParentsAndVariable: (GLVariable *) variable fromSet: (NSMutableSet *) reachableVariables
{
	if ( ![self isInternalVariable: variable] ) {
		return;
	}
	
	// This is a valid, reachable variable
	[reachableVariables removeObject: variable];
	
	// So now let's traverse a bit higher up
	GLVariableOperation *operation = variable.lastOperation;
	if (operation.operationType == kGLUnaryOperation)
	{
		GLUnaryOperation *aUnaryOperation = (GLUnaryOperation *) operation;
		[self removeParentsAndVariable: aUnaryOperation.operand fromSet: reachableVariables];
	}
	else if (operation.operationType == kGLBinaryOperation)
	{
		GLBinaryOperation *aBinaryOperation = (GLBinaryOperation *) operation;
		[self removeParentsAndVariable: aBinaryOperation.firstOperand fromSet: reachableVariables];
		[self removeParentsAndVariable: aBinaryOperation.secondOperand fromSet: reachableVariables];
	}
	else if (operation.operationType == kGLUnaryVectorOperation)
	{
		NSLog(@"We cannot handle unary vector operations at this time. This will fail.");
	}
}

- (GLVariable *) firstAnscestorOfVariable: (GLVariable *) variable byClearingSet: (NSMutableSet *) reachableVariables
{
	if ( ![self isInternalVariable: variable] ) {
		return nil;
	}
	
	GLVariableOperation *operation = variable.lastOperation;
	if ([reachableVariables containsObject: variable])
	{
		// Ha! We got one, now don't let the anyone else get at this variable or its parents.
		if (operation.operationType == kGLUnaryOperation) {
			GLUnaryOperation *aUnaryOperation = (GLUnaryOperation *) operation;
			[self removeParentsAndVariable: aUnaryOperation.operand fromSet: reachableVariables];
		}
		else if (operation.operationType == kGLBinaryOperation) {
			GLBinaryOperation *aBinaryOperation = (GLBinaryOperation *) operation;
			[self removeParentsAndVariable: aBinaryOperation.firstOperand fromSet: reachableVariables];
			[self removeParentsAndVariable: aBinaryOperation.secondOperand fromSet: reachableVariables];
		}
		
		[reachableVariables removeObject: variable];
		return variable;
	} 
	
	// We didn't find one, so move on up the tree
	if (operation.operationType == kGLUnaryOperation)
	{
		GLUnaryOperation *aUnaryOperation = (GLUnaryOperation *) operation;
		return [self firstAnscestorOfVariable: aUnaryOperation.operand byClearingSet: reachableVariables];
	}
	else if (operation.operationType == kGLBinaryOperation)
	{
		GLBinaryOperation *aBinaryOperation = (GLBinaryOperation *) operation;
		
		GLVariable *varA = [self firstAnscestorOfVariable: aBinaryOperation.firstOperand byClearingSet: reachableVariables];
		GLVariable *varB = [self firstAnscestorOfVariable: aBinaryOperation.secondOperand byClearingSet: reachableVariables];
		
		// The trick here is that if varB was a better match than varA, it will have a value.
		return varB ?: varA;
	}
	
	return nil;
}

- (GLVariable *) firstCommonAnscestorOfVariable: (GLVariable *) variable
{
    // First check to see if we've already mapped this.
    GLVariable *anscestor = [self.urChildUrParentMap objectForKey: variable];
    if (anscestor) {
        if ( (NSNull *) anscestor == [NSNull null]) {
            return nil;
        } else {
            return anscestor;
        }
    }
    
	NSMutableSet *intermediateVariables = [[NSMutableSet alloc] init];
	
    // Now check to see if it's even an internal variable
	if ( ![self isInternalVariable: variable] ) {
		return nil;
	}   
    
	// Okay, so let's find the first common anscestor of this variable.
	GLVariableOperation *operation = variable.lastOperation;
	if (operation.operationType == kGLUnaryOperation)
	{
		GLUnaryOperation *aUnaryOperation = (GLUnaryOperation *) operation;
		
		if ( ![self isInternalVariable: aUnaryOperation.operand] ) {
			return nil;
		}
		
		anscestor = aUnaryOperation.operand;
	}
	else if (operation.operationType == kGLBinaryOperation)
	{
		GLBinaryOperation *aBinaryOperation = (GLBinaryOperation *) operation;
		
		if ( ![self isInternalVariable: aBinaryOperation.firstOperand] ) {
			anscestor = aBinaryOperation.secondOperand;
		} else if ( ![self isInternalVariable: aBinaryOperation.secondOperand] ) {
			anscestor = aBinaryOperation.firstOperand;
		} else {
			NSMutableSet *reachableVariables = [[NSMutableSet alloc] init];
			[self addParentsAndVariable: aBinaryOperation.firstOperand toSet:reachableVariables];
			anscestor = [self firstAnscestorOfVariable: aBinaryOperation.secondOperand byClearingSet: reachableVariables];
			[intermediateVariables unionSet: reachableVariables];
			
			reachableVariables = [[NSMutableSet alloc] init];
			[self addParentsAndVariable: aBinaryOperation.secondOperand toSet:reachableVariables];
			anscestor = [self firstAnscestorOfVariable: aBinaryOperation.firstOperand byClearingSet: reachableVariables];
			[intermediateVariables unionSet: reachableVariables];
		}
		// At this point reachableVariables contains the collection of variables internal to the thread-split and thread-join.
		// Any variable in this set can't be used by thread-split variable (but can be by the thread-join, if it can do in-place).	
	}
	else if (operation.operationType == kGLUnaryVectorOperation)
	{
		return nil;
		NSLog(@"We cannot handle unary vector operations at this time. This will fail.");
	}
	
    if (anscestor) {
        [self.urChildUrParentMap setObject: anscestor forKey: variable];
		[self.urChildIntermediateVariableMap setObject: intermediateVariables forKey: variable];
    } else {
        [self.urChildUrParentMap setObject: [NSNull null] forKey: variable];
		[self.urChildIntermediateVariableMap setObject: intermediateVariables forKey: variable];
    }
	
	if (anscestor) {
		GLVariable *existingChild = [self.urParentUrChildMap objectForKey: anscestor];
		if (!existingChild) {
			[self.urParentUrChildMap setObject: variable forKey: anscestor];
			[self.urParentIntermediateVariableMap setObject: intermediateVariables forKey: anscestor];
		} else {
			if ( [intermediateVariables member: existingChild] )
			{	// If the child identified is actually an intermediate variable, then we found a larger set.
				[self.urParentUrChildMap setObject: variable forKey: anscestor];
				[self.urParentIntermediateVariableMap setObject: intermediateVariables forKey: anscestor];
			}
		}
		
//		NSLog(@"Parent: %@ (<0x%lx>) - Child: %@ (<0x%lx>) - Intermediates: %ld", NSStringFromClass([anscestor.lastOperation class]), (NSUInteger) anscestor, NSStringFromClass([variable.lastOperation class]), (NSUInteger) variable, intermediateVariables.count);

	}
    
	return anscestor;
}

- (void) mapAllGraphCyclesStartingWithVariable: (GLVariable *) variable
{
	if ( ![self isInternalVariable: variable] )
	{
		return;
	}
	else
	{
		[self firstCommonAnscestorOfVariable: variable];
		
		GLVariableOperation *operation = variable.lastOperation;
		if (operation.operationType == kGLUnaryOperation)
		{
			GLUnaryOperation *aUnaryOperation = (GLUnaryOperation *) operation;
			[self mapAllGraphCyclesStartingWithVariable: aUnaryOperation.operand];
		}
		else if (operation.operationType == kGLBinaryOperation)
		{
			GLBinaryOperation *aBinaryOperation = (GLBinaryOperation *) operation;
			[self mapAllGraphCyclesStartingWithVariable: aBinaryOperation.firstOperand];
			[self mapAllGraphCyclesStartingWithVariable: aBinaryOperation.secondOperand];
		}
	}
}

/************************************************************/
/*		Categorizing Reachable Variables					*/
/************************************************************/

#pragma mark -
#pragma mark Categorizing Reachable Variables
#pragma mark

@synthesize variableSharableVariablesMap;

- (void) mapAllReachableVariablesStartingWithVariable: (GLVariable *) variable
{
	if ( ![self isInternalVariable: variable] ) {
		return;
	}

	GLVariableOperation *operation = variable.lastOperation;
	
	// Have we mapped this one yet?
	if (![self.variableSharableVariablesMap objectForKey: variable])
	{
		// First look down:
		// A data buffer may be reused by the urChild and all its descendents...
		NSMutableSet *reachableVariables = [[NSMutableSet alloc] init];
		
		[reachableVariables addObject: variable];
		
		GLVariable *urChild = [self.urParentUrChildMap objectForKey: variable];
		if (urChild) {
			[self addVariableAndDecendents: urChild toSet: reachableVariables];
		}
		
		
		// ...but not by the urChild if it can't operate in place and its a direct descendent
		GLVariableOperation *urChildOperation = variable.lastOperation;
		if (urChildOperation && !urChildOperation.canOperateInPlace) {
			NSSet *children = [self childrenOfVariable: variable];
			if ( [children member: urChild] ) {
				[reachableVariables removeObject: urChild];
			}
		}
		
		// Now look up:
		// A data buffer may be reused by its reachable parents...
		GLVariableOperation *operation = variable.lastOperation;
		if (operation.operationType == kGLUnaryOperation) {
			GLUnaryOperation *aUnaryOperation = (GLUnaryOperation *) operation;
			[self addParentsAndVariable: aUnaryOperation.operand toSet:reachableVariables];
		}
		else if (operation.operationType == kGLBinaryOperation) {
			GLBinaryOperation *aBinaryOperation = (GLBinaryOperation *) operation;
			[self addParentsAndVariable: aBinaryOperation.firstOperand toSet:reachableVariables];
			[self addParentsAndVariable: aBinaryOperation.secondOperand toSet:reachableVariables];
		}
		
		// ...but not by an urParent to which we're an intermediate variable
		for (GLVariable *urChild in self.urChildIntermediateVariableMap) {
			NSSet *intermediateVariables = [self.urChildIntermediateVariableMap objectForKey: urChild];
			if ( [intermediateVariables member: variable] ) {
				GLVariable *urParent = [self.urChildUrParentMap objectForKey: urChild];
				if ( (NSNull *) urParent != [NSNull null] ) {
					[reachableVariables removeObject: urParent];
				}
			}
		}
		
		// and we may also not share with our parents if we can't operate in place.
		if (!operation.canOperateInPlace) {
			if (operation.operationType == kGLUnaryOperation) {
				GLUnaryOperation *aUnaryOperation = (GLUnaryOperation *) operation;
				if ([self isInternalVariable: aUnaryOperation.operand]) {
					[reachableVariables removeObject: aUnaryOperation.operand];
				}
			} else if (operation.operationType == kGLBinaryOperation) {
				GLBinaryOperation *aBinaryOperation = (GLBinaryOperation *) operation;
				if ([self isInternalVariable: aBinaryOperation.firstOperand]) {
					[reachableVariables removeObject: aBinaryOperation.firstOperand];
				}
				if ([self isInternalVariable: aBinaryOperation.secondOperand]) {
					[reachableVariables removeObject: aBinaryOperation.secondOperand];
				}
			}
		}
		
		[self.variableSharableVariablesMap setObject: reachableVariables forKey: variable];
	}
	
	// Walk up the chain.
	if (operation.operationType == kGLUnaryOperation)
	{
		GLUnaryOperation *aUnaryOperation = (GLUnaryOperation *) operation;
					
		[self mapAllReachableVariablesStartingWithVariable: aUnaryOperation.operand];
	}
	else if (operation.operationType == kGLBinaryOperation)
	{
		GLBinaryOperation *aBinaryOperation = (GLBinaryOperation *) operation;
		
		[self mapAllReachableVariablesStartingWithVariable: aBinaryOperation.firstOperand];
		[self mapAllReachableVariablesStartingWithVariable: aBinaryOperation.secondOperand];
	}
}


/************************************************/
/*		Diagnostic Functions                    */
/************************************************/

#pragma mark -
#pragma mark Diagnostic Functions
#pragma mark

- (BOOL) isInternalVariable: (GLVariable *) variable
{
	if ([self.topVariables containsObject: variable]) {
		return NO;
	} else if (variable.pendingOperations.count == 0) {
		return NO;
	} else if (variable.pendingOperations.count > 1) {
        NSLog(@"Too many pending operations for variable: %@", [variable description]);
		return NO;
	}
	
	return YES;
}

- (BOOL) canVariable: (GLVariable *) variable useDataBuffer: (NSMutableData *) dataBuffer
{
	// The data is explicitly sharable with this variable.
	NSSet *dataSharable = [self.dataSharableVariablesMap objectForKey: dataBuffer];
	
	if ( [dataSharable member: variable] ) {
		// The variables that this variable is allowed to share with---it only includes variables higher up the graph.
		NSSet *variableReachableSharableMap = [self.variableSharableVariablesMap objectForKey: variable];
		// The variables that are using the data
		NSSet *variablesUsingData = [self.dataVariableMap objectForKey: dataBuffer];
		
		return [variablesUsingData isSubsetOfSet: variableReachableSharableMap];
	}
	
	return NO;
}

// Returns all existing data buffers usable by the variable, if even they're too small.
- (NSHashTable *) existingDataBuffersAvailableToVariable: (GLVariable *) variable
{
	NSHashTable *availableDataBuffers = [NSHashTable hashTableWithOptions: NSHashTableObjectPointerPersonality];
	for (NSMutableData *data in self.dataVariableMap) {
		if ([self canVariable: variable useDataBuffer:data]) {
			[availableDataBuffers addObject: data];
		}
	}
	
	return availableDataBuffers;
}

// Returns the most appropriately sized existing data, if possible.
- (NSMutableData *) existingDataBufferOfAppropriateSizeForVariable: (GLVariable *) variable fromSet: (NSHashTable *) buffers
{	
	NSMutableData *bestData;
	for (NSMutableData *data in buffers) {
		if (data.length == variable.dataBytes) {
			bestData = data;
		} else if (data.length >= variable.dataBytes && !bestData) {
			bestData = data;
		}
	}
	
//	if (buffers.count && !bestData) {
//		NSLog(@"Damn, I totally would have taken it if it were big enough!");
//		bestData =[buffers anyObject];
//		[bestData setLength: variable.dataBytes];
//	}
	
	return bestData;
}

// Find all cycles.
// Reachable variables are either 1) explicitly shareable or 2) explicitly forbidden
// Can I use some existing variable?
// Are all the variables using that data reachable to me? Nope, I could be above them.
// If they are all reachable to me and explicitly shareable, I'm okay
// If I'm not explicitly forbidden, then it's okay.

/************************************************/
/*		Memory Buffers Offerings                */
/************************************************/

#pragma mark -
#pragma mark Memory Buffers Offerings
#pragma mark

- (void) assignDataBuffer: (NSMutableData *) dataBuffer toVariable: (GLVariable *) variable
{
	[self.variableDataMap setObject: dataBuffer forKey: variable];
	
//	NSLog(@"Initial Data buffer (<0x%lx>) - Allowed: %ld", (NSUInteger) dataBuffer,[[self.dataSharableVariableMap objectForKey: dataBuffer] count]);
	
	// Add this variable to the list of variables that is using this particular data buffer.
	NSMutableSet *variables = [self.dataVariableMap objectForKey: dataBuffer];
	if (!variables) {
		variables = [[NSMutableSet alloc] init];
		[self.dataVariableMap setObject: variables forKey: dataBuffer];
	}
	[variables addObject: variable];
	
	NSSet *additionalSharables = [self.variableSharableVariablesMap objectForKey: variable];
	
	// Now we need to add the explicitly shareable variables to this set.
	NSMutableSet *sharableVariables = [self.dataSharableVariablesMap objectForKey: dataBuffer];
	if (!sharableVariables) {
		sharableVariables = [[NSMutableSet alloc] init];
		[self.dataSharableVariablesMap setObject: sharableVariables forKey: dataBuffer];
		[sharableVariables unionSet: additionalSharables];
	} else {
		[sharableVariables intersectSet: additionalSharables];
	}

//	NSLog(@"%@ (<0x%lx>) - Allowed: %ld", NSStringFromClass([variable.lastOperation class]), (NSUInteger) variable, additionalSharables.count);
//	NSLog(@"Final Data buffer (<0x%lx>) - Allowed: %ld", (NSUInteger) dataBuffer,[[self.dataSharableVariableMap objectForKey: dataBuffer] count]);
	
}
/************************************************/
/*		Assigning Memory Buffers                */
/************************************************/

#pragma mark -
#pragma mark Assigning Memory Buffers
#pragma mark

@synthesize variableDataMap;
@synthesize dataVariableMap;

@synthesize dataSharableVariablesMap;

- (BOOL) assignMemoryBufferToParentsAndVariable: (GLVariable *) variable
{
	if ([self.topVariables containsObject: variable])
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
		// Now we work further up the tree and deal with the parents.
		GLVariableOperation *operation = variable.lastOperation;
		
		if ( ![self.variableDataMap objectForKey: variable] && ![self.bottomVariables containsObject: variable] )
		{	// It's not a top variable, bottom variable, and doesn't already have a buffer
			NSHashTable *buffers = [self existingDataBuffersAvailableToVariable: variable];
			NSMutableData *dataBuffer = [self existingDataBufferOfAppropriateSizeForVariable: variable fromSet: buffers];
			if (!dataBuffer)
			{
				static int count=0;
				count++;
				NSLog(@"Total internal memory buffers allocated: %d", count);
				
				// Here we fetch a data object, of appropriate size, for the variable.
				dataBuffer = [[GLMemoryPool sharedMemoryPool] dataWithLength: variable.dataBytes];
			}
			
			[self assignDataBuffer: dataBuffer toVariable: variable];
		}
		
		if (operation.operationType == kGLUnaryOperation)
		{
			GLUnaryOperation *aUnaryOperation = (GLUnaryOperation *) operation;
			return [self assignMemoryBufferToParentsAndVariable: aUnaryOperation.operand];
		}
		else if (operation.operationType == kGLBinaryOperation)
		{
			GLBinaryOperation *aBinaryOperation = (GLBinaryOperation *) operation;
			
			// Do NOT put these methods into one if (a && b) call, otherwise the second one may not get executed. They
			// need to execute every time!
			BOOL a = [self assignMemoryBufferToParentsAndVariable: aBinaryOperation.firstOperand];
			BOOL b = [self assignMemoryBufferToParentsAndVariable: aBinaryOperation.secondOperand];
			
			return a && b;
		}
		else if (operation.operationType == kGLUnaryVectorOperation)
		{
			NSLog(@"Cannot appropriately assign memory buffers to a unaryVectorOperation. Not implemented.");
			return NO;
		}
		
		return NO;
		
	}
	else
	{
		NSLog(@"Trying to assign memory buffers, but pending operations for this variable are greater than 1. This violates our current assumptions.");
		return NO;
	}
}

/************************************************/
/*		Brute Force Optimization                */
/************************************************/

#pragma mark -
#pragma mark Brute Force Optimization
#pragma mark

- (void) computeAllPossibleDataThreads
{
	NSMutableSet *internalVariables = [[NSMutableSet alloc] init];
	for (GLVariable *variable in self.bottomVariables) {
		[self addParentsAndVariable: variable toSet: internalVariables];
	}
	
	int i = 0;
	for (GLVariable *variable in internalVariables)
	{
		NSLog(@"Times through: %d", i); i++;
		NSMutableSet *added = [[NSMutableSet alloc] init];
		NSMutableSet *allowed = [[NSMutableSet alloc] init];
		NSMutableSet *remaining = [[NSMutableSet alloc] init];
		
		NSSet *variableReachableSharableMap = [self.variableSharableVariablesMap objectForKey: variable];
		[added addObject: variable];
		[allowed unionSet: variableReachableSharableMap];
		[remaining unionSet: internalVariables];
		[remaining removeObject: variable];
		
		[self variablesAdded: added allowed:allowed remaining:remaining completedThreads: [[NSMutableSet alloc] init]];
	}
}


- (void) variablesAdded: (NSMutableSet *) added allowed: (NSMutableSet *) allowed remaining: (NSMutableSet *) remaining completedThreads: (NSMutableSet *) master
{
	if (!remaining.count) {
		NSLog(@"Found a set. Total threads: %ld", master.count);
		return;
	}
	
	// The untapped remaining variables are ones that we may be able to add to the set.
	NSMutableSet *untappedRemaining = [[NSMutableSet alloc] initWithSet: allowed];
	[untappedRemaining minusSet: added];
	[untappedRemaining intersectSet: remaining];
	
	// This set will now contain ALL possible variables that can be added to this set.
	NSMutableSet *untappedRemainingAllowed = [[NSMutableSet alloc] init];
	for (GLVariable *variable in untappedRemaining) {
		NSSet *variableReachableSharableMap = [self.variableSharableVariablesMap objectForKey: variable];
		if ( [added isSubsetOfSet: variableReachableSharableMap] ) {
			[untappedRemainingAllowed addObject: variable];
		}
	}
	
	if (untappedRemainingAllowed.count)
	{
		for (GLVariable *variable in untappedRemainingAllowed) {
			NSSet *variableReachableSharableMap = [self.variableSharableVariablesMap objectForKey: variable];
			
			NSMutableSet *newAdded = [[NSMutableSet alloc] initWithSet: added];
			NSMutableSet *newAllowed = [[NSMutableSet alloc] initWithSet: allowed];
			NSMutableSet *newRemaining = [[NSMutableSet alloc] initWithSet: remaining];
			
			[newAdded addObject: variable];
			[newAllowed intersectSet: variableReachableSharableMap];
			[newRemaining removeObject: variable];
			
			[self variablesAdded: newAdded allowed: newAllowed remaining: newRemaining completedThreads: master];
		}
	}
	else
	{
		// There are no other variables that can be added to this particular thread
		NSMutableSet *newMaster = [[NSMutableSet alloc] initWithSet: master];
		[newMaster addObject: added];
		
		// So start a new thread, with the same master list
		for ( GLVariable *variable in remaining )
		{
			NSMutableSet *newAdded = [[NSMutableSet alloc] init];
			NSMutableSet *newAllowed = [[NSMutableSet alloc] init];
			NSMutableSet *newRemaining = [[NSMutableSet alloc] init];
			
			NSSet *variableReachableSharableMap = [self.variableSharableVariablesMap objectForKey: variable];
			[newAdded addObject: variable];
			[newAllowed unionSet: variableReachableSharableMap];
			[newRemaining unionSet: remaining];
			[newRemaining removeObject: variable];
			
			[self variablesAdded: newAdded allowed: newAllowed remaining: newRemaining completedThreads: newMaster];
		}
	}
}



@end
