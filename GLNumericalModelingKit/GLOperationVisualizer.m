//
//  GLOperationVisualizer.m
//  GLNumericalModelingKitCopy
//
//  Created by Jeffrey Early on 12/1/12.
//
//

// Good reference
// http://www.mactech.com/articles/mactech/Vol.25/25.01/IntroductiontoGraphviz/index.html

#import "GLOperationVisualizer.h"

@interface GLOperationVisualizer ()
- (BOOL) grabAllVariablesAndOperationsFromVariable: (GLVariable *) variable forTopVariables: (NSArray *) topVariables bottomVariables: (NSArray *) bottomVariables;
@end

@implementation GLOperationVisualizer

- (GLOperationVisualizer *) initWithTopVariables: (NSArray *) topVariables bottomVariables: (NSArray *) bottomVariables;
{
    if ((self=[super init])) {
        
        self.internalVariables = [[NSMutableArray alloc] init];
        self.otherVariables = [[NSMutableArray alloc] init];
        self.operations = [[NSMutableArray alloc] init];
        
        self.operationOperationEdges = [[NSMutableArray alloc] init];
        self.otherVariableOperationEdges = [[NSMutableArray alloc] init];
        self.internalVariableOperationEdges  = [[NSMutableArray alloc] init];
        self.operationOperationInputEdges = [[NSMutableArray alloc] init];
        
        BOOL success = YES;
        for (GLVariable *bottomVariable in bottomVariables) {
            success &=[self grabAllVariablesAndOperationsFromVariable: bottomVariable forTopVariables: topVariables bottomVariables: bottomVariables];
        }
        
        NSMutableSet *bottomOperations = [[NSMutableSet alloc] init];
        for (GLVariable *bottomVariable in bottomVariables) {
            if (bottomVariable.lastOperation){
                [bottomOperations addObject: bottomVariable.lastOperation];
                [self.otherVariableOperationEdges addObject: [NSString stringWithFormat: @"\t%ld -> %ld:in1:n [color=black, style=dotted];\n", (NSInteger) bottomVariable.lastOperation, (NSInteger) bottomVariable]];
            }
        }
        
        for (GLVariableOperation *operation in bottomOperations) {
            success &= [self mapPreliminariesWithOperation: operation forTopVariables: topVariables];
        }
        
    }
    return self;
}

- (NSString *) graphVisDescriptionFromOperation: (GLVariableOperation *) operation
{
    if (operation.operationType == kGLNullaryOperation) {
        return [NSString stringWithFormat: @"\t%ld [label=\"<op> %@\"];\n", (NSInteger) operation, operation.graphvisDescription];
    } else if (operation.operationType == kGLUnaryOperation) {
		return [NSString stringWithFormat: @"\t%ld [label=\"<in1> | <op> %@\"];\n", (NSInteger) operation, operation.graphvisDescription];
	} else if (operation.operationType == kGLBinaryOperation) {
		return [NSString stringWithFormat: @"\t%ld [label=\"<in1> | <op> %@ | <in2>\"];\n", (NSInteger) operation, operation.graphvisDescription];
	} else if (operation.operationType == kGLUnaryVectorOperation) {
		return [NSString stringWithFormat: @"\t%ld [label=\"<in1> | <op> %@\"];\n", (NSInteger) operation, operation.graphvisDescription];
	} else if (operation.operationType == kGLBinaryVectorOperation) {
		return [NSString stringWithFormat: @"\t%ld [label=\"<in1> | <op> %@ | <in2>\"];\n", (NSInteger) operation, operation.graphvisDescription];
	}
    
    return nil;
}

- (NSString *) graphvisDescription
{
    NSMutableString *graphDescription = [NSMutableString stringWithFormat: @"\ndigraph G {\n"];
    [graphDescription appendString: @"node [shape=ellipse,color=Red,fontname=Helvetica];\n\n"];
    for (GLVariable *variable in self.otherVariables) {
        [graphDescription appendFormat: @"\t%ld [label=\"%@\"];\n", (NSInteger) variable, variable.graphvisDescription];
    }
//    for (GLVariable *variable in self.internalVariables) {
//        [graphDescription appendFormat: @"\t%ld [label=\"%@\"]\n", (NSInteger) variable, variable.graphvisDescription];
//    }
    
    [graphDescription appendString: @"\n\nnode [shape=record,color=black,fontname=Helvetica];\n\n"];
    for (GLVariableOperation *operation in self.operations) {
        [graphDescription appendString: [self graphVisDescriptionFromOperation: operation]];
        //[graphDescription appendFormat: @"\t%ld [label=\"%@\"];\n", (NSInteger) operation, operation.graphvisDescription];
    }
    [graphDescription appendFormat: @"\tMasterOp [label=\"Master Operation\"];\n"];
    
    [graphDescription appendString: @"\n/*Edges defining the operation flow control*/\n"];
    for (NSString *descrip in self.operationOperationEdges) {
        [graphDescription appendString: descrip];
    }
    
    [graphDescription appendString: @"\n/*Top and bottom variable to operation input edges*/\n"];
    for (NSString *descrip in self.otherVariableOperationEdges) {
        [graphDescription appendString: descrip];
    }
    
    [graphDescription appendString: @"\n/*Operation output to operation input edges. Note that this shortcuts the intermediate variable*/\n"];
    for (NSString *descrip in self.operationOperationInputEdges) {
        [graphDescription appendString: descrip];
    }
    
    [graphDescription appendString: @"}\n"];
    
    return graphDescription;
}



- (BOOL) grabAllVariablesAndOperationsFromVariable: (GLVariable *) variable forTopVariables: (NSArray *) topVariables bottomVariables: (NSArray *) bottomVariables
{
	if ( [self.internalVariables containsObject: variable] || [self.otherVariables containsObject: variable] )
	{	// Already assigned a buffer to this variable, and therefore its parents as well.
		return YES;
	}
    else if ([topVariables containsObject: variable])
	{	// Top variables don't need buffers and we can exit the recusion.
        [self.otherVariables addObject: variable];
		return YES;
	}
	else if (variable.pendingOperations.count == 0)
	{
		// Precomputed variables are a dead end. It's values are fixed and already populated.
        [self.otherVariables addObject: variable];        
		return YES;
	}
	else if (variable.pendingOperations.count == 1)
	{
        if ([bottomVariables containsObject: variable]) {
            [self.otherVariables addObject: variable];
        } else {
            [self.internalVariables addObject: variable];
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
			return [self grabAllVariablesAndOperationsFromVariable: aUnaryOperation.operand forTopVariables: topVariables bottomVariables: bottomVariables];
		}
		else if (operation.operationType == kGLBinaryOperation)
		{
			GLBinaryOperation *aBinaryOperation = (GLBinaryOperation *) operation;
			
			// Do NOT put these methods into one if (a && b) call, otherwise the second one may not get executed. They
			// need to execute every time!
			BOOL a = [self grabAllVariablesAndOperationsFromVariable: aBinaryOperation.firstOperand forTopVariables: topVariables bottomVariables: bottomVariables];
			BOOL b = [self grabAllVariablesAndOperationsFromVariable: aBinaryOperation.secondOperand forTopVariables: topVariables bottomVariables: bottomVariables];
			
			return a && b;
		}
		else if (operation.operationType == kGLUnaryVectorOperation) {
			GLUnaryVectorOperation *aUnaryVectorOperation = (GLUnaryVectorOperation *) operation;
			BOOL success = YES;
			for (GLVariable *aVariable in aUnaryVectorOperation.operand) {
				success &= [self grabAllVariablesAndOperationsFromVariable: aVariable forTopVariables: topVariables bottomVariables: bottomVariables];
			}
			return success;
		} else if (operation.operationType == kGLBinaryVectorOperation) {
			GLBinaryVectorOperation *aBinaryVectorOperation = (GLBinaryVectorOperation *) operation;
			BOOL success = YES;
			for (GLVariable *aVariable in aBinaryVectorOperation.firstOperand) {
				success &= [self grabAllVariablesAndOperationsFromVariable: aVariable forTopVariables: topVariables bottomVariables: bottomVariables];
			}
			for (GLVariable *aVariable in aBinaryVectorOperation.secondOperand) {
				success &= [self grabAllVariablesAndOperationsFromVariable: aVariable forTopVariables: topVariables bottomVariables: bottomVariables];
			}
			return success;
		}
		
		return NO;
	} else {
		NSLog(@"Trying to record nodes, but pending operations for this variable are greater than 1. This violates our current assumptions.");
		return NO;
	}
}

- (NSString *) serialBlockGraphVisDesciptionFromOperation: (NSString *) fromOp toOperation: (NSString *) toOp
{
    return [NSString stringWithFormat: @"\t%@:op -> %@:op:n [color=black, style=bold];\n", fromOp, toOp];
}

- (NSString *) parallelBlockGraphVisDesciptionFromOperation: (NSString *) fromOp toOperation: (NSString *) toOp
{
    return [NSString stringWithFormat: @"\t%@:op -> %@:op:n [color=red, style=bold];\n", fromOp, toOp];
}

- (NSString *) parallelGroupGraphVisDesciptionFromOperation: (NSString *) fromOp toOperation: (NSString *) toOp
{
    return [NSString stringWithFormat: @"\t%@:op -> %@:op:n [color=blue, style=bold];\n", fromOp, toOp];
}

- (BOOL) mapPreliminariesWithOperation: (GLVariableOperation *) operation forTopVariables: (NSArray *) topVariables
{
	if ( [self.operations containsObject: operation] ) {
		return YES;
	} else {
		[self.operations addObject: operation];
	}
	
	// We now determine the responsibilities of each parent.
	// See the documentation for an explanation of the algorithm.
	BOOL success = NO;
	NSMutableArray *operands = [[NSMutableArray alloc] init];
	
    NSMutableArray *varOpDescriptions = [[NSMutableArray alloc] init];
    NSMutableArray *opOpDescriptions = [[NSMutableArray alloc] init];
	if (operation.operationType == kGLUnaryOperation) {
		GLUnaryOperation *aUnaryOperation = (GLUnaryOperation *) operation;
		[operands addObject: aUnaryOperation.operand];
        
        [varOpDescriptions addObject: [NSString stringWithFormat: @"\t%ld -> %ld:in1:n [color=black, style=dotted];\n", (NSInteger) aUnaryOperation.operand, (NSInteger) operation]];
        
        [opOpDescriptions addObject: [NSString stringWithFormat: @"\t%ld -> %ld:in1:n [color=black, style=dotted];\n", (NSInteger) aUnaryOperation.operand.lastOperation, (NSInteger) operation]];
	} else if (operation.operationType == kGLBinaryOperation) {
		GLBinaryOperation *aBinaryOperation = (GLBinaryOperation *) operation;
		[operands addObject: aBinaryOperation.firstOperand];
		[operands addObject: aBinaryOperation.secondOperand];
        
        [varOpDescriptions addObject: [NSString stringWithFormat: @"\t%ld -> %ld:in1:n [color=black, style=dotted];\n", (NSInteger) aBinaryOperation.firstOperand, (NSInteger) operation]];
        [varOpDescriptions addObject: [NSString stringWithFormat: @"\t%ld -> %ld:in2:n [color=black, style=dotted];\n", (NSInteger) aBinaryOperation.secondOperand, (NSInteger) operation]];
        
        [opOpDescriptions addObject: [NSString stringWithFormat: @"\t%ld -> %ld:in1:n [color=black, style=dotted];\n", (NSInteger) aBinaryOperation.firstOperand.lastOperation, (NSInteger) operation]];
        [opOpDescriptions addObject: [NSString stringWithFormat: @"\t%ld -> %ld:in2:n [color=black, style=dotted];\n", (NSInteger) aBinaryOperation.secondOperand.lastOperation, (NSInteger) operation]];
	} else if (operation.operationType == kGLUnaryVectorOperation) {
		GLUnaryVectorOperation *aUnaryVectorOperation = (GLUnaryVectorOperation *) operation;
		[operands addObjectsFromArray: aUnaryVectorOperation.operand];
        for (GLVariable *operandVar in aUnaryVectorOperation.operand) {
            [varOpDescriptions addObject: [NSString stringWithFormat: @"\t%ld -> %ld:in1:n [color=black, style=dotted];\n", (NSInteger) operandVar, (NSInteger) operation]];
            [opOpDescriptions addObject: [NSString stringWithFormat: @"\t%ld -> %ld:in1:n [color=black, style=dotted];\n", (NSInteger) operandVar.lastOperation, (NSInteger) operation]];
        }
	} else if (operation.operationType == kGLBinaryVectorOperation) {
		GLBinaryVectorOperation *aBinarVectorOperation = (GLBinaryVectorOperation *) operation;
		[operands addObjectsFromArray: aBinarVectorOperation.firstOperand];
		[operands addObjectsFromArray: aBinarVectorOperation.secondOperand];
        for (GLVariable *operandVar in aBinarVectorOperation.firstOperand) {
            [varOpDescriptions addObject: [NSString stringWithFormat: @"\t%ld -> %ld:in1:n [color=black, style=dotted];\n", (NSInteger) operandVar, (NSInteger) operation]];
            [opOpDescriptions addObject: [NSString stringWithFormat: @"\t%ld -> %ld:in1:n [color=black, style=dotted];\n", (NSInteger) operandVar.lastOperation, (NSInteger) operation]];
        }
        for (GLVariable *operandVar in aBinarVectorOperation.secondOperand) {
            [varOpDescriptions addObject: [NSString stringWithFormat: @"\t%ld -> %ld:in2:n [color=black, style=dotted];\n", (NSInteger) operandVar, (NSInteger) operation]];
            [opOpDescriptions addObject: [NSString stringWithFormat: @"\t%ld -> %ld:in2:n [color=black, style=dotted];\n", (NSInteger) operandVar.lastOperation, (NSInteger) operation]];
        }
	}
	
	NSMutableArray *topVariableOperands = [[NSMutableArray alloc] init];
	NSMutableArray *precomputedVariableOperands = [[NSMutableArray alloc] init];
	NSMutableArray *otherVariableOperandOperations = [[NSMutableArray alloc] init];
    NSMutableArray *otherVariableOperands = [[NSMutableArray alloc] init];
	
	for (NSUInteger i=0; i<operands.count; i++) {
        GLVariable *aVariable = operands[i];
		if ( [topVariables containsObject: aVariable] ) {
			[topVariableOperands addObject: aVariable];
            [self.otherVariableOperationEdges addObject: varOpDescriptions[i]];
		} else if ( aVariable.pendingOperations.count == 0 ) {
			[precomputedVariableOperands addObject: aVariable];
            [self.otherVariableOperationEdges addObject: varOpDescriptions[i]];
		} else {
            
            [self.internalVariableOperationEdges addObject: varOpDescriptions[i]];
            [self.operationOperationInputEdges addObject: opOpDescriptions[i]];
            
			if ( ![otherVariableOperandOperations containsObject: aVariable.lastOperation] ) {
				[otherVariableOperandOperations addObject: aVariable.lastOperation];
                [otherVariableOperands addObject: aVariable];
			}
		}
	}
	
	if ( precomputedVariableOperands.count && !topVariableOperands.count && !otherVariableOperandOperations.count ) {
		//NSLog(@"This operation depends only on precomputed variables. We have not yet implemented the appropriate optimization to deal with this.");
        [self.operationOperationEdges addObject: [self serialBlockGraphVisDesciptionFromOperation: @"MasterOp" toOperation: [NSString stringWithFormat: @"%ld", (NSInteger) operation]]];
	} else if ( operation.operationType == kGLNullaryOperation ) {
		//[self incrementSerialBlockCountForOperation: (GLVariableOperation *) self];
        // This operation depends on nothing, so a 'master' needs to get it started
        [self.operationOperationEdges addObject: [self serialBlockGraphVisDesciptionFromOperation: @"MasterOp" toOperation: [NSString stringWithFormat: @"%ld", (NSInteger) operation]]];
	} else if ( topVariableOperands.count && !otherVariableOperandOperations.count ) {
		//[self incrementSerialBlockCountForOperation: (GLVariableOperation *) self];
        // This operation depends only on top variables, so a 'master' needs to get it started
        [self.operationOperationEdges addObject: [self serialBlockGraphVisDesciptionFromOperation: @"MasterOp" toOperation: [NSString stringWithFormat: @"%ld", (NSInteger) operation]]];
	} else if ( otherVariableOperandOperations.count == 1 ) {
        //[self incrementSerialBlockCountForOperation: [otherVariableOperandOperations objectAtIndex: 0]];
        // This operation depends only on a single operand, so that operand needs to get it started
        [self.operationOperationEdges addObject: [self serialBlockGraphVisDesciptionFromOperation: [NSString stringWithFormat: @"%ld", (NSInteger) otherVariableOperandOperations[0]] toOperation: [NSString stringWithFormat: @"%ld", (NSInteger) operation]]];
	} else if ( otherVariableOperandOperations.count > 1 ){
		//[self incrementParallelBlockCountForOperation: [otherVariableOperandOperations objectAtIndex: 0]];
        [self.operationOperationEdges addObject: [self parallelBlockGraphVisDesciptionFromOperation: [NSString stringWithFormat: @"%ld", (NSInteger) otherVariableOperandOperations[0]] toOperation: [NSString stringWithFormat: @"%ld", (NSInteger) operation]]];
		for (NSUInteger i=1; i<otherVariableOperandOperations.count; i++) {
			//[self incrementParallelGroupCountForOperation: [otherVariableOperandOperations objectAtIndex: i]];
            [self.operationOperationEdges addObject: [self parallelGroupGraphVisDesciptionFromOperation: [NSString stringWithFormat: @"%ld", (NSInteger) otherVariableOperandOperations[i]] toOperation: [NSString stringWithFormat: @"%ld", (NSInteger) operation]]];
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

@end
