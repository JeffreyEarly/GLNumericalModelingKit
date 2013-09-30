//
//  GLOperationVisualizer.h
//  GLNumericalModelingKitCopy
//
//  Created by Jeffrey Early on 12/1/12.
//
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLEquation.h>
#import <GLNumericalModelingKit/GLVariable.h>
#import <GLNumericalModelingKit/GLVariableOperations.h>

@interface GLOperationVisualizer : NSObject

// 1. Initialize with the variables at the top of the tree, and the bottom of the tree.
- (GLOperationVisualizer *) initWithTopVariables: (NSArray *) topVariables bottomVariables: (NSArray *) bottomVariables;

// These are straight stores of GLVariables and GLOperations
@property(readwrite, strong, nonatomic) NSMutableArray *internalVariables;
@property(readwrite, strong, nonatomic) NSMutableArray *otherVariables;
@property(readwrite, strong, nonatomic) NSMutableArray *operations;

// These are arrays of graphvis strings

// Operations are responsible for starting other operations. So these edges define the flow control of the program.
@property(readwrite, strong, nonatomic) NSMutableArray *operationOperationEdges;

// Top variable and precomputed variable connection to operations.
// Also includes the bottom operation connection to their variables.
@property(readwrite, strong, nonatomic) NSMutableArray *otherVariableOperationEdges;

// Same as above, but for internal variables only.
@property(readwrite, strong, nonatomic) NSMutableArray *internalVariableOperationEdges;

// This is a 'shortcut' description where internal variables are bypassed, and edges are simply
// drawn from an operation's output, to another operations input.
@property(readwrite, strong, nonatomic) NSMutableArray *operationOperationInputEdges;



@property(readonly) NSString *graphvisDescription;

@end
