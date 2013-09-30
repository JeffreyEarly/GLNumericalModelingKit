//
//  GLNetCDFFetchDataOperation.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 1/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import "GLVariableOperations.h"
#import "GLNetCDFVariable.h"
#import "GLLowLevelNetCDF.h"

// Critically, this operation does NOT depend on any other operations. It fetches from a NetCDF file in order
// to populate data buffer of the result variable---all the necessary information to do this is in the
// GLNetCDFVariable instance.
@interface GLNetCDFFetchDataOperation : GLVariableOperation

// This operation does not depend on the variable passed it. It will not keep a strong reference.
- (id) initWithNetCDFVariable: (GLNetCDFVariable *) variable indexRange: (NSArray *) ranges flatten: (BOOL) aFlag;
- (id) initWithResult: (GLVariable *) aResult netCDFVariable: (GLNetCDFVariable *) variable indexRange: (NSArray *) ranges flatten: (BOOL) aFlag;

@property(readwrite, strong, nonatomic) NSArray *theRanges;
@property(readwrite, assign, nonatomic) BOOL shouldFlatten;

@property(readwrite, strong) GLLowLevelNetCDF *file;
@property BOOL isComplex;
@property int variableID;
@property int imagpVariableID;

@property(readwrite, strong, nonatomic) GLVariable *result;

// To be called after the result and operands are set.
- (void) setupDependencies;

// To be called after the operation is complete.
- (void) tearDownDependencies;

@end
