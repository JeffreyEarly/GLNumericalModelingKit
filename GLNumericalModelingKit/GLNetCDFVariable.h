//
//  GLNetCDFVariable.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 10/21/11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#import <GLNumericalModelingKit/GLVariable.h>

// This variable is backed by a NetCDF file. The variable depends on a GLNetCDFFetchDataOperation which
// will pull the data from the NetCDF file, if necessary. It may be the case only a subdomain
// of the variable is requested, or perhaps the variable is only to be appended to and never read. In
// these cases the fetch data operation will (hopefully) never be called and this variable simply stands
// as a proxy for the data.

@class GLLowLevelNetCDF;
@interface GLNetCDFVariable : GLVariable <NSCopying>

// Initializes a NetCDF backed variable with the same properties as the existingVariable.
// The uniqueIDs match. You will still need to set the file, variableID and -setupDependency.
+ (id) variableWithVariable: (GLVariable *) existingVariable;

// He we override the readonly attribute to allow the uniqueID to be set based on an existing variable.
@property(readwrite, assign, nonatomic) NSUInteger uniqueID;

// The variable needs a file to fetch from, and the variableID to fetch.
@property(readwrite, strong, nonatomic) GLLowLevelNetCDF *file;
@property int variableID;
@property int imagpVariableID;

// This adds the appropriate GLNetCDFFetchDataOperation as a dependency.
// Must be called after the file and variableID are set.
- (void) setupDependency;

@end


@interface GLMutableNetCDFVariable : GLNetCDFVariable <NSCopying>

/************************************************/
/*		Dimension Gymnastics					*/
/************************************************/

#pragma mark -
#pragma mark Dimension Gymnastics
#pragma mark

// These operations can only be applied along an existing mutable dimension.
// The result of these operation is that the mutable dimension will increase in length,
// and the variable will add new data.

// Both variables must be n-dimensional.
// The receiver's dimensions must include aDim (at mutableDimensionIndex).
// The dimensions at indices other than mutableDimensionIndex must have the same number of points.
// The dimensions at mutableDimensionIndex do not need to match. If the first operand's dimension is evenly spaced
// then it will be extended, otherwise the values from the second operand's dimension will be added.
- (void) concatenateWithVariable: (GLVariable *) otherVariable alongDimensionAtIndex: (NSUInteger) mutableDimensionIndex;

// The receiver must be n-dimensional and the other variable must be (n-1)-dimensional.
// The (n-1) dimensions of the two variables must have the same number of points.
// If the mutableDimension is evenly spaced, then it will be extended to length pointIndex+1, if necessary.
// If the mutableDimension is not evenly spaced, then it must already have the correct value.
- (void) concatenateWithLowerDimensionalVariable: (GLVariable *) otherVariable alongDimensionAtIndex: (NSUInteger) mutableDimensionIndex toIndex: (NSUInteger) pointIndex;

@end