//
//  NetCDFFile.h
//  FPSimulator
//
//  Created by Jeffrey Early on 7/6/09.
//  Copyright 2009 __MyCompanyName__. All rights reserved.
//

#import <Cocoa/Cocoa.h>

// Need to make this thread-safe for adding variables.

@class GLDimension, GLFunction, GLLowLevelNetCDF, GLNetCDFVariable, GLEquation;
@interface GLNetCDFFile : NSObject

/************************************************/
/*		Initialization							*/
/************************************************/

#pragma mark -
#pragma mark Initialization
#pragma mark

// Open an existing NetCDF file, or create a new one if none exists.
- (GLNetCDFFile *) initWithURL: (NSURL *) anURL forEquation: (GLEquation *) anEquation;

// Same as above, but with the option of overwriting the file if one already exists.
- (GLNetCDFFile *) initWithURL: (NSURL *) anURL forEquation: (GLEquation *) anEquation overwriteExisting: (BOOL) shouldOverwrite;

@property(readonly, strong) NSURL *URL;

@property(readonly, strong) GLEquation *equation;

// Blocks the current thread until all operations on the NetCDF file are finished.
- (void) waitUntilAllOperationsAreFinished;

// Clear the variables, close down the file, and be done.
- (void) close;

// Defaults to YES, in order to save space.
@property BOOL shouldAlwaysWriteSinglePrecision;

/************************************************/
/*		Dimensions, Variables & Attributes		*/
/************************************************/

#pragma mark -
#pragma mark Dimensions, Variables & Attributes
#pragma mark

// The read-only versions of the three basic properties of a NetCDF file.

// Dimensions are strongly retained and therefore remain in memory as long
// as the NetCDF file is still around in memory.
@property(readonly, strong) NSArray *dimensions;
@property(readonly, strong) NSArray *staticDimensions;
@property(readonly, strong) NSArray *mutableDimensions;

// Variables are not strongly retained. The instances returned may be different each time,
// although they will test true for -isEqual:.
// The strategy here allows the variables to be deallocated if they're not being used.
@property(readonly, strong) NSArray *variables;
@property(readonly, strong) NSArray *staticVariables;
@property(readonly, strong) NSArray *mutableVariables;

@property(readonly, strong) NSDictionary *globalAttributes;
- (void) setGlobalAttribute: (id) anObject forKey: (id) aKey;

// Dimensions can be added explicitly, or they will be added automatically when a variable is added.
// Dimensions must have a name. Dimensions automatically create a 'coordinate variable' behind the scenes.
// The NetCDF file creates a strong reference to the instance.
- (void) addDimension: (GLDimension *) dimension;

// Adds the variable, and any dimensions necessary, in one giant chunk. Variables must have a name.
// No strong reference to the variable is made. A GLNetCDFVariable is returned. No data is fetched.
// You cannot add a variable that doesn't have any valid data, unless it's a mutable variable.
// The returned variable will be mutable if dimensions dicate.
- (id) addVariable: (GLFunction *) variable;

- (GLDimension *) dimensionWithName: (NSString *) name;
- (GLNetCDFVariable *) variableWithName: (NSString *) name;

// Defaults to 0 (no compression). Level 1 and 2 are usually sufficient, but will slow the write speed.
@property(readwrite, assign) NSUInteger variableCompressionLevel;

/************************************************/
/*		Internal Representation					*/
/************************************************/

#pragma mark -
#pragma mark Internal Representation
#pragma mark

@property(readonly, strong) GLLowLevelNetCDF *file;

@end
