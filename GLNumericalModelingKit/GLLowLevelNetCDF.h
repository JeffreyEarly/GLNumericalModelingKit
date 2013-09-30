//
//  GLLowLevelNetCDF.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 10/22/11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <netcdf.h>

// Keys of the fileProperties dictionary
NSString *GLNumberOfDimensionsKey;
NSString *GLNumberOfVariablesKey;
NSString *GLGlobalAttributesKey;
NSString *GLUnlimitedDimensionsArrayKey;

// Keys of the dimensionProperties dictionary
NSString *GLDimensionNameKey;
NSString *GLDimensionLengthKey;
NSString *GLDimensionIDKey;

// Keys of the variableProperties dictionary
NSString *GLVariableNameKey;
NSString *GLVariableTypeKey;
NSString *GLVariableIDKey;
NSString *GLVariableDimensionsArrayKey;
NSString *GLVariableAttributesKey;

// This is a *simple* objective-c wrapper to NetCDF's C interface that requires no external classes.
// GLNetCDFFile provides a higher level interface and may be more appropriate to use for some cases.
//
// NetCDF's C interface is supposedly not thread safe, but this class uses a serial queue on NetCDF calls.
// Thus, even if you message an instance from multiple threads, you should be safe.
//
// This is by no means a complete implementation of the NetCDF interface, but it should cover most use cases.
// If it doesn't, then it's quite easily extendable.
//
// This class should not be using not be usin a mutex lock. Instead, it should use a serial dispatch queue.
// This would allow write operations to be asynchronous.

@interface GLLowLevelNetCDF : NSObject

/************************************************/
/*		Initialization							*/
/************************************************/

#pragma mark -
#pragma mark Initialization
#pragma mark

- (void) openAtURL: (NSURL *) anURL withOptions: (int) options;
- (void) createAtURL: (NSURL *) anURL withOptions: (int) options;
- (void) waitUntilAllOperationsAreFinished;
- (void) close;

@property(readwrite, strong, nonatomic) NSURL *url;
@property(readwrite, assign, nonatomic) int fileID;

/************************************************/
/*		Inquiry									*/
/************************************************/

#pragma mark -
#pragma mark Inquiry
#pragma mark

// Returns a dictionary with the four file property keys.
- (NSDictionary *) fileProperties;

// Returns a dictionary with the two dimension property keys for the given dimension id.
// Dimensions are given IDs 0...numDimensions-1, where numDimensions can be found in the file properties.
- (NSDictionary *) propertiesOfDimensionWithID: (int) dimensionID;

// Returns a dictionary with the four variable property keys.
// Variables are given IDs 0...numVariables-1, where numVariables can be found in the file properties.
- (NSDictionary *) propertiesOfVariableWithID: (int) variableID;

/************************************************/
/*		Creation								*/
/************************************************/

#pragma mark -
#pragma mark Creation
#pragma mark

// The dictionary should contain objects that NSStrings or NSNumbers.
// Everything else will be ignored.
- (void) addGlobalAttributes: (NSDictionary *) attributesDictionary;

// Returns the ID of the new dimension. NetCDF considers names unique.
// nPoints can be unlimitedDimension.
- (int) addDimensionWithName: (NSString *) aName length: (int) nPoints;

// The dimensions should be an array of NSNumbers corresponding to the existing dimensionID.
// The name needs to be unique amongst variables.
// The attributes dictionary should contain objects that NSStrings or NSNumbers.
// Returns the ID of the new variable.
- (int) addVariableOfType: (nc_type) type withName: (NSString *) aName dimensions: (NSArray *) dimArray attributes: (NSDictionary *) attributesDictionary;

/************************************************/
/*		Add Data								*/
/************************************************/

#pragma mark -
#pragma mark Add Data
#pragma mark

// Simple set some float or double float data for a particular variable ID. You'd better have the length right.
- (void) setFloatData: (NSData *) data forVariableWithID: (int) variableID;
- (void) setDoubleData: (NSData *) data forVariableWithID: (int) variableID;

// All ranges count the number of elements, not the number of bytes!
- (void) writeFloatData: (NSData *) data toVariableWithID: (int) variableID atIndexRange: (NSArray *) ranges;
- (void) writeDoubleData: (NSData *) data toVariableWithID: (int) variableID atIndexRange: (NSArray *) ranges;

/************************************************/
/*		Get Data								*/
/************************************************/

#pragma mark -
#pragma mark Get Data
#pragma mark

- (void) readFloatVariableWithID: (int) variableID intoData: (NSMutableData *) data indexRange: (NSArray *) ranges;
- (void) readDoubleVariableWithID: (int) variableID intoData: (NSMutableData *) data indexRange: (NSArray *) ranges;

- (void) readFloatVariableWithID: (int) variableID intoBuffer: (void *) buffer indexRange: (NSArray *) ranges;
- (void) readDoubleVariableWithID: (int) variableID intoBuffer: (void *) buffer indexRange: (NSArray *) ranges;
@end
