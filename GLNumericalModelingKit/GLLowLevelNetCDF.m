//
//  GLLowLevelNetCDF.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 10/22/11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#import "GLLowLevelNetCDF.h"
#import <netcdf.h>

#define ERR(e) {NSLog(@"Error: %s\n", nc_strerror(e));}

NSString *GLNumberOfDimensionsKey = @"GLNumberOfDimensionsKey";
NSString *GLNumberOfVariablesKey = @"GLNumberOfVariablesKey";
NSString *GLGlobalAttributesKey = @"GLGlobalAttributesKey";
NSString *GLUnlimitedDimensionsArrayKey = @"GLUnlimitedDimensionsArrayKey";

NSString *GLDimensionNameKey = @"GLDimensionNameKey";
NSString *GLDimensionLengthKey = @"GLDimensionLengthKey";
NSString *GLDimensionIDKey = @"GLDimensionIDKey";

NSString *GLVariableNameKey = @"GLVariableNameKey";
NSString *GLVariableTypeKey = @"GLVariableTypeKey";
NSString *GLVariableIDKey = @"GLVariableIDKey";
NSString *GLVariableDimensionsArrayKey = @"GLVariableDimensionsArrayKey";
NSString *GLVariableAttributesKey = @"GLVariableAttributesKey";

@interface GLLowLevelNetCDF () {
	int _fileID;
}

@property dispatch_queue_t readwriteQueue;

- (void) setAttributes: (NSDictionary *) attributes forVariableID: (int) varID;
- (NSDictionary *) attributesForVariableID: (int) varID totalExpected: (int) nAttributes;

- (void) setShortAttribute: (short) aValue forVariableID: (int) varID withKey: (NSString *) aKey;
- (void) setFloatAttribute: (float) aValue forVariableID: (int) varID withKey: (NSString *) aKey;
- (void) setDoubleAttribute: (double) aValue forVariableID: (int) varID withKey: (NSString *) aKey;
- (void) setStringAttribute: (NSString *) aValue forVariableID: (int) varID withKey: (NSString *) aKey;

- (id) valueForVariableID: (int) varID name: (NSString *) name type: (nc_type) type ofLength: (size_t) length;

@end

@implementation GLLowLevelNetCDF

/************************************************/
/*		Initialization							*/
/************************************************/

#pragma mark -
#pragma mark Initialization
#pragma mark

- (id) init
{
	if ((self=[super init])) {
		self.readwriteQueue = dispatch_queue_create( "com.EarlyInnovations.NetCDFReadWriteQueue", NULL );
	}
	return self;
}

- (void) waitUntilAllOperationsAreFinished
{
	__block int trivial=0;
	dispatch_sync( self.readwriteQueue, ^{
		trivial++;
	});
}

- (void) dealloc
{
	[self close];
}

- (void) createAtURL: (NSURL *) anURL withOptions: (int) options
{
	dispatch_async(self.readwriteQueue, ^{
		self.url = anURL;
		
		// Create the file. NC_SHARE makes NetCDF think we're planning on sharing the file while writing it.
		// This has the effect of not buffering the data and writing it immediately. This is a good thing!
		int retval;
		int anID;
		if ((retval = nc_create([self.url.path cStringUsingEncoding: NSASCIIStringEncoding], options, &anID))) ERR(retval);
		self.fileID = anID;
		
		// Exit define mode.
		if ((retval = nc_enddef(self.fileID))) ERR(retval);
	});
}


- (void) openAtURL: (NSURL *) anURL withOptions: (int) options
{
	dispatch_async(self.readwriteQueue, ^{
		self.url = anURL;
		
		int status;
		int anID;
		if ((status = nc_open([self.url.path cStringUsingEncoding: NSASCIIStringEncoding], NC_NOWRITE, &anID))) ERR(status);
		self.fileID = anID;
	});
}

- (void) close
{
	dispatch_sync(self.readwriteQueue, ^{
		if (self.fileID) {
			int retval;
			if ((retval = nc_close(self.fileID))) ERR(retval);
			self.fileID = 0;
		}
	});
}

@synthesize url;
@synthesize fileID=_fileID;
@synthesize readwriteQueue;

/************************************************/
/*		Inquiry									*/
/************************************************/

#pragma mark -
#pragma mark Inquiry
#pragma mark

- (NSDictionary *) fileProperties
{
	__block NSMutableDictionary *dictionary;
	
	dispatch_sync(self.readwriteQueue, ^{
		int status, ndims, nvars, nGlobalAttributes;
		if ((status = nc_inq( self.fileID, &ndims, &nvars, &nGlobalAttributes, NULL))) ERR(status);
		
		int nunlimdimsp;
		if ((status = nc_inq_unlimdims(self.fileID, &nunlimdimsp, NULL))) ERR(status);
		
		NSMutableArray *unlimitedDimensions = [NSMutableArray array];
		if (nunlimdimsp) {
			int *unlimdimidsp = malloc(nunlimdimsp*sizeof(int));
			if ((status = nc_inq_unlimdims(self.fileID, NULL, unlimdimidsp))) ERR(status);
			for ( int i = 0; i < nunlimdimsp; i++) {
				[unlimitedDimensions addObject: [NSNumber numberWithInt: unlimdimidsp[i]]];
			}
			free(unlimdimidsp);
		}
		
		NSDictionary *globalAttributes = [self attributesForVariableID: NC_GLOBAL totalExpected: nGlobalAttributes];
		
		dictionary = [NSMutableDictionary dictionary];
		[dictionary setObject: [NSNumber numberWithInt: ndims] forKey: GLNumberOfDimensionsKey];
		[dictionary setObject: [NSNumber numberWithInt: nvars] forKey: GLNumberOfVariablesKey];
		[dictionary setObject: globalAttributes forKey: GLGlobalAttributesKey];
		[dictionary setObject: unlimitedDimensions forKey: GLUnlimitedDimensionsArrayKey];
	});
	
	return dictionary;
}

- (NSDictionary *) propertiesOfDimensionWithID: (int) dimID
{
	__block NSMutableDictionary *dictionary;
	
	dispatch_sync(self.readwriteQueue, ^{
		char name[NC_MAX_NAME+1];
		size_t length;
		int status;
		
		if ((status = nc_inq_dim(self.fileID, dimID, name, &length))) ERR(status);
		NSString *nameString = [[NSString alloc] initWithCString: name encoding: NSASCIIStringEncoding];
		
		dictionary = [NSMutableDictionary dictionary];
		[dictionary setObject: nameString forKey: GLDimensionNameKey];
		[dictionary setObject: [NSNumber numberWithUnsignedLong: length] forKey: GLDimensionLengthKey];
		[dictionary setObject: [NSNumber numberWithInt: dimID] forKey: GLDimensionIDKey];
	});
	
	return dictionary;
}

- (NSDictionary *) propertiesOfVariableWithID: (int) varID
{
	__block NSMutableDictionary *dictionary;
	
	dispatch_sync(self.readwriteQueue, ^{
		char name[NC_MAX_NAME+1];
		nc_type varType;
		int status, nVarDims, varDims[NC_MAX_VAR_DIMS], nAttributes;
		
		if ((status = nc_inq_var(self.fileID, varID, name, &varType, &nVarDims, varDims, &nAttributes))) ERR(status);
		
		NSMutableArray *variableDimensions = [[NSMutableArray alloc] init];
		for ( int i=0; i < nVarDims; i++ ) {
			[variableDimensions addObject: [NSNumber numberWithInt: varDims[i]]];
		}
		
		NSString *nameString = [[NSString alloc] initWithCString: name encoding: NSASCIIStringEncoding];
		NSDictionary *attributesDictionary = [self attributesForVariableID: varID totalExpected: nAttributes];
		
		dictionary = [NSMutableDictionary dictionary];
		[dictionary setObject: nameString forKey: GLVariableNameKey];
		[dictionary setObject: variableDimensions forKey: GLVariableDimensionsArrayKey];
		[dictionary setObject: [NSNumber numberWithInt: varType] forKey: GLVariableTypeKey];
		[dictionary setObject: attributesDictionary forKey: GLVariableAttributesKey];
		[dictionary setObject: [NSNumber numberWithInt: varID] forKey: GLVariableIDKey];
	});
	
	return dictionary;
}

/************************************************/
/*		Creation								*/
/************************************************/

#pragma mark -
#pragma mark Creation
#pragma mark

- (void) addGlobalAttributes: (NSDictionary *) attributesDictionary
{
	dispatch_async(self.readwriteQueue, ^{
		int retval;
		nc_redef(self.fileID);
		[self setAttributes: attributesDictionary forVariableID: NC_GLOBAL];
		if ((retval = nc_enddef(self.fileID))) ERR(retval);	
	});
}

- (int) addDimensionWithName: (NSString *) aName length: (int) nPoints;
{
	__block int dimID;
	dispatch_sync(self.readwriteQueue, ^{
		nc_redef(self.fileID);
		
		int retval;
		if ((retval = nc_def_dim(self.fileID, [aName cStringUsingEncoding: NSASCIIStringEncoding], nPoints, &dimID))) ERR(retval);
		
		if ((retval = nc_enddef(self.fileID))) ERR(retval);		
	});
	return dimID;
}

- (int) addVariableOfType: (nc_type) type withName: (NSString *) aName dimensions: (NSArray *) dimArray compressionLevel: (NSUInteger) compLevel attributes: (NSDictionary *) attributesDictionary;
{
	__block int varID;
	dispatch_sync(self.readwriteQueue, ^{
		nc_redef(self.fileID);
		
		int retval, i=0;
		
		int *dimensionIDs = malloc( dimArray.count * sizeof(int) );
		for ( NSNumber *aDim in dimArray ) {
			dimensionIDs[i] = [aDim intValue];
			i++;
		}
		
		if ((retval = nc_def_var(self.fileID, [aName cStringUsingEncoding: NSASCIIStringEncoding], type, (int) dimArray.count, dimensionIDs, &varID))) ERR(retval);
		free(dimensionIDs);
		
		if (compLevel) {
			if ((retval = nc_def_var_deflate(self.fileID, varID, 0, 1, (int) compLevel))) ERR(retval);
		}
		
		[self setAttributes: attributesDictionary forVariableID: varID];
		
		if ((retval = nc_enddef(self.fileID))) ERR(retval);
	});
	return varID;
}

/************************************************/
/*		Add Data								*/
/************************************************/

#pragma mark -
#pragma mark Add Data
#pragma mark

- (void) setFloatData: (NSData *) data forVariableWithID: (int) variableID
{
	int fileID = self.fileID;
	dispatch_async(self.readwriteQueue, ^{
		int retval;
		if ((retval = nc_put_var_float( fileID,  variableID, data.bytes ))) ERR(retval);
	});
}

- (void) setDoubleData: (NSData *) data forVariableWithID: (int) variableID
{
	int fileID = self.fileID;
	dispatch_async(self.readwriteQueue, ^{
		int retval;
		if ((retval = nc_put_var_double( fileID,  variableID, data.bytes ))) ERR(retval);
	});
}

- (void) writeFloatData: (NSData *) data toVariableWithID: (int) variableID atIndexRange: (NSArray *) ranges
{
	int fileID = self.fileID;
	dispatch_async(self.readwriteQueue, ^{
		int retval;
		
		size_t *count = malloc( ranges.count * sizeof(size_t ));
		size_t *start = calloc( ranges.count, sizeof(size_t ));
		int i=0;
		for (NSValue *value in ranges) {
			NSRange range = [value rangeValue];
			start[i] = range.location;
			count[i] = range.length;
			i++;
		}
		
		if ((retval = nc_put_vara_float( fileID, variableID, start, count, data.bytes ))) ERR(retval);
		
		if ((retval = nc_sync( fileID )) ) ERR(retval);
		
		free(count); free(start);
		//NSLog(@"Wrote some float data for varID: %d at location: %lu.", variableID, [ranges[0] rangeValue].location );
//		NSMutableString *desc = [NSMutableString stringWithString: @"Added data at "];
//		for (NSValue *value in ranges) {
//			NSRange range = [value rangeValue];
//			[desc appendFormat: @"(%lu, %lu) ", range.location, range.length];
//		}
//		NSLog(@"%@", desc);
	});
}

- (void) writeDoubleData: (NSData *) data toVariableWithID: (int) variableID atIndexRange: (NSArray *) ranges
{
	int fileID = self.fileID;
	dispatch_async(self.readwriteQueue, ^{
		int retval;
		
		size_t *count = malloc( ranges.count * sizeof(size_t ));
		size_t *start = calloc( ranges.count, sizeof(size_t ));
		int i=0;
		for (NSValue *value in ranges) {
			NSRange range = [value rangeValue];
			start[i] = range.location;
			count[i] = range.length;
			i++;
		}
		
		if ((retval = nc_put_vara_double( fileID, variableID, start, count, data.bytes ))) ERR(retval);
		
		if ((retval = nc_sync( fileID )) ) ERR(retval);
		
		free(count); free(start);
	});
}

/************************************************/
/*		Get Data								*/
/************************************************/

#pragma mark -
#pragma mark Get Data
#pragma mark

- (void) readFloatVariableWithID: (int) variableID intoData: (NSMutableData *) theData indexRange: (NSArray *) ranges
{
	__block NSMutableData *data = theData;
	int fileID = self.fileID;
	dispatch_sync(self.readwriteQueue, ^{
		size_t *count = malloc( ranges.count * sizeof(size_t ));
		size_t *start = calloc( ranges.count, sizeof(size_t ));
		int i=0;
		for (NSValue *value in ranges) {
			NSRange range = [value rangeValue];
			start[i] = range.location;
			count[i] = range.length;
			i++;
		}
		
		nc_get_vara_float( fileID, variableID, start, count, data.mutableBytes);
		
		free(count); free(start);
	});
}

- (void) readDoubleVariableWithID: (int) variableID intoData: (NSMutableData *) theData indexRange: (NSArray *) ranges
{
	__block NSMutableData *data = theData;
	int fileID = self.fileID;
	dispatch_sync(self.readwriteQueue, ^{
		size_t *count = malloc( ranges.count * sizeof(size_t ));
		size_t *start = calloc( ranges.count, sizeof(size_t ));
		int i=0;
		for (NSValue *value in ranges) {
			NSRange range = [value rangeValue];
			start[i] = range.location;
			count[i] = range.length;
			i++;
		}
		
		nc_get_vara_double( fileID, variableID, start, count, data.mutableBytes);
		
		free(count); free(start);
	});
}

- (void) readFloatVariableWithID: (int) variableID intoBuffer: (void *) buffer indexRange: (NSArray *) ranges
{
	int fileID = self.fileID;
	dispatch_sync(self.readwriteQueue, ^{
		size_t *count = malloc( ranges.count * sizeof(size_t ));
		size_t *start = calloc( ranges.count, sizeof(size_t ));
		int i=0;
		for (NSValue *value in ranges) {
			NSRange range = [value rangeValue];
			start[i] = range.location;
			count[i] = range.length;
			i++;
		}
		
		nc_get_vara_float( fileID, variableID, start, count, buffer);
		
		free(count); free(start);
	});
}

- (void) readDoubleVariableWithID: (int) variableID intoBuffer: (void *) buffer indexRange: (NSArray *) ranges
{
	int fileID = self.fileID;
	dispatch_sync(self.readwriteQueue, ^{
		size_t *count = malloc( ranges.count * sizeof(size_t ));
		size_t *start = calloc( ranges.count, sizeof(size_t ));
		int i=0;
		for (NSValue *value in ranges) {
			NSRange range = [value rangeValue];
			start[i] = range.location;
			count[i] = range.length;
			i++;
		}
		
		nc_get_vara_double( fileID, variableID, start, count, buffer);
		
		free(count); free(start);
	});
}

/************************************************/
/*		Private: Attribute Dictionaries			*/
/************************************************/

#pragma mark -
#pragma mark Private: Attribute Dictionaries
#pragma mark

// These private functions should NOT lock/syncronize
// The should NOT change netcdf modes
// Both of those functions should be done by higher level methods.

- (void) setAttributes: (NSDictionary *) attributes forVariableID: (int) varID
{
	for (NSString *key in attributes) {
		id object = [attributes objectForKey: key];
		if ( [[object class] isSubclassOfClass: [NSString class]] )
		{
			[self setStringAttribute: object forVariableID: varID withKey: key];
		}
		else if ( [[object class] isSubclassOfClass: [NSNumber class]] )
		{
			NSNumber *number = object;
			if (strcmp(number.objCType, "c") == 0 ||
				strcmp(number.objCType, "C") == 0 ||
				strcmp(number.objCType, "S") == 0)
			{
				[self setShortAttribute: number.shortValue forVariableID: varID withKey:key];
			}
			else if (strcmp(number.objCType, "f") == 0)
			{
				[self setFloatAttribute: number.floatValue forVariableID: varID withKey:key];
			}
			else
			{
				[self setDoubleAttribute: number.doubleValue forVariableID: varID withKey:key];
			}
		}
	}
}

- (NSDictionary *) attributesForVariableID: (int) varID totalExpected: (int) nAttributes
{
	int status;
	NSMutableDictionary *attributesDictionary = [[NSMutableDictionary alloc] init];
	for ( int i=0; i < nAttributes; i++ ) {
		char attributeName[NC_MAX_NAME+1];
		nc_type attributeType;
		size_t attributeLength;
		if ((status = nc_inq_attname( self.fileID, varID, i, attributeName))) ERR(status);
		if ((status = nc_inq_att(self.fileID, varID, attributeName, &attributeType, &attributeLength))) ERR(status);
		
		NSString *attributeNameString = [[NSString alloc] initWithCString: attributeName encoding: NSASCIIStringEncoding];
		
		id value = [self valueForVariableID: varID name: attributeNameString type: attributeType ofLength: attributeLength];
		if (value && attributeNameString) [attributesDictionary setValue: value forKey: attributeNameString];
	}
	return attributesDictionary;
}

/************************************************/
/*		Private: Attribute Accessors			*/
/************************************************/

#pragma mark -
#pragma mark Private: Attribute Accessors
#pragma mark

- (void) setShortAttribute: (short) aValue forVariableID: (int) varID withKey: (NSString *) aKey
{
	int retval;
	if ((retval = nc_put_att_short(self.fileID, varID, [aKey cStringUsingEncoding: NSASCIIStringEncoding], NC_SHORT, 1, &aValue))) ERR(retval);
}

- (void) setFloatAttribute: (float) aValue forVariableID: (int) varID withKey: (NSString *) aKey
{
	int retval;
	if ((retval = nc_put_att_float(self.fileID, varID, [aKey cStringUsingEncoding: NSASCIIStringEncoding], NC_FLOAT, 1, &aValue))) ERR(retval);
}

- (void) setDoubleAttribute: (double) aValue forVariableID: (int) varID withKey: (NSString *) aKey
{
	int retval;
	if ((retval = nc_put_att_double(self.fileID, varID, [aKey cStringUsingEncoding: NSASCIIStringEncoding], NC_DOUBLE, 1, &aValue))) ERR(retval);
}

- (void) setStringAttribute: (NSString *) aValue forVariableID: (int) varID withKey: (NSString *) aKey
{
	int retval;
	const char *text = [aValue cStringUsingEncoding: NSASCIIStringEncoding];
	size_t length = [aValue lengthOfBytesUsingEncoding: NSASCIIStringEncoding];
	if ((retval = nc_put_att_text(self.fileID, varID, [aKey cStringUsingEncoding: NSASCIIStringEncoding], length, text))) ERR(retval);
}

- (id) valueForVariableID: (int) varID name: (NSString *) name type: (nc_type) type ofLength: (size_t) length
{
	if (length == 0) return nil;
	
	id value;
	
	int status;
	if ( type == NC_CHAR )
	{
		char *text = malloc( (length+1) * sizeof(char));
		if ((status = nc_get_att_text(self.fileID, varID, [name cStringUsingEncoding: NSASCIIStringEncoding], text))) ERR(status);
		value = [NSString stringWithCString: text encoding: NSASCIIStringEncoding];
		value = [value substringWithRange: NSMakeRange(0, length)];
		free(text);
	}
	else if (type == NC_SHORT)
	{
		SInt16 *array = malloc( length * sizeof(SInt16) );
		if ((status = nc_get_att_short(self.fileID, varID, [name cStringUsingEncoding: NSASCIIStringEncoding], array))) ERR(status);
		if (length == 1) {
			value = [NSNumber numberWithShort: array[0]];
		} else {
			NSMutableArray *objectArray = [[NSMutableArray alloc] init];
			for ( int i=0; i<length; i++) {
				[objectArray addObject: [NSNumber numberWithShort: array[i]]];
			}
			value = objectArray;
		}
		free(array);
	}
	else if (type == NC_INT)
	{
		int *array = malloc( length * sizeof(int) );
		if ((status = nc_get_att_int(self.fileID, varID, [name cStringUsingEncoding: NSASCIIStringEncoding], array))) ERR(status);
		if (length == 1) {
			value = [NSNumber numberWithInt: array[0]];
		} else {
			NSMutableArray *objectArray = [[NSMutableArray alloc] init];
			for ( int i=0; i<length; i++) {
				[objectArray addObject: [NSNumber numberWithInt: array[i]]];
			}
			value = objectArray;
		}
		free(array);
	}
	else if (type == NC_FLOAT)
	{
		float *floats = malloc( length * sizeof(float) );
		if ((status = nc_get_att_float(self.fileID, varID, [name cStringUsingEncoding: NSASCIIStringEncoding], floats))) ERR(status);
		if (length == 1) {
			value = [NSNumber numberWithFloat: *floats];
		} else {
			NSMutableArray *floatsArray = [[NSMutableArray alloc] init];
			for ( int i=0; i<length; i++) {
				[floatsArray addObject: [NSNumber numberWithFloat: floats[i]]];
			}
			value = floatsArray;
		}
		free(floats);
	}
	else if (type == NC_DOUBLE)
	{
		double *floats = malloc( length * sizeof(double) );
		if ((status = nc_get_att_double(self.fileID, varID, [name cStringUsingEncoding: NSASCIIStringEncoding], floats))) ERR(status);
		if (length == 1) {
			value = [NSNumber numberWithDouble: *floats];
		} else {
			NSMutableArray *floatsArray = [[NSMutableArray alloc] init];
			for ( int i=0; i<length; i++) {
				[floatsArray addObject: [NSNumber numberWithDouble: floats[i]]];
			}
			value = floatsArray;
		}
		free(floats);
	}
	
	return value;
}

@end
