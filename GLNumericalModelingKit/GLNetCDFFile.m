//
//  NetCDFFile.m
//  FPSimulator
//
//  Created by Jeffrey Early on 7/6/09.
//  Copyright 2009 __MyCompanyName__. All rights reserved.
//

#import "GLNetCDFFile.h"

#import "GLDimension.h"
#import "GLFunction.h"
#import "GLLowLevelNetCDF.h"
#import "GLEquation.h"
#import "GLNetCDFVariable.h"

static NSString *GLNetCDFSchemaVersionKey = @"GLNetCDFSchemaVersion";
static NSString *GLNetCDFSchemaDimensionIDKey = @"dimensionID";
static NSString *GLNetCDFSchemaIsCoordinateVariableKey = @"isCoordinateVariable";
static NSString *GLNetCDFSchemaIsPeridiocKey = @"isPeriodic";
static NSString *GLNetCDFSchemaMutableKey = @"isMutable";
static NSString *GLNetCDFSchemaIsEvenlySampledKey = @"isEvenlySampled";
static NSString *GLNetCDFSchemaDomainMinimumKey = @"domainMin";
static NSString *GLNetCDFSchemaSampleIntervalKey = @"sampleInterval";
static NSString *GLNetCDFSchemaIsFrequencyDomainKey = @"isFrequencyDomain";
static NSString *GLNetCDFSchemaIsComplexKey = @"isComplex";
static NSString *GLNetCDFSchemaProperNameKey = @"properName";
static NSString *GLNetCDFSchemaIsRealPartKey = @"isRealPart";
static NSString *GLNetCDFSchemaIsImaginaryPartKey = @"isImaginaryPart";
static NSString *GLNetCDFSchemaUnitsKey = @"units";
static NSString *GLNetCDFSchemaUniqueVariableIDKey = @"uniqueVariableID";

static NSString *GLNetCDFSchemaBasisFunctionKey = @"basisFunction";
static NSString *GLNetCDFSchemaGridTypeKey = @"gridType";
static NSString *GLNetCDFSchemaDomainLengthKey = @"domainLength";

static NSInteger GLCurrentNetCDFSchemaVersion = 12;

/************************************************/
/*		Private Methods 						*/
/************************************************/

#pragma mark -
#pragma mark Private Methods
#pragma mark

@interface GLNetCDFFile ()
@property(readwrite, strong) NSURL *URL;
@property(readwrite, strong) GLEquation *equation;

@property(readwrite, strong) NSDictionary *globalAttributes;

@property(readwrite, strong) GLLowLevelNetCDF *file;

@property(readwrite, strong) NSMapTable *dimensionDimensionIDMapTable;
@property(readwrite, strong) NSMapTable *dimensionVariableIDMapTable;

- (void) loadFromFileUsingGLNetCDFSchema: (NSDictionary *) fileProperties;
- (void) loadFromFileWithUnknownSource: (NSDictionary *) fileProperties;
- (GLDimension *) createDimensionFromProperties: (NSDictionary *) dimProperties variables: (NSMutableArray *) vars;
- (GLFunction *) createVariableFromProperties: (NSDictionary *) varProperties;

@end

@implementation GLNetCDFFile
{
	NSMutableDictionary *_globalAttributes;
	NSMutableArray *_dimensions;
	NSMutableArray *_staticDimensions;
	NSMutableArray *_mutableDimensions;
	NSMutableArray *_variables;
	NSMutableArray *_staticVariables;
	NSMutableArray *_mutableVariables;
}

/************************************************/
/*		Initialization							*/
/************************************************/

#pragma mark -
#pragma mark Initialization
#pragma mark

@synthesize file;
@synthesize URL;
@synthesize equation;

- (GLNetCDFFile *) initWithURL: (NSURL *) anURL forEquation: (GLEquation *) anEquation;
{
	return [self initWithURL: anURL forEquation: anEquation overwriteExisting: NO];
}

- (GLNetCDFFile *) initWithURL: (NSURL *) anURL forEquation: (GLEquation *) anEquation overwriteExisting: (BOOL) shouldOverwrite
{
	if (!anURL) {
		NSLog(@"GLNetCDFFile.m: attempting to initialize a file without an url!");
		return nil;
	} else if (!anEquation) {
		NSLog(@"GLNetCDFFile.m: attempting to initialize a file without an equation!");
	}
	
	if ( (self=[super init]) )
	{
		self.URL = anURL;
		self.equation = anEquation;
		
		_dimensions = [NSMutableArray array];
		_staticDimensions = [NSMutableArray array];
		_mutableDimensions = [NSMutableArray array];
		
		_variables = [NSMutableArray array];
		_staticVariables = [NSMutableArray array];
		_mutableVariables = [NSMutableArray array];
		
		_globalAttributes = [NSMutableDictionary dictionary];
		
		// We strongly hold on to the dimensions, but only do pointer comparisons, not isEqual.
		NSPointerFunctionsOptions options = NSPointerFunctionsStrongMemory | NSPointerFunctionsObjectPointerPersonality;
		self.dimensionDimensionIDMapTable = [NSMapTable mapTableWithKeyOptions:options valueOptions:options];
		self.dimensionVariableIDMapTable = [NSMapTable mapTableWithKeyOptions:options valueOptions:options];
		
		self.file = [[GLLowLevelNetCDF alloc] init];
		
		NSFileManager *fileManager = [[NSFileManager alloc] init];
		if ( shouldOverwrite || ![fileManager fileExistsAtPath: anURL.path] ) {
			[self.file createAtURL: self.URL withOptions: 0];
			[self setGlobalAttribute: [NSNumber numberWithInteger: GLCurrentNetCDFSchemaVersion] forKey: GLNetCDFSchemaVersionKey];
		} else {
			[self.file openAtURL: self.URL withOptions: 0];
			
			NSDictionary *fileProperties = self.file.fileProperties;
			self.globalAttributes = [[fileProperties objectForKey: GLGlobalAttributesKey] mutableCopy];
			if ([self.globalAttributes objectForKey: GLNetCDFSchemaVersionKey]) {
				[self loadFromFileUsingGLNetCDFSchema: fileProperties];
			} else {
				[self loadFromFileWithUnknownSource: fileProperties];
			}
		}
	}
	
	return self;
}

- (void) waitUntilAllOperationsAreFinished
{
	[self.file waitUntilAllOperationsAreFinished];
}

- (void) close
{
	for (GLMutableDimension *dimension in self.mutableDimensions) {
		[self stopObservingDimension: dimension];
	}
	[self.file close];
}

- (void) loadFromFileUsingGLNetCDFSchema: (NSDictionary *) fileProperties
{
	// If everything really is kosher, we can attach the variables directly.	
	int numVariables = [[fileProperties objectForKey: GLNumberOfVariablesKey] intValue];
	NSMutableArray *vars = [NSMutableArray array];
	for (int i=0; i<numVariables; i++)
	{
		[vars addObject: [self.file propertiesOfVariableWithID: i]];
	}
	
	// Find the coordinate variables, and build dimensions from them.
	for ( NSDictionary *variableProperties in vars )
	{
		NSDictionary *variableAttributes = variableProperties[GLVariableAttributesKey];
		
		if ( [variableAttributes[GLNetCDFSchemaIsCoordinateVariableKey] boolValue] )
		{
			int dimensionID = [variableAttributes[GLNetCDFSchemaDimensionIDKey] intValue];
			int variableID = [variableProperties[GLVariableIDKey] intValue];
			
			NSDictionary *dimensionProperties = [self.file propertiesOfDimensionWithID: dimensionID];
			int dimPoints = [[dimensionProperties objectForKey: GLDimensionLengthKey] intValue];
			GLGridType grid = [variableAttributes[GLNetCDFSchemaGridTypeKey] unsignedIntegerValue];
            
			GLDimension * dimension;
			Class GLDimensionClass = [variableAttributes[GLNetCDFSchemaMutableKey] boolValue] ? [GLMutableDimension class] : [GLDimension class];
            
            if (grid == kGLUnevenGrid)
            {
                NSArray *indexRange = [NSArray arrayWithObject: [NSValue valueWithRange: NSMakeRange(0, dimPoints)]];
				NSMutableData *buffer = [[NSMutableData alloc] initWithLength: dimPoints*sizeof(GLFloat)];
				if ( sizeof(GLFloat) == sizeof(double) ) {
					[self.file readDoubleVariableWithID: variableID intoData: buffer indexRange: indexRange];
				} else {
					[self.file readFloatVariableWithID: variableID intoData: buffer indexRange: indexRange];
				}
				dimension = [[GLDimension alloc] initWithNPoints: dimPoints values: buffer];
            }
            else
            {
                GLFloat domainLength = [variableAttributes[GLNetCDFSchemaDomainLengthKey] doubleValue];
				GLFloat domainMin = [variableAttributes[GLNetCDFSchemaDomainMinimumKey] doubleValue];
				
				dimension = [[GLDimensionClass alloc] initDimensionWithGrid: grid nPoints: dimPoints domainMin: domainMin length: domainLength];
            }
            
			if ( variableAttributes[GLNetCDFSchemaBasisFunctionKey] ) {
				dimension.basisFunction = [variableAttributes[GLNetCDFSchemaBasisFunctionKey] unsignedIntegerValue];
			}

			if (!dimension) {
                [NSException raise: @"GLNetCDFReadVariableException" format:@"GLNetCDFFile.m: something went horribly wrong. Unable to build a dimension from a known schema."];
				continue;
			}
			dimension.name = variableProperties[GLVariableNameKey];
			dimension.units = variableAttributes[GLNetCDFSchemaUnitsKey];
			
			[self.dimensionDimensionIDMapTable setObject: [NSNumber numberWithInt: dimensionID] forKey: dimension];
			[self.dimensionVariableIDMapTable setObject: [NSNumber numberWithInteger: variableID] forKey: dimension];
			
			[_dimensions addObject: dimension];
			if (dimension.isMutable) {
				[self startObservingDimension: (GLMutableDimension *)dimension];
				[_mutableDimensions addObject: dimension];
			} else {
				[_staticDimensions addObject: dimension];
			}
		}
	}
	
	// Now we build the variables that are NOT coordinate variables.
	NSMutableArray *alreadyCreated = [NSMutableArray array];
	for ( NSDictionary *variableProperties in vars )
	{
		NSDictionary *variableAttributes = variableProperties[GLVariableAttributesKey];
		
		if ( ![variableAttributes[GLNetCDFSchemaIsCoordinateVariableKey] boolValue] && ![alreadyCreated containsObject: variableProperties] )
		{
			int variableID = [variableProperties[GLVariableIDKey] intValue];
			
			NSArray *dimensionIDs = variableProperties[GLVariableDimensionsArrayKey];
			NSMutableArray *dimensions = [NSMutableArray array];
			
			for ( NSNumber *dimensionID in dimensionIDs) {
				for ( GLDimension *dimension in self.dimensionDimensionIDMapTable ) {
					if ([dimensionID isEqualToNumber: [self.dimensionDimensionIDMapTable objectForKey: dimension]]) {
						[dimensions addObject: dimension];
					}
				}
			}
			
			BOOL isComplex = [variableAttributes[GLNetCDFSchemaIsComplexKey] boolValue];
			GLDataFormat format = isComplex ? kGLSplitComplexDataFormat : kGLRealDataFormat;
			GLNetCDFVariable *netcdfVariable = [GLNetCDFVariable functionOfType: format withDimensions: dimensions forEquation: self.equation];
			netcdfVariable.file = self.file;
			netcdfVariable.name = variableProperties[GLVariableNameKey];
			netcdfVariable.units = variableAttributes[GLNetCDFSchemaUnitsKey];
			netcdfVariable.uniqueID = [variableAttributes[GLNetCDFSchemaUniqueVariableIDKey] unsignedIntegerValue];
			[netcdfVariable.metadata addEntriesFromDictionary: [variableProperties objectForKey: GLVariableAttributesKey]];
			
			if (isComplex)
			{
				netcdfVariable.name = variableAttributes[GLNetCDFSchemaProperNameKey];
				int otherID=-1;
				for ( NSDictionary *otherVariableProperties in vars ) {
					NSDictionary *otherVariableAttributes = variableProperties[GLVariableAttributesKey];
					if (otherVariableProperties != variableProperties && netcdfVariable.uniqueID == [otherVariableAttributes[GLNetCDFSchemaUniqueVariableIDKey] unsignedIntegerValue]) {
						otherID = (int) [vars indexOfObject: otherVariableProperties];
					}
				}
				
				if (otherID != -1) {
					if ( [[variableProperties objectForKey: GLNetCDFSchemaIsRealPartKey] boolValue] ) {
						netcdfVariable.variableID = variableID;
						netcdfVariable.imagpVariableID = otherID;
					} else {
						netcdfVariable.variableID = otherID;
						netcdfVariable.imagpVariableID = variableID;
					}
				} else {
					NSLog(@"GLNetCDFFile.m: unable to locate the other half of the complex variable.");
					netcdfVariable = nil;
				}				
			}
			else
			{
				netcdfVariable.variableID = variableID;
			}
			
			if (!netcdfVariable) {
				NSLog(@"GLNetCDFFile.m: something went horribly wrong. Unable to build a variable from a known schema.");
				continue;
			}
			
			[netcdfVariable setupDependency];
			[_variables addObject: netcdfVariable];
			if (netcdfVariable.isMutable) {
				[_mutableVariables addObject: netcdfVariable];
			} else {
				[_staticVariables addObject: netcdfVariable];
			}
		}
	}
}

- (void) loadFromFileWithUnknownSource: (NSDictionary *) fileProperties
{
	
	int numDimensions = [[fileProperties objectForKey: GLNumberOfDimensionsKey] intValue];
	int numVariables = [[fileProperties objectForKey: GLNumberOfVariablesKey] intValue];
	
	// Read all the dimension properties
	NSMutableArray *dims = [NSMutableArray array];
	for (int i=0; i<numDimensions; i++)
	{
		[dims addObject: [self.file propertiesOfDimensionWithID: i]];
	}
	
	// Read all the variable properties
	NSMutableArray *vars = [NSMutableArray array];
	for (int i=0; i<numVariables; i++)
	{
		[vars addObject: [self.file propertiesOfVariableWithID: i]];
	}
	
	// Create the dimensions
	for ( NSDictionary *dimProperties in dims )
	{
		GLDimension *aDimension = [self createDimensionFromProperties: dimProperties variables: vars];
		
		if (aDimension)
		{
			[self.dimensionDimensionIDMapTable setObject: [dimProperties objectForKey: GLDimensionIDKey] forKey: aDimension];
			[_dimensions addObject:aDimension];
			
			if (aDimension.isMutable) {
				[self startObservingDimension: (GLMutableDimension *)aDimension];
				[_mutableDimensions addObject: aDimension];
			} else {
				[_staticDimensions addObject:aDimension];
			}
		}
	}
	
	// Finally, create the variables
	for ( NSDictionary *varProperties in vars )
	{
		GLFunction *aVariable = [self createVariableFromProperties: varProperties];
		
		if (aVariable) {
			[_variables addObject:aVariable];
			if (aVariable.isMutable) {
				[_mutableVariables addObject: aVariable];
			} else {
				[_staticVariables addObject: aVariable];
			}
		}
	}
}

// Does the best it can to find the associated coordinate variable, but if not, it will still return a valid dimension.
// If it *does* find a coordinate variable, it removes it from vars array.
- (GLDimension *) createDimensionFromProperties: (NSDictionary *) dimProperties variables: (NSMutableArray *) vars
{
	int dimPoints = [[dimProperties objectForKey: GLDimensionLengthKey] intValue];
	NSString *dimName = [dimProperties objectForKey: GLDimensionNameKey];
	int dimID = [[dimProperties objectForKey: GLDimensionIDKey] intValue];
	
	// Can we find an associated coordinate variable?
//	NSDictionary *bestGuess = nil;
	NSDictionary *coordVar = nil;
	for ( NSDictionary *varProperties in vars )
	{
		NSArray *varArray = [varProperties objectForKey: GLVariableDimensionsArrayKey];
		
		// This variable had better only have one dimension
		if ( varArray.count != 1) continue;
		
		// The one dimension had better be THIS dimension
		if ( [[varArray lastObject] intValue] != dimID ) continue;
		
		// We'll use this as the coordinate variable if we can't find one where the names match.
//		bestGuess = varProperties;
		
		// And it must have the same name as the dimension
		if ( [[varProperties objectForKey: GLVariableNameKey] isEqualToString: dimName] ) {
			coordVar = varProperties;
			break;
		}
	}
//	if (!coordVar && bestGuess) {
//		coordVar = bestGuess;
//	}
	
	// Is the dimension unlimited?
	BOOL isUnlimited = NO;
	for ( NSNumber *unlimited in [self.globalAttributes objectForKey: GLUnlimitedDimensionsArrayKey] )
	{
		if (unlimited.intValue == [[dimProperties objectForKey: GLDimensionIDKey] intValue]) {
			isUnlimited = YES;
			break;
		}
	}
	
	// Now go create the appropriate versions of the dimensions
	Class GLDimensionClass = isUnlimited ? [GLMutableDimension class] : [GLDimension class];
	GLDimension *dimension;
	if (coordVar)
	{
		// Read the values in from the coordinate variable
		int varID = [[coordVar objectForKey: GLVariableIDKey] intValue];
		NSArray *indexRange = [NSArray arrayWithObject: [NSValue valueWithRange: NSMakeRange(0, dimPoints)]];
		NSMutableData *buffer = [[NSMutableData alloc] initWithLength: dimPoints*sizeof(GLFloat)];
		if ( sizeof(GLFloat) == sizeof(double) ) {
			[self.file readDoubleVariableWithID: varID intoData: buffer indexRange: indexRange];
		} else {
			[self.file readFloatVariableWithID: varID intoData: buffer indexRange: indexRange];
		}
		
		// Dump the coordinate variable from future consideration
		[vars removeObject: coordVar];
		
		dimension = [[GLDimensionClass alloc] initWithNPoints:dimPoints values:buffer];
		dimension.name = dimName;
		[self.dimensionVariableIDMapTable setObject: [NSNumber numberWithInteger: varID] forKey: dimension];
	} else {
        dimension = [[GLDimensionClass alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: dimPoints domainMin: 0.0 length: 1.0*(dimPoints-1)];
	}
	
	[self.dimensionDimensionIDMapTable setObject: [NSNumber numberWithInt: dimID] forKey: dimension];
	
	return dimension;
}

// This is called after all the dimensions have been created.
- (GLFunction *) createVariableFromProperties: (NSDictionary *) variableProperties
{
	
	// Create an array of the dimension
	NSMutableArray *dimensions = [NSMutableArray array];
	NSArray *dimensionIDs = [variableProperties objectForKey: GLVariableDimensionsArrayKey];
	
	for ( NSNumber *dimensionID in dimensionIDs) {
		for ( GLDimension *dimension in self.dimensionDimensionIDMapTable ) {
			if ([dimensionID isEqualToNumber: [self.dimensionDimensionIDMapTable objectForKey: dimension]]) {
				[dimensions addObject: dimension];
			}
		}
	}
	
	GLNetCDFVariable *netcdfVariable = [GLNetCDFVariable functionOfType: kGLRealDataFormat withDimensions: dimensions forEquation: self.equation];
	netcdfVariable.file = self.file;
	netcdfVariable.name = [variableProperties objectForKey: GLVariableNameKey];
	[netcdfVariable.metadata addEntriesFromDictionary: [variableProperties objectForKey: GLVariableAttributesKey]];
	netcdfVariable.variableID = [[variableProperties objectForKey: GLVariableIDKey] intValue];
	[netcdfVariable setupDependency];
	
	return netcdfVariable;
}

const static NSString *MutableDimensionContext = @"com.EarlyInnovations.MutableDimensionContext";
- (void) startObservingDimension: (GLMutableDimension *) mutableDimension
{
	[mutableDimension addObserver: self forKeyPath: @"nPoints" options: (NSKeyValueObservingOptionNew|NSKeyValueObservingOptionOld) context: &MutableDimensionContext];
}

- (void) stopObservingDimension: (GLMutableDimension *) mutableDimension
{
	[mutableDimension removeObserver:self forKeyPath:@"nPoints" context: &MutableDimensionContext];
}

- (void) observeValueForKeyPath:(NSString *)keyPath ofObject:(id)object change:(NSDictionary *)change context:(void *)context
{
	if (context == &MutableDimensionContext)
	{
		GLMutableDimension *mutableDimension = object;
		NSNumber *oldNPoints = [change objectForKey: NSKeyValueChangeOldKey];
		NSNumber *newNPoints = [change objectForKey: NSKeyValueChangeNewKey];
		NSNumber *variableIDNumber = [self.dimensionVariableIDMapTable objectForKey: mutableDimension];

		if (oldNPoints && newNPoints && variableIDNumber) {
			NSData *theData = mutableDimension.data;
			NSRange byteRange = NSMakeRange(oldNPoints.intValue * sizeof(GLFloat), (newNPoints.intValue-oldNPoints.intValue)*sizeof(GLFloat));
			NSRange indexRange = NSMakeRange(byteRange.location/sizeof(GLFloat), byteRange.length/sizeof(GLFloat));
			NSArray *array = [NSArray arrayWithObject: [NSValue valueWithRange: indexRange]];
			if (NSMaxRange(byteRange) <= theData.length) {
				NSData *subdata = [theData subdataWithRange: byteRange];
				if (sizeof(GLFloat)==sizeof(float)) {
					[self.file writeFloatData: subdata toVariableWithID: [variableIDNumber intValue] atIndexRange: array];
				} else {
					[self.file writeDoubleData: subdata toVariableWithID: [variableIDNumber intValue] atIndexRange: array];
				}
			}
		} else {
			NSLog(@"Ack! GLNetCDFFile.m unable to retrieve old and new points of mutable dimension");
		}
	}
}

/************************************************/
/*		Dimensions, Variables & Attributes		*/
/************************************************/

#pragma mark -
#pragma mark Dimensions, Variables & Attributes
#pragma mark

@synthesize dimensions=_dimensions;
@synthesize staticDimensions=_staticDimensions;
@synthesize mutableDimensions=_mutableDimensions;

@synthesize variables=_variables;
@synthesize staticVariables=_staticVariables;
@synthesize mutableVariables=_mutableVariables;

- (NSArray *) variables {
	NSMutableArray *array = [NSMutableArray arrayWithCapacity: _variables.count];
	for (GLNetCDFVariable *var in _variables) {
		[array addObject: [var copy]];
	}
	return array;
}

- (NSArray *) staticVariables {
	NSMutableArray *array = [NSMutableArray arrayWithCapacity: _staticVariables.count];
	for (GLNetCDFVariable *var in _staticVariables) {
		[array addObject: [var copy]];
	}
	return array;
}

- (NSArray *) mutableVariables {
	NSMutableArray *array = [NSMutableArray arrayWithCapacity: _mutableVariables.count];
	for (GLNetCDFVariable *var in _mutableVariables) {
		[array addObject: [var copy]];
	}
	return array;
}

- (GLDimension *) dimensionWithName: (NSString *) name {
	for (GLDimension *dim in self.dimensions) {
		if ( [dim.name isEqualToString: name] ) return dim;
	}
	return nil;
}

- (GLNetCDFVariable *) variableWithName: (NSString *) name {
	for (GLNetCDFVariable *var in _variables) {
		if ( [var.name isEqualToString: name] ) return [var copy];
	}
	return nil;
}

@synthesize globalAttributes=_globalAttributes;
- (void) setGlobalAttribute: (id) anObject forKey: (id) aKey {
	[self.file addGlobalAttributes: [NSDictionary dictionaryWithObject: anObject forKey: aKey]];
	[_globalAttributes setObject: anObject forKey: aKey];
}

@synthesize dimensionDimensionIDMapTable;
@synthesize dimensionVariableIDMapTable;

- (void) addDimension:(GLDimension *)dimension
{
	if ([self.dimensionDimensionIDMapTable objectForKey: dimension]) {
		return;
	}
	
	if (!dimension.name) {
		NSLog(@"GLNetCDFFile.m: attempting to add a dimension without a name!");
		return;
	}
		
	int dimensionID = [self.file addDimensionWithName: dimension.name length: dimension.isMutable ? NC_UNLIMITED : (int) dimension.nPoints];
	[self.dimensionDimensionIDMapTable setObject: [NSNumber numberWithInt: dimensionID] forKey: dimension];
	
	nc_type type = (sizeof(GLFloat) == 4) ? NC_FLOAT : NC_DOUBLE;
	NSMutableDictionary *properties = [NSMutableDictionary dictionary];
	[properties setObject: [NSNumber numberWithInteger: GLCurrentNetCDFSchemaVersion] forKey: GLNetCDFSchemaVersionKey];
	[properties setObject: [NSNumber numberWithInt: dimensionID] forKey: GLNetCDFSchemaDimensionIDKey];
	[properties setObject: [NSNumber numberWithBool: YES] forKey: GLNetCDFSchemaIsCoordinateVariableKey];
	[properties setObject: [NSNumber numberWithBool: dimension.isMutable] forKey: GLNetCDFSchemaMutableKey];
    [properties setObject: @(dimension.domainLength) forKey: GLNetCDFSchemaDomainLengthKey];
	if (dimension.isEvenlySampled) {
		[properties setObject: [NSNumber numberWithDouble: dimension.sampleInterval] forKey: GLNetCDFSchemaSampleIntervalKey];
		[properties setObject: [NSNumber numberWithDouble: dimension.domainMin] forKey: GLNetCDFSchemaDomainMinimumKey];
	}
	[properties setObject: [NSNumber numberWithBool: dimension.isFrequencyDomain] forKey: GLNetCDFSchemaIsFrequencyDomainKey];
	[properties setObject: dimension.units ?: @"" forKey: GLNetCDFSchemaUnitsKey];
	[properties setObject: @(dimension.basisFunction) forKey: GLNetCDFSchemaBasisFunctionKey];
	[properties setObject: @(dimension.gridType) forKey: GLNetCDFSchemaGridTypeKey];
	
	NSInteger variableID = [self.file addVariableOfType: type withName: dimension.name dimensions: [NSArray arrayWithObject: [NSNumber numberWithInt: dimensionID]] attributes: properties];
	[self.dimensionVariableIDMapTable setObject: [NSNumber numberWithInteger: variableID] forKey: dimension];
	NSArray *array = [NSArray arrayWithObject: [NSValue valueWithRange: NSMakeRange(0, dimension.nPoints)]];
	if (sizeof(GLFloat)==sizeof(float)) {
		[self.file writeFloatData: dimension.data toVariableWithID: (int) variableID atIndexRange: array];
	} else {
		[self.file writeDoubleData: dimension.data toVariableWithID: variableID atIndexRange: array];
	}
	
	[_dimensions addObject: dimension];
	if (dimension.isMutable) {
		[self startObservingDimension: (GLMutableDimension *)dimension];
		[_mutableDimensions addObject: dimension];
	} else {
		[_staticDimensions addObject: dimension];
	}
}

- (GLNetCDFVariable *) addVariable:(GLFunction *)variable
{
	if (!variable.name) {
		NSLog(@"GLNetCDFFile.m: attempting to add a variable without a name!");
		return nil;
	} else if (!variable.isMutable && !variable.hasData) {
		[variable solve];
//		NSLog(@"GLNetCDFFile.m: attempting to add a static variable without any data!");
//		return nil;
	}
		
	for (GLDimension *dimension in variable.dimensions) {
		[self addDimension: dimension];
	}
	
	// Create an array of the dimension IDs
	NSMutableArray *dimensionIDs = [NSMutableArray array];
	for (GLDimension *dimension in variable.dimensions) {
		NSNumber *dimensionID = [self.dimensionDimensionIDMapTable objectForKey: dimension];
		if (dimensionID) {
			[dimensionIDs addObject: dimensionID];
		} else {
			NSLog(@"GLNetCDFFile.m: Attempting to add variable. Unable to get dimension ID for %@!", [dimension description]);
			return nil;
		}
	}
	
	nc_type type = (sizeof(GLFloat) == 4) ? NC_FLOAT : NC_DOUBLE;
	
	NSMutableDictionary *properties = [NSMutableDictionary dictionary];
	[properties addEntriesFromDictionary: variable.metadata];
	[properties setObject: [NSNumber numberWithInteger: GLCurrentNetCDFSchemaVersion] forKey: GLNetCDFSchemaVersionKey];
	[properties setObject: [NSNumber numberWithUnsignedInteger: variable.uniqueID] forKey: GLNetCDFSchemaUniqueVariableIDKey];
	[properties setObject: [NSNumber numberWithBool: NO] forKey: GLNetCDFSchemaIsCoordinateVariableKey];
	[properties setObject: [NSNumber numberWithBool: variable.isMutable] forKey: GLNetCDFSchemaMutableKey];
	[properties setObject: [NSNumber numberWithBool: variable.isFrequencyDomain] forKey: GLNetCDFSchemaIsFrequencyDomainKey];
	[properties setObject: [NSNumber numberWithBool: variable.isComplex] forKey: GLNetCDFSchemaIsComplexKey];
	[properties setObject: variable.units ?: @"" forKey: GLNetCDFSchemaUnitsKey];
	
	NSString *name = variable.isComplex ? [NSString stringWithFormat: @"%@_realp", variable.name] : variable.name;
	if (variable.isComplex) {
		[properties setObject: [NSNumber numberWithBool: YES] forKey: GLNetCDFSchemaIsRealPartKey];
		[properties setObject: [NSNumber numberWithBool: NO] forKey: GLNetCDFSchemaIsImaginaryPartKey];
		[properties setObject: variable.name forKey: GLNetCDFSchemaProperNameKey];
	}
	int variableID = [self.file addVariableOfType: type withName: name dimensions: dimensionIDs attributes: properties];
	
	// First we add the real part to the file.
	[self.equation solveForVariable: variable];
	if (variable.hasData) {
		NSData *data;
		if (variable.isComplex) {
			GLSplitComplex split = variable.splitComplex;
			data = [NSData dataWithBytes: split.realp length: (variable.nDataPoints)*(sizeof(GLFloat))];
		} else {
			data = variable.data;
		}
		if (type == NC_FLOAT) {
			[self.file setFloatData: data forVariableWithID: variableID];
		} else {
			[self.file setDoubleData: data forVariableWithID: variableID];
		}
	}
	
	// Now we deal with the complex part.
	int imagpVariableID = -1;
	if (variable.isComplex && variable.hasData) {
		name = [NSString stringWithFormat: @"%@_imagp", variable.name];
		[properties setObject: [NSNumber numberWithBool: NO] forKey: GLNetCDFSchemaIsRealPartKey];
		[properties setObject: [NSNumber numberWithBool: YES] forKey: GLNetCDFSchemaIsImaginaryPartKey];
		imagpVariableID = [self.file addVariableOfType: type withName: name dimensions: dimensionIDs attributes: properties];
		
		NSData *data = [NSData dataWithBytes: variable.splitComplex.imagp length: (variable.nDataPoints)*(sizeof(GLFloat))];
		if (type == NC_FLOAT) {
			[self.file setFloatData: data forVariableWithID: imagpVariableID];
		} else {
			[self.file setDoubleData: data forVariableWithID: imagpVariableID];
		}
	}
	
	GLNetCDFVariable *netcdfVariable = [GLNetCDFVariable variableWithVariable: variable];
	netcdfVariable.file = self.file;
	netcdfVariable.variableID = variableID;
	netcdfVariable.imagpVariableID = imagpVariableID;
	[netcdfVariable setupDependency];
	
	[_variables addObject: netcdfVariable];
	if (netcdfVariable.isMutable) {
		[_mutableVariables addObject: netcdfVariable];
	} else {
		[_staticVariables addObject: netcdfVariable];
	}
	
	return netcdfVariable;
}


/************************************************/
/*		Internal Representation					*/
/************************************************/

#pragma mark -
#pragma mark Internal Representation
#pragma mark

@end
