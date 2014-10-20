//
//  GLDimension.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 4/21/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import "GLDimension.h"

@interface GLDimension ()
{
    // These are made public so that the subclass NSMutableDimension can use them directly.
    @public
    GLFloat _sampleInterval;
    GLFloat _domainMin;
    GLFloat _domainLength;
    NSUInteger _nPoints;
    NSUInteger _dataBytes;
    
    @private
    NSMapTable *_scaledDimensionMapTable;
}

- (void) populateEvenlySampledValues;

// If I don't set this as strong-strong, we could easily lose the spatialDomainDimension.
// If I do set this as strong-strong, we have a retain cycle. What to do!
@property(readwrite, strong) GLDimension *frequencyDomainDimension;
@property(readwrite, weak) GLDimension *spatialDomainDimension;
@property(readwrite, nonatomic) BOOL isMutable;
@property(readwrite, strong, nonatomic) NSMapTable *rangeSubdimensionMap;

@end

static NSMapTable *dimensionMapTableMap = nil;
static NSMapTable *transformSpatialDimensionMap = nil;

static NSString *GLDimensionNameKey = @"GLDimensionNameKey";
static NSString *GLDimensionUnitsKey = @"GLDimensionUnitsKey";
static NSString *GLDimensionNPointsKey = @"GLDimensionNPointsKey";
static NSString *GLDimensionDomainMinKey = @"GLDimensionDomainMinKey";
static NSString *GLDimensionDomainLengthKey = @"GLDimensionDomainLengthKey";
static NSString *GLDimensionGridTypeKey = @"GLDimensionGridTypeKey";
static NSString *GLDimensionIsEvenlySampledKey = @"GLDimensionIsEvenlySampledKey";
static NSString *GLDimensionSampleIntervalKey = @"GLDimensionSampleIntervalKey";
static NSString *GLDimensionIsMutableKey = @"GLDimensionIsMutableKey";
static NSString *GLDimensionBasisFunctionKey = @"GLDimensionBasisFunctionKey";
static NSString *GLDimensionDifferentiationBasisKey = @"GLDimensionDifferentiationBasisKey";
static NSString *GLDimensionIsStrictlyPositiveKey = @"GLDimensionIsStrictlyPositiveKey";
static NSString *GLDimensionNDataBytesKey = @"GLDimensionNDataBytesKey";
static NSString *GLDimensionDataKey = @"GLDimensionDataKey";
static NSString *GLDimensionMapTableKey = @"GLDimensionMapTableKey";
static NSString *GLDimensionSpatialDimensionKey = @"GLDimensionSpatialDimensionKey";

@implementation GLDimension


- (void)encodeWithCoder:(NSCoder *)coder
{
    if (self.name) [coder encodeObject:self.name forKey:GLDimensionNameKey];
    if (self.units) [coder encodeObject:self.units forKey:GLDimensionUnitsKey];
    
    [coder encodeObject: @(self.nPoints) forKey:GLDimensionNPointsKey];
    [coder encodeObject: @(self.domainMin) forKey:GLDimensionDomainMinKey];
    [coder encodeObject: @(self.domainLength) forKey:GLDimensionDomainLengthKey];
    [coder encodeObject: @(self.gridType) forKey:GLDimensionGridTypeKey];
    [coder encodeObject: @(self.isEvenlySampled) forKey:GLDimensionIsEvenlySampledKey];
    [coder encodeObject: @(self.sampleInterval) forKey:GLDimensionSampleIntervalKey];
    [coder encodeObject: @(self.isMutable) forKey:GLDimensionIsMutableKey];
    [coder encodeObject: @(self.basisFunction) forKey:GLDimensionBasisFunctionKey];
    [coder encodeObject: @(self.differentiationBasis) forKey:GLDimensionDifferentiationBasisKey];
    [coder encodeObject: @(self.isStrictlyPositive) forKey:GLDimensionIsStrictlyPositiveKey];
    
    [coder encodeObject: @(self.dataBytes) forKey:GLDimensionNDataBytesKey];
    [coder encodeObject: self.data forKey:GLDimensionDataKey];
    
    if (self.basisFunction == kGLDeltaBasis) {
        NSMapTable *dimensionMapTable = [GLDimension transformMapForDimension: self];
        if (dimensionMapTable) [coder encodeObject: dimensionMapTable forKey:GLDimensionMapTableKey];
    } else {
        GLDimension *spatialDimension = [GLDimension spatialDimensionForTransformDimension: self];
        if (spatialDimension) [coder encodeObject: spatialDimension forKey:GLDimensionSpatialDimensionKey];
    }
}

- (id)initWithCoder:(NSCoder *)decoder
{
    if ((self=[super init])) {
        _name = [decoder decodeObjectForKey: GLDimensionNameKey];
        _units = [decoder decodeObjectForKey: GLDimensionUnitsKey];
        
        _nPoints = [[decoder decodeObjectForKey: GLDimensionNPointsKey] unsignedIntegerValue];
        _domainMin = [[decoder decodeObjectForKey: GLDimensionDomainMinKey] doubleValue];
        _domainLength = [[decoder decodeObjectForKey: GLDimensionDomainLengthKey] doubleValue];
        _gridType = [[decoder decodeObjectForKey: GLDimensionGridTypeKey] unsignedIntegerValue];
        _isEvenlySampled = [[decoder decodeObjectForKey: GLDimensionIsEvenlySampledKey] boolValue];
        _sampleInterval = [[decoder decodeObjectForKey: GLDimensionSampleIntervalKey] doubleValue];
        _isMutable = [[decoder decodeObjectForKey: GLDimensionIsMutableKey] boolValue];
        _basisFunction = [[decoder decodeObjectForKey: GLDimensionBasisFunctionKey] unsignedIntegerValue];
        _differentiationBasis = [[decoder decodeObjectForKey: GLDimensionDifferentiationBasisKey] unsignedIntegerValue];
        _isStrictlyPositive = [[decoder decodeObjectForKey: GLDimensionIsStrictlyPositiveKey] boolValue];
        
        _dataBytes = [[decoder decodeObjectForKey: GLDimensionNDataBytesKey] unsignedIntegerValue];
        _data = [decoder decodeObjectForKey: GLDimensionDataKey];
        
        NSMapTable *dimensionMapTable = [decoder decodeObjectForKey: GLDimensionMapTableKey];
        if (dimensionMapTable) {
            NSMapTable *existingDimensionMapTable = [GLDimension transformMapForDimension: self];
            for (id key in dimensionMapTable) {
                [existingDimensionMapTable setObject: [dimensionMapTable objectForKey: key] forKey: key];
            }
        }
        
        GLDimension *spatialDimension = [decoder decodeObjectForKey: GLDimensionSpatialDimensionKey];
        if (spatialDimension) {
            [GLDimension setSpatialDimension:spatialDimension forTransformDimension:self];
        }
    }
    return self;
}

/************************************************/
/*		Related Dimensions						*/
/************************************************/

#pragma mark -
#pragma mark Related Dimensions
#pragma mark

// These methods are used to store dimensions that are related to each other, either as transformed dimensions or subdimensions.

// Returns a map table for the requested spatial dimension which maps to transformed dimensions
+ (NSMapTable *) transformMapForDimension: (GLDimension *) aDimension
{
	if (!dimensionMapTableMap) {
		dimensionMapTableMap = [NSMapTable mapTableWithKeyOptions: NSMapTableObjectPointerPersonality valueOptions:NSMapTableObjectPointerPersonality];
	}
	
	NSMapTable *dimensionMapTable = [dimensionMapTableMap objectForKey: aDimension];
	if (!dimensionMapTable) {
		dimensionMapTable = [NSMapTable mapTableWithKeyOptions: NSMapTableObjectPointerPersonality valueOptions:NSMapTableObjectPointerPersonality];
		[dimensionMapTableMap setObject: dimensionMapTable forKey: aDimension];
	}
	return dimensionMapTable;
}

+ (GLDimension *) transformOfDimension: (GLDimension *) existingDimension transformedToBasis: (GLBasisFunction) basis strictlyPositive: (BOOL) positive
{
	NSMapTable *transformMap = [self transformMapForDimension: existingDimension];
	
	NSUInteger key = (basis << 4);
	key |= positive;
	
	return [transformMap objectForKey: [NSNumber numberWithUnsignedInteger: key]];
}

+ (void) setTransform: (GLDimension *) tDim ofDimension: (GLDimension *) existingDimension transformedToBasis: (GLBasisFunction) basis strictlyPositive: (BOOL) positive
{
	NSMapTable *transformMap = [self transformMapForDimension: existingDimension];
	
	NSUInteger key = (basis << 4);
	key |= positive;
	
	[transformMap setObject: tDim forKey:[NSNumber numberWithUnsignedInteger: key]];
}



+ (void) setSpatialDimension: (GLDimension *) sDim forTransformDimension: (GLDimension *) tDim
{
	if (!transformSpatialDimensionMap) {
		transformSpatialDimensionMap = [NSMapTable mapTableWithKeyOptions: NSMapTableObjectPointerPersonality valueOptions:NSMapTableObjectPointerPersonality];
	}
    [transformSpatialDimensionMap setObject: sDim forKey: tDim];
}

+ (GLDimension *) spatialDimensionForTransformDimension: (GLDimension *) tDim
{
	if (!transformSpatialDimensionMap) {
		transformSpatialDimensionMap = [NSMapTable mapTableWithKeyOptions: NSMapTableObjectPointerPersonality valueOptions:NSMapTableObjectPointerPersonality];
	}
    return [transformSpatialDimensionMap objectForKey: tDim];
}

/************************************************/
/*		Convenience Methods						*/
/************************************************/

#pragma mark -
#pragma mark Convenience Methods
#pragma mark

+ (NSArray *) dimensionsForRealFunctionWithDimensions: (NSArray *) dimensions transformedToBasis: (NSArray *) basis
{
	if ( basis.count == 1 && dimensions.count > 1) {
        NSMutableArray *array = [NSMutableArray array];
        for (NSUInteger i=0; i < dimensions.count; i++) {
            [array addObject: basis.lastObject];
        }
        basis = array;
    }
	
	
	NSMutableArray *transformedDimensions = [NSMutableArray array];
	NSUInteger lastExponentialIndex = NSNotFound;
	for (NSUInteger i=0; i < dimensions.count; i++) {
		if ([basis[i] unsignedIntegerValue] == kGLExponentialBasis) {
			lastExponentialIndex = i;
		}
	}
    
	for (NSUInteger i=0; i < dimensions.count; i++) {
		BOOL strictlyPositive = ( lastExponentialIndex == i || [basis[i] unsignedIntegerValue] > 1) ? YES : NO;
		GLDimension *dim = [[GLDimension alloc] initAsDimension: dimensions[i] transformedToBasis: [basis[i] unsignedIntegerValue] strictlyPositive: strictlyPositive];
		[transformedDimensions addObject: dim];
	}
	
	return transformedDimensions;
}

+ (NSArray *) dimensionsForComplexFunctionWithDimensions: (NSArray *) dimensions transformedToBasis: (NSArray *) basis
{
	if ( basis.count == 1 && dimensions.count > 1) {
        NSMutableArray *array = [NSMutableArray array];
        for (NSUInteger i=0; i < dimensions.count; i++) {
            [array addObject: basis.lastObject];
        }
        basis = array;
    }
	
	
	NSMutableArray *transformedDimensions = [NSMutableArray array];
	NSUInteger lastExponentialIndex = NSNotFound;
	for (NSUInteger i=0; i < dimensions.count; i++) {
		BOOL strictlyPositive = ( lastExponentialIndex == i || [basis[i] unsignedIntegerValue] > 1) ? YES : NO;
		GLDimension *dim = [[GLDimension alloc] initAsDimension: dimensions[i] transformedToBasis: [basis[i] unsignedIntegerValue] strictlyPositive: strictlyPositive];
		[transformedDimensions addObject: dim];
	}
	
	return transformedDimensions;
}

+ (GLDimension *) dimensionXWithNPoints: (NSUInteger) numPoints length: (GLFloat) theLength;{
	GLDimension *aDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:numPoints domainMin:0 length:theLength];
	aDim.name = @"x";
	aDim.units = @"unitless";	
	aDim.basisFunction = kGLDeltaBasis;
	return aDim;
}
+ (GLDimension *) dimensionYWithNPoints: (NSUInteger) numPoints length: (GLFloat) theLength;{
	GLDimension *aDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:numPoints domainMin:0 length:theLength];
	aDim.name = @"y";
	aDim.units = @"unitless";
	aDim.basisFunction = kGLDeltaBasis;
	return aDim;
}
+ (GLDimension *) dimensionZWithNPoints: (NSUInteger) numPoints length: (GLFloat) theLength;{
	GLDimension *aDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:numPoints domainMin:0 length:theLength];
	aDim.name = @"z";
	aDim.units = @"unitless";
	aDim.basisFunction = kGLDeltaBasis;
	return aDim;
}
+ (GLMutableDimension *) dimensionT {
	GLMutableDimension *aDim = [[GLMutableDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints:1 domainMin:0 length:0];
	aDim.name = @"t";
	aDim.units = @"unitless";
	aDim.basisFunction = kGLDeltaBasis;
	return aDim;
}

/************************************************/
/*		Initialization							*/
/************************************************/

#pragma mark -
#pragma mark Initialization
#pragma mark

- (GLDimension *) initDimensionWithGrid: (GLGridType) gridType nPoints: (NSUInteger) numPoints domainMin: (GLFloat) theMin length: (GLFloat) theLength
{
	if ((self = [super init]))
	{
		_gridType = gridType;
		_nPoints = numPoints;
		_domainMin = theMin;
		_domainLength = theLength;
		_basisFunction = kGLDeltaBasis;
		_isMutable = NO;
        _isStrictlyPositive = ( _domainMin >= 0 );

		_dataBytes = _nPoints*sizeof(GLFloat);
		_data =[NSMutableData dataWithLength: _dataBytes];
		
		GLFloat *f = self.data.mutableBytes;
		
		if (gridType == kGLEndpointGrid)
		{
			_isEvenlySampled = YES;
			_sampleInterval = self.nPoints > 1 ? _domainLength / ( (double) (self.nPoints-1)) : 0;
			for (NSUInteger i=0; i<self.nPoints; i++) {
				f[i] = _sampleInterval * ( (GLFloat) i) + _domainMin;
			}
			_differentiationBasis = kGLDeltaBasis;
		}
		else if (gridType == kGLInteriorGrid)
		{
			_isEvenlySampled = YES;
			_sampleInterval = self.nPoints > 1 ? _domainLength / ( (double) (self.nPoints)) : 0;
			for (NSUInteger i=0; i<self.nPoints; i++) {
				f[i] = _sampleInterval * ( (GLFloat) i + 0.5) + _domainMin;
			}
			_differentiationBasis = kGLCosineBasis;
		}
		else if (gridType == kGLPeriodicGrid)
		{
			_isEvenlySampled = YES;
			_sampleInterval = self.nPoints > 1 ? _domainLength / ( (double) (self.nPoints)) : 0;
			for (NSUInteger i=0; i<self.nPoints; i++) {
				f[i] = _sampleInterval * ( (GLFloat) i ) + _domainMin;
			}
			_differentiationBasis = kGLExponentialBasis;
		}
		else if (gridType == kGLChebyshevEndpointGrid)
		{
			_isEvenlySampled = NO;
			_sampleInterval = self.nPoints > 1 ? M_PI / ( (double) (self.nPoints-1)) : 0;
			for (NSUInteger i=0; i<self.nPoints; i++) {
				f[i] = 0.5 * _domainLength * (cos( _sampleInterval * ( (GLFloat) i ) ) + 1.0) + _domainMin;
			}
			_differentiationBasis = kGLChebyshevBasis;
		}
		else if (gridType == kGLChebyshevInteriorGrid)
		{
			_isEvenlySampled = NO;
			_sampleInterval = self.nPoints > 1 ? M_PI / ( (double) (self.nPoints)) : 0;
			for (NSUInteger i=0; i<self.nPoints; i++) {
				f[i] = 0.5 * _domainLength * (cos( _sampleInterval * ( (GLFloat) i + 0.5 ) ) + 1.0) + _domainMin;
			}
			_differentiationBasis = kGLChebyshevBasis;
		}
	}
	
	return self;
}

- (GLDimension *) initWithNPoints: (NSUInteger) numPoints values: (NSData *) values
{
	if (numPoints < 1) {
        [NSException raise: @"BadFormatException" format:@"Cannot initialize an unevenly sampled dimension with less than 1 point."];
	}
	
	if ((self = [super init]))
	{
		GLFloat *fFrom = (GLFloat *) values.bytes;
		
		_isEvenlySampled = NO;
		_gridType = kGLEndpointGrid;
		_nPoints = numPoints;
		_domainMin = fFrom[0];
		_domainLength = fFrom[numPoints-1] - fFrom[0];
		_basisFunction = kGLDeltaBasis;
		_isMutable = NO;
		
		_dataBytes = self.nPoints*sizeof(GLFloat);
		_data =[NSMutableData dataWithLength: _dataBytes];
		
		GLFloat *f = self.data.mutableBytes;
		for (NSUInteger i=0; i < self.nPoints; i++) {
			f[i] = fFrom[i];
		}
	}
	
	return self;
}

- (GLDimension *) initWithPoints: (NSArray *) pointsArray
{
	if (pointsArray.count < 1) {
        [NSException raise: @"BadFormatException" format:@"Cannot initialize an unevenly sampled dimension with less than 1 point."];
	}
	
	if ((self = [super init]))
	{
		_isEvenlySampled = NO;
		_gridType = kGLEndpointGrid;
		_nPoints = pointsArray.count;
		_domainMin = [[pointsArray objectAtIndex:0] doubleValue];
		_domainLength = [[pointsArray lastObject] doubleValue] - [[pointsArray objectAtIndex:0] doubleValue];
		_basisFunction = kGLDeltaBasis;
		_isMutable = NO;
		
		_dataBytes = self.nPoints*sizeof(GLFloat);
		_data =[NSMutableData dataWithLength: _dataBytes];
		
		GLFloat *f = self.data.mutableBytes;
		NSUInteger i=0;
		for ( NSNumber *point in pointsArray) {
			f[i]=[point doubleValue];
			i++;
		}
	}
	
	return self;
}

// The logic here sucks. It's unnecessarily complicated.
- (GLDimension *) initAsDimension: (GLDimension *) existingDimension transformedToBasis: (GLBasisFunction) basis strictlyPositive: (BOOL) isPositive
{
	if (existingDimension.basisFunction == basis) {
		return existingDimension;
	}
	
	// Identify the spatial dimension
	GLDimension *spatialDimension;
	if (existingDimension.basisFunction == kGLDeltaBasis) {
		spatialDimension = existingDimension;
	} else {
		spatialDimension = [GLDimension spatialDimensionForTransformDimension: existingDimension];
        if (!spatialDimension) {
            [NSException raise: @"DimensionException" format: @"Current assumptions suggest this should already exist."];
        }
	}
	
	// If that's the dimension requested, return it.
	if (basis == kGLDeltaBasis) {
		return spatialDimension;
	}
	
	// See if we've already done this match.
	GLDimension *alreadyCreated = [GLDimension transformOfDimension: existingDimension transformedToBasis: basis strictlyPositive: isPositive];
	if (alreadyCreated) {
		return alreadyCreated;
	}
	
	// If we're switching between the DCT-I and DST-I basis, we need a *different* spatial dimension
	if ( existingDimension.basisFunction == kGLDiscreteCosineTransformIBasis && basis == kGLDiscreteSineTransformIBasis ) {
        [NSException raise: @"DeprecationException" format:@"This functionality might be deprecated."];
//		GLFloat aMin = spatialDimension.domainMin;
//		GLFloat aSampleInterval = spatialDimension.domainLength / ( (GLFloat) spatialDimension.nPoints+1);
//		NSString *aName = spatialDimension.name;
//		NSString *aUnits = spatialDimension.units;
//        ;
//		spatialDimension = [[GLDimension alloc] initPeriodicDimension: NO nPoints: spatialDimension.nPoints domainMin: aMin + aSampleInterval sampleInterval: aSampleInterval];
//		spatialDimension.name = aName;
//		spatialDimension.units = aUnits;
	} else if ( existingDimension.basisFunction == kGLDiscreteSineTransformIBasis && basis == kGLDiscreteCosineTransformIBasis ) {
        [NSException raise: @"DeprecationException" format:@"This functionality might be deprecated."];
//		NSString *aName = spatialDimension.name;
//		NSString *aUnits = spatialDimension.units;
//		spatialDimension = [[GLDimension alloc] initPeriodicDimension: NO nPoints: spatialDimension.nPoints domainMin: spatialDimension.domainMin-spatialDimension.sampleInterval length: spatialDimension.domainLength+2*spatialDimension.sampleInterval];
//		spatialDimension.name = aName;
//		spatialDimension.units = aUnits;
	}
	
	if ((self = [super init]))
	{		
		NSString *inverseName = nil;
		if ( [spatialDimension.name isEqualToString: @"x"] ) {
			inverseName = @"k";
		}
		else if ( [spatialDimension.name isEqualToString: @"y"] ) {
			inverseName = @"l";
		}
		else if ( [spatialDimension.name isEqualToString: @"z"] ) {
			inverseName = @"m";
		}
		else if ( [spatialDimension.name isEqualToString: @"t"] ) {
			inverseName = @"omega";
		} else {
			inverseName = [NSString stringWithFormat: @"%@_bar", self.name];
		}
		
		NSString *inverseUnits = spatialDimension.units ? [NSString stringWithFormat: @"cycles per %@", self.units] : @"cycles";
		
		_isEvenlySampled = YES;
		_gridType = kGLEndpointGrid;
		_isMutable = NO;
		_differentiationBasis = basis;
		
		self.basisFunction = basis;
		self.isStrictlyPositive = isPositive;
        
        _domainMin = 0.0;
        
        // The user may have incorrectly decided whether or not the existingDimension was periodic.
        // These particular choices of basis *demand* choice--so we correct the mistake.
		if (self.basisFunction == kGLExponentialBasis)
		{
			if (!self.isStrictlyPositive)
			{	// 0,..,fc-1,-fc,..,-1
				if (spatialDimension.isPeriodic) {
                    _domainLength = ((GLFloat) spatialDimension.nPoints)/spatialDimension.domainLength;
                } else {
                    _domainLength = ((GLFloat) spatialDimension.nPoints)/(spatialDimension.domainLength+spatialDimension.sampleInterval);
                }
				_nPoints = spatialDimension.nPoints;
				_domainMin = -_domainLength/2;
				_sampleInterval = _domainLength / ( (double) (self.nPoints));
			}
			else
			{	// 0..fc --- negative frequencies ignored because they're found with the conjugate of the function.
				
                if (spatialDimension.isPeriodic) {
                    _domainLength = ((GLFloat) spatialDimension.nPoints/2 + 1)/spatialDimension.domainLength;
                } else {
                    _domainLength = ((GLFloat) spatialDimension.nPoints/2 + 1)/(spatialDimension.domainLength+spatialDimension.sampleInterval);
                }
				_nPoints = spatialDimension.nPoints/2 + 1;
				_sampleInterval = _domainLength / ( (double) (self.nPoints));
			}
		}
		else if (self.basisFunction == kGLDiscreteCosineTransformIBasis )
		{
            // 0..fc -- negative frequencies ignored, as they're found with the even symmetric part of the function
            if (spatialDimension.isPeriodic) {
                _domainLength = ((GLFloat) spatialDimension.nPoints-1)/(2*(spatialDimension.domainLength-spatialDimension.sampleInterval));
            } else {
                _domainLength = ((GLFloat) spatialDimension.nPoints-1)/(2*spatialDimension.domainLength);
            }
            _nPoints = spatialDimension.nPoints;
            _sampleInterval = _domainLength / ( (double) (self.nPoints-1));
		}
        else if (self.basisFunction == kGLCosineBasis)
        {
            // 0..fc -- negative frequencies ignored, as they're found with the 
            // Same as the cosine transform above, but the actual domain length should have been treated as one point longer, as if it were periodic.
            if (spatialDimension.gridType == kGLInteriorGrid) {
                _domainLength = ((GLFloat) spatialDimension.nPoints-1)/(2*(spatialDimension.domainLength));
            } else {
                NSLog(@"Inappropriate grid type!");
                if (spatialDimension.isPeriodic) {
                    _domainLength = ((GLFloat) spatialDimension.nPoints-1)/(2*(spatialDimension.domainLength));
                } else {
                    _domainLength = ((GLFloat) spatialDimension.nPoints-1)/(2*(spatialDimension.domainLength+spatialDimension.sampleInterval));
                }
            }
            
            _nPoints = spatialDimension.nPoints;
            _sampleInterval = _domainLength / ( (double) (self.nPoints-1));
        }
		else if (self.basisFunction == kGLDiscreteSineTransformIBasis)
		{
            // 0..fc -- negative frequencies ignored, as they're found with the 
            if (spatialDimension.isPeriodic) {
                _domainLength = ((GLFloat) spatialDimension.nPoints)/(2*(spatialDimension.domainLength+spatialDimension.sampleInterval));
            } else {
                _domainLength = ((GLFloat) spatialDimension.nPoints)/(2*(spatialDimension.domainLength+2*spatialDimension.sampleInterval));
            }
            _nPoints = spatialDimension.nPoints;
            _sampleInterval = _domainLength / ( (double) (self.nPoints));
			_domainMin = _sampleInterval;
			
		}
        else if (self.basisFunction == kGLSineBasis)
        {
            // 0..fc -- negative frequencies ignored, as they're found with the 
            // Same as the cosine transform above, but the actual domain length should have been treated as one point longer, as if it were periodic.
            if (spatialDimension.gridType == kGLInteriorGrid) {
                _domainLength = ((GLFloat) spatialDimension.nPoints-1)/(2*(spatialDimension.domainLength));
            } else {
                NSLog(@"Inappropriate grid type!");
                if (spatialDimension.isPeriodic) {
                    _domainLength = ((GLFloat) spatialDimension.nPoints-1)/(2*(spatialDimension.domainLength));
                } else {
                    _domainLength = ((GLFloat) spatialDimension.nPoints-1)/(2*(spatialDimension.domainLength+spatialDimension.sampleInterval));
                }
            }
            _nPoints = spatialDimension.nPoints;
            _sampleInterval = _domainLength / ( (double) (self.nPoints-1));
			_domainMin = _sampleInterval;
        }
		else if (self.basisFunction == kGLChebyshevBasis )
		{
            _domainLength = 2.0*((GLFloat) spatialDimension.nPoints-1)/(spatialDimension.domainLength);
            _nPoints = spatialDimension.nPoints;
            _sampleInterval = _domainLength / ( (double) (self.nPoints-1));
		}
		else {
			[NSException raise: @"CaseNotImplemented" format: @"This case has not been implemented"];
		}
		
		_dataBytes = self.nPoints*sizeof(GLFloat);
		_data =[NSMutableData dataWithLength: _dataBytes];
		[self populateEvenlySampledValues];
		
		self.name = inverseName;
		self.units = inverseUnits;
	}
	
	[GLDimension setTransform: self ofDimension:spatialDimension transformedToBasis:basis strictlyPositive: isPositive];
    [GLDimension setSpatialDimension: spatialDimension forTransformDimension:self];
	if ( spatialDimension != existingDimension) {
		[GLDimension setTransform: self ofDimension:existingDimension transformedToBasis:basis strictlyPositive: isPositive];
	}
	
	return self;
}

- (GLDimension *) initAsFourierTransformOfDimension: (GLDimension *) existingDimension
{
	if (existingDimension.isFrequencyDomain) {
		if (existingDimension.spatialDomainDimension) {
			return existingDimension.spatialDomainDimension;
		} else {
			NSLog(@"Error in GLDimension -initAsFourierTransformOfDimension. By assumption the space/time domain dimension should exist");
		}
	} else {
		if (existingDimension.frequencyDomainDimension) {
			return existingDimension.frequencyDomainDimension;
		}
	}
	
	if ((self = [super init]))
	{		
		NSString *inverseName = nil;
		if ( [existingDimension.name isEqualToString: @"x"] ) {
			inverseName = @"k";
		}
		else if ( [existingDimension.name isEqualToString: @"y"] ) {
			inverseName = @"l";
		}
		else if ( [existingDimension.name isEqualToString: @"z"] ) {
			inverseName = @"m";
		}
		else if ( [existingDimension.name isEqualToString: @"t"] ) {
			inverseName = @"omega";
		} else {
			inverseName = [NSString stringWithFormat: @"%@_bar", self.name];
		}
		
		NSString *inverseUnits = existingDimension.units ? [NSString stringWithFormat: @"cycles per %@", self.units] : @"cycles";
		
		_isEvenlySampled = YES;
		_nPoints = existingDimension.nPoints;
		_domainMin = 0.0;
		_domainLength = ((GLFloat) existingDimension.nPoints)/existingDimension.domainLength;
		_isMutable = NO;
		_basisFunction = kGLExponentialBasis;
		_isStrictlyPositive = NO;
		
		_dataBytes = self.nPoints*sizeof(GLFloat);
		_data =[NSMutableData dataWithLength: _dataBytes];
		[self populateEvenlySampledValues];
		
		self.name = inverseName;
		self.units = inverseUnits;
		
		existingDimension.frequencyDomainDimension = self;
		self.spatialDomainDimension = existingDimension;
	}
	
	return self;
	
}

- (GLDimension *) initWithDimension: (GLDimension *) existingDim scaledBy: (GLFloat) scale translatedBy: (GLFloat) delta withUnits: (NSString *) newUnits
{
    if ((self=[super init])) {
        _nPoints = existingDim.nPoints;
        _gridType = existingDim.gridType;
        _domainMin = (existingDim.domainMin)*scale+delta;
        _domainLength = (existingDim.domainLength)*scale;
        _sampleInterval = (existingDim.sampleInterval)*scale;
        
        _isEvenlySampled = existingDim.isEvenlySampled;
        _basisFunction = existingDim.basisFunction;
        _isStrictlyPositive = existingDim.isStrictlyPositive;
        _differentiationBasis = existingDim.differentiationBasis;
        
        _dataBytes = self.nPoints*sizeof(GLFloat);
		_data =[NSMutableData dataWithLength: _dataBytes];
        if (self.isEvenlySampled) {
            [self populateEvenlySampledValues];
        } else {
            GLFloat *f = self.data.mutableBytes;
            GLFloat *g = existingDim.data.mutableBytes;
            for (NSUInteger i=0; i<self.nPoints; i++) {
                f[i] = scale*g[i];
            }
        }
        
        self.name = existingDim.name;
        self.units = existingDim.units;
    }
    return self;
}

/************************************************/
/*		Populating point values					*/
/************************************************/

#pragma mark -
#pragma mark Populating point values
#pragma mark

- (GLFloat) evenlySampledValueAtIndex: (NSUInteger) index
{
	if (self.basisFunction == kGLExponentialBasis && self.isStrictlyPositive == NO) {
		if (index < self.nPoints/2) {
			return ( self.domainLength * ( (GLFloat) index) / ((GLFloat) self.nPoints)); 
		} else {
			return ( self.domainLength * ( ((GLFloat) index - (GLFloat) self.nPoints)) / ((GLFloat) self.nPoints));
		}
	} else {
		return ( self.sampleInterval * ( (GLFloat) index) + self.domainMin);
	}
}

- (void) populateEvenlySampledValues
{
	GLFloat *f = self.data.mutableBytes;
	for (NSUInteger i=0; i<self.nPoints; i++) {
		f[i] = [self evenlySampledValueAtIndex: i];
	}
}

/************************************************/
/*		Metadata								*/
/************************************************/

#pragma mark -
#pragma mark Metadata
#pragma mark

/************************************************/
/*		Essential Properties					*/
/************************************************/

#pragma mark -
#pragma mark Essential Properties
#pragma mark

- (BOOL) isPeriodic {
	return (_gridType == kGLPeriodicGrid);
}

- (BOOL) isFrequencyDomain {
	return self.basisFunction > 0 ? YES : NO;
}

@synthesize domainLength=_domainLength;
@synthesize domainMin=_domainMin;
@synthesize sampleInterval=_sampleInterval;
@synthesize nPoints=_nPoints;
@synthesize dataBytes=_dataBytes;

/************************************************/
/*		Value									*/
/************************************************/

#pragma mark -
#pragma mark Value
#pragma mark

- (GLFloat) valueAtIndex: (NSUInteger) index
{
	GLFloat *f = self.data.mutableBytes;
	if (index < self.nPoints) {
		return f[index];
	} else {
		NSLog(@"GLDimension requested index out of bounds!");
	}
    return 0.0;
}

- (NSArray *) points
{
	NSMutableArray *array = [NSMutableArray array];
	for (NSUInteger i=0; i<self.nPoints; i++) {
		[array addObject: [NSNumber numberWithDouble: [self valueAtIndex: i]]];
	}
	return array;
}

/************************************************/
/*		Derived Dimensions						*/
/************************************************/

#pragma mark -
#pragma mark Derived Dimensions
#pragma mark

- (GLDimension *) fourierTransformedDimension
{
	if (self.isFrequencyDomain) {
        return [[GLDimension alloc] initAsDimension: self transformedToBasis: kGLDeltaBasis strictlyPositive: NO];
    } else {
        return [[GLDimension alloc] initAsDimension: self transformedToBasis: kGLExponentialBasis strictlyPositive: NO];
    }
    
    //return [[GLDimension alloc] initAsFourierTransformOfDimension: self];
}

- (GLDimension *) subdimensionWithRange: (NSRange) range
{
	if (self.rangeSubdimensionMap) {
		for (NSValue *value in self.rangeSubdimensionMap) {
			if ( NSEqualRanges( value.rangeValue, range)) {
				return [self.rangeSubdimensionMap objectForKey: value];
			}
		}
	}
	
	if (range.location > self.nPoints - 1) return nil;
	if ((range.location + range.length) > self.nPoints) return nil;
	if (range.location == 0 && range.length == self.nPoints) return self;
	
	GLFloat minValue = [self valueAtIndex: range.location];
	GLFloat maxValue = [self valueAtIndex: range.location + range.length-1];
	
	GLDimension *aDim;
	if (self.isEvenlySampled) {
		aDim= [[GLDimension alloc] initDimensionWithGrid: self.gridType nPoints: range.length domainMin: minValue length: maxValue-minValue];
	} else {
		aDim = [[GLDimension alloc] initWithPoints: [self.points objectsAtIndexes: [NSIndexSet indexSetWithIndexesInRange: range]]];
	}
	aDim.name = self.name;
	aDim.units = self.units;
	
	if (!self.rangeSubdimensionMap) {
		self.rangeSubdimensionMap = [NSMapTable strongToStrongObjectsMapTable];
	}
	[self.rangeSubdimensionMap setObject: aDim forKey: [NSValue valueWithRange: range]];
	
	return aDim;
}



- (void) setDimension: (GLDimension *) dim forScale: (GLFloat) scale translation: (GLFloat) delta units: (NSString *) theUnits
{
	if (!_scaledDimensionMapTable) {
		_scaledDimensionMapTable = [NSMapTable mapTableWithKeyOptions: NSMapTableObjectPointerPersonality valueOptions:NSMapTableObjectPointerPersonality];
	}
	
	NSArray *keys = @[@(scale),@(delta), theUnits];
	[_scaledDimensionMapTable setObject: dim forKey: keys];
}

- (GLDimension *) dimensionForScale: (GLFloat) scale translation: (GLFloat) delta units: (NSString *) theUnits
{
	if (!_scaledDimensionMapTable) {
		_scaledDimensionMapTable = [NSMapTable mapTableWithKeyOptions: NSMapTableObjectPointerPersonality valueOptions:NSMapTableObjectPointerPersonality];
	}
	
	for (NSArray *keys in _scaledDimensionMapTable) {
		if ([keys[0] floatValue] == scale && [keys[1] floatValue] == delta && [keys[2] isEqualToString: theUnits]) {
			return [_scaledDimensionMapTable objectForKey: keys];
		}
	}
	return nil;
}

- (GLDimension *) scaledBy: (GLFloat) scale translatedBy: (GLFloat) delta withUnits: (NSString *) newUnits
{
	GLDimension *transformedDim = [self dimensionForScale: scale translation: delta units: newUnits];
	if (!transformedDim) {
		transformedDim = [[GLDimension alloc] initWithDimension: self scaledBy:scale translatedBy:delta withUnits:newUnits];
		[self setDimension: transformedDim forScale:scale translation:delta units: newUnits];
	}
	return transformedDim;
}

/************************************************/
/*		Equality & Comparison					*/
/************************************************/

#pragma mark -
#pragma mark Equality & Comparison
#pragma mark

- (BOOL) isEqualToDimension: (GLDimension *) otherDimension
{
	// Most basic check
	if (self == otherDimension) return YES;
	
	// Next most basic
	if (self.nPoints != otherDimension.nPoints ||
		self.isPeriodic != otherDimension.isPeriodic ||
		self.isEvenlySampled != otherDimension.isEvenlySampled ||
		self.basisFunction != otherDimension.basisFunction ||
        self.isStrictlyPositive != otherDimension.isStrictlyPositive ||
		![self.name isEqualToString: otherDimension.name] ||
		![self.units isEqualToString: otherDimension.units]) {
		return NO;
	}
	
	if (self.isEvenlySampled) {
		return (fabs(self.sampleInterval-otherDimension.sampleInterval) < FLT_EPSILON) ? YES : NO;
	} else {
		GLFloat *f = self.data.mutableBytes;
		GLFloat *fOther = otherDimension.data.mutableBytes;
		for (NSUInteger i=0; i<self.nPoints; i++) {
			if ( fabs(f[i]-fOther[i]) > FLT_EPSILON ) {
				return NO;
			}
		}
		return YES;
	}
}

- (BOOL) isEqual:(id)object
{
	if (self == object) return YES;
	if (![[object class] isSubclassOfClass: [GLDimension class]] ) {
		return NO;
	}
	return [self isEqualToDimension: (GLDimension *) object];
}


/************************************************/
/*		Other									*/
/************************************************/

#pragma mark -
#pragma mark Other
#pragma mark

+ (NSNumber *) uniqueKeyForDimensions: (NSArray *) dimensions
{	
	NSUInteger key = 0;
	for ( GLDimension *aDim in dimensions ) {
		key += aDim.hash;
	}
	return [NSNumber numberWithInteger: key];
}

+ (NSUInteger) strideOfDimensionAtIndex: (NSUInteger) index inArray: (NSArray *) dimensionArray
{
    NSUInteger stride = 1;
    for (NSUInteger i=dimensionArray.count-1; i>index; i--) {
        stride *= [dimensionArray[i] nPoints];
    }
    return stride;
}

+ (NSArray *) rangesFromIndexString:  (NSString *) indexString usingDimensions: (NSArray *) dimensions
{
    indexString = [indexString stringByTrimmingCharactersInSet:[NSCharacterSet whitespaceCharacterSet]];
    NSArray *exploded = [indexString componentsSeparatedByString:@","];
    if (exploded.count != dimensions.count) {
        [NSException raise: @"DimensionsMismatch" format: @"You have specified the incorrect number of dimensions"];
    }
    NSMutableArray *ranges = [NSMutableArray array];
    for ( NSUInteger i=0; i<dimensions.count; i++ ) {
        GLDimension *dim = dimensions[i];
        NSString *string = exploded[i];
        
        if ( [string isEqualToString: @":"] ) {
            ranges[i] = [NSValue valueWithRange: NSMakeRange(0, dim.nPoints)];
        } else {
            NSArray *broken = [string componentsSeparatedByString: @":"];
            if (broken.count == 1) {
                NSInteger a = [broken[0] isEqualToString: @"end"] ? dim.nPoints-1 : [broken[0] integerValue];
                if ( a < 0 || a > dim.nPoints-1) {
                    [NSException raise: @"DimensionsNonsensical" format: @"Index out of range."];
                }
				ranges[i] = [NSValue valueWithRange: NSMakeRange(a,1)];
            } else if (broken.count == 2){
                NSInteger a = [broken[0] isEqualToString: @"end"] ? dim.nPoints-1 : [broken[0] integerValue];
                if ( a < 0 || a > dim.nPoints-1) {
                    [NSException raise: @"DimensionsNonsensical" format: @"Index out of range."];
                }
                
                NSInteger b = [broken[1] isEqualToString: @"end"] ? dim.nPoints-1 : [broken[1] integerValue];
                if ( b < 0 || b > dim.nPoints-1) {
                    [NSException raise: @"DimensionsNonsensical" format: @"Index out of range."];
                }
                
                ranges[i] = [NSValue valueWithRange: NSMakeRange(a, b-a+1)];
            } else {
				[NSException raise: @"DimensionsNonsensical" format: @"I don't understand the dimensions you're asking for."];
			}
        }
    }
    
    return ranges;
}

/************************************************/
/*		Superclass Overrides					*/
/************************************************/

#pragma mark -
#pragma mark Superclass Overrides
#pragma mark

- (NSUInteger)hash {
    NSUInteger hash = 0;
	
    hash += self.name.hash;
    hash += self.units.hash;
	
	hash += 11 * abs(self.domainMin*100); 
	hash += 23 * abs(self.domainLength*100); 
	hash += 31 * abs( 31 - __builtin_clz( (unsigned int) self.nPoints ));
	
	hash += 41 * self.isFrequencyDomain ? 1231:1237;
	
    return hash;
}

- (NSString *) description
{
	NSString *basis;
	if (self.basisFunction == kGLDeltaBasis) {
		basis = @"delta";
	} else if (self.basisFunction == kGLExponentialBasis) {
		basis = @"exponential";
	}
	else if (self.basisFunction == kGLDiscreteCosineTransformIBasis) {
		basis = @"cosine";
	}
	else if (self.basisFunction == kGLDiscreteSineTransformIBasis) {
		basis = @"sine";
	}
	else if (self.basisFunction == kGLCosineBasis) {
		basis = @"cosine half-shift";
	}
	else if (self.basisFunction == kGLSineBasis) {
		basis = @"sine half-shift";
	}
	return [NSString stringWithFormat: @"%@ <0x%lx> (%@: %lu points, %@ basis)", NSStringFromClass([self class]), (NSUInteger)self, self.name, self.nPoints, basis];
}

- (NSString *) graphvisDescription
{
	NSString *basis;
	if (self.basisFunction == kGLDeltaBasis) {
		basis = @"delta";
	} else if (self.basisFunction == kGLExponentialBasis) {
		basis = @"exponential";
	}
	else if (self.basisFunction == kGLDiscreteCosineTransformIBasis) {
		basis = @"cosine";
	}
	else if (self.basisFunction == kGLDiscreteSineTransformIBasis) {
		basis = @"sine";
	}
	else if (self.basisFunction == kGLCosineBasis) {
		basis = @"cosine half-shift";
	}
	else if (self.basisFunction == kGLSineBasis) {
		basis = @"sine half-shift";
	}
	return [NSString stringWithFormat: @"%@: %lu points, %@ basis", self.name, self.nPoints, basis];
}

@end





/************************************************/
/*		GLMutableDimension						*/
/************************************************/

#pragma mark -
#pragma mark GLMutableDimension
#pragma mark

@implementation GLMutableDimension

- (BOOL) isMutable {
	return YES;
}


/************************************************/
/*		Changing Dimensional Length				*/
/************************************************/

#pragma mark -
#pragma mark Changing Dimensional Length
#pragma mark

+ (BOOL)automaticallyNotifiesObserversForKey:(NSString *)theKey {
    
    BOOL automatic = NO;
    if ([theKey isEqualToString:@"domainMin"]) {
        automatic = NO;
    } else if ([theKey isEqualToString:@"domainLength"]) {
        automatic = NO;
    } else if ([theKey isEqualToString:@"isPeriodic"]) {
        automatic = NO;
    } else if ([theKey isEqualToString:@"sampleInterval"]) {
        automatic = NO;
    } else if ([theKey isEqualToString:@"nPoints"]) {
        automatic = NO;
    }  else {
        automatic=[super automaticallyNotifiesObserversForKey:theKey];
    }
    return automatic;
}

- (void) setDomainMin:(GLFloat) min
{
	_domainMin = min;
	[self populateEvenlySampledValues];
}

- (void) setDomainLength:(GLFloat) length
{
	[self willChangeValueForKey: @"sampleInterval"];
	_domainLength = length;
	if (self.isPeriodic)
	{
		if (_nPoints > 0) {
			_sampleInterval = _domainLength / ( (double) _nPoints);
		}
	}
	else
	{
		if (_nPoints > 1) {
			_sampleInterval = _domainLength / ( (double) (_nPoints-1));
		} else {
			_sampleInterval = 0;
		}
	}
	[self didChangeValueForKey: @"sampleInterval"];
	
	[self populateEvenlySampledValues];
}


- (void) setSampleInterval:(GLFloat)anInterval
{
	[self willChangeValueForKey: @"domainLength"];
	[self willChangeValueForKey: @"sampleInterval"];
	_sampleInterval = anInterval;
	if (self.isPeriodic)
	{
		_domainLength = ( (double) _nPoints) * _sampleInterval;
	}
	else
	{
		if (_nPoints == 0) {
			_domainLength = 0;
		} else {
			_domainLength = ( (double) (_nPoints-1)) * _sampleInterval;
		}
	}
	[self didChangeValueForKey: @"sampleInterval"];
	[self didChangeValueForKey: @"domainLength"];
	
	[self populateEvenlySampledValues];
}


//- (void) setIsPeriodic:(BOOL)periodicity
//{
//    [self willChangeValueForKey: @"domainLength"];
//	[self willChangeValueForKey: @"isPeriodic"];
//	isPeriodic = periodicity;
//	if (isPeriodic)
//	{
//		_domainLength = ( (double) _nPoints) * sampleInterval;
//	}
//	else
//	{
//		if (_nPoints == 0) {
//			_domainLength = 0;
//		} else {
//			_domainLength = ( (double) (_nPoints-1)) * sampleInterval;
//		}
//	}
//	[self didChangeValueForKey: @"isPeriodic"];
//	[self didChangeValueForKey: @"domainLength"];
//	
//	[self populateEvenlySampledValues];
//}
//
//- (BOOL) isPeriodic {
//    return [super isPeriodic];
//}


/************************************************/
/*		Increasing Number of Points				*/
/************************************************/

#pragma mark -
#pragma mark Increasing Number of Points
#pragma mark

- (void) setNPoints:(NSUInteger)n
{
	if (_nPoints == n) {
		return;
	}
	[self willChangeValueForKey: @"nPoints"];
	[self willChangeValueForKey: @"domainLength"];
	_nPoints = n;
    _dataBytes = _nPoints*sizeof(GLFloat);
    [self.data setLength: _dataBytes];
	if (self.isPeriodic)
	{
		_domainLength = ( (double) _nPoints) * _sampleInterval;
	}
	else
	{
		if (_nPoints == 0) {
			_domainLength = 0;
		} else {
			_domainLength = ( (double) (_nPoints-1)) * _sampleInterval;
		}
	}
	[self didChangeValueForKey: @"domainLength"];
	[self didChangeValueForKey: @"nPoints"];
	
	[self populateEvenlySampledValues];
} 

- (NSUInteger) nPoints {
    return [super nPoints];
}


- (void) addPoint: (NSNumber *) aPoint
{
    [self addPointsFromArray: [NSArray arrayWithObject: aPoint]];
}

- (void) addPointsFromArray: (NSArray *) pointsArray
{
    [self willChangeValueForKey: @"nPoints"];
    [self willChangeValueForKey: @"domainLength"];
    _nPoints = _nPoints+pointsArray.count;
    _dataBytes = _nPoints*sizeof(GLFloat);
    [self.data setLength: _dataBytes];
    
    GLFloat *f = self.data.mutableBytes;
    NSUInteger i=_nPoints-pointsArray.count;
    for ( NSNumber *point in pointsArray) {
        f[i]= (GLFloat) [point doubleValue];
        i++;
    }
    
	_domainLength = f[_nPoints-1] - f[0];
	[self didChangeValueForKey: @"domainLength"];
    [self didChangeValueForKey: @"nPoints"];
}

//- (void) insertPoint:(NSNumber *) point atIndex: (NSUInteger) index
//{
//	
//}
//
//- (void) insertObjects:(NSArray *) points atIndexes: (NSIndexSet *) indexes
//{
//	
//}


@end
