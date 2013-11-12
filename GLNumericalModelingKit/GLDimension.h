//
//  GLDimension.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 4/21/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/Precision.h>

// Dimensional policies:
// LOOSE	A 2x6 matrix can be added to a 4x3 matrix because they have the same data lengths.
// NORMAL	Dimensional lengths must match.
// STRICT	Dimensional units, and colocation points must align.

// Implementation notes
// --------------------
//
// Dimensions are independent of all other classes, including the GLMemoryPool. It is assummed that
// dimensions will generally be small enough, and created infrequently enough, that there is not an
// advantage to using the memory pool. They are always memory backed with their own memory buffer.
//
// It is similarly assumed that an individual dimension may be associated with multiple variables 
// and multiple NetCDF files.

typedef NS_ENUM(NSUInteger, GLGridType) {
	kGLEndpointGrid = 0,            // A good grid for finite differencing
	kGLInteriorGrid = 1,            // A good grid for transforming to a sine or cosine basis
	kGLPeriodicGrid = 2,            // A good grid for transforming to an exponential basis.
	kGLChebyshevEndpointGrid = 3,   // aka, 'extrema' or 'Lobatto' grid
	kGLChebyshevInteriorGrid = 4,   // aka, 'roots' or 'Gauss' grid
    kGLUnevenGrid = 10              // Some other non-standard grid.
};

typedef NS_ENUM(NSUInteger, GLBasisFunction) {
	kGLDeltaBasis = 0,			// Spatial domain
    kGLExponentialBasis = 1,	// Frequency domain
	kGLCosineBasis = 2,			// Frequency domain
	kGLSineBasis = 3,			// Frequency domain
	kGLCosineHalfShiftBasis = 4,// Frequency domain
	kGLSineHalfShiftBasis = 5,	// Frequency domain
	kGLChebyshevBasis = 6
};

@class GLMutableDimension;
@interface GLDimension : NSObject

/************************************************/
/*		Convenience Methods						*/
/************************************************/

#pragma mark -
#pragma mark Convenience Methods
#pragma mark

// Pre-named, evenly spaced, periodic, min=0, length=1.
+ (GLDimension *) dimensionXWithNPoints: (NSUInteger) numPoints length: (GLFloat) theLength;
+ (GLDimension *) dimensionYWithNPoints: (NSUInteger) numPoints length: (GLFloat) theLength;
+ (GLDimension *) dimensionZWithNPoints: (NSUInteger) numPoints length: (GLFloat) theLength;

// Pre-named, evenly spaced, non-periodic, nPoints=1, sample interval 1.
+ (GLMutableDimension *) dimensionT;

/************************************************/
/*		Initialization							*/
/************************************************/

#pragma mark -
#pragma mark Initialization
#pragma mark

// Dimensions should be held onto for as long as they're used. In other words,
// don't just reinitialize a new one if you lost track

// Methods for creating *evenly spaced* dimensions.
- (GLDimension *) initDimensionWithGrid: (GLGridType) gridType nPoints: (NSUInteger) numPoints domainMin: (GLFloat) theMin length: (GLFloat) theLength;

// Create the corresponding transformed dimension in a different basis and optionally restrict the dimension to positive points (useful for common symmetries).
- (GLDimension *) initAsDimension: (GLDimension *) existingDimension transformedToBasis: (GLBasisFunction) basis strictlyPositive: (BOOL) isPositive;

// Method for create the fourier transformed version of this dimension.
// This is a fully invertible process, the two related dimensions hold on to each other.
- (GLDimension *) initAsFourierTransformOfDimension: (GLDimension *) existingDimension;

- (GLDimension *) initWithDimension: (GLDimension *) existingDim scaledBy: (GLFloat) scale translatedBy: (GLFloat) delta withUnits: (NSString *) newUnits;

// Methods for creating *unevenly spaced* dimensions. Evenly spaced dimensions may benefit from
// some optimizations.
- (GLDimension *) initWithNPoints: (NSUInteger) numPoints values: (NSData *) values;
- (GLDimension *) initWithPoints: (NSArray *) pointsArray;

/************************************************/
/*		Metadata								*/
/************************************************/

#pragma mark -
#pragma mark Metadata
#pragma mark

// Name of the dimension, e.g. "x".
@property(readwrite, copy, nonatomic) NSString *name;

// Units of the dimension, e.g. "meters"
@property(readwrite, copy, nonatomic) NSString *units;

/************************************************/
/*		Essential Properties					*/
/************************************************/

#pragma mark -
#pragma mark Essential Properties
#pragma mark

// This properties cannot be changed after initialization.
// If you need to vary the dimension, use a GLMutableDimension

// The number of discrete points, e.g. 128
@property(readonly, nonatomic) NSUInteger nPoints;

// The dimensional starting point and size of the dimension, e.g. min = 0 meters & length = 100 meters.
@property(readonly, nonatomic) GLFloat domainMin;
@property(readonly, nonatomic) GLFloat domainLength;

@property(readwrite, nonatomic) GLGridType gridType;

// Changing the periodicity will change the domain length because a periodic domain
// is one sample interval longer than a non-periodic domain.
@property(readonly, nonatomic) BOOL isPeriodic;

// This cannot be changed after initialization.
@property(readonly, nonatomic) BOOL isEvenlySampled;

// This is only valid if the dimension is evenly sampled.
// Change the sampleInterval will change the domainLength, and viceversa.
@property(readonly, nonatomic) GLFloat sampleInterval;

// You can only initialize dimensions in the time/space domain, so only by fetching the
// fourier transformated dimension is this ever set to YES.
@property(readonly, nonatomic) BOOL isFrequencyDomain;

@property(readonly, nonatomic) BOOL isMutable;

@property(readwrite, nonatomic) GLBasisFunction basisFunction;

/// The default basis which should be used to take a derivative. This can only be changed for spatial dimensions.
@property(readwrite, nonatomic) GLBasisFunction differentiationBasis;

@property(readwrite, nonatomic) BOOL isStrictlyPositive;

/************************************************/
/*		Value									*/
/************************************************/

#pragma mark -
#pragma mark Value
#pragma mark

// Returns the dimensional value at a given index.
- (GLFloat) valueAtIndex: (NSUInteger) index;

// An array of nPoints of GLFloats
@property(readonly, strong, nonatomic) NSMutableData *data;
@property(readonly, nonatomic) NSUInteger dataBytes;

// An array of NSNumbers, generated upon request.
@property(readonly, nonatomic) NSArray *points;

/************************************************/
/*		Derived Dimensions						*/
/************************************************/

#pragma mark -
#pragma mark Derived Dimensions
#pragma mark

// Returns the fourier transformed version of this dimesion.
// This is a fully invertible process, the two related dimensions hold on to each other.
@property(readonly) GLDimension *fourierTransformedDimension;

// This will return nil if the range is out of bounds.
// It will return self if the range includes the whole dimension.
- (GLDimension *) subdimensionWithRange: (NSRange) range;

// The dimension will be scaled first, then translated.
- (GLDimension *) scaledBy: (GLFloat) scale translatedBy: (GLFloat) delta withUnits: (NSString *) newUnits;

/************************************************/
/*		Equality & Comparison					*/
/************************************************/

#pragma mark -
#pragma mark Equality & Comparison
#pragma mark

// Performs a comprehensive check, if necessary.
- (BOOL) isEqualToDimension: (GLDimension *) otherDimension;

/************************************************/
/*		Other									*/
/************************************************/

#pragma mark -
#pragma mark Other
#pragma mark

// Return a unique NSNumber for an ordered array of dimensions.
// The number is unique only up to the number of points in each dimensions.
+ (NSNumber *) uniqueKeyForDimensions: (NSArray *) dimensions;

+ (NSUInteger) strideOfDimensionAtIndex: (NSUInteger) index inArray: (NSArray *) dimensionArray;

// Convert a Matlab style index string into an array of ranges. String should be @"start:end,start:end,..."
// Legal values are integers, :, and end.
+ (NSArray *) rangesFromIndexString:  (NSString *) indexString usingDimensions: (NSArray *) dimensions;

@property(readonly) NSString *graphvisDescription;

@end


/************************************************/
/*		GLMutableDimension						*/
/************************************************/

#pragma mark -
#pragma mark GLMutableDimension
#pragma mark

@interface GLMutableDimension : GLDimension

/************************************************/
/*		Changing Dimensional Length				*/
/************************************************/

#pragma mark -
#pragma mark Changing Dimensional Length
#pragma mark

// Changing these parameters will changes the values at the different points,
// but will not effect the number of points.

// The dimensional starting point and size of the dimension, e.g. min = 0 meters & length = 100 meters.
@property(readwrite, nonatomic) GLFloat domainMin;
@property(readwrite, nonatomic) GLFloat domainLength;

// Changing the periodicity will change the domain length because a periodic domain
// is one sample interval longer than a non-periodic domain.
@property(readwrite, nonatomic) BOOL isPeriodic;

// This is only valid if the dimension is evenly sampled.
// Change the sampleInterval will change the domainLength, and viceversa.
@property(readwrite, nonatomic) GLFloat sampleInterval;

/************************************************/
/*		Increasing Number of Points				*/
/************************************************/

#pragma mark -
#pragma mark Increasing Number of Points
#pragma mark

// Changing this value holds the sampleInterval fixed, and varies the domainLength.
// If your dimension is evenly sampled, this is how you change its length.
@property(readwrite, nonatomic) NSUInteger nPoints;

- (void) addPoint: (NSNumber *) aPoint;

- (void) addPointsFromArray: (NSArray *) points;

//- (void) insertPoint:(NSNumber *) point atIndex: (NSUInteger) index;
//
//- (void) insertObjects:(NSArray *) points atIndexes: (NSIndexSet *) indexes;

@end

