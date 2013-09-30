//
//  GLFourierTransformPool.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 4/22/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>

@class GLFourierTransform;
@interface GLFourierTransformPool : NSObject

// Access the shared instance which maintains the pool of available data objects
+ (id)sharedFourierTransformPool;

- (id) initWithMaxTransforms: (NSUInteger) maxT forClass: (Class) aClass;

// Returns a transform object that responds to:
//	- (void) transform: (GLFloat *) f forward: (GLSplitComplex *) fbar;
//	- (void) transform: (GLSplitComplex *) fbar inverse: (GLFloat *) f;
- (GLFourierTransform *) transformWithDimensions: (NSArray *) dimensions;

// Returns a unique key for the dimensions. All dimensions that have this key can use the same fourier transform.
- (NSNumber *) keyForDimensions: (NSArray *) dimensions;

// This converts the key back to an array of GLDimension objects. Only the number of points in these dimensions
// will be recovered from the key---no other information.
- (NSArray *) dimensionsForKey: (NSNumber *) number;

// Request a fourier transform directly, based on the key.
- (GLFourierTransform *) transformForKey: (NSNumber *) key;

// Return the transform so that it can be used elsewhere.
- (void) returnTransform: (id) transform;

// Destroy the entire pool.
- (void) destroy;

/*******************************************************/
//			Private
/*******************************************************/

@property dispatch_queue_t poolQueue;

// The key of this dictionary is an NSNumber corresponding to the length of the data.
// The dictionary returns an array of NSMutableData objects that have been initialized to that length.
@property(strong) NSMutableDictionary *bytesArrayDictionary;

@property dispatch_semaphore_t semaphore;
//@property Class protoClass;

@end
