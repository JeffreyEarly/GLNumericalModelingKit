//
//  GLFourierTransformPool.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 4/22/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import "GLFourierTransformPool.h"
#import "GLFourierTransform.h"
#import "GLMatrixFourierTransform.h"
#import "GLOpenCLFourierTransform.h"
#import "GLDimension.h"

@implementation GLFourierTransformPool

+ (id)sharedFourierTransformPool
{
    static dispatch_once_t pred;
    static GLFourierTransformPool *transformPool = nil;
    
    dispatch_once(&pred, ^{ transformPool = [[self alloc] initWithMaxTransforms: 10 forClass: nil]; });
    return transformPool;
}

@synthesize poolQueue;
@synthesize semaphore;
@synthesize bytesArrayDictionary;

- (id) init
{
    if ((self=[super init]))
    {
        self.poolQueue = dispatch_queue_create( "com.EarlyInnovations.SharedFourierTransformPoolQueue", NULL );
        self.bytesArrayDictionary = [NSMutableDictionary dictionary];
    }
    return self;
}

- (id) initWithMaxTransforms: (NSUInteger) maxT forClass: (Class) aClass
{
	if ((self=[super init]))
    {
        self.poolQueue = dispatch_queue_create( "com.EarlyInnovations.SharedFourierTransformPoolQueue", NULL );
        self.bytesArrayDictionary = [NSMutableDictionary dictionary];
		self.semaphore = dispatch_semaphore_create(maxT);
    }
    return self;
}

- (NSNumber *) keyForDimensions: (NSArray *) dimensions
{	
	NSUInteger key = 0;
	NSInteger i = 0;
	for ( GLDimension *aDim in dimensions ) {
		// a is now the log2(nPoints)
		NSUInteger a = 31 - __builtin_clz( (unsigned int) aDim.nPoints );
		key |= a << (i*8);
		i++;
	}
	return [NSNumber numberWithUnsignedInteger: key];
}

- (NSArray *) dimensionsForKey: (NSNumber *) number
{
	NSMutableArray *array = [NSMutableArray array];
	
	NSUInteger key = [number unsignedIntegerValue];
	NSUInteger mask = 0xFF;
	
	while (key & mask) {
		unsigned long int theLog = key & mask;
		NSUInteger nPoints = 1;
		nPoints = nPoints << theLog;
		[array addObject: [GLDimension dimensionXWithNPoints: nPoints length: 1.0]];
		key = key >> 8;
	}
	
	return array;
}

- (GLFourierTransform *) transformWithDimensions: (NSArray *) dimensions
{
    __block GLFourierTransform *transform = nil;
    __block NSHashTable *transformArray = nil;
    
	if (self.semaphore) {
		dispatch_semaphore_wait(self.semaphore, DISPATCH_TIME_FOREVER);
	}
    dispatch_sync( self.poolQueue, ^{
        transformArray = [self.bytesArrayDictionary objectForKey: [self keyForDimensions: dimensions]];
        if (!transformArray) {
            transformArray = [NSHashTable hashTableWithOptions: NSHashTableObjectPointerPersonality];
//			[transformArray addObject: [[GLMatrixFourierTransform alloc] initWithDimensions: dimensions]];
//			[transformArray addObject: [[GLMatrixFourierTransform alloc] initWithDimensions: dimensions]];
//			[transformArray addObject: [[GLOpenCLFourierTransform alloc] initWithDimensions: dimensions]];
            [self.bytesArrayDictionary setObject: transformArray forKey: [self keyForDimensions: dimensions]];
        }
        
        transform = [transformArray anyObject];
        if (!transform) {
            transform = [GLFourierTransform fourierTransformWithDimensions: dimensions];
        } else {
			[transformArray removeObject: transform];
		}
    });
    

	
    return transform;
}

- (GLFourierTransform *) transformForKey: (NSNumber *) key
{
    __block GLFourierTransform *transform = nil;
    __block NSHashTable *transformArray = nil;
    
	if (self.semaphore) {
		dispatch_semaphore_wait(self.semaphore, DISPATCH_TIME_FOREVER);
	}
    dispatch_sync( self.poolQueue, ^{
        transformArray = [self.bytesArrayDictionary objectForKey: key];
        if (!transformArray) {
            transformArray = [NSHashTable hashTableWithOptions: NSHashTableObjectPointerPersonality];
//			[transformArray addObject: [[GLMatrixFourierTransform alloc] initWithDimensions: [self dimensionsForKey:key]]];
//			[transformArray addObject: [[GLMatrixFourierTransform alloc] initWithDimensions: [self dimensionsForKey:key]]];
//			[transformArray addObject: [[GLOpenCLFourierTransform alloc] initWithDimensions: [self dimensionsForKey:key]]];
            [self.bytesArrayDictionary setObject: transformArray forKey: key];
        }
        
        transform = [transformArray anyObject];
        if (!transform) {
            transform = [GLFourierTransform fourierTransformWithDimensions: [self dimensionsForKey:key]];
        } else {
			[transformArray removeObject: transform];
		}
    });
    
    return transform;
}

- (void) returnTransform: (id) transform;
{
    __block NSHashTable *dataArray = nil;
    
    dispatch_sync( self.poolQueue, ^{
        dataArray = [self.bytesArrayDictionary objectForKey: [self keyForDimensions: [transform dimensions]]];
        if (!dataArray) {
            dataArray = [NSHashTable hashTableWithOptions: NSHashTableObjectPointerPersonality];
            [self.bytesArrayDictionary setObject: dataArray forKey: [self keyForDimensions: [transform dimensions]]];
        }
        
        [dataArray addObject: transform];
    });
	
	if (self.semaphore) {
		dispatch_semaphore_signal(self.semaphore);
	}
}

- (void) destroy
{
    dispatch_sync( self.poolQueue, ^{
        self.bytesArrayDictionary = [NSMutableDictionary dictionary];
    });
}

@end
