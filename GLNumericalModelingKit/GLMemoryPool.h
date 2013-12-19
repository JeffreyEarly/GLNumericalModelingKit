//
//  GLMemoryPool.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 4/18/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>

// You don't have to use

@interface GLMemoryPool : NSObject

// Access the shared instance which maintains the pool of available data objects
+ (id)sharedMemoryPool;

/// Request an empty (but not necessarily zeroed) data object of a particular size in bytes.
- (NSMutableData *) dataWithLength: (NSUInteger) numBytes;

// Return the data object so that it can be used elsewhere.
- (void) returnData: (NSData *) data;

// Destroy the entire pool.
- (void) destroy;


/*******************************************************/
//			Private
/*******************************************************/

@property dispatch_queue_t poolQueue;

// The key of this dictionary is an NSNumber corresponding to the length of the data.
// The dictionary returns an array of NSMutableData objects that have been initialized to that length.
@property(strong) NSMutableDictionary *bytesArrayDictionary;

@end
