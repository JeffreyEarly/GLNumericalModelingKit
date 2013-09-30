//
//  GLMemoryPool.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 4/18/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import "GLMemoryPool.h"
#import <fftw3.h>

#define FFT_IS_ALIGNED(p, alignSize)		((((intptr_t)(p)) & (alignSize-1)) == 0)

//@interface GLMutableData : NSMutableData
//
//@end
//
//@implementation GLMutableData
//
//- (NSMutableData *) initWithCapacity: (NSUInteger) length
//{
//	return [
//}
//
//-(void) dealloc
//{
//	NSLog(@"Why oh why should I dealloc");
//}
//
//@end


@implementation GLMemoryPool

@synthesize poolQueue;
@synthesize bytesArrayDictionary;

+ (id)sharedMemoryPool
{
    static dispatch_once_t pred;
    static GLMemoryPool *memoryPool = nil;
    
    dispatch_once(&pred, ^{ memoryPool = [[self alloc] init]; });
    return memoryPool;
}


- (id) init
{
    if ((self=[super init]))
    {
        self.poolQueue = dispatch_queue_create( "com.EarlyInnovations.SharedMemoryPoolQueue", NULL );
        self.bytesArrayDictionary = [NSMutableDictionary dictionary];
    }
    return self;
}

- (NSData *) dataWithLength: (NSUInteger) numBytes
{
    __block NSData *data;
    __block NSHashTable *dataArray;
    
    dispatch_sync( self.poolQueue, ^{
        dataArray = [self.bytesArrayDictionary objectForKey: [NSNumber numberWithUnsignedInteger: numBytes]];
        if (!dataArray) {
            dataArray = [NSHashTable hashTableWithOptions: NSHashTableObjectPointerPersonality];
            [self.bytesArrayDictionary setObject: dataArray forKey: [NSNumber numberWithUnsignedInteger: numBytes]];
        }
        
		data = [dataArray anyObject];
        if (!data) {
            data = [NSMutableData dataWithLength: numBytes];
			
//			void *f = fftwf_malloc(numBytes);
//			data = [NSMutableData dataWithBytesNoCopy: f length:numBytes freeWhenDone:YES];
			
			if (!FFT_IS_ALIGNED(data.bytes, 16)) {
				NSLog(@"Data is not aligned");
			}
        } else {
			[dataArray removeObject: data];
		}
    });
    
    return data;
}

- (void) returnData: (NSData *) data
{    
    dispatch_async( self.poolQueue, ^{
         NSHashTable *dataArray = [self.bytesArrayDictionary objectForKey: [NSNumber numberWithUnsignedInteger: data.length]];
        if (!dataArray) {
            dataArray = [NSHashTable hashTableWithOptions: NSHashTableObjectPointerPersonality];
            [self.bytesArrayDictionary setObject: dataArray forKey: [NSNumber numberWithUnsignedInteger: data.length]];
        }
        
        [dataArray addObject: data];
    });
}

- (void) destroy
{
    dispatch_sync( self.poolQueue, ^{
        self.bytesArrayDictionary = [NSMutableDictionary dictionary];
    });
}

@end
