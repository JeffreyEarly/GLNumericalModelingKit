//
//  GLBuffer.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 11/6/13.
//  Copyright (c) 2013 Jeffrey J. Early. All rights reserved.
//

#import "GLBuffer.h"

@implementation GLBuffer
- (GLBuffer *) initWithLength: (NSUInteger) numBytes {
	if ((self=[super init])) {
		self.numBytes = numBytes;
	}
	return self;
}
@end
