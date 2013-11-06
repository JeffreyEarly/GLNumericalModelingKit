//
//  GLBuffer.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 11/6/13.
//  Copyright (c) 2013 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface GLBuffer : NSObject
- (GLBuffer *) initWithLength: (NSUInteger) numBytes;
@property NSUInteger numBytes;
@end
