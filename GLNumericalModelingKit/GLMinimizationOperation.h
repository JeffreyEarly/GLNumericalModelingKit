//
//  GLMinimizationOperation.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 11/14/13.
//  Copyright (c) 2013 Jeffrey J. Early. All rights reserved.
//

#import <GLNumericalModelingKit/GLNumericalModelingKit.h>

typedef GLScalar * (^yFromX)(NSArray *);

@interface GLMinimizationOperation : GLVariableOperation

- (GLMinimizationOperation *) initAtPoint: (NSArray *) startPoint withDeltas: (NSArray *) deltas forFunction: (yFromX) functionBlock;

@end
