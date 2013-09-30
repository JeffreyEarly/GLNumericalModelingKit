//
//  GLOpenCLFourierTransform.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 5/15/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import "GLFourierTransform.h"
#include <OpenCL/opencl.h>

@interface GLOpenCLFourierTransform : GLFourierTransform

- (GLFourierTransform *) initWithDimensions: (NSArray *) theDimensions;

@property(readwrite, strong, nonatomic) NSArray *dimensions;
@property(readwrite, strong, nonatomic) NSMutableData *zerosBuffer;
@property(readwrite, assign, nonatomic) NSUInteger dataPoints;

@property(readwrite, assign, nonatomic) cl_device_id device_id;
@property(readwrite, assign, nonatomic) cl_context context;
@property(readwrite, assign, nonatomic) cl_command_queue queue;

@property(readwrite, assign, nonatomic) cl_mem data_in_real;
@property(readwrite, assign, nonatomic) cl_mem data_in_imag;
@property(readwrite, assign, nonatomic) cl_mem data_out_real;
@property(readwrite, assign, nonatomic) cl_mem data_out_imag;

@end
