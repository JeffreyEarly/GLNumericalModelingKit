//
//  GLOpenCLFourierTransform.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 5/15/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import "GLOpenCLFourierTransform.h"
#import "GLDimension.h"
#import "GLMemoryPool.h"

#import "clFFT.h"

@interface GLOpenCLFourierTransform ()
@property(readwrite, assign, nonatomic) clFFT_Plan plan;
@end

@implementation GLOpenCLFourierTransform
{
	NSUInteger _dataPoints;
	cl_device_id _device_id;
	cl_context _context;
	cl_command_queue _queue;
	
	clFFT_Plan _plan;
	
	cl_mem _data_in_real;
	cl_mem _data_in_imag;
	cl_mem _data_out_real;
	cl_mem _data_out_imag;
}
@synthesize dimensions;
@synthesize zerosBuffer;
@synthesize dataPoints=_dataPoints;
@synthesize device_id=_device_id;
@synthesize context=_context;
@synthesize queue=_queue;
@synthesize data_in_real=_data_in_real;
@synthesize data_in_imag=_data_in_imag;
@synthesize data_out_real=_data_out_real;
@synthesize data_out_imag=_data_out_imag;
@synthesize plan=_plan;


cl_device_type getGlobalDeviceType()
{
	char *force_cpu = getenv( "CL_DEVICE_TYPE" );
	if( force_cpu != NULL )
	{
		if( strcmp( force_cpu, "gpu" ) == 0 || strcmp( force_cpu, "CL_DEVICE_TYPE_GPU" ) == 0 )
			return CL_DEVICE_TYPE_GPU;
		else if( strcmp( force_cpu, "cpu" ) == 0 || strcmp( force_cpu, "CL_DEVICE_TYPE_CPU" ) == 0 )
			return CL_DEVICE_TYPE_CPU;
		else if( strcmp( force_cpu, "accelerator" ) == 0 || strcmp( force_cpu, "CL_DEVICE_TYPE_ACCELERATOR" ) == 0 )
			return CL_DEVICE_TYPE_ACCELERATOR;
		else if( strcmp( force_cpu, "CL_DEVICE_TYPE_DEFAULT" ) == 0 )
			return CL_DEVICE_TYPE_DEFAULT;
	}
	// default
	return CL_DEVICE_TYPE_GPU;
}

- (GLFourierTransform *) initWithDimensions: (NSArray *) theDimensions
{
	if ( sizeof(GLFloat) != sizeof(float) ) {
		NSLog(@"OpenCL FFT only supports single precision floating point.");
		return nil;
	}
	
	// Apparently this FFT routine *only* works on GPUs. So check that we have a GPU.
	cl_device_type device_type = getGlobalDeviceType();	
	if(device_type != CL_DEVICE_TYPE_GPU)  {
		NSLog(@"Test only supported on DEVICE_TYPE_GPU\n");
		return nil;
	}
	
	// Now get the list of available devices.
	cl_device_id device_ids[16];
	unsigned int num_devices;
	cl_int err = clGetDeviceIDs(NULL, device_type, sizeof(device_ids), device_ids, &num_devices);
	if(err) {		
		NSLog(@"clGetComputeDevice failed");
		return nil;
	}
	
	// Loop through the devices and grab the first one that we can compute with.
	cl_device_id device_id;
	for(unsigned int i = 0; i < num_devices; i++)
	{
	    cl_bool available;
	    err = clGetDeviceInfo(device_ids[i], CL_DEVICE_AVAILABLE, sizeof(cl_bool), &available, NULL);
	    if(err) {
			NSLog(@"Cannot check device availability of device # %d\n", i);
	    }
	    
	    if(available) {
	        device_id = device_ids[i];
	        break;
	    }
	    else {
	        char name[200];
	        err = clGetDeviceInfo(device_ids[i], CL_DEVICE_NAME, sizeof(name), name, NULL);
	        if(err == CL_SUCCESS) {
				NSLog(@"Device %s not available for compute\n", name);
	        }
	        else {
				NSLog(@"Device # %d not available for compute\n", i);
	        }
	    }
	}
	
	// Hmmm, didn't find one available, so we need to bail.
	if (!device_id) {
		NSLog(@"None of the devices available for compute ... aborting");
		return nil;
	}
	
	cl_context context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
	if(!context || err) {
		NSLog(@"clCreateContext failed");
		return nil;
	}
	
    cl_command_queue queue = clCreateCommandQueue(context, device_id, 0, &err);
    if(!queue || err) {
        NSLog(@"clCreateCommandQueue() failed.\n");
		clReleaseContext(context);
		return nil;
    } 
	
	cl_ulong gMemSize;
	err = clGetDeviceInfo(device_id, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &gMemSize, NULL);
	if(err) {
		NSLog(@"Failed to get global mem size\n");
		clReleaseContext(context);
		clReleaseCommandQueue(queue);
		return nil;
	}
	gMemSize /= (1024*1024);
	
	// Now check the memory requirements.
	unsigned int dims[3];
	clFFT_Dim3 clDims = { 1, 1, 1 };
	NSUInteger dataPoints = 1;
	NSUInteger i=0;
	for ( GLDimension *aDim in theDimensions ) {
		dataPoints = dataPoints * aDim.nPoints;
		vDSP_Length logDim = 31 - __builtin_clz( (unsigned int) aDim.nPoints );
		dims[i] = (unsigned int) logDim;
		i++;
	}
	
	// A temporary hack since we only support two dimensions anyway.
	// This FFT implementation indexes everything differently than we do.
	clFFT_Dimension clFFTDim;
	if (theDimensions.count) {
		clFFTDim = clFFT_1D;
		clDims.x = (unsigned int)[[theDimensions objectAtIndex: 1] nPoints];
		if (theDimensions.count > 1) {
			clFFTDim = clFFT_2D;
			clDims.y = (unsigned int)[[theDimensions objectAtIndex: 0] nPoints];
			if (theDimensions.count > 2) {
				clFFTDim = clFFT_3D;
				clDims.z = (unsigned int)[[theDimensions objectAtIndex: 2] nPoints];
			}
		}
	}
	
	// We're assuming out-of-place FFT and a "batchSize" of 1.
	cl_ulong memReq = 3;
	memReq *= clDims.x*clDims.y*clDims.z*sizeof(clFFT_Complex);
	memReq = memReq/1024/1024;
	if(memReq >= gMemSize) {
		NSLog(@"FFT memory requires exceeds device memory");
		return nil;
	}
	
	if ((self=[super init])) {
		self.dimensions = theDimensions;
		self.dataPoints = dataPoints;
		self.device_id = device_id;
		self.context = context;
		self.queue = queue;
		
		self.zerosBuffer = [[GLMemoryPool sharedMemoryPool] dataWithLength: self.dataPoints*sizeof(GLFloat)];
		self.plan = clFFT_CreatePlan( self.context, clDims, clFFTDim, clFFT_SplitComplexFormat, &err );
		if(!self.plan || err)  {
			NSLog(@"clFFT_CreatePlan failed");
			return nil;
		}
				
	}
	return self;
}

- (void) transform: (GLFloat *) f forward: (GLSplitComplex *) fbar
{
	// Do an out-of-place FFT
	GLSplitComplex split;
	split.realp = f;
	split.imagp = self.zerosBuffer.mutableBytes;
	vGL_vclr( split.imagp, 1, self.dataPoints );
	
	cl_int err;
	cl_mem data_in_real = clCreateBuffer(self.context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, self.dataPoints*sizeof(float), split.realp, &err);
	cl_mem data_in_imag = clCreateBuffer(self.context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, self.dataPoints*sizeof(float), split.imagp, &err);
	cl_mem data_out_real = clCreateBuffer(self.context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, self.dataPoints*sizeof(float), fbar->realp, &err);
	cl_mem data_out_imag = clCreateBuffer(self.context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, self.dataPoints*sizeof(float), fbar->imagp, &err);
	
	if(!data_in_real || !data_in_imag || !data_out_real || !data_out_imag || err)  {
		NSLog(@"clCreateBuffer failed\n");
		return;
	}
	
	err |= clFFT_ExecutePlannar(self.queue, self.plan, 1, clFFT_Forward, data_in_real, data_in_imag, data_out_real, data_out_imag, 0, NULL, NULL);
	
	err |= clEnqueueReadBuffer(self.queue, data_out_real, CL_TRUE, 0, self.dataPoints*sizeof(float), fbar->realp, 0, NULL, NULL);
	err |= clEnqueueReadBuffer(self.queue, data_out_imag, CL_TRUE, 0, self.dataPoints*sizeof(float), fbar->imagp, 0, NULL, NULL);
	
	clReleaseMemObject(data_in_real);
	clReleaseMemObject(data_in_imag);
	clReleaseMemObject(data_out_real);
	clReleaseMemObject(data_out_imag);
	
	GLFloat inversesize = 1.0 / ( (GLFloat) self.dataPoints );
	vGL_vsmul( fbar->realp, 1, &inversesize, fbar->realp, 1, self.dataPoints );
	vGL_vsmul( fbar->imagp, 1, &inversesize, fbar->imagp, 1, self.dataPoints );
}

- (void) transform: (GLSplitComplex *) fbar inverse: (GLFloat *) f
{
	GLSplitComplex split;
	split.realp = f;
	split.imagp = self.zerosBuffer.mutableBytes;
	
	cl_int err;
	cl_mem data_in_real = clCreateBuffer(self.context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, self.dataPoints*sizeof(float), fbar->realp, &err);
	cl_mem data_in_imag = clCreateBuffer(self.context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, self.dataPoints*sizeof(float), fbar->imagp, &err);
	cl_mem data_out_real = clCreateBuffer(self.context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, self.dataPoints*sizeof(float), split.realp, &err);
	cl_mem data_out_imag = clCreateBuffer(self.context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, self.dataPoints*sizeof(float), split.imagp, &err);
	
	if(!data_in_real || !data_in_imag || !data_out_real || !data_out_imag || err)  {
		NSLog(@"clCreateBuffer failed\n");
		return;
	}
	
	err |= clFFT_ExecutePlannar(self.queue, self.plan, 1, clFFT_Inverse, data_in_real, data_in_imag, data_out_real, data_out_imag, 0, NULL, NULL);
	
	err |= clEnqueueReadBuffer(self.queue, data_out_real, CL_TRUE, 0, self.dataPoints*sizeof(float), split.realp, 0, NULL, NULL);
	
	clReleaseMemObject(data_in_real);
	clReleaseMemObject(data_in_imag);
	clReleaseMemObject(data_out_real);
	clReleaseMemObject(data_out_imag);
	
//	GLFloat inversesize = 1.0 / ( (GLFloat) self.dataPoints );
//	vGL_vsmul( f, 1, &inversesize, f, 1, self.dataPoints );
}

- (void) dealloc
{
	clFFT_DestroyPlan(self.plan);
	clReleaseContext(self.context);
	clReleaseCommandQueue(self.queue);
}

@end
