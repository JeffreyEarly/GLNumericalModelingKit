//
//  GLDSPFourierTransform.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 10/11/11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#import "GLDSPFourierTransform.h"
#import "GLDimension.h"
#import "GLMemoryPool.h"

@implementation GLDSPFourierTransform

@synthesize dimensions;
@synthesize zerosBuffer;
@synthesize setup;
@synthesize logDimensions;
@synthesize dataPoints;

- (GLFourierTransform *) initWithDimensions: (NSArray *) theDimensions
{
	if ((self=[super init])) {
		self.dimensions = theDimensions;
		self.logDimensions = [NSMutableArray array];
		
		vDSP_Length maxLogDim = 0;
		self.dataPoints = 1;
		for ( GLDimension *aDim in self.dimensions ) {
			self.dataPoints = self.dataPoints * aDim.nPoints;
			vDSP_Length logDim = 31 - __builtin_clz( (unsigned int) aDim.nPoints );
			[self.logDimensions addObject: [NSNumber numberWithUnsignedLong: logDim]];
			if (logDim > maxLogDim) maxLogDim = logDim;
		}
		
		self.setup = vGL_create_fftsetup( maxLogDim, kFFTRadix2 );
		
		self.zerosBuffer = [[GLMemoryPool sharedMemoryPool] dataWithLength: self.dataPoints*sizeof(GLFloat)];
		
		
	}
	return self;
}

- (void) transform: (GLFloat *) f forward: (GLSplitComplex *) fbar
{
	if (self.dimensions.count != 2) {
		[NSException raise: @"BadFFTDimensionsException" format: @"Cannot transform anything other than two dimensions!"];
        return;
	}
	// Do an out-of-place FFT with a garbage buffer
	GLSplitComplex split;
	split.realp = f;
	split.imagp = self.zerosBuffer.mutableBytes;
	vGL_vclr( split.imagp, 1, self.dataPoints );
	
	vDSP_Length logX = [[self.logDimensions objectAtIndex: 0] unsignedLongValue];
	vDSP_Length logY = [[self.logDimensions objectAtIndex: 1] unsignedLongValue];
	vGL_fft2d_zop(self.setup, &split, 1, 0, fbar, 1, 0, logY, logX, FFT_FORWARD);
	
	// Fix the forward transformation scaling factor
	// Yes, the *forward*, not the inverse. If we don't fix it here, then the fourier transform
	// is dependent on the number of discrete points, which isn't helpful.
	GLFloat inversesize = 1.0 / ( (GLFloat) self.dataPoints );
	vGL_vsmul( fbar->realp, 1, &inversesize, fbar->realp, 1, self.dataPoints );
	vGL_vsmul( fbar->imagp, 1, &inversesize, fbar->imagp, 1, self.dataPoints );
}

- (void) transform: (GLSplitComplex *) fbar inverse: (GLFloat *) f
{
	if (self.dimensions.count != 2) {
		[NSException raise: @"BadFFTDimensionsException" format: @"Cannot transform anything other than two dimensions!"];
        return;
	}
	
	GLSplitComplex split;
	split.realp = f;
	split.imagp = self.zerosBuffer.mutableBytes;
	
	vDSP_Length logX = [[self.logDimensions objectAtIndex: 0] unsignedLongValue];
	vDSP_Length logY = [[self.logDimensions objectAtIndex: 1] unsignedLongValue];
	vGL_fft2d_zop(self.setup, fbar, 1, 0, &split, 1, 0, logY, logX, FFT_INVERSE);
}

- (void) dealloc
{
	[[GLMemoryPool sharedMemoryPool] returnData: self.zerosBuffer];
	vGL_destroy_fftsetup( setup );	
}

@end
