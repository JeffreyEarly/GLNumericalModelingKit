//
//  GLFFTWFourierTransform.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 5/16/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import "GLFFTWFourierTransform.h"
#import "GLDimension.h"
#import "GLMemoryPool.h"
#import <fftw3.h>

@interface GLFFTWFourierTransform ()
@property(readwrite, assign, nonatomic) fftw_iodim *iodims;
@end

@implementation GLFFTWFourierTransform

@synthesize dimensions;
@synthesize zerosBuffer;
@synthesize iodims;
@synthesize dataPoints;

- (GLFourierTransform *) initWithDimensions: (NSArray *) theDimensions
{
	if ( theDimensions.count > 3 ) {
		[NSException raise: @"NumberOfDimensionsUnsupportedException" format: @"MatrixFFT only supports up to two dimensions."];
		return nil;
	}
	
	if ((self=[super init])) {
		self.dimensions = theDimensions;
		
		self.iodims = malloc(3*sizeof(fftw_iodim));
		self.dataPoints = 1;
		fftw_iodim *dims = self.iodims;
		for (NSInteger i=self.dimensions.count-1; i >= 0; i--)
		{
			GLDimension *aDim = [self.dimensions objectAtIndex: i];
			if (i == self.dimensions.count-1 ) {
				dims[i].n = (int)aDim.nPoints;
				dims[i].is = 1;
				dims[i].os = 1;
			} else {
				dims[i].n = (int)aDim.nPoints;
				dims[i].is = dims[i+1].is*dims[i+1].n;
				dims[i].os = dims[i+1].is*dims[i+1].n;
			}
			
			self.dataPoints = self.dataPoints * aDim.nPoints;
		}
		
		self.zerosBuffer = [[GLMemoryPool sharedMemoryPool] dataWithLength: self.dataPoints*sizeof(GLFloat)];
		
	}
	return self;
}

// -->One *should* store the plans and use the new-array execute functions: http://www.fftw.org/doc/New_002darray-Execute-Functions.html
//
// The problem is that storing the plans requires that the separation between the real and imaginary
// memory buffers remain the same---but we're using a fixed zero buffer. So we're hosed.
//
// --> One *should* use real-to-complex transforms when possible and get a 2x speed up.
//
// But, this has the drawback of overwriting the input data, with no workaround on multidimensions.

- (void) transform: (GLFloat *) f forward: (GLSplitComplex *) fbar
{
	// Do an out-of-place FFT
	GLSplitComplex split;
	split.realp = f;
	split.imagp = self.zerosBuffer.mutableBytes;
	vGL_vclr( split.imagp, 1, self.dataPoints );
	
	vGL_fftw_plan plan = vGL_fftw_plan_guru_split_dft((int)self.dimensions.count, self.iodims, 0, NULL, split.realp, split.imagp, fbar->realp, fbar->imagp, FFTW_ESTIMATE);
	vGL_fftw_execute(plan);
	vGL_fftw_destroy_plan(plan);
		
	// Fix the forward transformation scaling factor
	// Yes, the *forward*, not the inverse. If we don't fix it here, then the fourier transform
	// is dependent on the number of discrete points, which isn't helpful.
	GLFloat inversesize = 1.0 / ( (GLFloat) self.dataPoints );
	vGL_vsmul( fbar->realp, 1, &inversesize, fbar->realp, 1, self.dataPoints );
	vGL_vsmul( fbar->imagp, 1, &inversesize, fbar->imagp, 1, self.dataPoints );
}

- (void) transform: (GLSplitComplex *) fbar inverse: (GLFloat *) f
{
	GLSplitComplex split;
	split.realp = f;
	split.imagp = self.zerosBuffer.mutableBytes;
	
	vGL_fftw_plan plan = vGL_fftw_plan_guru_split_dft((int)self.dimensions.count, self.iodims, 0, NULL, split.imagp, split.realp, fbar->imagp, fbar->realp, FFTW_ESTIMATE);
	vGL_fftw_execute(plan);
	vGL_fftw_destroy_plan(plan);	
}

- (void) dealloc
{
	[[GLMemoryPool sharedMemoryPool] returnData: self.zerosBuffer];
}


@end
