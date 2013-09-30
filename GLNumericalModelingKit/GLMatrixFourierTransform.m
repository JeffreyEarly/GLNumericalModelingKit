//
//  GLMatrixFourierTransform.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 10/11/11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#import "GLMatrixFourierTransform.h"

#import "GLDimension.h"
#import "GLMemoryPool.h"

#include <libMatrixFFT/fftUtils.h>
#include <libMatrixFFT/MatrixFFT.h>
#include <libMatrixFFT/vdspUtils.h>
#include <libMatrixFFT/complexBufUtils.h>

@interface GLMatrixFourierTransform ()
@property(readwrite, assign, nonatomic) MatrixFFTPlan mfftPlan;
@end

@implementation GLMatrixFourierTransform

@synthesize dimensions;
@synthesize zerosBuffer;
@synthesize mfftPlan;
@synthesize dataPoints;

- (GLFourierTransform *) initWithDimensions: (NSArray *) theDimensions
{
	if ( theDimensions.count > 2 ) {
		[NSException raise: @"NumberOfDimensionsUnsupportedException" format: @"MatrixFFT only supports up to two dimensions."];
		return nil;
	}
	
	if ((self=[super init])) {
		dimensions = theDimensions;
		
		unsigned int dims[2];
		dataPoints = 1;
		unsigned int i=0;
		for ( GLDimension *aDim in self.dimensions ) {
			dataPoints = dataPoints * aDim.nPoints;
			vDSP_Length logDim = 31 - __builtin_clz( (unsigned int) aDim.nPoints );
			dims[i] = (unsigned int) logDim;
			i++;
		}

		// Create the plan
		mfftCreatePlan(i, dims, false, 0, 0, &mfftPlan);
		
		zerosBuffer = [[GLMemoryPool sharedMemoryPool] dataWithLength: self.dataPoints*sizeof(GLFloat)];
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
	
	mfftExecute(mfftPlan, MEF_TransposeOutput, true, (FFTComplex *) &split, (FFTComplex *) fbar);
	
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
	
	mfftExecute(mfftPlan, MEF_TransposeInput, false, (FFTComplex *) fbar, (FFTComplex *) &split);
}

- (void) dealloc
{
	[[GLMemoryPool sharedMemoryPool] returnData: self.zerosBuffer];
}

@end
