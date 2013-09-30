/*	File: mvfyConvolve.cpp
	
	Description:
		Verify proper operation of fftConvolve() and its supporting functions.
	
	Copyright:
		Copyright (C) 2008 Apple Inc.  All rights reserved.
	
	Disclaimer:
		IMPORTANT:  This Apple software is supplied to you by Apple
		Computer, Inc. ("Apple") in consideration of your agreement to
		the following terms, and your use, installation, modification
		or redistribution of this Apple software constitutes acceptance
		of these terms.  If you do not agree with these terms, please
		do not use, install, modify or redistribute this Apple
		software.

		In consideration of your agreement to abide by the following
		terms, and subject to these terms, Apple grants you a personal,
		non-exclusive license, under Apple’s copyrights in this
		original Apple software (the "Apple Software"), to use,
		reproduce, modify and redistribute the Apple Software, with or
		without modifications, in source and/or binary forms; provided
		that if you redistribute the Apple Software in its entirety and
		without modifications, you must retain this notice and the
		following text and disclaimers in all such redistributions of
		the Apple Software.  Neither the name, trademarks, service
		marks or logos of Apple Computer, Inc. may be used to endorse
		or promote products derived from the Apple Software without
		specific prior written permission from Apple.  Except as
		expressly stated in this notice, no other rights or licenses,
		express or implied, are granted by Apple herein, including but
		not limited to any patent rights that may be infringed by your
		derivative works or by other works in which the Apple Software
		may be incorporated.

		The Apple Software is provided by Apple on an "AS IS" basis.
		APPLE MAKES NO WARRANTIES, EXPRESS OR IMPLIED, INCLUDING
		WITHOUT LIMITATION THE IMPLIED WARRANTIES OF NON-INFRINGEMENT,
		MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, REGARDING
		THE APPLE SOFTWARE OR ITS USE AND OPERATION ALONE OR IN
		COMBINATION WITH YOUR PRODUCTS.

		IN NO EVENT SHALL APPLE BE LIABLE FOR ANY SPECIAL, INDIRECT,
		INCIDENTAL OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
		TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
		DATA, OR PROFITS; OR BUSINESS INTERRUPTION) ARISING IN ANY WAY
		OUT OF THE USE, REPRODUCTION, MODIFICATION AND/OR DISTRIBUTION
		OF THE APPLE SOFTWARE, HOWEVER CAUSED AND WHETHER UNDER THEORY
		OF CONTRACT, TORT (INCLUDING NEGLIGENCE), STRICT LIABILITY OR
		OTHERWISE, EVEN IF APPLE HAS BEEN ADVISED OF THE POSSIBILITY OF
		SUCH DAMAGE.
*/

/*
 * Created 12/23/2008. 
 * Copyright 2008 by Apple, Inc. 
 *
 * This verifier is a little awkward in that the vImage convolution routines,
 * which we use as a reference, only work on single precision floats. When
 * we are running in double precision (per the compile-time configuration in
 * MatrixFFTConfig.h), we have to do a bunch of shuffling and conversion of the 
 * original signals and of the vImage results.
 */
 
#include <stdlib.h>
#include <strings.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <Accelerate/Accelerate.h>
#include <CoreFoundation/CoreFoundation.h>
#include <libMatrixFFT/MatrixFFT.h>
#include <libMatrixFFT/fftConvolve.h>
#include <libMatrixFFT/complexBufUtils.h>
#include <libMatrixFFT/fftUtils.h>
#include <libMatrixFFT/devRandom.h>

#define KERNEL_SIZE_DEF		5
#define IMAGE_ROWS_DEF		8
#define IMAGE_COLS_DEF		4

static void usage(char **argv)
{
	printf("usage: %s [options]\n", argv[0]);
	printf("Options:\n");
	printf("  -k kernelSize    -- kernel size (per side); default %u\n", KERNEL_SIZE_DEF);
	printf("  -r imageRows     -- number of rows in image; default %u\n", IMAGE_ROWS_DEF);
	printf("  -c imageCols     -- number of columns in image; default %u\n", IMAGE_COLS_DEF);
	printf("  -C               -- test correlation; default is convolution\n");
	printf("  -i               -- incrementing image data; default is random\n");
	printf("  -v               -- custom trivial test kernel\n");
	printf("  -d               -- dump data\n");
	printf("  -b               -- no banner\n");
	printf("  -T numThreads    -- default is # of CPU cores\n");
	exit(1);
}

/*
 * For regression test purposes, here is a set of observed max delta and RMSE
 * values for a range of 'n'. If one of these is exceeded on a subsequent run, you
 * probably broke something because it used to work this well. 
 * These are very rough approximations - when using random data you get different
 * maxDelta and RMSE every run. 
 */
typedef struct {
	double maxDelta;
	double rmse;
} FFTErrors;

#define DISABLE_FFT_ERR		0
#if		DISABLE_FFT_ERR	
static FFTErrors fftErrors[1];
static unsigned fftNumErrors = 0;
#else	/* DISABLE_FFT_ERR */

#if FFT_DOUBLE_PREC

/* 
 * The index info this array is log2(total FFT size)
 */
static const FFTErrors fftErrors[] = 
{
	/* 0 */		{ 0.0,		0.0		},
	/* 1 */		{ 7.0e-07,  6.0e-07 },
	/* 2 */		{ 8.0e-07,	4.0e-07	},
	/* 3 */		{ 8.0e-07,	4.0e-07	},
	/* 4 */		{ 8.0e-07,	4.0e-07	},
	/* 5 */		{ 8.0e-07,	4.0e-07	},
	/* 6 */		{ 8.0e-07,	4.0e-07	},
	/* 7 */		{ 8.0e-07,	4.0e-07	},
	/* 8 */		{ 8.0e-07,	4.0e-07	},
	/* 9 */		{ 8.0e-07,	4.0e-07	},
	/* 10 */	{ 8.0e-07,	4.0e-07	},
	/* 11 */	{ 8.0e-07,	4.0e-07	},
	/* 12 */	{ 8.0e-07,	4.0e-07	},
	/* 13 */	{ 8.0e-07,	4.0e-07	},
	/* 14 */	{ 8.0e-07,	4.0e-07	},
	/* 15 */	{ 8.0e-07,	4.0e-07	},
	/* 16 */	{ 8.0e-07,	4.0e-07	},
	/* 17 */	{ 8.0e-07,	4.0e-07	},
	/* 18 */	{ 8.0e-07,	4.0e-07	},
	/* 19 */	{ 8.0e-07,	4.0e-07	},
	/* 20 */	{ 8.0e-07,	4.0e-07	},
	/* 21 */	{ 8.0e-07,	4.0e-07	},
	/* 22 */	{ 8.0e-07,	4.0e-07	},
	/* 23 */	{ 8.0e-07,	4.0e-07	},
	/* 24 */	{ 8.0e-07,	4.0e-07	},
	/* 25 */	{ 8.0e-07,	4.0e-07	},
};
#else	/* single precision */

static const FFTErrors fftErrors[] = 
{
	/* 0 */		{ 0.0,		0.0		},
	/* 1 */		{ 7.0e-07,  6.0e-07 },
	/* 2 */		{ 8.0e-07,	4.0e-07	},
	/* 3 */		{ 8.0e-07,	4.0e-07	},
	/* 4 */		{ 8.0e-07,	4.0e-07	},
	/* 5 */		{ 8.0e-07,	4.0e-07	},
	/* 6 */		{ 8.0e-07,	4.0e-07	},
	/* 7 */		{ 8.0e-07,	4.0e-07	},
	/* 8 */		{ 8.0e-07,	4.0e-07	},
	/* 9 */		{ 8.0e-07,	4.0e-07	},
	/* 10 */	{ 8.0e-07,	4.0e-07	},
	/* 11 */	{ 8.0e-07,	4.0e-07	},
	/* 12 */	{ 8.0e-07,	4.0e-07	},
	/* 13 */	{ 8.0e-07,	4.0e-07	},
	/* 14 */	{ 8.0e-07,	4.0e-07	},
	/* 15 */	{ 8.0e-07,	4.0e-07	},
	/* 16 */	{ 8.0e-07,	4.0e-07	},
	/* 17 */	{ 8.0e-07,	4.0e-07	},
	/* 18 */	{ 8.0e-07,	4.0e-07	},
	/* 19 */	{ 8.0e-07,	4.0e-07	},
	/* 20 */	{ 8.0e-07,	4.0e-07	},
	/* 21 */	{ 8.0e-07,	4.0e-07	},
	/* 22 */	{ 8.0e-07,	4.0e-07	},
	/* 23 */	{ 8.0e-07,	4.0e-07	},
	/* 24 */	{ 8.0e-07,	4.0e-07	},
	/* 25 */	{ 8.0e-07,	4.0e-07	},
};
#endif	/* FFT_DOUBLE_PREC */

static unsigned fftNumErrors = sizeof(fftErrors) / sizeof(fftErrors[0]);

#endif	/* DISABLE_FFT_ERR */

/*	
 * Fill an FFTFloat array with the real values of a Gaussian distribution:
 *
 * G(u,v) = (1 / 2 * pi * sigma^2) * e ^ -(u^2 + v^2) / (2 * sigma^2)
 *
 * See http://en.wikipedia.org/wiki/Gaussian_blur for more info.
 *
 * The results are normalized so that the sum of all of the G(u,v) values we
 * write is 1.0. Thus during the calculations we can ignore the 
 * (1 / 2 * pi * sigma^2) component since it's constant and would get
 * normalized out at the end. 
 *
 * Caller specifies both the size of a side of the buffer, which must be odd,
 * and the Gaussian radius. Normally a Gaussian kernel of radius r fits nicely 
 * in a square 6r+1 on a side, but the caller can specify any arbitrary sizes 
 * (though very small values like a bufSize < 4 will yield not-very-useful 
 * results).
 */
static void fillKernelGaussian(
	unsigned	radius,				/* a.k.a sigma */
	FFTFloat	*buf, 
	size_t		bufSize,			/* on a side */
	bool		dumpKernel)
{
	int			u, v;		
	FFTFloat	posSquared;			/* u^2 + v^2 */
	FFTFloat	sum = 0.0;			/* for accumulating normalization factor */
	FFTFloat	*fp = buf;
	
	/* denominator of expf() argument */
	FFTFloat expDenom = 2.0 * radius * radius;
	
	int maxPos = (int)(bufSize / 2);	/* trucated, it's odd */
	
	for(v=-maxPos; v<=maxPos; v++) {
		
		/* 
		 * Conceptually we're going from -maxPos <= y <= +maxPos, with y=-maxPos going to the 
		 * first output value.  
		 */
		
		/* invariant, of course, in the inner loop */
		FFTFloat vSquare = v * v;
		
		for(u=-maxPos; u<=maxPos; u++) {
			
			/* posSquared := u^2 + v^2 */
			posSquared = u*u + vSquare;
			
			FFTFloat exponent = -posSquared / expDenom;
			FFTFloat pointval = expf(exponent);
			
			sum += pointval;
			*fp++ = pointval;
		}
	}
	
	if(dumpKernel) {
		fftDumpBuf("Gaussian kernel", buf, bufSize, bufSize);
	}
	
	/* normalize */
	float normalizer = 1.0/sum;
	fp = buf;
	
	for (size_t r=0; r<bufSize; r++) {
		for (size_t c=0; c<bufSize; c++) {
			*fp++ *= normalizer;
		}
	}

	if(dumpKernel) {
		fftDumpBuf("Gaussian kernel, normalized", buf, bufSize, bufSize);
	}
}

/* 
 * Copy between native FFTFloat and float. 
 */
#if		FFT_DOUBLE_PREC
static void copyFFTFloatToFloat(
	const FFTFloat *ff, 
	float *f, 
	size_t numSamples)
{
	for(size_t dex=0; dex<numSamples; dex++) {
		*f++ = *ff++;
	}
}

#endif	/* FFT_DOUBLE_PREC */

/* 
 * For a real array in FFTComplex form, perform the "reverse" op in
 * preparation for a correlation:
 *
 *     y[row, col] swaps with y[(height-row) mod height, (width - col) mod width]
 *
 * This is slightly tricky because of the storage of the reals in 
 * FFTComplex form; odd column indices are in the imaginary 
 * component, and even column indices are in the real component. 
 */
static void reverseReals(
	FFTComplex *y,
	size_t fftRows,
	size_t fftCols)
{
	size_t complexCols = fftCols >> 1;
	
	/* the number of rows and columns we process - half size in each dim */
	size_t halfRows = fftRows >> 1;
	size_t halfCols = complexCols;		/* same number, different meaning */
	for(size_t row=0; row<halfRows+1; row++) {
		size_t swapRow = (row == 0) ? 0 : (fftRows - row);
		size_t numCols;
		if((row == 0) || (row == halfRows)) {
			/* these two rows reflect on themselves, avoid dup swaps... */
			numCols = halfCols;
		}
		else {
			/* other rows operate on all columns */
			numCols = fftCols;
		}
		for(size_t col=0; col<numCols; col++) {
			size_t swapCol		  = (col == 0) ? 0 : (fftCols - col);
			size_t colComplex	  = col     >> 1;
			size_t swapColComplex = swapCol >> 1;
			bool colOdd			  = (col & 0x01);
			bool swapColOdd		  = (swapCol & 0x01);
			
			/* pointers to the two floats we swap */
			FFTFloat *x1;
			FFTFloat *x2;
			
			/* offsets in complex of x1 and x2 */
			size_t x1Off = (row * complexCols) + colComplex;
			size_t x2Off = (swapRow * complexCols) + swapColComplex;
			
			#if FFT_SPLIT_COMPLEX
			x1 = colOdd     ? (y->imag + x1Off) : (y->real + x1Off);
			x2 = swapColOdd ? (y->imag + x2Off) : (y->real + x2Off);
			#else
			x1 = colOdd     ? &(y[x1Off].imag) : &(y[x1Off].real);
			x2 = swapColOdd ? &(y[x2Off].imag) : &(y[x2Off].real);
			#endif
			
			FFTFloat tmp = *x1;
			*x1 = *x2;
			*x2 = tmp;
		}
	}
}

#define TEST_KERNEL_SIZE	3

static void genTestKernel(
	FFTFloat *buf)
{
	*buf++ = 0.0; *buf++ = 0.0;  *buf++ = 0.0;
	*buf++ = 0.0; *buf++ = 1.00; *buf++ = 0.0;
	*buf++ = 0.0; *buf++ = 0.0;  *buf++ = 0.0;
}

int main(int argc, char **argv)
{
	size_t kernelSize = KERNEL_SIZE_DEF;
	size_t imageRows = IMAGE_ROWS_DEF;
	size_t imageCols = IMAGE_COLS_DEF;
	bool incrData = false;
	bool dumpData = false;
	bool printBanner = true;
	bool trivialKernel = false;
	bool logTiming = false;
	unsigned numThreads = 0;
	bool doCorrelate = false;
	
	int arg;
	while ((arg = getopt(argc, argv, "k:r:c:idbvtT:Ch")) != -1) {
		switch (arg) {
			case 'k':
				kernelSize = atoi(optarg);
				break;
			case 'r':
				imageRows = atoi(optarg);
				break;
			case 'c':
				imageCols = atoi(optarg);
				break;
			case 'i':
				incrData = true;
				break;
			case 'd':
				dumpData = true;
				break;
			case 'b':
				printBanner = false;
				break;
			case 't':
				logTiming = true;
				break;
			case 'v':
				trivialKernel = true;
				kernelSize = TEST_KERNEL_SIZE;
				break;
			case 'T':
				numThreads = atoi(optarg);
				break;
			case 'C':
				doCorrelate = true;
				break;
			case 'h':
			default:
				usage(argv);
		}
	}
	if(optind != argc) {
		usage(argv);
	}
	if((kernelSize & 0x01) == 0) {
		printf("***kernel size must be odd\n");
		exit(1);
	}
	
	/* 
	 * Min buffer size, before rounding up to powers of two, is
	 * (imageRows + kernelSize - 1) by (imageCols + kernelSize - 1)
	 */
	size_t fftRows;
	size_t fftCols;
	unsigned log2FftRows;
	unsigned log2FftCols;
	
	fftRows = fftRoundNumSamples(imageRows + kernelSize - 1, &log2FftRows);
	fftCols = fftRoundNumSamples(imageCols + kernelSize - 1, &log2FftCols);
	
	size_t totalFftSamples = fftRows * fftCols;
	size_t imageSamples = imageRows * imageCols;
	
	/*
	 * Alloc and init FFTComplexes, FFTFloat, and float arrays. 
	 *  -- one FFTFloat for image
	 *  -- one FFTFloat for kernel
	 *  -- one FFTComplex for split/padded image
	 *  -- one FFTComplex for split/padded kernel
	 *  -- one FFTComplex for fftConvolve() output 
	 *  -- one float for vImage image, if we're double precision
	 *  -- one float for vImage kernel, if we're double
	 *  -- one float for vImage output
	 */
	FFTFloat *image = NULL;
	FFTFloat *kernel = NULL;
	/* aligned buffers */
	FFTComplex *splitImage  = NULL;
	FFTComplex *splitKernel = NULL;
	FFTComplex *splitResult = NULL;
	/* to-be-freed buffers */
	FFTComplex *splitImageFree  = NULL;
	FFTComplex *splitKernelFree = NULL;
	FFTComplex *splitResultFree = NULL;
	
    float *vImageImage = NULL;
	float *vImageKernel = NULL;
	float *vImageResult = NULL;
	
	/* Raw image and kernel */
	image  = (FFTFloat *)malloc(sizeof(FFTFloat) * imageSamples);
	kernel = (FFTFloat *)malloc(sizeof(FFTFloat) * kernelSize * kernelSize);
	
	/* vImage stuff */
	size_t mallocSize   = sizeof(float) * imageSamples;
	vImageResult = (float *)malloc(mallocSize);
	#if		FFT_DOUBLE_PREC
	vImageImage  = (float *)malloc(mallocSize);
	vImageKernel = (float *)malloc(mallocSize);
	#else
	/* If single precision, vImage can use image, kernel directly */
	vImageImage  = image;
	vImageKernel = kernel;
	#endif	/* FFT_DOUBLE_PREC */
	
	/* split buffers */
	size_t numComplex = totalFftSamples / 2;
	splitImage  = fftAllocComplexArrayAlign(numComplex, FFT_MEM_ALIGNMENT, &splitImageFree);
	splitKernel = fftAllocComplexArrayAlign(numComplex, FFT_MEM_ALIGNMENT, &splitKernelFree);
	splitResult = fftAllocComplexArrayAlign(numComplex, FFT_MEM_ALIGNMENT, &splitResultFree);
	
	if(!image || !kernel || 
	   !vImageImage || !vImageKernel || !vImageResult ||
	   !splitImage  || !splitKernel  || !splitResult) {
		printf("***Malloc failure for numFftSamples = %lu\n", (unsigned long)totalFftSamples);
		exit(1);
	}
	
	genConstComplex(splitResult, numComplex, 0.0);
	
	/*
	 * Note: Throughout this test, we are attemtping to verify proper operation 
	 * of every function in fftConvolve.cpp, not to run really fast. Thus
	 * we place the image and kernel in FFTFloat arrays and use fftConvCopyKernel()
	 * and fftConvCopyImage() to reformat them. In the real world we'd probably
	 * place image and kernel directly in FFTComplex buffers. 
	 */
	if(incrData) {
		genIncrFloat(image, imageSamples);
	}
	else {
		genRandFloat(image, imageSamples);
	}
	
	/* 
	 * Figure out a reasonable Gaussian radius to fit in the specified 
	 * kernel size. Normally, a Gaussian kernel of radius r goes into
	 * a square 6r+1 samples on a side, so a kernel n samples on a side
	 * can fit a Gaussian of radius (n-1)/6. 
	 */
	unsigned radius = (kernelSize - 1) / 6;
	if(radius == 0) {
		radius = 1;
	}
	if(trivialKernel) {
		genTestKernel(kernel);
	}
	else {
		fillKernelGaussian(radius, kernel, kernelSize, dumpData);
	}
	
	/* 
	 * If we're built with double precision, we have to copy over
	 * image and kernel to the vImage float arrays.
	 */
	#if		FFT_DOUBLE_PREC
	copyFFTFloatToFloat(kernel, vImageKernel, kernelSize * kernelSize);
	copyFFTFloatToFloat(image, vImageImage, imageSamples);
	#endif	/* FFT_DOUBLE_PREC */
	
	if(dumpData) {
		fftDumpBuf("Original signal", image, imageRows, imageCols);
	}
	
	/* Set up a MatrixFFTPlan */
	MatrixFFTPlan mfftPlan = NULL; 
	MFFTReturn mrtn;
	unsigned n[2] = {log2FftRows, log2FftCols};
	mrtn = mfftCreatePlan(2, n, true, 0, numThreads, &mfftPlan);
	if(mrtn) {
		mfftPrintErrInfo("mfftCreatePlan", mrtn);
		exit(1);
	}
	
	/* 
	 * Convert image and kernel from FFTFloat arrays to split buffers. 
	 */
	fftConvCopyKernel(kernel, splitKernel, kernelSize, log2FftRows, log2FftCols, true);
	fftConvCopyImage(image,   splitImage,  imageRows,  imageCols,   log2FftRows, log2FftCols, true);
	
	if(dumpData) {
		fftDumpMatrixReal("Split Kernel", splitKernel, fftRows, fftCols);
		fftDumpMatrixReal("Split Image",  splitImage, fftRows, fftCols);
	}
	
	/* Convolve via vImage */
	vImage_Buffer imageBuf;
	vImage_Buffer dstBuf;
	
	imageBuf.data     = vImageImage;
	imageBuf.height   = imageRows;
	imageBuf.width    = imageCols;
	imageBuf.rowBytes = imageCols * sizeof(float);
	
	dstBuf.data       = vImageResult;
	dstBuf.height     = imageRows;
	dstBuf.width      = imageCols;
	dstBuf.rowBytes   = imageCols * sizeof(float);
	
	vImage_Error verr;
	verr = vImageConvolve_PlanarF(&imageBuf, &dstBuf, 
		NULL,		// temp buffer - forget it 
		0, 0,							// region of interest
		vImageKernel, kernelSize, kernelSize, 
		/* background fill 0.0 */
		0.0, kvImageBackgroundColorFill);
	if(verr) {
		printf("***vImageConvolve_PlanarF returned %ld\n", (long)verr);
		exit(1);
	}
	
	/* convolve or correlate via MatrixFFT */
	if(doCorrelate) {
		/*
		 * To test this, we assume that
		 *
		 * (x correlate y) == (x convolve (reverse y)))
		 *
		 * And, more useful here, 
		 *
		 * (x correlate (reverse y)) == (x convolve y)
		 *
		 * Since we're comparing against the vDSP output, namely (x convolve y),
		 * we use the second form above, which means we just do a swap of 
		 * y, which here is splitImage, and involve the correlation option 
		 * in fftConvolve().
		 *
		 * First, the reverse:
		 *
		 *     y[row, col] swaps with y[(height-row) mod height, (width - col) mod width]
		 */
		reverseReals(splitImage, fftRows, fftCols);
		if(dumpData) {
			fftDumpMatrixReal("split image reversed", splitImage, fftRows, fftCols);
		}

		/* now the conjugate-style convolution */
		mrtn = fftConvolve(mfftPlan, splitKernel, splitImage, 
			log2FftRows, log2FftCols,
			true,		// conjugate
			splitResult);
		if(mrtn) {
			mfftPrintErrInfo("fftConvolve", mrtn);
			exit(1);
		}
	}
	else {
		/* convolve via MatrixFFT */
		mrtn = fftConvolve(mfftPlan, splitKernel, splitImage, 
			log2FftRows, log2FftCols,
			false,		// conjugate
			splitResult);
		if(mrtn) {
			mfftPrintErrInfo("fftConvolve", mrtn);
			exit(1);
		}
	}
	
	/* trimmed result --> splitImage */
	fftConvRetrieveImage(splitResult, splitImage, log2FftRows, log2FftCols,
		imageRows, imageCols);

	size_t complexCols = (imageCols + 1) >> 1;
	if(dumpData) {
		fftDumpMatrixRect("fftConvolve result", splitImage, imageRows, complexCols);
	}

	/* for display and compare, convert vImage result to split format in splitKernel */
	fftCopyFloatToSplit(vImageResult, splitKernel, imageSamples);
	if(dumpData) {
		fftDumpMatrixRect("vImage result", splitKernel, imageRows, complexCols);
	}
	
	double maxDelta;
	double rmse;
	fftCompareRealAmplitudes(splitImage, splitKernel, imageSamples, 
		&maxDelta, &rmse);
	int ourRtn = 0;
	
	/* 
	 * Skip these checks for incrementing data, which leads to much larger
	 * maxDelta and RMSE
	 */
	unsigned log2FftSize = log2FftRows + log2FftCols;
	if(log2FftSize >= fftNumErrors) {
		printf("***Can't compare output, need bigger table\n");
		ourRtn = 1;
	}
	else if(!incrData) {
		static const FFTErrors *fftErr = &fftErrors[log2FftSize];
		if(maxDelta > fftErr->maxDelta) {
			printf("***maxDelta exceeded: expect %.3e, got %.3e\n",
				fftErr->maxDelta, maxDelta);
			ourRtn = 1;
		}
		if(rmse > fftErr->rmse) {
			printf("***rmseForward exceeded: expect %.3e, got %.3e\n",
				fftErr->rmse, rmse);
			ourRtn = 1;
		}
		
		/* test the test */
		if(rmse > maxDelta) {
			printf("*** Impossible: rmse > maxDelta\n");
			ourRtn = 1;
		}
	}
		
	if(printBanner) {
		fftPrintTestBanner(doCorrelate ? "Two-dimension Correlation" : "Two-dimension Convolve",
		    "MatrixFFT", 
			FFT_DOUBLE_PREC ? true : false, 
			incrData ? "Incrementing" : "Random", NULL, 0, numThreads);
		printf("\n");
		printf("  image R | image C | Kernel | FFT R | FFT C |  maxDelta |   RMSE   \n");
		printf(" ---------+---------+--------+-------+-------+-----------+----------\n");
	}
	printf("   %5u     %5u   %5u     %5u   %5u   %.3e   %.3e\n", 
		(unsigned)imageRows, (unsigned)imageCols, (unsigned)kernelSize,
		(unsigned)fftRows, (unsigned)fftCols,
		maxDelta, rmse);
			
	COND_FREE(image);
	COND_FREE(kernel);
	fftFreeComplexArrayAlign(splitImage, splitImageFree);
	fftFreeComplexArrayAlign(splitKernel, splitKernelFree);
	fftFreeComplexArrayAlign(splitResult, splitResultFree);
	
	#if		FFT_DOUBLE_PREC
	COND_FREE(vImageImage);
	COND_FREE(vImageKernel);
	COND_FREE(vImageResult);
	#endif
	
	mfftFreePlan(mfftPlan);
	return ourRtn;
}
