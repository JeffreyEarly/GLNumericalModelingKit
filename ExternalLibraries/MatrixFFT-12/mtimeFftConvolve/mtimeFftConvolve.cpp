/*	File: mtimeFftConvolve.cpp
	
	Description:
		Measure performance of MatrixFFT-based convolution.
	
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
 */
 
#include <stdlib.h>
#include <strings.h>
#include <stdio.h>
#include <unistd.h>
#include <libMatrixFFT/fftUtils.h>
#include <libMatrixFFT/fftConvolve.h>
#include <libMatrixFFT/MatrixFFT.h>
#include <libMatrixFFT/complexBufUtils.h>
#include <libMatrixFFT/fftUtils.h>
#include <Accelerate/Accelerate.h>

#define FC_IMAGE_ROWS_DEF	300		
#define FC_IMAGE_COLS_DEF	400

/* 
 * Kernel is a square of dimension k x k where k is odd.
 * If min != max we iterate for(k=min; k<=max; k+=incr).
 */
#define FC_KERNEL_SIZE_MIN	3	
#define FC_KERNEL_SIZE_MAX	17
#define FC_KERNEL_SIZE_INCR	2

static void usage(char **argv)
{
	printf("usage: %s [options]\n", argv[0]);
	printf("Options:\n");
	printf("  -r imageRows      -- image rows; default %u\n", FC_IMAGE_ROWS_DEF);
	printf("  -c imageCols      -- image columns; default %u\n", FC_IMAGE_COLS_DEF);
	printf("  -k kernelSizeMin  -- min kernel size, must be odd, default %u\n", 
		FC_KERNEL_SIZE_MIN);
	printf("  -K kernelSizeMax  -- max kernel size, must be odd, default %u\n", 
		FC_KERNEL_SIZE_MAX);
	printf("  -i kernelSizeIncr -- kernel size increment, must be even, default %u\n",
		FC_KERNEL_SIZE_INCR);
	printf("  -p                -- precalculate kernel FFT\n");
	printf("  -u                -- user time; default is wall time\n");
	printf("  -C                -- output in Crandall format; default is table\n");
	printf("  -t                -- trivial data for kernel and image\n");
	printf("  -n                -- no banner\n");
	printf("  -T numThreads     -- default is # of CPU cores\n");
	printf("  -d                -- dump buffers\n");
	exit(1);
}

/*
 * Parameters for one test run. All memory preallocated and written with 
 * test pattern before doTest() called. 
 * For clarity, the doTest() routine is given raw image and kernel in FFTFloat array
 * format; doTest() copies these to power-of-2-sized FFTComplex buffers as
 * would be done in the real world. These copies are not part of the timed operation.
 */
typedef struct {
	/* inputs */
	FFTFloat			*image;				/* raw image */
	FFTFloat			*kernel;			/* raw kernel */
	FFTComplex			*splitImage;		/* image, split format, 2^n x 2^m sized */
	FFTComplex			*splitKernel;		/* kernel, split format, 2^n x 2^m sized */
	FFTComplex			*result;

	size_t				kernelSize;			/* 'k', as in k x k */
	size_t				imageRows;	
	size_t				imageCols;	
	MatrixFFTPlan		mfftPlan;
	
	bool				precalcKernelFft;	/* when true, FFT(kernel) is not included in 
											 *   timing */
	bool				wallTime;
	bool				dumpBufs;			/* dump data for debug */
	/* output */
	double				runTime;			/* in seconds */
} FCParams;


/* 
 * Core test routine; performs one convolution. All data buffers allocated by caller. 
 * Returns nonzero on error. 
 */
static int doTest(FCParams *fcp)
{
	/*
	 * Calculate the size of the FFT. Each dimension is the power of 2 at least as large
	 * as (kernelSize + imageSize - 1). 
	 */
	unsigned log2FftRows;
	unsigned log2FftCols;
	
	size_t fftRows = fftRoundNumSamples(fcp->imageRows + fcp->kernelSize - 1, &log2FftRows);
	size_t fftCols = fftRoundNumSamples(fcp->imageCols + fcp->kernelSize - 1, &log2FftCols);

	/* 
	 * Convert image and kernel from FFTFloat arrays to split buffers. 
	 */
	fftConvCopyKernel(fcp->kernel, fcp->splitKernel, fcp->kernelSize, 
		log2FftRows, log2FftCols, true);
	fftConvCopyImage(fcp->image, fcp->splitImage, fcp->imageRows, fcp->imageCols, 
		log2FftRows, log2FftCols, true);
		
	/* 
	 * Optionally pre-calculate kernel FFT (in-place) outside of the timed section.
	 */
	MFFTReturn mrtn;
	if(fcp->precalcKernelFft) {
		mrtn = mfftExecute(fcp->mfftPlan, 0, true, fcp->splitKernel, fcp->splitKernel);
		if(mrtn) {
			mfftPrintErrInfo("mfftExecute", mrtn);
			return -1;
		}
	}

	/* Begin timed operation */
	double startTime = fftGetTime(fcp->wallTime);

	/* FFT(image), in-place */
	mrtn = mfftExecute(fcp->mfftPlan, 0, true, fcp->splitImage, fcp->splitImage);
	if(mrtn) {
		mfftPrintErrInfo("mfftExecute", mrtn);
		return -1;
	}

	/* optionally, FFT(kernel) */
	if(!fcp->precalcKernelFft) {
		mrtn = mfftExecute(fcp->mfftPlan, 0, true, fcp->splitKernel, fcp->splitKernel);
		if(mrtn) {
			mfftPrintErrInfo("mfftExecute", mrtn);
			return -1;
		}
	}

	/* Obtain format for the dyadic mul */
	MFFTFormat format;
	mrtn = mfftNativeFormat(fcp->mfftPlan, false, &format);
	if(mrtn) {
		return mrtn;
	}

	/* complex multiply of the two FFTs */
	fftConvDyadicMul(fcp->splitImage, fcp->splitKernel, fcp->result, log2FftRows, log2FftCols, format);
	
	if(fcp->dumpBufs) {
		fftDumpBuf("raw image", fcp->image, fcp->imageRows, fcp->imageCols);
		fftDumpBuf("raw kernel", fcp->kernel, fcp->kernelSize, fcp->kernelSize);
		fftDumpMatrixReal("FFT(image)", fcp->splitImage, fftRows, fftCols);
		fftDumpMatrixReal("FFT(kernel)", fcp->splitKernel, fftRows, fftCols);
		fftDumpMatrixReal("dyadic mul", fcp->result, fftRows, fftCols);
	}
	
	/* inverse FFT in place */
	mrtn = mfftExecute(fcp->mfftPlan, MEF_NormOutput, false, fcp->result, fcp->result);
	if(mrtn) {
		mfftPrintErrInfo("mfftExecute(inverse)", mrtn);
		return -1;
	}
	
	/*
	 * Scale once more. The MEF_NormOutput flag to the inverse FFT caused
	 * a scale by 1/2N, but we need to scale once more, by 0.5, since we
	 * did two forward FFTs (each of which produced an output 2x the actual
	 * values) and one inverse.
	 */
	fftScaleComplex(fcp->result, 0.5, (fftRows * fftCols) / 2);
	
	double endTime = fftGetTime(fcp->wallTime);

	if(fcp->dumpBufs) {
		fftDumpMatrixReal("result", fcp->result, fftRows, fftCols);
	}
	
	fcp->runTime = endTime - startTime;
	return 0;
}

/* 
 * Generate a string for exporting results to external tools for analysis.
 */
static unsigned crandCtr = 0;

static char *appendCrandallFormat(
	vImagePixelCount imageRows,
	vImagePixelCount imageCols,
	vImagePixelCount kernelSize,
	double runTime,
	char *crandallStr)				/* optional existing entry */
{
	/* nSqkSq = R * C * m^2 */
	unsigned long long rckSq = imageRows * imageCols;
	rckSq *= (unsigned long long)kernelSize;
	rckSq *= (unsigned long long)kernelSize;
	
	char newEntry[200];
	
	snprintf(newEntry, 200, "{%f, %llu}, ", runTime, rckSq);
	if((++crandCtr % 4) == 0) {
		strcat(newEntry, "\n");
	}
	
	unsigned newLen = strlen(newEntry) + 1;
	if(crandallStr) {
		newLen += strlen(crandallStr);
	}
	char *newStr = (char *)malloc(newLen);
	if(crandallStr) {
		strcpy(newStr, crandallStr);
		free(crandallStr);
	}
	else {
		newStr[0] = '\0';
	}
	strcat(newStr, newEntry);
	return newStr;
}

int main(int argc, char **argv)
{
	FCParams fcParams;
	memset(&fcParams, 0, sizeof(fcParams));
	fcParams.imageRows = FC_IMAGE_ROWS_DEF;
	fcParams.imageCols = FC_IMAGE_COLS_DEF;
	
	size_t kernelSizeMin = FC_KERNEL_SIZE_MIN;
	size_t kernelSizeMax = FC_KERNEL_SIZE_MAX;
	size_t kernelSizeIncr = FC_KERNEL_SIZE_INCR;
	bool crandallFormat = false;
	bool printBanner = true;
	bool trivialData = false;
	const char *optStr = NULL;
	unsigned numThreads = 0;
	fcParams.wallTime = true;
	
	int arg;
	while ((arg = getopt(argc, argv, "r:c:k:K:i:puCntdT:h")) != -1) {
		switch (arg) {
			case 'r':
				fcParams.imageRows = atoi(optarg);
				break;
			case 'c':
				fcParams.imageCols = atoi(optarg);
				break;
			case 'k':
				kernelSizeMin = atoi(optarg);
				break;
			case 'K':
				kernelSizeMax = atoi(optarg);
				break;
			case 'i':
				kernelSizeIncr = atoi(optarg);
				break;
			case 'p':
				fcParams.precalcKernelFft = true;
				optStr = "Precalculate kernel FFT";
				break;
			case 'u':
				fcParams.wallTime = false;
				break;
			case 'C':
				crandallFormat = true;
				break;
			case 'n':
				printBanner = false;
				break;
			case 'd':
				fcParams.dumpBufs = true;
				break;
			case 't':
				trivialData = true;
				break;
			case 'T':
				numThreads = atoi(optarg);
				break;
			case 'h':
				usage(argv);
		}
	}
	if(optind != argc) {
		usage(argv);
	}
	
	/* validate inputs */
	if(!(kernelSizeMin & 0x01) || !(kernelSizeMax & 0x01)) {
		printf("***kernelSize min and max (%lu, %lu) must be odd\n",
			(unsigned long)kernelSizeMin, (unsigned long)kernelSizeMax);
		exit(1);
	}
	if(kernelSizeIncr & 0x01) {
		printf("***kernelSizeIncr (%u) must be even\n", (unsigned)kernelSizeIncr);
		exit(1);
	}
	if(kernelSizeMax < kernelSizeMin) {
		printf("***kernelSizeMax must be >= kernelSizeMin\n");
		exit(1);
	}
	
	unsigned actThreads = numThreads;
	if(actThreads == 0) {
		/* Find out how many are actually going to be used */
		MFFTReturn mrtn = mfftNumThreads(NULL, &actThreads);
		if(mrtn) {
			mfftPrintErrInfo("mfftNumThreads", mrtn);
			exit(1);
		}
	}

	if(printBanner) {
		fftPrintTestBanner("Two-dimension Real Convolve", "MatrixFFT", 
			FFT_DOUBLE_PREC ? true : false,
			"Random", optStr, 0, actThreads);
		printf("\n");
		printf("Image Rows | Image Cols | Kernel size | Convolution time (s)\n");
		printf("-----------+------------+-------------+---------------------\n");
	}
	
	size_t imageSize = fcParams.imageRows * fcParams.imageCols;
	size_t maxKernelSize = kernelSizeMax * kernelSizeMax;
	
	/*
	 * Calculate the max size of the FFT. Each dimension is the power of 2 at least as large
	 * as (kernelSize + imageSize - 1). 
	 */
	unsigned maxFftRows;
	unsigned maxFftCols;
	unsigned log2FftRows;
	unsigned log2FftCols;
	maxFftRows = fftRoundNumSamples(fcParams.imageRows + kernelSizeMax - 1, &log2FftRows);
	maxFftCols = fftRoundNumSamples(fcParams.imageCols + kernelSizeMax - 1, &log2FftCols);

	/* to-be-freed pointers */
	FFTComplex *splitImageFree;
	FFTComplex *splitKernelFree;
	FFTComplex *resultFree;

	/* allocate buffers big enough for kernelSizeMax */
	fcParams.image       = (FFTFloat *)malloc(imageSize * sizeof(FFTFloat));
	fcParams.kernel      = (FFTFloat *)malloc(maxKernelSize * sizeof(FFTFloat));
	/* size of FFT buffers in complex elements */
	size_t fftSizeCompl  = maxFftRows * maxFftCols / 2;
	fcParams.splitImage  = fftAllocComplexArrayAlign(fftSizeCompl, FFT_MEM_ALIGNMENT, &splitImageFree);
	fcParams.splitKernel = fftAllocComplexArrayAlign(fftSizeCompl, FFT_MEM_ALIGNMENT, &splitKernelFree);
	fcParams.result      = fftAllocComplexArrayAlign(fftSizeCompl, FFT_MEM_ALIGNMENT, &resultFree);
	
	if(!fcParams.image        || !fcParams.kernel	   || 
		!fcParams.splitImage  || !fcParams.splitKernel || 
		!fcParams.result) {
		printf("***malloc failure (maxFftRows %lu maxFftCols %lu)\n",
			(unsigned long)maxFftRows, (unsigned long)maxFftCols);
		exit(1);
	}
	
	/* MatrixFFTPlan */
	unsigned n[2] = { log2FftRows, log2FftCols };
	MFFTReturn mrtn = mfftCreatePlan(2, n, true, 0, numThreads, &fcParams.mfftPlan);
	if(mrtn) {
		mfftPrintErrInfo("mfftCreatePlan", mrtn);
		exit(1);
	}

	if(trivialData) {
		/* 
		 * trival data for debug:
		 * -- image = incrementing 
		 * -- kernel = one 1.0 (result in an FFT of all 1.0, actually 2.0 pre-scale)
		 */
		genIncrFloat(fcParams.image, imageSize);
		memset(fcParams.kernel, 0, maxKernelSize * sizeof(FFTFloat));
		fcParams.kernel[0] = 1.0;
	}
	else {
		/* both image and kernel = random */
		genRandFloat(fcParams.image, imageSize);
		genRandFloat(fcParams.kernel, maxKernelSize);
	}
	
	/* here we go */
	char *crandallStr = NULL;

	if(crandallFormat) {
		crandallStr = strdup("\n/* {time(seconds), n^2 * m^2} */\n");
	}
	
	size_t kernelSize;
	for(kernelSize=kernelSizeMin; kernelSize<=kernelSizeMax; kernelSize+=kernelSizeIncr) {
		fcParams.kernelSize = kernelSize;
		int irtn = doTest(&fcParams);
		if(irtn) {
			printf("***test failure; aborting\n");
			exit(1);
		}
		
		printf( "%8lu   | %8lu   | %8lu    | %10.5f\n", 
			(unsigned long)fcParams.imageRows, (unsigned long)fcParams.imageCols, 
			(unsigned long)kernelSize,
			fcParams.runTime);
			
		if(crandallFormat) {
			/* save this result by appending it to the string we'll output when finished */
			crandallStr = appendCrandallFormat(fcParams.imageRows, fcParams.imageCols, 
				kernelSize, fcParams.runTime, crandallStr);
		}
	}
	
	if(crandallFormat) {
		printf("\n%s\n", crandallStr);
	}
	
	/* clean up */
	COND_FREE(fcParams.image);
	COND_FREE(fcParams.kernel);
	fftFreeComplexArrayAlign(fcParams.splitImage, splitImageFree);
	fftFreeComplexArrayAlign(fcParams.splitKernel, splitKernelFree);
	fftFreeComplexArrayAlign(fcParams.result, resultFree);
	mfftFreePlan(fcParams.mfftPlan);
	return 0;
}
