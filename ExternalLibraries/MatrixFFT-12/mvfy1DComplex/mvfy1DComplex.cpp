/*	File: mvfy1DComplex.cpp 
	
	Description:
		Verify proper operation of 1-D complex FFT.
	
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
 * Created 12/17/2008. 
 * Copyright 2008 by Apple, Inc. 
 */
 
#include <stdlib.h>
#include <strings.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <Accelerate/Accelerate.h>
#include <CoreFoundation/CoreFoundation.h>
#include <libMatrixFFT/fftUtils.h>
#include <libMatrixFFT/complexBufUtils.h>
#include <libMatrixFFT/MatrixFFT.h>
#include <libMatrixFFT/vdspUtils.h>
#include <libMatrixFFT/src/fftPlatformConf.h>       /* private SPI */

#define MIN_SIZE_DEF		64
#define MAX_SIZE_DEF		(128 * 1024)

static void usage(char **argv)
{
	printf("usage: %s [options]\n", argv[0]);
	printf("Options:\n");
	printf("  -s minSize         -- minimum size; default is %u\n", MIN_SIZE_DEF);
	printf("  -S maxSize         -- maximum size; default is %u\n", MAX_SIZE_DEF);
	printf("  -o                 -- out-of-place; default is in-place\n");
	printf("  -m                 -- manual output transpose; default is lib\n");     
	printf("  -i                 -- incrementing data; default is random\n"); 
	printf("  -c                 -- chirp signal; default is random\n");
	printf("  -T numThreads      -- default is # of CPU cores\n");
	printf("  -N                 -- manual normalization (default is lib)\n");
	printf("  -O offset          -- signed rectangle offset; default is 0\n");
    printf("  -H                 -- set MCF_HintTranspose\n");
    printf("  -d on|off          -- force raw vDSP on/off (default is per config)\n");
	printf("  -D                 -- dump data\n");
	printf("  -u                 -- continue on error\n");
	printf("  -p                 -- pause for MallocDebug\n");
    printf("  -b                 -- no banner\n");
	exit(1);
}


/*
 * For regression test purposes, here is a set of observed max delta and RMSE
 * values for a range of 'n' where the total size of the arry is 2^n.
 * If one of these is exceeded on a subsequent run, you probably broke something 
 * because it used to work this well. 
 * These are very rough approximations - when using random data you get different
 * maxDelta and RMSE every run. 
 */
typedef struct {
	double mxe;
	double rmse;
} FFTErrors;

#pragma mark --- Allowed errors ---

#define DISABLE_FFT_ERR		0
#if		DISABLE_FFT_ERR	
static FFTErrors fftErrorsForward[1];
static FFTErrors roundTripErrors[1];
static unsigned fftNumErrors = 0;
#else	/* DISABLE_FFT_ERR */

#if FFT_DOUBLE_PREC

/* 
 * Forward errors. The index info this array is 'n'.
 */
static const FFTErrors fftErrorsForward[] = 
{
	/* 0 */		{ 0.0,		0.0		},
	/* 1 */		{ 3.0e-15,  1.0e-15 },
	/* 2 */		{ 1.0e-14,	1.0e-15	},
	/* 3 */		{ 9.0e-15,	3.0e-15 },
	/* 4 */		{ 2.0e-14,	5.0e-15 },
	/* 5 */		{ 6.0e-14,	2.0e-14 },
	/* 6 */		{ 3.0e-13,	2.5e-14 },
	/* 7 */		{ 9.0e-13,	5.8e-14 },
	/* 8 */		{ 4.0e-12,	1.3e-13 },
	/* 9 */		{ 2.0e-11,	3.1e-13 },
	/* 10 */	{ 7.0e-11,	7.0e-13 },
	/* 11 */	{ 3.0e-10,	1.4e-12 },
	/* 12 */	{ 7.3e-10,	2.0e-12 },
	/* 13 */	{ 7.3e-10,	2.0e-12 },
	/* 14 */	{ 7.3e-10,	2.0e-12 },
	/* 15 */	{ 7.3e-10,	2.0e-12 },
	/* 16 */	{ 7.3e-10,	3.0e-12 },
	/* 17 */	{ 7.3e-10,	4.0e-12 },
	/* 18 */	{ 7.3e-10,	6.0e-12 },
	/* 19 */	{ 8.2e-10,	1.0e-11 },
	/* 20 */	{ 4.0e-9,	2.0e-11 },
	/* 21 */	{ 8.0e-9,	3.0e-11 },
	/* 22 */	{ 2.0e-8,	8.0e-11 },
	/* 23 */	{ 3.0e-8, 	1.0e-10 },
	/* 24 */	{ 6.4e-8,	2.0e-10 },
	/* 25 */	{ 6.4e-8,	4.0e-10 },
	/* 26 */	{ 6.4e-8,	1.0e-9 },
	/* 27 */	{ 6.4e-8,	2.0e-9 },
	/* 28 */	{ 6.4e-8,	4.0e-9 },
	/* 29 */	{ 6.4e-8,	1.0e-8 },
	/* 30 */	{ 6.4e-8,	2.0e-8 },
	/* 31 */	{ 6.4e-8,	4.0e-8 },
	/* 32 */	{ 6.4e-8,	1.0e-7 },
};

/* Round trip errors */
static const FFTErrors roundTripErrors[] = 
{
	/* 0 */		{ 0.0,		0.0		},
	/* 1 */		{ 3.0e-15,  1.0e-15 },
	/* 2 */		{ 3.0e-15,  1.0e-15	},
	/* 3 */		{ 3.0e-15,  1.0e-15 },
	/* 4 */		{ 3.0e-15,  1.0e-15 },
	/* 5 */		{ 3.0e-15,  1.0e-15 },
	/* 6 */		{ 3.0e-15,  1.0e-15 },
	/* 7 */		{ 3.0e-15,  1.0e-15 },
	/* 8 */		{ 3.0e-15,  1.0e-15 },
	/* 9 */		{ 3.0e-15,  1.0e-15 },
	/* 10 */	{ 3.0e-15,  1.0e-15 },
	/* 11 */	{ 3.0e-15,  1.0e-15 },
	/* 12 */	{ 8.0e-15,  1.0e-15 },
	/* 13 */	{ 3.0e-14,  1.0e-14 },
	/* 14 */	{ 3.0e-14,  1.0e-14 },
	/* 15 */	{ 3.0e-14,  1.0e-14 },
	/* 16 */	{ 3.0e-14,  1.0e-14 },
	/* 17 */	{ 3.0e-14,  1.0e-14 },
	/* 18 */	{ 3.0e-14,  1.0e-14 },
	/* 19 */	{ 3.0e-14,  1.0e-14 },
	/* 20 */	{ 3.0e-13,  1.0e-13 },
	/* 21 */	{ 3.0e-13,  1.0e-13 },
	/* 22 */	{ 3.0e-13,  1.0e-13 },
	/* 23 */	{ 3.0e-13,  1.0e-13 },
	/* 24 */	{ 3.0e-13,  1.0e-13 },
	/* 25 */	{ 3.0e-13,  1.0e-13 },
	/* 26 */	{ 3.0e-13,  1.0e-13 },
	/* 27 */	{ 3.0e-13,  1.0e-13 },
	/* 28 */	{ 3.0e-13,  1.0e-13 },
	/* 29 */	{ 3.0e-13,  1.0e-13 },
	/* 30 */	{ 3.0e-13,  1.0e-13 },
	/* 31 */	{ 3.0e-13,  1.0e-13 },
	/* 32 */	{ 3.0e-13,  1.0e-13 },
};

#else	/* single precision */

/* forward */
static const FFTErrors fftErrorsForward[] = 
{
	/* 0 */		{ 0.0,		0.0		},
	/* 1 */		{ 7.0e-07,  6.0e-07 },
	/* 2 */		{ 8.0e-07,	4.0e-07	},
	/* 3 */		{ 6.0e-06,	2.0e-07 },
	/* 4 */		{ 2.0e-05,	2.0e-06 },
	/* 5 */		{ 4.0e-06,	2.0e-06 },
	/* 6 */		{ 4.0e-05,	1.0e-05 },
	/* 7 */		{ 2.0e-05,	1.0e-05 },
	/* 8 */		{ 3.0e-05,	1.0e-05 },
	/* 9 */		{ 6.0e-05,	1.0e-05 },
	/* 10 */	{ 8.0e-05,	5.0e-05 },
	/* 11 */	{ 9.0e-05,	5.0e-05 },
	/* 12 */	{ 4.0e-04,	6.0e-05 },
	/* 13 */	{ 4.0e-04,	6.0e-05 },
	/* 14 */	{ 1.0e-03,	1.0e-04 },
	/* 15 */	{ 3.0e-03,	1.0e-04 },
	/* 16 */	{ 2.0e-02,	8.0e-04 },
	/* 17 */	{ 4.0e-02,	2.0e-03 },
	/* 18 */	{ 1.0e-01,	4.0e-03 },
	/* 19 */	{ 2.0e-01,	4.0e-03 },
	/* 20 */	{ 3.0e-01,	8.0e-03 },
	/* 21 */	{ 4.0e-01,	8.0e-03 },
	/* 22 */	{ 6.0e-01,	4.0e-02 },
	/* 23 */	{ 7.0e-01,	4.0e-02 },
	/* 24 */	{ 9.0e-01,	6.0e-02 },	
	/* 25 */	{ 2.0e-00,	6.0e-02 },
	/* 26 */	{ 4.0e-00,	6.0e-02 },
	/* 27 */	{ 8.0e-00,	6.0e-02 },
	/* 28 */	{ 2.0e+01,	6.0e-02 },
	/* 29 */	{ 3.0e+01,	6.0e-02 },
	/* 30 */	{ 4.0e+01,	6.0e-02 },
	/* 31 */	{ 5.0e+01,	6.0e-02 },
	/* 32 */	{ 6.0e+01,	6.0e-02 },
};

/* round trip */
static const FFTErrors roundTripErrors[] = 
{
	/* 0 */		{ 0.0,		0.0		},
	/* 1 */		{ 5.0e-06,  5.0e-07 },
	/* 2 */		{ 5.0e-06,  5.0e-07	},
	/* 3 */		{ 5.0e-06,  5.0e-07},
	/* 4 */		{ 5.0e-06,  5.0e-07 },
	/* 5 */		{ 5.0e-06,  5.0e-07 },
	/* 6 */		{ 5.0e-06,  5.0e-07 },
	/* 7 */		{ 5.0e-06,  5.0e-07 },
	/* 8 */		{ 5.0e-06,  5.0e-07 },
	/* 9 */		{ 5.0e-06,  5.0e-07 },
	/* 10 */	{ 5.0e-06,  5.0e-07 },
	/* 11 */	{ 5.0e-06,  5.0e-07 },
	/* 12 */	{ 5.0e-06,  5.0e-07 },
	/* 13 */	{ 5.0e-06,  5.0e-07 },
	/* 14 */	{ 5.0e-06,  5.0e-07 },
	/* 15 */	{ 5.0e-06,  8.0e-07 },
	/* 16 */	{ 5.0e-06,  8.0e-07 },
	/* 17 */	{ 1.0e-05,  1.0e-06 },
	/* 18 */	{ 1.0e-05,  1.0e-06 },
	/* 19 */	{ 1.0e-05,  1.0e-06 },
	/* 20 */	{ 1.0e-05,  1.0e-06 },
	/* 21 */	{ 1.0e-05,  1.0e-06 },
	/* 22 */	{ 1.0e-05,  1.0e-06 },
	/* 23 */	{ 1.0e-05,  1.0e-06 },
	/* 24 */	{ 1.0e-05,  5.0e-06 },	
	/* 25 */	{ 5.0e-05,  5.0e-06 },
	/* 26 */	{ 5.0e-05,  5.0e-06 },
	/* 27 */	{ 5.0e-05,  5.0e-06 },
	/* 28 */	{ 5.0e-05,  5.0e-06 },
	/* 29 */	{ 5.0e-05,  5.0e-06 },
	/* 30 */	{ 5.0e-05,  5.0e-06 },
	/* 31 */	{ 5.0e-05,  5.0e-06 },
	/* 32 */	{ 5.0e-05,  5.0e-06 },
};

#endif	/* FFT_DOUBLE_PREC */

static unsigned fftNumErrors = sizeof(fftErrorsForward) / sizeof(fftErrorsForward[0]);

#endif	/* DISABLE_FFT_ERR */

typedef enum {
	ID_Random,
	ID_Increment,
	ID_Chirp
} InputData;

#pragma mark --- run test on one size ---

typedef struct {
	unsigned		log2N;
	bool			inPlace;
	InputData		inputData;
	bool			dumpData;
	bool			manualTranspose;
    bool            hintTranspose;
	bool			doPause;
	unsigned		numThreads;				/* 0 here means default as in API */
	bool			manualNorm;
} TestParams;

static int doTest(TestParams *tp)
{
	int ourRtn = 0;
	unsigned dimArray = tp->log2N;
	MatrixFFTPlan mfftPlan = NULL;
	MFFTReturn mrtn;

	vDSPComplex zBufIn = {NULL, NULL};				/* reference signal */
	vDSPComplex zBufOut = {NULL, NULL};				/* vDSP output */
	FFTComplex *fftSrc = NULL;						/* aligned MatrixFFT input */
	FFTComplex *fftDst = NULL;						/* aligned MatrixFFT output if OOP */
	FFTComplex *fftSrcFree = NULL;					/* to-be-freed MatrixFFT input */
	FFTComplex *fftDstFree = NULL;					/* to-be-freed MatrixFFT output if OOP */
	FFTComplex *actFftDst = NULL;					/* either fftSrc or fftDst */
	
	FFT_Setup dspSetup = NULL;

	size_t N = (size_t)1 << tp->log2N;
	
	FFTFloat normFact = 1.0 / (FFTFloat)N;	
	
	size_t zeroSize;
	uint32_t flags = 0;

    if(tp->hintTranspose) {
        flags |= MCF_HintTranspose;
    }
    
	/* These nonzero values are markers for "we didn't measure" */
	double maxDeltaForward = 7.7e7;
	double rmseForward = 7.7e7;
	double maxDeltaReverse;
	double rmseReverse;
	
	mrtn = mfftCreatePlan(1, &dimArray, false, flags, tp->numThreads, &mfftPlan);
	if(mrtn) {
		mfftPrintErrInfo("mfftCreatePlan", mrtn);
		return -1;
	}

	/* These are for displaying data */
	size_t numRows;
	size_t numCols;
	mrtn = mfftRectangle(mfftPlan, &numRows, &numCols);
	if(mrtn) {
		mfftPrintErrInfo("mfftRectangle", mrtn);
		return -1;
	}
			
	dspSetup = FFTCreateSetup(tp->log2N);
	if(dspSetup == NULL) {
		printf("***malloc failure on FFTCreateSetup()\n");
		exit(1);
	}

	/* Alloc & init buffers */
	if(fftAllocDSPComplex(&zBufIn, N)) {
		printf("***Malloc failure for numSamples = %llu\n", 
			(unsigned long long)N);
		return 1;
	}
	if(fftAllocDSPComplex(&zBufOut, N)) {
		printf("***Malloc failure for numSamples = %llu\n", 
			(unsigned long long)N);
		ourRtn = -1;
		goto errOut;
	}
	fftSrc = fftAllocComplexArrayAlign(N, FFT_MEM_ALIGNMENT, &fftSrcFree);
	if(fftSrc == NULL) {
		printf("***Malloc failure for totalSamples = %llu\n", 
			(unsigned long long)N);
		ourRtn = -1;
		goto errOut;
	} 
	if(tp->inPlace) {
		actFftDst = fftSrc;
	}
	else {
		fftDst = fftAllocComplexArrayAlign(N, FFT_MEM_ALIGNMENT, &fftDstFree);
		if(fftDst == NULL) {
			printf("***Malloc failure for totalSamples = %llu\n", 
				(unsigned long long)N);
			ourRtn = -1;
			goto errOut;
		} 
		actFftDst = fftDst;
	}
	
	switch(tp->inputData) {
		case ID_Increment:
			genIncrComplexDSP(&zBufIn, N);
			break;
		case ID_Random:
			genRandComplexDSP(&zBufIn, N);
			break;
		case ID_Chirp:
			fftGenChirp1DDSP(&zBufIn, N);
			break;
	}
	fftCopyFromDSP(&zBufIn, fftSrc, N);
	
	zeroSize = N * sizeof(FFTFloat);
	memset(zBufOut.realp, 0, zeroSize);
	memset(zBufOut.imagp, 0, zeroSize);
	if(!tp->inPlace) {
		genConstComplex(fftDst, N, 0.0);
	}
	
	if(tp->dumpData) {
		fftDump2DDSPComplex("Original", &zBufIn, numRows, numCols);
	}
	
	if(tp->doPause) {
		printf("Pausing after setup, before FFTs: CR to continue: ");
		fflush(stdout);
		getchar();
	}

	/* FFT via vDSP */
	FFTComplex1dOP(dspSetup, &zBufIn, &zBufOut, tp->log2N, FFT_FORWARD);
	if(tp->dumpData) {
		fftDump2DDSPComplex("vDSP", &zBufOut, numRows, numCols);
		fftDump2DDSPComplex("vDSP", &zBufOut, numCols, numRows);
	}

	/* FFT by MatrixFFT */
	flags = 0;
	if(!tp->manualTranspose) {
		flags |= MEF_TransposeOutput;
	}
	mrtn = mfftExecute(mfftPlan, flags, true, fftSrc, actFftDst);
	if(mrtn) {
		mfftPrintErrInfo("mfftExecute", mrtn);
		ourRtn = -1;
		goto errOut;
	}
	if(tp->manualTranspose) {
		if(tp->dumpData) {
			/* it's in column order here */
			fftDump2DComplex("mfft pre-transpose", actFftDst, numRows, numCols);
			fftDump2DComplex("mfft pre-transpose", actFftDst, numCols, numRows);
		}
		mrtn = mfftTranspose(mfftPlan, true, false, actFftDst, actFftDst);
		if(mrtn) {
			mfftPrintErrInfo("mfftTranspose", mrtn);
			ourRtn = -1;
			goto errOut;
		}
	}
	if(tp->dumpData) {
		fftDump2DComplex("MatrixFFT forward output", actFftDst, numRows, numCols);
	}

	fftCompareAmplitudesWithDSP(actFftDst, &zBufOut, N, 
			&maxDeltaForward, &rmseForward);
	
	/* 
	 * Skip these checks for incrementing data, which leads to much larger
	 * maxDelta and RMSE.
	 * Also skip for outTransType == OTT_None, where we know we don't match the 
	 * standard FFT Output. 
	 */
	if(tp->log2N >= fftNumErrors) {
		printf("***Can't compare output, need bigger table\n");
		ourRtn = 1;
	}
	else if(tp->inputData != ID_Increment) {
		const FFTErrors *fftErr = &fftErrorsForward[tp->log2N];
		if(maxDeltaForward > fftErr->mxe) {
			printf("***MXE Forward exceeded: expect %.3e, got %.3e\n",
				fftErr->mxe, maxDeltaForward);			
			ourRtn = 1;
		}
		if(rmseForward > fftErr->rmse) {
			printf("***rmse Forward exceeded: expect %.3e, got %.3e\n",
				fftErr->rmse, rmseForward);
			ourRtn = 1;
		}
		
		/* test the test */
		if(rmseForward > maxDeltaForward) {
			printf("*** Impossible: rmseForward > maxDeltaForward\n");
			ourRtn = 1;
		}
	}
	
	/* inverse FFT */
	if(tp->manualTranspose) {
		mrtn = mfftTranspose(mfftPlan, false, true, actFftDst, actFftDst);
		if(mrtn) {
			mfftPrintErrInfo("mfftTranspose", mrtn);
			ourRtn = -1;
			goto errOut;
		}
		flags = 0;
	}
	else {
		flags = MEF_TransposeInput;
	}
	if(!tp->manualNorm) {
		flags |= MEF_NormOutput;
	}

	mrtn = mfftExecute(mfftPlan, flags, false, actFftDst, fftSrc);
	if(mrtn) {
		mfftPrintErrInfo("mfftExecute", mrtn);
		ourRtn = -1;
		goto errOut;
	}

	if(tp->dumpData) {
		fftDump2DComplex("MatrixFFT reverse", fftSrc, numRows, numCols);
	}

	if(tp->manualNorm) {
		/* 
		 * scale by 1/N
		 */
		fftScaleComplex(fftSrc, normFact, N);
		if(tp->dumpData) {
			fftDump2DComplex("post scale", fftSrc, numRows, numCols);
		}
	}
	
	fftCompareAmplitudesWithDSP(fftSrc, &zBufIn, N, 
		&maxDeltaReverse, &rmseReverse);
	
	/* 
	 * Skip these checks for incrementing data, which leads to much larger
	 * maxDelta and RMSE
	 */
	if((tp->inputData != ID_Increment) && (tp->log2N < fftNumErrors)) {
		const FFTErrors *fftErr = &roundTripErrors[tp->log2N];
		if(maxDeltaReverse > fftErr->mxe) {
			printf("***MXE Reverse exceeded: expect %.3e, got %.3e\n",
				fftErr->mxe, maxDeltaReverse);
			ourRtn = 1;
		}
		if(rmseReverse > fftErr->rmse) {
			printf("***rmse Reverse exceeded: expect %.3e, got %.3e\n",
				fftErr->rmse, rmseReverse);
			ourRtn = 1;
		}
		
		/* test the test */
		if(rmseReverse > maxDeltaReverse) {
			printf("*** Impossible: rmseReverse > maxDeltaReverse\n");
			ourRtn = 1;
		}
	}
	
	if(tp->log2N < 10) {
		printf("    2^%u ", tp->log2N);
	}
	else {
		printf("    2^%u", tp->log2N);
	}
	printf("      %.1e    %.1e    %.1e     %.1e\n", 
		maxDeltaForward, rmseForward, maxDeltaReverse, rmseReverse);

errOut:
	fftFreeDSPComplex(&zBufIn);
	fftFreeDSPComplex(&zBufOut);
	fftFreeComplexArrayAlign(fftSrc, fftSrcFree);
	if(!tp->inPlace) {
		fftFreeComplexArrayAlign(fftDst, fftDstFree);
	}
	if(mfftPlan) {
		mfftFreePlan(mfftPlan);
	}
	if(dspSetup) {
		FFTFreeSetup(dspSetup);
	}
	return ourRtn;
}

#pragma mark --- main() ---

int main(int argc, char **argv)
{
	char optStr[200];
	TestParams tp;
	size_t minSize = MIN_SIZE_DEF;
	size_t maxSize = MAX_SIZE_DEF;
	bool continueOnErr = false;
	int rectOffset = 0;
    bool rectOffsetSpecd = false;
	const char *dataFormStr = "Random";
	bool printBanner = true;
	MFFT_ForceVdsp forceVdsp = MF_Default;
    
	optStr[0] = '\0';
	memset(&tp, 0, sizeof(tp));
	tp.inPlace = true;
	tp.manualTranspose = false;
	tp.inputData = ID_Random;
	tp.dumpData = false;
	
	int arg;
	while ((arg = getopt(argc, argv, "s:S:omicDpCT:O:NHbd:vh")) != -1) {
		switch (arg) {
			case 's':
				minSize = fftParseStringRep(optarg);
				break;
			case 'S':
				maxSize = fftParseStringRep(optarg);
				break;
			case 'o':
				tp.inPlace = false;
				break;
			case 'm':
				tp.manualTranspose = true;
				appendOptStr(optStr, "Manual transpose");
				break;
			case 'i':
				tp.inputData = ID_Increment;
				dataFormStr = "Incrementing";
				break;
			case 'c':
				tp.inputData = ID_Chirp;
				dataFormStr = "1-D Chirp";
				break;
			case 'D':
				tp.dumpData = true;
				break;
			case 'p':
				tp.doPause = true;
				break;
			case 'u':
				continueOnErr = true;
				break;
			case 'T':
				tp.numThreads = atoi(optarg);
				break;
			case 'O':
				rectOffset = atoi(optarg);
                rectOffsetSpecd = true;
				break;
			case 'N':
				tp.manualNorm = true;
				appendOptStr(optStr, "Manual normalize");
				break;
            case 'H':
                tp.hintTranspose = true;
				appendOptStr(optStr, "HintTranspose");
                break;
            case 'b':
                printBanner = false;
                break;
            case 'd':
                if(!strcmp(optarg, "on")) {
                    forceVdsp = MF_ForceVdsp;
                }
                else if(!strcmp(optarg, "off")) {
                    forceVdsp = MF_ForceMfft;
                }
                else {
                    printf("***Bad forceVdsp value\n");
                    usage(argv);
                }
                break;
			case 'h':
			default:
				usage(argv);
		}
	}
	if(optind != argc) {
		usage(argv);
	}
	if(!tp.inPlace) {
		appendOptStr(optStr, "Out-of-place");
	}
	else {
		appendOptStr(optStr, "In-place");
	}

	unsigned actThreads = tp.numThreads;
	if(actThreads == 0) {
		/* Find out how many are actually going to be used */
		MFFTReturn mrtn = mfftNumThreads(NULL, &actThreads);
		if(mrtn) {
			mfftPrintErrInfo("mfftNumThreads", mrtn);
			exit(1);
		}
	}
	if(rectOffsetSpecd) {
		mfftSetRectangleOffset(rectOffset, true);
		char str[100];
		sprintf(str, "rectOffset=%d", rectOffset);
		appendOptStr(optStr, str);
	}
	if(forceVdsp != MF_Default) {
        mfftSetForceVdsp(forceVdsp);
		appendOptStr(optStr, (forceVdsp == MF_ForceVdsp) ? "Force vDSP" : "Force !vDSP");
    }
	if(tp.doPause) {
		printf("Pausing at top; CR to continue: ");
		fflush(stdout);
		getchar();
	}

    if(printBanner) {
        fftPrintTestBanner("One-dimension complex", "MatrixFFT", 
            FFT_DOUBLE_PREC ? true : false, 
            dataFormStr, optStr, 0, actThreads);
        printf("\n");
        printf("   samples  | MXE fwd | RMSE fwd | MXE round | RMSE round\n");
        printf(" -----------+---------+----------+-----------+-----------\n");
    }
    
	if(!fftIsPowerOfTwo(minSize, &tp.log2N)) {
		printf("***size must be power of two\n");
		exit(1);
	}

	int ourRtn = 0;
	do {
		ourRtn = doTest(&tp);
		if(ourRtn && !continueOnErr) {
			break;
		}
		if(tp.doPause) {
			printf("Pausing at end of loop; CR to continue: ");
			fflush(stdout);
			getchar();
		}
		tp.log2N++;
	} while(((size_t)1 << tp.log2N) <= maxSize);

	if(tp.doPause) {
		printf("Pausing after test, before frees for malloc debug; CR to continue: ");
		fflush(stdout);
		getchar();
	}
	
	/* clean up here */
	return ourRtn;
}
