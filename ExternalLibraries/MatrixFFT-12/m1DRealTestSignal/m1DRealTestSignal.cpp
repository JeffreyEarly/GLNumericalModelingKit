/*	File: m1DRealTestSignal.cpp 
	
	Description:
		Verify proper operation of 1-D real FFT using a 1-D real test signal. 
	
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
#include <libMatrixFFT/real1DTestSignal.h>

#define MIN_SIZE_DEF		64
#define MAX_SIZE_DEF		(128 * 1024)

static void usage(char **argv)
{
	printf("usage: %s [options]\n", argv[0]);
	printf("Options:\n");
	printf("  -s minSize         -- minimum size; default is %u\n", MIN_SIZE_DEF);
	printf("  -S maxSize         -- maximum size; default is %u\n", MAX_SIZE_DEF);
	printf("  -T numThreads      -- default is # of CPU cores\n");
	printf("  -d                 -- dump data\n");
	printf("  -p                 -- pause for MallocDebug\n");
	printf("  -v                 -- verbose\n");
    printf("  -b                 -- no banner\n");
	exit(1);
}

/* Callback function for test signal functions */
static int callbackSeen = 0;

static void fftCallback(void *arg, unsigned percent)
{
	if(callbackSeen) {
		for(unsigned dex=0; dex<4; dex++) {
			putchar('\b'); 
		}
	}
	else {
		printf("Percent complete: ");
		callbackSeen = 1;
	}
	printf("%3d%%", percent);
	fflush(stdout);
}

#pragma mark --- run test on one size ---

typedef struct {
	unsigned		log2N;
	bool			dumpData;
	bool			doPause;
	unsigned		numThreads;				/* 0 here means default as in API */
	bool			verbose;			
} TestParams;

static int doTest(TestParams *tp)
{
	int ourRtn = 0;
	unsigned dimArray = tp->log2N;
	MatrixFFTPlan mfftPlan = NULL;
	MFFTReturn mrtn;

	/* 
	 * To minimize memory usage this always runs in place, and both forward and reverse
	 * FFTs are analyzed in-place as well.
	 */
	FFTComplex *fftBuf = NULL;	
	FFTComplex *fftBufFree = NULL;	
	size_t N = (size_t)1 << tp->log2N;				/* total real samples */
	size_t complexN = N >> 1;						/* total complex samples */
	
	uint32_t flags = 0;

	/* These nonzero values are markers for "we didn't measure" */
	double mxeForward  = 0.0;
	double rmseForward = 0.0;
	double mxeReverse  = 0.0;
	double rmseReverse = 0.0;
	
	if(tp->verbose) {
		printf("...creating MatrixFFTPlan...\n");
	}
	mrtn = mfftCreatePlan(1, &dimArray, true, flags, tp->numThreads, &mfftPlan);
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
		
	if(tp->verbose) {
		printf("...allocating data buffer...\n");
	}
	fftBuf = fftAllocComplexArrayAlign(complexN, FFT_MEM_ALIGNMENT, &fftBufFree);
	if(fftBuf == NULL) {
		printf("***Malloc failure for totalSamples = %llu\n", 
			(unsigned long long)complexN);
		ourRtn = -1;
		goto errOut;
	} 

	/* Generate test signal */
	if(tp->verbose) {
		printf("...generating test signal; ");
	}
	ourRtn = fftGenReal1DTestSignal(fftBuf, N, tp->verbose ? fftCallback : NULL, NULL);
	if(tp->verbose) {
		putchar('\n');
		callbackSeen = 0;
	}
	if(ourRtn) {
		printf("***Error generating test signal\n");
		goto errOut;
	}
	
	if(tp->dumpData) {
		fftDumpMatrixReal("Original", fftBuf, numRows, numCols);
	}
	
	if(tp->doPause) {
		printf("Pausing after setup, before FFTs: CR to continue: ");
		fflush(stdout);
		getchar();
	}

	/* forward FFT */
	if(tp->verbose) {
		printf("...performing forward FFT\n");
	}
	mrtn = mfftExecute(mfftPlan, 0, true, fftBuf, fftBuf);
	if(mrtn) {
		mfftPrintErrInfo("mfftExecute", mrtn);
		ourRtn = -1;
		goto errOut;
	}
	if(tp->dumpData) {
		fftDumpMatrixReal("Forward FFT output", fftBuf, numRows, numCols);
	}

	/* analyze forward FFT */
	if(tp->verbose) {
		printf("...analyzing forward FFT; ");
	}
	ourRtn = fftAnalyzeReal1DTestSignalFFT(mfftPlan, fftBuf, N, ST1D_Column, 
		tp->verbose ? fftCallback : NULL, NULL,
		&mxeForward, &rmseForward);
	if(tp->verbose) {
		putchar('\n');
		callbackSeen = 0;
	}
	if(ourRtn) {
		printf("***Error analyzing forward FFT\n");
		goto errOut;
	}
	
	/* inverse FFT - library normalizes */
	if(tp->verbose) {
		printf("...performing inverse FFT\n");
	}
	flags = MEF_NormOutput;
	mrtn = mfftExecute(mfftPlan, flags, false, fftBuf, fftBuf);
	if(mrtn) {
		mfftPrintErrInfo("mfftExecute", mrtn);
		ourRtn = -1;
		goto errOut;
	}

	if(tp->dumpData) {
		fftDumpMatrixReal("Inverse FFT output", fftBuf, numRows, numCols);
	}

	/* analyze round trip */
	if(tp->verbose) {
		printf("...analyzing inverse FFT; ");
	}
	ourRtn = fftAnalyzeReal1DTestSignal(fftBuf, N, tp->verbose ? fftCallback : NULL, NULL,
		&mxeReverse, &rmseReverse);
	if(tp->verbose) {
		putchar('\n');
		callbackSeen = 0;
	}
	if(ourRtn) {
		printf("***Error analyzing inverse FFT\n");
		goto errOut;
	}
	
	if(tp->log2N < 10) {
		printf("    2^%u ", tp->log2N);
	}
	else {
		printf("    2^%u", tp->log2N);
	}
	printf("        %.1e     %.1e     %.1e    %.1e\n", 
		rmseReverse, mxeReverse, rmseForward, mxeForward);

errOut:
	fftFreeComplexArrayAlign(fftBuf, fftBufFree);
	if(mfftPlan) {
		mfftFreePlan(mfftPlan);
	}
	return ourRtn;
}

#pragma mark --- main() ---

int main(int argc, char **argv)
{
	TestParams tp;
	size_t minSize = MIN_SIZE_DEF;
	size_t maxSize = MAX_SIZE_DEF;
	bool printBanner = true;
    
	memset(&tp, 0, sizeof(tp));
	tp.dumpData = false;
	
	int arg;
	while ((arg = getopt(argc, argv, "s:S:T:dphvb")) != -1) {
		switch (arg) {
			case 's':
				minSize = fftParseStringRep(optarg);
				break;
			case 'S':
				maxSize = fftParseStringRep(optarg);
				break;
			case 'd':
				tp.dumpData = true;
				break;
			case 'p':
				tp.doPause = true;
				break;
			case 'T':
				tp.numThreads = atoi(optarg);
				break;
			case 'v':
				tp.verbose = true;
				break;
            case 'b':
                printBanner = false;
                break;
			case 'h':
			default:
				usage(argv);
		}
	}
	if(optind != argc) {
		usage(argv);
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
	
	if(tp.doPause) {
		printf("Pausing at top; CR to continue: ");
		fflush(stdout);
		getchar();
	}

    if(printBanner) {
        fftPrintTestBanner("One-dimension real", "MatrixFFT", 
            FFT_DOUBLE_PREC ? true : false, 
            "1-D Real Test Signal", NULL, 0, actThreads);
        printf("\n");
        printf("   samples  | RMSE round | MXE round | RMSE fwd |  MXE fwd  \n");
        printf(" -----------+------------+-----------+----------+-----------\n");
    }
    
	if(!fftIsPowerOfTwo(minSize, &tp.log2N)) {
		printf("***size must be power of two\n");
		exit(1);
	}

	int ourRtn = 0;
	do {
		ourRtn = doTest(&tp);
		if(ourRtn) {
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
