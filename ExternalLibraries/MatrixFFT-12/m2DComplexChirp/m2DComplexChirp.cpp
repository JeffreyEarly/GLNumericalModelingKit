/*	File: m2DComplexChirp.cpp 
	
	Description:
		Verify proper operation of 2-D complex FFT using chirp signal. 
	
	Copyright:
		Copyright (C) 2009 Apple Inc.  All rights reserved.
	
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
 * Created 01/12/2009. 
 * Copyright 2009 by Apple, Inc. 
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
#include <libMatrixFFT/complexChirpSignal.h>

#define MIN_ROW_SIZE_DEF		32
#define MAX_ROW_SIZE_DEF		(4 * 1024)

static void usage(char **argv)
{
	printf("usage: %s [options]\n", argv[0]);
	printf("Options:\n");
	printf("  -s minRowSize      -- minimum row size; default is %u\n", MIN_ROW_SIZE_DEF);
	printf("  -S maxRowSize      -- maximum row size; default is %u\n", MAX_ROW_SIZE_DEF);
    printf("  -c colSize         -- one run, minRowSize x colSize\n");
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
	unsigned		log2NumRows;
	unsigned		log2NumCols;
	bool			dumpData;
	bool			doPause;
	unsigned		numThreads;				/* 0 here means default as in API */
	bool			verbose;			
} TestParams;

static int doTest(TestParams *tp)
{
	int ourRtn = 0;
	unsigned dimArray[2] = {tp->log2NumRows, tp->log2NumCols};
	MatrixFFTPlan mfftPlan = NULL;
	MFFTReturn mrtn;
	char rowStr[5];
	char colStr[5];
	unsigned log2Total = tp->log2NumRows + tp->log2NumCols;

	/* 
	 * To minimize memory usage this always runs in place, and both forward and reverse
	 * FFTs are analyzed in-place as well.
	 */
	FFTComplex *fftBuf = NULL;	
	size_t numRows = (size_t)(1) << tp->log2NumRows;
	size_t numCols = (size_t)(1) << tp->log2NumCols;
	size_t numSamples = numRows * numCols;			/* total complex samples */
	
	uint32_t flags = 0;

	double che  = 0.0;
	double mxeReverse  = 0.0;
	double rmseReverse = 0.0;
	
	if(tp->verbose) {
		printf("...creating MatrixFFTPlan...\n");
	}
	mrtn = mfftCreatePlan(2, dimArray, false, flags, tp->numThreads, &mfftPlan);
	if(mrtn) {
		mfftPrintErrInfo("mfftCreatePlan", mrtn);
		return -1;
	}

	if(tp->verbose) {
		printf("...allocating data buffer...\n");
	}
	fftBuf = fftAllocComplexArray(numSamples);
	if(fftBuf == NULL) {
		printf("***Malloc failure for numSamples = %llu\n", 
			(unsigned long long)numSamples);
		ourRtn = -1;
		goto errOut;
	} 

	/* Generate test signal */
	if(tp->verbose) {
		printf("...generating test signal; ");
	}
	fftGenChirp2D(fftBuf, numCols, numRows, tp->verbose ? fftCallback : NULL, NULL);
	if(tp->verbose) {
		putchar('\n');
		callbackSeen = 0;
	}
	
	if(tp->dumpData) {
		fftDumpMatrixRect("Original", fftBuf, numRows, numCols);
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
		fftDumpMatrixRect("Forward FFT output", fftBuf, numRows, numCols);
	}

	/* analyze forward FFT */
	if(tp->verbose) {
		printf("...analyzing forward FFT; ");
	}
	che = fftMaxChirpError(fftBuf, numSamples, 
		tp->verbose ? fftCallback : NULL, NULL);
	if(tp->verbose) {
		putchar('\n');
		callbackSeen = 0;
	}
	
	/* inverse FFT - have the library normalize */
	flags = MEF_NormOutput;
	if(tp->verbose) {
		printf("...performing inverse FFT\n");
	}
	mrtn = mfftExecute(mfftPlan, flags, false, fftBuf, fftBuf);
	if(mrtn) {
		mfftPrintErrInfo("mfftExecute", mrtn);
		ourRtn = -1;
		goto errOut;
	}

	if(tp->dumpData) {
		fftDumpMatrixRect("Inverse FFT output", fftBuf, numRows, numCols);
	}

	/* analyze round trip */
	if(tp->verbose) {
		printf("...analyzing inverse FFT; ");
	}
	fftAnalyzeChirp2D(fftBuf, numCols, numRows, 
		tp->verbose ? fftCallback : NULL, NULL,
		&mxeReverse, &rmseReverse);
	if(tp->verbose) {
		putchar('\n');
		callbackSeen = 0;
	}
	
	sprintf(rowStr, "2^%u", tp->log2NumRows);
	sprintf(colStr, "2^%u", tp->log2NumCols);
	printf("  %-4s   %-4s ", rowStr, colStr);
	if(log2Total < 10) {
		printf("     2^%u", log2Total);
	}
	else {
		printf("    2^%u", log2Total);
	}
	printf("     %.1e      %.1e    %.1e\n", 
		rmseReverse, mxeReverse, che);

errOut:
	fftFreeComplexArray(fftBuf);
	if(mfftPlan) {
		mfftFreePlan(mfftPlan);
	}
	return ourRtn;
}

#pragma mark --- main() ---

int main(int argc, char **argv)
{
	TestParams tp;
	size_t minRowSize = MIN_ROW_SIZE_DEF;
	size_t maxRowSize = MAX_ROW_SIZE_DEF;
	bool printBanner = true;
	size_t colSizeSpec = 0;
	
	memset(&tp, 0, sizeof(tp));
	tp.dumpData = false;
	
	int arg;
	while ((arg = getopt(argc, argv, "s:S:c:T:dpbhv")) != -1) {
		switch (arg) {
			case 's':
				minRowSize = fftParseStringRep(optarg);
				break;
			case 'S':
				maxRowSize = fftParseStringRep(optarg);
				break;
			case 'c':
				colSizeSpec = fftParseStringRep(optarg);
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
        fftPrintTestBanner("Two-dimension complex", "MatrixFFT", 
            FFT_DOUBLE_PREC ? true : false, 
            "2-D Complex Chirp Signal", NULL, 0, actThreads);
        printf("\n");
        printf("  rows | cols | samples | RMSE round | MXE round |   CHE   \n");
        printf(" ------+------+---------+------------+-----------+----------\n");
    }
    
	if(!fftIsPowerOfTwo(minRowSize, &tp.log2NumCols)) {
		printf("***size must be power of two\n");
		exit(1);
	}
	if(colSizeSpec != 0) {
        if(!fftIsPowerOfTwo(colSizeSpec, &tp.log2NumRows)) {
            printf("***Size must be a power of 2\n");
            exit(1);
        }
    }
    else {
        tp.log2NumRows = tp.log2NumCols;
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
		if(minRowSize == maxRowSize) {
			/* infer: user just wants this one run */
			break;
		}
		if(colSizeSpec != 0) {
			/* just this one */
			break;
		}
		if(tp.log2NumRows == tp.log2NumCols) {
			tp.log2NumRows++;
		}
		else {
			tp.log2NumCols++;
		}
	} while(((size_t)1 << tp.log2NumCols) <= maxRowSize);

	if(tp.doPause) {
		printf("Pausing after test, before frees for malloc debug; CR to continue: ");
		fflush(stdout);
		getchar();
	}
	
	/* clean up here */
	return ourRtn;
}
