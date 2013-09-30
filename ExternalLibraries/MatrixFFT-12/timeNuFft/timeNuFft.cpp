/*	File: timeNuFft.cpp 
	
	Description:
		Measure timing of one-dimensional Nonuniform FFT
	
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
 * Created 4/17/2009. 
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
#include <libMatrixFFT/NUFFT.h>
#include <libMatrixFFT/vdspUtils.h>

#define FFT_MIN_SIZE_DEF		32
#define FFT_MAX_SIZE_DEF		(1024 * 1024)
#define LOOPS_DEF				10
#define BITS_DEF                16

static void usage(char **argv)
{
	printf("usage: %s [options]\n", argv[0]);
	printf("Options:\n");
	printf("  -s minFftSize      -- minimum size in samples; default is %u\n", FFT_MIN_SIZE_DEF);
	printf("  -S maxFftSize      -- maximum size in samples; default is %u\n", FFT_MAX_SIZE_DEF);
    printf("  -b bits            -- bits of accuracy; default is %u\n", BITS_DEF);
    printf("  -r                 -- test reference NUFFT; default is optimized\n");
    printf("  -i                 -- in-place; default is OOP\n");
	printf("  -u                 -- user time (default is wall time)\n");
	printf("  -l loops           -- default = %u\n", LOOPS_DEF);
	printf("  -T maxThreads      -- default is # of host cores\n");
	printf("  -p                 -- pause for Shark/MallocDebug\n");
	printf("  -P                 -- pause for Shark/MallocDebug INSIDE the test loop\n");
	printf("  -D                 -- dump buffers\n");
	printf("  -A                 -- allow error (i.e. abort but return OK status)\n");
	printf("  -B                 -- no banner\n");
	printf("  -v                 -- verbose\n");
	exit(1);
}

typedef struct {
	size_t				numSamples;
	bool				doublePrec;
    bool                testReference;
    bool                inPlace;
    bool                testPause;
    unsigned            bits;
	bool				verbose;
	bool				wallTime;
	unsigned			loops;
	bool				dumpBuffers;
	unsigned			numThreads;			/* 0 here means default as in API */
} TestParams;

static int doTest(
	TestParams *tp)
{
	int ourRtn = 0;
	MFFTReturn mrtn;
	char *numSampStr = NULL;
	
	if(tp->testPause) {
		fpurge(stdin);
		printf("Pausing at top of test before any mallocs\n");
        printf("CR to proceed: ");
		getchar();
	}

    /*
     * Optimized NUFFT modifies the input buffer. To avoid input degenerating into
     * zeroes we have to essentially stop the clock, restore input, and restart the clock
     * after every loop. Also do this if running reference implementation in-place.
     */
    if(tp->inPlace && !tp->testReference) {
        printf("***Hey! Optimized NUFFT can't run in-place!\n");
        return -1;
    }
    
    bool restoreInput = !tp->testReference || tp->inPlace;
    
	/* 
	 * Create plan.
	 */
	if(tp->verbose) {
		printf("...setting up plan\n");
	}
	NUFFTPlan nuFftPlan = NULL;
	
	unsigned log2NumSamples;
	if(!fftIsPowerOfTwo(tp->numSamples, &log2NumSamples)) {
		printf("***NUFFT operate on powers of 2 only\n");
		return -1;
	}
	
    FFTComplex *rawSrc = NULL;                      /* if restoreInput */
	FFTComplex *fftSrc = NULL;						/* aligned input */
	FFTComplex *fftDst = NULL;						/* aligned output if OOP */
	FFTComplex *fftSrcFree = NULL;					/* to-be-freed input */
	FFTComplex *fftDstFree = NULL;					/* to-be-freed output if OOP */
    FFTFloat *tau = NULL;                           /* nuFFT tau */
	size_t totalSamples = tp->numSamples;
	
	double setupStart = 0.0;
	double setupEnd = 0.0;
	
	double startTime = 0.0;
	double endTime = 0.0;
	double elapsedTime;
    double msPerLoop;
    uint32_t implFlag = tp->testReference ? NF_Reference : NF_Optimized;
	    
    /*
     * If restoreInput is true, accumulate the time spent copying from 
     * rawSrc to fftSrc and subtract that from the total at the end of 
     * the loop.
     */
    double restoreTotal = 0.0;

	/* Create plan */
	if(tp->verbose) {
		setupStart = fftGetTime(tp->wallTime);
	}	
    mrtn = nufftCreatePlan(1, &log2NumSamples, tp->bits, implFlag, tp->numThreads, &nuFftPlan);
    if(mrtn) {
        mfftPrintErrInfo("nufftCreatePlan", mrtn);
        ourRtn = -1;
        goto errOut;
    }
	if(tp->verbose) {
		setupEnd = fftGetTime(tp->wallTime);
	}
	
	/*
	 * Alloc and init buffers. 
	 */
	if(tp->verbose) {
		printf("...setting up buffers\n");
	}
    fftSrc = fftAllocComplexArrayAlign(totalSamples, FFT_MEM_ALIGNMENT, &fftSrcFree);
    if(fftSrc == NULL) {
        printf("***Malloc failure for totalSamples = %llu\n", 
            (unsigned long long)totalSamples);
        ourRtn = -1;
        goto errOut;
    } 
    if(!tp->inPlace) {
        fftDst = fftAllocComplexArrayAlign(totalSamples, FFT_MEM_ALIGNMENT, &fftDstFree);
        if(fftDst == NULL) {
            printf("***Malloc failure for totalSamples = %llu", 
                (unsigned long long)totalSamples);
            ourRtn = -1;
            goto errOut;
        } 
    }
    else {
        fftDst = fftSrc;
    }
    genRandComplex(fftSrc, totalSamples);
    
    if(restoreInput) {
        rawSrc = fftAllocComplexArray(totalSamples);
        fftCopyComplexArray(fftSrc, rawSrc, totalSamples);
    }
    
    tau = (FFTFloat *)malloc(totalSamples * sizeof(FFTFloat));
    if(tau == NULL) {
		printf("***Malloc failure for totalSamples = %llu\n", 
			(unsigned long long)totalSamples);
		ourRtn = -1;
		goto errOut;
    }
    fftGenTau(tau, totalSamples);
	
	if(tp->dumpBuffers) {
		fftDump1DComplex("Starting time domain", fftSrc, totalSamples);
		fftDump1DFloat("tau", tau, totalSamples);
	}

	if(tp->verbose) {
		printf("...performing %llu element NUFFT\n", (unsigned long long)totalSamples);
	}
	
	startTime = fftGetTime(tp->wallTime);

    for(unsigned loop=0; loop<tp->loops; loop++) {
        mrtn = nuFftExecute(nuFftPlan, implFlag, fftSrc, tau, fftDst);
        if(mrtn) {
            mfftPrintErrInfo("nuFftExecute", mrtn);  
            ourRtn = -1;
            goto errOut;
        }
        if(restoreInput) {
            double restoreStart = fftGetTime(tp->wallTime);
            fftCopyComplexArray(rawSrc, fftSrc, totalSamples);
            double restoreEnd   = fftGetTime(tp->wallTime);
            restoreTotal += (restoreEnd - restoreStart);
        }

    }
	
	endTime = fftGetTime(tp->wallTime) - restoreTotal;

	if(tp->dumpBuffers) {
		fftDump1DComplex("Ending time domain", fftSrc, totalSamples);
	}
		
	elapsedTime = endTime - startTime;
	msPerLoop = 1000.0 * elapsedTime / (double)tp->loops;
    
	numSampStr = fftStringRepPow2(tp->numSamples);
	printf(" %8s | %9.3f     | %5.2f\n", numSampStr, elapsedTime, msPerLoop);

	if(tp->testPause) {
		fpurge(stdin);
		printf("Pausing at end of test before any frees\n");
        printf("CR to proceed: ");
		getchar();
	}

errOut:
	COND_FREE(numSampStr);
	fftFreeComplexArrayAlign(fftSrc, fftSrcFree);
	fftFreeComplexArrayAlign(fftDst, fftDstFree);
    fftFreeComplexArray(rawSrc);
    COND_FREE(tau);
    
	if(nuFftPlan) {
		nufftFreePlan(nuFftPlan);
	}

	return ourRtn;	
}
	
int main(int argc, char **argv)
{
	char optStr[200];
	TestParams tp;
	bool doPause = false;
	size_t minSize = FFT_MIN_SIZE_DEF;
	size_t maxSize = FFT_MAX_SIZE_DEF;
	int ourRtn = 0;
	bool allowError = false;
	bool printBanner = true;
	
	memset(&tp, 0, sizeof(tp));
	optStr[0] = '\0';
	
	tp.wallTime = true;
	tp.doublePrec = FFT_DOUBLE_PREC ? true : false;
	tp.loops = LOOPS_DEF;
	tp.bits = BITS_DEF;
    tp.testReference = false;
    
	int arg;
	while ((arg = getopt(argc, argv, "s:S:b:iul:T:pPDABvrh")) != -1) {
		switch (arg) {
			case 's':
				minSize = fftParseStringRep(optarg);
				break;
			case 'S':
				maxSize = fftParseStringRep(optarg);
				break;
			case 'v':
				tp.verbose = true;
				break;
			case 'u':
				tp.wallTime = false;
				break;
			case 'l':
				tp.loops = atoi(optarg);
				break;
			case 'p':
				doPause = true;
				break;
			case 'P':
				tp.testPause = true;
				break;
			case 'D':
				tp.dumpBuffers = true;
				break;
			case 'i':
				tp.inPlace = true;
				break;
			case 'T':
				tp.numThreads = atoi(optarg);
				break;
			case 'A':
				allowError = true;
				break;
			case 'B':
				printBanner = false;
				break;
			case 'b':
				tp.bits = atoi(optarg);
				break;
            case 'r':
                tp.testReference = true;
                break;
			case 'h':
			default:
				usage(argv);
		}
	}
	if(optind != argc) {
		usage(argv);
	}
	if((maxSize == FFT_MAX_SIZE_DEF) && (minSize > maxSize)) {
		/* allow specification of single large size easily */
		maxSize = minSize;
	}
	if(minSize > maxSize) {
		printf("***maxSize must be greater than or equal to minSize\n");
		exit(1);
	}
	if(tp.inPlace) {
		appendOptStr(optStr, "In-place");
	}
	else {
		appendOptStr(optStr, "Out-of-place");
	}
    
    char bitStr[100];
    sprintf(bitStr, "bits=%u", tp.bits);
	appendOptStr(optStr, bitStr);
	
	unsigned actThreads = tp.numThreads;
	if(actThreads == 0) {
		/* Find out how many are actually going to be used */
		MFFTReturn mrtn = mfftNumThreads(NULL, &actThreads);
		if(mrtn) {
			mfftPrintErrInfo("mfftNumThreads", mrtn);
			exit(1);
		}
	}

	if(printBanner) {
		fftPrintTestBanner("1-D NUFFT", 
            tp.testReference ? "Reference NUFFT" : "Optimized NUFFT", 
            tp.doublePrec, "Random", 
			optStr, tp.loops, actThreads);
		printf("\n");
		
		printf("  Samples | %s time (s) | ms/loop \n", tp.wallTime ? "Wall" : "User");
		printf(" ---------+---------------+----------\n");
	}
	
	if(doPause) {
		fpurge(stdin);
		printf("Pausing at top of loop\n");
        printf("CR to proceed: ");
		getchar();
	}
		
	for(tp.numSamples=minSize; tp.numSamples<=maxSize; tp.numSamples <<= 1) {
		ourRtn = doTest(&tp);
		if(ourRtn) {
			break;
		}
		if(doPause) {
			fpurge(stdin);
			printf("Pausing at end of loop\n");
            printf("CR to proceed: ");
			getchar();
		}
	}
	
	return allowError ? 0 : ourRtn;
}
