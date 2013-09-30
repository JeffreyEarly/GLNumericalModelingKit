/*	File: mtime1DComplex.cpp 
	
	Description:
		Measure timing of one-dimensional complex FFT 
	
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

#define FFT_MIN_SIZE_DEF		32
#define FFT_MAX_SIZE_DEF		(1024 * 1024)
#define LOOPS_DEF				10

static void usage(char **argv)
{
	printf("usage: %s [options]\n", argv[0]);
	printf("Options:\n");
	printf("  -s minFftSize      -- minimum size in samples; default is %u\n", FFT_MIN_SIZE_DEF);
	printf("  -S maxFftSize      -- maximum size in samples; default is %u\n", FFT_MAX_SIZE_DEF);
	printf("  -V                 -- vDSP; default is MatrixFFT\n");
	printf("  -o                 -- out-of-place; default is in-place\n");
	printf("  -u                 -- user time (default is wall time)\n");
	printf("  -l loops           -- default = %u\n", LOOPS_DEF);
	printf("  -k                 -- skip Matrix output transpositions\n");
	printf("  -1 (one)           -- don't time first loop\n");
	printf("  -T maxThreads      -- default is # of host cores\n");
	printf("  -O offset          -- signed rectangle offset (default is platform-specific)\n");
    printf("  -P stripeSize      -- column stripe size (default is platform-specific)\n");
	printf("  -r                 -- override vDSP limits\n");
	printf("  -N                 -- manual normalization (default is lib)\n");
    printf("  -H                 -- set MCF_HintTranspose\n");
	printf("  -p                 -- pause for Shark/MallocDebug\n");
	printf("  -D                 -- dump buffers\n");
    printf("  -d on|off          -- force raw vDSP on/off (default is per config)\n");
	printf("  -A                 -- allow error (i.e. abort but return OK status)\n");
    printf("  -a                 -- show all results, not just best\n");
	printf("  -B                 -- no banner\n");
	printf("  -v                 -- verbose\n");
	exit(1);
}

/* Which implementation */
typedef enum {
	FW_vDSP,
	FW_Matrix
} FFT_Which;

typedef struct {
	size_t				numSamples;
	bool				doublePrec;
	bool				verbose;
	bool				wallTime;
	unsigned			loops;
	FFT_Which			whichFFT;
	bool				skipTranspose;
    bool                hintTranspose;
	bool				dumpBuffers;
	bool				outOfPlace;
	bool				dontTimeFirstLoop;
	unsigned			numThreads;			/* 0 here means default as in API */
	bool				display2D;			/* show size as rows x cols */
	bool				overrideLimits;
	bool				manualNorm;

    /*
     * bestTime is cumulative fastest (shortest) elapsed time of all 
     * runs at current size. Test updates bestTime if the current run 
     * is faster than the previous fastest. 
     * Results are only displayed if lastRun is true (fastest time 
     * displayed) or displayAll is true (display current time regardless).
     */
    double              bestTime;
    bool                lastRun;
    bool                displayAll;
} TestParams;

static void makeStringRep(
	size_t n,
	char *str)
{
	unsigned log2n;
	fftIsPowerOfTwo(n, &log2n);
	sprintf(str, "2^%u", log2n);
}

static int doTest(
	TestParams *tp)
{
	int ourRtn = 0;
	unsigned actLoops;
	MFFTReturn mrtn;
	char *numSampStr = NULL;
	
	if((tp->whichFFT == FW_vDSP) && 
	   (tp->numSamples > VDSP_MAX_1D_COMPLEX) &&
	   !tp->overrideLimits) {
	   
		/*
		 * This facilitates large signal testing in the scripts without
		 * having to be limited by the vDSP limits.
		 */
		numSampStr = fftStringRepPow2(tp->numSamples);
		printf(" %8s |     0.000     | 0.000\n", numSampStr);
		free(numSampStr);
		return 0;
	}
	
	/* 
	 * Create plans.
	 */
	if(tp->verbose) {
		printf("...setting up plan\n");
	}
	FFT_Setup fftSetup = NULL;
	MatrixFFTPlan mfftPlan = NULL;
	
	unsigned log2NumSamples;
	if(!fftIsPowerOfTwo(tp->numSamples, &log2NumSamples)) {
		printf("***vDSP and MatrixFFT operate on powers of 2 only\n");
		return -1;
	}
	
    /* temp test */
    unsigned loops = tp->loops;
    if(log2NumSamples <= 22) {
        //loops <<= 1;
    } 
    /* end temp test */
    
	/* for vDSP */
	vDSPComplex zBufIn = {NULL, NULL};				
	vDSPComplex zBufOut = {NULL, NULL};				/* out of place */
	
	/* MatrixFFT */
	FFTComplex *fftSrc = NULL;						/* aligned MatrixFFT input */
	FFTComplex *fftDst = NULL;						/* aligned MatrixFFT output if OOP */
	FFTComplex *fftSrcFree = NULL;					/* to-be-freed MatrixFFT input */
	FFTComplex *fftDstFree = NULL;					/* to-be-freed MatrixFFT output if OOP */
	size_t totalSamples = tp->numSamples;
	
	double setupStart = 0.0;
	double setupEnd = 0.0;
	
	double startTime = 0.0;
	double endTime = 0.0;
	FFTFloat normFactor = 1.0 / tp->numSamples;
	double elapsedTime;
	double CTGs;
	uint32_t flagsForward = 0;
	uint32_t flagsReverse = 0;
	uint32_t createFlags = 0;
	double ops = 5.0 * (double)totalSamples * (double)log2NumSamples * 2.0 * (double)loops;
	
	if(!tp->manualNorm) {
		flagsReverse = MEF_NormOutput;
	}
    if(tp->hintTranspose) {
        createFlags |= MCF_HintTranspose;
    }

	/* Create plans */
	if(tp->verbose) {
		setupStart = fftGetTime(tp->wallTime);
	}	
	switch(tp->whichFFT) {
		case FW_vDSP:
			/* Straight vDSP */
			fftSetup = FFTCreateSetup(log2NumSamples);
			if (fftSetup == NULL) {
				printf("***Error: unable to create FFT setup.\n");
				return -1;
			}
			break;
			
		case FW_Matrix:
			mrtn = mfftCreatePlan(1, &log2NumSamples, false, createFlags, tp->numThreads, &mfftPlan);
			if(mrtn) {
				mfftPrintErrInfo("mfftCreatePlan", mrtn);
				ourRtn = -1;
				goto errOut;
			}
			break;
	}

	if(tp->verbose) {
		setupEnd = fftGetTime(tp->wallTime);
	}
	
	/*
	 * Alloc and init SplitBuffers. 
	 */
	if(tp->verbose) {
		printf("...setting up buffers\n");
	}
	switch(tp->whichFFT) {
		case FW_vDSP:
			if(fftAllocDSPComplex(&zBufIn, tp->numSamples)) {
				printf("***Malloc failure for numSamples = %llu\n", 
					(unsigned long long)tp->numSamples);
				ourRtn = -1;
				goto errOut;
			}
			if(tp->outOfPlace) {
				if(fftAllocDSPComplex(&zBufOut, tp->numSamples)) {
					printf("***Malloc failure for numSamples = %llu", 
						(unsigned long long)tp->numSamples);
					ourRtn = -1;
					goto errOut;
				}
			}
			genRandComplexDSP(&zBufIn, tp->numSamples);
            fftFlushComplexDSP(&zBufIn, tp->numSamples);
			break;

		case FW_Matrix:
		{
			fftSrc = fftAllocComplexArrayAlign(totalSamples, FFT_MEM_ALIGNMENT, &fftSrcFree);
			if(fftSrc == NULL) {
				printf("***Malloc failure for totalSamples = %llu\n", 
					(unsigned long long)totalSamples);
				ourRtn = -1;
				goto errOut;
			} 
			if(tp->outOfPlace) {
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
            fftFlushComplex(fftSrc, totalSamples);
			break;
		}
	}
	
	if(tp->dumpBuffers) {
		if(tp->whichFFT == FW_vDSP) {
			fftDump1DDSPComplex("Starting time domain", &zBufIn, tp->numSamples);
		}
		else {
			fftDump1DComplex("Starting time domain", fftSrc, totalSamples);
		}
	}

	/* MatrixFFT only */
	if(!tp->skipTranspose) {
		flagsForward |= MEF_TransposeOutput;
		flagsReverse |= MEF_TransposeInput;
	}
	
	if(tp->verbose) {
		printf("...performing %llu element FFT\n", (unsigned long long)totalSamples);
	}
	
	startTime = fftGetTime(tp->wallTime);

	actLoops = loops;
	if(tp->dontTimeFirstLoop) {
		actLoops++;
	}
	switch(tp->whichFFT) {
		case FW_vDSP:
			for(unsigned loop=0; loop<actLoops; loop++) {
				if(tp->outOfPlace) {
					/* complex -- 1-d -- out of place */
					FFTComplex1dOP(fftSetup, &zBufIn, &zBufOut, log2NumSamples, FFT_FORWARD);
					FFTComplex1dOP(fftSetup, &zBufOut, &zBufIn, log2NumSamples, FFT_INVERSE);
				}
				else {
					/* complex -- 1-d -- double -- in place */
					FFTComplex1d(fftSetup, &zBufIn, log2NumSamples, FFT_FORWARD);
					FFTComplex1d(fftSetup, &zBufIn, log2NumSamples, FFT_INVERSE);
				}
				FFTvScale(zBufIn.realp, normFactor, tp->numSamples);	
				FFTvScale(zBufIn.imagp, normFactor, tp->numSamples);	
				if(tp->dontTimeFirstLoop && (loop == 0)) {
					/* restart timer */
					startTime = fftGetTime(tp->wallTime);
				}
			}
			break;
		
		case FW_Matrix:
			for(unsigned loop=0; loop<actLoops; loop++) {
				mrtn = mfftExecute(mfftPlan, flagsForward, true, fftSrc, fftDst);
				if(mrtn) {
					mfftPrintErrInfo("mfftExecute", mrtn);
					ourRtn = -1;
					goto errOut;
				}
				mrtn = mfftExecute(mfftPlan, flagsReverse, false, fftDst, fftSrc);
				if(mrtn) {
					mfftPrintErrInfo("mfftExecute", mrtn);
					ourRtn = -1;
					goto errOut;
				}
				
				if(tp->manualNorm) {
					/* 
					 * scale by 1/N
					 */
					fftScaleComplex(fftSrc, normFactor, totalSamples);
				}
				
				if(tp->dontTimeFirstLoop && (loop == 0)) {
					/* restart timer */
					startTime = fftGetTime(tp->wallTime);
				}
			}
			break;
	}

	endTime = fftGetTime(tp->wallTime);

	if(tp->dumpBuffers) {
		if(tp->whichFFT == FW_vDSP) {
			fftDump1DDSPComplex("Ending time domain", &zBufIn, tp->numSamples);
		}
		else {
			fftDump1DComplex("Ending time domain", fftSrc, totalSamples);
		}
	}
		
	elapsedTime = endTime - startTime;
    
    if(!tp->displayAll) {
        /*
         * Accumulate best time
         */
        if((tp->bestTime == 0.0) ||             // first time thru
           (tp->bestTime > elapsedTime)) {      // new best
            tp->bestTime = elapsedTime;
        }
        if(!tp->lastRun) {
            /* We're done, no display this time thru */
            goto errOut;
        }
        
        /* Last run: display cumulative best */
        elapsedTime = tp->bestTime;
    }
    
	CTGs = (ops / elapsedTime) / 1.0e+9;
	
	if(tp->verbose) {
		printf("complexity = 5 * totalSamples * lg(n) * 2 * loops\n");
		printf("           = 5 * %llu * %u * 2 * %u\n",
				(unsigned long long)totalSamples, log2NumSamples, tp->loops);
		printf("           = %.1f\n", ops);
		printf("Setup time = %.4f s\n", setupEnd - setupStart);
	}
	
	if(tp->display2D) {
		size_t numRows;
		size_t numCols;
		mfftRectangle(mfftPlan, &numRows, &numCols);
		char rowStr[20];
		char colStr[20];
		makeStringRep(numRows, rowStr);
		makeStringRep(numCols, colStr);
		printf(" %4s x %4s | %9.3f     | %5.3f\n", rowStr, colStr, elapsedTime, CTGs);
	}
	else {
		numSampStr = fftStringRepPow2(tp->numSamples);
		printf(" %8s | %9.3f     | %5.3f\n", numSampStr, elapsedTime, CTGs);
	}
errOut:
	COND_FREE(numSampStr);
	fftFreeComplexArrayAlign(fftSrc, fftSrcFree);
	if(tp->outOfPlace && (fftDst != NULL)) {
		fftFreeComplexArrayAlign(fftDst, fftDstFree);
	}

	COND_FREE(zBufIn.realp);
	COND_FREE(zBufIn.imagp);
	COND_FREE(zBufOut.realp);
	COND_FREE(zBufOut.imagp);
	
	if(fftSetup) {
		FFTFreeSetup(fftSetup);
	}

	if(mfftPlan) {
		mfftFreePlan(mfftPlan);
	}
	return ourRtn;	
}
	
int main(int argc, char **argv)
{
	const char *implStr = "MatrixFFT";
	char optStr[200];
	TestParams tp;
	bool doPause = false;
	size_t minSize = FFT_MIN_SIZE_DEF;
	size_t maxSize = FFT_MAX_SIZE_DEF;
	int ourRtn = 0;
	bool allowError = false;
	int rectOffset = 0;
    bool rectOffsetSpecd = false;
	bool printBanner = true;
	size_t columnStripeSize = 0;
	MFFT_ForceVdsp forceVdsp = MF_Default;
    
	memset(&tp, 0, sizeof(tp));
	optStr[0] = '\0';
	
	tp.wallTime = true;
	tp.doublePrec = FFT_DOUBLE_PREC ? true : false;
	tp.loops = LOOPS_DEF;
	tp.whichFFT = FW_Matrix;
	
	int arg;
	while ((arg = getopt(argc, argv, "s:S:vul:VkpDo1T:O:AaBrNHP:d:h")) != -1) {
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
			case 'V':
				tp.whichFFT = FW_vDSP;
				implStr = "Accelerate";
				break;
			case 'k':
				tp.skipTranspose = true;
				break;
			case 'p':
				doPause = true;
				break;
			case 'D':
				tp.dumpBuffers = true;
				break;
			case 'o':
				tp.outOfPlace = true;
				break;
			case '1':
				tp.dontTimeFirstLoop = true;
				appendOptStr(optStr, "Don't time first loop");
				break;
			case 'T':
				tp.numThreads = atoi(optarg);
				break;
			case 'O':
				rectOffset = atoi(optarg);
                rectOffsetSpecd = true;
				tp.display2D = true;
				break;
			case 'A':
				allowError = true;
				break;
			case 'B':
				printBanner = false;
				break;
			case 'r':
				tp.overrideLimits = true;
				break;
			case 'N':
				tp.manualNorm = true;
				appendOptStr(optStr, "Manual normalize");
				break;
            case 'H':
                tp.hintTranspose = true;
				appendOptStr(optStr, "HintTranspose");
                break;
            case 'P':
                columnStripeSize = strtol(optarg, NULL, 0);
                break;
            case 'a':
                tp.displayAll = true;
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
	if((maxSize == FFT_MAX_SIZE_DEF) && (minSize > maxSize)) {
		/* allow specification of single large size easily */
		maxSize = minSize;
	}
	if(minSize > maxSize) {
		printf("***maxSize must be greater than or equal to minSize\n");
		exit(1);
	}
	if(tp.outOfPlace) {
		appendOptStr(optStr, "Out-of-place");
	}
	else {
		appendOptStr(optStr, "In-place");
	}
	if(tp.skipTranspose) {
		if(tp.whichFFT != FW_Matrix) {
			printf("***transpose options only valid for MatrixFFT\n");
			exit(1);
		}
		appendOptStr(optStr, "Skip transposition");
	}
	
	unsigned actThreads = tp.numThreads;
	if((actThreads == 0) && (tp.whichFFT != FW_vDSP)) {
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
	if(columnStripeSize) {
        /* Note this is always complex stripe size, real not used here */
		mfftSetColumnStripeSize(false, columnStripeSize);
		char str[100];
		sprintf(str, "columnStripeSize=%llu", (unsigned long long)columnStripeSize);
		appendOptStr(optStr, str);
	}
	if(forceVdsp != MF_Default) {
        mfftSetForceVdsp(forceVdsp);
		appendOptStr(optStr, (forceVdsp == MF_ForceVdsp) ? "Force vDSP" : "Force !vDSP");
    }

	if(printBanner) {
		fftPrintTestBanner("One-dimension complex", implStr, tp.doublePrec, "Random", 
			optStr, tp.loops, actThreads);
		printf("\n");
		
		if(tp.display2D) {
			printf("    Samples  | %s time (s) |  CTGs \n", tp.wallTime ? "Wall" : "User");
			printf(" ------------+---------------+--------\n");
		}
		else {
			printf("  Samples | %s time (s) |  CTGs \n", tp.wallTime ? "Wall" : "User");
			printf(" ---------+---------------+--------\n");
		}
	}
	
	if(doPause) {
		fpurge(stdin);
		printf("Pausing at top of loop; CR to proceed: ");
		getchar();
	}
		
	for(tp.numSamples=minSize; tp.numSamples<=maxSize; tp.numSamples <<= 1) {
        unsigned numIter = fftIterationsForSize(tp.numSamples);
        
        tp.bestTime = 0.0;
        tp.lastRun = false;
        for(unsigned iter=1; iter<=numIter; iter++) {
            if(iter == numIter) {
                tp.lastRun = true;
            }
            ourRtn = doTest(&tp);
            if(ourRtn) {
                break;
            }
            if(doPause) {
                fpurge(stdin);
                printf("Pausing at end of loop; CR to proceed: ");
                getchar();
            }
        }
	}
	
	return allowError ? 0 : ourRtn;
}
