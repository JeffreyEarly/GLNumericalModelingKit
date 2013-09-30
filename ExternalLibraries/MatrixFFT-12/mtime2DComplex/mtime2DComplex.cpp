/*	File: mtime2DComplex.cpp 
	
	Description:
		Measure timing of two-dimensional complex FFT 
	
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

#define LOOPS_DEF				10
#define MIN_ROW_SIZE_DEF		32
#define MAX_ROW_SIZE_DEF		(4 * 1024)

static void usage(char **argv)
{
	printf("usage: %s [options]\n", argv[0]);
	printf("Options:\n");
	printf("  -s minRowSize      -- minimum row size; default is %u\n", MIN_ROW_SIZE_DEF);
	printf("  -S maxRowSize      -- maximum row size; default is %u\n", MAX_ROW_SIZE_DEF);
	printf("  -c colSize         -- one run, minRowSize x colSize\n");
	printf("  -V                 -- vDSP; default is MatrixFFT\n");
	printf("  -o                 -- out-of-place; default is in-place\n");
	printf("  -f                 -- forward only\n");
	printf("  -Q                 -- square only\n");
	printf("  -u                 -- user time (default is wall time)\n");
	printf("  -l loops           -- default = %u\n", LOOPS_DEF);
	printf("  -k                 -- skip Matrix output transpositions\n");
	printf("  -1 (one)           -- don't time first loop\n");
	printf("  -T maxThreads      -- default is # of host cores\n");
    printf("  -P stripeSize      -- column stripe size (default is platform-specific)\n");
	printf("  -N                 -- manual normalization (default is lib)\n");
	printf("  -r                 -- override vDSP limits\n");
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
	size_t				numRows;
	size_t				numCols;
	bool				doublePrec;
	bool				verbose;
	bool				wallTime;
	unsigned			loops;
	FFT_Which			whichFFT;
	bool				skipTranspose;
	bool				dumpBuffers;
	bool				outOfPlace;
	bool				dontTimeFirstLoop;
	bool				forwardOnly;
	unsigned			numThreads;		 	/* 0 here means default as in API */
	bool				overrideLimits;
	bool				manualNorm;
    
    /*
     * bestTime is cumulative fastest (shortest) elapsed time of all 
     * runs at current size. Test updates bestTime if the current run 
     * is faster than the previous festest. 
     * Results are only displayed if lastRun is true (fastest time 
     * displayed) or displayAll is true (display current time regardless).
     */
    double              bestTime;
    bool                lastRun;
    bool                displayAll;
} TestParams;

static int doTest(
	TestParams *tp)
{
	int ourRtn = 0;
	unsigned actLoops;
	size_t total2DSamples = tp->numRows * tp->numCols;
	unsigned dims[2];
	
	unsigned log2TotalSamples;
	unsigned log2NumRows;
	unsigned log2NumCols;
	if(!fftIsPowerOfTwo(tp->numRows, &log2NumRows) ||
	   !fftIsPowerOfTwo(tp->numCols, &log2NumCols) ) {
		printf("***We only operate on powers of 2 only\n");
		return -1;
	}
	log2TotalSamples = log2NumRows + log2NumCols;
	unsigned maxLog = max(log2NumCols, log2NumRows);    /* for vDSP setup */
	
	if((tp->whichFFT == FW_vDSP) && 
	   (total2DSamples > VDSP_MAX_2D_COMPLEX) &&
	   !tp->overrideLimits) {
	   
		/*
		 * This facilitates large signal testing in the scripts without
		 * having to be limited by the vDSP limits.
		 */
		printf("   2^%-2u |  2^%-2u  |   2^%-2u  |     0.000     |  0.000\n",
			log2NumCols, log2NumRows, log2TotalSamples);
		return 0;
	}

	FFT_Setup fftSetup = NULL;
	MatrixFFTPlan mfftPlan = NULL;
	MFFTReturn mrtn;
	
	/* for vDSP */
	vDSPComplex zBufIn = {NULL, NULL};				
	vDSPComplex zBufOut = {NULL, NULL};				/* out of place */
	
	/* MatrixFFT */
	FFTComplex *fftSrc = NULL;						/* aligned MatrixFFT input */
	FFTComplex *fftDst = NULL;						/* aligned MatrixFFT output if OOP */
	FFTComplex *fftSrcFree = NULL;					/* to-be-freed MatrixFFT input */
	FFTComplex *fftDstFree = NULL;					/* to-be-freed MatrixFFT output if OOP */
	FFTComplex *actFftDst = NULL;
	
	double setupStart = 0.0;
	double setupEnd = 0.0;
	
	double startTime = 0.0;
	double endTime = 0.0;
	FFTFloat normFactor = 1.0 / total2DSamples;
	double elapsedTime;
	double ops = 5.0 * (double)total2DSamples * log2TotalSamples;
	if(tp->forwardOnly) { 
		ops *= tp->loops;
	}
	else {
		ops *= (2.0 * tp->loops);
		/* else we only execute one loop, forward only */
	}
	double CTGs;
	uint32_t flagsForward = 0;
	uint32_t flagsReverse = 0;
	uint32_t flagsCreate = 0;
	
	if(!tp->manualNorm) {
		flagsReverse = MEF_NormOutput;
	}

	/* 
	 * Create plans.
	 */
	if(tp->verbose) {
		printf("...setting up plans\n");
	}
	if(tp->verbose) {
		setupStart = fftGetTime(tp->wallTime);
	}
	
	switch(tp->whichFFT) {
		case FW_vDSP:
			/* Straight vDSP */
			fftSetup = FFTCreateSetup(maxLog);
			if (fftSetup == NULL) {
				printf("***Error: unable to create FFT setup.\n");
				return -1;
			}
			break;
		case FW_Matrix:
			dims[0] = log2NumRows;
			dims[1] = log2NumCols;
			mrtn = mfftCreatePlan(2, dims, false, flagsCreate, tp->numThreads, &mfftPlan);
			if(mrtn) {
				mfftPrintErrInfo("mfftCreatePlan", mrtn);
				return -1;
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
			if(fftAllocDSPComplex(&zBufIn, total2DSamples)) {
				printf("***Malloc failure for total2DSamples = %llu\n", 
					(unsigned long long)total2DSamples);
				ourRtn = -1;
				goto errOut;
			}
			if(tp->outOfPlace) {
				if(fftAllocDSPComplex(&zBufOut, total2DSamples)) {
					printf("***Malloc failure for total2DSamples = %llu", 
						(unsigned long long)total2DSamples);
					ourRtn = -1;
					goto errOut;
				}
			}
			genRandComplexDSP(&zBufIn, total2DSamples);
            fftFlushComplexDSP(&zBufIn, total2DSamples);
			break;

		case FW_Matrix:
		{
			fftSrc = fftAllocComplexArrayAlign(total2DSamples, FFT_MEM_ALIGNMENT, &fftSrcFree);
			if(fftSrc == NULL) {
				printf("***Malloc failure for total2DSamples = %llu\n",
					(unsigned long long)total2DSamples);
				ourRtn = -1;
				goto errOut;
			} 
			if(tp->outOfPlace) {
				fftDst = fftAllocComplexArrayAlign(total2DSamples, FFT_MEM_ALIGNMENT, &fftDstFree);
				if(fftDst == NULL) {
					printf("***Malloc failure for total2DSamples = %llu", 
						(unsigned long long)total2DSamples);
					ourRtn = -1;
					goto errOut;
				} 
			}			
			genRandComplex(fftSrc, total2DSamples);
            fftFlushComplex(fftSrc, total2DSamples);
		}
	}
	
	if(tp->dumpBuffers) {
		if(tp->whichFFT == FW_vDSP) {
			fftDump2DDSPComplex("Starting time domain", &zBufIn, tp->numRows, tp->numCols);
		}
		else {
			fftDump2DComplex("Starting time domain", fftSrc, tp->numRows, tp->numCols);
		}
	}

	/* MatrixFFT only */
	if(!tp->skipTranspose) {
		flagsForward |= MEF_TransposeOutput;
		flagsReverse |= MEF_TransposeInput;
	}

	if(tp->verbose) {
		printf("...performing %llu element FFT\n", 
			(unsigned long long)total2DSamples);
	}
	
	startTime = fftGetTime(tp->wallTime);
	actLoops = tp->loops;
	if(tp->dontTimeFirstLoop) {
		actLoops++;
	}
	switch(tp->whichFFT) {
		case FW_vDSP:
			for(unsigned loop=0; loop<actLoops; loop++) {
				if(tp->outOfPlace) {
					/* complex -- 2-d -- out of place */
					FFTComplex2dOP(fftSetup, &zBufIn, &zBufOut, log2NumCols, log2NumRows, FFT_FORWARD);
					if(!tp->forwardOnly) {
						FFTComplex2dOP(fftSetup, &zBufOut, &zBufIn, log2NumCols, log2NumRows, FFT_INVERSE);
					}
				}
				else {
					/* complex -- 2-d -- in place */
					FFTComplex2d(fftSetup, &zBufIn, log2NumCols, log2NumRows, FFT_FORWARD);
					if(!tp->forwardOnly) {
						FFTComplex2d(fftSetup, &zBufIn, log2NumCols, log2NumRows, FFT_INVERSE);
					}
				}
				if(!tp->forwardOnly) {
					FFTvScale(zBufIn.realp, normFactor, total2DSamples);	
					FFTvScale(zBufIn.imagp, normFactor, total2DSamples);	
				}
				if(tp->dontTimeFirstLoop && (loop == 0)) {
					/* restart timer */
					startTime = fftGetTime(tp->wallTime);
				}
			}
			break;
			
		case FW_Matrix:
			actFftDst = tp->outOfPlace ? fftDst : fftSrc;
			for(unsigned loop=0; loop<actLoops; loop++) {
				mrtn = mfftExecute(mfftPlan, flagsForward, true, fftSrc, actFftDst);
				if(mrtn) {
					mfftPrintErrInfo("mfftExecute", mrtn);
					ourRtn = -1;
					goto errOut;
				}
				if(!tp->forwardOnly) {
					mrtn = mfftExecute(mfftPlan, flagsReverse, false, actFftDst, fftSrc);
					if(mrtn) {
						mfftPrintErrInfo("mfftExecute", mrtn);
						ourRtn = -1;
						goto errOut;
					}
					
					if(tp->manualNorm) {
						/* 
						 * scale by 1/N
						 */
						fftScaleComplex(fftSrc, normFactor, total2DSamples);
					}
				}
				if(tp->dontTimeFirstLoop && (loop == 0)) {
					/* restart timer */
					startTime = fftGetTime(tp->wallTime);
				}
			}
			break;
	}
	
	endTime = fftGetTime(tp->wallTime);
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
	
	if(tp->dumpBuffers) {
		if(tp->whichFFT == FW_vDSP) {
			fftDump2DDSPComplex("Ending time domain", &zBufIn, tp->numRows, tp->numCols);
		}
		else {
			fftDump2DComplex("Ending time domain", fftSrc, tp->numRows, tp->numCols);
		}
	}
	

	if(tp->verbose) {
		printf("complexity = 5 * totalSamples * lg(n) * 2 * loops\n");
		printf("           = 5 * %llu * %u * 2 * %u\n",
			(unsigned long long)total2DSamples, log2TotalSamples, tp->loops);
		printf("           = %.1f\n", ops);
		printf("Setup time = %.4f s\n", setupEnd - setupStart);
	}
	printf("   2^%-2u |  2^%-2u  |   2^%-2u  | %9.3f     | %6.3f\n",
		log2NumCols, log2NumRows, 
		log2TotalSamples,
		elapsedTime,
		CTGs);
	
errOut:
	
	fftFreeComplexArrayAlign(fftSrc, fftSrcFree);
	if(tp->outOfPlace) {
		fftFreeComplexArrayAlign(fftDst, fftDstFree);
	}
	
	COND_FREE(zBufIn.realp);
	COND_FREE(zBufIn.imagp);
	if(tp->outOfPlace) {
		COND_FREE(zBufOut.realp);
		COND_FREE(zBufOut.imagp);
	}
	
	if(fftSetup) {
		FFTFreeSetup(fftSetup);
	}
	if(mfftPlan)  {
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
	size_t minRowSize = MIN_ROW_SIZE_DEF;
	size_t maxRowSize = MAX_ROW_SIZE_DEF;
	int ourRtn = 0;
	bool allowError = false;
	size_t colSizeSpec = 0;
	bool squareOnly = false;
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
	while ((arg = getopt(argc, argv, "s:S:c:Vkoul:1pDAfT:NQrP:Bad:vh")) != -1) {
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
			case 'V':
				tp.whichFFT = FW_vDSP;
				implStr = "Accelerate";
				break;
			case 'k':
				tp.skipTranspose = true;
				break;
			case 'o':
				tp.outOfPlace = true;
				break;
			case 'u':
				tp.wallTime = false;
				break;
			case 'l':
				tp.loops = atoi(optarg);
				break;
			case '1':
				tp.dontTimeFirstLoop = true;
				appendOptStr(optStr, "Don't time first loop");
				break;
			case 'p':
				doPause = true;
				break;
			case 'f':
				tp.forwardOnly = true;
				break;
			case 'D':
				tp.dumpBuffers = true;
				break;
			case 'v':
				tp.verbose = true;
				break;
			case 'A':
				allowError = true;
				break;
			case 'T':
				tp.numThreads = atoi(optarg);
				break;
			case 'Q':
				squareOnly = true;
				break;
			case 'r':
				tp.overrideLimits = true;
				break;
			case 'N':
				tp.manualNorm = true;
				appendOptStr(optStr, "Manual normalize");
				break;
            case 'P':
                columnStripeSize = strtol(optarg, NULL, 0);
                break;
			case 'B':
				printBanner = false;
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
	if(columnStripeSize) {
        /* Note this is always complex stripe size, real not used here */
		mfftSetColumnStripeSize(false, columnStripeSize);
		char str[100];
		sprintf(str, "columnStripeSize=%llu", (unsigned long long)columnStripeSize);
		appendOptStr(optStr, str);
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
	if(forceVdsp != MF_Default) {
        mfftSetForceVdsp(forceVdsp);
		appendOptStr(optStr, (forceVdsp == MF_ForceVdsp) ? "Force vDSP" : "Force !vDSP");
    }
		
	if(printBanner) {
        fftPrintTestBanner("Two-dimension complex", implStr, tp.doublePrec, "Random", 
            optStr, tp.loops, actThreads);
        printf("\n");
        
        printf("  Width | Height | Samples | %s time (s) |  CTGs \n",
            tp.wallTime ? "Wall" : "User");
        printf(" -------+--------+---------+---------------+---------\n");
    }
    
	if(doPause) {
		fpurge(stdin);
		printf("Pausing at top of loop; CR to proceed: ");
		getchar();
	}
	
	tp.numCols = minRowSize;
	if(colSizeSpec != 0) {
		tp.numRows = colSizeSpec;
	}
	else {
		tp.numRows = tp.numCols;
	}
	
	do {
        unsigned numIter = fftIterationsForSize(tp.numRows * tp.numCols);
        
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
                printf("Pausing at end of loop; CR to continue: ");
                fflush(stdout);
                getchar();
            }
        }
		if(minRowSize == maxRowSize) {
			/* infer: user just wants this one run */
			break;
		}
		if(colSizeSpec != 0) {
			/* just this one */
			break;
		}
		if(tp.numRows == tp.numCols) {
			tp.numRows <<= 1;
			if(squareOnly) {
				tp.numCols <<= 1;
			}
		}
		else {
			tp.numCols <<= 1;
		}
	} while(tp.numCols <= maxRowSize);
	
	return allowError ? 0 : ourRtn;
}
