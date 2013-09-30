/*	File: vfyNuFft.cpp 
	
	Description:
		Verify proper operation of 1-D nonuniform FFT.
	
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
 * Created 04/10/2009. 
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
#include <libMatrixFFT/devRandom.h>
#include <libMatrixFFT/src/PolyComplex.h>       /* SPI */

#define MIN_SIZE_DEF		64
#define MAX_SIZE_DEF		(128 * 1024)
#define BITS_DEF            16

static void usage(char **argv)
{
	printf("usage: %s [options]\n", argv[0]);
	printf("Options:\n");
	printf("  -s minSize         -- minimum size; default is %u\n", MIN_SIZE_DEF);
	printf("  -S maxSize         -- maximum size; default is %u\n", MAX_SIZE_DEF);
    printf("  -b bits            -- bits of accuracy; default is %u\n", BITS_DEF);
    printf("  -u                 -- uniform tau; default is random\n");
	printf("  -i                 -- incrementing data; default is random\n"); 
    printf("  -m                 -- compare against FFT; default is NU DFT\n");
    printf("  -r                 -- test reference NUFFT; default is optimized\n");
	printf("  -T numThreads      -- default is # of CPU cores\n");
	printf("  -d                 -- dump data\n");
	printf("  -U                 -- continue on error\n");
	printf("  -p                 -- pause for MallocDebug\n");
	exit(1);
}

#pragma mark --- Private typedefs and routines ---

/* signal type */
typedef enum {
	ID_Random,
	ID_Increment
} InputData;

typedef struct {
	unsigned		log2N;
    unsigned        bits;
	InputData		inputData;
    bool            randomTau;
    bool            testReference;
    bool            compareFFT;
	bool			dumpData;
	bool			doPause;
	unsigned		numThreads;				/* 0 here means default as in API */
} TestParams;

/*
 * Given an NUFFT output and a reference output (usually NU DFT, but it can be FFT),
 * determine the error epsilon
 *
 *      max_k |X_k - X_k(exact)|)/sum_{j=0}^{D-1} |x_j|
 */
static double nuFftEpsilon(
    FFTComplex *nuFft,
    FFTComplex *refFft,
    size_t N)
{
    /*
     * Calculate both running sum of the NUFFT and the max error over
     * the entire signal.
     */
    double currMax = 0.0;           // square of max error
    double nuFftSum = 0.0;          // sum X[j] 
    PolyComplex X(nuFft);           // signal being verified
    PolyComplex ref(refFft);        // reference
    
    for(size_t dex=0; dex<N; dex++) {
        /* amplitude^2 of X[j] - ref[j] */
		double d = X.real() - ref.real();
		double deltaSquare = d * d;
		d = X.imag() - ref.imag();
		deltaSquare += (d * d);
		if(deltaSquare > currMax) {
            /* new max amplitude; store its square */
			currMax = deltaSquare;
		}
        
        /* accumulate: amplitude of current X */
        nuFftSum += sqrt((X.real() * X.real()) + (X.imag() * X.imag()));
        
        ++X;
        ++ref;
    }
    
    /* max amplitude */
    double maxAmp = sqrt(currMax);
    
    /* epsilon := max amplitude / sum */
    return maxAmp / nuFftSum;
}

#pragma mark --- run test on one size ---

static int doTest(TestParams *tp)
{
	int ourRtn = 0;
	unsigned dimArray = tp->log2N;
    
    /* One plan for MatrixFFT (optional), one for NUFFT */
	MatrixFFTPlan mfftPlan = NULL;
    NUFFTPlan nuFftPlan = NULL;
	MFFTReturn mrtn;

    FFTComplex *srcArray = NULL;            // nuFFT source, aligned
    FFTComplex *srcArrayFree = NULL;        // " " to free
    FFTComplex *dstArray = NULL;            // nuFFT destination, aligned
    FFTComplex *dstArrayFree = NULL;        // " " to free
    FFTComplex *mftArray = NULL;            // MatrixFFT src/dest, aligned
    FFTComplex *mftArrayFree = NULL;        // " " to free
    FFTComplex *nuDftSrcArray = NULL;       // nuDFT in
    FFTComplex *nuDftDstArray = NULL;       // nuDFT out
    FFTFloat *tau = NULL;                   // nuFFT tau
	double epsilon = 0.0;
    double allowError = 1.0 / ((double)(1 << tp->bits));
    uint32_t implFlag = tp->testReference ? NF_Reference : NF_Optimized;
    uint32_t planFlag = tp->compareFFT ? implFlag : (implFlag | NF_Discrete);

	size_t N = (size_t)1 << tp->log2N;
	
    /*
     * Create plans.
     */
    if(tp->compareFFT) {
        /* here's our reference implementation */
        mrtn = mfftCreatePlan(1, &dimArray, false, 0, tp->numThreads, &mfftPlan);
        if(mrtn) {
            mfftPrintErrInfo("mfftCreatePlan", mrtn);
            return -1;
        }
    }
    mrtn = nufftCreatePlan(1, &dimArray, tp->bits, planFlag, tp->numThreads, &nuFftPlan);
	if(mrtn) {
		mfftPrintErrInfo("nufftCreatePlan", mrtn);
		return -1;
	}
    
	/* 
     * Alloc & init buffers.
     */
	srcArray = fftAllocComplexArrayAlign(N, FFT_MEM_ALIGNMENT, &srcArrayFree);
	if(srcArray == NULL) {
		printf("***Malloc failure for totalSamples = %llu\n", 
			(unsigned long long)N);
		ourRtn = -1;
		goto errOut;
	} 
	dstArray = fftAllocComplexArrayAlign(N, FFT_MEM_ALIGNMENT, &dstArrayFree);
	if(dstArray == NULL) {
		printf("***Malloc failure for totalSamples = %llu\n", 
			(unsigned long long)N);
		ourRtn = -1;
		goto errOut;
	} 
    if(tp->compareFFT) {
        /* well-aligned reference destination */
        mftArray = fftAllocComplexArrayAlign(N, FFT_MEM_ALIGNMENT, &mftArrayFree);
        if(mftArray == NULL) {
            printf("***Malloc failure for totalSamples = %llu\n", 
                (unsigned long long)N);
            ourRtn = -1;
            goto errOut;
        } 
    }
    else {
        nuDftSrcArray = fftAllocComplexArray(N);
        if(nuDftSrcArray == NULL) {
            printf("***Malloc failure for totalSamples = %llu\n", 
                (unsigned long long)N);
            ourRtn = -1;
            goto errOut;
        } 
        nuDftDstArray = fftAllocComplexArray(N);
        if(nuDftDstArray == NULL) {
            printf("***Malloc failure for totalSamples = %llu\n", 
                (unsigned long long)N);
            ourRtn = -1;
            goto errOut;
        } 
    }
        
    tau = (FFTFloat *)malloc(N * sizeof(FFTFloat));
    if(tau == NULL) {
		printf("***Malloc failure for totalSamples = %llu\n", 
			(unsigned long long)N);
		ourRtn = -1;
		goto errOut;
    }
    
	switch(tp->inputData) {
		case ID_Increment:
			genIncrComplex(srcArray, N);
			break;
		case ID_Random:
			genRandComplex(srcArray, N);
			break;
	}
    if(tp->compareFFT) {
        fftCopyComplexArray(srcArray, mftArray, N);
    }
    else {
        fftCopyComplexArray(srcArray, nuDftSrcArray, N);
    }
	genConstComplex(dstArray, N, 0.0);
	
    /*
     * Initialize w array. 
     */
    if(tp->randomTau) {
        fftGenTau(tau, N);
    }
    else {
        for(size_t tauDex=0; tauDex<N; tauDex++) {
            tau[tauDex] = (FFTFloat)tauDex;
        }
    }
    
	if(tp->dumpData) {
		fftDump1DComplex("Original", srcArray, N);
		fftDump1DFloat("tau", tau, N);
	}
	
	if(tp->doPause) {
		printf("Pausing after setup, before FFTs: CR to continue: ");
		fflush(stdout);
		getchar();
	}

    /*
     * Op to test: nonuniform FFT.
     */
    mrtn = nuFftExecute(nuFftPlan, implFlag, srcArray, tau, dstArray);
	if(mrtn) {
		mfftPrintErrInfo("nuFftExecute", mrtn);
		ourRtn = -1;
		goto errOut;
	}
	if(tp->dumpData) {
		fftDump1DComplex("NUFFT forward output", dstArray, N);
	}
    
    if(tp->compareFFT) {
        /* 
         * FFT by MatrixFFT - always transpose so output is in row order.
         */
        mrtn = mfftExecute(mfftPlan, MEF_TransposeOutput, true, mftArray, mftArray);
        if(mrtn) {
            mfftPrintErrInfo("mfftExecute", mrtn);
            ourRtn = -1;
            goto errOut;
        }
        if(tp->dumpData) {
            fftDump1DComplex("MatrixFFT forward output", mftArray, N);
        }
        epsilon = nuFftEpsilon(dstArray, mftArray, N);
    }
    else {
        /* 
         * Nonuniform DFT. Slow, but it's the ultimate reference - and it's the 
         * *only* reference when using nonuniform w.
         */
        mrtn = nuFftExecute(nuFftPlan, NF_Discrete, nuDftSrcArray, tau, nuDftDstArray);
        if(mrtn) {
            mfftPrintErrInfo("nuFftExecute(NF_Discrete)", mrtn);
            ourRtn = -1;
            goto errOut;
        }
        if(tp->dumpData) {
            fftDump1DComplex("nuDFT forward output", nuDftDstArray, N);
        }
        epsilon = nuFftEpsilon(dstArray, nuDftDstArray, N);
    }

    if(epsilon > allowError) {
        printf("***Max error exceeded compared to %s\n",
            tp->compareFFT ? "FFT" : "nuDFT");
        printf("bits              : %u\n", tp->bits);
        printf("max allowed error : %.3e\n", allowError);
        printf("epsilon           : %.3e\n", epsilon);
        ourRtn = -1;
    }

	if(tp->log2N < 10) {
		printf("    2^%u ", tp->log2N);
	}
	else {
		printf("    2^%u", tp->log2N);
	}
	printf("       %.1e   %.1e\n", 
		allowError, epsilon);

errOut:
	fftFreeComplexArrayAlign(srcArray, srcArrayFree);
	fftFreeComplexArrayAlign(dstArray, dstArrayFree);
	fftFreeComplexArrayAlign(mftArray, mftArrayFree);
    fftFreeComplexArray(nuDftDstArray);
    fftFreeComplexArray(nuDftSrcArray); 
    COND_FREE(tau);
	if(mfftPlan) {
		mfftFreePlan(mfftPlan);
	}
	if(nuFftPlan) {
		nufftFreePlan(nuFftPlan);
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
	const char *dataFormStr = "Random";
	
	optStr[0] = '\0';
	memset(&tp, 0, sizeof(tp));
	tp.inputData = ID_Random;
	tp.dumpData = false;
	tp.bits = BITS_DEF;
    tp.randomTau = true;
    tp.testReference = false;
    
	int arg;
	while ((arg = getopt(argc, argv, "Us:S:b:idpuT:mrh")) != -1) {
		switch (arg) {
			case 's':
				minSize = fftParseStringRep(optarg);
				break;
			case 'S':
				maxSize = fftParseStringRep(optarg);
				break;
            case 'b':
                tp.bits = atoi(optarg);
                break;
			case 'i':
				tp.inputData = ID_Increment;
				dataFormStr = "Incrementing";
				break;
			case 'd':
				tp.dumpData = true;
				break;
			case 'p':
				tp.doPause = true;
				break;
			case 'U':
				continueOnErr = true;
				break;
			case 'T':
				tp.numThreads = atoi(optarg);
				break;
            case 'u':
                tp.randomTau = false;
                break;
            case 'm':
                tp.compareFFT = true;
                appendOptStr(optStr, "Compare to MFFT");
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
    
    if(tp.randomTau && tp.compareFFT) {
        /*
         * Allow this just as a sanity check during development - it WILL fail. It
         * should!
         */
        printf("***WARNING: randomTau and compareFFT both asserted, expect error!\n");
    }
    
    char bitStr[100];
    sprintf(bitStr, "bits=%u", tp.bits);
	appendOptStr(optStr, bitStr);
    if(tp.randomTau) {
        appendOptStr(optStr, "Random Tau");
    }
    else {
        appendOptStr(optStr, "Regular Tau");
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

	fftPrintTestBanner("1-D NUFFT", 
        tp.testReference ? "Reference NUFFT" : "Optimized NUFFT", 
		FFT_DOUBLE_PREC ? true : false, 
		dataFormStr, optStr, 0, actThreads);
	printf("\n");
	printf("   samples  | 1/2^bits | epsilon \n");
	printf(" -----------+----------+---------\n");

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
