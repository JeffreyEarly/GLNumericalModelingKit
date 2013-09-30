/*
 * Copyright (c) 2009 Apple Computer, Inc. All Rights Reserved.
 * 
 * real1D.c.cpp - simple 1-D real-to-complex FFT example.
 *
 * Created Oct. 29 2009. 
 */
 
#include <stdlib.h>
#include <strings.h>
#include <stdio.h>
#include <unistd.h>
#include <libMatrixFFT/MatrixFFT.h>
#include <libMatrixFFT/complexBufUtils.h>
#include <libMatrixFFT/fftUtils.h>

#define SAMPLES_DEF			64

static void usage(char **argv)
{
	printf("usage: %s [options]\n", argv[0]);
	printf("Options:\n");
	printf("  -s samples      -- # real samples; default %u\n", SAMPLES_DEF);
	exit(1);
}

int main(int argc, char **argv)
{
	size_t          numRealSamples = SAMPLES_DEF;
    size_t          numComplexSamples = 0;
	unsigned        n;                     // log2(numSamples)
    FFTComplex      *alignedBuf = NULL;
    FFTComplex      *freeBuf = NULL;
    size_t          dex;
    MFFTReturn      mrtn;
    MatrixFFTPlan   mfftPlan = NULL;
    int             ourRtn = -1;
    
	extern char     *optarg;
	int             arg;
    
	while ((arg = getopt(argc, argv, "s:h")) != -1) {
		switch (arg) {
			case 's':
				numRealSamples = strtol(optarg, NULL, 0);
				break;
			default:
			case 'h':
				usage(argv);
		}
	}
	if(optind != argc) {
		usage(argv);
	}
	
    /*
     * MatrixFFT requires a size which is a power of two 
     */
    if(!fftIsPowerOfTwo(numRealSamples, &n)) {
        printf("***FFT size must be a power of 2.\n");
        exit(1);
    }
    numComplexSamples = numRealSamples >> 1;
    
    /*
     * Allocate a well-aligned complex array.
     * Note we don't have to know what the complex format is here.
     */
    alignedBuf = fftAllocComplexArrayAlign(numComplexSamples, FFT_MEM_ALIGNMENT, &freeBuf);
    if(alignedBuf == NULL) {
        printf("***Memory allocation error\n");
        exit(1);
    }
    /* subsequent errors to errOut: */
    
    /*
     * Initialize alignedBuf with incrementing data.
     * -- Two different code paths here, for the two complex formats.
     * -- Using the macros at the end of complexBufUtils.h would make this simpler;
     *    we're illustrating what's actually going on here.
     * -- Note that this code is independent of the precision for which
     *    libMatrixFFT is currently configured since we use the FFTFloat 
     *    type. 
     */
    #if     FFT_SPLIT_COMPLEX
        for(dex=0; dex<numRealSamples; ) {
            size_t bufDex = dex/2;
            alignedBuf->real[bufDex] = (FFTFloat)dex++;
            alignedBuf->imag[bufDex] = (FFTFloat)dex++;
        }
    #else   /* !FFT_SPLIT_COMPLEX */
        for(dex=0; dex<numRealSamples; ) {
            size_t bufDex = dex/2;
            alignedBuf[bufDex].real = (FFTFloat)dex++;
            alignedBuf[bufDex].imag = (FFTFloat)dex++;
        }
    #endif
    
    fftDump1DComplex("Original data", alignedBuf, numComplexSamples);
    
    /*
     * Set up a MatrixFFTPlan.
     * This can be used for multiple FFTs as long as they are the same size and
     * type (real/complex). 
     */
    mrtn = mfftCreatePlan(1,        // dims
                &n,
                true,               // real
                0,                  // flags
                0,                  // numThreads - 0 = default
                &mfftPlan);
    if(mrtn) {
        mfftPrintErrInfo("mfftCreatePlan", mrtn);
        ourRtn = -1;
        goto errOut;
    }
    
    /* 
     * Execute the forward FFT in-place.
     * MEF_TransposeOutput --> output in row order
     */
    mrtn = mfftExecute(mfftPlan,
        MEF_TransposeOutput,
        true,                   // forward
        alignedBuf, alignedBuf);
    if(mrtn) {
        mfftPrintErrInfo("mfftExecute(forward)", mrtn);
        ourRtn = -1;
        goto errOut;
    }

    fftDump1DComplex("Forward FFT output", alignedBuf, numComplexSamples);
 
    /* 
     * Execute the inverse FFT in-place.
     * MEF_TransposeInput --> input is currently in row order
     * MEF_NormOutput     --> normalize output
     */
    mrtn = mfftExecute(mfftPlan,
        MEF_TransposeInput | MEF_NormOutput,
        false,                   // forward
        alignedBuf, alignedBuf);
    if(mrtn) {
        mfftPrintErrInfo("mfftExecute(inverse)", mrtn);
        ourRtn = -1;
        goto errOut;
    }

    fftDump1DComplex("Inverse FFT output", alignedBuf, numComplexSamples);

errOut:
    if(alignedBuf != NULL) {
        fftFreeComplexArrayAlign(alignedBuf, freeBuf);
    }
    if(mfftPlan != NULL) {
        mfftFreePlan(mfftPlan);
    }
	return ourRtn;
}
