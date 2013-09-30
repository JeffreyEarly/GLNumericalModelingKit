/*
 * Copyright (c) 2009 Apple Computer, Inc. All Rights Reserved.
 * 
 * complex1D.c.cpp - simple 2-D complex FFT example.
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

#define ROWS_DEF			4
#define COLS_DEF            4

static void usage(char **argv)
{
	printf("usage: %s [options]\n", argv[0]);
	printf("Options:\n");
	printf("  -r rows      -- # rows; default %u\n", ROWS_DEF);
	printf("  -c cols      -- # columns; default %u\n", COLS_DEF);
	exit(1);
}

int main(int argc, char **argv)
{
    size_t          numRows = ROWS_DEF;
    size_t          numCols = COLS_DEF;
	size_t          totalSamples = 0;
	unsigned        n[2];                     // log2(numRows), log2(numCols)
    FFTComplex      *alignedBuf = NULL;
    FFTComplex      *freeBuf = NULL;
    size_t          dex;
    MFFTReturn      mrtn;
    MatrixFFTPlan   mfftPlan = NULL;
    int             ourRtn = -1;
    
	extern char     *optarg;
	int             arg;
    
	while ((arg = getopt(argc, argv, "r:c::h")) != -1) {
		switch (arg) {
			case 'r':
				numRows = strtol(optarg, NULL, 0);
				break;
			case 'c':
				numCols = strtol(optarg, NULL, 0);
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
     * MatrixFFT requires sizes which are powers of two 
     */
    if(!fftIsPowerOfTwo(numRows, &n[0])) {
        printf("***FFT size must be a power of 2.\n");
        exit(1);
    }
    if(!fftIsPowerOfTwo(numCols, &n[1])) {
        printf("***FFT size must be a power of 2.\n");
        exit(1);
    }
    
    /*
     * Allocate a well-aligned complex array.
     * Note we don't have to know what the complex format is here.
     */
    totalSamples = numRows * numCols;
    alignedBuf = fftAllocComplexArrayAlign(totalSamples, FFT_MEM_ALIGNMENT, &freeBuf);
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
        for(dex=0; dex<totalSamples*2; ) {
            size_t bufDex = dex/2;
            alignedBuf->real[bufDex] = (FFTFloat)dex++;
            alignedBuf->imag[bufDex] = (FFTFloat)dex++;
        }
    #else   /* !FFT_SPLIT_COMPLEX */
        for(dex=0; dex<totalSamples*2; ) {
            size_t bufDex = dex/2;
            alignedBuf[bufDex].real = (FFTFloat)dex++;
            alignedBuf[bufDex].imag = (FFTFloat)dex++;
        }
    #endif
    
    fftDump2DComplex("Original data", alignedBuf, numRows, numCols);
    
    /*
     * Set up a MatrixFFTPlan.
     * This can be used for multiple FFTs as long as they are the same size and
     * type (real/complex). 
     */
    mrtn = mfftCreatePlan(2,        // dims
                n,
                false,              // real
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

    fftDump2DComplex("Forward FFT output", alignedBuf, numRows, numCols);
 
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

    fftDump2DComplex("Inverse FFT output", alignedBuf, numRows, numCols);

errOut:
    if(alignedBuf != NULL) {
        fftFreeComplexArrayAlign(alignedBuf, freeBuf);
    }
    if(mfftPlan != NULL) {
        mfftFreePlan(mfftPlan);
    }
	return ourRtn;
}
