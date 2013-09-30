/*	File:  mfftAutoConf.cpp

	Description:
		Platform-dependent configuration utility for MatrixFFT library.
	
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
 * mfftAutoConf.cpp - MatrixFFT configuration utility.
 *
 * Created Nov. 6 2009. 
 *
 * See Configuration.txt in the root directory of the MatrixFFT project for 
 * information on using this tool.
 */
 
#include <stdlib.h>
#include <strings.h>
#include <stdio.h>
#include <unistd.h>
#include <libMatrixFFT/MatrixFFT.h>
#include <libMatrixFFT/fftUtils.h>
#include <libMatrixFFT/complexBufUtils.h>
#include <libMatrixFFT/src/fftPlatformConf.h>   /* private SPI */
#include <CoreFoundation/CoreFoundation.h>
#include "getHwInfo.h"

#define MIN_SIZE_DEF        128
#define MIN_OFFSET          (-4)
#define MAX_OFFSET          3
#define LOOPS_TINY          2
#define LOOPS_QUICK         4
#define LOOPS_NORM          16
#define LOOPS_LARGE         64

/*
 * This defines the indent of the output we generate 
 */
#define INDENT  "      "


static void usage(char **argv)
{
	printf("usage: %s [options]\n", argv[0]);
	printf("Options:\n");
    printf("  -n minSize      -- minimum 1D real FFT size; K, M, 2^n OK; default %u\n", MIN_SIZE_DEF);
    printf("  -x maxSize      -- maximum 1D real FFT size; K, M, 2^n OK; default per HW config\n");
    printf("  -q              -- quick (Small loop counts and sizes)\n");
    printf("  -v              -- verbose\n");
    printf("  -V              -- really verbose\n");
    printf("  -1              -- (one) 1-dimension only\n");
    printf("  -2              -- 2-dimension only\n");
    printf("  -c              -- complex only\n");
    printf("  -r              -- real only\n");
    printf("  -M              -- min MatrixFFT sizes only\n");
    printf("  -i              -- just print CPU/memory info as C comment\n");
    
	exit(1);
}

#pragma mark --- Utility functions ---

/* 
 * Determine the max size for a 1-D real FFT for this platform.
 * Returns nonzero on error.
 * What we do is:
 *  -- Divide RAM size by two
 *  -- Find the largest power-of-two number of FFTFloats that fit in that size
 */
static int maxFftSize(
    uint64_t    ramSize,
    size_t      *maxFftSize)        /* RETURNED */
{
    if((sizeof(int *) == 4) && (ramSize > 0x100000000ULL)) {
        printf("**************************************************************************\n");
        printf("***WARNING: MatrixFFT is currently built as a 32 bit binary.\n");
        printf("            There is more than 4 GB of RAM on this machine, but this\n");
        printf("            configuration will be limited to using only 4 GB.\n");
        printf("**************************************************************************\n");
        ramSize = 0x100000000ULL;
    }
    size_t availRAM = ramSize >> 1;
    size_t maxFloats = availRAM / sizeof(FFTFloat);
    
    unsigned exp;
    fftIsPowerOfTwo(maxFloats, &exp);
    
    *maxFftSize = ((size_t)1) << exp;
    return 0;
}

static unsigned calcLoops(
    bool quick,
    unsigned n,
    bool twoDim)    /* 2-D runs faster, use more loops */
{
    if(n > 27) {
        /* enormous, just one round-trip will do it */
        return 1;
    }
    
    unsigned loops = 0;
    if(n > 24) {
        loops = LOOPS_TINY;
    }
    else if(n > 20) {
        loops = LOOPS_QUICK;
    }
    else if(n > 16) {
        loops = LOOPS_NORM;
    }
    else {
        loops = LOOPS_LARGE;
    }
    if(twoDim && (n == 26)) {
        loops <<= 1;
    }
    if(quick) {
        loops >>= 1;
    }   
    return loops;
}

static void printStripeSizeTable(
    const char  *arrayName,      // e.g. "fftStripeSize2DReal"
    size_t      *sizeArray,
    unsigned    numSizes)
{
    printf("%sstatic size_t %s[] = \n", INDENT, arrayName);
    printf("%s       { ", INDENT);
    for(unsigned dex=0; dex<numSizes; dex++) {
        printf("%llu", (unsigned long long)(sizeArray[dex]));
        if(dex < (numSizes - 1)) {
            printf(", ");
        }
    }
    printf("};\n");
}

static void showInfo(HWInfo *hwInfo)
{
    printf("/*\n");
    printf(" * CPU                     : %s\n", hwInfo->cpuInfo);

    char *ramStr = fftStringRep(hwInfo->RAMSize);
    printf(" * RAM                     : %sB\n", ramStr);
    free(ramStr);
    
    printf(" * Physical CPUs           : %u\n", (unsigned)hwInfo->numPhysicalCPUs);
    printf(" * Logical  CPUs           : %u\n", (unsigned)hwInfo->numVirtualCPUs);
    
    printf(" * L1 data cache size      : %llu\n", (unsigned long long)hwInfo->cacheInfo[0].cacheSize);
    for(unsigned dex=2; dex<=MAX_NUM_CACHES; dex++) {
        uint64_t cacheSize = hwInfo->cacheInfo[dex-1].cacheSize;
        if(cacheSize != 0) {
            printf(" * L%u cache size           : %llu\n", dex, (unsigned long long)cacheSize);
        }
    }
    
    for(unsigned dex=1; dex<=MAX_NUM_CACHES; dex++) {
        uint32_t procsPerCache = hwInfo->cacheInfo[dex-1].processorsPerCache;
        if(procsPerCache != 0) {
            printf(" * L%u processors per cache : %lu\n", dex, (unsigned long)procsPerCache);
        }
    }
    printf(" */\n\n");
}

/*
 * Print a preprocessor expression denoting the current configuration.
 */
static void configAsCppString()
{
    printf("((FFT_DOUBLE_PREC == %u) && (FFT_SPLIT_COMPLEX == %u))",
        FFT_DOUBLE_PREC, FFT_SPLIT_COMPLEX);
}

#pragma mark --- Rectangle offset and stripe size calclulation ---

/*
 * Basic FFT parameters
 */
typedef struct {
    unsigned    dims;
    unsigned    n[2];
    bool        isReal;
    FFTComplex  *buf;
    unsigned    loops;
} FFTParams;

typedef enum {
    V_Normal,
    V_Verbose,
    V_Noisy
} Verbosity;

/*
 * -- Set up a plan for specified FFT type and column stripe size
 * -- time the execution of 'loops' pairs of FFT;
 * -- if elapsed time < current fastest, replace fasted time & stripe size with current
 */
static int doFftForStripeSize(
    FFTParams *params,
    size_t newStripeSize,
    CFAbsoluteTime *fastest,        // IN/OUT
    size_t *fastestStripeSize,      // IN/OUT
    CFAbsoluteTime *slowest,        // IN/OUT
    Verbosity verbosity)
{
    uint32_t fwdFlags = 0;
    uint32_t invFlags = MEF_NormOutput;
    MatrixFFTPlan mfftPlan;
    MFFTReturn mrtn;
    
    mfftSetColumnStripeSize(params->isReal, newStripeSize);
    mrtn = mfftCreatePlan(params->dims, params->n, params->isReal, 0, 0, &mfftPlan);
    if(mrtn) {
        mfftPrintErrInfo("mfftCreatePlan", mrtn);
        return -1;
    }
    
    CFAbsoluteTime startTime, endTime, elapsed;
    
    if(verbosity == V_Noisy) {
        printf("+++++++++ stripeSize %llu\n", (unsigned long long)newStripeSize);
    }
    startTime = CFAbsoluteTimeGetCurrent();
    for(unsigned loop=0; loop<params->loops; loop++) {
        mrtn = mfftExecute(mfftPlan, fwdFlags, true, params->buf, params->buf);
        if(mrtn) {
            mfftPrintErrInfo("mfftExecute(fwd)", mrtn);
            break;
        }
        mrtn = mfftExecute(mfftPlan, invFlags, false, params->buf, params->buf);
        if(mrtn) {
            mfftPrintErrInfo("mfftExecute(inv)", mrtn);
            break;
        }
    }
    endTime = CFAbsoluteTimeGetCurrent();
    elapsed = endTime - startTime;
    if(verbosity == V_Noisy) {
        printf("          elapsed %.2e\n", elapsed);
    }
    if((*fastest == 0.0) ||     // i.e. nothing measured yet
       (*fastest > elapsed)) {
        *fastest = elapsed;
        *fastestStripeSize = newStripeSize;
    }
    if(elapsed > *slowest) {
        *slowest = elapsed;
    }
    mfftFreePlan(mfftPlan);
    return (mrtn == MR_Success) ? 0 : 1;
}

static int doStripeSizes(
    FFTParams *params,
    CFAbsoluteTime *fastest,        // IN/OUT
    size_t *fastestStripeSize,      // IN/OUT
    CFAbsoluteTime *slowest,        // IN/OUT
    Verbosity verbosity)
{
    int ourRtn = 0;
    for(unsigned stripeSize=1; stripeSize<16; stripeSize++) {
        ourRtn = doFftForStripeSize(params, stripeSize, fastest, 
            fastestStripeSize, slowest, verbosity);
        if(ourRtn) {
            return ourRtn;
        }
    }
    for(unsigned stripeSize=16; 
                 stripeSize<=(8*FFT_SUBMATRIX_ATOM); 
                 stripeSize+=FFT_SUBMATRIX_ATOM) {
        ourRtn = doFftForStripeSize(params, stripeSize, fastest, 
            fastestStripeSize, slowest, verbosity);
        if(ourRtn) {
            return ourRtn;
        }
    }
    return 0;
}

/* generate rectangle offsets and stripeSize1DComplex[] */
static int gen1DComplex(
    unsigned log2MinSize, 
    unsigned log2MaxSize, 
    int *rectOffset,                // RETURNED
    size_t *stripeSize1DComplex,    // RETURNED
    bool quick, 
    Verbosity verbosity)
{
    /*
     * One-time only allocation of all the memory we need.
     */
    size_t maxComplex = ((size_t)1) << log2MaxSize;
    FFTComplex *alignBuf = NULL;
    FFTComplex *freeBuf = NULL;
    alignBuf = fftAllocComplexArrayAlign(maxComplex, FFT_MEM_ALIGNMENT, &freeBuf);
    if(alignBuf == NULL) {
        printf("***gen1DComplex: malloc failure; maxComplex %llu\n",
            (unsigned long long)maxComplex);
        return -1;
    }
    
    FFTParams params;
    params.dims = 1;
    params.isReal = false;
    params.buf = alignBuf;
    
    int ourRtn = 0;
    
    for(unsigned n=log2MinSize; n<=log2MaxSize; n++) {
        printf("+++ timing 1D complex N=2^%u\n", n);
        
        params.n[0] = n;
        params.loops = calcLoops(quick, n, false);
        
        CFAbsoluteTime totalFastest = 0.0;
        CFAbsoluteTime slowestTime = 0.0;
        int fastestOffset = 0;
        size_t fastestStripeSize = 0;
        
        /* 
         * Re-init here to make sure memory is resident
         */
        size_t numComplex = ((size_t)1) << n;
        genRandComplex(alignBuf, numComplex);
        
        /*
         * Clamp offset for small signals 
         */
        int actMinOff = MIN_OFFSET;
        int actMaxOff = MAX_OFFSET;
        unsigned nOver2Minus1 = (n >> 1) - 1;
        if(abs(actMinOff) >= nOver2Minus1) {
            actMinOff = -((int)nOver2Minus1) + 1;
        }
        if(actMaxOff >= nOver2Minus1) {
            actMaxOff = nOver2Minus1 - 1;
        }
        for(int offset=actMinOff; offset<=actMaxOff; offset++) {
            if(verbosity != V_Normal) {
                printf("++++++ timing offset=%d\n", offset);
            }   
            mfftSetRectangleOffset(offset, true);
            CFAbsoluteTime fastestForOffset = 0.0;
            size_t stripeSizeForOffset = 0;
            
            ourRtn = doStripeSizes(&params, &fastestForOffset, &stripeSizeForOffset,
                &slowestTime, verbosity);
            if(ourRtn) {
                break;
            }
            
            if(verbosity == V_Noisy) {
                printf("++++++ offset=%d best time %.2e at size %llu\n", 
                    offset, fastestForOffset, (unsigned long long)stripeSizeForOffset);
            }   
            if((totalFastest == 0.0) || // first time thru
               (totalFastest > fastestForOffset)) {
                /* 
                 * New fastest - i.e. the fastest stripe size for THIS offset is the 
                 * fastest we've seen so far 
                 */
                totalFastest = fastestForOffset;
                fastestStripeSize = stripeSizeForOffset;
                fastestOffset = offset;
            }
        }
        if(ourRtn) {
            break;
        }
        
        if(verbosity != V_Normal) {
            printf("+++ 1DComplex(%u): offset %d  stripeSize %llu  fastest %.2e  slowest %.2e\n",
                n, fastestOffset, (unsigned long long)fastestStripeSize,
                totalFastest, slowestTime);
        }
        rectOffset[n-log2MinSize] = fastestOffset;
        stripeSize1DComplex[n-log2MinSize] = fastestStripeSize;
    }
    
    fftFreeComplexArrayAlign(alignBuf, freeBuf);
    return ourRtn;
}

#define HALF_RECT_ENABLE    0

static int gen2DComplex(
    unsigned log2MinSize,           // log2(numRows)
    unsigned log2MaxSize, 
    size_t *stripeSize2DComplex,    // RETURNED
    bool quick, 
    Verbosity verbosity)
{   
    /*
     * One-time only allocation of all the memory we need.
     */
    size_t maxComplex = ((size_t)1) << (log2MaxSize * 2);
    FFTComplex *alignBuf = NULL;
    FFTComplex *freeBuf = NULL;
    alignBuf = fftAllocComplexArrayAlign(maxComplex, FFT_MEM_ALIGNMENT, &freeBuf);
    if(alignBuf == NULL) {
        printf("***gen2DComplex: malloc failure; maxComplex %llu\n",
            (unsigned long long)maxComplex);
        return -1;
    }
    
    FFTParams params;
    params.dims = 2;
    params.isReal = false;
    params.buf = alignBuf;
    
    int ourRtn = 0;
    
    for(unsigned nNumRows=log2MinSize; nNumRows<=log2MaxSize; nNumRows++) {
        printf("+++ timing 2D complex numRows=2^%u\n", nNumRows);
        
        params.n[0] = nNumRows;
        
        /* 
         * NOTE WELL: a 2-D complex FFT with numRows==numCols, with a sufficiently
         * large signal, is impervious to modification of column strips size because
         * the FFT_2D_SQUARE_ENABLE mechanism (see fft2DComplex.cpp) does an 
         * in-place (square) transpose of the whole signal instead of column peels.
         * So we *never* want a square signal here. 
         */
        /* small signal, square */
        unsigned nNumCols = nNumRows-1;
        #if HALF_RECT_ENABLE
        if(nNumRows > 12) {
            /* 
             * We really don't need to test numCols == numRows - everything we're
             * interested in is solely determined by the column size. Let's cut
             * the total signal size in half. 
             */
            nNumCols = nNumRows >> 1;
        }
        #endif  /* HALF_RECT_ENABLE */
        
        params.n[1] = nNumCols;
        unsigned totalN = nNumRows + nNumCols;
        params.loops = calcLoops(quick, totalN, true);
        
        /* 
         * Re-init here to make sure memory is resident
         */
        size_t numComplex = ((size_t)1) << totalN;
        genRandComplex(alignBuf, numComplex);
        
        CFAbsoluteTime fastestTime = 0.0;
        CFAbsoluteTime slowestTime = 0.0;
        size_t fastestStripeSize = 0;
        
        ourRtn = doStripeSizes(&params, &fastestTime, &fastestStripeSize, &slowestTime, verbosity);
        if(ourRtn) {
            break;
        }
            
        if(verbosity != V_Normal) {
            printf("+++ 2DComplex(%u): stripeSize %llu  fastest %.2e  slowest %.2e\n",
                nNumRows, (unsigned long long)fastestStripeSize,
                fastestTime, slowestTime);
        }
        stripeSize2DComplex[nNumRows-log2MinSize] = fastestStripeSize;
    }
    
    fftFreeComplexArrayAlign(alignBuf, freeBuf);
    return ourRtn;
}

static int gen2DReal(
    unsigned log2MinSize,           // log2(numRows)
    unsigned log2MaxSize, 
    size_t *stripeSize2DReal,       // RETURNED
    bool quick, 
    Verbosity verbosity)
{   
    /*
     * One-time only allocation of all the memory we need.
     */
    size_t maxComplex = ((size_t)1) << ((log2MaxSize * 2) - 1);
    FFTComplex *alignBuf = NULL;
    FFTComplex *freeBuf = NULL;
    alignBuf = fftAllocComplexArrayAlign(maxComplex, FFT_MEM_ALIGNMENT, &freeBuf);
    if(alignBuf == NULL) {
        printf("***gen2DReal: malloc failure; maxComplex %llu\n",
            (unsigned long long)maxComplex);
        return -1;
    }
    
    FFTParams params;
    params.dims = 2;
    params.isReal = true;
    params.buf = alignBuf;
    
    int ourRtn = 0;
    
    for(unsigned nNumRows=log2MinSize; nNumRows<=log2MaxSize; nNumRows++) {
        printf("+++ timing 2D real    numRows=2^%u\n", nNumRows);
        
        params.n[0] = nNumRows;
        
        /* small signal, square */
        unsigned nNumCols = nNumRows;
        
        #if HALF_RECT_ENABLE
        if(nNumRows > 12) {
            /* 
             * We really don't need to test numCols == numRows - everything we're
             * interested in is solely determined by the column size. Let's cut
             * the total signal size in half. 
             */
            nNumCols = nNumRows >> 1;
        }
        #endif/* HALF_RECT_ENABLE */
        
        params.n[1] = nNumCols;
        unsigned totalN = nNumRows + nNumCols;
        params.loops = calcLoops(quick, totalN, true);
        
        /* 
         * Re-init here to make sure memory is resident
         */
        size_t numComplex = ((size_t)1) << (totalN - 1);
        genRandComplex(alignBuf, numComplex);
        
        CFAbsoluteTime fastestTime = 0.0;
        CFAbsoluteTime slowestTime = 0.0;
        size_t fastestStripeSize = 0;
        
        ourRtn = doStripeSizes(&params, &fastestTime, &fastestStripeSize, &slowestTime, verbosity);
        if(ourRtn) {
            break;
        }
            
        if(verbosity != V_Normal) {
            printf("+++ 2DReal(%u): stripeSize %llu  fastest %.2e  slowest %.2e\n",
                nNumRows, (unsigned long long)fastestStripeSize,
                fastestTime, slowestTime);
        }
        stripeSize2DReal[nNumRows-log2MinSize] = fastestStripeSize;
    }
    
    fftFreeComplexArrayAlign(alignBuf, freeBuf);
    return ourRtn;
}

#pragma mark --- calculation of min sizes for MatrixFFT ---
 
/* 
 * These constants define the log2 of the starting sizes for
 * the 4 searches. 2-D values are log2(numRows).
 */
#define COMPLEX_1D_START    5
#define COMPLEX_2D_START    4
#define REAL_1D_START       4
#define REAL_2D_START       4

/* 
 * For specified FFT:
 * -- set real and complex stripe sizes if provided;
 * -- allocate & init complex buffer;
 * -- set MF_ForceVdsp;
 * -- time 'loops' round-trip FFTs;
 * -- set MF_ForceMfft;
 * -- time 'loops' round-trip FFTs;
 * -- if MF_ForceMfft is faster, *mfftFaster is true, else it's false;
 */
static int compareFFTs(
    FFTParams *params,
    size_t complexStripeSize,
    size_t realStripeSize,
    int rectOffset,                 // valid iff rectOffsetSpecd
    bool rectOffsetSpecd,
    Verbosity verbosity,
    bool *mfftFaster)               // RETURNED
{
    /* Optional pre-calculated stripe and rectangle offset specification */
    if(complexStripeSize != 0) {
        mfftSetColumnStripeSize(false, complexStripeSize);
    }
    if(realStripeSize != 0) {
        mfftSetColumnStripeSize(true, realStripeSize);
    }
    mfftSetRectangleOffset(rectOffset, rectOffsetSpecd);
    
    if(verbosity == V_Noisy) {
        printf("+++++++++ MF_ForceVdsp(%u) complexStripe %u realStripe %u rectOff %d rectSpecd %s\n", 
            params->n[0] + params->n[1], 
            (unsigned)complexStripeSize,
            (unsigned)realStripeSize,
            rectOffset, rectOffsetSpecd ? "true" : "false");
    }

    /* allocate buffers for this op */
    size_t numComplex = ((size_t)1) << params->n[0];
    if(params->dims == 2) {
        numComplex *= (((size_t)1) << params->n[1]);
    }
    if(params->isReal) {
        numComplex >>= 1;
    }
    FFTComplex *alignBuf = NULL;
    FFTComplex *freeBuf = NULL;
    alignBuf = fftAllocComplexArrayAlign(numComplex, FFT_MEM_ALIGNMENT, &freeBuf);
    if(alignBuf == NULL) {
        printf("***compareFFTs: malloc failure; numComplex %llu\n",
            (unsigned long long)numComplex);
        return -1;
    }
    params->buf = alignBuf;
    genRandComplex(alignBuf, numComplex);
    
    /* subsequent errors to errOut: */
    int ourRtn = 0;
    
    CFAbsoluteTime vdspTime = 0.0;
    CFAbsoluteTime mfftTime = 0.0;
    CFAbsoluteTime startTime;
    CFAbsoluteTime endTime;
    
    uint32_t fwdFlags = 0;
    uint32_t invFlags = MEF_NormOutput;
    MatrixFFTPlan mfftPlan;
    MFFTReturn mrtn;

    /*
     * Time raw vDSP mode
     */
    mfftSetForceVdsp(MF_ForceVdsp);
    mrtn = mfftCreatePlan(params->dims, params->n, params->isReal, 0, 0, &mfftPlan);
    if(mrtn) {
        mfftPrintErrInfo("mfftCreatePlan", mrtn);
        ourRtn = -1;
        goto errOut;
    }
    startTime = CFAbsoluteTimeGetCurrent();
    for(unsigned loop=0; loop<params->loops; loop++) {
        mrtn = mfftExecute(mfftPlan, fwdFlags, true, params->buf, params->buf);
        if(mrtn) {
            mfftPrintErrInfo("mfftExecute(fwd)", mrtn);
            ourRtn = -1;
            goto errOut;
        }
        mrtn = mfftExecute(mfftPlan, invFlags, false, params->buf, params->buf);
        if(mrtn) {
            mfftPrintErrInfo("mfftExecute(inv)", mrtn);
            ourRtn = -1;
            goto errOut;
        }
    }
    endTime = CFAbsoluteTimeGetCurrent();
    vdspTime = endTime - startTime;
    if(verbosity == V_Noisy) {
        printf("          MF_ForceVdsp elapsed %.2e\n", vdspTime);
    }
    mfftFreePlan(mfftPlan);
    
    /*
     * Time MatrixFFT-only mode
     */
    mfftSetForceVdsp(MF_ForceMfft);
    mrtn = mfftCreatePlan(params->dims, params->n, params->isReal, 0, 0, &mfftPlan);
    if(mrtn) {
        mfftPrintErrInfo("mfftCreatePlan", mrtn);
        ourRtn = -1;
        goto errOut;
    }
    startTime = CFAbsoluteTimeGetCurrent();
    for(unsigned loop=0; loop<params->loops; loop++) {
        mrtn = mfftExecute(mfftPlan, fwdFlags, true, params->buf, params->buf);
        if(mrtn) {
            mfftPrintErrInfo("mfftExecute(fwd)", mrtn);
            ourRtn = -1;
            goto errOut;
        }
        mrtn = mfftExecute(mfftPlan, invFlags, false, params->buf, params->buf);
        if(mrtn) {
            mfftPrintErrInfo("mfftExecute(inv)", mrtn);
            ourRtn = -1;
            goto errOut;
        }
    }
    endTime = CFAbsoluteTimeGetCurrent();
    mfftTime = endTime - startTime;
    if(verbosity == V_Noisy) {
        printf("          MF_ForceMfft elapsed %.2e\n", mfftTime);
    }
    mfftFreePlan(mfftPlan);
    
    *mfftFaster = (mfftTime < vdspTime);
errOut:
    fftFreeComplexArrayAlign(alignBuf, freeBuf);
    return ourRtn;
}

/* 
 * Generate MFFT_1D_COMPLEX_MIN_SIZE using the arrays generated 
 * in gen1DComplex(). If this is a standalone MIN_SIZE run, pass
 * NULL for the generated rectOffset and stripeSize1DComplex arrays,
 * and we'll just use the current configuration for those.
 */
static int gen1DComplexMinSize(
    /* 
     * These two determine the bounds of the complex rectOffset and 
     * stripeSize arrays. Typically we start *below* log2MinSize.
     */
    unsigned log2MinSize,      
    unsigned log2MaxSize, 
    const int *rectOffset, 
    const size_t *stripeSize1DComplex,
    unsigned *min1DComplexSize,    // RETURNED
    Verbosity verbosity)
{
    bool rectOffSpecd = (rectOffset != NULL);
    FFTParams params;
    
    params.dims = 1;
    params.isReal = false;
    params.n[1] = 0;
    
    for(unsigned n=COMPLEX_1D_START; n<=log2MaxSize; n++) {
        int rectOff = 0;
        size_t complexStripeSize = 0;
        if(rectOffSpecd) {
            if(n < log2MinSize) {
                rectOff = rectOffset[0];
            }
            else {
                rectOff = rectOffset[n - log2MinSize];
            }
        }
        if(stripeSize1DComplex != NULL) {
            if(n < log2MinSize) {
                complexStripeSize = stripeSize1DComplex[0];
            }
            else {
                complexStripeSize = stripeSize1DComplex[n - log2MinSize];
            }
        }
        bool mfftFaster = false;
        int rtn;
        
        params.n[0] = n;
        params.loops = calcLoops(false, n, false);
        rtn = compareFFTs(&params,
            complexStripeSize,
            0,              // realStripeSize
            rectOff, rectOffSpecd,
            verbosity,
            &mfftFaster);
        if(rtn) {
            return rtn;
        }
        if(mfftFaster) {
            if(verbosity != V_Normal) {
                printf("++++++ 1D complex minSize=%u\n", n);
            }  
            *min1DComplexSize = n;
            return 0;
        }
    }
    
    /*
     * I'd be surprised if we ever reach here, other than if 
     * we're running with a very restricted max size
     */
    if(verbosity != V_Normal) {
        printf("++++++ 1D complex loop exhausted; minSize=%u\n", log2MaxSize);
    }  
    *min1DComplexSize = log2MaxSize;
    return 0;
}

/* 
 * Generate MFFT_2D_COMPLEX_MIN_SIZE using the arrays generated 
 * in gen2DComplex(). If this is a standalone MIN_SIZE run, pass
 * NULL for the generated stripeSize2DComplex array, and we'll just 
 * use the current configuration for it.
 */
static int gen2DComplexMinSize(
    /* 
     * These two determine the bounds of the complex stripeSize 
     * array. Typically we start *below* log2MinSize.
     * The sizes are numRows, not the total signal size.
     */
    unsigned log2MinSize,      
    unsigned log2MaxSize, 
    const size_t *stripeSize2DComplex,
    unsigned *min2DComplexSize,    // RETURNED
    Verbosity verbosity)
{
    FFTParams params;
    
    params.dims = 2;
    params.isReal = false;
    
    unsigned nRow = COMPLEX_2D_START;
    unsigned nCol = COMPLEX_2D_START;
    
    do { 
        unsigned nTotal = nRow + nCol;
        size_t complexStripeSize = 0;
        if(stripeSize2DComplex != NULL) {
            if(nRow < log2MinSize) {
                complexStripeSize = stripeSize2DComplex[0];
            }
            else {
                complexStripeSize = stripeSize2DComplex[nRow - log2MinSize];
            }
        }
        bool mfftFaster = false;
        int rtn;
        
        params.n[0] = nRow;
        params.n[1] = nCol;
        params.loops = calcLoops(false, nTotal, true);
        rtn = compareFFTs(&params,
            complexStripeSize,
            0,              // realStripeSize
            0, false,       // rectOff, rectOffSpecd
            verbosity,
            &mfftFaster);
        if(rtn) {
            return rtn;
        }
        if(mfftFaster) {
            if(verbosity != V_Normal) {
                printf("++++++ 2D complex minSize=%u\n", nTotal);
            }  
            *min2DComplexSize = nTotal;
            return 0;
        }
        
        if(nRow == nCol) {
            nRow++;
        }
        else {
            nCol++;
        }
    } while(nRow <=log2MaxSize);
    
    /*
     * I'd be surprised if we ever reach here, other than if 
     * we're running with a very restricted max size
     */
    if(verbosity != V_Normal) {
        printf("++++++ 2D complex loop exhausted; minSize=%u\n", 2 * log2MaxSize);
    }  
    *min2DComplexSize = 2* log2MaxSize;
    return 0;
}

/* 
 * Generate MFFT_1D_REAL_MIN_SIZE using the arrays generated 
 * in gen1DComplex(). If this is a standalone MIN_SIZE run, pass
 * NULL for the generated rectOffset and stripeSize1DComplex arrays,
 * and we'll just use the current configuration for those.
 */
static int gen1DRealMinSize(
    unsigned log2MaxSizeReal,
    /* 
     * These two determine the bounds of the complex rectOffset and 
     * stripeSize arrays. Typically we start *below* log2MinSize.
     */
    unsigned log2MinSizeCompl,      
    unsigned log2MaxSizeCompl, 
    const int *rectOffset, 
    const size_t *stripeSize1DComplex,
    unsigned *min1DRealSize,    // RETURNED
    Verbosity verbosity)
{
    bool rectOffSpecd = (rectOffset != NULL);
    FFTParams params;
    
    params.dims = 1;
    params.isReal = true;
    params.n[1] = 0;
    
    for(unsigned n=REAL_1D_START; n<=log2MaxSizeReal; n++) {
        unsigned nComplex = n - 1;
        int rectOff = 0;
        size_t complexStripeSize = 0;
        if(rectOffSpecd) {
            if(nComplex < log2MinSizeCompl) {
                rectOff = rectOffset[0];
            }
            else {
                rectOff = rectOffset[nComplex - log2MinSizeCompl];
            }
        }
        if(stripeSize1DComplex != NULL) {
            if(nComplex < log2MinSizeCompl) {
                complexStripeSize = stripeSize1DComplex[0];
            }
            else {
                complexStripeSize = stripeSize1DComplex[nComplex - log2MinSizeCompl];
            }
        }
        bool mfftFaster = false;
        int rtn;
        
        params.n[0] = n;
        params.loops = calcLoops(false, n, false);
        rtn = compareFFTs(&params,
            complexStripeSize,
            0,              // realStripeSize
            rectOff, rectOffSpecd,
            verbosity,
            &mfftFaster);
        if(rtn) {
            return rtn;
        }
        if(mfftFaster) {
            if(verbosity != V_Normal) {
                printf("++++++ 1D real minSize=%u\n", n);
            }  
            *min1DRealSize = n;
            return 0;
        }
    }
    
    /*
     * I'd be surprised if we ever reach here, other than if 
     * we're running with a very restricted max size
     */
    if(verbosity != V_Normal) {
        printf("++++++ 1D real loop exhausted; minSize=%u\n", log2MaxSizeReal);
    }  
    *min1DRealSize = log2MaxSizeReal;
    return 0;
}

/* 
 * Generate MFFT_2D_REAL_MIN_SIZE using the arrays generated 
 * in gen2DReal(). If this is a standalone MIN_SIZE run, pass
 * NULL for the generated stripeSize2DReal array, and we'll just 
 * use the current configuration for it.
 */
static int gen2DRealMinSize(
    /* 
     * These two determine the bounds of the stripeSize 
     * array. Typically we start *below* log2MinSize.
     * The sizes are numRows, not the total signal size.
     */
    unsigned log2MinSize,      
    unsigned log2MaxSize, 
    const size_t *stripeSize2DReal,
    unsigned *min2DRealSize,    // RETURNED
    Verbosity verbosity)
{
    FFTParams params;
    
    params.dims = 2;
    params.isReal = true;
    
    unsigned nRow = REAL_2D_START;
    unsigned nCol = REAL_2D_START;
    
    do { 
        unsigned nTotal = nRow + nCol;
        size_t realStripeSize = 0;
        if(stripeSize2DReal != NULL) {
            if(nRow < log2MinSize) {
                realStripeSize = stripeSize2DReal[0];
            }
            else {
                realStripeSize = stripeSize2DReal[nRow - log2MinSize];
            }
        }
        bool mfftFaster = false;
        int rtn;
        
        params.n[0] = nRow;
        params.n[1] = nCol;
        params.loops = calcLoops(false, nTotal, true);
        rtn = compareFFTs(&params,
            0,              // complexStripeSize,
            realStripeSize,
            0, false,       // rectOff, rectOffSpecd
            verbosity,
            &mfftFaster);
        if(rtn) {
            return rtn;
        }
        if(mfftFaster) {
            if(verbosity != V_Normal) {
                printf("++++++ 2D real minSize=%u\n", nTotal);
            }  
            *min2DRealSize = nTotal;
            return 0;
        }
        
        if(nRow == nCol) {
            nRow++;
        }
        else {
            nCol++;
        }
    } while(nRow <=log2MaxSize);
    
    /*
     * I'd be surprised if we ever reach here, other than if 
     * we're running with a very restricted max size
     */
    if(verbosity != V_Normal) {
        printf("++++++ 2D real loop exhausted; minSize=%u\n", 2 * log2MaxSize);
    }  
    *min2DRealSize = 2* log2MaxSize;
    return 0;
}


#pragma mark --- main() ---

#define MAX_ARRAY_SIZE      32

int main(int argc, char **argv)
{
	size_t minSize1DReal = MIN_SIZE_DEF;
    size_t maxSize1DReal = 0;
    Verbosity verbosity = V_Normal;
    bool quick = false;
    HWInfo hwInfo;
    bool showHwInfo = false;
    
    /* 
     * enables for debugging one specific table 
     * Note 1-D real only is impossible - it requires 1-D complex
     */
    bool oneDimEnable = true;
    bool twoDimEnable = true;
    bool realEnable = true;
    bool complexEnable = true;
    bool minMfftEnable = true;
    
    if(getHwInfo(&hwInfo)) {
        printf("***Error obtianing CPU/cache info\n");
        exit(1);
    }

	int arg;
	while ((arg = getopt(argc, argv, "n:x:vV12criqiMh")) != -1) {
		switch (arg) {
			case 'n':
				minSize1DReal = fftParseStringRep(optarg);
				break;
			case 'x':
				maxSize1DReal = fftParseStringRep(optarg);
				break;
			case 'v':
				verbosity = V_Verbose;
				break;
			case 'V':
				verbosity = V_Noisy;
				break;
            case 'q':
                quick = true;
                break;
            case '1':
                twoDimEnable = false;
                break;
            case '2':
                oneDimEnable = false;
                break;
            case 'c':
                realEnable = false;
                break;
            case 'r':
                complexEnable = false;
                break;
            case 'M':
                complexEnable = false;
                realEnable = false;
                oneDimEnable = false;
                twoDimEnable = false;
                break;
            case 'i':
                showHwInfo = true;
                break;
			default:
			case 'h':
				usage(argv);
		}
	}
	if(optind != argc) {
		usage(argv);
	}
    if(showHwInfo) {
        showInfo(&hwInfo);
        return 0;
    }
    if(maxSize1DReal == 0) {
        /* User hasn't specified max size; infer from RAM size. */
        if(maxFftSize(hwInfo.RAMSize, &maxSize1DReal)) {
            printf("***Error determining max FFT size\n");
            exit(1);
        }
    }
    
    /* 
     * Note: even though we don'e calculate stripe size for 1-D real (because
     * that op doesn't do anything with stripes), we specify sizes here in terms
     * if 1-D real because that's the easiest boundary to visualize.
     */
    unsigned log2Min1DReal;
    unsigned log2Max1DReal;
    if(!fftIsPowerOfTwo(minSize1DReal, &log2Min1DReal)) {
        printf("***min and max sizes must be powers of 2\n");
        exit(1);
    }
    if(!fftIsPowerOfTwo(maxSize1DReal, &log2Max1DReal)) {
        printf("***min and max sizes must be powers of 2\n");
        exit(1);
    }
    if(minSize1DReal > maxSize1DReal) {
        printf("***minSize must be less than maxSize\n");
        exit(1);
    }
    
    /*
     * -- minSize1DReal and maxSize1DReal are the literal bounds of the 1-D real FFT table.
     * -- log2Min1DReal and log2Max1DReal are log2 of the 1-D min and max values. The actual
     *    table is indexed via n-log2Min1DReal (i.e. log2Min1DReal is the first element in
     *    the table).
     * -- 1-D complex table is needed to calculate 1-D real table since 1-D real FFT first
     *    performs a 1-D complex FFT. When calculating entries for 1-D real of signal size
     *    n we set the app-specified offset and complex stripe size to element n-1 in the 
     *    1-D complex tables. 
     * -- The bounds for the 1-D complex table are log2*1DComplex = log2*1DReal-1
     * -- Lookups for 1-D signals - real and complex - are via log2 of the total signal size.
     * -- Lookups for 2-D signals are via log2(numRows).  
     * -- Bounds for 2-D real: log2*2DReal = log2*1DReal/2
     * -- Bounds for 2-D complex: log2*2DComplex = log2*2DReal-1
     * -- rectOffset table is ONLY used for 1-D complex, so it has the same bounds and 
     *    indexing as the 1-D complex table. 
     */
    unsigned log2Min2DReal    = log2Min1DReal >> 1;
    unsigned log2Max2DReal    = log2Max1DReal >> 1;
    unsigned log2Min1DComplex = log2Min1DReal - 1;
    unsigned log2Max1DComplex = log2Max1DReal - 1;
    unsigned log2Min2DComplex = log2Min2DReal - 1;
    unsigned log2Max2DComplex = log2Max2DReal - 1;
    
    /* The tables whose generation are the reason we exist */
    int rectOffset[MAX_ARRAY_SIZE];
    size_t stripeSize1DComplex[MAX_ARRAY_SIZE];
    size_t stripeSize2DComplex[MAX_ARRAY_SIZE];
    size_t stripeSize1DReal[MAX_ARRAY_SIZE];
    size_t stripeSize2DReal[MAX_ARRAY_SIZE];
   
    memset(rectOffset, 0, sizeof(rectOffset));
    memset(stripeSize1DComplex, 0, sizeof(stripeSize1DComplex));
    memset(stripeSize2DComplex, 0, sizeof(stripeSize2DComplex));
    memset(stripeSize1DReal, 0, sizeof(stripeSize1DReal));
    memset(stripeSize2DReal, 0, sizeof(stripeSize2DReal));
    
    printf("\nGenerating configuration tables for the following maximum FFT sizes:\n");
    if(oneDimEnable) {
        if(complexEnable) {
            printf("1-D complex   : 2^%u elements\n", log2Max1DComplex);
        }
        if(realEnable) {
            printf("1-D real      : 2^%u elements\n", log2Max1DReal);
        }
    }
    if(twoDimEnable) {
        if(complexEnable) {
            printf("2-D complex   : 2^%u rows\n", log2Max2DComplex);
        }
        if(realEnable) {
            printf("2-D real      : 2^%u rows\n", log2Max2DReal);
        }
    }
    printf("Configuration : %s precision %s complex\n\n",
        FFT_DOUBLE_PREC   ? "Double" : "Single", 
        FFT_SPLIT_COMPLEX ? "Split" : "Interleaved"); 
    
    if(oneDimEnable) {
        if(complexEnable) {
            /* This generates rectangle offsets and stripeSize1DComplex[] */
            if(gen1DComplex(log2Min1DComplex, log2Max1DComplex, 
                        rectOffset, stripeSize1DComplex, 
                        quick, verbosity)) {
                exit(1);
            }
        }
    }
    
    if(twoDimEnable) {
        if(complexEnable) {
            if(gen2DComplex(log2Min2DComplex, log2Max2DComplex, stripeSize2DComplex, quick, verbosity)) {
                exit(1);
            }
        }
        if(realEnable) {
            if(gen2DReal(log2Min2DReal, log2Max2DReal, stripeSize2DReal, quick, verbosity)) {
                exit(1);
            }
        }
    }
    
    /* 
     * Min MatrixFFT sizes
     */
    unsigned minMfftComplex1D = 0;
    unsigned minMfftComplex2D = 0;
    unsigned minMfftReal1D = 0;
    unsigned minMfftReal2D = 0; 
    
    if(minMfftEnable) {
        if(gen1DComplexMinSize(log2Min1DComplex, log2Max1DComplex,
                (complexEnable && oneDimEnable) ? rectOffset : NULL,
                (complexEnable && oneDimEnable) ? stripeSize1DComplex : NULL,
                &minMfftComplex1D, verbosity)) {
            exit(1);
        }
        if(gen2DComplexMinSize(log2Min2DComplex, log2Max2DComplex,
                (complexEnable && twoDimEnable) ? stripeSize2DComplex : NULL,
                &minMfftComplex2D, verbosity)) {
            exit(1);
        }
        /* this one uses stripe sizes and rectOffset from the 1-D complex run */
        if(gen1DRealMinSize(log2Max1DReal, log2Min1DComplex, log2Max1DComplex,
                (complexEnable && oneDimEnable) ? rectOffset : NULL,
                (complexEnable && oneDimEnable) ? stripeSize1DComplex : NULL,
                &minMfftReal1D, verbosity)) {
            exit(1);
        }
        if(gen2DRealMinSize(log2Min2DReal, log2Max2DReal,
                (realEnable && twoDimEnable) ? stripeSize2DReal : NULL,
                &minMfftReal2D, verbosity)) {
            exit(1);
        }
    }
    
    /* 
     * Output compilable C code of all arrays as enabled
     */
    printf("\n\n");
    printf("   #if   ");
    configAsCppString();
    putchar('\n');
    
    printf("\n%s/* %s precision %s complex */\n\n",
        INDENT, 
        FFT_DOUBLE_PREC   ? "Double" : "Single", 
        FFT_SPLIT_COMPLEX ? "Split" : "Interleaved");
    
    if(realEnable) {
        if(twoDimEnable) {
            printf("%s#define CONFIG_TABLE_2D_REAL_START  %2u   ", INDENT, log2Min2DReal);
            printf("/* log2(n) of first entry of 2-D real table */\n");
            printf("%s#define CONFIG_TABLE_2D_REAL_END    %2u   ", INDENT, log2Max2DReal);
            printf("/* log2(n) of last  entry of 2-D real table */\n");
        }
    }
    if(complexEnable) {
        if(oneDimEnable) {
            printf("%s#define CONFIG_TABLE_1D_CMPL_START  %2u   ", INDENT, log2Min1DComplex);
            printf("/* log2(n) of first entry of 1-D complex table */\n");
            printf("%s#define CONFIG_TABLE_1D_CMPL_END    %2u   ", INDENT, log2Max1DComplex);
            printf("/* log2(n) of last  entry of 1-D complex table */\n");
        }
        if(twoDimEnable) {
            printf("%s#define CONFIG_TABLE_2D_CMPL_START  %2u   ", INDENT, log2Min2DComplex);
            printf("/* log2(n) of first entry of 2-D complex table */\n");
            printf("%s#define CONFIG_TABLE_2D_CMPL_END    %2u   ", INDENT, log2Max2DComplex);
            printf("/* log2(n) of last  entry of 2-D complex table */\n");
        }
    }
    if(realEnable || complexEnable) {
        putchar('\n');
    }
    
    unsigned numSizes2DReal    = log2Max2DReal    - log2Min2DReal    + 1;
    unsigned numSizes1DComplex = log2Max1DComplex - log2Min1DComplex + 1;
    unsigned numSizes2DComplex = log2Max2DComplex - log2Min2DComplex + 1;
    
    if(oneDimEnable) {
        printf("%sstatic int mfftRectOffsetTable[] = \n", INDENT);
        printf("%s       { ", INDENT);
        for(unsigned dex=0; dex<numSizes1DComplex; dex++) {
            printf("%d", rectOffset[dex]);
            if(dex < (numSizes1DComplex - 1)) {
                printf(", ");
            }
        }
        printf("};\n\n");
    }
    
    if(complexEnable) {
        if(oneDimEnable) {
            printStripeSizeTable("fftStripeSize1DComplex", stripeSize1DComplex, numSizes1DComplex);
        }
        if(twoDimEnable) {
            printStripeSizeTable("fftStripeSize2DComplex", stripeSize2DComplex, numSizes2DComplex);
        }
    }
    if(realEnable) {
        if(twoDimEnable) {
            printStripeSizeTable("fftStripeSize2DReal",    stripeSize2DReal, numSizes2DReal);
        }
    }
    
    if(minMfftEnable) {
        if(realEnable || complexEnable) {
            printf("\n");
        }
        printf("%s#define MFFT_1D_COMPLEX_MIN_SIZE    %2u\n", INDENT, minMfftComplex1D);
        printf("%s#define MFFT_2D_COMPLEX_MIN_SIZE    %2u\n", INDENT, minMfftComplex2D);
        printf("%s#define MFFT_1D_REAL_MIN_SIZE       %2u\n", INDENT, minMfftReal1D);
        printf("%s#define MFFT_2D_REAL_MIN_SIZE       %2u\n\n", INDENT, minMfftReal2D);
    }

    printf("   #endif   /* ");
    configAsCppString();
    printf(" */\n\n");

	return 0;
}
