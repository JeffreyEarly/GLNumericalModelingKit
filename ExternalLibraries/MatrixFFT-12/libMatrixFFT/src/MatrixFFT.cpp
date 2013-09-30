/*	File: MatrixFFT.cpp  
	
	Description:
		Public MatrixFFT functions.
	
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
 * MatrixFFT.cpp - Public MatrixFFT functions. 
 *			  
 * Created 12/16/2008. 
 * Copyright 2008 by Apple, Inc. 
 */

#include <libMatrixFFT/MatrixFFT.h>
#include <libMatrixFFT/fftUtils.h>
#include <libMatrixFFT/complexBufUtils.h>
#include "MatrixFFTPlan.h"
#include "fftSinCos.h"
#include "fftPriv.h"
#include "ThreadPool.h"
#include "fftTranspose.h"
#include "fftThreadOps.h"
#include "FFTEngineVDSP.h"
#include "fftDebug.h"

#define LOG_STRIPE_BUF_SIZE		0

/*
 * Determine L2 cache available for each FFTEngine. 
 * If the L2 cache per core is a power of 2, each engine gets half of that. 
 * Else each engine gets the largest power of 2 smaller than L2 cache size.
 */
static uint32_t calcL2CachePerEngine()
{
	uint64_t cachePerCore = l2CachePerCore();
	
	if(cachePerCore == 0) {
		/* sysctl error, punt */
		return 0;
	}
	if(cachePerCore > ((uint32_t) -1)) {
		/* Really? over 4 GB of L2 cache per core? */
		dbprintf("l2CachePerEngine: cachePerCore > uint32_t\n");
		cachePerCore = 0x8000000UL;
	}
	unsigned pwrTwo;
	bool isPwrTwo = fftIsPowerOfTwo((size_t)cachePerCore, &pwrTwo);
	if(isPwrTwo) {
		return cachePerCore >> 1;
	}
	else {
		return (uint32_t)1 << pwrTwo;
	}
}

/*
 * As of 23 Nov. 2009, this flag hasn't been seen to make much of a
 * difference either way, but it's still conceptually the right thing to 
 * do since our threads don't share anything. 
 */
#define THREAD_POOL_OPTIONS TPO_SeparateAffinity

/* 
 * Common MatrixFFTPlan creation. 
 */
MFFTReturn mfftPlanCreateCommon(
	unsigned			dims,			// 1 or 2
	FFTSineTableType	sinTableType,
	bool				real,
	uint32_t			optFlags,
	/* rectangle dimensions, even for 1-D */
	unsigned			nRow,			// two dimension : N = 2 ^ (nRow + nCol)
	unsigned			nCol,
	unsigned			numThreads,
	size_t				sinPeriod,		// only for STT_Optimized
    uint32_t            options,        // MFFT_PLAN_OPT_*
	MatrixFFTPlan		*mfftPlanRtn)	// RETURNED
{
	MatrixFFTPlan mfftPlan = (MatrixFFTPlan)malloc(sizeof(*mfftPlan));
	if(mfftPlan == NULL) {
		return MR_Memory;
	}
	/* subsequent errors to errOut: */
	MFFTReturn mrtn = MR_Success;
	memset(mfftPlan, 0, sizeof(*mfftPlan));
	
	mfftPlan->nRow = nRow;
	mfftPlan->nCol = nCol;
	mfftPlan->numRows = (size_t)1 << nRow;
	mfftPlan->numCols = (size_t)1 << nCol;
	mfftPlan->N = (size_t)1 << (nRow + nCol);
	mfftPlan->realFft = real;
	if(real) {
		mfftPlan->numColsComplex = mfftPlan->numCols >> 1;
		mfftPlan->nColComplex = nCol - 1;
	}
	else {
		mfftPlan->numColsComplex = mfftPlan->numCols;
		mfftPlan->nColComplex = nCol;
	}
	mfftPlan->sinTableType = sinTableType;

	if(!(options & MFFT_PLAN_OPT_NO_VDSP)) {
		unsigned maxLog = max(nRow, nCol);
		mfftPlan->vdspSetup = FFTCreateSetup(maxLog);
		if(mfftPlan->vdspSetup == NULL) {
			printf("***Error creating vDSP setup for log2(N) = %u\n", maxLog);
			return MR_Memory;
		}
	}
	
	if(numThreads == 0) {
		numThreads = numCpuCores();
		if(numThreads == 0) {
			/* Something went wrong with numCpuCores()... */
			numThreads = 1;
		}
	}
	
	mfftPlan->l2CachePerEngine = calcL2CachePerEngine();
	
	if(!(options & MFFT_PLAN_OPT_NO_ENGINES)) {
        /* 
         * Additional FFTEngines - e.g. OpenCL-based - can be added here. 
         * For now we just use numThreads vDSP-based engines.
         */
        unsigned numEngines = mfftPlan->numHostEngines = numThreads;
        mfftPlan->fftEngines = (FFTEngine **)malloc(numEngines * sizeof(FFTEngine *));
        if(mfftPlan->fftEngines == NULL) {
            mrtn = MR_Memory;
            goto errOut;
        }
        try {
            float opRatio = 1.0 / (float)(numEngines);
            for(unsigned engineDex=0; engineDex<numEngines; engineDex++) {
                FFTEngineVDSP *hostEngine = new FFTEngineVDSP(mfftPlan, dims, 
                    mfftPlan->numRows, mfftPlan->numCols, real);
                mfftPlan->fftEngines[engineDex] = hostEngine;
                sprintf(hostEngine->mEngineName, "vDSP%u", engineDex);
                for(unsigned opDex=0; opDex<FEO_NumOps; opDex++) {
                    hostEngine->mOpRatios[opDex] = opRatio;
                }
            }
        }
        catch(FFTError cerr) {
            printf("***mfftPlanCreateCommon: FFTEngine creation error\n");
            mrtn = cerr.err();
            goto errOut;
        }
        catch(...) {
            printf("***mfftPlanCreateCommon: Unknown exception\n");
            mrtn = MR_Internal;
            goto errOut;
        }

        #if		LOG_STRIPE_BUF_SIZE
        {
            unsigned long stripeSize = mfftPlan->fftEngines[0]->stripeBufSize() * 
                2 * sizeof(FFTFloat);
            printf("==== stripeBufSize %lu bytes\n", stripeSize);
        }
        #endif	/* LOG_STRIPE_BUF_SIZE */
	}
    
	/* One pthread per engine - but avoid creating unnecessary singleton pthread */
	if(!(options & MFFT_PLAN_OPT_NO_THREADS)) {
		if(numThreads > 1) {
			mrtn = tpThreadInit(&mfftPlan->threadPool, numThreads, 
				0,		/* no extra - nothing done in main() when using threads */
				sizeof(TP_TaskUnion_U),
                THREAD_POOL_OPTIONS);
			if(mrtn) {
				goto errOut;
			}
		}
	}
	/* always need this to get threadPoolP->numThreads */
	mfftPlan->threadPoolP = &mfftPlan->threadPool;
	
	switch(sinTableType) {
		case STT_None:
		case STT_External:
			break;
		default:
			mfftPlan->sinCosNMask  = mfftPlan->N - 1;
			mfftPlan->piOver2      = mfftPlan->N >> 2;
			mfftPlan->pi           = mfftPlan->N >> 1;
			mfftPlan->threePiOver2 = mfftPlan->piOver2 * 3;
			break;
	}
	
	switch(sinTableType) {
		case STT_None:
		case STT_External:
			/* we're done */
			break;
		case STT_Standard:
			if(fftSetupSineNorm(mfftPlan)) {
				mrtn = MR_Memory;
			}
			break;
		case STT_Optimized:
			mfftPlan->sinPeriod = sinPeriod;
			if(fftSetupSineOpt(mfftPlan)) {
				mrtn = MR_Memory;
			}
			break;
	}
	
	mfftPlan->addNormFact = 0.0;
	
errOut:
	if(mrtn) {
		mfftFreePlan(mfftPlan);
	}
	else {
		*mfftPlanRtn = mfftPlan;
	}
	return mrtn;
}

/* Generic NULL transpose */
MFFTReturn fftNullTranspose(
	MatrixFFTPlan		mfftPlan,
	bool				forward,		
	bool				input, 
	const FFTComplex	*inBuf,
	FFTComplex			*outBuf)
{
	if(inBuf == outBuf) {
		return MR_Success;
	}
	size_t toMove = mfftPlan->numRows * mfftPlan->numCols * sizeof(FFTComplex);
	memmove(outBuf, inBuf, toMove);
	return MR_Success;
}

#pragma mark --- Public API ---

/* See public header for documentation on these functions */

MFFTReturn mfftCreatePlan(
	unsigned			dims,		
	unsigned			*n,			
	bool				real,
	uint32_t			optFlags,
	unsigned			numThreads,
	MatrixFFTPlan		*mfftPlanRtn)		/* RETURNED */
{
	if(mfftPlanRtn == NULL) {
		return MR_IllegalArg;
	}
	
	/*
	 * Parse the args into an alg-specific configuration. 
	 */
	MFFTReturn mrtn = MR_Success;
	MatrixFFTPlan mfftPlan = NULL;
	
	if(real) {
		switch(dims) {
			case 1:
				/* 1-D Real */
				mrtn = mfftCreatePlan1DReal(n, optFlags, numThreads, &mfftPlan);
				break;
			case 2:
				/* 2-D Real */
				mrtn = mfftCreatePlan2DReal(n, optFlags, numThreads, &mfftPlan);
				break;
			default:
				printf("***mfftCreatePlan: unsupported dimension\n");
				return MR_Unsupported;
		}
	}
	else {
		/* Complex */
		switch(dims) {
			case 1:
				/* 1-D complex */
				mrtn = mfftCreatePlan1DComplex(n, optFlags, numThreads, false, false, false, &mfftPlan);
				break;
			case 2:
				/* 2-D complex */
				mrtn = mfftCreatePlan2DComplex(n, optFlags, numThreads, &mfftPlan);
				break;
			default:
				printf("***mfftCreatePlan: unsupported dimension\n");
				return MR_Unsupported;
		}
	}
	if(mrtn) {
		goto errOut;
	}

	RFASSERT(mfftPlan->execFcn != NULL);
	RFASSERT(mfftPlan->transFcn != NULL);
	
errOut:
	if(mrtn) {
		if(mfftPlan != NULL) {
			mfftFreePlan(mfftPlan);
		}
	}
	else {
		*mfftPlanRtn = mfftPlan;
	}
	return mrtn;
}

/*
 * Free a MatrixFFTPlan.
 */
void mfftFreePlan(MatrixFFTPlan mfftPlan)
{
	if(mfftPlan == NULL) {
		return;
	}
	tpThreadShutdown(&mfftPlan->threadPool);
	fftFreeComplexArrayAlign(mfftPlan->auxBuf, mfftPlan->auxBufFree);
	COND_FREE(mfftPlan->sine);
	COND_FREE(mfftPlan->cosine);
	
	for(unsigned dex=0; dex<mfftPlan->numHostEngines; dex++) {
		delete mfftPlan->fftEngines[dex];
	}
	COND_FREE(mfftPlan->fftEngines);
	
	if(mfftPlan->vdspSetup) {
		FFTFreeSetup(mfftPlan->vdspSetup);
		mfftPlan->vdspSetup = NULL;
	}
	#if	!FFT_SPLIT_COMPLEX
	fftFreeDSPComplex(&mfftPlan->vdspBufFree);
	#endif
	
	if(mfftPlan->subPlan != NULL) {
		mfftFreePlan(mfftPlan->subPlan);
		mfftPlan->subPlan = NULL;
	}
    
	free(mfftPlan);
}

/* 
 * Execute FFT for the specified MatrixFFTPlan.
 */
MFFTReturn mfftExecute(
	MatrixFFTPlan		mfftPlan,
	uint32_t			optFlags,
	bool				forward,		
	FFTComplex			*inBuf,
	FFTComplex			*outBuf)
{
	if(mfftPlan == NULL) {
		return MR_IllegalArg;
	}
	if(mfftPlan->execFcn == NULL) {
		return MR_Internal;
	}
	return mfftPlan->execFcn(mfftPlan, optFlags, forward, inBuf, outBuf);
}

/* 
 * Obtain optimal native format for the input to, or the output from,
 * the forward FFTs for the specified MatrixFFTPlan. 
 */
MFFTReturn mfftNativeFormat(
	MatrixFFTPlan		mfftPlan,
	bool				input,
	MFFTFormat			*format)		/* RETURNED */
{
	if((mfftPlan == NULL) || (format == NULL)) {
		return MR_IllegalArg;
	}
	if(input) {
		*format = mfftPlan->inputFormat;
	}
	else {
		*format =  mfftPlan->outputFormat;
	}
	return MR_Success;
}

/*
 * Return the number of threads in use by a MatrixFFTPlan.
 */
extern MFFTReturn mfftNumThreads(
	MatrixFFTPlan		mfftPlan,
	unsigned			*numThreads)	/* RETURNED */
{
	if(numThreads == NULL) {
		return MR_IllegalArg;
	}
	if(mfftPlan == NULL) {
		/* return the default */
		*numThreads = numCpuCores();
	}
	else {
		RFASSERT(mfftPlan->numHostEngines == mfftPlan->threadPoolP->numThreads);
		*numThreads = mfftPlan->numHostEngines;
	}
	return MR_Success;
}

/*
 * Obtain the number of rows and columns used in a MatrixFFTPlan.
 */
MFFTReturn mfftRectangle(
	MatrixFFTPlan		mfftPlan,
	size_t				*numRows,		/* RETURNED */
	size_t				*numCols)		/* RETURNED */
{
	if((mfftPlan == NULL) || (numRows == NULL) || (numCols == NULL)) {
		return MR_IllegalArg;
	}
	*numRows = mfftPlan->numRows;
	*numCols = mfftPlan->numCols;
	return MR_Success;
}

/*
 * Transpose a signal to/from normal row order.
 */
MFFTReturn mfftTranspose(
	MatrixFFTPlan		mfftPlan,
	bool				forward,		
	bool				input, 
	const FFTComplex	*inBuf,
	FFTComplex			*outBuf)
{
	if(mfftPlan == NULL) {
		return MR_IllegalArg;
	}
	if(mfftPlan->transFcn == NULL) {
		return MR_Internal;
	}
	return mfftPlan->transFcn(mfftPlan, forward, input, inBuf, outBuf);
}

/* 
 * Obtain C-string version of a MFFTReturn.
 */
const char *mfftErrStr(MFFTReturn mrtn)
{
	static char unk[200];
	
	switch(mrtn) {
		case MR_Success: return "MR_Success";
		case MR_Unsupported: return "MR_Unsupported";
		case MR_VDSP: return "MR_VDSP";
		case MR_Memory: return "MR_Memory";
		case MR_IllegalArg: return "MR_IllegalArg";
		case MR_Internal: return "MR_Internal";
		default:
			/* not thread safe */
			sprintf(unk, "Unknown error <%u>", (unsigned)mrtn);
			return unk;
	}
}

/* dump error info to stdout */
void mfftPrintErrInfo(
	const char *func, 
	MFFTReturn mrtn)
{
	printf("%s returned %s\n", func, mfftErrStr(mrtn));
}

