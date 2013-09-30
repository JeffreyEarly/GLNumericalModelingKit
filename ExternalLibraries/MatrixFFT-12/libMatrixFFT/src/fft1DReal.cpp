/*	File: fft1DReal.cpp 
	
	Description:
		1-D Real-signal FFT via matrix decomposition.
	
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
 * fft1DComplex.cpp - 1-D Real-signal FFT via matrix decomposition. 
 *			  
 * Created 1/5/2009. 
 * Copyright 2009 by Apple, Inc. 
 */

#include <libMatrixFFT/MatrixFFT.h>
#include <libMatrixFFT/fftUtils.h>
#include <libMatrixFFT/complexBufUtils.h>
#include "MatrixFFTPlan.h"
#include "fftPriv.h"
#include "fftTranspose.h"
#include "ThreadPool.h"
#include "fftSinCos.h"
#include "fftEngineDispatch.h"
#include "fftDebug.h"
#include "fftPlatformConf.h"

#pragma mark --- Forward 1-D Real FFT ---

/* 
 * Execute a forward 1-D real-signal FFT. 
 * In-place or out-of-place (if in-place, inBuf == outBuf). 
 */
static MFFTReturn fftExec1DRealFwd(
	MatrixFFTPlan		mfftPlan,
	uint32_t			optFlags,
	FFTComplex			*inBuf,
	FFTComplex			*outBuf)
{
	dumpRectReal("Exec1DRealFwd input", inBuf, mfftPlan->numRows, mfftPlan->numCols);

	/*
	 * 1-D complex FFT, output in column order (but with same numRows and numCols
	 * as original signal).
	 */
	MFFTReturn mrtn;
	if(mfftPlan->subPlan == NULL) {
		return MR_Internal;
	}
	mrtn = mfftExecute(mfftPlan->subPlan, 0, true, inBuf, outBuf);
	if(mrtn) {
		return mrtn;
	}
	
	dumpMatrixRect("Exec1DRealFwd 1-D complex output", outBuf, mfftPlan->numRows, 
		mfftPlan->numColsComplex);
		
	/* 
	 * Twist, in-place on outBuf.
	 */
	mrtn = fe1DRealTwist(mfftPlan, true, outBuf, outBuf);
	if(mrtn) {
		return mrtn;
	}
	dumpMatrixRect("Exec1DRealFwd 1-D twist output", outBuf, mfftPlan->numRows,
		mfftPlan->numColsComplex);
	
	/*
	 * Optional frequency-domain transpose.
	 */
	if(optFlags & MEF_TransposeOutput) {
		mrtn = mfftPlan->transFcn(mfftPlan, true, false, outBuf, outBuf);
		dumpMatrixRect("Exec1DRealFwd transposed output", outBuf, mfftPlan->numRows,
			mfftPlan->numColsComplex);
	}
	return mrtn;
}

#pragma mark --- Inverse 1-D Real FFT ---

/* 
 * Execute a reverse 1-D real-signal FFT. 
 * In-place or out-of-place (if in-place, inBuf == outBuf). 
 */
static MFFTReturn fftExec1DRealInv(
	MatrixFFTPlan		mfftPlan,
	uint32_t			optFlags,
	FFTComplex			*inBuf,
	FFTComplex			*outBuf)
{
	if(mfftPlan->subPlan == NULL) {
		return MR_Internal;
	}
	MFFTReturn mrtn;

	/*
	 * Optional transpose to column order.
	 */
	if(optFlags & MEF_TransposeInput) {
		dumpMatrixRect("Exec1DRealInv pre-transpose", inBuf, 
			mfftPlan->numRows, mfftPlan->numColsComplex);
		mrtn = mfftPlan->transFcn(mfftPlan, false, true, inBuf, inBuf);
		if(mrtn) {
			return mrtn;
		}
	}
	dumpMatrixRect("Exec1DRealInv input", inBuf, 
		mfftPlan->numRows, mfftPlan->numColsComplex);
		
	/*
	 * Twist, in-place on inBuf.
	 */
	mrtn = fe1DRealTwist(mfftPlan, false, inBuf, inBuf);
	if(mrtn) {
		return mrtn;
	}
	dumpMatrixRect("Exec1DRealInv twist output", outBuf, mfftPlan->numRows,
		mfftPlan->numColsComplex);
		
	/*
	 * 1-D Complex FFT, inBuf-->outBuf. Pass along MEF_NormOutput flag.
	 */
	mrtn = mfftExecute(mfftPlan->subPlan, 
		optFlags & MEF_NormOutput, 
		false, inBuf, outBuf);
	if(mrtn) {
		return mrtn;
	}
	dumpRectReal("Exec1DRealInv output", outBuf, mfftPlan->numRows,	mfftPlan->numCols);
	return mrtn;
}

#pragma mark --- MatrixFFTPlan callouts ---

static MFFTReturn fftExec1DReal(
	MatrixFFTPlan		mfftPlan,
	uint32_t			optFlags,
	bool				forward,		
	FFTComplex			*inBuf,
	FFTComplex			*outBuf)
{
	if(forward) {
		return fftExec1DRealFwd(mfftPlan, optFlags, inBuf, outBuf);
	}
	else {
		return fftExec1DRealInv(mfftPlan, optFlags, inBuf, outBuf);
	}
}

#pragma mark --- 1-D Real FFT via vDSP ---

#if     FFT_SPLIT_COMPLEX

/*
 * Split complex form, this is pretty easy...
 */
static MFFTReturn fftExec1DRealVdsp(
	MatrixFFTPlan		mfftPlan,
	uint32_t			optFlags,
	bool				forward,		
	FFTComplex			*inBuf,
	FFTComplex			*outBuf)
{
    vDSPComplex vbufIn;
    vDSPComplex vbufOut;
    
    fftComplexToVDSP(inBuf, &vbufIn);
    fftComplexToVDSP(outBuf, &vbufOut);

    bool inPlace = (inBuf == outBuf) || 
                   ( (inBuf->real == outBuf->real) && 
                     (inBuf->imag == outBuf->imag));
                    
	FFTDirection fdir = forward ? FFT_FORWARD : FFT_INVERSE;
    if(inPlace) {
        FFTReal1d(mfftPlan->vdspSetup, &vbufIn, mfftPlan->nCol, fdir);
    }
    else {
        FFTReal1dOP(mfftPlan->vdspSetup, &vbufIn, &vbufOut, mfftPlan->nCol, fdir);
    }
    if(!forward && (optFlags & MEF_NormOutput)) {
        /* scale by 1/2N */
        FFTFloat normFact = 0.5 / (FFTFloat)(mfftPlan->numCols);
		fftScaleComplex(outBuf, normFact, mfftPlan->numCols >> 1);
    }
    return MR_Success;
}

#else   /* FFT_SPLIT_COMPLEX */

/* 
 * Interleaved, we have to convert to/from vDSP's native split format.
 */
static MFFTReturn fftExec1DRealVdsp(
	MatrixFFTPlan		mfftPlan,
	uint32_t			optFlags,
	bool				forward,		
	FFTComplex			*inBuf,
	FFTComplex			*outBuf)
{
    RFASSERT(mfftPlan->vdspBuf.realp != NULL);
	FFTDirection fdir = forward ? FFT_FORWARD : FFT_INVERSE;
    size_t numComplex = mfftPlan->numCols >> 1;
    fftIntToVDSP(inBuf, &mfftPlan->vdspBuf, numComplex);
    FFTReal1d(mfftPlan->vdspSetup, &mfftPlan->vdspBuf, mfftPlan->nCol, fdir);
    fftVDSPToInt(&mfftPlan->vdspBuf, outBuf, numComplex);
    if(!forward && (optFlags & MEF_NormOutput)) {
        /* scale by 1/2N */
        FFTFloat normFact = 0.5 / (FFTFloat)(mfftPlan->numCols);
		fftScaleComplex(outBuf, normFact, numComplex);
    }
    return MR_Success;
}

#endif  /* FFT_SPLIT_COMPLEX */

#pragma mark --- config-specific MatrixFFTPlan init ---

/*
 * Set up a MatrixFFTPlan for raw vDSP-based FFTs.
 */
static MFFTReturn mfftCreatePlan1DRealVdsp(
	unsigned			*n,			
	uint32_t			optFlags,
	MatrixFFTPlan		*mfftPlanRtn)
{
    MatrixFFTPlan mfftPlan = NULL;
    MFFTReturn mrtn = mfftPlanCreateCommon(1, 
            STT_None,   // no sine table
            true,       // real
            optFlags, 
            0, n[0],    // nRows, nCols
            0,          // numThreads, actually will be ignored
            0,          // sinPeriod
            MFFT_PLAN_OPT_NO_THREADS | MFFT_PLAN_OPT_NO_ENGINES,
            &mfftPlan);
    if(mrtn) {
        return mrtn;
    }

    mfftPlan->inputFormat = MF_RowOrder;
    mfftPlan->transFcn    = fftNullTranspose;
    mfftPlan->execFcn	  = fftExec1DRealVdsp;

    #if		!FFT_SPLIT_COMPLEX
    /* We need a vDSP buffer for conversion from interleaved */
    if(fftAllocDSPComplexAlign(&mfftPlan->vdspBuf, 
            mfftPlan->numCols >> 1, VECTOR_SIZE, 
            &mfftPlan->vdspBufFree)) {
        mfftFreePlan(mfftPlan);
        return MR_Memory;
    }
    #endif

    *mfftPlanRtn = mfftPlan;
    return MR_Success;
}

MFFTReturn mfftCreatePlan1DReal(
	unsigned			*n,			
	uint32_t			optFlags,
	unsigned			numThreads,
	MatrixFFTPlan		*mfftPlanRtn)
{
    if(mfftShouldUseVdsp(n[0], true, 1)) {
        /* Set up a plan for vDSP only */
        return mfftCreatePlan1DRealVdsp(n, optFlags, mfftPlanRtn);
    }

	/*
	 * First, set up a subplan for performing a 1-D complex FFT.
	 * We're going to set up a fully populated sine table, larger than
	 * the one the 1-D complex needs (by a factor of 2), so disable
	 * sin table creation here.
	 */
	MatrixFFTPlan subPlan = NULL;
	MFFTReturn mrtn;
	unsigned nComplex = (*n)-1;
	
	mrtn = mfftCreatePlan1DComplex(&nComplex, 
        optFlags & MCF_HintTranspose,       // pass this along 
        numThreads,
		true,		// externalSineTable
		true,		// disableThreads
        true,       // expect column order output
		&subPlan);
	if(mrtn) {
		return mrtn;
	}
	/* subsequent errors to errOut: */
	
	/*
	 * The subplan's numRows and numCols are in complex elements; our coordinates
	 * are in reals, using the same layout as the subplan. 
	 */
	unsigned nRow = subPlan->nRow;
	unsigned nCol = subPlan->nCol + 1;
	
	MatrixFFTPlan mfftPlan = NULL;
	
	/*
	 * We do not need a vDSP plan; all vDSP ops are performed by the 
	 * subplan. All we need is a sin table for the 1-D real twist. 
	 */
	mrtn = mfftPlanCreateCommon(1, STT_Standard,
		true, optFlags, nRow, nCol, numThreads,
		0,				// sinPeriod,
        MFFT_PLAN_OPT_NO_VDSP,
		&mfftPlan);
	if(mrtn) {
		goto errOut;
	}

	mfftPlan->inputFormat   = MF_RowOrder;
	mfftPlan->outputFormat  = MF_ColumnOrder;
	if(subPlan->nRow == subPlan->nCol) {
		mfftPlan->transFcn  = fftTranspose2DComplexSquare;
	}
	else {
		mfftPlan->transFcn  = fftTranspose1DComplexRect;
	}
	mfftPlan->execFcn	    = fftExec1DReal;
	
	/* Connect our plan to our 1-D complex subplan */
	mfftPlan->subPlan = subPlan;
	subPlan->externSineTable = mfftPlan;
	subPlan->sineShift = 1;
	subPlan->threadPoolP = &mfftPlan->threadPool;
	
	/*
	 * Normalizing done by the subplan. 
	 * 1D complex normFactor = 1 / Ncomplex
	 * 1D Real normFactor    = 1 / (2 * Nreal)
	 * 1 / (2 * Nreal) = (1 / Ncomplex) * 0.25
	 */
	subPlan->addNormFact = 0.25;
	
	*mfftPlanRtn = mfftPlan;
	return MR_Success;

errOut:
	if(subPlan != NULL) {
		mfftFreePlan(subPlan);
	}
	if(mfftPlan != NULL) {
		mfftFreePlan(mfftPlan);
	}
	return mrtn;
}
