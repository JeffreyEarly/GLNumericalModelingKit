/*	File: fft2DComplex.cpp 
	
	Description:
		2-D Complex FFT functions.
	
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
 * fft2DComplex.cpp - 2-D Complex FFT functions.
 *			  
 * Created 9/24/2008. 
 * Copyright 2008 by Apple, Inc. 
 */

#include <libMatrixFFT/MatrixFFT.h>
#include <libMatrixFFT/fftUtils.h>
#include <libMatrixFFT/complexBufUtils.h>
#include "MatrixFFTPlan.h"
#include "fftPriv.h"
#include "fftTranspose.h"
#include "ThreadPool.h"
#include "fftEngineDispatch.h"
#include "fftDebug.h"
#include "fftPlatformConf.h"

/* 
 * For square 2-D FFTs at least 2^FFT_2D_ROW_TRANS_SIZE elements on a side,
 * perform full square transpose and do the FFT-columns as rows, leaving
 * the result in column order.
 */
#if	FFT_DOUBLE_PREC
	#define FFT_2D_ROW_TRANS_SIZE	10
#else
	#if	FFT_SPLIT_COMPLEX
		#define FFT_2D_ROW_TRANS_SIZE	10
	#else
		#define FFT_2D_ROW_TRANS_SIZE	14
	#endif
#endif

/* Turn this off to disable square transpose option */
#define FFT_2D_SQUARE_ENABLE        1

#pragma mark --- 2-D complex FFT with column FFTs in place ---

/* 
 * Execute a 2-D complex FFT. 
 * In-place or out-of-place (if in-place, inBuf == outBuf). 
 */
static MFFTReturn fftExecute2DComplex(
	MatrixFFTPlan		mfftPlan,
	uint32_t			optFlags,
	bool				forward,		
	FFTComplex			*inBuf,
	FFTComplex			*outBuf)
{
	MFFTReturn mrtn = MR_Success;

	dumpMatrixRect("fftExecute2DComplex input", inBuf, mfftPlan->numRows, mfftPlan->numCols);
	if(forward) {
		mrtn = feFftAllRows(mfftPlan, forward, 0, inBuf, outBuf, 0.0);
		if(mrtn) {
			return mrtn;
		}
		dumpMatrixRect("fftExecute2DComplex(fwd) ceFftAllRows output", outBuf, 
			mfftPlan->numRows, mfftPlan->numCols);
		mrtn = feFftColumns(mfftPlan, forward, 0, outBuf, outBuf, 0.0);
	}
	else {
		FFTFloat normFact = 0.0;
		
		if(optFlags & MEF_NormOutput) {
			normFact = 1.0 / ((FFTFloat)(mfftPlan->numRows * mfftPlan->numCols));
		}
		mrtn = feFftColumns(mfftPlan, forward, 0, inBuf, outBuf, 0.0);
		dumpMatrixRect("fftExecute2DComplex(inv) ceFftAllCols output", outBuf, 
			mfftPlan->numRows, mfftPlan->numCols);
		if(mrtn) {
			return mrtn;
		}
		mrtn = feFftAllRows(mfftPlan, forward, 0, outBuf, outBuf, normFact);			
	}
	dumpMatrixRect("fftExecute2DComplex output", outBuf, mfftPlan->numRows, mfftPlan->numCols);
	return mrtn;
}

#pragma mark --- 2-D complex FFT with internal transpose ---

/* 
 * Execute a 2-D square complex FFT, with internal transpose; 
 * frequency domain data is in column order.
 */
static MFFTReturn fftExecute2DComplexSquare(
	MatrixFFTPlan		mfftPlan,
	uint32_t			optFlags,
	bool				forward,		
	FFTComplex			*inBuf,
	FFTComplex			*outBuf)
{
	MFFTReturn mrtn = MR_Success;

	dumpMatrixRect("fftExecute2DComplexSquare input", inBuf, mfftPlan->numRows, mfftPlan->numCols);
	if(forward) {
		mrtn = feFftAllRows(mfftPlan, forward, 0, inBuf, outBuf, 0.0);
		if(mrtn) {
			return mrtn;
		}
		dumpMatrixRect("fftExecute2DComplexSquare(fwd) ceFftAllRows output", outBuf, 
			mfftPlan->numRows, mfftPlan->numCols);
		mrtn = fftTransposeSquare(mfftPlan, outBuf, mfftPlan->numRows);
		if(mrtn) {
			return mrtn;
		}
		dumpMatrixRect("fftExecute2DComplexSquare(fwd) transpose output", outBuf, 
			mfftPlan->numRows, mfftPlan->numCols);
		mrtn = feFftAllRows(mfftPlan, forward, 0, outBuf, outBuf, 0.0);
		if(mrtn) {
			return mrtn;
		}
		if(optFlags & MEF_TransposeOutput) {
			dumpMatrixRect("fftExecute2DComplexSquare(fwd) output pre-transpose", outBuf, 
				mfftPlan->numRows, mfftPlan->numCols);
			mrtn = fftTransposeSquare(mfftPlan, outBuf, mfftPlan->numRows);
		}
	}
	else {
		FFTFloat normFact = 0.0;
		
		if(optFlags & MEF_NormOutput) {
			normFact = 1.0 / ((FFTFloat)(mfftPlan->numRows * mfftPlan->numCols));
		}

		if(optFlags & MEF_TransposeInput) {
			mrtn = fftTransposeSquare(mfftPlan, inBuf, mfftPlan->numRows);
			dumpMatrixRect("fftExecute2DComplexSquare(inv) input post-transpose", outBuf, 
				mfftPlan->numRows, mfftPlan->numCols);
		}
		if(mrtn) {
			return mrtn;
		}
		mrtn = feFftAllRows(mfftPlan, forward, 0, inBuf, outBuf, 0.0);
		dumpMatrixRect("fftExecute2DComplexSquare(inv) feFftAllRows output", outBuf, 
			mfftPlan->numRows, mfftPlan->numCols);
		if(mrtn) {
			return mrtn;
		}
		mrtn = fftTransposeSquare(mfftPlan, outBuf, mfftPlan->numRows);
		if(mrtn) {
			return mrtn;
		}
		dumpMatrixRect("fftExecute2DComplexSquare(inv) transpose output", outBuf, 
			mfftPlan->numRows, mfftPlan->numCols);
		mrtn = feFftAllRows(mfftPlan, forward, 0, outBuf, outBuf, normFact);			
	}
	dumpMatrixRect("fftExecute2DComplexSquare output", outBuf, mfftPlan->numRows, mfftPlan->numCols);
	return mrtn;
}

#pragma mark --- generic 2-D transpose functions ---

/* 
 * Perform 2-D complex square transpose. 
 */
MFFTReturn fftTranspose2DComplexSquare(
	MatrixFFTPlan		mfftPlan,
	bool				forward,		
	bool				input, 
	const FFTComplex	*inBuf,
	FFTComplex			*outBuf)
{
	if(forward == input) {
		/* nop */
		return MR_Success;
	}
	
	/*
	 * There are two flavors; in-place is faster.
	 */
	RFASSERT(mfftPlan->numRows == mfftPlan->numColsComplex);
	
	if(inBuf == outBuf) {
		fftTransposeSquare(mfftPlan, outBuf, mfftPlan->numRows);
	}
	else {
		fftTransposeOP(mfftPlan, inBuf, outBuf, mfftPlan->numRows, mfftPlan->numRows);
	}
	
	return MR_Success;
}

#pragma mark --- 2-D Complex FFT via vDSP ---

#if     FFT_SPLIT_COMPLEX

/*
 * Split complex form, this is pretty easy...
 */
static MFFTReturn fftExec2DComplexVdsp(
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
        FFTComplex2d(mfftPlan->vdspSetup, &vbufIn, mfftPlan->nCol, mfftPlan->nRow, fdir);
    }
    else {
        FFTComplex2dOP(mfftPlan->vdspSetup, &vbufIn, &vbufOut, mfftPlan->nCol, mfftPlan->nRow, fdir);
    }
    if(!forward && (optFlags & MEF_NormOutput)) {
        size_t totalSamples = mfftPlan->numCols * mfftPlan->numRows;
        FFTFloat normFact = 1.0 / (FFTFloat)totalSamples;
		fftScaleComplex(outBuf, normFact, totalSamples);
    }
    return MR_Success;
}

#else   /* FFT_SPLIT_COMPLEX */

/* 
 * Interleaved, we have to convert to/from vDSP's native split format.
 */
static MFFTReturn fftExec2DComplexVdsp(
	MatrixFFTPlan		mfftPlan,
	uint32_t			optFlags,
	bool				forward,		
	FFTComplex			*inBuf,
	FFTComplex			*outBuf)
{
    RFASSERT(mfftPlan->vdspBuf.realp != NULL);
	FFTDirection fdir = forward ? FFT_FORWARD : FFT_INVERSE;
    size_t totalSamples = mfftPlan->numCols * mfftPlan->numRows;
    fftIntToVDSP(inBuf, &mfftPlan->vdspBuf, totalSamples);
    FFTComplex2d(mfftPlan->vdspSetup, &mfftPlan->vdspBuf, mfftPlan->nCol, mfftPlan->nRow, fdir);
    fftVDSPToInt(&mfftPlan->vdspBuf, outBuf, totalSamples);
    if(!forward && (optFlags & MEF_NormOutput)) {
        FFTFloat normFact = 1.0 / (FFTFloat)totalSamples;
		fftScaleComplex(outBuf, normFact, totalSamples);
    }
    return MR_Success;
}

#endif  /* FFT_SPLIT_COMPLEX */

#pragma mark --- config-specific MatrixFFTPlan init ---

/*
 * Set up a MatrixFFTPlan for raw vDSP-based FFTs.
 */
static MFFTReturn mfftCreatePlan2DComplexVdsp(
	unsigned			*n,			
	uint32_t			optFlags,
	MatrixFFTPlan		*mfftPlanRtn)
{
    MatrixFFTPlan mfftPlan = NULL;
    MFFTReturn mrtn = mfftPlanCreateCommon(2, 
            STT_None,   // no sine table
            false,      // real
            optFlags, 
            n[0], n[1], // nRows, nCols
            0,          // numThreads, actually will be ignored
            0,          // sinPeriod
            MFFT_PLAN_OPT_NO_THREADS | MFFT_PLAN_OPT_NO_ENGINES,
            &mfftPlan);
    if(mrtn) {
        return mrtn;
    }

    mfftPlan->inputFormat = MF_RowOrder;
    mfftPlan->transFcn    = fftNullTranspose;
    mfftPlan->execFcn	  = fftExec2DComplexVdsp;

    #if		!FFT_SPLIT_COMPLEX
    /* We need a vDSP buffer for conversion from interleaved */
    if(fftAllocDSPComplexAlign(&mfftPlan->vdspBuf, 
            mfftPlan->numCols * mfftPlan->numRows, VECTOR_SIZE, 
            &mfftPlan->vdspBufFree)) {
        mfftFreePlan(mfftPlan);
        return MR_Memory;
    }
    #endif

    *mfftPlanRtn = mfftPlan;
    return MR_Success;
}

MFFTReturn mfftCreatePlan2DComplex(
	unsigned			*n,			
	uint32_t			optFlags,
	unsigned			numThreads,
	MatrixFFTPlan		*mfftPlanRtn)
{
    if(mfftShouldUseVdsp(n[0] + n[1], false, 2)) {
        /* Set up a plan for vDSP only */
        return mfftCreatePlan2DComplexVdsp(n, optFlags, mfftPlanRtn);
    }

	MatrixFFTPlan mfftPlan = NULL;
	MFFTReturn mrtn = mfftPlanCreateCommon(2, STT_None, 
		false, optFlags, n[0], n[1], numThreads,
		0,				// sinPeriod not used
		MFFT_PLAN_OPT_NONE,	
		&mfftPlan);
	if(mrtn) {
		return mrtn;
	}
	
	if(FFT_2D_SQUARE_ENABLE &&
       (mfftPlan->nRow == mfftPlan->nCol) && 
       (mfftPlan->nRow >= FFT_2D_ROW_TRANS_SIZE)) {
		/* Large square, with internal transpose and output in column order */
		mfftPlan->inputFormat  = MF_RowOrder;
		mfftPlan->outputFormat = MF_ColumnOrder;
		mfftPlan->execFcn      = fftExecute2DComplexSquare;	
		mfftPlan->transFcn	   = fftTranspose2DComplexSquare;
	}
	else {
		/* FFT columns in place, input and output in row order */
		mfftPlan->inputFormat  = MF_RowOrder;
		mfftPlan->outputFormat = MF_RowOrder;
		mfftPlan->execFcn      = fftExecute2DComplex;	
		mfftPlan->transFcn	   = fftNullTranspose;
	}
	
	*mfftPlanRtn = mfftPlan;
	return MR_Success;
}

