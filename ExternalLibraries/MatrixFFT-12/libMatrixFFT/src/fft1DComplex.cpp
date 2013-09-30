/*	File: fft1DComplex.cpp 
	
	Description:
		1-D Complex FFT via matrix decomposition.
	
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
 * fft1DComplex.cpp - 1-D Complex FFT via matrix decomposition. 
 *			  
 * Created 9/24/2008. 
 * Copyright 2008 by Apple, Inc. 
 */

#include <libMatrixFFT/MatrixFFT.h>
#include <libMatrixFFT/fftUtils.h>
#include <libMatrixFFT/complexBufUtils.h>
#include <libMatrixFFT/src/fftPlatformConf.h>       /* private SPI */
#include "MatrixFFTPlan.h"
#include "fftPriv.h"
#include "fftTranspose.h"
#include "fftSinCos.h"
#include "fftEngineDispatch.h"
#include "fftDebug.h"
#include "fftPlatformConf.h"

#pragma mark --- Forward 1-D Complex FFT ---

/* 
 * Execute a forward 1-D complex FFT. 
 * In-place or out-of-place (if in-place, inBuf == outBuf). 
 */
static MFFTReturn fftExec1DComplexFwd(
	MatrixFFTPlan		mfftPlan,
	uint32_t			optFlags,
	FFTComplex			*inBuf,
	FFTComplex			*outBuf)
{
	/* 
	 * input in row order;
	 * FFT each column in place;
	 * <<no transpose>>
	 * perform the twist in place;
	 * FFT each row in place;
	 * output is in column order;
	 */
	 
	dumpMatrixRect("Exec1DComplexFwd input", inBuf, 
		mfftPlan->numRows, mfftPlan->numCols);
	
	/* FFT each column inBuf --> outBuf */
	MFFTReturn mrtn = feFftColumns(mfftPlan, true, 0, inBuf, outBuf, 0.0);
	if(mrtn) {
		return mrtn;
	}

	dumpMatrixRect("Exec1DComplexFwd after FFT each column", outBuf, 
		mfftPlan->numRows, mfftPlan->numCols);
	
	/* Twist and FFT each row, a stripe at a time */
	mrtn = feTwistPlusFftRow(mfftPlan, outBuf, outBuf);
	if(mrtn) {
		return mrtn;
		
	}
	dumpMatrixRect("Exec1DComplexFwd after twist and FFT each row", outBuf, 
		mfftPlan->numRows, mfftPlan->numCols);
	
	if(optFlags & MEF_TransposeOutput) {
		mrtn = mfftPlan->transFcn(mfftPlan, true, false, outBuf, outBuf);
	}
	return mrtn;
}

#pragma mark --- Inverse 1-D Complex FFT ---

/* 
 * Execute a reverse 1-D complex FFT. 
 * In-place or out-of-place (if in-place, inBuf == outBuf). 
 */
static MFFTReturn fftExec1DComplexInv(
	MatrixFFTPlan		mfftPlan,
	uint32_t			optFlags,
	FFTComplex			*inBuf,
	FFTComplex			*outBuf)
{
	/* 
	 * input in column order;
	 * FFT(inv) each row inBuf-->outBuf;
	 * perform the inverse twist in place;
	 * FFT(Inv) each column in place;
	 * output is in row order;
	 */
	MFFTReturn mrtn = MR_Success;
	if(optFlags & MEF_TransposeInput) {
		dumpMatrixRect("Exec1DComplexInv pre-transpose", inBuf, 
			mfftPlan->numRows, mfftPlan->numCols);
		mrtn = mfftPlan->transFcn(mfftPlan, false, true, inBuf, inBuf);
		if(mrtn) {
			return mrtn;
		}
	}
	dumpMatrixRect("Exec1DComplexInv input", inBuf, 
		mfftPlan->numRows, mfftPlan->numCols);

	/* FFT and twist each row, a stripe at a time */
	mrtn = feFftRowPlusInvTwist(mfftPlan, inBuf, outBuf);
	if(mrtn) {
		return mrtn;
	}
	dumpMatrixRect("Exec1DComplexInv after FFT and twist each row", outBuf, 
		mfftPlan->numRows, mfftPlan->numCols);
	
	/* FFT each column with optional normalization */
	FFTFloat normFact = 0.0;
	if(optFlags & MEF_NormOutput) {
		normFact = 1.0 / (FFTFloat)(mfftPlan->numRows * mfftPlan->numCols);
		if(mfftPlan->addNormFact != 0.0) {
			normFact *= mfftPlan->addNormFact;
		}
	}
	
	mrtn = feFftColumns(mfftPlan, false, 0, outBuf, outBuf, normFact);
	if(mrtn) {
		return mrtn;
	}

	dumpMatrixRect("Exec1DComplexInv after FFT each column", outBuf, 
		mfftPlan->numRows, mfftPlan->numCols);
	
	return mrtn;
}

#pragma mark --- 1-D nonsquare transpose ---

/* 
 * Perform 1-D complex nonsquare transpose. 
 * This has to be out of place; we'll malloc the auxBuf if
 * we have to (and overwrite its contents if it's already present). 
 *
 * We have to swap rows and columns here when transposing.
 * Although the output of the forward FFT is considered to be in column 
 * order, for things to work right, the transpose back to row order has 
 * to use the original coordinates for the source.
 *
 * For the record here is what happens during a 1-D FFT of 32 elements:
 *
 *	-- initial rectangle: numRows=4  numCols=8
 *	-- input to forward FFT is in normal row order
 *	-- FFT each column; totalCols=8
 *  -- Twist & FFT each row in fft1DFwdTwistFft. Still numRows=4, numCols=8. 
 *	-- Here's the oddity. Although we're still looking at data with 4 rows
 *	   and 8 columns, we have to transpose it, from a rows=4 x cols=8 
 *     matrix to a rows=8, cols=4 matrix to match vDSP/output. 
 *
 * For the inverse, assume the data is in "row order" i.e. it matches other
 * FFT implementations. 
 * 
 *	-- data comes in as rows=8 x cols=4 rectangle; we transpose it to 
 *	   rows=4 x cols=8.
 *	-- ALL subsequent ops operate on those coordinates. 
 */
MFFTReturn fftTranspose1DComplexRect(
	MatrixFFTPlan		mfftPlan,
	bool				forward,		
	bool				input, 
	const FFTComplex	*inBuf,
	FFTComplex			*outBuf)
{
	if(forward == input) {
		/* nop, but provide for OOP anyway */
		return fftNullTranspose(mfftPlan, forward, input, inBuf, outBuf);;
	}

	size_t srcRows;
	size_t srcCols;
	
	if(forward) {
		/*
		 * Data is in column order, viewing expected output with 
		 * rows and columns swapped (i.e. the input - column order data -
		 * has the same coordinates as the original time-domain data).
		 */
		srcRows = mfftPlan->numRows;
		srcCols = mfftPlan->numColsComplex;
	}
	else {
		/*
		 * Converting from row order to column order: input has rows and 
		 * columns swapped.
		 */
		srcRows = mfftPlan->numColsComplex;
		srcCols = mfftPlan->numRows;
	}
	
	FFTComplex *dst = outBuf;
	MFFTReturn mrtn = MR_Success;
	size_t numElts = srcRows * srcCols;
	
	if(inBuf == outBuf) {
		mrtn = fftAllocAuxBuf(mfftPlan, numElts);
		if(mrtn) {
			return mrtn;
		}
		dst = mfftPlan->auxBuf;
	}
	fftTransposeOP(mfftPlan, inBuf, dst, srcRows, srcCols);
	if(inBuf == outBuf) {
		fftCopyComplexArray(dst, outBuf, numElts);
	}
	return MR_Success;
}

static MFFTReturn fftExec1DComplex(
	MatrixFFTPlan		mfftPlan,
	uint32_t			optFlags,
	bool				forward,		
	FFTComplex			*inBuf,
	FFTComplex			*outBuf)
{
	if(forward) {
		return fftExec1DComplexFwd(mfftPlan, optFlags, inBuf, outBuf);
	}
	else {
		return fftExec1DComplexInv(mfftPlan, optFlags, inBuf, outBuf);
	}
}

#pragma mark --- 1-D Complex FFT via vDSP ---

#if     FFT_SPLIT_COMPLEX

/*
 * Split complex form, this is pretty easy...
 */
static MFFTReturn fftExec1DComplexVdsp(
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
        FFTComplex1d(mfftPlan->vdspSetup, &vbufIn, mfftPlan->nCol, fdir);
    }
    else {
        FFTComplex1dOP(mfftPlan->vdspSetup, &vbufIn, &vbufOut, mfftPlan->nCol, fdir);
    }
    if(!forward && (optFlags & MEF_NormOutput)) {
        FFTFloat normFact = 1.0 / (FFTFloat)(mfftPlan->numCols);
		fftScaleComplex(outBuf, normFact, mfftPlan->numCols);
    }
    return MR_Success;
}

#else   /* FFT_SPLIT_COMPLEX */

/* 
 * Interleaved, we have to convert to/from vDSP's native split format.
 */
static MFFTReturn fftExec1DComplexVdsp(
	MatrixFFTPlan		mfftPlan,
	uint32_t			optFlags,
	bool				forward,		
	FFTComplex			*inBuf,
	FFTComplex			*outBuf)
{
    RFASSERT(mfftPlan->vdspBuf.realp != NULL);
	FFTDirection fdir = forward ? FFT_FORWARD : FFT_INVERSE;
    fftIntToVDSP(inBuf, &mfftPlan->vdspBuf, mfftPlan->numCols);
    FFTComplex1d(mfftPlan->vdspSetup, &mfftPlan->vdspBuf, mfftPlan->nCol, fdir);
    fftVDSPToInt(&mfftPlan->vdspBuf, outBuf, mfftPlan->numCols);
    if(!forward && (optFlags & MEF_NormOutput)) {
        FFTFloat normFact = 1.0 / (FFTFloat)(mfftPlan->numCols);
		fftScaleComplex(outBuf, normFact, mfftPlan->numCols);
    }
    return MR_Success;
}

#endif  /* FFT_SPLIT_COMPLEX */

#pragma mark --- config-specific MatrixFFTPlan init ---

/*
 * Set up a MatrixFFTPlan for raw vDSP-based FFTs.
 */
static MFFTReturn mfftCreatePlan1DComplexVdsp(
	unsigned			*n,			
	uint32_t			optFlags,
	MatrixFFTPlan		*mfftPlanRtn)
{
    MatrixFFTPlan mfftPlan = NULL;
    MFFTReturn mrtn = mfftPlanCreateCommon(1, 
            STT_None,   // no sine table
            false,      // real
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
    mfftPlan->execFcn	  = fftExec1DComplexVdsp;

    #if		!FFT_SPLIT_COMPLEX
    /* We need a vDSP buffer for conversion from interleaved */
    if(fftAllocDSPComplexAlign(&mfftPlan->vdspBuf, 
            mfftPlan->numCols, VECTOR_SIZE, 
            &mfftPlan->vdspBufFree)) {
        mfftFreePlan(mfftPlan);
        return MR_Memory;
    }
    #endif

    *mfftPlanRtn = mfftPlan;
    return MR_Success;
}

/* 
 * This is also called by other FFT plan setups, thus we have some peculiar options.
 *
 * When disableThreads is true, caller has set up a thread pool elsewhere
 * and will "hook up" our threadPoolP to that external thread pool. 
 *
 * When externalSineTable is true we don't have to set up a sine lookup
 * table; caller will connect our externSineTable to its MatrixFFTPlan.
 *
 * When expectColOut is true, caller expects to see the output of our
 * forward FFT in column output, so we should *not* attempt to use
 * vDSP FFTs for small signals (since the output of that would be in
 * row order).
 */
MFFTReturn mfftCreatePlan1DComplex(
	unsigned			*n,			
	uint32_t			optFlags,
	unsigned			numThreads,
	bool				externalSineTable,	/* don't create sine tables */
	bool				disableThreads,		/* skip ThreadPool init */
    bool                expectColOut,       /* Don't optimize by using vDSP */
	MatrixFFTPlan		*mfftPlanRtn)
{
    if(!expectColOut) {
        if(mfftShouldUseVdsp(n[0], false, 1)) {
            /* Set up a plan for vDSP only */
            return mfftCreatePlan1DComplexVdsp(n, optFlags, mfftPlanRtn);
        }
    }
    
	/* Divvy up into a configuration-dependent rectangle */
	unsigned nRows;
	unsigned nCols;
	mfftMakeRectangle(n[0], optFlags, &nRows, &nCols);
	
	MatrixFFTPlan mfftPlan = NULL;
	
	FFTSineTableType tableType = STT_Standard;
	size_t sinPeriod = 0;
	size_t numRows = (size_t)1 << nRows;
	if(externalSineTable) {
		/* Currently used when 1-D real uses us as a subplan */
		tableType = STT_External;
	}
	else {
		if((numRows > FFT_FLOATS_PER_VECTOR) && !FFT_FWD_TWIST_COMPLEX_SMALL && FFT_INTEL) {
			/* go for the optimized sin table */
			tableType = STT_Optimized;
			sinPeriod = FFT_SIN_RECALC_COMPLEX;
		}
	}
	MFFTReturn mrtn = mfftPlanCreateCommon(1, tableType,
		false, optFlags, nRows, nCols, numThreads,
		sinPeriod,
		disableThreads ? MFFT_PLAN_OPT_NO_THREADS : MFFT_PLAN_OPT_NONE,
		&mfftPlan);
	if(mrtn) {
		return mrtn;
	}

	mfftPlan->inputFormat   = MF_RowOrder;
	mfftPlan->outputFormat  = MF_ColumnOrder;
	if(nRows == nCols) {
		mfftPlan->transFcn  = fftTranspose2DComplexSquare;
	}
	else {
		mfftPlan->transFcn  = fftTranspose1DComplexRect;
	}
	mfftPlan->execFcn	    = fftExec1DComplex;

	*mfftPlanRtn = mfftPlan;
	return MR_Success;
}

