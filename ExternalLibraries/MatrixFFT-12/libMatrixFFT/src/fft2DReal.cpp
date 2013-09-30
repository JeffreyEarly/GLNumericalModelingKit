/*	File: fft2DReal.cpp 
	
	Description:
		2-D Real-signal FFT functions.
	
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
 * fft2DReal.cpp - 2-D Real-signal FFT functions.
 *			  
 * Created 12/19/2008. 
 * Copyright 2008 by Apple, Inc. 
 */

#include <libMatrixFFT/MatrixFFT.h>
#include <libMatrixFFT/fftUtils.h>
#include <libMatrixFFT/complexBufUtils.h>
#include "MatrixFFTPlan.h"
#include "fftPriv.h"
#include "fftTranspose.h"
#include "fftEngineDispatch.h"
#include "fftDebug.h"
#include "fftPlatformConf.h"

/* 
 * For square 2-D FFTs at least 2^FFT_2D_ROW_TRANS_SIZE rows,
 * perform full square transpose and do the FFT-columns as rows, leaving
 * the result in column order.
 */
#if 0
/* force all sizes to same method for debug */
#define FFT_2D_ROW_TRANS_SIZE		100
#else
#if	FFT_DOUBLE_PREC
	#if		FFT_SPLIT_COMPLEX
		/* split double */
		#define FFT_2D_ROW_TRANS_SIZE	8
	#else
		/* interleaved double */
		#define FFT_2D_ROW_TRANS_SIZE	8
	#endif
#else
	#if		FFT_SPLIT_COMPLEX
		/* split single */
		#define FFT_2D_ROW_TRANS_SIZE	10
	#else
		/* interleaved single */
		#define FFT_2D_ROW_TRANS_SIZE	9
	#endif
#endif
#endif	/* 1 */

/* Turn this off to disable square transpose option */
#define FFT_2D_SQUARE_ENABLE        1

#if		FFT_SPLIT_COMPLEX

#pragma mark --- Routines to copy to/from column 0 as complexes ---

/* 
 * Copy the reals and then imaginaries from column 0 of inBuf as complexes to 
 * outBuf.
 */
void fftCopyFromColumn0(
	const FFTComplex *inBuf,
	FFTComplex *outBuf,
	size_t rowSize,			/* size of *inp row in FFTComplexes */
	size_t numComplex)		/* # of complex elements to move to each output row */
{
	FFTFloat *realp = outBuf->real;
	FFTFloat *imagp = outBuf->imag;
	
	/* 
	 * Reals in inBuf column 0 - DC components - first set of complex in output
	 */
	const FFTFloat *inp = inBuf->real;
	for(size_t dex=0; dex<numComplex; dex++) {
		*realp++ = *inp;
		inp     += rowSize;
		*imagp++ = *inp;
		inp     += rowSize;
	}
	
	/* 
	 * Imaginaries in inBuf column 0 - Nyquist components - 2nd set of complex in output
	 */
	inp = inBuf->imag;
	for(size_t dex=0; dex<numComplex; dex++) {
		*realp++ = *inp;
		inp     += rowSize;
		*imagp++ = *inp;
		inp     += rowSize;
	}

}

/*
 * Copy two sets of complexes from inBuf to the real and then the imaginary
 * components of column 0 of outBuf.
 */
void fftCopyToColumn0(
	const FFTComplex *inBuf,
	FFTComplex *outBuf,
	size_t rowSize,			/* size of *outp in FFTComplexes */
	size_t numComplex)		/* # of complex elements to move from each input row */
{
	FFTFloat *realp = inBuf->real;
	FFTFloat *imagp = inBuf->imag;
	
	/* 
	 * Reals in inBuf column 0 - DC components - first set of complex in input
	 */
	FFTFloat *outp = outBuf->real;
	for(size_t dex=0; dex<numComplex; dex++) {
		*outp = *realp++;
		outp += rowSize;
		*outp = *imagp++;
		outp += rowSize;
	}
	
	/* 
	 * Imaginaries in inBuf column 0 - Nyquist components - 2nd set of complex in input
	 */
	outp = outBuf->imag;
	for(size_t dex=0; dex<numComplex; dex++) {
		*outp = *realp++;
		outp += rowSize;
		*outp = *imagp++;
		outp += rowSize;
	}
}

#else	/* !FFT_SPLIT_COMPLEX */

/* 
 * Copy the reals and then imaginaries from column 0 of inBuf as complexes to 
 * outBuf.
 */
void fftCopyFromColumn0(
	const FFTComplex *inBuf,
	FFTComplex *outBuf,			/* two rows of numComplex elements */
	size_t rowSize,				/* size of *inBuf row in FFTComplexes */
	size_t numComplex)			/* # of complex elements to move to each output row */
{
	FFTComplex *dcOut  = outBuf;				// dest of all reals in inBuf
	FFTComplex *nyqOut = outBuf + numComplex;	// dest of all imaginaries in inBuf
	
	for(size_t dex=0; dex<numComplex; dex++) {
		dcOut->real  = inBuf->real;
		nyqOut->real = inBuf->imag;
		inBuf += rowSize;
		dcOut->imag  = inBuf->real;
		nyqOut->imag = inBuf->imag;
		inBuf += rowSize;
		dcOut++;
		nyqOut++;
	}
}

/*
 * Copy two sets of complexes from inBuf to the real and then the imaginary
 * components of column 0 of outBuf.
 */
void fftCopyToColumn0(
	const FFTComplex *inBuf,	/* two rows of numComplex elements */
	FFTComplex *outBuf,
	size_t rowSize,				/* size of *outBuf row in FFTComplexes */
	size_t numComplex)			/* # of complex elements to move from each input row */
{
	const FFTComplex *dcIn  = inBuf;					// source of of all reals in outBuf
	const FFTComplex *nyqIn = inBuf + numComplex;		// dest of all imaginaries in inBuf
	
	for(size_t dex=0; dex<numComplex; dex++) {
		outBuf->real = dcIn->real;
		outBuf->imag = nyqIn->real;
		outBuf += rowSize;
		outBuf->real = dcIn->imag;
		outBuf->imag = nyqIn->imag;
		outBuf += rowSize;
		dcIn++;
		nyqIn++;
	}
	
}
#endif	/* FFT_SPLIT_COMPLEX */

/* 
 * Special case for real-signal FFT of the first two columns. Since this operates on
 * half-column-size ops, and since it's small relative to the overall
 * signal, we don't clutter up the fftEngineDispatch with this case;
 * we just use vDSP directly here. 
 */
#if		FFT_SPLIT_COMPLEX
static void fftRealColumns(
	MatrixFFTPlan		mfftPlan,
	bool				forward,
	FFTComplex			*buf)				// always in place here - auxBuf
{
	size_t nRow = mfftPlan->nRow;			// N for the FFT
	size_t rowSize = 1 << (nRow-1);			// in complexes
	vDSPComplex vBuf;
	fftComplexToVDSP(buf, &vBuf);
	FFTDirection fdir = forward ? FFT_FORWARD : FFT_INVERSE;
	FFTReal1d(mfftPlan->vdspSetup, &vBuf, nRow, fdir);
	vBuf.realp += rowSize;
	vBuf.imagp += rowSize;
	FFTReal1d(mfftPlan->vdspSetup, &vBuf, nRow, fdir);
}

#else	/* !FFT_SPLIT_COMPLEX */

/* interleaved version */

static void fftRealColumns(
	MatrixFFTPlan		mfftPlan,
	bool				forward,
	FFTComplex			*buf)		// always in place here - auxBuf
{
	size_t nRow = mfftPlan->nRow;			// N for the FFT
	size_t rowSize = 1 << (nRow-1);			// in complexes
	FFTDirection fdir = forward ? FFT_FORWARD : FFT_INVERSE;
	RFASSERT(mfftPlan->vdspBuf.realp != NULL);
	
	for(unsigned dex=0; dex<2; dex++) {
		fftIntToVDSP(buf, &mfftPlan->vdspBuf, rowSize);
		FFTReal1d(mfftPlan->vdspSetup, &mfftPlan->vdspBuf, nRow, fdir);
		fftVDSPToInt(&mfftPlan->vdspBuf, buf, rowSize);
		buf += rowSize;
	}
}

#endif	/* FFT_SPLIT_COMPLEX */

#pragma mark --- 2-D real FFT, output in row order ---

/* 
 * Execute a 2-D Real-signal FFT with output in row order.
 * In-place or out-of-place (if in-place, inBuf == outBuf). 
 *
 * Scaling notes
 * -------------
 * This function produces output in the same form and with the same values
 * as the corresponding vDSP functions.
 *
 * Forward: output is 2x the actual value, the same as with vDSP. However we have
 * to scale the output of the 2nd real FFT on the DC and Nyquist columns by 0.5 since
 * that data is processed by two 1-D real-signal FFTs (which is the source of 
 * the 2x in the output). Thus the output of the 2nd real FFT is scaled by 0.5.
 *
 * Reverse output is (numRows * numCols) x the actual value, same as with vDSP.
 * We scale by 0.5 at the end to match vDSP.
 */
static MFFTReturn fftExecute2DReal(
	MatrixFFTPlan		mfftPlan,
	uint32_t			optFlags,
	bool				forward,		
	FFTComplex			*inBuf,
	FFTComplex			*outBuf)
{
	MFFTReturn mrtn = MR_Success;
	FFTComplex  *colSrc = inBuf;
	size_t totalElts = mfftPlan->numRows * mfftPlan->numCols;
	FFTFloat normFact = 0.0;
	
	if(optFlags & MEF_NormOutput) {
		normFact = 0.5 / (FFTFloat)totalElts;
	}

	dprintf(forward ? "===fftExecute2DReal forward top\n" : "===fftExecute2DReal inverse top\n");	
	dumpRectReal("fftExecute2DReal input", inBuf, mfftPlan->numRows, mfftPlan->numCols);
	
	/* 
	 * Alloc column-sized aux buf for FFTing the first two columns
	 */
	mrtn = fftAllocAuxBuf(mfftPlan, mfftPlan->numRows);
	if(mrtn) {
		return mrtn;
	}
	
	if(forward) {
		/* FFT all rows, real-->complex */
		mrtn = feFftAllRowsReal(mfftPlan, forward, inBuf, outBuf, 0.0);
		if(mrtn) {
			return mrtn;
		}
		colSrc = outBuf;
		dumpRectReal("fftExecute2DReal fft rows output", outBuf, mfftPlan->numRows, mfftPlan->numCols);
	}
	
	/*
	 * Here's the odd part. During a forward FFT, the first colummn in colSrc contains 
	 * the DC and Nyquist components of the output of feFftAllRowsReal in its
	 * real and imaginary components respectively. We need to perform a real-signal
	 * FFT on those values. Transpose over to auxBuf as two sets of complexes, 
	 * perform the FFT (on two elements, each of which is size numRows/2), and 
	 * transpose back to original position. 
	 * This works the same for forward and inverse FFTs. 
	 * For inverse out-of-place FFTs, here's where we switch from inBuf to 
	 * outBuf. For forward FFTs we already switched to outBuf. 
	 */
	size_t colSizeComplex = mfftPlan->numRows >> 1;
	size_t rowSizeComplex = mfftPlan->numCols >> 1;
	fftCopyFromColumn0(colSrc, mfftPlan->auxBuf, rowSizeComplex, colSizeComplex);
	dumpRectReal("fftExecute2DReal fromCol0 transpose", mfftPlan->auxBuf, 2, mfftPlan->numRows);
	fftRealColumns(mfftPlan, forward, mfftPlan->auxBuf);
	dumpRectReal("fftExecute2DReal col0 FFT output prescale", mfftPlan->auxBuf, 2, mfftPlan->numRows);
	
	/*
	 * Scaling - see comments above.
	 */
	if(forward) {
		fftScaleComplex(mfftPlan->auxBuf, 0.5, mfftPlan->numRows);
		dumpRectReal("fftExecute2DReal col0 FFT output postscale", mfftPlan->auxBuf, 2, mfftPlan->numRows);
	}
	fftCopyToColumn0(mfftPlan->auxBuf, outBuf, rowSizeComplex, colSizeComplex);
	dumpRectReal("fftExecute2DReal copyToCol0 output postscale", outBuf, mfftPlan->numRows,
		mfftPlan->numCols);
	
	/* 
	 * Now complex FFT of all remaining columns. 
	 * Note feFftColumns() handles any alignment problems here. 
	 */
	mrtn = feFftColumns(mfftPlan, forward, 
		1,		// startCol - in COMPLEX elts
		colSrc, outBuf,
		0.0);	// no normalization here
	if(mrtn) {
		return mrtn;
	}
	dumpRectReal("fftExecute2DReal feFftColumns output", outBuf, mfftPlan->numRows, mfftPlan->numCols);

	if(!forward) {
		/* FFT all rows, in-place in outBuf */
		mrtn = feFftAllRowsReal(mfftPlan, forward, outBuf, outBuf, normFact);
		dumpRectReal("fftExecute2DReal inverse FFT rows", outBuf, mfftPlan->numRows, mfftPlan->numCols);
		
		/* scale by 0.5 to match vDSP output */
		// done by feFftAllRowsReal()... fftScaleThr(mfftPlan, outBuf, 0.5, totalComplexElts);
	}
	
	dumpRectReal("fftExecute2DReal output", outBuf, mfftPlan->numRows, mfftPlan->numCols);
	return mrtn;
}

#pragma mark --- 2-D real FFT, output in column order ---

/*
 * Custom-order transpose for square 2-D real.
 */
static MFFTReturn fft2DRealTranspose(
	MatrixFFTPlan		mfftPlan,
	bool				forward,		
	bool				input, 
	const FFTComplex	*inBuf,
	FFTComplex			*outBuf)
{
	if(forward == input) {
		return MR_Success;
	}

	size_t colSizeComplex = mfftPlan->numRows >> 1;
	size_t rowSizeComplex = mfftPlan->numCols >> 1;
	RFASSERT(mfftPlan->numRows == mfftPlan->numColsComplex);
	MFFTReturn mrtn;
	
	/* 
	 * Alloc column-sized aux buf for saving the first two columns
	 */
	mrtn = fftAllocAuxBuf(mfftPlan, mfftPlan->numRows);
	if(mrtn) {
		return mrtn;
	}

	if(forward) {
		/*
		 * Row 0 contains real-signal FFTs of the DC and Nyquist components
		 * of the per-row real-signal FFTs. Save a copy...
		 */
		fftCopyComplexArray(inBuf, mfftPlan->auxBuf, rowSizeComplex);
		
		/* Normal square transpose on the whole array - might be OOP! */
		mrtn = fftTranspose2DComplexSquare(mfftPlan, true, false, inBuf, outBuf);
		if(mrtn) {
			return mrtn;
		}
		
		/* Transpose auxBuf to column 0 of output */
		fftCopyToColumn0(mfftPlan->auxBuf, outBuf, rowSizeComplex, colSizeComplex);
	}
	else {
		/*
		 * Input to inverse FFT. First save a transposed copy of column 0.
		 */
		fftCopyFromColumn0(inBuf, mfftPlan->auxBuf, rowSizeComplex, colSizeComplex);

		/* Normal square transpose on the whole array - might be OOP! */
		mrtn = fftTranspose2DComplexSquare(mfftPlan, false, true, inBuf, outBuf);
		if(mrtn) {
			return mrtn;
		}
		
		/* Copy column 0 to row 0 */
		fftCopyComplexArray(mfftPlan->auxBuf, outBuf, rowSizeComplex);
	}
	return MR_Success;
}

/* 
 * Execute a 2-D Real-signal FFT with output in column order.
 * In-place or out-of-place (if in-place, inBuf == outBuf). 
 */
static MFFTReturn fftExecute2DRealSquare(
	MatrixFFTPlan		mfftPlan,
	uint32_t			optFlags,
	bool				forward,		
	FFTComplex			*inBuf,
	FFTComplex			*outBuf)
{
	MFFTReturn mrtn = MR_Success;
	RFASSERT(mfftPlan->numRows == mfftPlan->numColsComplex);
	
	dprintf(forward ? "===fftExecute2DRealSquare forward top\n" : "===fftExecute2DReal inverse top\n");	
	dumpRectReal("fftExecute2DRealSquare input", inBuf, mfftPlan->numRows, mfftPlan->numCols);

	size_t totalElts = mfftPlan->numRows * mfftPlan->numCols;
	FFTFloat normFact = 0.0;
	
	if(optFlags & MEF_NormOutput) {
		normFact = 0.5 / (FFTFloat)totalElts;
	}

	/* 
	 * Alloc column-sized aux buf for FFTing the first two columns
	 */
	mrtn = fftAllocAuxBuf(mfftPlan, mfftPlan->numRows);
	if(mrtn) {
		return mrtn;
	}
	size_t colSizeComplex = mfftPlan->numRows >> 1;
	size_t rowSizeComplex = mfftPlan->numCols >> 1;

	if(forward) {
		/* FFT all rows, real-->complex */
		mrtn = feFftAllRowsReal(mfftPlan, forward, inBuf, outBuf, 0.0);
		if(mrtn) {
			return mrtn;
		}
		dumpRectReal("fftExecute2DRealSquare fft rows output", outBuf, mfftPlan->numRows, mfftPlan->numCols);
		
		/* 
		 * Transpose column 0, containing DC and Nyquist components of the FFT-rows op,
		 * to aux buf. 
		 */
		fftCopyFromColumn0(outBuf, mfftPlan->auxBuf, rowSizeComplex, colSizeComplex);
		dumpRectReal("fftExecute2DRealSquare fromCol0 transpose", mfftPlan->auxBuf, 2, mfftPlan->numRows);
		
		/* 
		 * Real-FFT those two half-rows. Keep result in auxBuf.
		 */
		fftRealColumns(mfftPlan, true, mfftPlan->auxBuf);
		dumpRectReal("fftExecute2DRealSquare col0 FFT output prescale", mfftPlan->auxBuf, 2, 
			mfftPlan->numRows);

		/*
		 * Scaling - see comments above in fftExecute2DReal().
		 */
		fftScaleComplex(mfftPlan->auxBuf, 0.5, mfftPlan->numRows);
		dumpRectReal("fftExecute2DRealSquare col0 FFT output postscale", mfftPlan->auxBuf, 2, 
			mfftPlan->numRows);

		/* 
		 * Transpose the whole signal as one square - including the first column, because it's 
		 * faster to avoid the special unaligned case if we didn't transpose it. We're going to 
		 * overwrite the result of the transposed first column with auxBuf after the transpose.
		 */
		mrtn = fftTransposeSquare(mfftPlan, outBuf, mfftPlan->numRows);
		if(mrtn) {
			return mrtn;
		}
		dumpRectReal("fftExecute2DRealSquare int transpose output", outBuf, 
			mfftPlan->numRows, mfftPlan->numCols);

		/* 
		 * Complex FFT all rows after the first one. 
		 */
		mrtn = feFftAllRows(mfftPlan, true, 1, outBuf, outBuf, 0.0);
		dumpRectReal("fftExecute2DRealSquare fft cols output", outBuf, mfftPlan->numRows, 
			mfftPlan->numCols);
		
		
		/* Now copy over the two half-rows to the first output row. */
		fftCopyComplexArray(mfftPlan->auxBuf, outBuf, rowSizeComplex);
		
		if(optFlags & MEF_TransposeOutput) {
			dumpRectReal("fftExecute2DRealSquare(fwd) output pre-transpose", outBuf, 
				mfftPlan->numRows, mfftPlan->numCols);
			mrtn = fft2DRealTranspose(mfftPlan, true, false, outBuf, outBuf);
		}
	}
	else {
		/*
		 * Optionally transpose input from row order to custom. 
		 */
		if(optFlags & MEF_TransposeInput) {
			mrtn = fft2DRealTranspose(mfftPlan, false, true, inBuf, inBuf);
			dumpRectReal("fftExecute2DRealSquare input post-transpose", outBuf, 
				mfftPlan->numRows, mfftPlan->numCols);
		}
		
		/*
		 * Copy the first row to auxBuf. Inverse real FFT in place in auxBuf, save
		 * result for later. 
		 */
		fftCopyComplexArray(inBuf, mfftPlan->auxBuf, rowSizeComplex);
		fftRealColumns(mfftPlan, false, mfftPlan->auxBuf);
		
		/* 
		 * Complex FFT all rows after the first one. 
		 */
		mrtn = feFftAllRows(mfftPlan, false, 1, inBuf, outBuf, 0.0);
		dumpRectReal("fftExecute2DRealSquare fft cols output", outBuf, mfftPlan->numRows, 
			mfftPlan->numCols);
		
		/*
		 * Transpose the whole square, including the first row, which we'll overwrite as 
		 * the first column with data from auxBuf.
		 */
		mrtn = fftTransposeSquare(mfftPlan, outBuf, mfftPlan->numRows);
		if(mrtn) {
			return mrtn;
		}
		dumpRectReal("fftExecute2DRealSquare after transpose", outBuf, mfftPlan->numRows, 
			mfftPlan->numCols);
		
		/*
		 * Drop in the DC and Nyquist components to the first column. 
		 */
		fftCopyToColumn0(mfftPlan->auxBuf, outBuf, rowSizeComplex, colSizeComplex);
		
		/*
		 * Finally, real-FFT all rows. Result is in row order.
		 */
		mrtn = feFftAllRowsReal(mfftPlan, false, outBuf, outBuf, normFact);
		if(mrtn) {
			return mrtn;
		}
		
		/* scale by 0.5 to match vDSP output */
		//fftScaleThr(mfftPlan, outBuf, 0.5, mfftPlan->numRows * mfftPlan->numColsComplex);
	}
	dumpRectReal("fftExecute2DRealSquare output", outBuf, 
		mfftPlan->numRows, mfftPlan->numCols);
	return MR_Success;
}

#pragma mark --- 2-D Real FFT via vDSP ---

#if     FFT_SPLIT_COMPLEX

/*
 * Split complex form, this is pretty easy...
 */
static MFFTReturn fftExec2DRealVdsp(
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
        FFTReal2d(mfftPlan->vdspSetup, &vbufIn, mfftPlan->nCol, mfftPlan->nRow, fdir);
    }
    else {
        FFTReal2dOP(mfftPlan->vdspSetup, &vbufIn, &vbufOut, mfftPlan->nCol, mfftPlan->nRow, fdir);
    }
    if(!forward && (optFlags & MEF_NormOutput)) {
        /* scale by 1/2N */
        size_t totalRealSamples = mfftPlan->numCols * mfftPlan->numRows;
        FFTFloat normFact = 0.5 / (FFTFloat)totalRealSamples;
		fftScaleComplex(outBuf, normFact, totalRealSamples >> 1);
    }
    return MR_Success;
}

#else   /* FFT_SPLIT_COMPLEX */

/* 
 * Interleaved, we have to convert to/from vDSP's native split format.
 */
static MFFTReturn fftExec2DRealVdsp(
	MatrixFFTPlan		mfftPlan,
	uint32_t			optFlags,
	bool				forward,		
	FFTComplex			*inBuf,
	FFTComplex			*outBuf)
{
    RFASSERT(mfftPlan->vdspBuf.realp != NULL);
	FFTDirection fdir = forward ? FFT_FORWARD : FFT_INVERSE;
    size_t totalRealSamples = mfftPlan->numCols * mfftPlan->numRows;
    size_t totalComplexSamples = totalRealSamples >> 1;
    fftIntToVDSP(inBuf, &mfftPlan->vdspBuf, totalComplexSamples);
    FFTReal2d(mfftPlan->vdspSetup, &mfftPlan->vdspBuf, mfftPlan->nCol, mfftPlan->nRow, fdir);
    fftVDSPToInt(&mfftPlan->vdspBuf, outBuf, totalComplexSamples);
    if(!forward && (optFlags & MEF_NormOutput)) {
        /* scale by 1/2N */
        FFTFloat normFact = 0.5 / (FFTFloat)totalRealSamples;
		fftScaleComplex(outBuf, normFact, totalComplexSamples);
    }
    return MR_Success;
}

#endif  /* FFT_SPLIT_COMPLEX */

#pragma mark --- 2-D Real MatrixFFTPlan init ---

/*
 * Set up a MatrixFFTPlan for raw vDSP-based FFTs.
 */
static MFFTReturn mfftCreatePlan2DRealVdsp(
	unsigned			*n,			
	uint32_t			optFlags,
	MatrixFFTPlan		*mfftPlanRtn)
{
    MatrixFFTPlan mfftPlan = NULL;
    MFFTReturn mrtn = mfftPlanCreateCommon(2, 
            STT_None,   // no sine table
            true,       // real
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
    mfftPlan->execFcn	  = fftExec2DRealVdsp;

    #if		!FFT_SPLIT_COMPLEX
    /* We need a vDSP buffer for conversion from interleaved */
    if(fftAllocDSPComplexAlign(&mfftPlan->vdspBuf, 
            (mfftPlan->numCols * mfftPlan->numRows) >> 1, VECTOR_SIZE, 
            &mfftPlan->vdspBufFree)) {
        mfftFreePlan(mfftPlan);
        return MR_Memory;
    }
    #endif

    *mfftPlanRtn = mfftPlan;
    return MR_Success;
}

MFFTReturn mfftCreatePlan2DReal(
	unsigned			*n,			
	uint32_t			optFlags,
	unsigned			numThreads,
	MatrixFFTPlan		*mfftPlanRtn)
{
    if(mfftShouldUseVdsp(n[0] + n[1], true, 2)) {
        /* Set up a plan for vDSP only */
        return mfftCreatePlan2DRealVdsp(n, optFlags, mfftPlanRtn);
    }

	MatrixFFTPlan mfftPlan = NULL;
	MFFTReturn mrtn = mfftPlanCreateCommon(2, STT_None, 
		true, optFlags, n[0], n[1], numThreads,
		0,				// sinPeriod not used
		MFFT_PLAN_OPT_NONE,	
		&mfftPlan);
	if(mrtn) {
		return mrtn;
	}
	
	if(FFT_2D_SQUARE_ENABLE &&
      (mfftPlan->nRow == mfftPlan->nColComplex) && 
      (mfftPlan->nRow >= FFT_2D_ROW_TRANS_SIZE)) {
		/* Large square, with internal transpose and output in custom order */
		mfftPlan->inputFormat  = MF_RowOrder;
		mfftPlan->outputFormat = MF_CustomOrder;
		mfftPlan->execFcn      = fftExecute2DRealSquare;	
		mfftPlan->transFcn	   = fft2DRealTranspose;
	}
	else {
		/* FFT columns in place, input and output in row order */
		mfftPlan->inputFormat  = MF_RowOrder;
		mfftPlan->outputFormat = MF_RowOrder;
		mfftPlan->execFcn      = fftExecute2DReal;	
		mfftPlan->transFcn	   = fftNullTranspose;
	}
	
	#if		!FFT_SPLIT_COMPLEX
	/* We need a vDSP buffer for fftRealColumns() */
	if(fftAllocDSPComplexAlign(&mfftPlan->vdspBuf, mfftPlan->numRows, VECTOR_SIZE, 
			&mfftPlan->vdspBufFree)) {
		mfftFreePlan(mfftPlan);
		return MR_Memory;
	}
	#endif
	
	*mfftPlanRtn = mfftPlan;
	return MR_Success;
}
