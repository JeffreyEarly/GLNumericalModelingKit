/*	File: FFTEngineVDSP.cpp  
	
	Description:
		Accelerate.framework-based FFTEngine.
	
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
 * FFTEngineVDSP.cpp - Accelerate.framework-based FFTEngine. 
 *			  
 * Created 12/16/2008. 
 * Copyright 2008 by Apple, Inc. 
 */

#include "FFTEngineVDSP.h"
#include "fftPriv.h"
#include "fftTranspose.h"
#include "fftDebug.h"
#include <libMatrixFFT/fftUtils.h>
#include "fftPlatformConf.h"

#ifdef	DEBUG
/* sanity check, slower */
#define SAFE_CAST(thing, type)	dynamic_cast<type>(thing)
#else
/* faster */
#define SAFE_CAST(thing, type)	static_cast<type>(thing)
#endif	/* DEBUG */


/* 
 * Stripe size for "depth-first" FFT and twist.
 * This doesn't affect any memory allocation; it's the quanta we break up
 * Twist/FFT ops into. It's always best to keep this to the minimum possible
 * value of 1, to maximize caching of the data being processed across the 
 * two ops.
 */
#define VDSP_DEPTH_FIRST_SIZE		1

/* Set to 1 to check buffer alignment */
#define HOST_CHECK_ALIGN			0

#if		HOST_CHECK_ALIGN

static void inline fftCheckAlign(size_t numElts, void *p)
{
	if((1 << numElts) < FFT_FLOATS_PER_SUBMATRIX) {
		return;
	}
	intptr_t ip = (intptr_t)p;
	if(ip & (FFT_FLOATS_PER_SUBMATRIX - 1)) {
		printf("***Hey! Unaligned ptr (%p) with numElts %llu\n",
			p, (unsigned long long)numElts);
		exit(1);
	}
}

#else	/* HOST_CHECK_ALIGN */

#define fftCheckAlign(n, p)

#endif	/* HOST_CHECK_ALIGN */

#if		FFT_SPLIT_COMPLEX

/* Split complex version */
MFFTReturn FFTEngineVDSP::fftRowsAtom(
	bool				forward,		
	FFTComplex			*inBuf,
	FFTComplex			*outBuf,
	FFTFloat			normFact,
	size_t				startRow, 
	size_t				rowsToProcess,
	unsigned			log2Size)		/* can be numRows() when called by fftColumnsAtom() */
{
	size_t rowSize = 1 << log2Size;
	size_t offset  = rowSize * startRow;
	vDSPComplex vinBuf;
	vDSPComplex voutBuf;
	bool inPlace = (inBuf == outBuf) || (inBuf->real == outBuf->real);
	bool doScale = (normFact != 0.0);

	fftComplexToVDSP(inBuf, &vinBuf);
	if(!inPlace) {
		fftComplexToVDSP(outBuf, &voutBuf);
	}
	vinBuf.realp += offset;
	vinBuf.imagp += offset;
	if(!inPlace) {
		voutBuf.realp += offset;
		voutBuf.imagp += offset;
	}
	FFTDirection fdir = forward ? FFT_FORWARD : FFT_INVERSE;
		
	for(size_t dex=0; dex<rowsToProcess; dex++) {
		if(inPlace) {
			FFTComplex1d(mFftSetup, &vinBuf, log2Size, fdir);
			if(doScale) {
				fftScaleDSPComplex(&vinBuf, normFact, rowSize);
			}
		}
		else {
			FFTComplex1dOP(mFftSetup, &vinBuf, &voutBuf, log2Size, fdir);
			if(doScale) {
				fftScaleDSPComplex(&voutBuf, normFact, rowSize);
			}
			voutBuf.realp += rowSize;
			voutBuf.imagp += rowSize;
		}
		vinBuf.realp += rowSize;
		vinBuf.imagp += rowSize;
	}
	
	return MR_Success;
}

#else /* Interleaved version */

MFFTReturn FFTEngineVDSP::fftRowsAtom(
	bool				forward,		
	FFTComplex			*inBuf,
	FFTComplex			*outBuf,
	FFTFloat			normFact,
	size_t				startRow, 
	size_t				rowsToProcess,
	unsigned			log2Size)		/* can be numRows() when called by fftColumnsAtom() */
{
	size_t rowSize      = 1 << log2Size;
	size_t offset		= rowSize * startRow;
	inBuf += offset;
	outBuf += offset;
	FFTDirection fdir   = forward ? FFT_FORWARD : FFT_INVERSE;
	bool doScale		= (normFact != 0.0);
	
	RFASSERT(mVdspBuf.realp != NULL);
	
	for(size_t dex=0; dex<rowsToProcess; dex++) {
		fftIntToVDSP(inBuf, &mVdspBuf, rowSize);
		FFTComplex1d(mFftSetup, &mVdspBuf, log2Size, fdir);
		if(doScale) {
			fftScaleDSPComplex(&mVdspBuf, normFact, rowSize);
		}
		fftVDSPToInt(&mVdspBuf, outBuf, rowSize);
		inBuf  += rowSize;
		outBuf += rowSize;
	}
	
	return MR_Success;
}

#endif	/* FFT_SPLIT_COMPLEX */
	
#pragma mark --- Public SPI ---

FFTEngineVDSP::FFTEngineVDSP(
	MatrixFFTPlan		_mfftPlan,
	unsigned			dims,
	size_t				_numRows,
	size_t				_numCols,
	bool				_real)
		:
		FFTEngine(_mfftPlan, _numRows, _numCols, _real),
		mFftSetup(_mfftPlan->vdspSetup),
		mLog2NumRows(0), mLog2NumCols(0), mLog2NumColsComplex(0),
		mColStripeSize(0)
		#if		!FFT_SPLIT_COMPLEX
		, mVdspBuf(), mVdspBufFree()
		#endif
{
	fftIsPowerOfTwo(_numRows, &mLog2NumRows);
	fftIsPowerOfTwo(_numCols, &mLog2NumCols);
	if(_real) {
		mLog2NumColsComplex = mLog2NumCols - 1;
	}
	else {
		mLog2NumColsComplex = mLog2NumCols;
	}
	
    /* 
     * mfftColumnStripeSize takes log2(total signal size) for 1-D and logs(numRows) for 2-D
     */
    unsigned l2n = (dims == 1) ? (mLog2NumRows + mLog2NumCols) : mLog2NumRows;
    mColStripeSize = mfftColumnStripeSize(l2n, _real, dims, _mfftPlan->l2CachePerEngine);
	size_t stripeSize = _numRows * mColStripeSize;
	
	/* stripe buffer, aligned to submatrix size */
	MFFTReturn mrtn = allocStripeBufs(1, stripeSize, CACHE_LINE_SIZE);
	if(mrtn) {
		FFTError::throwMe(mrtn);
	}
	fftAssertBufAlign(mStripeBuf[0], CACHE_LINE_SIZE);
	
	#if		!FFT_SPLIT_COMPLEX
	/* 
	 * Buffer for converting between native interleaved and vDSP format. 
	 * Size is larger of (numRows, numCols). Vector-aligned. 
	 */
	size_t numElts = max(_numRows, _numCols);
	if(fftAllocDSPComplexAlign(&mVdspBuf, numElts, VECTOR_SIZE, &mVdspBufFree)) {
		FFTError::throwMe(MR_Memory);
	}
	RFASSERT(FFT_IS_ALIGNED(mVdspBuf.realp, VECTOR_SIZE));
	RFASSERT(FFT_IS_ALIGNED(mVdspBuf.imagp, VECTOR_SIZE));
	#endif	/* FFT_SPLIT_COMPLEX */
	
}

FFTEngineVDSP::~FFTEngineVDSP()
{
	#if		!FFT_SPLIT_COMPLEX
	fftFreeDSPComplex(&mVdspBufFree);
	#endif	/* FFT_SPLIT_COMPLEX */
}
	
/* basic "FFT a contiguous set of rows" */
MFFTReturn FFTEngineVDSP::fftRows(
	bool				forward,
	FFTComplex			*inBuf,
	FFTComplex			*outBuf,
	FFTFloat			normFact,
	size_t				startRow, 
	size_t				rowsToProcess)
{
	return fftRowsAtom(forward, inBuf, outBuf, normFact, startRow, rowsToProcess, mLog2NumColsComplex);
}
	
/* 1-D forward complex twist, then FFT a contiguous set of rows */
MFFTReturn FFTEngineVDSP::twistPlusFftRow(
	FFTComplex			*inBuf,
	FFTComplex			*outBuf,
	size_t				startRow, 
	size_t				rowsToProcess)
{
	size_t numStripes  = (rowsToProcess + VDSP_DEPTH_FIRST_SIZE - 1) / VDSP_DEPTH_FIRST_SIZE;
	size_t rowsToGo = rowsToProcess;
	MFFTReturn mrtn = MR_Internal;
	size_t currRow = startRow;
	
	for(size_t stripe=0; stripe<numStripes; stripe++) {
		size_t thisRows = VDSP_DEPTH_FIRST_SIZE;
		if(thisRows > rowsToGo) {
			thisRows = rowsToGo;
		}
		mrtn = fft1DTwist(mfftPlan(), inBuf, true, numRows(), numCols(),
			currRow, thisRows);
		if(mrtn) {
			return mrtn;
		}
		mrtn = fftRowsAtom(true, inBuf, outBuf, 0.0, currRow, thisRows, mLog2NumCols);
		if(mrtn) {
			return mrtn;
		}
		
		currRow  += thisRows;
		rowsToGo -= thisRows;
	}
	return mrtn;
}
	
/* (inverse) FFT, then inverse complex twist a contiguous set of rows */
MFFTReturn FFTEngineVDSP::fftRowPlusInvTwist(
	FFTComplex			*inBuf,
	FFTComplex			*outBuf,
	size_t				startRow, 
	size_t				rowsToProcess)
{
	size_t numStripes  = (rowsToProcess + VDSP_DEPTH_FIRST_SIZE - 1) / VDSP_DEPTH_FIRST_SIZE;
	size_t rowsToGo = rowsToProcess;
	MFFTReturn mrtn = MR_Internal;
	size_t currRow = startRow;
	
	for(size_t stripe=0; stripe<numStripes; stripe++) {
		size_t thisRows = VDSP_DEPTH_FIRST_SIZE;
		if(thisRows > rowsToGo) {
			thisRows = rowsToGo;
		}
		mrtn = fftRowsAtom(false, inBuf, outBuf, 0.0, currRow, thisRows, mLog2NumCols);
		if(mrtn) {
			return mrtn;
		}
		mrtn = fft1DTwist(mfftPlan(), outBuf, false, numRows(), numCols(),
			currRow, thisRows);
		if(mrtn) {
			return mrtn;
		}
		
		currRow  += thisRows;
		rowsToGo -= thisRows;
	}
	return mrtn;
}
	
/* 
 * FFT a contiguous set of columns
 *
 * for each mColStripeSize-wide stripe {
 *		peel a column to a stripeBuf
 *		FFT the rows in stripeBuf
 *		optionally normalize
 *		unpeel from stripeBuf to outBuf
 * }
 */
MFFTReturn FFTEngineVDSP::fftColumns(
	bool				forward,
	FFTComplex			*inBuf,
	FFTComplex			*outBuf,
	FFTFloat			normFact,
	size_t				startCol, 
	size_t				colsToProcess)
{
	RFASSERT(mNumStripeBufs >= 1);
	FFTComplex *stripeBuf = mStripeBuf[0];
	RFASSERT(stripeBuf != NULL);
	size_t numStripes  = (colsToProcess + mColStripeSize - 1) / mColStripeSize;
	size_t colsToGo = colsToProcess;
	MFFTReturn mrtn = MR_Internal;
	size_t currCol = startCol;
	
	for(size_t stripe=0; stripe<numStripes; stripe++) {
		size_t thisCols = mColStripeSize;
		if(thisCols > colsToGo) {
			thisCols = colsToGo;
		}
		mrtn = fftTransposeFromColumn(inBuf, stripeBuf, numRows(), numColsComplex(), currCol, thisCols);
		if(mrtn) {
			break;
		}
		mrtn = fftRowsAtom(forward, stripeBuf, stripeBuf, normFact, 0, thisCols, mLog2NumRows);
		if(mrtn) {
			break;
		}
		mrtn = fftTransposeToColumn(stripeBuf, outBuf, numRows(), numColsComplex(), currCol, thisCols);
		if(mrtn) {
			break;
		}
		currCol += thisCols;
		colsToGo -= thisCols;
		if(colsToGo == 0) {
			break;
		}
	}
	return mrtn;
}

/* 
 * real-to-complex FFT a contiguous set of rows
 * Note: we do the optional scaling here, within the inner loops,
 * so that the fftScaleDSPComplex operates on the row while it's in 
 * cache (hopefully even in L1 cache). 
 */
#if		FFT_SPLIT_COMPLEX

MFFTReturn FFTEngineVDSP::fftRowsReal(
	bool				forward,
	FFTComplex			*inBuf,
	FFTComplex			*outBuf,
	FFTFloat			normFact,
	size_t				startRow, 
	size_t				rowsToProcess)
{
	size_t nCol = mfftPlan()->nCol;
	size_t rowSizeComplex = numColsComplex();
	size_t offsetComplex = rowSizeComplex * startRow;
	bool doScale = (normFact != 0.0);
	
	vDSPComplex vinBuf;
	vDSPComplex voutBuf;
	bool inPlace = (inBuf == outBuf) || (inBuf->real == outBuf->real);
	fftComplexToVDSP(inBuf, &vinBuf);
	if(!inPlace) {
		fftComplexToVDSP(outBuf, &voutBuf);
	}
	vinBuf.realp += offsetComplex;
	vinBuf.imagp += offsetComplex;
	if(!inPlace) {
		voutBuf.realp += offsetComplex;
		voutBuf.imagp += offsetComplex;
	}
	FFTDirection fdir = forward ? FFT_FORWARD : FFT_INVERSE;
		
	for(size_t dex=0; dex<rowsToProcess; dex++) {
		if(inPlace) {
			FFTReal1d(mFftSetup, &vinBuf, nCol, fdir);
			if(doScale) {
				fftScaleDSPComplex(&vinBuf, normFact, rowSizeComplex);
			}
		}
		else {
			FFTReal1dOP(mFftSetup, &vinBuf, &voutBuf, nCol, fdir);
			if(doScale) {
				fftScaleDSPComplex(&voutBuf, normFact, rowSizeComplex);
			}
			voutBuf.realp += rowSizeComplex;
			voutBuf.imagp += rowSizeComplex;
		}
		vinBuf.realp += rowSizeComplex;
		vinBuf.imagp += rowSizeComplex;
	}
	
	return MR_Success;

}

#else	/* !FFT_SPLIT_COMPLEX */

/* Interleaved version */

MFFTReturn FFTEngineVDSP::fftRowsReal(
	bool				forward,
	FFTComplex			*inBuf,
	FFTComplex			*outBuf,
	FFTFloat			normFact,
	size_t				startRow, 
	size_t				rowsToProcess)
{
	size_t nCol = mfftPlan()->nCol;
	size_t rowSizeComplex = numColsComplex();
	size_t offsetComplex = rowSizeComplex * startRow;
	inBuf += offsetComplex;
	outBuf += offsetComplex;
	FFTDirection fdir = forward ? FFT_FORWARD : FFT_INVERSE;
	bool doScale = (normFact != 0.0);
	
	RFASSERT(mVdspBuf.realp != NULL);
	
	for(size_t dex=0; dex<rowsToProcess; dex++) {
		fftIntToVDSP(inBuf, &mVdspBuf, rowSizeComplex);
		FFTReal1d(mFftSetup, &mVdspBuf, nCol, fdir);
		if(doScale) {
			fftScaleDSPComplex(&mVdspBuf, normFact, rowSizeComplex);
		}
		fftVDSPToInt(&mVdspBuf, outBuf, rowSizeComplex);
		inBuf  += rowSizeComplex;
		outBuf += rowSizeComplex;
	}
	
	return MR_Success;
}

#endif	/* FFT_SPLIT_COMPLEX */

MFFTReturn FFTEngineVDSP::realTwist(
	bool				forward,
	FFTComplex			*inBuf,
	FFTComplex			*outBuf,
	size_t				startRow, 
	size_t				rowsToProcess)
{
	RFASSERT(startRow <= ((numRows() >> 1) + 1));
	return fft1DRealTwist(mfftPlan(), inBuf, outBuf, forward,
		numRows(), numColsComplex(), startRow, rowsToProcess);
}

/* determine number of rows or columns need to process given op */
size_t FFTEngineVDSP::totalEltsForOp(
	FFTEngineOp			op)
{
	switch(op) {
		case FEO_FftCols:
			return numColsComplex();
		case FEO_RealTwist:
			/*
			 * We have to have a bit of detailed knowledge of the algorithm here.
			 * This twist operates on both ends of the array, with a special case
			 * for row 0 and the final row to be processed as row numRows/2. 
			 */
			return (numRows() >> 1) + 1;
		default:
			return numRows();
	}
}

/*
 * When nonzero, align FEO_FftCols to column stripe size, not submatrix size.
 */
#define VDSP_ALIGN_TO_STRIPE 0

/*
 * Determine optimal alignment for specified op. 
 */
size_t FFTEngineVDSP::alignmentForOp(
	FFTEngineOp			op)
{
	switch(op) {
		case FEO_FftRows:
		case FEO_FftRowsReal:
		case FEO_RealTwist:
			/* 
			 * These don't really need an alignment since we just go 
			 * one row at a time. However give it a nontrivial value
			 * to avoid needlessly breaking up small ops among multiple 
			 * engines (and hence threads).
			 */
			return 4;
		case FEO_TwistFft:
		case FEO_FftInvTwist:
			/* see comments where this is defined */
			return VDSP_DEPTH_FIRST_SIZE;
		case FEO_FftCols:
			/* column stripe size */
            #if     VDSP_ALIGN_TO_STRIPE
			return mColStripeSize;
            #else   /* !VDSP_ALIGN_TO_STRIPE */
			/* 
			 * Although we are capable of mColStripeSize columns at a time, 
			 * to enable full multithreading for smaller ops, we just
			 * need an alignment of one submatrix. 
			 */
			#if FFT_SPLIT_COMPLEX
			return FFT_FLOATS_PER_SUBMATRIX;
			#else
			return FFT_COMPLEX_PER_SUBMATRIX;
			#endif
            #endif  /* VDSP_ALIGN_TO_STRIPE */
		default:
			RFASSERT(0);
			printf("***FFTEngineVDSP::alignmentForOp error\n");
			return 0;
	}
}

