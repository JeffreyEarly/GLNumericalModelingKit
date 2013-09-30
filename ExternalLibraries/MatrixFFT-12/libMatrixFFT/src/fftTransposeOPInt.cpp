/*	File: fftTransposeOPInt.cpp 
	
	Description:
		Out-of-place transpose, interleaved complex version.
	
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
 * Copyright (c) 2008 Apple Computer, Inc. All Rights Reserved.
 * 
 * fftTransposeOPInt.cpp - out-of-place transpose, interleaved complex version
 *
 * Created Sep. 26 2008. 
 */

#include "MatrixFFTPlan.h"
#include "fftTranspose.h"
#include "ThreadPool.h"
#include "fftPriv.h"
#include "fftDebug.h"
#include "fftIntel.h"

#if		!FFT_SPLIT_COMPLEX

#if		FFT_INTEL

#pragma mark --- Intel precision-specific optimizations ---

/* Force 'basic' operation for debug */
#define TRANS_STRIPE_FORCE_BASIC	0

#if	FFT_DOUBLE_PREC

#define DOUBLE_STREAM_ENABLE	1

/*
 * Transpose one submatrix (FFT_COMPLEX_PER_SUBMATRIX on a side)
 * out of place. 
 * Intel, double precision, interleaved complex. 
 * This is quite trivial since an FFTComplex is the same size as an 
 * FFTVector. There is no shuffling of components between FFTVectors,
 * it's a straight memory move. 
 */

#if !TRANS_STRIPE_FORCE_BASIC

/* 
 * Normal memory store version, when the transposed data needs
 * to remain in cache, e.g. for transposing a column to a stripe buffer.
 */
static void fftOPSubTrans_store(
	const FFTComplex	*_src,
	FFTComplex			*_dst,
	size_t				srcRowSize,		// src, in FFTComplex, a.k.a. src numCols
	size_t				dstRowSize)		// dst, in FFTComplex, a.k.a. dst numCols
{
	/* rowDex, colDex refer to _src */
	for(size_t rowDex=0; rowDex<FFT_COMPLEX_PER_SUBMATRIX; rowDex++) {
		const FFTVector *invp = (const FFTVector *)(_src + (srcRowSize * rowDex));
		FFTVector *outvp = (FFTVector *)(_dst + rowDex);
		for(size_t colDex=0; colDex<FFT_COMPLEX_PER_SUBMATRIX; colDex++) {
			register FFTVector tmp = *invp++;
			*outvp = tmp;
			outvp += dstRowSize;
		}
	}
}

#endif  /* TRANS_STRIPE_FORCE_BASIC */

#if		DOUBLE_STREAM_ENABLE
/* _mm_stream_pd version, used for transposing from a stripe buffer to columns. */
static void fftOPSubTrans(
	  const FFTComplex	*_src,
	  FFTComplex		*_dst,
	  size_t			srcRowSize,		// src, in FFTComplex, a.k.a. src numCols
	  size_t			dstRowSize)		// dst, in FFTComplex, a.k.a. dst numCols
{
	/* rowDex, colDex refer to _src */
	for(size_t colDex=0; colDex<FFT_COMPLEX_PER_SUBMATRIX; colDex++) {

		const FFTVector *invp = (const FFTVector *)(_src + colDex);
		FFTVector		*outvp = (FFTVector *)(_dst + colDex*dstRowSize);
		
		for(size_t rowDex=0; rowDex<FFT_COMPLEX_PER_SUBMATRIX; rowDex++) {

			register FFTVector tmp = *invp;
			
			_mm_stream_pd((double*)outvp, tmp);
		
			outvp += 1;
			invp += srcRowSize;
		}
	}
}

#else
#define fftOPSubTrans fftOPSubTrans_store
#endif

#else	/* !FFT_DOUBLE_PREC */

#if     !TRANS_STRIPE_FORCE_BASIC

/* 
 * Single precision, _mm_store_pd version. Used when the transposed data needs
 * to remain in cache, e.g. for transposing a column to a stripe buffer. 
 *
 * Transpose one submatrix (FFT_COMPLEX_PER_SUBMATRIX on a side)
 * out of place. 
 * Based on Jason Klivington's implementation in the RowFFT project.  
 * This is a bit of a hack which assumes/knows that the size of 
 * a FFTComplex is the same as the size of a double. Since we're
 * just moving data, not computing anything with it, this is OK...
 */
static void fftOPSubTrans_store(
	const FFTComplex	*_src,
	FFTComplex			*_dst,
	size_t				srcRowSize,		// src, in FFTComplex, a.k.a. src numCols
	size_t				dstRowSize)		// dst, in FFTComplex, a.k.a. dst numCols
{
	double *src = (double *)_src;
	double *dst = (double *)_dst;
	
	dumpSub("fftOPSubTrans_store start", _src, srcRowSize);

	/* 
	 * row and col refer to coordinates in src 
	 * row size of dst is dstRowSize
	 */
	unsigned curcol;
	
	for(curcol=0; curcol<FFT_COMPLEX_PER_SUBMATRIX; curcol+=2) {
		__m128d vin1;
		__m128d vin2;
		__m128d vin3;
		__m128d vin4;
		__m128d vin5;
		__m128d vin6;
		__m128d vin7;
		__m128d vin8;

		__m128d vOut_row1_1;
		__m128d vOut_row1_2;
		__m128d vOut_row1_3;
		__m128d vOut_row1_4;
		__m128d vOut_row2_1;
		__m128d vOut_row2_2;
		__m128d vOut_row2_3;
		__m128d vOut_row2_4;

		const double *pIn = src + curcol;
		double *pOut = dst + curcol*dstRowSize;
		
		// load in two columns from src at curcol
		vin1 = _mm_load_pd(pIn+0*srcRowSize);
		vin2 = _mm_load_pd(pIn+1*srcRowSize);
		vin3 = _mm_load_pd(pIn+2*srcRowSize);
		vin4 = _mm_load_pd(pIn+3*srcRowSize);
		vin5 = _mm_load_pd(pIn+4*srcRowSize);
		vin6 = _mm_load_pd(pIn+5*srcRowSize);
		vin7 = _mm_load_pd(pIn+6*srcRowSize);
		vin8 = _mm_load_pd(pIn+7*srcRowSize);

		///////////////////////////////////////////////
		// transpose for first row out
		
		vOut_row1_1 = _mm_unpacklo_pd(vin1, vin2);
		vOut_row1_2 = _mm_unpacklo_pd(vin3, vin4);
		vOut_row1_3 = _mm_unpacklo_pd(vin5, vin6);
		vOut_row1_4 = _mm_unpacklo_pd(vin7, vin8);
		
		_mm_store_pd(pOut+(0*FFT_COMPLEX_PER_VECTOR), vOut_row1_1);
		_mm_store_pd(pOut+(1*FFT_COMPLEX_PER_VECTOR), vOut_row1_2);
		_mm_store_pd(pOut+(2*FFT_COMPLEX_PER_VECTOR), vOut_row1_3);
		_mm_store_pd(pOut+(3*FFT_COMPLEX_PER_VECTOR), vOut_row1_4);

		///////////////////////////////////////////////
		// transpose for second row out
		pOut += dstRowSize;
		
		vOut_row2_1 = _mm_unpackhi_pd(vin1, vin2);
		vOut_row2_2 = _mm_unpackhi_pd(vin3, vin4);
		vOut_row2_3 = _mm_unpackhi_pd(vin5, vin6);
		vOut_row2_4 = _mm_unpackhi_pd(vin7, vin8);
		
		_mm_store_pd(pOut+(0*FFT_COMPLEX_PER_VECTOR), vOut_row2_1);
		_mm_store_pd(pOut+(1*FFT_COMPLEX_PER_VECTOR), vOut_row2_2);
		_mm_store_pd(pOut+(2*FFT_COMPLEX_PER_VECTOR), vOut_row2_3);
		_mm_store_pd(pOut+(3*FFT_COMPLEX_PER_VECTOR), vOut_row2_4);
	}
	
	dumpSub("fftOPSubTrans_store end", _dst, dstRowSize);
}

#endif  /* TRANS_STRIPE_FORCE_BASIC */

/* 
 * Intel single precision, _mm_stream_pd version, used for transposing
 * from a stripe buffer to columns. 
 */
static void fftOPSubTrans(
  const FFTComplex	*_src,
  FFTComplex		*_dst,
  size_t			srcRowSize,		// src, in FFTComplex, a.k.a. src numCols
  size_t			dstRowSize)		// dst, in FFTComplex, a.k.a. dst numCols
{
	double *src = (double *)_src;
	double *dst = (double *)_dst;
	
	dumpSub("fftOPSubTrans start", _src, srcRowSize);
	
	/* 
	 * row and col refer to coordinates in src 
	 * row size of dst is dstRowSize
	 */
	unsigned curcol;
	
	for(curcol=0; curcol<FFT_COMPLEX_PER_SUBMATRIX; curcol+=2) {
		__m128d vin1;
		__m128d vin2;
		__m128d vin3;
		__m128d vin4;
		__m128d vin5;
		__m128d vin6;
		__m128d vin7;
		__m128d vin8;
		
		__m128d vOut_row1_1;
		__m128d vOut_row1_2;
		__m128d vOut_row1_3;
		__m128d vOut_row1_4;
		__m128d vOut_row2_1;
		__m128d vOut_row2_2;
		__m128d vOut_row2_3;
		__m128d vOut_row2_4;
		
		const double *pIn = src + curcol;
		double *pOut = dst + curcol*dstRowSize;
		
		// load in two columns from src at curcol
		vin1 = _mm_load_pd(pIn+0*srcRowSize);
		vin2 = _mm_load_pd(pIn+1*srcRowSize);
		vin3 = _mm_load_pd(pIn+2*srcRowSize);
		vin4 = _mm_load_pd(pIn+3*srcRowSize);
		vin5 = _mm_load_pd(pIn+4*srcRowSize);
		vin6 = _mm_load_pd(pIn+5*srcRowSize);
		vin7 = _mm_load_pd(pIn+6*srcRowSize);
		vin8 = _mm_load_pd(pIn+7*srcRowSize);
		
		///////////////////////////////////////////////
		// transpose for first row out
		
		vOut_row1_1 = _mm_unpacklo_pd(vin1, vin2);
		vOut_row1_2 = _mm_unpacklo_pd(vin3, vin4);
		vOut_row1_3 = _mm_unpacklo_pd(vin5, vin6);
		vOut_row1_4 = _mm_unpacklo_pd(vin7, vin8);
		
		_mm_stream_pd(pOut+(0*FFT_COMPLEX_PER_VECTOR), vOut_row1_1);
		_mm_stream_pd(pOut+(1*FFT_COMPLEX_PER_VECTOR), vOut_row1_2);
		_mm_stream_pd(pOut+(2*FFT_COMPLEX_PER_VECTOR), vOut_row1_3);
		_mm_stream_pd(pOut+(3*FFT_COMPLEX_PER_VECTOR), vOut_row1_4);
		
		///////////////////////////////////////////////
		// transpose for second row out
		pOut += dstRowSize;
		
		vOut_row2_1 = _mm_unpackhi_pd(vin1, vin2);
		vOut_row2_2 = _mm_unpackhi_pd(vin3, vin4);
		vOut_row2_3 = _mm_unpackhi_pd(vin5, vin6);
		vOut_row2_4 = _mm_unpackhi_pd(vin7, vin8);
		
		_mm_stream_pd(pOut+(0*FFT_COMPLEX_PER_VECTOR), vOut_row2_1);
		_mm_stream_pd(pOut+(1*FFT_COMPLEX_PER_VECTOR), vOut_row2_2);
		_mm_stream_pd(pOut+(2*FFT_COMPLEX_PER_VECTOR), vOut_row2_3);
		_mm_stream_pd(pOut+(3*FFT_COMPLEX_PER_VECTOR), vOut_row2_4);
	}
	
	dumpSub("fftOPSubTrans end", _dst, dstRowSize);
}



#endif	/* !FFT_DOUBLE_PREC */

#pragma mark --- Intel precision-independent optimizations ---

/* 
 * Optimized out-of-place transpose using submatrices. The size of the 
 * matrix must be >= FFT_COMPLEX_PER_SUBMATRIX in each dimension and 
 * also be a multiple of FFT_COMPLEX_PER_SUBMATRIX in each dimension.
 */
static MFFTReturn fftOPSwapRowsColsOpt(
	const FFTComplex	*src,
	FFTComplex			*dst,
	/* total size of src */
	size_t				numRows,
	size_t				numCols,
	/* portion to move - must be FFT_COMPLEX_PER_SUBMATRIX aligned */
	size_t				startRow,
	size_t				rowsToMove)
{
	size_t numRowSubs = rowsToMove / FFT_COMPLEX_PER_SUBMATRIX;
	size_t numColSubs = numCols / FFT_COMPLEX_PER_SUBMATRIX;
	size_t row = startRow;
	unsigned colSubIncr = FFT_COMPLEX_PER_SUBMATRIX * numRows;

	tpThreadDebug("fftOPSwapRowsColsOpt top: startRow %lu rowsToMove %lu\n", 
		(unsigned long)startRow, (unsigned long)rowsToMove);

	RFASSERT((numRows % FFT_COMPLEX_PER_SUBMATRIX) == 0);
	RFASSERT((numCols % FFT_COMPLEX_PER_SUBMATRIX) == 0);
	RFASSERT((startRow % FFT_COMPLEX_PER_SUBMATRIX) == 0);
	RFASSERT((rowsToMove % FFT_COMPLEX_PER_SUBMATRIX) == 0);
	
	for(size_t rowSubDex=0; rowSubDex<numRowSubs; rowSubDex++) {
	
		const FFTComplex *srcSub = src + (row * numCols);
		FFTComplex *dstSub = dst + row;
		
		for(size_t colSubDex=0; colSubDex<numColSubs; colSubDex++) {
			fftOPSubTrans(srcSub, dstSub, numCols, numRows);
			srcSub += FFT_COMPLEX_PER_SUBMATRIX;
			dstSub += colSubIncr;
		}
		row += FFT_COMPLEX_PER_SUBMATRIX;
	}
	return MR_Success;
}

/* Called out from thread module */
static MFFTReturn transOPThr(TP_TaskUnion *u)
{
	TP_TransposeInt *ttop  = &u->transI;
	fftOPSwapRowsColsOpt(ttop->src, ttop->dst, 
		ttop->numRows, ttop->numCols,
		ttop->startRow, ttop->rowsToMove);
	return MR_Success;
}

#pragma mark --- Threaded out-of-place transpose ---

/* 
 * Threaded out-of-place transpose using submatrices. The size of the 
 * matrix must be >= FFT_COMPLEX_PER_SUBMATRIX in each dimension and 
 * also be a multiple of FFT_COMPLEX_PER_SUBMATRIX in each dimension.
 */
static MFFTReturn fftOPSwapRowsColsThr(
	MatrixFFTPlan		mfftPlan,
	const FFTComplex	*src,
	FFTComplex			*dst,
	size_t				numRows,		// rows in src
	size_t				numCols)		// columns in src
{
	tpThreadDebug("fftOPSwapRowsColsThr top: numRows %lu numCols %lu src %p\n",
		(unsigned long)numRows, (unsigned long)numCols, src);

	MFFTReturn ourRtn = MR_Success;
	size_t numRowSubs = numRows / FFT_COMPLEX_PER_SUBMATRIX;
	RFASSERT(mfftPlan->threadPoolP->numThreads < numRowSubs);

	unsigned numThreads = mfftPlan->threadPoolP->numThreads;
	if(numThreads > numRowSubs) {
		numThreads = numRowSubs;
	}
	
	/* submatrices per thread, we (main) get remainder */
	size_t subsPerThread = numRowSubs / numThreads;
	
	/* rows per thread (ours might be bigger) */
	size_t rowsPerThread = subsPerThread * FFT_COMPLEX_PER_SUBMATRIX;

	/* set up a TP_TransposeInt for each thread */	
	for(unsigned dex=0; dex<numThreads; dex++) {
		TP_PerThread *pt   = &mfftPlan->threadPoolP->perThread[dex];
		TP_Task *task      = &pt->task;
		TP_TransposeInt *ttop = &task->u->transI;

		task->op        = TPO_App;
		task->threadFcn  = transOPThr;
		ttop->src        = src;
		ttop->dst        = dst;
		ttop->numRows    = numRows;
		ttop->numCols    = numCols;
		ttop->startRow   = dex * subsPerThread * FFT_COMPLEX_PER_SUBMATRIX;
		ttop->rowsToMove = rowsPerThread;
	}
	
	/* last one might have more rows */
	TP_TransposeInt *ttop = &mfftPlan->threadPoolP->perThread[numThreads-1].task.u->transI;
	ttop->rowsToMove = numRows - ((numThreads - 1) * rowsPerThread);

	/* GO */
	ourRtn = tpThreadDispatch(mfftPlan->threadPoolP, numThreads);
	if(ourRtn) {
		return ourRtn;
	}
	ourRtn = tpThreadFinish(mfftPlan->threadPoolP, numThreads);
	tpThreadDebug("fftOPSwapRowsColsThr end; status %d\n", ourRtn);
	return ourRtn;
}


#endif	/* FFT_INTEL */

/* 
 * Brute force OOP transpose, for small or non-Intel ops
 */
static MFFTReturn fftTransposeOPSmall(
	const FFTComplex	*srcBuf,
	FFTComplex			*dstBuf,
	size_t				srcRows,
	size_t				srcCols)
{
	for(size_t row=0; row<srcRows; row++) {
		/*
		 * ip : source : moves across a row in src
		 * op : dest   : moves down a column in dst 
		 */
		const FFTComplex *ip = srcBuf + (row * srcCols);	
		FFTComplex *op = dstBuf + row;
		for(size_t col=0; col<srcCols; col++) {
			*op = *ip++;
			op += srcRows;
		}
	}
	return MR_Success;
}

#if (FFT_INTEL && !TRANS_STRIPE_FORCE_BASIC) 

#pragma mark --- Intel-optimized columnar stripe transpose routines ---

/*
 * Transpose a columnar stripe to a stripe buffer
 */
static MFFTReturn fftTransposeFromColumnOpt(
	const FFTComplex	*src,		/* numRows x numCols */
	FFTComplex			*dst,		/* colsToMove rows x numRows cols */
	/* total size of src */
	size_t				numRows,
	size_t				numCols,
	/* portion to move - in src - must be FFT_COMPLEX_PER_SUBMATRIX aligned */
	size_t				startCol,
	size_t				colsToMove)
{
	size_t numRowSubs = numRows    / FFT_COMPLEX_PER_SUBMATRIX;
	size_t numColSubs = colsToMove / FFT_COMPLEX_PER_SUBMATRIX;
	size_t srcRow = 0;
	unsigned colSubIncr = FFT_COMPLEX_PER_SUBMATRIX * numRows;

	tpThreadDebug("fftOPTransposeFromColumn top: startCol %lu colsToMove %lu\n", 
		(unsigned long)startCol, (unsigned long)colsToMove);

	RFASSERT((numRows % FFT_COMPLEX_PER_SUBMATRIX) == 0);
	RFASSERT((numCols % FFT_COMPLEX_PER_SUBMATRIX) == 0);
	RFASSERT((startCol % FFT_COMPLEX_PER_SUBMATRIX) == 0);
	RFASSERT((colsToMove % FFT_COMPLEX_PER_SUBMATRIX) == 0);
	
	for(size_t rowSubDex=0; rowSubDex<numRowSubs; rowSubDex++) {
	
		const FFTComplex *srcSub = src + (srcRow * numCols) + startCol;
		FFTComplex       *dstSub = dst + srcRow;
		
		for(size_t colSubDex=0; colSubDex<numColSubs; colSubDex++) {
			fftOPSubTrans_store(srcSub, dstSub, numCols, numRows);
			srcSub += FFT_COMPLEX_PER_SUBMATRIX;	// to right, one submatrix
			dstSub += colSubIncr;				// down, one submatrix
		}
		srcRow += FFT_COMPLEX_PER_SUBMATRIX;
	}
	return MR_Success;
}

/*
 * Transpose a stripe buffer to a columnar stripe
 */
static MFFTReturn fftTransposeToColumnOpt(
	const FFTComplex	*src,		/* colsToMove rows x numRows cols  */
	FFTComplex			*dst,		/* numRows x numCols */
	/* total size of dst */
	size_t				numRows,
	size_t				numCols,
	/* portion to move - in dst - must be FFT_COMPLEX_PER_SUBMATRIX aligned */
	size_t				startCol,
	size_t				colsToMove)
{
	size_t numRowSubs = numRows    / FFT_COMPLEX_PER_SUBMATRIX;
	size_t numColSubs = colsToMove / FFT_COMPLEX_PER_SUBMATRIX;
	size_t dstRow = 0;
	unsigned colSubIncr = FFT_COMPLEX_PER_SUBMATRIX * numRows;

	tpThreadDebug("fftOPTransposeFromColumn top: startCol %lu colsToMove %lu\n", 
		(unsigned long)startCol, (unsigned long)colsToMove);

	RFASSERT((numRows % FFT_COMPLEX_PER_SUBMATRIX) == 0);
	RFASSERT((numCols % FFT_COMPLEX_PER_SUBMATRIX) == 0);
	RFASSERT((startCol % FFT_COMPLEX_PER_SUBMATRIX) == 0);
	RFASSERT((colsToMove % FFT_COMPLEX_PER_SUBMATRIX) == 0);
	
	for(size_t rowSubDex=0; rowSubDex<numRowSubs; rowSubDex++) {
	
		FFTComplex		 *dstSub = dst + (dstRow * numCols) + startCol;
		const FFTComplex *srcSub = src + (rowSubDex * FFT_COMPLEX_PER_SUBMATRIX);
		
		for(size_t colSubDex=0; colSubDex<numColSubs; colSubDex++) {
			fftOPSubTrans(srcSub, dstSub, numRows, numCols);
			dstSub += FFT_COMPLEX_PER_SUBMATRIX;	// to right, one submatrix
			srcSub += colSubIncr;				// down, one submatrix
		}
		dstRow += FFT_COMPLEX_PER_SUBMATRIX;
	}
	return MR_Success;
}

#endif	/* FFT_INTEL */

/*
 * Brute force implementation of "transpose from columnar stripe"
 */
static MFFTReturn fftTransposeFromColumnBasic(
	const FFTComplex	*src,		/* numRows x numCols */
	FFTComplex			*dst,		/* colsToMove rows x numRows cols */
	/* total size of src */
	size_t				numRows,
	size_t				numCols,
	/* portion to move - in src - must be FFT_COMPLEX_PER_SUBMATRIX aligned */
	size_t				startCol,
	size_t				colsToMove)
{
	for(size_t srcRow=0; srcRow<numRows; srcRow++) {
		const FFTComplex *inp  = src + (srcRow * numCols) + startCol;
		FFTComplex       *outp = dst + srcRow;
		for(size_t colDex=0; colDex<colsToMove; colDex++) {
			*outp = *inp++;
			outp += numRows;
		}
	}
	return MR_Success;
}

/*
 * Brute force implementation of "transpose to columnar stripe"
 */
static MFFTReturn fftTransposeToColumnBasic(
	const FFTComplex	*src,		/* colsToMove rows x numRows cols  */
	FFTComplex			*dst,		/* numRows x numCols */
	/* total size of dst */
	size_t				numRows,
	size_t				numCols,
	/* portion to move - in dst - must be FFT_COMPLEX_PER_SUBMATRIX aligned */
	size_t				startCol,
	size_t				colsToMove)
{
	for(size_t dstRow=0; dstRow<numRows; dstRow++) {
		FFTComplex       *outp = dst + (dstRow * numCols) + startCol;
		const FFTComplex *inp  = src + dstRow;
		for(size_t colDex=0; colDex<colsToMove; colDex++) {
			*outp++ = *inp;
			inp += numRows;
		}
	}
	return MR_Success;
}

#pragma mark --- Public API ---

MFFTReturn fftTransposeOP(
	MatrixFFTPlan		mfftPlan,
	const FFTComplex	*srcBuf,
	FFTComplex			*dstBuf,
	size_t				srcRows,
	size_t				srcCols)
{
	#if	!FFT_INTEL
	/* For now, use brute force */
	return fftTransposeOPSmall(srcBuf, dstBuf, srcRows, srcCols);
	#else
	
	if((srcRows >= FFT_COMPLEX_PER_SUBMATRIX) &&
	   (srcCols >= FFT_COMPLEX_PER_SUBMATRIX)) {
		unsigned numThreads = mfftPlan->threadPoolP->numThreads;		/* # of pthreads */
		if((numThreads > 1) &&											/* threading enabled */
		   ((srcRows / FFT_COMPLEX_PER_SUBMATRIX) > numThreads))	{	/* more subs then threads */
		   return fftOPSwapRowsColsThr(mfftPlan, srcBuf, dstBuf, srcRows, srcCols);
		}
		else {
			return fftOPSwapRowsColsOpt(srcBuf, dstBuf, srcRows, srcCols, 0, srcRows);
		}
	}
	else {
		return fftTransposeOPSmall(srcBuf, dstBuf, srcRows, srcCols);
	}
	#endif	/* !FFT_INTEL */
	
}

MFFTReturn fftTransposeFromColumn(
	const FFTComplex	*src,		/* numRows x numCols */
	FFTComplex			*dst,		/* colsToMove rows x numRows cols */
	/* total size of src */
	size_t				numRows,
	size_t				numCols,
	/* portion to move - in src */
	size_t				startCol,
	size_t				colsToMove)
{
	#if		(FFT_INTEL && !TRANS_STRIPE_FORCE_BASIC)
	if((numRows >= FFT_COMPLEX_PER_SUBMATRIX) &&
	   (numCols >= FFT_COMPLEX_PER_SUBMATRIX) &&
	   FFT_IS_ALIGNED(colsToMove, FFT_COMPLEX_PER_SUBMATRIX) &&
	   FFT_IS_ALIGNED(startCol, FFT_COMPLEX_PER_SUBMATRIX)) {
		return fftTransposeFromColumnOpt(src, dst, numRows, numCols, startCol, colsToMove);
	}
	/* else drop thru and use the brute force implementation */
	#endif	/* FFT_INTEL) */

	return fftTransposeFromColumnBasic(src, dst, numRows, numCols, startCol, colsToMove);
}

MFFTReturn fftTransposeToColumn(
	const FFTComplex	*src,		/* colsToMove rows x numRows cols  */
	FFTComplex			*dst,		/* numRows x numCols */
	/* total size of dst */
	size_t				numRows,
	size_t				numCols,
	/* portion to move - in dst */
	size_t				startCol,
	size_t				colsToMove)
{
	#if		(FFT_INTEL && !TRANS_STRIPE_FORCE_BASIC)
	if((numRows >= FFT_COMPLEX_PER_SUBMATRIX) &&
	   (numCols >= FFT_COMPLEX_PER_SUBMATRIX) &&
	   FFT_IS_ALIGNED(colsToMove, FFT_COMPLEX_PER_SUBMATRIX) &&
	   FFT_IS_ALIGNED(startCol, FFT_COMPLEX_PER_SUBMATRIX)) {
		return fftTransposeToColumnOpt(src, dst, numRows, numCols, startCol, colsToMove);
	}
	/* else drop thru and use the brute force implementation */
	#endif	/* !FFT_INTEL */

	return fftTransposeToColumnBasic(src, dst, numRows, numCols, startCol, colsToMove);
}

#endif	/* !FFT_SPLIT_COMPLEX */
