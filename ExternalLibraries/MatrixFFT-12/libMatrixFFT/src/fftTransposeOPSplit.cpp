/*	File: fftTransposeOPSplit.cpp  
	
	Description:
		Out-of-place transpose, Split complex version.
	
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
 * fftTransposeOPSplit.cpp - out-of-place transpose, Split complex version. 
 *
 * Created Sep. 26 2008. 
 */

#include "MatrixFFTPlan.h"
#include "fftTranspose.h"
#include "ThreadPool.h"
#include "fftPriv.h"
#include "fftDebug.h"
#include "fftIntel.h"
#include "fftPlatformConf.h"

#include <libMatrixFFT/fftUtils.h>

/* This entire file is #if FFT_SPLIT_COMPLEX */
#if		FFT_SPLIT_COMPLEX

#define FULL_VECTOR_TRANSPOSE	1
#define VECTOR_TRANSPOSE_ENABLE	(FULL_VECTOR_TRANSPOSE && FFT_INTEL)

/* Possible platform-specific vector disabled; currently unused. */
#define TRANS_STRIPE_DISABLE_VECTOR     0

/* Force 'basic' operation for debug */
#define TRANS_STRIPE_BASIC_DEBUG    0

#define TRANS_STRIPE_FORCE_BASIC	(TRANS_STRIPE_DISABLE_VECTOR || TRANS_STRIPE_BASIC_DEBUG)


#if		VECTOR_TRANSPOSE_ENABLE

#if		FFT_DOUBLE_PREC

/* Intel Double Precision */

/* 
 * Transpose one submatrix (FFT_FLOATS_PER_SUBMATRIX on a side)
 * out of place, _mm_stream_pd version. 
 * Intel; Double Precision; Split Complex only. 
 * Used when transposing from a stripe buffer to columns; the store to 
 * the actual column buffer bypasses the cache since in this op we're
 * done with the data.
 */
static void  fftOPSubTrans(
	const FFTFloat	*src,
	FFTFloat		*dst,
	size_t			rowSize,		// src, in FFTFloat, a.k.a. numCols
	size_t			colSize)		// src, in FFTFloat, a.k.a. numRows
{
	dumpSub("fftOPSubTrans start", src, rowSize);

	/* 
	 * row and col refer to coordinates in src 
	 * row size of dst is colSize
	 */
	unsigned curcol;
	
	for(curcol=0; curcol<FFT_FLOATS_PER_SUBMATRIX; curcol+=2) {
		FFTVector vin1;
		FFTVector vin2;
		FFTVector vin3;
		FFTVector vin4;
		FFTVector vin5;
		FFTVector vin6;
		FFTVector vin7;
		FFTVector vin8;

		FFTVector vOut_row1_1;
		FFTVector vOut_row1_2;
		FFTVector vOut_row1_3;
		FFTVector vOut_row1_4;
		FFTVector vOut_row2_1;
		FFTVector vOut_row2_2;
		FFTVector vOut_row2_3;
		FFTVector vOut_row2_4;

		const double *pIn = src + curcol;
		double *pOut = dst + curcol*colSize;
		
		// load in two columns from src at curcol
		vin1 = _mm_load_pd(pIn+0*rowSize);
		vin2 = _mm_load_pd(pIn+1*rowSize);
		vin3 = _mm_load_pd(pIn+2*rowSize);
		vin4 = _mm_load_pd(pIn+3*rowSize);
		vin5 = _mm_load_pd(pIn+4*rowSize);
		vin6 = _mm_load_pd(pIn+5*rowSize);
		vin7 = _mm_load_pd(pIn+6*rowSize);
		vin8 = _mm_load_pd(pIn+7*rowSize);

		///////////////////////////////////////////////
		// transpose for first row out
		
		vOut_row1_1 = _mm_unpacklo_pd(vin1, vin2);
		vOut_row1_2 = _mm_unpacklo_pd(vin3, vin4);
		vOut_row1_3 = _mm_unpacklo_pd(vin5, vin6);
		vOut_row1_4 = _mm_unpacklo_pd(vin7, vin8);
		
		_mm_stream_pd(pOut+(0*FFT_FLOATS_PER_VECTOR), vOut_row1_1);
		_mm_stream_pd(pOut+(1*FFT_FLOATS_PER_VECTOR), vOut_row1_2);
		_mm_stream_pd(pOut+(2*FFT_FLOATS_PER_VECTOR), vOut_row1_3);
		_mm_stream_pd(pOut+(3*FFT_FLOATS_PER_VECTOR), vOut_row1_4);

		///////////////////////////////////////////////
		// transpose for second row out
		pOut += colSize;
		
		vOut_row2_1 = _mm_unpackhi_pd(vin1, vin2);
		vOut_row2_2 = _mm_unpackhi_pd(vin3, vin4);
		vOut_row2_3 = _mm_unpackhi_pd(vin5, vin6);
		vOut_row2_4 = _mm_unpackhi_pd(vin7, vin8);
		
		_mm_stream_pd(pOut+(0*FFT_FLOATS_PER_VECTOR), vOut_row2_1);
		_mm_stream_pd(pOut+(1*FFT_FLOATS_PER_VECTOR), vOut_row2_2);
		_mm_stream_pd(pOut+(2*FFT_FLOATS_PER_VECTOR), vOut_row2_3);
		_mm_stream_pd(pOut+(3*FFT_FLOATS_PER_VECTOR), vOut_row2_4);
	}
	
	dumpSub("fftOPSubTrans end", dst, colSize);
}

/* 
 * Transpose one submatrix (FFT_FLOATS_PER_SUBMATRIX on a side)
 * out of place, straight memory write version. 
 * Intel; Double Precision; Split Complex only. 
 * Used when transposing from a column to a stripe buffer; data remains in cache.
 * the actual column buffer bypasses the cache since in this op we're
 * done with the data.
 */
#define writeVect(p, v)     (*(FFTVector *)(p)) = v

static void  fftOPSubTrans_store(
	const FFTFloat	*src,
	FFTFloat		*dst,
	size_t			rowSize,		// src, in FFTFloat, a.k.a. numCols
	size_t			colSize)		// src, in FFTFloat, a.k.a. numRows
{
	dumpSub("fftOPSubTrans_store start", src, rowSize);

	/* 
	 * row and col refer to coordinates in src 
	 * row size of dst is colSize
	 */
	unsigned curcol;
	
	for(curcol=0; curcol<FFT_FLOATS_PER_SUBMATRIX; curcol+=2) {
		FFTVector vin1;
		FFTVector vin2;
		FFTVector vin3;
		FFTVector vin4;
		FFTVector vin5;
		FFTVector vin6;
		FFTVector vin7;
		FFTVector vin8;

		FFTVector vOut_row1_1;
		FFTVector vOut_row1_2;
		FFTVector vOut_row1_3;
		FFTVector vOut_row1_4;
		FFTVector vOut_row2_1;
		FFTVector vOut_row2_2;
		FFTVector vOut_row2_3;
		FFTVector vOut_row2_4;

		const double *pIn = src + curcol;
		double *pOut = dst + curcol*colSize;
		
		// load in two columns from src at curcol
		vin1 = _mm_load_pd(pIn+0*rowSize);
		vin2 = _mm_load_pd(pIn+1*rowSize);
		vin3 = _mm_load_pd(pIn+2*rowSize);
		vin4 = _mm_load_pd(pIn+3*rowSize);
		vin5 = _mm_load_pd(pIn+4*rowSize);
		vin6 = _mm_load_pd(pIn+5*rowSize);
		vin7 = _mm_load_pd(pIn+6*rowSize);
		vin8 = _mm_load_pd(pIn+7*rowSize);

		///////////////////////////////////////////////
		// transpose for first row out
		
		vOut_row1_1 = _mm_unpacklo_pd(vin1, vin2);
		vOut_row1_2 = _mm_unpacklo_pd(vin3, vin4);
		vOut_row1_3 = _mm_unpacklo_pd(vin5, vin6);
		vOut_row1_4 = _mm_unpacklo_pd(vin7, vin8);
		
		writeVect(pOut+(0*FFT_FLOATS_PER_VECTOR), vOut_row1_1);
		writeVect(pOut+(1*FFT_FLOATS_PER_VECTOR), vOut_row1_2);
		writeVect(pOut+(2*FFT_FLOATS_PER_VECTOR), vOut_row1_3);
		writeVect(pOut+(3*FFT_FLOATS_PER_VECTOR), vOut_row1_4);

		///////////////////////////////////////////////
		// transpose for second row out
		pOut += colSize;
		
		vOut_row2_1 = _mm_unpackhi_pd(vin1, vin2);
		vOut_row2_2 = _mm_unpackhi_pd(vin3, vin4);
		vOut_row2_3 = _mm_unpackhi_pd(vin5, vin6);
		vOut_row2_4 = _mm_unpackhi_pd(vin7, vin8);
		
		writeVect(pOut+(0*FFT_FLOATS_PER_VECTOR), vOut_row2_1);
		writeVect(pOut+(1*FFT_FLOATS_PER_VECTOR), vOut_row2_2);
		writeVect(pOut+(2*FFT_FLOATS_PER_VECTOR), vOut_row2_3);
		writeVect(pOut+(3*FFT_FLOATS_PER_VECTOR), vOut_row2_4);
	}
	
	dumpSub("fftOPSubTrans_store end", dst, colSize);
}

#else	/* !FFT_DOUBLE_PREC */

/* Intel Single Precision */

#define TRANSPOSE_4x4_VEC(in, out) \
do {\
	vFloat	temp1, temp2, temp3, temp4;\
	\
	temp1  = _mm_unpacklo_ps( in##1, in##3 );\
	temp2  = _mm_unpacklo_ps( in##2, in##4 );\
	temp3  = _mm_unpackhi_ps( in##1, in##3 );\
	temp4  = _mm_unpackhi_ps( in##2, in##4 );\
	out##1 = _mm_unpacklo_ps( temp1, temp2 );\
	out##2 = _mm_unpackhi_ps( temp1, temp2 );\
	out##3 = _mm_unpacklo_ps( temp3, temp4 );\
	out##4 = _mm_unpackhi_ps( temp3, temp4 );\
} while (0)

#define WRITE_VECT_NORM		1

#if		WRITE_VECT_NORM
#define writeVect(p, v)		*((FFTVector *)p) = v
#else
#define writeVect(p, v)		_mm_stream_ps(p, v)
#endif  /* WRITE_VECT_NORM */

/* 
 * Transpose one submatrix (FFT_FLOATS_PER_SUBMATRIX on a side)
 * out of place. 
 * Intel; Single Precision; Split Complex only. 
 */
static void fftOPSubTrans(
	const FFTFloat		*src,
	FFTFloat			*dst,
	size_t				srcRowSize,		// src, in FFTFloats, a.k.a. src numCols
	size_t				dstRowSize)		// dst, in FFTFloats, a.k.a. dst numCols
{
	/*
	 * Each iteration of the outer loop moves 4 columns from src to 4 rows of
	 * dst.
	 */
	for(size_t col=0; col<FFT_FLOATS_PER_SUBMATRIX; col+=FFT_FLOATS_PER_VECTOR) {
		const FFTFloat *inP = src + col;
		FFTFloat *outRowStart = dst + (col * dstRowSize);
		FFTVector vin1, vin2, vin3, vin4;
		FFTVector vout1, vout2, vout3, vout4;

		for(size_t row=0; row<FFT_FLOATS_PER_SUBMATRIX; row+=FFT_FLOATS_PER_VECTOR) {
			/* read 4 rows, 4 columns wide */
			vin1 = _mm_load_ps(inP);
			inP += srcRowSize;
			vin2 = _mm_load_ps(inP);
			inP += srcRowSize;
			vin3 = _mm_load_ps(inP);
			inP += srcRowSize;
			vin4 = _mm_load_ps(inP);
			inP += srcRowSize;
			
			/* Transpose */
			TRANSPOSE_4x4_VEC(vin, vout);
			
			/* write 4 columns of 4 rows */
			FFTFloat *outP = outRowStart + row;
			writeVect(outP, vout1);
			outP += dstRowSize;
			writeVect(outP, vout2);
			outP += dstRowSize;
			writeVect(outP, vout3);
			outP += dstRowSize;
			writeVect(outP, vout4);
		}
	}
}

#define fftOPSubTrans_store(s, d, r, c)     fftOPSubTrans(s, d, r, c)

#endif	/* !FFT_DOUBLE_PREC */

#else	/* !VECTOR_TRANSPOSE_ENABLE */

/* 
 * Platform- and precision-independent routine to copy FFT_FLOATS_PER_VECTOR
 * continguous FFTFloats from *row to consecutive columns at *col.
 */
static FFT_INLINE void swapRowCol(
	const FFTFloat *row,
	FFTFloat *col,
	size_t rowSize)		// in floats
{
	for(unsigned dex=0; dex<FFT_FLOATS_PER_VECTOR; dex++) {
		*col = *row++;
		col += rowSize;
	}
}


/* 
 * Transpose one submatrix (FFT_FLOATS_PER_SUBMATRIX on a side)
 * out of place. 
 */
static void fftOPSubTrans(
	const FFTFloat		*src,
	FFTFloat			*dst,
	size_t				srcRowSize,		// src, in FFTFloats, a.k.a. src numCols
	size_t				dstRowSize)		// dst, in FFTFloats, a.k.a. dst numCols
{
	dumpSub("fftOPSubTrans start", src, srcRowSize);
	unsigned colIncr = FFT_FLOATS_PER_VECTOR * dstRowSize;
	
	/* row and col refer to coordinates in src */
	for(unsigned row=0; row<FFT_FLOATS_PER_SUBMATRIX; row++) {
		/*
		 * ip : source : moves across a row in src
		 * op : dest   : moves down a column in dst, by its rowSize, which is 'our' 
		 *               colSize 
		 */
		const FFTFloat *ip = src + (row * srcRowSize);
		fft_prefetch(ip);
		FFTFloat *op = dst + row;

		for(unsigned col=0; col<FFT_FLOATS_PER_SUBMATRIX; col+=FFT_FLOATS_PER_VECTOR) {
			swapRowCol(ip, op, dstRowSize);
			ip += FFT_FLOATS_PER_VECTOR;
			op += colIncr;
		}
	}
	dumpSub("fftOPSubTrans end", dst, dstRowSize);
}

#define fftOPSubTrans_store(s, d, r, c)     fftOPSubTrans(s, d, r, c)

#endif	/* !VECTOR_TRANSPOSE_ENABLE */

/* 
 * Optimized out-of-place transpose using submatrices. The size of the 
 * matrix must be >= FFT_FLOATS_PER_SUBMATRIX in each dimension and 
 * also be a multiple of FFT_FLOATS_PER_SUBMATRIX in each dimension.
 */
static MFFTReturn fftOPSwapRowsColsOpt(
	const FFTFloat		*src,
	FFTFloat			*dst,
	/* total size of src */
	size_t				numRows,
	size_t				numCols,
	/* portion to move - must be FFT_FLOATS_PER_SUBMATRIX aligned */
	size_t				startRow,
	size_t				rowsToMove)
{
	size_t numRowSubs = rowsToMove / FFT_FLOATS_PER_SUBMATRIX;
	size_t numColSubs = numCols / FFT_FLOATS_PER_SUBMATRIX;
	size_t row = startRow;
	unsigned colSubIncr = FFT_FLOATS_PER_SUBMATRIX * numRows;

	dumpSubSwapPrint("fftOPSwapRowsColsOpt top: startRow %lu rowsToMove %lu src %p dst %p\n", 
		(unsigned long)startRow, (unsigned long)rowsToMove,
		src, dst);

	RFASSERT((numRows % FFT_FLOATS_PER_SUBMATRIX) == 0);
	RFASSERT((numCols % FFT_FLOATS_PER_SUBMATRIX) == 0);
	RFASSERT((startRow % FFT_FLOATS_PER_SUBMATRIX) == 0);
	RFASSERT((rowsToMove % FFT_FLOATS_PER_SUBMATRIX) == 0);
	
	for(size_t rowSubDex=0; rowSubDex<numRowSubs; rowSubDex++) {
	
		const FFTFloat *srcSub = src + (row * numCols);
		FFTFloat *dstSub = dst + row;
		
		for(size_t colSubDex=0; colSubDex<numColSubs; colSubDex++) {
			fftOPSubTrans(srcSub, dstSub, numCols, numRows);
			srcSub += FFT_FLOATS_PER_SUBMATRIX;
			dstSub += colSubIncr;
		}
		row += FFT_FLOATS_PER_SUBMATRIX;
	}
	return MR_Success;
}

/* Called out from thread module */
static MFFTReturn transOPThr(TP_TaskUnion *u)
{
	TP_TransposeSplit *ttop  = &u->transS;
	fftOPSwapRowsColsOpt(ttop->src, ttop->dst, 
		ttop->numRows, ttop->numCols,
		ttop->startRow, ttop->rowsToMove);
	return MR_Success;
}

/* 
 * Threaded out-of-place transpose using submatrices. The size of the 
 * matrix must be >= FFT_FLOATS_PER_SUBMATRIX in each dimension and 
 * also be a multiple of FFT_FLOATS_PER_SUBMATRIX in each dimension.
 */
static MFFTReturn fftOPSwapRowsColsThr(
	MatrixFFTPlan		mfftPlan,
	const FFTFloat		*src,
	FFTFloat			*dst,
	size_t				numRows,		// rows in src
	size_t				numCols)		// columns in src
{
	tpThreadDebug("fftOPSwapRowsColsThr top: numRows %lu numCols %lu src %p\n",
		(unsigned long)numRows, (unsigned long)numCols, src);

	MFFTReturn ourRtn = MR_Success;
	size_t numRowSubs = numRows / FFT_FLOATS_PER_SUBMATRIX;
	RFASSERT(mfftPlan->threadPoolP->numThreads < numRowSubs);

	unsigned numThreads = mfftPlan->threadPoolP->numThreads;
	if(numThreads > numRowSubs) {
		numThreads = numRowSubs;
	}
	
	/* submatrices per thread, we (main) get remainder */
	size_t subsPerThread = numRowSubs / numThreads;
	
	/* rows per thread (ours might be bigger) */
	size_t rowsPerThread = subsPerThread * FFT_FLOATS_PER_SUBMATRIX;

	/* set up a TP_TransposeSplit for each thread */	
	for(unsigned dex=0; dex<numThreads; dex++) {
		TP_PerThread *pt   = &mfftPlan->threadPoolP->perThread[dex];
		TP_Task *task      = &pt->task;
		TP_TransposeSplit *ttop = &task->u->transS;

		task->op		 = TPO_App;
		task->threadFcn  = transOPThr;
		ttop->src		 = src;
		ttop->dst		 = dst;
		ttop->numRows	 = numRows;
		ttop->numCols	 = numCols;
		ttop->startRow	 = dex * subsPerThread * FFT_FLOATS_PER_SUBMATRIX;
		ttop->rowsToMove = rowsPerThread;
	}
	
	/* last one might have more rows */
	TP_TransposeSplit *ttop = &mfftPlan->threadPoolP->perThread[numThreads-1].task.u->transS;
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

/* Force brute-force for debug */
#define OP_TRANSPOSE_BASIC					0

/* 
 * Brute force OOP transpose, for small or non-Intel ops
 */
static MFFTReturn fftTransposeOPSmall(
	const FFTFloat		*srcBuf,
	FFTFloat			*dstBuf,
	size_t				srcRows,
	size_t				srcCols)
{
	for(size_t row=0; row<srcRows; row++) {
		/*
		 * ip : source : moves across a row in src
		 * op : des    : moves down a column in dst 
		 */
		const FFTFloat *ip = srcBuf + (row * srcCols);	
		FFTFloat *op = dstBuf + row;
		for(size_t col=0; col<srcCols; col++) {
			*op = *ip++;
			op += srcRows;
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
	if(!VECTOR_TRANSPOSE_ENABLE || OP_TRANSPOSE_BASIC) {
		/* For now, use brute force */
		fftTransposeOPSmall(srcBuf->real, dstBuf->real, srcRows, srcCols);
		return fftTransposeOPSmall(srcBuf->imag, dstBuf->imag, srcRows, srcCols);
	}
	
	MFFTReturn mrtn = MR_Success;
	if((srcRows >= FFT_FLOATS_PER_SUBMATRIX) &&
	   (srcCols >= FFT_FLOATS_PER_SUBMATRIX)) {
		unsigned numThreads = mfftPlan->threadPoolP->numThreads;	/* # of pthreads */
		if((numThreads > 1) &&										/* threading enabled */
		   ((srcRows / FFT_FLOATS_PER_SUBMATRIX) > numThreads))	{	/* more subs then threads */
			mrtn = fftOPSwapRowsColsThr(mfftPlan, srcBuf->real, dstBuf->real, srcRows, srcCols);
			if(mrtn) {
				return mrtn;
			}
			return fftOPSwapRowsColsThr(mfftPlan, srcBuf->imag, dstBuf->imag, srcRows, srcCols);
		}
		else {
			mrtn = fftOPSwapRowsColsOpt(srcBuf->real, dstBuf->real, srcRows, srcCols, 0, srcRows);
			if(mrtn) {
				return mrtn;
			}
			return fftOPSwapRowsColsOpt(srcBuf->imag, dstBuf->imag, srcRows, srcCols, 0, srcRows);
		}
	}
	else {
		mrtn = fftTransposeOPSmall(srcBuf->real, dstBuf->real, srcRows, srcCols);
		if(mrtn) {
			return mrtn;
		}
		return fftTransposeOPSmall(srcBuf->imag, dstBuf->imag, srcRows, srcCols);
	}
}

#if !TRANS_STRIPE_FORCE_BASIC

/*
 * Transpose a columnar stripe to a stripe buffer
 */
static void fftTransposeFromColumnOpt(
	const FFTFloat		*src,		/* numRows x numCols */
	FFTFloat			*dst,		/* colsToMove rows x numRows cols */
	/* total size of src */
	size_t				numRows,
	size_t				numCols,
	/* portion to move - in src - must be FFT_FLOATS_PER_SUBMATRIX aligned */
	size_t				startCol,
	size_t				colsToMove)
{
	size_t numRowSubs = numRows    / FFT_FLOATS_PER_SUBMATRIX;
	size_t numColSubs = colsToMove / FFT_FLOATS_PER_SUBMATRIX;
	size_t srcRow = 0;
	size_t colSubIncr = FFT_FLOATS_PER_SUBMATRIX * numRows;

	dumpSubSwapPrint("fftOPTransposeFromColumnOpt top: startCol %lu colsToMove %lu floatsPerSub %u\n", 
		(unsigned long)startCol, (unsigned long)colsToMove, (unsigned)FFT_FLOATS_PER_SUBMATRIX);

	RFASSERT((numRows    % FFT_FLOATS_PER_SUBMATRIX) == 0);
	RFASSERT((numCols    % FFT_FLOATS_PER_SUBMATRIX) == 0);
	RFASSERT((startCol   % FFT_FLOATS_PER_SUBMATRIX) == 0);
	RFASSERT((colsToMove % FFT_FLOATS_PER_SUBMATRIX) == 0);

	for(size_t rowSubDex=0; rowSubDex<numRowSubs; rowSubDex++) {
		const FFTFloat *srcSub = src + (srcRow * numCols) + startCol;
		FFTFloat *dstSub = dst + srcRow;
		
		for(size_t colSubDex=0; colSubDex<numColSubs; colSubDex++) {
			fftOPSubTrans_store(srcSub, dstSub, numCols, numRows);
			srcSub += FFT_FLOATS_PER_SUBMATRIX;		// to right, one submatrix
			dstSub += colSubIncr;					// down, one submatrix
		}
		srcRow += FFT_FLOATS_PER_SUBMATRIX;
	}
}

/*
 * Transpose a stripe buffer to a columnar stripe
 */
static void fftTransposeToColumnOpt(
	const FFTFloat		*src,		/* colsToMove rows x numRows cols  */
	FFTFloat			*dst,		/* numRows x numCols */
	/* total size of dst */
	size_t				numRows,
	size_t				numCols,
	/* portion to move - in dst - must be FFT_FLOATS_PER_SUBMATRIX aligned */
	size_t				startCol,
	size_t				colsToMove)
{
	size_t numRowSubs = numRows    / FFT_FLOATS_PER_SUBMATRIX;
	size_t numColSubs = colsToMove / FFT_FLOATS_PER_SUBMATRIX;
	size_t dstRow = 0;
	size_t colSubIncr = FFT_FLOATS_PER_SUBMATRIX * numRows;

	dumpSubSwapPrint("fftOPTransposeFromColumn top: startCol %lu colsToMove %lu floatsPerSub %u\n", 
		(unsigned long)startCol, (unsigned long)colsToMove, (unsigned)FFT_FLOATS_PER_SUBMATRIX);

	RFASSERT((numRows    % FFT_FLOATS_PER_SUBMATRIX) == 0);
	RFASSERT((numCols    % FFT_FLOATS_PER_SUBMATRIX) == 0);
	RFASSERT((startCol   % FFT_FLOATS_PER_SUBMATRIX) == 0);
	RFASSERT((colsToMove % FFT_FLOATS_PER_SUBMATRIX) == 0);

	for(size_t rowSubDex=0; rowSubDex<numRowSubs; rowSubDex++) {
	
		const FFTFloat *srcSub = src + (rowSubDex * FFT_FLOATS_PER_SUBMATRIX);
		FFTFloat *dstSub = dst + (dstRow * numCols) + startCol;

		for(size_t colSubDex=0; colSubDex<numColSubs; colSubDex++) {
			fftOPSubTrans(srcSub, dstSub, numRows, numCols);
			dstSub += FFT_FLOATS_PER_SUBMATRIX;		// to right, one submatrix
			srcSub += colSubIncr;					// down, one submatrix
		}
		dstRow += FFT_FLOATS_PER_SUBMATRIX;
	}
}

#endif	/* TRANS_STRIPE_FORCE_BASIC */

/*
 * Brute force implementations of stripe transpose
 */
static MFFTReturn fftTransposeFromColumnBasic(
	const FFTComplex	*src,		/* numRows x numCols */
	FFTComplex			*dst,		/* colsToMove rows x numRows cols */
	/* total size of src */
	size_t				numRows,
	size_t				numCols,
	/* portion to move - in src - must be FFT_FLOATS_PER_SUBMATRIX aligned */
	size_t				startCol,
	size_t				colsToMove)
{
	const FFTFloat *realIn;
	const FFTFloat *imagIn;
	FFTFloat *realOut;
	FFTFloat *imagOut;
	
	for(size_t srcRow=0; srcRow<numRows; srcRow++) {
		realIn  = src->real + (srcRow * numCols) + startCol;
		imagIn  = src->imag + (srcRow * numCols) + startCol;
		realOut = dst->real + srcRow;
		imagOut = dst->imag + srcRow;
		for(size_t colDex=0; colDex<colsToMove; colDex++) {
			*realOut = *realIn++;
			*imagOut = *imagIn++;
			realOut += numRows;
			imagOut += numRows;
		}
	}
	return MR_Success;
}

/*
 * Transpose a stripe buffer to a columnar stripe
 */
static MFFTReturn fftTransposeToColumnBasic(
	const FFTComplex	*src,		/* colsToMove rows x numRows cols  */
	FFTComplex			*dst,		/* numRows x numCols */
	/* total size of dst */
	size_t				numRows,
	size_t				numCols,
	/* portion to move - in dst - must be FFT_FLOATS_PER_SUBMATRIX aligned */
	size_t				startCol,
	size_t				colsToMove)
{
	const FFTFloat *realIn;
	const FFTFloat *imagIn;
	FFTFloat *realOut;
	FFTFloat *imagOut;

	for(size_t dstRow=0; dstRow<numRows; dstRow++) {
		realOut = dst->real + (dstRow * numCols) + startCol;
		imagOut = dst->imag + (dstRow * numCols) + startCol;
		realIn  = src->real + dstRow;
		imagIn  = src->imag + dstRow;
		for(size_t colDex=0; colDex<colsToMove; colDex++) {
			*realOut++ = *realIn;
			*imagOut++ = *imagIn;
			realIn += numRows;
			imagIn += numRows;
		}
	}
	return MR_Success;
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
	#if		!TRANS_STRIPE_FORCE_BASIC
	if((numRows     >= FFT_FLOATS_PER_SUBMATRIX) &&
	   (numCols     >= FFT_FLOATS_PER_SUBMATRIX) &&
	   (colsToMove  >= FFT_FLOATS_PER_SUBMATRIX) &&
	   ((colsToMove & (FFT_FLOATS_PER_SUBMATRIX - 1)) == 0) &&
	   ((startCol   & (FFT_FLOATS_PER_SUBMATRIX - 1)) == 0)) {
		fftTransposeFromColumnOpt(src->real, dst->real, numRows, numCols, startCol, colsToMove);
		fftTransposeFromColumnOpt(src->imag, dst->imag, numRows, numCols, startCol, colsToMove);
		return MR_Success;
	}
	/* else drop thru and use the brute force implementation */
	#endif	/* TRANS_STRIPE_FORCE_BASIC */

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
	#if		!TRANS_STRIPE_FORCE_BASIC
	if((numRows     >= FFT_FLOATS_PER_SUBMATRIX) &&
	   (numCols     >= FFT_FLOATS_PER_SUBMATRIX) &&
	   (colsToMove  >= FFT_FLOATS_PER_SUBMATRIX) &&
	   ((colsToMove & (FFT_FLOATS_PER_SUBMATRIX - 1)) == 0) &&
	   ((startCol   & (FFT_FLOATS_PER_SUBMATRIX - 1)) == 0)) {
		fftTransposeToColumnOpt(src->real, dst->real, numRows, numCols, startCol, colsToMove);
		fftTransposeToColumnOpt(src->imag, dst->imag, numRows, numCols, startCol, colsToMove);
		return MR_Success;
	}
	/* else drop thru and use the brute force implementation */
	#endif	/* TRANS_STRIPE_FORCE_BASIC */

	return fftTransposeToColumnBasic(src, dst, numRows, numCols, startCol, colsToMove);
}

#endif	/* FFT_SPLIT_COMPLEX */

