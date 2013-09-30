/*	File: fftTransposeIPInt.cpp  
	
	Description:
		Square in-place transpose, interleaved version.
	
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
 * fftTransposeIPInt.cpp - square in-place transpose, interleaved version. 
 *
 * Created Sep. 26 2008. 
 */

/*
 * Submatrices - the ones we work on at the lowest level - are squares whose sides
 * are cache line sized. Row/column transposition is performed on one submatrix
 * at a time for in-place transpositions (the ones on the x=y diagonal of the whole
 * 2D matrix), or two submatrices - one above and one below the x=y diagonal for
 * other operations. The idea is that the submatrices we're operating on fit entirely
 * within cache (hopefully L1 cache), avoiding the cache misses normally taken
 * when processing an entire column of the whole 2D matrix.
 */

#include "fftTranspose.h"
#include "ThreadPool.h"
#include "fftPriv.h"
#include "fftDebug.h"
#include "fftIntel.h"

#if		!FFT_SPLIT_COMPLEX

#if FFT_INTEL 

#pragma mark --- Intel precision-dependent optimizations ---

#if	FFT_DOUBLE_PREC

/*
 * Transpose one submatrix in place.
 * Intel, double precision, interleaved complex. 
 * Used for submatrices containing the x=y diagonal.
 * The submatrix's diagonal is untouched. 
 */
static void fftSubTransInPlace(
	FFTComplex *sub,
	size_t rowSize)			/* offset between rows, in FFTComplex */
{
	dumpSub("fftSubTransInPlace start", sub, rowSize);

	for(unsigned row=0; row<FFT_COMPLEX_PER_SUBMATRIX; row++) {
		/* 
		 * start row pointer at [row, row+1]
		 * start column pointer at [row+1, row]
		 */
		FFTVector *rowPtr = (FFTVector *)(sub + (row * rowSize) + row + 1);
		fft_prefetch(rowPtr);
		/* down one row, back one column */
		FFTVector *colPtr = rowPtr + rowSize - 1;
		fft_prefetch(colPtr);
		
		/*
		 * Note everything is vector-aligned since the basic element in the 
		 * array is vector-sized; no special case for elements between the 
		 * diagonal and the next vector-aligned element as is done in the 
		 * single precision case. Proceed with the swap.
		 */
		for(size_t col=row+1; col<FFT_COMPLEX_PER_SUBMATRIX; col++) {
			FFTVector vTmp = *rowPtr;
			*rowPtr++ = *colPtr;
			*colPtr = vTmp;
			colPtr += rowSize;
		}
	}
	dumpSub("fftSubTransInPlace end", sub, rowSize);
}

/*
 * Transpose rows and columns between two submatrices. The two submatrices
 * are from symmtric positions relative to the x=y diagonal in the original 
 * array; rows in sub1 become columns in sub2 and vice versa.
 * Intel, double precision, interleaved complex. 
 * This is really easy since an FFTComplex is the same size as an FFTVector;
 * there is no shuffling of vector components. 
 */
static void fftSubSwap(
	FFTComplex *sub1,
	FFTComplex *sub2,
	size_t rowSize)			/* offset between rows, in FFTComplex */
{
	dumpSub("fftSubSwap start sub1", sub1, rowSize);
	dumpSub("fftSubSwap start sub2", sub2, rowSize);

	/* rowDex, colDex refer to sub1 */
	for(size_t rowDex=0; rowDex<FFT_COMPLEX_PER_SUBMATRIX; rowDex++) {
		FFTVector *rowPtr = (FFTVector *)(sub1 + (rowSize * rowDex));
		FFTVector *colPtr = (FFTVector *)(sub2 + rowDex);
		for(size_t colDex=0; colDex<FFT_COMPLEX_PER_SUBMATRIX; colDex++) {
			FFTVector vTmp = *rowPtr;
			*rowPtr++ = *colPtr;
			*colPtr = vTmp;
			colPtr += rowSize;
		}
	}

	dumpSub("fftSubSwap end sub1", sub1, rowSize);
	dumpSub("fftSubSwap end sub2", sub2, rowSize);
}

#else	/* !FFT_DOUBLE_PREC */

/*
 * Optimized, vector-based implementation for Intel, single-precision only.
 */
 
/* mask for detecting FFTComplex counts which are not vector aligned */
#define COMPLEX_VECTOR_MASK	0x01

/* 
 * Swap 2 contiguous FFTComplex's at *row with 2 FFTComplex's. 
 * in consecutive rows at *col.
 * This is a bit of a hack which assumes/knows that the size of 
 * a FFTComplex is the same as the size of a double. Since we're
 * just moving data, not computing anything with it, this is OK...
 */
static void swapRowCol(
	FFTComplex *_row,
	FFTComplex *_col,
	size_t rowSize)		// in doubles
{
	double *row = (double *)_row;
	double *col = (double *)_col;
	__m128d cval = _mm_set_pd(col[rowSize], col[0]); 
	__m128d *rd = (__m128d *)row;
	__m128d rval = *rd;
	*rd = cval;
	_mm_storel_pd(col, rval);
	_mm_storeh_pd(col+rowSize, rval); 
}

/*
 * Transpose one submatrix in place.
 * Used for submatrices containing the x=y diagonal.
 * The submatrix's diagonal is untouched. 
 */
static void fftSubTransInPlace(
	FFTComplex *sub,
	size_t rowSize)			/* offset between rows, in FFTComplex */
{
	dumpSub("fftSubTransInPlace start", sub, rowSize);
	size_t colIncr = FFT_COMPLEX_PER_VECTOR * rowSize;
	
	for(unsigned row=0; row<FFT_COMPLEX_PER_SUBMATRIX; row++) {
		/* 
		 * start row pointer at [row, row+1]
		 * start column pointer at [row+1, row]
		 */
		FFTComplex *rowPtr = sub + (row * rowSize) + row + 1;
		fft_prefetch(rowPtr);
		/* down one row, back one column */
		FFTComplex *colPtr = rowPtr + rowSize - 1;
		fft_prefetch(colPtr);
		
		/* 
		 * Do some swaps outside of the loop on rows on which the current
		 * rowPtr is not vector-aligned so inner loop processes a vector-sized
		 * number of elements (and works on vector-aligned pointers). 
		 */
		FFTComplex tmp;
		size_t col = row + 1;
		size_t fillerRows = (FFT_COMPLEX_PER_VECTOR - col) & COMPLEX_VECTOR_MASK;
		
		for(unsigned dex=0; dex<fillerRows; dex++) {
			tmp       = *rowPtr;
			*rowPtr++ = *colPtr;
			*colPtr   = tmp;
			colPtr   += rowSize;
			col++;
		}

		for(; col<FFT_COMPLEX_PER_SUBMATRIX; col+=FFT_COMPLEX_PER_VECTOR) {
			swapRowCol(rowPtr, colPtr, rowSize);
			rowPtr += FFT_COMPLEX_PER_VECTOR;
			colPtr += colIncr;
		}
	}
	dumpSub("fftSubTransInPlace end", sub, rowSize);
}

/*
 * Transpose rows and columns between two submatrices. The two submatrices
 * are from symmtric positions relative to the x=y diagonal in the original 
 * array; rows in sub1 become columns in sub2 and vice versa.
 */
static void fftSubSwap(
	FFTComplex *sub1,
	FFTComplex *sub2,
	size_t rowSize)			/* offset between rows, in FFTComplex */
{
	dumpSub("fftSubSwap start sub1", sub1, rowSize);
	dumpSub("fftSubSwap start sub2", sub2, rowSize);
	unsigned colIncr = FFT_COMPLEX_PER_VECTOR * rowSize;

	/* row and col refer to coordinates in sub1 */
	for(unsigned row=0; row<FFT_COMPLEX_PER_SUBMATRIX; row++) {
		FFTComplex *row1Ptr = sub1 + (row * rowSize);
		fft_prefetch(row1Ptr);
		FFTComplex *col2Ptr = sub2 + row;

		fft_prefetch(col2Ptr);
		for(size_t col=0; col<FFT_COMPLEX_PER_SUBMATRIX; col+=FFT_COMPLEX_PER_VECTOR) {
			swapRowCol(row1Ptr, col2Ptr, rowSize);
			row1Ptr += FFT_COMPLEX_PER_VECTOR;
			col2Ptr += colIncr;
		}
	}
	dumpSub("fftSubSwap end sub1", sub1, rowSize);
	dumpSub("fftSubSwap end sub2", sub2, rowSize);
}

#endif	/* FFT_DOUBLE_PREC */

#pragma mark --- Intel precision-independent optimizations ---

/*
 * Basic top-level row/column swapper for a single square array of FFTFloats.
 * Size of the array and the portion to transform must be power-of-two multiples
 * of FFT_COMPLEX_PER_SUBMATRIX.
 */
static MFFTReturn fftSwapRowsColumnsOpt(
	FFTComplex		*buf,
	size_t			rowSize,		/* offset between rows, in FFTComplex */
	/* portion to move - must be FFT_COMPLEX_PER_SUBMATRIX aligned */
	size_t			startRow,
	size_t			rowsToMove)
{
	size_t numRowSubs = rowsToMove / FFT_COMPLEX_PER_SUBMATRIX;
	size_t numColSubs = rowSize / FFT_COMPLEX_PER_SUBMATRIX;
	size_t row = startRow;
	
	tpThreadDebug("fftSwapRowsColumnsOpt top: startRow %lu rowsToMove %lu\n", 
		(unsigned long)startRow, (unsigned long)rowsToMove);
	RFASSERT((numRowSubs * FFT_COMPLEX_PER_SUBMATRIX) == rowsToMove);
	RFASSERT((numColSubs * FFT_COMPLEX_PER_SUBMATRIX) == rowSize);
	
	/*
	 * -- base is the upper left corner of the submatrix on the diagonal we start 
	 *    with at the top of each outer loop.
	 * -- sub1 is the upper left corner of submatrices to the right of base.
	 * -- sub2 is the upper left corner of submatrices below base.
	 *
	 * All pointer arithmetic is in FFTComplex.
	 * rowSubDex and colSubDex are indices into an array of submatrices.
	 */
	size_t startRowSub = startRow / FFT_COMPLEX_PER_SUBMATRIX;
	RFASSERT((startRowSub * FFT_COMPLEX_PER_SUBMATRIX) == startRow);
	size_t endRowSub = startRowSub + numRowSubs;
	size_t sub2Incr = FFT_COMPLEX_PER_SUBMATRIX * rowSize;

	/* this points to the submatrix on the diagonal */
	FFTComplex *base = buf + (row * rowSize) + row;
	
	/* increment between successive base ptrs */
	size_t baseIncr = FFT_COMPLEX_PER_SUBMATRIX * (rowSize + 1);
	
	for(size_t rowSubDex=startRowSub; 
	    rowSubDex<endRowSub; 
		rowSubDex++, row+=FFT_COMPLEX_PER_SUBMATRIX) {
		
		/* swap the submatrix on the diagonal in place */
		fftSubTransInPlace(base, rowSize);
		
		/* now swap between [x,y] and [y,x] */
		/* sub1 is 'n' submatrices to the right of base */
		FFTComplex *sub1 = base + FFT_COMPLEX_PER_SUBMATRIX;
		/* sub2 is 'n' submatrices down from base */
		FFTComplex *sub2 = base + sub2Incr;
		for(size_t colSubDex=rowSubDex+1; colSubDex<numColSubs; colSubDex++) { 
			fftSubSwap(sub1, sub2, rowSize);
			sub1 += FFT_COMPLEX_PER_SUBMATRIX;
			sub2 += sub2Incr;
		}
		base += baseIncr;
	}
	RFASSERT(row == (startRow + rowsToMove));
	return MR_Success;
}

/* Called out from thread module */
static MFFTReturn swapRowsColsThr(TP_TaskUnion *u)
{
	TP_TransposeInt *ttop  = &u->transI;
	fftSwapRowsColumnsOpt(ttop->dst, ttop->numRows, ttop->startRow, ttop->rowsToMove);
	return MR_Success;
}

#pragma mark --- Square threaded transpose ---

/* 
 * Transpose, square, threaded. 
 *
 * Breaking up a square into equal sized chunks of work, one for each thread, 
 * is not as simple as dividing up the rows as in other threaded routines. Here,
 * for example, see the unit square broken up into 4 segments, each of which is 
 * a portion of a square which in this case has its origin at the lower
 * right (the last row, actually the last submatrix, to get processed). 
 * To break up the total work into two equal chunks, we'd want a square
 * which is half the area of the total square; such a square has a side
 * of length sqrt(1/2), so that the small square has area 0.5. 
 * Similarly, a square with an area of 1/4 of the total square has sides 
 * of length sqrt(1/4). 
 *
 *   .................... <--- 1.0
 *   .                  .
 *   .   ...............x <--- sqrt(3/4)
 *   .   .              .
 *   .   .   ...........x <--- sqrt(2/4)
 *   .   .   .          .
 *   .   .   .   .......x <--- sqrt(1/4)
 *   .   .   .   .      .
 *   .................... <--- 0.0
 *
 * Since we actually start with row 0 at the top, what we really want for 
 * the above numRows coefficients is (1 - sqrt(3/4)), (1 - sqrt(1/2)), etc.
 *
 * To avoid calculating sqrt() functions every time we come through this routine,
 * the above cofficients are calculated once in fftRowThreadInit(). 
 */
 
/* minimum subs/numThreads to use squareCoeff algorithm */
#define MIN_SUBS_PER_THREAD_COEFF	16

static MFFTReturn fftSwapRowsColumnsThr(
	MatrixFFTPlan mfftPlan,
	FFTComplex	*buf,
	size_t		numRows)			/* offset between rows, in FFTComplex */
{
	tpThreadDebug("fftSwapRowsColumnsThr top: numRows %lu\n", (unsigned long)numRows);

	size_t numRowSubs = numRows / FFT_COMPLEX_PER_SUBMATRIX;
	RFASSERT((numRowSubs * FFT_COMPLEX_PER_SUBMATRIX) == numRows);
	RFASSERT(mfftPlan->threadPoolP->numThreads > 0);
	MFFTReturn ourRtn = MR_Success;
	
	unsigned numThreads = mfftPlan->threadPoolP->numThreads;
	if(numThreads > numRowSubs) {
		numThreads = numRowSubs;
	}
	
	/* 
	 * starting submatrix and submatrices per thread
	 */
	size_t numSubs[numThreads];
	size_t startSub[numThreads];
	
	if(numRowSubs < (numThreads * MIN_SUBS_PER_THREAD_COEFF)) {
		/* too small to use the squareCoeff algorithm; divide equally */
		size_t subsPerThread = numRowSubs / numThreads;
		for(unsigned dex=0; dex<numThreads; dex++) {
			numSubs[dex]  = subsPerThread;
			startSub[dex] = dex * subsPerThread;
		}
		/* last one, it might have more or less rows */
		numSubs[numThreads - 1] = numRowSubs - (subsPerThread * (numThreads - 1));
	}
	else {
		for(unsigned dex=0; dex<numThreads; dex++) {
			float squareCoeff = mfftPlan->threadPoolP->perThread[dex].squareCoeff;
			RFASSERT((squareCoeff >= 0.0) && (squareCoeff < 1.0));
			startSub[dex] = (size_t)((float)numRowSubs * squareCoeff);
		}
		size_t totalSubs = 0;
		for(unsigned dex=0; dex<(numThreads-1); dex++) {
			numSubs[dex] = startSub[dex+1] - startSub[dex];
			totalSubs += numSubs[dex];
		}
		numSubs[numThreads-1] = numRowSubs - totalSubs;
	}
	
	/* set up a TP_TransposeInt for each thread */	
	for(unsigned dex=0; dex<numThreads; dex++) {
		TP_PerThread *pt		= &mfftPlan->threadPoolP->perThread[dex];
		TP_Task *task			= &pt->task;
		TP_TransposeInt *ttop	= &task->u->transI;

		task->op         = TPO_App;
		task->threadFcn  = swapRowsColsThr;
		ttop->src        = NULL;		/* not used, make sure of it */
		ttop->dst        = buf;
		ttop->numRows    = numRows;
		ttop->numCols    = 0;		/* not used */
		ttop->startRow   = startSub[dex] * FFT_COMPLEX_PER_SUBMATRIX;
		ttop->rowsToMove = numSubs[dex] * FFT_COMPLEX_PER_SUBMATRIX;
	}
	
	/* GO */
	ourRtn = tpThreadDispatch(mfftPlan->threadPoolP, numThreads);
	if(ourRtn) {
		return ourRtn;
	}
	ourRtn = tpThreadFinish(mfftPlan->threadPoolP, numThreads);
	tpThreadDebug("fftSwapRowsColumnsThr end\n");
	return ourRtn;
}

#endif	/* FFT_INTEL */

/* Brute force manual implementation */
static MFFTReturn fftTransposeSquareSmall(
	FFTComplex			*buf,
	size_t				size)		/* size of each side in complex elements */
{
	for(size_t row=0; row<size; row++) {
		/* Start at (row, row+1) */
		FFTComplex *top = buf + (row * size) + row + 1;
		/* bottom is one down and one to the left */
		FFTComplex *bot = top + size - 1;
		for(size_t col=row+1; col<size; col++) {
			FFTComplex tmp = *top;
			*top++ = *bot;
			*bot   = tmp;
			bot   += size;
		}
	}
	return MR_Success;
}

#pragma mark --- Public API ---

/* minimum submatrices for threading */
#define THREAD_MIN_SUBS							2

MFFTReturn fftTransposeSquare(
	MatrixFFTPlan		mfftPlan,
	FFTComplex			*buf,
	size_t				numRows)		/* size of each side in complex elements */
{
	#if	!FFT_INTEL
	/* For now, use brute force */
	return fftTransposeSquareSmall(buf, numRows);
	#else
	
	if(numRows < FFT_COMPLEX_PER_SUBMATRIX) {
		/* too small to use submatrices */
		return fftTransposeSquareSmall(buf, numRows);
	}
	
	size_t numSubs = numRows / FFT_COMPLEX_PER_SUBMATRIX;
	if((mfftPlan->threadPoolP->numThreads > 1) && (numSubs >= THREAD_MIN_SUBS)) {
		return fftSwapRowsColumnsThr(mfftPlan, buf, numRows);
	}

	/* single threaded, submatrix optimization */
	return fftSwapRowsColumnsOpt(buf, numRows, 0, numRows);
	
	#endif	/* FFT_INTEL */
}

#endif	/* !FFT_SPLIT_COMPLEX */
