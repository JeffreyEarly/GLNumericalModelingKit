/*	File: fftTransposeIPSplit.cpp  
	
	Description:
		Square in-place transpose, split complex version.
	
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
 * fftTransposeIPSplit.cpp - square in-place transpose, split complex version
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

#if		FFT_SPLIT_COMPLEX

/* mask for detecting FFTFloat counts which are not vector aligned */
#if		FFT_DOUBLE_PREC
#define FFT_FLOAT_VECTOR_MASK	0x01
#else
#define FFT_FLOAT_VECTOR_MASK	0x03
#endif

#if	FFT_INTEL 

#pragma mark --- Intel-optimized swapRowCol() routines ---

#if	FFT_DOUBLE_PREC

/* 
 * Swap 2 (FFT_FLOATS_PER_VECTOR) contiguous doubles at *row with 2 doubles 
 * in consecutive rows at *col.
 */
static void swapRowCol(
	double *row,
	double *col,
	size_t rowSize)		// in doubles
{
	FFTVector cval = _mm_set_pd(col[rowSize], col[0]); 
	FFTVector *rd = (FFTVector *)row;
	FFTVector rval = *rd;
	*rd = cval;
	_mm_storel_pd(col, rval);
	_mm_storeh_pd(col+rowSize, rval); 
}

#else	/* !FFT_DOUBLE_PREC */


/* 
 * Swap 4 (FFT_FLOATS_PER_VECTOR) contiguous floats at *row with 4 floats 
 * in consecutive columns at *col.
 * Intel, single precision only.
 */
static void swapRowCol(
	float *row,
	float *col,
	size_t rowSize)		// in floats
{
	FFTVectUnion rowIn, colIn;
	
	/* read 4 floats from row at once */
	rowIn.v = *(__m128 *)row;
	
	/* assemble 4 floats from col */
	float *uvf = colIn.f;
	float *cp  = col;
	*uvf++ = *cp; cp += rowSize;
	*uvf++ = *cp; cp += rowSize;
	*uvf++ = *cp; cp += rowSize;
	*uvf   = *cp;
	
	/* assembled col --> row */
	*(__m128 *)row = colIn.v;
	
	/* disassemble row --> col */
	cp  = col;
	uvf = rowIn.f;
	*cp = *uvf++; cp += rowSize;
	*cp = *uvf++; cp += rowSize;
	*cp = *uvf++; cp += rowSize;
	*cp = *uvf;
}

#endif	/* FFT_DOUBLE_PREC */

#else	/* !FFT_INTEL */

#pragma mark --- Platform-independent swapRowCol() ---

/*
 * Brute force, platform independent routine to swap FFT_FLOATS_PER_VECTOR
 * contiguous FFTFloats at *row with FFT_FLOATS_PER_VECTOR FFTFloats 
 * in consecutive columns at *col.
 */
static FFT_INLINE void swapRowCol(
	FFTFloat *row,
	FFTFloat *col,
	size_t rowSize)		// in floats
{
	for(unsigned dex=0; dex<FFT_FLOATS_PER_VECTOR; dex++) {
		FFTFloat tmp = *row;
		*row++ = *col;
		*col = tmp;
		col += rowSize;
	}
}

#endif	/* !FFT_INTEL */

#pragma mark --- Platform-independent optimized transform support ---

/*
 * Transpose one submatrix in place.
 * Used for submatrices containing the x=y diagonal.
 * The submatrix's diagonal is untouched. 
 */
static void fftSubTransInPlace(
	FFTFloat *sub,
	size_t rowSize)			/* offset between rows, in FFTComplex */
{
	dumpSub("fftSubTransInPlace start", sub, rowSize);
	RFASSERT(FFT_IS_ALIGNED(sub, FFT_FLOATS_PER_SUBMATRIX));
	RFASSERT(FFT_IS_ALIGNED(rowSize, FFT_FLOATS_PER_SUBMATRIX));
	
	size_t colIncr = FFT_FLOATS_PER_VECTOR * rowSize;
	
	for(unsigned row=0; row<FFT_FLOATS_PER_SUBMATRIX; row++) {
		/* 
		 * start row pointer at [row, row+1]
		 * start column pointer at [row+1, row]
		 */
		FFTFloat *rowPtr = sub + (row * rowSize) + row + 1;
		fft_prefetch(rowPtr);
		/* down one row, back one column */
		FFTFloat *colPtr = rowPtr + rowSize - 1;
		fft_prefetch(colPtr);
		
		size_t col = row + 1;

		/* 
		 * Do some swaps outside of the loop on rows on which the current
		 * rowPtr is not vector-aligned so inner loop processes a vector-sized
		 * number of elements (and works on vector-aligned pointers). 
		 * Note: this is needed even though the *inputs* to this routine
		 * are well aligned (and they always are); this covers the elements
		 * next to the diagonal. 
		 */
		FFTFloat tmp;
		size_t fillerRows = (FFT_FLOATS_PER_VECTOR - col) & FFT_FLOAT_VECTOR_MASK;
		
		for(unsigned dex=0; dex<fillerRows; dex++) {
			tmp       = *rowPtr;
			*rowPtr++ = *colPtr;
			*colPtr   = tmp;
			colPtr   += rowSize;
			col++;
		}
		
		for(; col<FFT_FLOATS_PER_SUBMATRIX; col+=FFT_FLOATS_PER_VECTOR) {
			swapRowCol(rowPtr, colPtr, rowSize);
			rowPtr += FFT_FLOATS_PER_VECTOR;
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
	FFTFloat *sub1,
	FFTFloat *sub2,
	size_t rowSize)			/* offset between rows, in FFTComplex */
{
	dumpSub("fftSubSwap start sub1", sub1, rowSize);
	dumpSub("fftSubSwap start sub2", sub2, rowSize);
	RFASSERT(FFT_IS_ALIGNED(sub1, FFT_FLOATS_PER_SUBMATRIX));
	RFASSERT(FFT_IS_ALIGNED(sub2, FFT_FLOATS_PER_SUBMATRIX));
	RFASSERT(FFT_IS_ALIGNED(rowSize, FFT_FLOATS_PER_SUBMATRIX));

	unsigned colIncr = FFT_FLOATS_PER_VECTOR * rowSize;

	/* row and col refer to coordinates in sub1 */
	for(unsigned row=0; row<FFT_FLOATS_PER_SUBMATRIX; row++) {
		FFTFloat *row1Ptr = sub1 + (row * rowSize);
		fft_prefetch(row1Ptr);
		FFTFloat *col2Ptr = sub2 + row;

		fft_prefetch(col2Ptr);
		for(size_t col=0; col<FFT_FLOATS_PER_SUBMATRIX; col+=FFT_FLOATS_PER_VECTOR) {
			swapRowCol(row1Ptr, col2Ptr, rowSize);
			row1Ptr += FFT_FLOATS_PER_VECTOR;
			col2Ptr += colIncr;
		}
	}
	dumpSub("fftSubSwap end sub1", sub1, rowSize);
	dumpSub("fftSubSwap end sub2", sub2, rowSize);
}

/*
 * Basic top-level row/column swapper for a single square array of FFTFloats.
 * Size of the array and the portion to transform must be power-of-two multiples
 * of FFT_FLOATS_PER_SUBMATRIX.
 */
static MFFTReturn fftSwapRowsColumnsOpt(
	FFTFloat		*buf,
	size_t			rowSize,		/* offset between rows, in FFTComplex */
	/* portion to move - must be FFT_FLOATS_PER_SUBMATRIX aligned */
	size_t			startRow,
	size_t			rowsToMove)
{
	size_t numRowSubs = rowsToMove / FFT_FLOATS_PER_SUBMATRIX;
	size_t numColSubs = rowSize / FFT_FLOATS_PER_SUBMATRIX;
	size_t row = startRow;
	
	tpThreadDebug("fftSwapRowsColumnsOpt top: startRow %lu rowsToMove %lu\n", 
		(unsigned long)startRow, (unsigned long)rowsToMove);
	RFASSERT((numRowSubs * FFT_FLOATS_PER_SUBMATRIX) == rowsToMove);
	RFASSERT((numColSubs * FFT_FLOATS_PER_SUBMATRIX) == rowSize);
	
	/*
	 * -- base is the upper left corner of the submatrix on the diagonal we start 
	 *    with at the top of each outer loop.
	 * -- sub1 is the upper left corner of submatrices to the right of base.
	 * -- sub2 is the upper left corner of submatrices below base.
	 *
	 * All pointer arithmetic is in FFTComplex.
	 * rowSubDex and colSubDex are indices into an array of submatrices.
	 */
	size_t startRowSub = startRow / FFT_FLOATS_PER_SUBMATRIX;
	RFASSERT((startRowSub * FFT_FLOATS_PER_SUBMATRIX) == startRow);
	size_t endRowSub = startRowSub + numRowSubs;
	size_t sub2Incr = FFT_FLOATS_PER_SUBMATRIX * rowSize;

	/* this points to the submatrix on the diagonal */
	FFTFloat *base = buf + (row * rowSize) + row;
	
	/* increment between successive base ptrs */
	size_t baseIncr = FFT_FLOATS_PER_SUBMATRIX * (rowSize + 1);
	
	for(size_t rowSubDex=startRowSub; 
	    rowSubDex<endRowSub; 
		rowSubDex++, row+=FFT_FLOATS_PER_SUBMATRIX) {
		
		/* swap the submatrix on the diagonal in place */
		fftSubTransInPlace(base, rowSize);
		
		/* now swap between [x,y] and [y,x] */
		/* sub1 is 'n' submatrices to the right of base */
		FFTFloat *sub1 = base + FFT_FLOATS_PER_SUBMATRIX;
		/* sub2 is 'n' submatrices down from base */
		FFTFloat *sub2 = base + sub2Incr;
		for(size_t colSubDex=rowSubDex+1; colSubDex<numColSubs; colSubDex++) { 
			fftSubSwap(sub1, sub2, rowSize);
			sub1 += FFT_FLOATS_PER_SUBMATRIX;
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
	TP_TransposeSplit *ttop  = &u->transS;
	fftSwapRowsColumnsOpt(ttop->dst, ttop->numRows, ttop->startRow, ttop->rowsToMove);
	return MR_Success;
}

#pragma mark --- Threaded square transform ---

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
	FFTFloat	*buf,
	size_t		numRows)			/* offset between rows, in FFTComplex */
{
	tpThreadDebug("fftSwapRowsColumnsThr top: numRows %lu\n", (unsigned long)numRows);

	size_t numRowSubs = numRows / FFT_FLOATS_PER_SUBMATRIX;
	RFASSERT((numRowSubs * FFT_FLOATS_PER_SUBMATRIX) == numRows);
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
	
	/* set up a TP_TransposeSplit for each thread */	
	for(unsigned dex=0; dex<numThreads; dex++) {
		TP_PerThread *pt   = &mfftPlan->threadPoolP->perThread[dex];
		TP_Task *task      = &pt->task;
		TP_TransposeSplit *ttop = &task->u->transS;

		task->op         = TPO_App;
		task->threadFcn  = swapRowsColsThr;
		ttop->src        = NULL;		/* not used, make sure of it */
		ttop->dst        = buf;
		ttop->numRows    = numRows;
		ttop->numCols    = 0;		/* not used */
		ttop->startRow   = startSub[dex] * FFT_FLOATS_PER_SUBMATRIX;
		ttop->rowsToMove = numSubs[dex] * FFT_FLOATS_PER_SUBMATRIX;
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

/* Brute force manual implementation */
static MFFTReturn fftTransposeSquareSmall(
	FFTFloat			*buf,
	size_t				rowSize)		/* size of each side in complex elements */
{
	for(size_t row=0; row<rowSize; row++) {
		/* start rowPtr at (row, row+1) */
		FFTFloat *rowPtr = buf + (row * rowSize) + row + 1; 
		/* start colPtr one row down, one column back */
		FFTFloat *colPtr = rowPtr + rowSize - 1;
		
		for(size_t col=row+1; col<rowSize; col++) {
			FFTFloat tmp = *rowPtr;
			*rowPtr++    = *colPtr;
			*colPtr      = tmp;
			colPtr      += rowSize;
		}
	}
	return MR_Success;
}

#pragma mark --- Public API ---

/* minimum sumbtarices for threading */
#define THREAD_MIN_SUBS							2

/* Force brute-force for debug */
#define SQUARE_TRANSPOSE_BASIC					0

MFFTReturn fftTransposeSquare(
	MatrixFFTPlan		mfftPlan,
	FFTComplex			*buf,
	size_t				numRows)		/* size of each side in complex elements */
{
	if(!FFT_INTEL || SQUARE_TRANSPOSE_BASIC) {
		/* For now, use brute force */
		fftTransposeSquareSmall(buf->real, numRows);
		return fftTransposeSquareSmall(buf->imag, numRows);
	}
	
	if(numRows < FFT_FLOATS_PER_SUBMATRIX) {
		/* too small to use submatrices */
		fftTransposeSquareSmall(buf->real, numRows);
		return fftTransposeSquareSmall(buf->imag, numRows);
	}
	
	size_t numSubs = numRows / FFT_FLOATS_PER_SUBMATRIX;
	MFFTReturn mrtn;
	if((mfftPlan->threadPoolP->numThreads > 0) && (numSubs >= THREAD_MIN_SUBS)) {
		mrtn = fftSwapRowsColumnsThr(mfftPlan, buf->imag, numRows);
		if(mrtn) {
			return mrtn;
		}
		return fftSwapRowsColumnsThr(mfftPlan, buf->real, numRows);
	}

	/* single threaded, submatrix optimization */
	mrtn = fftSwapRowsColumnsOpt(buf->real, numRows, 0, numRows);
	if(mrtn) {
		return mrtn;
	}
	return fftSwapRowsColumnsOpt(buf->imag, numRows, 0, numRows);
}

#endif	/* FFT_SPLIT_COMPLEX */
