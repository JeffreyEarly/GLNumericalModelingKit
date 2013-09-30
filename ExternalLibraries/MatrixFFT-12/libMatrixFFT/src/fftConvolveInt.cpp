/*	File: fftConvolveSplit.cpp
	
	Description:
		MatrixFFT-based convolution: interleaved-complex routines
	
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
 * This is an implementation of Algortihm 4 in "Large-scale FFTs and convolutions on 
 * Apple hardware", available here:
 * http://images.apple.com/acg/pdf/FFTapps.pdf
 *
 * Created 12/23/2008. 
 * Copyright 2008 by Apple, Inc. 
 */
 
#include <libMatrixFFT/fftConvolve.h>
#include <strings.h>
#include <cstdio>
#include <libMatrixFFT/complexBufUtils.h>
#include "fftPriv.h"
#include <stdio.h>

#define FFT_CONV_DUMP		0
#if		FFT_CONV_DUMP
#define dumpConvFloat(title, buf, r, c)			fftDumpBuf(title, buf, r, c)
#define dumpConvSplit(title, buf, r, c)			fftDumpMatrixReal(title, buf, r, c)
#else
#define dumpConvFloat(title, buf, r, c)
#define dumpConvSplit(title, buf, r, c)
#endif	/* FFT_CONV_DUMP */

#ifdef	DEBUG
#define debugZero(b, n)		genConstComplex(b, n, 0.0)
#else
#define debugZero(b, n)
#endif


#if		!FFT_SPLIT_COMPLEX

/* 
 * Copy a square, odd-sized kernel to a power-of-two-sized FFTComplex.
 *
 * Example, a 5x5 kernel like this
 *
 *       00 01 02 03 04
 *       10 11 12 13 14
 *       20 21 22 23 24
 *       30 31 32 33 34
 *       40 41 42 43 44
 *
 * Gets placed in the FFTComplex like this:
 *
 *       22 23 24   ...    20 21
 *       32 33 34   ...    30 31
 *       42 43 44   ...    40 41
 *         ...      ...     ...
 *       02 03 04   ...    00 01
 *       12 13 14   ...    10 11
 */
void fftConvCopyKernel(
	const FFTFloat			*kernel,
	FFTComplex				*dst,
	size_t					kernelWidth,
	unsigned				log2dstNumRows,
	unsigned				log2dstNumCols,
	bool					doZero)
{
	/* destination dimensions in floats */
   	size_t dstNumRows = (size_t)1 << log2dstNumRows;
	size_t dstNumCols = (size_t)1 << log2dstNumCols;
	/* in complexes */
	size_t numColsComplex = dstNumCols / 2;
	
	dumpConvFloat("fftConvCopyKernel source", kernel, kernelWidth, kernelWidth);

	/* the number of floats per row in each of {dst->real, dst->imag} */

	#ifdef	DEBUG
	if((kernelWidth > dstNumRows) || (kernelWidth > dstNumCols)) { 
		printf("***fftConvCopyKernel: overflow\n");
		return;
	}
	if((kernelWidth & 0x01) == 0) {
		printf("***fftConvCopyKernel: even size kernel\n");
		return;
	}
	#endif
	
	if(doZero) {
		genConstComplex(dst, numColsComplex * dstNumRows, 0.0);
	}
	
	/* degenerate, trivial (though legal) case */
	if(kernelWidth == 1) {
		dst->real = *kernel;
		dumpConvSplit("fftConvCopyKernel dest", dst, dstNumRows, dstNumCols);
		return;
	}
	
	/* even and odd "half kernel" sizes */
	size_t nkOver2 = kernelWidth / 2;
	size_t nkOver2_p1 = nkOver2 + 1;
	FFTComplex *dstP;
	
	/* top: (kernelWidth / 2) + 1 rows */
	for(size_t row=0; row<nkOver2_p1; row++) {
		size_t offset = row * numColsComplex;
		dstP = dst + offset;
		
		/* 
		 * left: (kernelWidth / 2) + 1 elements, starting from
		 * kernel[nkOver2+row, nkOver2]
		 */
		const FFTFloat *rowStart = kernel + ((nkOver2 + row) * kernelWidth);
		const FFTFloat *src = rowStart + nkOver2;
		for(size_t dex=0; dex<nkOver2_p1; dex+=2) {
			dstP->real = *src++;
			if(dex == (nkOver2_p1 - 1)) {
				break;
			}
			dstP->imag = *src++;
			dstP++;
		}
		
		/* 
		 * Right: (kernelWidth / 2) elements, starting from
		 * kernel[nkOver2+row, 0]
		 * Go backwards, it's easier...
		 * start from kernel[nkOver2+row, nkOver2-1]
		 */
		src = rowStart + nkOver2 - 1;
		offset = ((row + 1) * numColsComplex) - 1;
		dstP = dst + offset;
		for(size_t dex=0; dex<nkOver2; dex+=2) {
			dstP->imag = *src--;
			if(dex == (nkOver2 - 1)) {
				break;
			}
			dstP->real = *src--;
			dstP--;
		}
	}
	
	/* bottom: (kernelWidth / 2) rows */
	size_t dstRow = dstNumRows - nkOver2;
	for(size_t kernRow=0; kernRow<nkOver2; kernRow++, dstRow++) {
		size_t offset = dstRow * numColsComplex;
		dstP = dst + offset;
		
		/* 
		 * left: (kernelWidth / 2) + 1 elements, starting from
		 * kernel[kernRow, nkOver2]
		 */
		const FFTFloat *rowStart = kernel + (kernRow * kernelWidth);
		const FFTFloat *src = rowStart + nkOver2;
		for(size_t dex=0; dex<nkOver2_p1; dex+=2) {
			dstP->real = *src++;
			if(dex == (nkOver2_p1 - 1)) {
				break;
			}
			dstP->imag = *src++;
			dstP++;
		}
		
		/* 
		 * Right: (kernelWidth / 2) elements, starting from
		 * kernel[nkOver2_p1+kernRow, 0]
		 * Go backwards, it's easier...
		 * start from kernel[kernRow, nkOver2-1]
		 */
		src = rowStart + nkOver2 - 1;
		offset = ((dstRow + 1) * numColsComplex) - 1;
		dstP = dst + offset;
		for(size_t dex=0; dex<nkOver2; dex+=2) {
			dstP->imag = *src--;
			if(dex == (nkOver2 - 1)) {
				break;
			}
			dstP->real = *src--;
			dstP--;
		}
	}
	dumpConvSplit("fftConvCopyKernel dest", dst, dstNumRows, dstNumCols);
}

/* 
 * If the destination for fftConvCopyImage() is this size or smaller,
 * we just bzero the whole dst buffer on entry instead of doing it
 * line-by-line in the main loop.
 */
#define FCT_DST_SIZE_THRESHHOLD		(8 * 1024)

/*
 * Copy an arbitrarily sized image to a power-of-two-sized FFTComplex.
 */
void fftConvCopyImage(
	const FFTFloat			*image,
	FFTComplex				*dst,
	size_t					imageNumRows,
	size_t					imageNumCols,
	unsigned				log2dstNumRows,
	unsigned				log2dstNumCols,
	bool					doZero)
{
	/* destination dimensions in floats */
   	size_t dstNumRows = (size_t)1 << log2dstNumRows;
	size_t dstNumCols = (size_t)1 << log2dstNumCols;
	/* in complexes */
	size_t numColsComplex = dstNumCols / 2;
	
	dumpConvFloat("fftConvCopyImage source", image, imageNumRows, imageNumCols);	

	#ifdef	DEBUG
	if((imageNumRows > dstNumRows) || (imageNumCols > dstNumCols)) {
		printf("***fftConvCopyImage: overflow\n");
		return;
	}
	#endif
	
	bool isZeroed = !doZero;
	size_t dstTotalFloats = dstNumRows * dstNumCols;
	if(doZero && (dstTotalFloats <= FCT_DST_SIZE_THRESHHOLD)) {
		genConstComplex(dst, dstTotalFloats >> 1, 0.0);
		isZeroed = true;
	}
	
	/* 
	 * number of (full) complexes to zero at end of row in each of dst.
	 * We handle odd-size image rows inside the loop.
	 */
	size_t remLengthComplex = (dstNumCols - imageNumCols) / 2;
	const FFTFloat *src = image;
	FFTComplex *dstP;
	
	for(size_t row=0; row<imageNumRows; row++) {
		size_t offset = row * numColsComplex;
		dstP = dst + offset;
		
		for(size_t col=0; col<imageNumCols; col+=2) {
			dstP->real = *src++;
			if(col == (imageNumCols - 1)) {
				/* odd size source, zero last imaginary */
				dstP->imag = 0.0;
				dstP++;
				break;
			}
			dstP->imag = *src++;
			dstP++;
		}
		if(!isZeroed  && (remLengthComplex != 0)) {
			genConstComplex(dstP, remLengthComplex, 0.0);
		}
	}
	
	/* now the last remaining rows */
	if(!isZeroed && (dstNumRows != imageNumRows)) {
		size_t offset = imageNumRows * numColsComplex;
		dstP = dst + offset;
		remLengthComplex = (dstNumRows - imageNumRows) * numColsComplex;
		genConstComplex(dstP, remLengthComplex, 0.0);
	}
	dumpConvSplit("fftConvCopyImage dest", dst, dstNumRows, dstNumCols);
}

/*
 * Copy an image from a power-of-two-sized FFTComplex.
 */
void fftConvRetrieveImage(
	const FFTComplex		*image,
	FFTComplex				*dst,
	unsigned				log2SrcNumRows,
	unsigned				log2SrcNumCols,
	size_t					dstNumRows,
	size_t					dstNumCols)
{
	/* 
	 * row sizes in complex.
	 * Round up the destination size - OK? We'll be copying beyond the technical
	 * end of the destination row...
	 */
	size_t srcRowSize = (size_t)1 << (log2SrcNumCols - 1);
	size_t dstRowSize = (dstNumCols + 1) >> 1;
	
	for(size_t row=0; row<dstNumRows; row++) {
		fftCopyComplexArray(image, dst, dstRowSize);
		image += srcRowSize;
		dst   += dstRowSize;
	}
}

/* 
 * Multiply two complexes c1 and c2; result goes to co.
 */
static inline void fftConvMulComplex(
	const FFTComplex *c1,
	const FFTComplex *c2,
	FFTComplex *co)
{
	co->real = (c1->real * c2->real) - (c1->imag * c2->imag);
	co->imag = (c1->real * c2->imag) + (c1->imag * c2->real);
}

/* 
 * Multiply two complexes c1 and c2 using the complex conjugate
 * of c2.
 */
static inline void fftConvMulComplexConj(
	const FFTComplex *c1,
	const FFTComplex *c2,
	FFTComplex *co)
{
	co->real = (c1->real * c2->real) + (c1->imag * c2->imag);
	co->imag = (c1->imag * c2->real) - (c1->real * c2->imag);
}

/* 
 * Dyadic mul of the outputs of a 1-D FFT.  
 * The first two elements are reals, the remainder are complex. 
 * If conj is true we use the complex conjugate of c2.
 */
static inline void fftConvMul1DFFT(
	const FFTComplex	*c1,
	const FFTComplex	*c2,
	FFTComplex			*co,
	bool				conj,
	size_t				n)				// row size in complex elements
{
	co->real = c1->real * c2->real;		/* real DC */
	co->imag = c1->imag * c2->imag;		/* real Nyquist */
	c1++;
	c2++;
	co++;
	
	/* leave conj test out of the loop... */
	if(conj) {
		for(size_t dex=1; dex<n; dex++) {
			fftConvMulComplexConj(c1++, c2++, co++);
		}
	}
	else {
		for(size_t dex=1; dex<n; dex++) {
			fftConvMulComplex(c1++, c2++, co++);
		}
	}
}

/* 
 * Perform dyadic multiply of the output of two forward 2-D real FFT ops.
 * Result goes to dst; inputs are not modified. 
 * If conj is true we use the complex conjugates of src2.
 */
MFFTReturn fftConvDyadicMulCom(
	const FFTComplex	*src1,
	const FFTComplex	*src2,
	FFTComplex			*dst,
	unsigned			log2NumRows,
	unsigned			log2NumCols,
	bool				conj)
{
	/*
	 * Sizes in reals. 
	 */
	size_t numRows = (size_t)1 << log2NumRows;
	size_t numCols = (size_t)1 << log2NumCols;

	/* 
	 * Sizes in complexes.
	 */
	size_t numColsComplex = numCols / 2;
	debugZero(dst, numRows * numColsComplex);

	/* 
	 * The first columns contain the result of 1-D real-signal FFT on the 
	 * DC and Nyquist values from the first row-by-row FFT. 
	 *
	 * Here's what we're going to do:
	 *
	 * -- transpose that first column to auxBufN
	 * -- perform requisite computation on auxBufN
	 * -- transpose auxBuf result back to column 0
	 * -- perform normal complex mltiply on remainder of each row
	 */
	 
	/* 
	 * 1. alloc 3 aux buffers for operating on column 0.
	 *    These don't have to be well aligned, there is no vectorized code
	 *    that will use these. 
	 */
	size_t dcNyqComplexes = numRows / 2;
	FFTComplex *auxBuf1 = fftAllocComplexArray(numRows);
	if(auxBuf1 == NULL) {
		return MR_Memory;
	}
	FFTComplex *auxBuf2 = fftAllocComplexArray(numRows);
	if(auxBuf2 == NULL) {
		return MR_Memory;
	}
	FFTComplex *auxBufDst = fftAllocComplexArray(numRows);
	if(auxBufDst == NULL) {
		return MR_Memory;
	}
	
	/* 
	 * 2. Transpose and multiply the 1-D FFT outputs from column 0. 
	 */ 
	fftCopyFromColumn0(src1, auxBuf1, numColsComplex, dcNyqComplexes);
	fftCopyFromColumn0(src2, auxBuf2, numColsComplex, dcNyqComplexes);
	fftConvMul1DFFT(auxBuf1, auxBuf2, auxBufDst, conj, dcNyqComplexes);
	fftConvMul1DFFT(auxBuf1+dcNyqComplexes, auxBuf2+dcNyqComplexes, 
		auxBufDst+dcNyqComplexes, conj, dcNyqComplexes);
		
	/* 
	 * 3. Transpose result to destination column 0.
	 */
	fftCopyToColumn0(auxBufDst, dst, numColsComplex, dcNyqComplexes);

	/* 
	 * 4. Normal complex multiply of remainder.
	 *    We don't have anything like vDSP_zvmul(), which operates on split complex
	 *    only, so we use brute force. 
	 */
	const FFTComplex *c1 = src1 + 1;
	const FFTComplex *c2 = src2 + 1;
	FFTComplex *co = dst + 1;

	if(conj) {
		for(size_t row=0; row<numRows; row++) {
			for(size_t col=1; col<numColsComplex; col++) {
				fftConvMulComplexConj(c1++, c2++, co++);
			}
			/* pointers now at start of next row; skip first element */
			c1++;
			c2++;
			co++;
		}
	}
	else {
		for(size_t row=0; row<numRows; row++) {
			for(size_t col=1; col<numColsComplex; col++) {
				fftConvMulComplex(c1++, c2++, co++);
			}
			/* pointers now at start of next row; skip first element */
			c1++;
			c2++;
			co++;
		}
	}
	
	fftFreeComplexArray(auxBuf1);
	fftFreeComplexArray(auxBuf2);
	fftFreeComplexArray(auxBufDst);
	
	return MR_Success;
}

/* 
 * Perform dyadic multiply of the output of two forward 2-D real FFT ops
 * when the FFT output is in MF_CustomOrder.
 * Result goes to dst; inputs are not modified. 
 * If conj is true we use the complex conjugates of src2.
 */
MFFTReturn fftConvDyadicMulCustomCom(
	const FFTComplex	*src1,
	const FFTComplex	*src2,
	FFTComplex			*dst,
	unsigned			log2NumRows,
	unsigned			log2NumCols,
	bool				conj)
{
	/*
	 * Sizes in reals. 
	 */
	size_t numRows = (size_t)1 << log2NumRows;
	size_t numCols = (size_t)1 << log2NumCols;

	/* 
	 * Sizes in complexes.
	 */
	size_t numColsComplex = numCols / 2;
	debugZero(dst, numRows * numColsComplex);

	/*
	 * This is actually quite a bit simpler than the version which operates on 
	 * row-order data. In this version the first *row* of src1 and src2 contain
	 * the real-signal FFT on the real and Nyquist values from the first
	 * row-by-row real-signal FFT. No transpose or temp storage needed;
	 * we just use fftConvMul1DFFT() to process those in-place and do 
	 * the rest of the array in bulk.
	 */
	if(log2NumRows != (log2NumCols - 1)) {
		printf("***fftConvDyadicMulCustomCom: not a square signal\n");
		return MR_IllegalArg;
	}
	
	size_t dcNyqComplexes = numColsComplex / 2;
	
	/* 
	 * 1. Multiply the 1-D FFT outputs from row 0. 
	 */ 
	fftConvMul1DFFT(src1, src2, dst, conj, dcNyqComplexes);
	fftConvMul1DFFT(src1+dcNyqComplexes, src2+dcNyqComplexes, dst+dcNyqComplexes, conj, dcNyqComplexes);
		
	/* 
	 * 4. Normal complex multiply of remainder.
	 *    We don't have anything like vDSP_zvmul(), which operates on split complex
	 *    only, so we use brute force. 
	 */
	const FFTComplex *c1 = src1 + numColsComplex;
	const FFTComplex *c2 = src2 + numColsComplex;
	FFTComplex       *co = dst  + numColsComplex;

	if(conj) {
		for(size_t row=1; row<numRows; row++) {
			for(size_t col=0; col<numColsComplex; col++) {
				fftConvMulComplexConj(c1++, c2++, co++);
			}
		}
	}
	else {
		for(size_t row=1; row<numRows; row++) {
			for(size_t col=0; col<numColsComplex; col++) {
				fftConvMulComplex(c1++, c2++, co++);
			}
		}
	}
	
	return MR_Success;
}


MFFTReturn fftConvDyadicMul(
	const FFTComplex		*src1,
	const FFTComplex		*src2,
	FFTComplex				*dst,
	unsigned				log2NumRows,
	unsigned				log2NumCols,
	MFFTFormat				format)
{
	switch(format) {
		case MF_RowOrder:
			return fftConvDyadicMulCom(src1, src2, dst, log2NumRows, log2NumCols, false);
		case MF_CustomOrder:
			return fftConvDyadicMulCustomCom(src1, src2, dst, log2NumRows, log2NumCols, false);
		default:
			return MR_IllegalArg;
	}
}
	
MFFTReturn fftConvDyadicMulConj(
	const FFTComplex		*src1,
	const FFTComplex		*src2,
	FFTComplex				*dst,
	unsigned				log2NumRows,
	unsigned				log2NumCols,
	MFFTFormat				format)
{
	switch(format) {
		case MF_RowOrder:
			return fftConvDyadicMulCom(src1, src2, dst, log2NumRows, log2NumCols, true);
		case MF_CustomOrder:
			return fftConvDyadicMulCustomCom(src1, src2, dst, log2NumRows, log2NumCols, true);
		default:
			return MR_IllegalArg;
	}
}


#endif	/* !FFT_SPLIT_COMPLEX */
