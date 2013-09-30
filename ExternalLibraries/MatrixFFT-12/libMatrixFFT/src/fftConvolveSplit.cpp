/*	File: fftConvolveSplit.cpp
	
	Description:
		MatrixFFT-based convolution: common API and split-complex routines
	
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
#include <libMatrixFFT/complexBufUtils.h>
#include "fftPriv.h"

/*
 * Note: in this module (as elsewhere in the MatrixFFT library), the origin of 
 * 2-D data is in the upper left corner when we talk about up/down/left/right.
 */
 
#define FFT_CONV_DUMP		0
#if		FFT_CONV_DUMP
#define dumpConvFloat(title, buf, r, c)			fftDumpBuf(title, buf, r, c)
#define dumpConvSplit(title, buf, r, c)			fftDumpMatrixReal(title, buf, r, c)
#else
#define dumpConvFloat(title, buf, r, c)
#define dumpConvSplit(title, buf, r, c)
#endif	/* FFT_CONV_DUMP */

#if		FFT_SPLIT_COMPLEX

#ifdef	DEBUG
#define debugZero(b, n)		genConstComplex(b, n, 0.0)
#else
#define debugZero(b, n)
#endif

/* Indirection for precision-independent vDSP_zvmul*() */
#if		FFT_DOUBLE_PREC

/* normal complex vector mul */
#define FFTZvmul(c1, c2, num, dst)		vDSP_zvmulD(c1, 1, c2, 1, dst, 1, num, 1)

/* complex conjugate vector mul */
#define FFTZvmulConj(c1, c2, num, dst)	vDSP_zvcmulD(c2, 1, c1, 1, dst, 1, num)

#else	/* !FFT_DOUBLE_PREC */

#define FFTZvmul(c1, c2, num, dst)		vDSP_zvmul(c1, 1, c2, 1, dst, 1, num, 1)
#define FFTZvmulConj(c1, c2, num, dst)	vDSP_zvcmul(c2, 1, c1, 1, dst, 1, num)

#endif	/* FFT_DOUBLE_PREC */

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
   	size_t dstNumRows = (size_t)1 << log2dstNumRows;
	size_t dstNumCols = (size_t)1 << log2dstNumCols;
	
	dumpConvFloat("fftConvCopyKernel source", kernel, kernelWidth, kernelWidth);

	/* the number of floats per row in each of {dst->real, dst->imag} */
	size_t dstNumRowFloats = dstNumCols / 2;

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
		size_t bufSize = dstNumRowFloats * dstNumRows * sizeof(FFTFloat);
		memset(dst->real, 0, bufSize);
		memset(dst->imag, 0, bufSize);
	}
	
	/* degenerate, trivial (though legal) case */
	if(kernelWidth == 1) {
		*dst->real = *kernel;
		dumpConvSplit("fftConvCopyKernel dest", dst, dstNumRows, dstNumCols);
		return;
	}
	
	/* even and odd "half kernel" sizes */
	size_t nkOver2 = kernelWidth / 2;
	size_t nkOver2_p1 = nkOver2 + 1;
	
	/* top: (kernelWidth / 2) + 1 rows */
	for(size_t row=0; row<nkOver2_p1; row++) {
		size_t offset = row * dstNumRowFloats;
		FFTFloat *realp = dst->real + offset;
		FFTFloat *imagp = dst->imag + offset;
		
		/* 
		 * left: (kernelWidth / 2) + 1 elements, starting from
		 * kernel[nkOver2+row, nkOver2]
		 */
		const FFTFloat *rowStart = kernel + ((nkOver2 + row) * kernelWidth);
		const FFTFloat *src = rowStart + nkOver2;
		for(size_t dex=0; dex<nkOver2_p1; dex+=2) {
			*realp++ = *src++;
			if(dex == (nkOver2_p1 - 1)) {
				break;
			}
			*imagp++ = *src++;
		}
		
		/* 
		 * Right: (kernelWidth / 2) elements, starting from
		 * kernel[nkOver2+row, 0]
		 * Go backwards, it's easier...
		 * start from kernel[nkOver2+row, nkOver2-1]
		 */
		src = rowStart + nkOver2 - 1;
		offset = ((row + 1) * dstNumRowFloats) - 1;
		realp = dst->real + offset;
		imagp = dst->imag + offset;
		for(size_t dex=0; dex<nkOver2; dex+=2) {
			*imagp-- = *src--;
			if(dex == (nkOver2 - 1)) {
				break;
			}
			*realp-- = *src--;
		}
	}
	
	/* bottom: (kernelWidth / 2) rows */
	size_t dstRow = dstNumRows - nkOver2;
	for(size_t kernRow=0; kernRow<nkOver2; kernRow++, dstRow++) {
		size_t offset = dstRow * dstNumRowFloats;
		FFTFloat *realp = dst->real + offset;
		FFTFloat *imagp = dst->imag + offset;
		
		/* 
		 * left: (kernelWidth / 2) + 1 elements, starting from
		 * kernel[kernRow, nkOver2]
		 */
		const FFTFloat *rowStart = kernel + (kernRow * kernelWidth);
		const FFTFloat *src = rowStart + nkOver2;
		for(size_t dex=0; dex<nkOver2_p1; dex+=2) {
			*realp++ = *src++;
			if(dex == (nkOver2_p1 - 1)) {
				break;
			}
			*imagp++ = *src++;
		}
		
		/* 
		 * Right: (kernelWidth / 2) elements, starting from
		 * kernel[nkOver2_p1+kernRow, 0]
		 * Go backwards, it's easier...
		 * start from kernel[kernRow, nkOver2-1]
		 */
		src = rowStart + nkOver2 - 1;
		offset = ((dstRow + 1) * dstNumRowFloats) - 1;
		realp = dst->real + offset;
		imagp = dst->imag + offset;
		for(size_t dex=0; dex<nkOver2; dex+=2) {
			*imagp-- = *src--;
			if(dex == (nkOver2 - 1)) {
				break;
			}
			*realp-- = *src--;
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
	size_t dstNumRows = (size_t)1 << log2dstNumRows;
	size_t dstNumCols = (size_t)1 << log2dstNumCols;
	
	dumpConvFloat("fftConvCopyImage source", image, imageNumRows, imageNumCols);

	/* the number of floats per row in each of {dst->realp, dst->imagp} */
	size_t dstNumRowFloats = dstNumCols / 2;
	
	#ifdef	DEBUG
	if((imageNumRows > dstNumRows) || (imageNumCols > dstNumCols)) {
		printf("***fftConvCopyImage: overflow\n");
		return;
	}
	#endif
	
	bool isZeroed = !doZero;
	if(doZero && ((dstNumRows * dstNumCols) <= FCT_DST_SIZE_THRESHHOLD)) {
		size_t bufSize = dstNumRowFloats * dstNumRows * sizeof(FFTFloat);
		memset(dst->real, 0, bufSize);
		memset(dst->imag, 0, bufSize);
		isZeroed = true;
	}
	
	/* number of floats to zero at end of row in each of {dst->realp, dst->imagp} */
	size_t remLengthFloats = (dstNumCols - imageNumCols) / 2;
	size_t remLengthBytes = remLengthFloats * sizeof(FFTFloat);
	const FFTFloat *src = image;
	
	for(size_t row=0; row<imageNumRows; row++) {
		size_t offset = row * dstNumRowFloats;
		FFTFloat *realp = dst->real + offset;
		FFTFloat *imagp = dst->imag + offset;
		
		for(size_t col=0; col<imageNumCols; col+=2) {
			*realp++ = *src++;
			if(col == (imageNumCols - 1)) {
				/* odd size source, zero last imaginary */
				*imagp++ = 0.0;
				break;
			}
			*imagp++ = *src++;
		}
		if(!isZeroed) {
			memset(realp, 0, remLengthBytes);
			memset(imagp, 0, remLengthBytes);
		}
	}
	
	/* now the last remaining rows */
	if(!isZeroed && (dstNumRows != imageNumRows)) {
		size_t offset = imageNumRows * dstNumRowFloats;
		FFTFloat *realp = dst->real + offset;
		FFTFloat *imagp = dst->imag + offset;

		remLengthBytes = (dstNumRows - imageNumRows) * dstNumRowFloats * sizeof(FFTFloat);
		memset(realp, 0, remLengthBytes);
		memset(imagp, 0, remLengthBytes);
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
	const FFTFloat *srcReal = image->real;
	const FFTFloat *srcImag = image->imag;
	FFTFloat *dstReal = dst->real; 
	FFTFloat *dstImag = dst->imag;
	
	size_t srcNumCols = (size_t)1 << log2SrcNumCols;
	size_t srcRowFloats = srcNumCols / 2;
	size_t toMove = dstNumCols * sizeof(FFTFloat);
	size_t dstRowFloats = dstNumCols / 2;
	
	for(size_t dex=0; dex<dstNumRows; dex++) {
		memmove(dstReal, srcReal, toMove);
		memmove(dstImag, srcImag, toMove);
		
		srcReal += srcRowFloats;
		srcImag += srcRowFloats;
		dstReal += dstRowFloats;
		dstImag += dstRowFloats;
	}
}
	
/* 
 * Multiply two complexes {r1, i1} and {r2, i2}
 * Result goes to {oR, oI}
 */
static inline void fftConvMulComplex(
	const FFTFloat r1,
	const FFTFloat i1,
	const FFTFloat r2,
	const FFTFloat i2,
	FFTFloat *oR,
	FFTFloat *oI)
{
	*oR = (r1 * r2) - (i1 * i2);
	*oI = (r1 * i2) + (i1 * r2);
}

/* 
 * Multiply two complexes {r1, i1} and {r2, i2} using the complex conjugate
 * of {r2, i2}.
 * Result goes to {oR, oI}
 */
static inline void fftConvMulComplexConj(
	const FFTFloat r1,
	const FFTFloat i1,
	const FFTFloat r2,
	const FFTFloat i2,
	FFTFloat *oR,
	FFTFloat *oI)
{
	*oR = (r1 * r2) + (i1 * i2);
	*oI = (i1 * r2) - (r1 * i2);
}


/* 
 * Dyadic mul of the outputs of a 1-D real FFT. The first two elements
 * are reals, the remainder are complex. 
 * If conj is true we use the complex conjugate of {real2, imag2}.
 */
static inline void fftConvMul1DFFT(
	const FFTFloat	*real1,
	const FFTFloat	*imag1,
	const FFTFloat	*real2,
	const FFTFloat	*imag2,
	FFTFloat		*dstReal,
	FFTFloat		*dstImag,
	bool			conj,
	size_t			n)		/* number of real elements in each of real* and imag* */
{
	*dstReal++ = *real1++ * *real2++;	/* real DC */
	*dstImag++ = *imag1++ * *imag2++;	/* real Nyquist */
	
	/* leave conj test out of the loop... */
	if(conj) {
		for(size_t dex=1; dex<n; dex++) {
			fftConvMulComplexConj(*real1++, *imag1++, *real2++, *imag2++, dstReal++, dstImag++);
		}
	}
	else {
		for(size_t dex=1; dex<n; dex++) {
			fftConvMulComplex(*real1++, *imag1++, *real2++, *imag2++, dstReal++, dstImag++);
		}
	}
}

/* 
 * When true, do the bulk of the complex multiplies with vDSP_zvmul*() 
 */
#define COMPLEX_MUL_VIA_VDSP	1

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
	size_t totalComplex = numRows * numColsComplex;
	debugZero(dst, totalComplex);

	/* 
	 * The first columns contain the result of 1-D real-signal FFT on the 
	 * DC and Nyquist values from the first row-by-row FFT. 
	 *
	 * Here's what we're going to do:
	 *
	 * -- transpose that first column to auxBufN
	 * -- perform requisite computation on auxBufN
	 * -- blindly do normal complex multiply on the entire arrays, including
	 *    the first columns (which is bogus but it'll be faster than skipping
	 *    it)
	 * -- transpose auxBufN back to column 0
	 */
	 
	/* 
	 * 1. alloc 3 aux buffers for operating on column 0.
	 *    No alignment needed, as no vectorized code is going to work on these. 
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
	
	FFTFloat *real1 = auxBuf1->real;
	FFTFloat *imag1 = auxBuf1->imag;
	FFTFloat *real2 = auxBuf2->real;
	FFTFloat *imag2 = auxBuf2->imag;
	FFTFloat *dstReal = auxBufDst->real;
	FFTFloat *dstImag = auxBufDst->imag;

	fftConvMul1DFFT(real1, imag1, real2, imag2, dstReal, dstImag, conj, dcNyqComplexes);
	real1 += dcNyqComplexes;
	imag1 += dcNyqComplexes;
	real2 += dcNyqComplexes;
	imag2 += dcNyqComplexes;
	dstReal += dcNyqComplexes;
	dstImag += dcNyqComplexes;
	fftConvMul1DFFT(real1, imag1, real2, imag2, dstReal, dstImag, conj, dcNyqComplexes);

	/* 
	 * 3. Normal complex multiply of the entire arrays including the bogus
	 *    column 0. 
	 */
	#if COMPLEX_MUL_VIA_VDSP
	vDSPComplex vBuf1, vBuf2, vBufDst;
	fftComplexToVDSP(src1, &vBuf1);
	fftComplexToVDSP(src2, &vBuf2);
	fftComplexToVDSP(dst, &vBufDst);
	if(conj) {
		FFTZvmulConj(&vBuf1, &vBuf2, totalComplex, &vBufDst);
	}
	else {
		FFTZvmul(&vBuf1, &vBuf2, totalComplex, &vBufDst);
	}
	
	#else	/* COMPLEX_MUL_VIA_VDSP */

	real1 = src1->real;
	imag1 = src1->imag;
	real2 = src2->real;
	imag2 = src2->imag;
	dstReal = dst->real;
	dstImag = dst->imag;
	
	if(conj) {
		for(size_t dex=0; dex<totalComplex; dex++) {
			fftConvMulComplexConj(*real1++, *imag1++, *real2++, *imag2++, dstReal++, dstImag++);
		}
	}
	else {
		for(size_t dex=0; dex<totalComplex; dex++) {
			fftConvMulComplex(*real1++, *imag1++, *real2++, *imag2++, dstReal++, dstImag++);
		}
	}
	#endif	/* COMPLEX_MUL_VIA_VDSP */
	
	/* 
	 * 4. Restore the proper column0 values. 
	 */
	fftCopyToColumn0(auxBufDst, dst, numColsComplex, dcNyqComplexes);
	
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
	size_t totalComplex = numRows * numColsComplex;
	debugZero(dst, totalComplex);

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
	size_t dcNyqComplexes = numRows / 2;
	
	/*
	 * 1. Multiply the 1-D FFT outputs from row 0. 
	 */
	FFTFloat *real1 = src1->real;
	FFTFloat *imag1 = src1->imag;
	FFTFloat *real2 = src2->real;
	FFTFloat *imag2 = src2->imag;
	FFTFloat *dstReal = dst->real;
	FFTFloat *dstImag = dst->imag;

	fftConvMul1DFFT(real1, imag1, real2, imag2, dstReal, dstImag, conj, dcNyqComplexes);
	real1 += dcNyqComplexes;
	imag1 += dcNyqComplexes;
	real2 += dcNyqComplexes;
	imag2 += dcNyqComplexes;
	dstReal += dcNyqComplexes;
	dstImag += dcNyqComplexes;
	fftConvMul1DFFT(real1, imag1, real2, imag2, dstReal, dstImag, conj, dcNyqComplexes);


	/* 
	 * 2. Bulk complex multiple of remainder.
	 */
	real1 += dcNyqComplexes;
	imag1 += dcNyqComplexes;
	real2 += dcNyqComplexes;
	imag2 += dcNyqComplexes;
	dstReal += dcNyqComplexes;
	dstImag += dcNyqComplexes;
	size_t remComplexes = totalComplex - numColsComplex;
	
	vDSPComplex vBuf1 = {real1, imag1};
	vDSPComplex vBuf2 = {real2, imag2};
	vDSPComplex vBufDst = {dstReal, dstImag};
	if(conj) {
		FFTZvmulConj(&vBuf1, &vBuf2, remComplexes, &vBufDst);
	}
	else {
		FFTZvmul(&vBuf1, &vBuf2, remComplexes, &vBufDst);
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

#endif	/* FFT_SPLIT_COMPLEX */
 
/* 
 * High-level interface, independent of complex format. 
 *
 * -- 2-D real FFT buf1 in place
 * -- 2-D real FFT buf2 in place
 * -- dyadic multiply the result, normal (conj=false) or conjugate (conj=true)
 * -- inverse 2-D real FFT the result --> dst
 */
MFFTReturn fftConvolve(
	MatrixFFTPlan			mfftPlan,
	FFTComplex				*buf1,
	FFTComplex				*buf2,
	unsigned				log2NumRows,
	unsigned				log2NumCols,
	bool					conj,
	FFTComplex				*dst)
{
	size_t fftRows = (size_t)1 << log2NumRows;
	size_t fftCols = (size_t)1 << log2NumCols;
	
	/* forward FFT buf1, buf2 in place */
	MFFTReturn mrtn = mfftExecute(mfftPlan, 0, true, buf1, buf1);
	if(mrtn) {
		return mrtn;
	}
	dumpConvSplit("fftConvolve FFT(src1)", buf1, fftRows, fftCols);
	
	mrtn = mfftExecute(mfftPlan, 0, true, buf2, buf2);
	if(mrtn) {
		return mrtn;
	}
	dumpConvSplit("fftConvolve FFT(src2)", buf2, fftRows, fftCols);

	/* Obtain format for the dyadic mul */
	MFFTFormat format;
	mrtn = mfftNativeFormat(mfftPlan, false, &format);
	if(mrtn) {
		return mrtn;
	}
	
	if(format == MF_RowOrder) {
		mrtn = fftConvDyadicMulCom(buf1, buf2, dst, log2NumRows, log2NumCols, conj);
	}
	else {
		mrtn = fftConvDyadicMulCustomCom(buf1, buf2, dst, log2NumRows, log2NumCols, conj);
	}
	if(mrtn) {
		return mrtn;
	}
	dumpConvSplit("fftConvDyadicMul result", dst, fftRows, fftCols);
	
	mrtn = mfftExecute(mfftPlan, MEF_NormOutput, false, dst, dst);
	dumpConvSplit("fftConvolve invFFT", dst, fftRows, fftCols);

	/* 
	 * Scale. The output of one forward 2-D real FFT contains values
	 * that are 2x the actual FFT. When executing a reverse FFT 
	 * this 2x is handled via MEF_NormOutput, but we've done 2 forward 
	 * FFTs and multiplied the results together, followed by only one 
	 * reverse FFT, so we have to divide by two once more. 
	 */
	fftScaleComplex(dst, 0.5, (fftRows * fftCols) / 2);
	return mrtn;
}
	
