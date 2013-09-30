/*	File: fftConvolve.h
	
	Description:
		Public interface for MatrixFFT-based convolution.
	
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
 * fftConvolve.h - Public interface for MatrixFFT-based convolution. 
 *
 * Created 12/23/2008. 
 * Copyright 2008 by Apple, Inc. 
 */
 
#ifndef	_FFT_CONVOLVE_H_
#define _FFT_CONVOLVE_H_

#include <Accelerate/Accelerate.h>
#include <stdint.h>
#include <stdbool.h>
#include <libMatrixFFT/MatrixFFT.h>

#ifdef __cplusplus
extern "C" {
#endif

/* 
 * Copy a square, odd-sized kernel to a power-of-two-sized FFTComplex.
 * Destination buffer must be zeroed somewhere; we'll do the zero here 
 * (with a massive memset) if doZero is true. 
 */
extern void fftConvCopyKernel(
	const FFTFloat			*kernel,
	FFTComplex				*dst,
	size_t					kernelWidth,
	unsigned				log2dstNumRows,
	unsigned				log2dstNumCols,
	bool					doZero);
	
/*
 * Copy an arbitrarily sized real-signal image to a power-of-two-sized 
 * FFTSplitComplex. Destination buffer must be zeroed somewhere; we'll 
 * do the zero here if doZero is true. Doing the zero here is more 
 * efficient than a big memset() if the image size is close to the
 * destination size. 
 */
extern void fftConvCopyImage(
	const FFTFloat			*image,
	FFTComplex				*dst,
	size_t					imageNumRows,
	size_t					imageNumCols,
	unsigned				log2dstNumRows,
	unsigned				log2dstNumCols,
	bool					doZero);

/*
 * Copy an image from a power-of-two-sized FFTComplex to a packed
 * split buffer. Used to recover the output of the reverse FFT at the 
 * end of a convolution. 
 */
extern void fftConvRetrieveImage(
	const FFTComplex		*image,
	FFTComplex				*dst,
	unsigned				log2SrcNumRows,
	unsigned				log2SrcNumCols,
	size_t					dstNumRows,
	size_t					dstNumCols);
	
/* 
 * Perform dyadic multiply of the output of two 2-D real FFT ops.
 * For each complex value 0..N-1:
 *
 *    dst[n] = src1[n] * src2[n] 
 *
 * Result goes to dst; inputs are not modified. 
 * The log2NumRows and log2NumCols values match the n[0] and n[1] 
 * values passed to mfftCreatePlan().
 *
 * Format is per the MatrixFFTPlan which produced the FFTs on 
 * which this operates; it's equal to MF_RowOrder or MF_CustomOrder,
 * depending on signal dimensions and configuration. 
 */
extern MFFTReturn fftConvDyadicMul(
	const FFTComplex		*src1,
	const FFTComplex		*src2,
	FFTComplex				*dst,
	unsigned				log2NumRows,
	unsigned				log2NumCols,
	MFFTFormat				format);
	
/* 
 * Perform conjugate dyadic multiply of the output of two 2-D real FFT ops.
 * For each complex value 0..N-1:
 *
 *    dst[n] = src1[n] * (src2[n]^*)
 *
 * ...where src2[n]^* is the complex conjugate of src2[n].
 *
 * Result goes to dst; inputs are not modified. 
 * The log2NumRows and log2NumCols values match the n[0] and n[1] 
 * values passed to mfftCreatePlan().
 *
 * Format is per the MatrixFFTPlan which produced the FFTs on 
 * which this operates; it's equal to MF_RowOrder or MF_CustomOrder,
 * depending on signal dimensions and configuration. 
 */
extern MFFTReturn fftConvDyadicMulConj(
	const FFTComplex		*src1,
	const FFTComplex		*src2,
	FFTComplex				*dst,
	unsigned				log2NumRows,
	unsigned				log2NumCols,
	MFFTFormat				format);

/* 
 * High-level interface:
 * -- 2-D real FFT buf1 in place
 * -- 2-D real FFT buf2 in place
 * -- dyadic multiply the result, normal (conj=false) or conjugate (conj=true)
 * -- inverse 2-D real FFT the result --> dst
 * -- normalize dst
 *
 * The MatrixFFTPlan is obtained from mfftCreatePlan() for a 2-D real FFT
 * at least as large as the signals on which this function operates. 
 */
extern MFFTReturn fftConvolve(
	MatrixFFTPlan			mfftPlan,
	FFTComplex				*buf1,
	FFTComplex				*buf2,
	unsigned				log2NumRows,
	unsigned				log2NumCols,
	bool					conj,
	FFTComplex				*dst);
	
#ifdef __cplusplus
}
#endif

#endif	/* _FFT_CONVOLVE_H_ */
