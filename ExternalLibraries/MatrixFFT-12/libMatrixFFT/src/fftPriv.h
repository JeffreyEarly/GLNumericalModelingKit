/*	File:  fftPriv.h 

	Description:
		Common private utility functions for MatrixFFT library.
	
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
 * fftPriv.h - common private utility functions for MatrixFFT library.
 *			  
 * Note: the precision of the floats used in this module are determined at 
 * compile time via the MatrixFFTConfig.h module. 
 *
 * Created 9/25/2008. 
 * Copyright 2008 by Apple, Inc. 
 */
 
#ifndef	_FFT_PRIV_H_
#define _FFT_PRIV_H_

#include <stdint.h>
#include <libMatrixFFT/MatrixFFT.h>
#include "MatrixFFTPlan.h"
#include <Accelerate/Accelerate.h>

#ifdef __cplusplus
extern "C" {
#endif


#define FFT_SWAP_VALS(elt, v1, v2)	{	\
	elt tmp = v1;						\
	v1 = v2;							\
	v2 = tmp;							\
}

#pragma mark --- Platform-dependent submatrix support ---

/*
 * Platform-dependent submatrix support. 
 * See source files e.g. fftTransposeOP.cpp for info on submatrices. 
 */
#if defined( __i386__ ) || defined( __x86_64__ )

/* collate these two compiler-generated symbols into one */
#define FFT_INTEL	1
#else
#define FFT_INTEL	0
#endif

#if FFT_INTEL

/* 
 * Machine-dependent vector and cache line sizes.
 * WANRING: if CACHE_LINE_SIZE changes, update FFT_MEM_ALIGNMENT in the (public)
 * MatrixFFT.h. Technically, FFT_MEM_ALIGNMENT is the larger of (CACHE_LINE_SIZE,
 * VECTOR_SIZE), but we don't want these machine-dependent #defines visible in the 
 * public API.
 */
#define CACHE_LINE_SIZE			64 			/* in bytes */
#define VECTOR_SIZE				16			/* bytes */

/* enable prefetch */
/* as of 11/11/08 I see slightly better performance with this on, at least on 8-Core MacPro */
#define FFT_ENABLE_PREFETCH		1
#if		FFT_ENABLE_PREFETCH

#define FFT_PREFETCH_LEVEL		_MM_HINT_T0
#define fft_prefetch(a)			_mm_prefetch((const char *)a, FFT_PREFETCH_LEVEL)

#else	/* FFT_ENABLE_PREFETCH */

#define fft_prefetch(a)

#endif	/* FFT_ENABLE_PREFETCH */

#else /* !x86 */

/* 
 * These are here so subsequent code compiles on PPC.
 * I don't think these values are right but we don't have any code that
 * uses them for PPC.
 */

#define CACHE_LINE_SIZE			64 			/* in bytes */
#define VECTOR_SIZE				16			/* bytes */

#define fft_prefetch(a)

#endif	/* x86 */

/* 
 * Platform dependent constants, expressed here in a platform-independent manner
 * as long as CACHE_LINE_SIZE and VECTOR_SIZE have been defined in machine-dependent 
 * section above.
 */

/*
 * Machine-dependent vectors in a cache line.
 */
#define FFT_VECTORS_PER_SUBMATRIX	(CACHE_LINE_SIZE / VECTOR_SIZE)

/* 
 * Number of configuration-dependent FFTFloats in a machine-dependent vector.
 */
#define FFT_FLOATS_PER_VECTOR		(VECTOR_SIZE / sizeof(FFTFloat))

/* 
 * Number of configuration-dependent FFTFloats in one side of a submatrix.
 * A submatrix is a cache line squared.
 */
#define FFT_FLOATS_PER_SUBMATRIX	(CACHE_LINE_SIZE / sizeof(FFTFloat))

#if	!FFT_SPLIT_COMPLEX

/* 
 * FFT_COMPLEX_PER_{VECTOR,SUBMATRIX} is only meaningful for !FFT_SPLIT_COMPLEX.
 */
 
/* 
 * Number of configuration-dependent FFTComplex's in a machine-dependent vector.
 */
#define FFT_COMPLEX_PER_VECTOR		(VECTOR_SIZE / sizeof(FFTComplex))

/* 
 * Number of configuration-dependent FFTComplex's in a submatrix.
 */
#define FFT_COMPLEX_PER_SUBMATRIX	(CACHE_LINE_SIZE / sizeof(FFTComplex))

#endif	/* !FFT_SPLIT_COMPLEX */



/* 
 * Useful macros and types 
 */
#define CFRELEASE(cf)	if(cf != NULL) { CFRelease(cf); cf = NULL; }

/* 
 * Check alignment - alignSize must be a power of 2.
 * Works with integers and pointers.
 * Returns true if well-aligned.
 */
#define FFT_IS_ALIGNED(p, alignSize)		((((intptr_t)(p)) & (alignSize-1)) == 0)

/*
 * Align a value (pointer, integer) up to the next alignSize, which must
 * be a power of 2. You'll need to cast this as appropriate. 
 */
#define FFT_ALIGN(p, alignSize)		((((intptr_t)(p)) + alignSize - 1) & ~((intptr_t)alignSize - 1))

/* Enforce alignment of FFTComplex, any form - implemented as ASSERT for DEBUG builds */
#ifdef	DEBUG

extern void fftAssertBufAlign(
	FFTComplex	*buf,
	size_t		alignSize);
	
#else	/* DEBUG */

#define fftAssertBufAlign(buf, alignSize)

#endif	/* DEBUG */

#pragma mark --- Misc. private routines ---

/* 
 * Alloc or realloc mfftPlan->auxBuf. 
 */
extern MFFTReturn fftAllocAuxBuf(
	MatrixFFTPlan		mfftPlan,
	size_t				newSize);		/* in in FFTComplex elements */
	
/* 1-D Complex twist function, in fft1DTwist{Split,Int}.cpp */
extern MFFTReturn fft1DTwist(
	MatrixFFTPlan		mfftPlan,
	FFTComplex			*buf,
	bool				forward,
	size_t				numRows,
	size_t				numCols,
	size_t				startRow,
	size_t				rowsToProcess);

/* 1-D real twist function, in fft1DRealTwist.cpp */
extern MFFTReturn fft1DRealTwist(
	MatrixFFTPlan		mfftPlan,
	const FFTComplex	*src,
	FFTComplex			*dst,
	bool				forward,
	size_t				numRows,
	size_t				numCols,		/* in complex elements */
	size_t				startRow,
	size_t				rowsToProcess);

/*
 * Routines to copy to/from column 0 as complexes.
 * Implemented in fft2DReal.cpp.
 */

/* 
 * Copy the reals and then imaginaries from column 0 of inBuf as complexes to 
 * outBuf.
 */
extern void fftCopyFromColumn0(
	const FFTComplex *inBuf,
	FFTComplex *outBuf,			/* two rows of numComplex elements */
	size_t rowSize,				/* size of *inBuf row in FFTComplexes */
	size_t numComplex);			/* # of complex elements to move to each output row */

/*
 * Copy two sets of complexes from inBuf to the real and then the imaginary
 * components of column 0 of outBuf.
 */
extern void fftCopyToColumn0(
	const FFTComplex *inBuf,	/* two rows of numComplex elements */
	FFTComplex *outBuf,
	size_t rowSize,				/* size of *outBuf row in FFTComplexes */
	size_t numComplex);			/* # of complex elements to move from each input row */

/* Private dyadic mul, implementations in fftConvolve{Int,Split}.cpp */
extern MFFTReturn fftConvDyadicMulCom(
	const FFTComplex	*src1,
	const FFTComplex	*src2,
	FFTComplex			*dst,
	unsigned			log2NumRows,
	unsigned			log2NumCols,
	bool				conj);

/* for 2-D custom format */
extern MFFTReturn fftConvDyadicMulCustomCom(
	const FFTComplex	*src1,
	const FFTComplex	*src2,
	FFTComplex			*dst,
	unsigned			log2NumRows,
	unsigned			log2NumCols,
	bool				conj);

/* Multithreaded FFTvScale() */
extern MFFTReturn fftScaleThr(
	MatrixFFTPlan		mfftPlan,
	FFTComplex			*buf,
	FFTFloat			scaleFact,
	/* total size of buf */
	size_t				totalElts);

/* Force fft1DTwistSmall() to debug the rest */
#define FFT_FWD_TWIST_COMPLEX_SMALL		0

/* Sine reclaculation period */
#define FFT_SIN_RECALC_COMPLEX			(256)

#ifdef __cplusplus
}
#endif

#endif	/* _FFT_PRIV_H_ */
