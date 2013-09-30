/*	File: vdspUtils.h 
	
	Description:
		Common public utility functions for vDSP ops.
	
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
 * vdspUtils.h - common public utility functions for vDSP ops.
 *
 * Created 9/23/2008. 
 * Copyright 2008 by Apple, Inc. 
 */
 
#ifndef	_VDSP_UTILS_H_
#define _VDSP_UTILS_H_

#include <libMatrixFFT/fftPrecision.h>
#include <libMatrixFFT/complexBuf.h>
#include <Accelerate/Accelerate.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Redirected vDSP types and function calls */
#if		FFT_DOUBLE_PREC

typedef	DSPDoubleSplitComplex				vDSPComplex;	/* for vDSP ops */
typedef DSPDoubleComplex					vDSPComplexInt;
typedef FFTSetupD							FFT_Setup;

/* redirected vDSP function calls */
#define FFTCreateSetup(n)					vDSP_create_fftsetupD(n, FFT_RADIX2)
#define FFTFreeSetup(f)						vDSP_destroy_fftsetupD(f)
#define FFTReal2d(f, c, lgCol, lcRow, dir)	vDSP_fft2d_zripD(f, c, 1, 0, lgCol, lcRow, dir)
#define FFTReal1d(f, c, lgN, dir)			vDSP_fft_zripD(f, c, 1, lgN, dir)
#define FFTReal1dOP(f, ci, co, lgN, dir)	vDSP_fft_zropD(f, ci, 1, co, 1, lgN, dir)
#define FFTReal2dOP(f, ci, co, lgCol, lgRow, dir) \
                                            vDSP_fft2d_zropD(f, ci, 1, 0, co, 1, 0, lgCol, lgRow, dir)
#define FFTComplex1d(f, c, lgN, dir)		vDSP_fft_zipD(f, c, 1, lgN, dir)
#define FFTComplex2d(f, c, lgNc, lgNr, dir)	vDSP_fft2d_zipD(f, c, 1, 0, lgNc, lgNr, dir)
#define FFTComplex1dOP(f, ci, co, lgN, dir)			\
                                            vDSP_fft_zopD(f, ci, 1, co, 1, lgN, dir)
#define FFTComplex2dOP(f, ci, co, lgNc, lgNr, dir)	\
                                            vDSP_fft2d_zopD(f, ci, 1, 0, co, 1, 0, lgNc, lgNr, dir)
        
/* buf[x] *= scale, x = 0..size-1 */
#define FFTvScale(buf, scale, size)			vDSP_vsmulD(buf, 1, &scale, buf, 1, size)

/* Interleaved inp --> DSP split outp */
#define FFT_vDSP_ctoz(inp, outp, len)		vDSP_ctozD(inp, 2, outp, 1, len)

/* DSP split inp --> interleaved outp */
#define FFT_vDSP_ztoc(inp, outp, len)		vDSP_ztocD(inp, 1, outp, 2, len)

/* C[x] := A[x] * B[x], x = 0..N-1 */
#define FFT_vmul(A, B, C, N)                vDSP_vmulD(A, 1, B, 1, C, 1, N)

/* 
 * C[x] := A, x = 0..N-1
 * careful, the vDSP function takes a ptr to A as the first arg, but our macro doesn't 
 */
#define FFT_vfill(A, C, N)                  vDSP_vfillD(&A, C, 1, N)

#else	/* Single precision */

typedef	DSPSplitComplex						vDSPComplex;
typedef DSPComplex							vDSPComplexInt;
typedef FFTSetup							FFT_Setup;

#define FFTCreateSetup(n)					vDSP_create_fftsetup(n, FFT_RADIX2)
#define FFTFreeSetup(f)						vDSP_destroy_fftsetup(f)
#define FFTReal2d(f, c, lgCol, lgRow, dir)	vDSP_fft2d_zrip(f, c, 1, 0, lgCol, lgRow, dir)
#define FFTReal1d(f, c, lgN, dir)			vDSP_fft_zrip(f, c, 1, lgN, dir)
#define FFTReal1dOP(f, ci, co, lgN, dir)	vDSP_fft_zrop(f, ci, 1, co, 1, lgN, dir)
#define FFTReal2dOP(f, ci, co, lgCol, lgRow, dir) \
                                            vDSP_fft2d_zrop(f, ci, 1, 0, co, 1, 0, lgCol, lgRow, dir)
#define FFTComplex1d(f, c, lgN, dir)		vDSP_fft_zip(f, c, 1, lgN, dir)
#define FFTComplex2d(f, c, lgNc, lgNr, dir)	vDSP_fft2d_zip(f, c, 1, 0, lgNc, lgNr, dir)
#define FFTComplex1dOP(f, ci, co, lgN, dir)			\
                                            vDSP_fft_zop(f, ci, 1, co, 1, lgN, dir)
#define FFTComplex2dOP(f, ci, co, lgNc, lgNr, dir)	\
                                            vDSP_fft2d_zop(f, ci, 1, 0, co, 1, 0, lgNc, lgNr, dir)
#define FFTvScale(buf, scale, size)			vDSP_vsmul(buf, 1, &scale, buf, 1, size)
#define FFT_vDSP_ctoz(inp, outp, len)		vDSP_ctoz(inp, 2, outp, 1, len)
#define FFT_vDSP_ztoc(inp, outp, len)		vDSP_ztoc(inp, 1, outp, 2, len)
#define FFT_vmul(A, B, C, N)                vDSP_vmul(A, 1, B, 1, C, 1, N)
#define FFT_vfill(A, C, N)                  vDSP_vfill(&A, C, 1, N)

#endif	/* FFT_DOUBLE_PREC */

#pragma mark --- Buffer allocation ---

/*
 * Allocate an vDSPComplex with specified number of complex elements.
 * Returns nonzero on malloc failure. Caller must free the returned
 * vDSPComplex pointers.
 */
extern int fftAllocDSPComplex(
	vDSPComplex				*buf,
	size_t					numComplex);
	
/*
 * Allocate an vDSPComplex with specified number of complex elements,
 * aligned to specified alignSize, which must be a power of 2. 
 * Returns nonzero on malloc failure. 
 * Caller must free the contents of the returned freeBuf. 
 */
extern int fftAllocDSPComplexAlign(
	vDSPComplex				*buf,
	size_t					numComplex,
	unsigned				alignSize,
	vDSPComplex				*freeBuf);

extern void fftFreeDSPComplex(
	vDSPComplex				*buf);
	
/* copy from one vDSPComplex to another */
extern void fftCopyDSPComplex(
	const vDSPComplex		*src,
	vDSPComplex				*dst,
	size_t					numSamples);
	
/* Multiply each element in an vDSPComplex by specified scale factor */
extern void fftScaleDSPComplex(
	vDSPComplex *buf,
	FFTFloat scaleFact,
	size_t numSamples);		/* in each of realp, imag */

/* copy from float array (i.e. vImage format) to vDSPComplex */
void fftCopyFloatToDSPComplex(
	const float *src,
	vDSPComplex *dst,
	size_t numSamples);			/* total real-signal samples */
	
#pragma mark --- Buffer initialization ---

/* random complex data */
extern void genRandComplexDSP(
	vDSPComplex *zBuf,
	size_t numSamples);

/* incrementing complex data */
extern void genIncrComplexDSP(
	vDSPComplex *buf,
	size_t numSamples);			/* in each of {real,imag} */

/* Generate 1-dimension "chirp" signal */
extern void fftGenChirp1DDSP(
	vDSPComplex *buf,
	size_t numSamples);			/* in each of {real,imag} */

/* Generate 2-dimension "chirp" signal */
extern void fftGenChirp2DDSP(
	vDSPComplex *buf,
	size_t width,
	size_t height);

/* flush contents of an vDSPComplex from cache */
extern void fftFlushComplexDSP(
	vDSPComplex *buf,
	size_t numComplex);

#pragma mark --- Buffer comparison and analysis ---

/* 
 * Analyze amplitudes of two complex signals. Return max delta and RMSE between
 * the two.  
 */
extern void fftCompareAmplitudesDSP(
	const vDSPComplex *buf1,
	const vDSPComplex *buf2,
	size_t numSamples,			/* In each of real & imag */
	double *maxDelta,			/* RETURNED */
	double *rmse);				/* RETURNED */

/* 
 * Analyze amplitudes of two real signals. Return max delta and RMSE between
 * the two.  
 */
extern void fftCompareRealAmplitudesDSP(
	const vDSPComplex *buf1,
	const vDSPComplex *buf2,
	size_t numSamples,			/* total real samples */
	double *maxDelta,			/* RETURNED */
	double *rmse);				/* RETURNED */

/* 
 * Find min and max values of amplitudes of all the elements of an FFT output.
 */ 
extern void fftAnalyzeAmplitudesDSP(
	const vDSPComplex *zBuf,
	size_t numSamples,
	double *minAmpl,		/* RETURNED */
	double *maxAmpl);		/* RETURNED */

#pragma mark --- Buffer dump/display ---

/* 
 * Dump up to 8 rows and columns of a square vDSPComplex data to stdout.
 */
extern void fftDumpDSPMatrix(	
	const char				*title,
	const vDSPComplex		*buf,
	unsigned				n);				/* total size = 2 ^ 2n */

/*
 * Dump up to 8 rows and columns of vDSPComplex data to stdout.
 */
extern void fftDumpDSPMatrixRect(	
	const char				*title,
	const vDSPComplex		*buf,
	size_t					numRows,
	size_t					numCols);	

/*
 * Glue for using fftDumpDSPMatrixRect() to dump an array of reals.
 * The numCols argument for fftDumpDSPMatrixRect() routine is in complex items.
 */
#define fftDumpDSPMatrixReal(title, buf, r, c)		fftDumpDSPMatrixRect(title, buf, r, (c)/2)

/* 
 * Dump 1- or 2-dimension DSPSplit*Complex data to stdout.
 */
extern void fftDump1DDSPComplexFloat(	
	const char					*title,
	const DSPSplitComplex		*buf,
	size_t						numSamples);	/* max to print */
extern void fftDump1DDSPComplexDouble(	
	const char					*title,
	const DSPDoubleSplitComplex	*buf,
	size_t						numSamples);	/* max to print */
extern void fftDump1DDSPComplex(	
	const char					*title,
	const vDSPComplex			*buf,
	size_t						numSamples);	/* max to print */
extern void fftDump2DDSPComplexFloat(	
	const char					*title,
	const DSPSplitComplex		*buf,
	size_t						numRows,		/* max to print */
	size_t						numCols);
extern void fftDump2DDSPComplexDouble(	
	const char					*title,
	const DSPDoubleSplitComplex	*buf,
	size_t						numRows,		/* max to print */
	size_t						numCols);
extern void fftDump2DDSPComplex(	
	const char					*title,
	const vDSPComplex			*buf,
	size_t						numRows,		/* max to print */
	size_t						numCols);

#if     !FFT_SPLIT_COMPLEX

/* 
 * Routines to convert between interleaved complex format and vDSP's native
 * split format.
 * NOTE vDSPComplex pointers are expected to be vector-aligned. 
 */
extern void fftVDSPToInt(
	const vDSPComplex *inBuf,
	FFTComplex *outBuf,
	size_t num);
extern void fftIntToVDSP(
	const FFTComplex *inBuf,
	vDSPComplex *outBuf,
	size_t num);

#endif	/* FFT_SPLIT_COMPLEX */

#ifdef __cplusplus
}
#endif

#endif	/* _VDSP_UTILS_H_ */

