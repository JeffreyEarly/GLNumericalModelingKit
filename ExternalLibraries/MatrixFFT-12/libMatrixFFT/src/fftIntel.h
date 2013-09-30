/*	File:  fftIntel.h 

	Description:
		Private Intel-only utility functions for MatrixFFT library.
	
	Copyright:
		Copyright (C) 2009 Apple Inc.  All rights reserved.
	
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
 * Created 01/06/2009. 
 * Copyright 2009 by Apple, Inc. 
 */
 
#ifndef	_FFT_INTEL_H_
#define _FFT_INTEL_H_

#include "fftPriv.h"

#ifndef	FFT_INTEL
#error You must #define FFT_INTEL.
#endif

#if		FFT_INTEL

#ifdef __cplusplus
extern "C" {
#endif

/* 
 * For precision-independent vectorized Intel code, we extend
 * the redirected types and functions in fftPrecision.h like so.
 */
#if		FFT_DOUBLE_PREC

typedef __m128d					FFTVector;

#define FFTVectLoad(p)			_mm_load_pd(p)
#define FFTVectLoadRev(p)		_mm_loadr_pd(p)			/* load in reverse order */
#define FFTVectStore(p, v)		_mm_store_pd(p, v)
#define FFTVectStoreRev(p, v)	_mm_storer_pd(p, v)		/* store in reverse order */
#define FFTVectAdd(a, b)		_mm_add_pd(a, b)
#define FFTVectSub(a, b)		_mm_sub_pd(a, b)		/* result = a - b */
#define FFTVectMul(a, b)		_mm_mul_pd(a, b)
#define FFTVectSet1(a)			_mm_set1_pd(a)			/* set all to a */

#ifdef __llvm__
/* llvm doesn't have __builtin_ia32_shufpd */
#define __builtin_ia32_shufpd(v1, v2, b)  _mm_shuffle_pd(v1, v2, b)
#endif  /* __GNUC__ */
    
#else	/* !FFT_DOUBLE_PREC */

typedef __m128					FFTVector;

#define FFTVectLoad(p)			_mm_load_ps(p)
#define FFTVectLoadRev(p)		_mm_loadr_ps(p)			/* load in reverse order */
#define FFTVectStore(p, v)		_mm_store_ps(p, v)
#define FFTVectStoreRev(p, v)	_mm_storer_ps(p, v)		/* store in reverse order */
#define FFTVectAdd(a, b)		_mm_add_ps(a, b)
#define FFTVectSub(a, b)		_mm_sub_ps(a, b)
#define FFTVectMul(a, b)		_mm_mul_ps(a, b)
#define FFTVectSet1(a)			_mm_set1_ps(a)

#ifdef __llvm__
/* llvm doesn't have __builtin_ia32_shufps */
#define __builtin_ia32_shufps(v1, v2, b)  _mm_shuffle_ps(v1, v2, b)
#endif  /* __GNUC__ */

#endif	/* FFT_DOUBLE_PREC */

typedef union {
	FFTVector	v;
	FFTFloat	f[FFT_FLOATS_PER_VECTOR];
} FFTVectUnion;

#pragma mark --- Interleaved complex support --- */
 
#if		!FFT_SPLIT_COMPLEX

/*
 * Routines to load and store real/imaginary pairs to interleaved-format
 * FFTComplex* . Inlined for performance reasons; these are used
 * extensively.
 */
 
#if		FFT_DOUBLE_PREC

/* load two vectors - one real, one imaginary - from a FFTComplex * */
static inline void fftLoadComplex(
	const FFTComplex	*comp,
	FFTVector			&imag,
	FFTVector			&real)
{
	/* 
	 * Layout in memory is:
	 *		real0 imag0  real1 imag1  
	 * And we want the returned real and imag to be 
	 *		{ real0 real1 }
	 *		{ imag0 imag1 }
	 */
	const FFTVector *vp = (FFTVector *)comp;
	FFTVector lo = vp[0];
	FFTVector hi = vp[1];
	real = __builtin_ia32_shufpd(lo, hi, _MM_SHUFFLE2(0, 0));
	imag = __builtin_ia32_shufpd(lo, hi, _MM_SHUFFLE2(1, 1));
}

/* Store two vectors to a FFTComplex * */
static inline void fftStoreComplex(
	const FFTComplex	*comp,
	FFTVector			imag,
	FFTVector			real)
{
	/*
	 * real  = { real0 real1 }
	 * imag  = { imag0 imag1 }
	 * 
	 * And we want
	 * lo    = { real0 imag0  }
	 * hi    = { real1 imag1 }
	 */
	FFTVector lo = __builtin_ia32_shufpd(real, imag, _MM_SHUFFLE2(0, 0));
	FFTVector hi = __builtin_ia32_shufpd(real, imag, _MM_SHUFFLE2(1, 1));
	FFTVector *vp = (FFTVector *)comp;
	vp[0] = lo;
	vp[1] = hi;
}

/* load two vectors - one real, one imaginary - in reverse order */
static inline void fftLoadComplexRev(
	const FFTComplex	*comp,
	FFTVector			&imag,
	FFTVector			&real)
{
	/* 
	 * Layout in memory is:
	 *		real0 imag0  real1 imag1  
	 * And we want the returned real and imag to be 
	 *		{ real1 real0 }
	 *		{ imag1 imag0 }
	 */
	const FFTVector *vp = (FFTVector *)comp;
	FFTVector lo = vp[0];
	FFTVector hi = vp[1];
	real = __builtin_ia32_shufpd(hi, lo, _MM_SHUFFLE2(0, 0));
	imag = __builtin_ia32_shufpd(hi, lo, _MM_SHUFFLE2(1, 1));
}

/* Store two vectors to a FFTComplex * in reverse order */
static inline void fftStoreComplexRev(
	const FFTComplex	*comp,
	FFTVector			imag,
	FFTVector			real)
{
	/*
	 * real  = { real0 real1 }
	 * imag  = { imag0 imag1 }
	 * 
	 * And we want
	 * lo    = { real1 imag1  }
	 * hi    = { real0 imag0 }
	 */
	FFTVector lo = __builtin_ia32_shufpd(real, imag, _MM_SHUFFLE2(1, 1));
	FFTVector hi = __builtin_ia32_shufpd(real, imag, _MM_SHUFFLE2(0, 0));
	FFTVector *vp = (FFTVector *)comp;
	vp[0] = lo;
	vp[1] = hi;
}


#else	/* !FFT_DOUBLE_PREC */

/* load two vectors - one real, one imaginary - from a FFTComplex * */
static inline void fftLoadComplex(
	const FFTComplex	*comp,
	FFTVector			&imag,
	FFTVector			&real)
{
	/*
	 * Layout in memory is:
	 *		real0 imag0  real1 imag1  real2 imag2  real3 imag3
	 * Initially real and imag are
	 *		{ real0 imag0 real1 imag1 }
	 * And hi is
	 *		{ real2 imag2 real3 imag3 }
	 * And we want the returned real and imag to be
	 *		{ real0 real1 real2 real3 }
	 *		{ imag0 imag1 imag2 imag3 }
	 */
	FFTVector *vp = (FFTVector *)comp;
	imag = real = vp[0];
	FFTVector hi = vp[1];
	real = __builtin_ia32_shufps(real, hi, _MM_SHUFFLE(2, 0, 2, 0));
	imag = __builtin_ia32_shufps(imag, hi, _MM_SHUFFLE(3, 1, 3, 1));
}

/* Store two vectors to a FFTComplex * */
static inline void fftStoreComplex(
	const FFTComplex	*comp,
	FFTVector			imag,
	FFTVector			real)
{
	FFTVector lo = real;
	FFTVector hi = imag;
	
	/*
	 * lo    = { real0 real1 real2 real3 }
	 * hi    = { imag0 imag1 imag2 imag3 }
	 * And we want
	 * lo    = { real0 imag0 real1 imag1 }
	 * hi    = { real2 imag2 real3 imag3 }
	 * But we have to go thru this first
	 * tmpLo = { real0 real1 imag0 imag1 }
	 * tmpHi = { real2 real3 imag2 imag3 }
	 */
	FFTVector tmpLo = __builtin_ia32_shufps(lo, hi, _MM_SHUFFLE(1, 0, 1, 0));
	FFTVector tmpHi = __builtin_ia32_shufps(lo, hi, _MM_SHUFFLE(3, 2, 3, 2));
	lo = __builtin_ia32_shufps(tmpLo, tmpLo, _MM_SHUFFLE(3, 1, 2, 0));
	hi = __builtin_ia32_shufps(tmpHi, tmpHi, _MM_SHUFFLE(3, 1, 2, 0));
	FFTVector *vp = (FFTVector *)comp;
	vp[0] = lo;
	vp[1] = hi;
}

/* 
 * load two vectors - one real, one imaginary - from a FFTComplex * in 
 * reverse order
 */
static inline void fftLoadComplexRev(
	const FFTComplex	*comp,
	FFTVector			&imag,
	FFTVector			&real)
{
	/*
	 * Layout in memory is:
	 *		real0 imag0  real1 imag1  real2 imag2  real3 imag3
	 * Initially lo and tmp are
	 *		{ real0 imag0 real1 imag1 }
	 * And imag and real are
	 *		{ real2 imag2 real3 imag3 }
	 * And we want the returned real and imag to be
	 *		{ real3 real2 real1 real0 }
	 *		{ imag3 imag2 imag1 imag0 }
	 */
	FFTVector *vp = (FFTVector *)comp;
	FFTVector lo = vp[0];
	imag = real = vp[1];
	real = __builtin_ia32_shufps(real, lo, _MM_SHUFFLE(0, 2, 0, 2));
	imag = __builtin_ia32_shufps(imag, lo, _MM_SHUFFLE(1, 3, 1, 3));
}

/* Store two vectors to a FFTComplex * in reverse order */
static inline void fftStoreComplexRev(
	const FFTComplex	*comp,
	FFTVector			imag,
	FFTVector			real)
{
	FFTVector lo = real;
	FFTVector hi = imag;
	
	/*
	 * lo    = { real0 real1 real2 real3 }
	 * hi    = { imag0 imag1 imag2 imag3 }
	 * And we want
	 * lo    = { real3 imag3 real2 imag2 }
	 * hi    = { real1 imag1 real0 imag0 }
	 * But we have to go thru this first
	 * tmpLo = { real0 real1 imag0 imag1 }
	 * tmpHi = { real2 real3 imag2 imag3 }
	 */
	FFTVector tmpLo = __builtin_ia32_shufps(lo, hi, _MM_SHUFFLE(1, 0, 1, 0));
	FFTVector tmpHi = __builtin_ia32_shufps(lo, hi, _MM_SHUFFLE(3, 2, 3, 2));
	hi = __builtin_ia32_shufps(tmpLo, tmpLo, _MM_SHUFFLE(2, 0, 3, 1));
	lo = __builtin_ia32_shufps(tmpHi, tmpHi, _MM_SHUFFLE(2, 0, 3, 1));
	FFTVector *vp = (FFTVector *)comp;
	vp[0] = lo;
	vp[1] = hi;
}

#endif	/* FFT_DOUBLE_PREC */
#endif	/* !FFT_SPLIT_COMPLEX */

#pragma mark --- Vector sine/cosine ---

/*
 * Unfortunately there is only a single-precision version of a vDSP function
 * to take the sin and cos of a vector. We have to jerry-rig double precision
 * versions here. 
 */
#if     FFT_DOUBLE_PREC

static inline FFTVector fftVectCosine(FFTVector x)
{
    FFTVectUnion vi;
    FFTVectUnion vo;
    vi.v = x;
    double *d_i = vi.f;
    double *d_o = vo.f;
    
    /* We know this is double precision only, so unroll the expected loop */
    d_o[0] = cos(d_i[0]);
    d_o[1] = cos(d_i[1]);

    return vo.v;
}

static inline FFTVector fftVectSine(FFTVector x)
{
    FFTVectUnion vi;
    FFTVectUnion vo;
    vi.v = x;
    double *d_i = vi.f;
    double *d_o = vo.f;
    
    d_o[0] = sin(d_i[0]);
    d_o[1] = sin(d_i[1]);
    return vo.v;
}


#else   /* !FFT_DOUBLE_PREC */

static inline FFTVector fftVectCosine(FFTVector x)
{
    return vcosf(x);
}

static inline FFTVector fftVectSine(FFTVector x)
{
    return vsinf(x);
}

#endif  /* FFT_DOUBLE_PREC */

#ifdef __cplusplus
}
#endif

#endif	/* FFT_INTEL */

#endif	/* _FFT_INTEL_H_ */
