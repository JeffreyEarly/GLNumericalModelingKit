/*	File: complexBufUtils.h
	
	Description:
		Routines for manipulating FFTComplex buffers.
	
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
 * complexBufUtils.h - Routines for manipulating FFTComplex buffers.
 *			  
 * Created 12/16/2008. 
 * Copyright 2008 by Apple, Inc. 
 */
 
#ifndef	_COMPLEX_BUF_UTILS_H_
#define _COMPLEX_BUF_UTILS_H_

#include <libMatrixFFT/complexBuf.h>
#include <libMatrixFFT/vdspUtils.h>
#include <Accelerate/Accelerate.h>

#ifdef __cplusplus
extern "C" {
#endif

#pragma mark --- Buffer allocation ---

/*
 * Allocate an FFTComplex array with specified number of complex elements.
 * Returns NULL on malloc failure. Caller must free the returned pointer
 * via fftFreeComplexArray().
 */
extern FFTComplex *fftAllocComplexArray(
	size_t					numComplex);

extern void fftFreeComplexArray(
	FFTComplex				*buf);

/*
 * Allocate an array of well-aligned memory. 
 * Free the returned freePtr free().
 */
extern void *fftAllocAlign(
    size_t                  allocSize,      // to be allocated
    void                    **freePtr);     // to be freed
    
/*
 * Allocate well-aligned FFTComplex array. Returns the aligned array.
 * Free the returned freePtr via fftFreeComplexArrayAlign().
 * Specified alignSize must be a power of 2. 
 */
extern FFTComplex *fftAllocComplexArrayAlign(
	size_t					numComplex,
	size_t					alignSize,		// in bytes
	FFTComplex				**freePtr);		// to be freed

extern void fftFreeComplexArrayAlign(
	FFTComplex				*buf,
    FFTComplex              *freeBuf);


/* copy from one FFTComplex to another */
extern void fftCopyComplexArray(
	const FFTComplex		*src,
	FFTComplex				*dst,
	size_t					numSamples);
	
/* Init one buffer as specified offset (both real and imag) from another */
extern void fftComplexOffset(
	const FFTComplex		*src,
	size_t					offset,
	FFTComplex				*dst);

#if     FFT_SPLIT_COMPLEX
/* Add specified offset to {real,imag} */	
extern void fftComplexAddOffset(
	FFTComplex				*buf,
	size_t					offset);
#endif

/* Return true if two FFTComplexes refer to same memory, else return false */
extern bool fftComplexEquiv(
    const FFTComplex        *buf1,
    const FFTComplex        *buf2);
    
#pragma mark --- Signal initialization ---

/* generate random complex signal */
extern void genRandComplex(
	FFTComplex *buf,
	size_t numSamples);			/* in each of {real,imag} */

/* generate complex signal, incrementing data */
extern void genIncrComplex(
	FFTComplex *buf,
	size_t numSamples);			/* in each of {real,imag} */

/* generate random FFTFloat signal */
extern void genRandFloat(
	FFTFloat *buf,
	size_t numSamples);	

/* generate incrementing FFTFloat signal */
extern void genIncrFloat(
	FFTFloat *buf,
	size_t numSamples);	

/* generate constant signal */
extern void genConstComplex(
	FFTComplex *buf,
	size_t numSamples,
	FFTFloat f);
	
/* Multiply each element in an FFTComplex array by specified scale factor */
extern void fftScaleComplex(
	FFTComplex *buf,
	FFTFloat scaleFact,
	size_t numSamples);		/* in each of realp, imag */

/* Generate random tau array for NUFFT */
extern void fftGenTau(
    FFTFloat *tau,
    size_t N);

#pragma mark --- Copy between FFTComplex and vDSPComplex ---

extern void fftCopyToDSP(
	const FFTComplex *src,
	vDSPComplex *dst,
	size_t numComplex);
	
extern void fftCopyFromDSP(
	const vDSPComplex *src,
	FFTComplex *dst,
	size_t numComplex);

/* copy from float array (i.e. vImage format) to FFTComplex */
extern void fftCopyFloatToSplit(
	const float *src,
	FFTComplex *dst,
	size_t numSamples);			/* total real-signal samples */

#pragma mark --- Misc. Buffer manipulation ---

/* Multiply a complex array by a vector of FFTFloats */
extern void fftMulComplexVect(
    FFTComplex *cbuf,           /* cbuf[n] *= fbuf[n] */
    const FFTFloat *fbuf,
    size_t N);
    
/* flush contents of an FFTComplex from cache */
extern void fftFlushComplex(
	FFTComplex *buf,
	size_t numComplex);	
    
#pragma mark --- Buffer comparison and analysis ---

/* 
 * Analyze amplitudes of two complex signals. Return max delta and RMSE between
 * the two.  
 */
extern void fftCompareAmplitudes(
	const FFTComplex *buf1,
	const FFTComplex *buf2,
	size_t numSamples,			/* In each of real & imag */
	double *maxDelta,			/* RETURNED */
	double *rmse);				/* RETURNED */

/* 
 * Analyze amplitudes of two complex signals - one in MatrixFFT form, one in 
 * vDSP form. Return max delta and RMSE between the two.  
 */
extern void fftCompareAmplitudesWithDSP(
	const FFTComplex *buf1,
	const vDSPComplex *buf2,
	size_t numSamples,			/* In each of real & imag */
	double *maxDelta,			/* RETURNED */
	double *rmse);				/* RETURNED */

/* 
 * Analyze amplitudes of two real signals. Return max delta and RMSE between
 * the two.  
 */
void fftCompareRealAmplitudes(
	const FFTComplex *buf1,
	const FFTComplex *buf2,
	size_t numSamples,			/* total real samples */
	double *maxDelta,			/* RETURNED */
	double *rmse);				/* RETURNED */

/* 
 * Analyze amplitudes of two real signals - one in MatrixFFT form, one in 
 * vDSP form. Return max delta and RMSE between the two.  
 */
extern void fftCompareRealAmplitudesWithDSP(
	const FFTComplex *buf1,
	const vDSPComplex *buf2,
	size_t numSamples,			/* total real samples */
	double *maxDelta,			/* RETURNED */
	double *rmse);				/* RETURNED */

/* 
 * Find min and max values of amplitudes of all the elements of an FFT output.
 */ 
extern void fftAnalyzeAmplitudes(
	const FFTComplex *zBuf,
	size_t numSamples,
	double *minAmpl,		/* RETURNED */
	double *maxAmpl);		/* RETURNED */

#pragma mark --- Buffer dump/display ---

/* 
 * Dump up to 8 rows and columns of a square FFTComplex data to stdout.
 */
extern void fftDumpSquare(	
	const char				*title,
	const FFTComplex		*buf,
	unsigned				n);				/* total size = 2 ^ 2n */

/*
 * Dump up to 8 rows and columns of FFTComplex data to stdout.
 */
extern void fftDumpMatrixRect(	
	const char				*title,
	const FFTComplex		*buf,
	size_t					numRows,
	size_t					numCols);	

/*
 * Glue for using fftDumpMatrixRect() to dump an array of reals.
 * The numCols argument for fftDumpMatrixRect() routine is in complex items.
 */
#define fftDumpMatrixReal(title, buf, r, c)		fftDumpMatrixRect(title, buf, r, (c)/2)

/* dump general FFTFloat buffer to stdout */
extern void fftDumpBuf(	
	const char			*title,
	const FFTFloat		*buf,
	size_t				numRows,
	size_t				numCols);	

void fftDump1DFloat(	
    const char			*title,
    const FFTFloat		*buf,
    size_t				numSamples);
    
/* 
 * Dump 1- or 2-dimension FFTComplex data to stdout.
 */
extern void fftDump1DComplex(	
	const char					*title,
	const FFTComplex			*buf,
	size_t						numSamples);	/* max to print */
extern void fftDump2DComplex(	
	const char					*title,
	const FFTComplex			*buf,
	size_t						numRows,		/* max to print */
	size_t						numCols);

/* Dump an array of size_t's to stdout */
extern void fftDump1DSize(
    const char                  *title,
    const size_t                *s,
    size_t                      numSamples);
    
#pragma mark --- Format-independent Complex manipulation ---

/*
 * These macros can be used to write (some) code which manipulates complex
 * data in a manner which is independent of split/interleaved format. 
 * They use a local variable whose name is passed as the first argument to 
 * FC_DECL() or FC_DECL_SET(); those macros declare the local, and 
 * FC_DECL_SET initializes the local to a specified FFTComplex*.
 * Note the two versions of these declare quite different things; if 
 * you use the declared local outside the scope of these macros you really
 * need to understand what's going on here. 
 *
 * FC_SET assigns the specified FFTComplex* value to the specified local.
 *
 * The FC_OFF() macro adds an offset, in complex elements, to the declared
 * variable. This adds the specified offset to each of {real,imag} in the 
 * split complex case, and adds the specified offset to the local pointer 
 * in the interleaved complex case.  
 *
 * FC_REAL() and FC_IMAG() access the current FFTFloat components. Each
 * can be used as an lval or as a term in an rval. But do NOT attempt
 * a ++ or -- operation when using them.
 *
 * NOTE: C++ code can use the much better mechanism in src/PolyComplex.h 
 * to do all of this and much more. 
 */
#if		FFT_SPLIT_COMPLEX

#define FC_DECL(varName)				FFTComplex varName
#define FC_DECL_SET(varName, initVal)	FFTComplex varName = *initVal
#define FC_PTR(varName)                 (&varName)
#define FC_SET(varName, initVal)		{ varName = *initVal; }
#define FC_OFF(varName, off)			{ varName.real += off; varName.imag += off; }
#define FC_REAL(varName)				(*(varName.real))
#define FC_IMAG(varName)				(*(varName.imag))

#else	/* !FFT_SPLIT_COMPLEX */

#define FC_DECL(varName, initVal)		FFTComplex *varName
#define FC_DECL_SET(varName, initVal)	FFTComplex *varName = initVal
#define FC_PTR(varName)                 (varName)
#define FC_SET(varName, initVal)		{ varName = initVal; }
#define FC_OFF(varName, off)			{ varName += off; }
#define FC_REAL(varName)				(varName->real)
#define FC_IMAG(varName)				(varName->imag)

#endif	/* FFT_SPLIT_COMPLEX */

#ifdef __cplusplus
}
#endif

#endif	/* _COMPLEX_BUF_UTILS_H_ */
