/*	File: complexBufUtilsSplit.cpp 
	
	Description:
		FFTComplex buffer routines, split complex specific.
	
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
 * complexBufUtilsSplit.cpp - FFTComplex buffer routines, split complex specific.
 *			  
 * Created 11/20/2008. 
 * Copyright 2008 by Apple, Inc. 
 */
 

#include <libMatrixFFT/MatrixFFTConfig.h>

#if		FFT_SPLIT_COMPLEX

#include <libMatrixFFT/complexBufUtils.h>
#include <stdlib.h>
#include <strings.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <libMatrixFFT/devRandom.h>
#include <libMatrixFFT/fftUtils.h>
#include "fftPriv.h"

#pragma mark --- Buffer allocation ---

/*
 * Allocate an FFTComplex array with specified number of complex elements.
 * Returns NULL on malloc failure. Caller must free the returned
 * pointer via fftFreeComplexArray().
 */
FFTComplex *fftAllocComplexArray(
	size_t numComplex)
{
	FFTComplex *fftComplex = (FFTComplex *)malloc(sizeof(FFTComplex));
	if(fftComplex == NULL) {
		printf("***fftAllocComplexArray: malloc error\n");
		return NULL;
	}
	fftComplex->real = fftComplex->imag = NULL;
	if(numComplex == 0) {
		/* that's OK here, we're done */
		return fftComplex;
	}
	
	size_t bufSize = numComplex * sizeof(FFTFloat);
	fftComplex->real = (FFTFloat *)malloc(bufSize);
	fftComplex->imag = (FFTFloat *)malloc(bufSize);
	if((fftComplex->real == NULL) || (fftComplex->imag == NULL)) {
		printf("***fftAllocComplexArray: malloc error\n");
		COND_FREE(fftComplex->real);
		COND_FREE(fftComplex->imag);
		free(fftComplex);
		return NULL;
	}
	
	return fftComplex;
}

void fftFreeComplexArray(
	FFTComplex *buf)
{
	if(buf == NULL) {
		return;
	}
	COND_FREE(buf->real);
	COND_FREE(buf->imag);
	free(buf);
}


/*
 * Allocate well-aligned FFTComplex array. Returns the aligned array.
 * Free the returned freePtr via fftFreeComplexArrayAlign().
 * Specified alignSize must be a power of 2. 
 */
FFTComplex *fftAllocComplexArrayAlign(
	size_t					numComplex,
	size_t					alignSize,		// in bytes
	FFTComplex				**freePtr)		// to be freed
{
	*freePtr = NULL;
	
	/* alloc two empty FFTComplexes */
	FFTComplex *alignComplex = fftAllocComplexArray(0);
	FFTComplex *freeComplex  = fftAllocComplexArray(0);
	if((alignComplex == NULL) || (freeComplex == NULL)) {
		COND_FREE(freeComplex);
		COND_FREE(alignComplex);
		printf("***fftAllocComplexArrayAlign: malloc error\n");
		return NULL;
	}
	
	/* alloc unaligned pointers */
	size_t bufSize = (numComplex * sizeof(FFTFloat)) + alignSize;
	freeComplex->real = (FFTFloat *)malloc(bufSize);
	freeComplex->imag = (FFTFloat *)malloc(bufSize);
	if((freeComplex->real == NULL) || (freeComplex->imag == NULL)) {
		printf("***fftAllocComplexArrayAlign: malloc error\n");
		fftFreeComplexArray(freeComplex);
		fftFreeComplexArray(alignComplex);
		return NULL;
	}
	
	/* obtain aligned pointers */
	alignComplex->real = (FFTFloat *)FFT_ALIGN(freeComplex->real, alignSize);
	alignComplex->imag = (FFTFloat *)FFT_ALIGN(freeComplex->imag, alignSize);
	
	*freePtr = freeComplex;
	return alignComplex;
}

/*
 * This configuration has to free two FFTComplex's, plus the 
 * real/imag pointers in *alignBuf. 
 */
void fftFreeComplexArrayAlign(
	FFTComplex				*buf,
    FFTComplex              *freeBuf)
{
    fftFreeComplexArray(freeBuf);
    COND_FREE(buf);
}


/* copy from one FFTComplex to another */
void fftCopyComplexArray(
	const FFTComplex	*src,
	FFTComplex			*dst,
	size_t				numSamples)
{
	size_t bufSize = numSamples * sizeof(FFTFloat);
	memmove(dst->real, src->real, bufSize);
	memmove(dst->imag, src->imag, bufSize);
}

/* Init one buffer as specified offset (both real and imag) from another */
void fftComplexOffset(
	const FFTComplex		*src,
	size_t					offset,
	FFTComplex				*dst)
{
	*dst = *src;
	dst->real += offset;
	dst->imag += offset;
}
	
void fftComplexAddOffset(
	FFTComplex				*buf,
	size_t					offset)
{
	buf->real += offset;
	buf->imag += offset;
}

/* Multiply each element in an FFTComplex by specified scale factor */
void fftScaleComplex(
	FFTComplex *buf,
	FFTFloat scaleFact,
	size_t numSamples)		/* in each of realp, imag */
{
	FFTvScale(buf->real, scaleFact, numSamples);
	FFTvScale(buf->imag, scaleFact, numSamples);
}

/* Return true if two FFTComplexes refer to same memory, else return false */
bool fftComplexEquiv(
    const FFTComplex        *buf1,
    const FFTComplex        *buf2)
{
    if((buf1 == buf2) ||
       (buf1->real == buf2->real) ||
       (buf1->imag == buf2->imag)) {
        return true;
    }
    else {
        return false;
    }
}

#pragma mark --- Copy between FFTComplex and vDSPComplex ---

void fftCopyToDSP(
	const FFTComplex *src,
	vDSPComplex *dst,
	size_t numComplex)
{
	FFTComplex fftDst;
	vDSPComplexToFFT(dst, &fftDst);
	fftCopyComplexArray(src, &fftDst, numComplex);
}
	
void fftCopyFromDSP(
	const vDSPComplex *src,
	FFTComplex *dst,
	size_t numComplex)
{
	FFTComplex fftSrc;
	vDSPComplexToFFT(src, &fftSrc);
	fftCopyComplexArray(&fftSrc, dst, numComplex);
}

#pragma mark --- Misc. Buffer manipulation ---

/* generate constant signal */
void genConstComplex(
	FFTComplex *buf,
	size_t numSamples,
	FFTFloat f)
{ 
    FFT_vfill(f, buf->real, numSamples);
    FFT_vfill(f, buf->imag, numSamples);
}

/* Multiply a complex array by a vector of FFTFloats */
void fftMulComplexVect(
    FFTComplex *cbuf,           /* cbuf[n] *= fbuf[n] */
    const FFTFloat *fbuf,
    size_t N)
{
    if(FFT_IS_ALIGNED(cbuf->real, VECTOR_SIZE) &&
       FFT_IS_ALIGNED(cbuf->imag, VECTOR_SIZE) &&
       FFT_IS_ALIGNED(fbuf, VECTOR_SIZE)) {
        FFT_vmul(cbuf->real, fbuf, cbuf->real, N);
        FFT_vmul(cbuf->imag, fbuf, cbuf->imag, N);
    }
    else {
        /* brute force */
        FFTFloat *realp = cbuf->real;
        FFTFloat *imagp = cbuf->imag;
        for(size_t dex=0; dex<N; dex++) {
            FFTFloat f = *fbuf++;
            *realp *= f;
            *imagp *= f;
            realp++;
            imagp++;
        }
    }
}

/* flush contents of an FFTComplex from cache */
void fftFlushComplex(
	FFTComplex *buf,
	size_t numComplex)
{
    size_t numBytes = numComplex * sizeof(FFTFloat);
    fftFlushCache(buf->real, numBytes);
    fftFlushCache(buf->imag, numBytes);
    
}

#pragma mark --- Buffer comparison and analysis ---

/* 
 * Analyze amplitudes of two complex signals - one in MatrixFFT form, one in 
 * vDSP form. Return max delta and RMSE between the two.  
 */
void fftCompareAmplitudesWithDSP(
	const FFTComplex *buf1,
	const vDSPComplex *buf2,
	size_t numSamples,			/* In each of real & imag */
	double *maxDelta,			/* RETURNED */
	double *rmse)				/* RETURNED */
{
	FFTComplex fftBuf2;
	vDSPComplexToFFT(buf2, &fftBuf2);
	fftCompareAmplitudes(buf1, &fftBuf2, numSamples, maxDelta, rmse);
}

/* 
 * Analyze amplitudes of two real signals - one in MatrixFFT form, one in 
 * vDSP form. Return max delta and RMSE between the two.  
 */
void fftCompareRealAmplitudesWithDSP(
	const FFTComplex *buf1,
	const vDSPComplex *buf2,
	size_t numSamples,			/* total real samples */
	double *maxDelta,			/* RETURNED */
	double *rmse)				/* RETURNED */
{
	FFTComplex fftBuf2;
	vDSPComplexToFFT(buf2, &fftBuf2);
	fftCompareRealAmplitudes(buf1, &fftBuf2, numSamples, maxDelta, rmse);
}


#pragma mark --- Buffer dump/display ---

#define MAX_TO_PRINT		8
#define MAX_ROWS_TO_PRINT	16

extern void fftDumpSquare(	
	const char				*title,
	const FFTComplex		*buf,
	unsigned				n)				/* total size = 2 ^ 2n */
{
	if(title) {
		printf("%s:\n", title);
	}
	
	size_t x = (size_t)1 << n;
	size_t toPrint = x;
	const char *eolStr = "\n";
	if(toPrint > MAX_TO_PRINT) {
		toPrint = MAX_TO_PRINT;
		eolStr = " ...\n";
	}
	for(size_t row=0; row<toPrint; row++) {
		for(size_t col=0; col<toPrint; col++) {
			size_t off = (row * x) + col;
			printf("{%.2f,%.2f}, ", buf->real[off], buf->imag[off]);
		}
		printf("%s", eolStr);
	}
	if(x > MAX_TO_PRINT) {
		printf("%s", eolStr);
	}
}

void fftDumpMatrixRect(	
	const char				*title,
	const FFTComplex		*buf,
	size_t					numRows,
	size_t					numCols)	
{
	if(title) {
		printf("%s:\n", title);
	}
	
	const char *eolStr = "\n";
	const char *endStr = NULL;
	size_t rowsToPrint = numRows;
	if(rowsToPrint > MAX_ROWS_TO_PRINT) {
		rowsToPrint = MAX_ROWS_TO_PRINT;
		endStr = " ...\n";
	}
	size_t colsToPrint = numCols;
	if(colsToPrint > MAX_TO_PRINT) {
		colsToPrint = MAX_TO_PRINT;
		eolStr = " ...\n";
	}
	
	for(size_t row=0; row<rowsToPrint; row++) {
		for(size_t col=0; col<colsToPrint; col++) {
			size_t off = (row * numCols) + col;
			printf("{%.2f,%.2f}, ", buf->real[off], buf->imag[off]);
		}
		printf("%s", eolStr);
	}
	if(endStr) {
		printf("%s", endStr);
	}
}

/* 
 * Dump 1-dimension FFTComplex data to stdout.
 */
void fftDump1DComplex(	
	const char					*title,
	const FFTComplex			*buf,
	size_t						numSamples)	/* max to print */
{
	if(title) {
		printf("%s:\n", title);
	}
	size_t toPrint = MAX_TO_PRINT * MAX_TO_PRINT;
	if(toPrint > numSamples) {
		toPrint = numSamples;
	}
	for(size_t off=0; off<toPrint; off++) {
		printf("{%.2f,%.2f}, ", buf->real[off], buf->imag[off]);
		if((off % MAX_TO_PRINT) == (MAX_TO_PRINT - 1)) {
			printf("\n");
		}
	}
	printf("\n");
}

/* 
 * Dump 2-dimension FFTComplex data to stdout.
 */
extern void fftDump2DComplex(	
	const char					*title,
	const FFTComplex			*buf,
	size_t						numRows,		/* max to print */
	size_t						numCols)
{
	if(title) {
		printf("%s:\n", title);
	}
	size_t rowsToPrint = MAX_TO_PRINT;
	if(rowsToPrint > numRows) {
		rowsToPrint = numRows;
	}
	size_t colsToPrint = MAX_TO_PRINT;
	if(colsToPrint > numCols) {
		colsToPrint = numCols;
	}
	for(size_t row=0; row<rowsToPrint; row++) {
		for(size_t col=0; col<colsToPrint; col++) {
			size_t off = (row * numCols) + col;
			printf("{%.2f,%.2f}, ", buf->real[off], buf->imag[off]);
		}
		printf("\n");
	}
	printf("\n");
}


#endif	/* FFT_SPLIT_COMPLEX */
