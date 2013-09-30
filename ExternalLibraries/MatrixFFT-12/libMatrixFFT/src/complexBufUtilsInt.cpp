/*	File: complexBufUtilsInt.cpp 
	
	Description:
		FFTComplex buffer routines, interleaved version.
	
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
 * complexBufUtilsInt.cpp - FFTComplex buffer routines, interleaved version.
 *			  
 * Created 11/20/2008. 
 * Copyright 2008 by Apple, Inc. 
 */
 
#include <libMatrixFFT/MatrixFFTConfig.h>
#include <cstdio>
#include <cstring>

#if		!FFT_SPLIT_COMPLEX

#include <libMatrixFFT/complexBufUtils.h>
#include <libMatrixFFT/devRandom.h>
#include <libMatrixFFT/fftUtils.h>
#include "fftPriv.h"
#include <stdio.h>
#include <strings.h>

#pragma mark --- Buffer allocation ---

/*
 * Allocate an FFTComplex array with specified number of complex elements.
 * Returns NULL on malloc failure. Caller must free the returned
 * pointer via fftFreeComplexArray()> 
 */
FFTComplex *fftAllocComplexArray(
	size_t numComplex)
{
	size_t bufSize = numComplex * sizeof(FFTComplex);
	FFTComplex *fftComplex = (FFTComplex *)malloc(bufSize);
	if(fftComplex == NULL) {
		printf("***fftAllocComplexArray: malloc error\n");
	}
	return fftComplex;
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
	size_t bufSize = (numComplex * sizeof(FFTComplex)) + alignSize;	
	FFTComplex *freeComplex = (FFTComplex *)malloc(bufSize);
	if(freeComplex == NULL) {
		printf("***fftAllocComplexArrayAlign: malloc error\n");
	}
	*freePtr = freeComplex;
	return (FFTComplex *)FFT_ALIGN(freeComplex, alignSize);

}

void fftFreeComplexArray(
	FFTComplex *buf)
{
	COND_FREE(buf);
}

/* 
 * This config, fftFreeComplexArrayAlign is trivial since it's all one 
 * pointer...
 */
void fftFreeComplexArrayAlign(
	FFTComplex				*buf,
    FFTComplex              *freeBuf)
{
    fftFreeComplexArray(freeBuf);
}

/* copy from one FFTSplitComplex to another */
void fftCopyComplexArray(
	const FFTComplex	*src,
	FFTComplex			*dst,
	size_t				numSamples)
{
	size_t bufSize = numSamples * sizeof(FFTComplex);
	memmove(dst, src, bufSize);
}

/* Multiply each element in an FFTComplex by specified scale factor */
void fftScaleComplex(
	FFTComplex *buf,
	FFTFloat scaleFact,
	size_t numSamples)		/* in each of realp, imag */
{
	/* A bit of a kludge but I think it should be OK as long
	 * as the real part comes first in the FFTComplex struct */
	FFTvScale(&buf->real, scaleFact, numSamples << 1);
}

/* Return true if two FFTComplexes refer to same memory, else return false */
bool fftComplexEquiv(
    const FFTComplex        *buf1,
    const FFTComplex        *buf2)
{
    return (buf1 == buf2);
}

#pragma mark --- Copy between FFTComplex and vDSPComplex ---

void fftCopyToDSP(
	const FFTComplex *src,
	vDSPComplex *dst,
	size_t numComplex)
{
	FFTFloat *realp = dst->realp;
	FFTFloat *imagp = dst->imagp;
	
	for(size_t dex=0; dex<numComplex; dex++, src++) {
		*realp++ = src->real;
		*imagp++ = src->imag;
	}
}
	
void fftCopyFromDSP(
	const vDSPComplex *src,
	FFTComplex *dst,
	size_t numComplex)
{
	const FFTFloat *realp = src->realp;
	const FFTFloat *imagp = src->imagp;
	
	for(size_t dex=0; dex<numComplex; dex++, dst++) {
		dst->real = *realp++;
		dst->imag = *imagp++;
	}
}

#pragma mark --- Misc. Buffer manipulation ---

/* generate constant signal */
void genConstComplex(
	FFTComplex *buf,
	size_t numSamples,
	FFTFloat f)
{
    FFT_vfill(f, &buf->real, numSamples << 1);
}

/* Multiply a complex array by a vector of FFTFloats */
void fftMulComplexVect(
    FFTComplex *cbuf,           /* cbuf[n] *= fbuf[n] */
    const FFTFloat *fbuf,
    size_t N)
{
    for(size_t dex=0; dex<N; dex++) {
        FFTFloat f = *fbuf++;
        cbuf->real *= f;
        cbuf->imag *= f;
        cbuf++;
    }
}

/* flush contents of an FFTComplex from cache */
void fftFlushComplex(
	FFTComplex *buf,
	size_t numComplex)
{
    size_t numBytes = numComplex * sizeof(FFTComplex);
    fftFlushCache(buf, numBytes);
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
	double totalErr = 0.0;		/* running total of delta squared */
	double currMax = 0.0;		/* squared! */
	
	for(size_t dex=0; dex<numSamples; dex++) {
		double d = buf1->real - buf2->realp[dex];
		double deltaSquare = d * d;
		d = buf2->imagp[dex] - buf1->imag;
		deltaSquare += (d * d);
		if(deltaSquare > currMax) {
			currMax = deltaSquare;
		}
		totalErr += deltaSquare;
		buf1++;
	}
	*maxDelta = sqrt(currMax);
	*rmse = sqrt(totalErr / numSamples);
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
	double totalErr = 0.0;		/* running total of delta squared */
	double currMax = 0.0;		/* squared! */
	
	for(size_t dex=0; dex<numSamples/2; dex++) {
		double d = buf1->real - buf2->realp[dex];
		double deltaSquare = d * d;
		if(deltaSquare > currMax) {
			currMax = deltaSquare;
		}
		totalErr += deltaSquare;
		
		d = buf1->imag - buf2->imagp[dex];
		deltaSquare = d * d;
		if(deltaSquare > currMax) {
			currMax = deltaSquare;
		}
		totalErr += deltaSquare;
		buf1++;
	}
	*maxDelta = sqrt(currMax);
	*rmse = sqrt(totalErr / numSamples);
}

#pragma mark --- Buffer dump/display ---

#define MAX_TO_PRINT		8
#define MAX_ROWS_TO_PRINT	16

void fftDumpSquare(	
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
			printf("{%.2f,%.2f}, ", buf[off].real, buf[off].imag);
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
			printf("{%.2f,%.2f}, ", buf[off].real, buf[off].imag);
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
		printf("{%.2f,%.2f}, ", buf[off].real, buf[off].imag);
		if((off % MAX_TO_PRINT) == (MAX_TO_PRINT - 1)) {
			printf("\n");
		}
	}
	printf("\n");
}

/* 
 * Dump 2-dimension FFTComplex data to stdout.
 */
void fftDump2DComplex(	
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
			printf("{%.2f,%.2f}, ", buf[off].real, buf[off].imag);
		}
		printf("\n");
	}
	printf("\n");
}


#endif	/* FFT_SPLIT_COMPLEX */
