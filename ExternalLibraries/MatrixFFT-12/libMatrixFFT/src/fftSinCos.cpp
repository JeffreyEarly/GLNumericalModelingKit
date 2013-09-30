/*	File: fftSinCos.cpp  
	
	Description:
		Sine/cosine related functions for MatrixFFT module.
	
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
 * fftSinCos.cpp - sine/cosine related functions for MatrixFFT module.
 *			  
 * Created 9/25/2008. 
 * Copyright 2008 by Apple, Inc. 
 */
 
#include "fftSinCos.h"
#include "fftPriv.h"
#include "fftDebug.h"
#include <libMatrixFFT/fftUtils.h>
#include <math.h>

#define SIN_COS_DEBUG	0	/* dump all sin/cos values */
#if		SIN_COS_DEBUG
#define scprint(args...)	printf(args)
#else
#define scprint(args...)
#endif	/* SIN_COS_DEBUG */

/* 
 * Set up normal, fully-populated sin table.
 * Returns nonzero on error.
 */
int fftSetupSineNorm(
	MatrixFFTPlan mfftPlan)
{
	/* 
	 * Set up sine lookup table. Cosine is derived from via cos(x) = sin(x + pi/2).
	 * Size of lookup table = (pi/2 + 1) 
	 */
	size_t mallocSize = (size_t)sizeof(FFTFloat) * (mfftPlan->piOver2 + 1);
	mfftPlan->sine = (FFTFloat *)malloc(mallocSize);
	if(mfftPlan->sine == NULL) {
		printf("***fftRowCreateSetupCommon: malloc failure (2)\n");
		return -1;
	}
	
	/* 
	 * sin(2 * pi * x / N) for x = 0..(N/4)
	 * Note we DO include the value for N/4 or pi/2...
	 */
	FFTFloat *sinP = mfftPlan->sine;
	FFTFloat twoPiOverN = 2.0 * M_PI / (FFTFloat)mfftPlan->N;
	
	for(size_t x=0; x<=mfftPlan->piOver2; x++) {
		FFTFloat xn = (FFTFloat)x * twoPiOverN;
		*sinP++ = FFTSin(xn);
	}
	return 0;
}

/*
 * Optimized sine lookup table when using "angle incremnent" sine/cos calculation. 
 * At run time, for each row, the caller needs sin and cos values for the following
 * columns. AngleReset is the frequency at which sin and cos are "reset" to actual
 * lookup-table values. VectorSize is the number of elements in a vector. 
 *
 * 0, 1, ... vectorSize 
 * <then skip to...> 
 * (angleReset*vectorSize) ... (angleReset*vectorSize + vectorSize) 
 * <skip>
 * (2*angleReset*vectorSize) ... (2*angleReset*vectorSize + vectorSize)
 *
 * Actually, all of the elements of row 0 and column 0 have the same sin, namely 0, 
 * and the same cos, namely 1. So we skip the entire first row and column 0
 * of each subsequent row and handle those as special cases at the top of 
 * fftSinCosOpt(). Then for each row, we have vectorSize entries, each angleReset
 * columns, for a total of ceil(rowSize / angleReset) * vectorSize entries per
 * row (for row > 0). 
 * 
 * The sinPeriod value must be a power of two since we need to perform shifts
 * based on it. 
 */
int fftSetupSineOpt(
	MatrixFFTPlan mfftPlan)
{
	RFASSERT(mfftPlan->nRow >= 2);
	RFASSERT(mfftPlan->N != 0);
	if(!fftIsPowerOfTwo(mfftPlan->sinPeriod, &mfftPlan->sinPeriodShift)) {
		printf("***fftRowSetupSineOpt: sinPeriod must be a power of 2\n");
		return -1;
	}
	
	size_t numCols = mfftPlan->numCols;
	size_t numRows = mfftPlan->numRows;
	size_t vectorsPerRow = (numCols + mfftPlan->sinPeriod - 1) / mfftPlan->sinPeriod;
	mfftPlan->valuesPerRow = vectorsPerRow * FFT_FLOATS_PER_VECTOR;
	
	/* we don't allocate row 0 */
	size_t mallocSize = (size_t)sizeof(FFTFloat) * (numRows - 1) * mfftPlan->valuesPerRow;
	mdprint("...allocating %lu bytes for sin/cos\n", (unsigned long)(mallocSize * 2));

	mfftPlan->sine = (FFTFloat *)malloc(mallocSize);
	if(mfftPlan->sine == NULL) {
		printf("***fftRowSetupSineOpt: malloc failure\n");
		return -1;
	}
	mfftPlan->cosine = (FFTFloat *)malloc(mallocSize);
	if(mfftPlan->cosine == NULL) {
		printf("***fftRowSetupSineOpt: malloc failure\n");
		free(mfftPlan->sine);
		mfftPlan->sine = NULL;
		return -1;
	}
	FFTFloat *sinP = mfftPlan->sine;
	FFTFloat *cosP = mfftPlan->cosine;
	FFTFloat twoPiOverN = 2.0 * M_PI / (FFTFloat)mfftPlan->N;
	
	for(size_t row=1; row<numRows; row++) {
		/* first, col 1..FFT_FLOATS_PER_VECTOR */
		size_t col;
		size_t rc;
		for(col=1; col<=FFT_FLOATS_PER_VECTOR; col++) {
			rc = row * col;
			FFTFloat rcn = (FFTFloat)rc * twoPiOverN;
			*sinP++ = FFTSin(rcn);
			*cosP++ = FFTCos(rcn);
		}
		
		/* then, FFT_FLOATS_PER_VECTOR values each sinPeriod */
		col = mfftPlan->sinPeriod;
		while((col + FFT_FLOATS_PER_VECTOR) < numCols) {
			size_t thisCol = col;
			for(unsigned dex=0; dex<FFT_FLOATS_PER_VECTOR; dex++) {
				rc = row * thisCol++;
				FFTFloat rcn = (FFTFloat)rc * twoPiOverN;
				*sinP++ = FFTSin(rcn);
				*cosP++ = FFTCos(rcn);
			}
			col += mfftPlan->sinPeriod;
		}
	}
	return 0;
}


/* 
 * Assuming precalculated tables (via fftSetupSineNorm()), obtain
 * sin and cos of (2 * pi * a * b / N).
 *
 * fftSine() breaks up the 2*pi range (0..N) into 4 quadrants, each of which 
 * is expressed in terms of the lookup table's values which only represent
 * the first quadrant, 0..pi/2. 
 * Note that our table *does* have an entry for pi/2. 
 * 0 <= ab < N on entry.
 */
static FFT_INLINE FFTFloat fftSine(
	MatrixFFTPlan mfftPlan,
	size_t ab)
{
	if(ab <= mfftPlan->piOver2) {
		/* quad 0: raw table value */
		return mfftPlan->sine[ab];
	}
	else if(ab <= mfftPlan->pi) {
		/* quad 1: reflect around the line x = pi/2 */
		return mfftPlan->sine[mfftPlan->pi - ab];
	}
	/* remaining two quads are the additive inverse of the previous two */
	else if(ab <= mfftPlan->threePiOver2) {
		/* value = -sin(ab - pi) */
		return -(mfftPlan->sine[ab - mfftPlan->pi]);
	}
	else {
		return -(mfftPlan->sine[mfftPlan->N - ab]);
	}
}

/* 
 * Assuming precalculated tables (via fftSetupSineNorm()), obtain
 * sine and/or cosine of (2 * pi * ab / N).
 */
void fftSinCos(
	MatrixFFTPlan mfftPlan,
	size_t ab,
	FFTFloat *cosv,		/* optionally RETURNED */
	FFTFloat *sinv)		/* optionally RETURNED */
{
	RFASSERT(mfftPlan->sinPeriod == 0);
	
	if(mfftPlan->sineShift) {
		RFASSERT(mfftPlan->sinTableType == STT_External);
		/* caller's data is smaller than N */
		ab <<= mfftPlan->sineShift;
	}
	if(mfftPlan->sinTableType == STT_External) {
		/* switch over to owner's table */
		RFASSERT(mfftPlan->externSineTable != NULL);
		mfftPlan = mfftPlan->externSineTable;
	}
	/* how could ab ever be larger than N? */
	RFASSERT(ab < mfftPlan->N);
	size_t dex = ab & mfftPlan->sinCosNMask;

	if(sinv) {
		*sinv = fftSine(mfftPlan, dex);
		scprint("sin(%lu) = %f\n", (unsigned long)ab, *sinv);
	}
	
	if(cosv) {
		/* cos(x) = sin(x + pi/2) */
		*cosv = fftSine(mfftPlan, (dex + mfftPlan->piOver2) & mfftPlan->sinCosNMask);
		//scprint("cos(%lu) = %f\n", (unsigned long)ab, *cosv);
	}
}

/* 
 * Perform lookup of sin/cos using optimized lookup tables.
 */
void fftSinCosOpt(
	MatrixFFTPlan mfftPlan,
	size_t row,
	size_t col,
	FFTFloat *cosv,		/* optionally RETURNED */
	FFTFloat *sinv)		/* optionally RETURNED */
{
	RFASSERT(mfftPlan->sinPeriod != 0);
	
	/* shift not supported */
	RFASSERT(mfftPlan->sineShift == 0);
	if(mfftPlan->sineShift) {
		printf("***fftSinCosOpt: sineShift not supported\n");
		return;
	}
	
	if((row == 0) || (col == 0)) {
		if(sinv) {
			*sinv = 0;
			scprint("sin(%lu,%lu) = 0\n", (unsigned long)row, (unsigned long)col);
		}
		if(cosv) {
			*cosv = 1.0;
			//scprint("cos(%lu,%lu) = 1.0\n", (unsigned long)row, (unsigned long)col);
		}
		return;
	}
	
	/* tables start at row 1 */
	size_t rowOffset = ((row - 1) * mfftPlan->valuesPerRow);
	
	/* first FFT_FLOATS_PER_VECTOR values are sin/cos[1..FFT_FLOATS_PER_VECTOR] */
	size_t colOffset = 0;
	if(col <= FFT_FLOATS_PER_VECTOR) {
		colOffset = col - 1;
	}
	else {
		/* then, FFT_FLOATS_PER_VECTOR values each sinPeriod */
		size_t whichVect = col >> mfftPlan->sinPeriodShift;
		size_t vectStart = whichVect << mfftPlan->sinPeriodShift;
		colOffset = col - vectStart;
		
		/* the crucial assertion: we're being asked for properly aligned column values */
		RFASSERT(colOffset < FFT_FLOATS_PER_VECTOR);
		
		/* now index into the sparse table */
		colOffset += (whichVect * FFT_FLOATS_PER_VECTOR);
	}

	size_t totalOffset = rowOffset + colOffset;
	if(sinv) {
		*sinv = mfftPlan->sine[totalOffset];
		scprint("sin(%lu,%lu) = %f\n", (unsigned long)row, (unsigned long)col, *sinv);
	}
	
	if(cosv) {
		*cosv = mfftPlan->cosine[totalOffset];
		//scprint("cos(%lu,%lu) = %f\n", (unsigned long)row, (unsigned long)col, *cosv);
	}

}
