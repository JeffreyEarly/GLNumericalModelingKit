/*	File: complexChirpSignal.cpp 
	
	Description:
		Generate complex chirp signals, measure expected FFT outputs thereof.
	
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
 * complexChirpSignal.h - generate complex chirp signals, measure expected FFT outputs thereof. 
 *
 * Created 02/12/2009. 
 * Copyright 2009 by Apple, Inc. 
 */
 
#include <libMatrixFFT/MatrixFFT.h>
#include <libMatrixFFT/fftUtils.h>
#include <math.h>
#include <cstdio>
#include <libMatrixFFT/complexChirpSignal.h>
#include "PolyComplex.h"
#include <stdio.h>

/* callbacks occur this often for 1-D */
#define CALLBACK_INTERVAL_MASK_1D	((size_t)0xffff)

/* ditto for 2-D, used on the row counter */
#define CALLBACK_INTERVAL_MASK_2D	((size_t)0xff)

#pragma mark --- 1-D chirp signal functions ---

/* Generate 1-dimension "chirp" signal */
extern void fftGenChirp1D(
	FFTComplex		*buf,
	size_t			numSamples,			/* # of complex samples */
	fftCallbackFcn	callback,			/* optional */
	void			*callbackArg)		/* optional */
{
	double pi = M_PI;
	double oneOverN = 1.0 / (FFTFloat)numSamples;
	
	PolyComplex pc(buf);
	
	/* 
	 * All internal calculations in double precision - see comments in
	 * fftAnalyzeChirp1D(), below. 
	 */
	for(size_t j=0; j<numSamples; j++) {
		/* arg := pi * j^2 / N */
		double arg = j * j;
		arg *= pi;
		arg *= oneOverN;
		pc.real(cos(arg));
		pc.imag(sin(arg));
		++pc;
		if((callback != NULL) && ((j & CALLBACK_INTERVAL_MASK_1D) == 0)) {
			/* let client know how we're doing */
			unsigned percent = ((j + 1) * 100) / numSamples;
			callback(callbackArg, percent);
		}
	}
	if(callback != NULL) {
		callback(callbackArg, 100);
	}
}

/* 
 * Analyze round-trip FFT of 1-D chirp signal, comparing result
 * against expected values. Return max delta (mxe) and RMSE between
 * actual and expected. 
 */
void fftAnalyzeChirp1D(
	FFTComplex		*buf,
	size_t			numSamples,			/* # of complex samples */
	fftCallbackFcn	callback,			/* optional */
	void			*callbackArg,		/* optional */
	double			*mxe,				/* RETURNED */
	double			*rmse)				/* RETURNED */
{
	double pi = M_PI;
	double oneOverN = 1.0 / (FFTFloat)numSamples;
	double totalErr = 0.0;		/* running total of delta squared */
	double currMax = 0.0;		/* squared! */
	
	PolyComplex pc(buf);
	
	/*
	 * Do all calculations in double precision. You'd think that
	 *    float f = cos(arg)
	 * ...would give the same result as 
	 *    float f = cosf(arg)
	 *
	 * But no, they don't, and if we use the single-precision cos/sin 
	 * functions here, we get results that are off by several orders of magnitude
	 * for the single-precision configuration. 
	 */
	for(size_t j=0; j<numSamples; j++) {
		/* arg := pi * j^2 / N */
		double arg = j * j;
		arg *= pi;
		arg *= oneOverN;
		double r = cos(arg);
		double i = sin(arg);
		
		double d = pc.real() - r;
		double deltaSquare = d * d;
		d = pc.imag() - i;
		deltaSquare += (d * d);
		if(deltaSquare > currMax) {
			currMax = deltaSquare;
		}
		totalErr += deltaSquare;

		++pc;
		if((callback != NULL) && ((j & CALLBACK_INTERVAL_MASK_1D) == 0)) {
			/* let client know how we're doing */
			unsigned percent = ((j + 1) * 100) / numSamples;
			callback(callbackArg, percent);
		}
	}
	if(callback != NULL) {
		callback(callbackArg, 100);
	}
	*mxe = sqrt(currMax);
	*rmse = sqrt(totalErr / numSamples);
}

#pragma mark --- 2-D chirp signal functions ---

/* Generate 2-dimension "chirp" signal */
extern void fftGenChirp2D(
	FFTComplex		*buf,
	size_t			width,				/* in complex samples */
	size_t			height,				/* # of rows */
	fftCallbackFcn	callback,			/* optional */
	void			*callbackArg)		/* optional */
{
	double pi = M_PI;
	PolyComplex pc(buf);
	size_t totalSamples = width * height;
	
	for(size_t k=0; k<height; k++) {
		for(size_t j=0; j<width; j++) {
			/* arg := pi * (j^2/W + k^2/H) */
			double arg = ((double)(j * j)) / (double)width;
			arg += ((double)(k * k)) / (double)height;
			arg *= pi;
			pc.real(cos(arg));
			pc.imag(sin(arg));
			++pc;
		}
		if((callback != NULL) && ((k & CALLBACK_INTERVAL_MASK_2D) == 0)) {
			/* let client know how we're doing */
			unsigned percent = ((k + 1) * 100) / totalSamples;
			callback(callbackArg, percent);
		}
	}
	if(callback != NULL) {
		callback(callbackArg, 100);
	}
}

/* 
 * Analyze round-trip FFT of 2-D chirp signal, comparing result
 * against expected values. Return max delta (mxe) and RMSE between
 * actual and expected. 
 */
extern void fftAnalyzeChirp2D(
	FFTComplex		*buf,
	size_t			width,				/* in complex samples */
	size_t			height,				/* # of rows */
	fftCallbackFcn	callback,			/* optional */
	void			*callbackArg,		/* optional */
	double			*mxe,				/* RETURNED */
	double			*rmse)				/* RETURNED */
{
	double pi = M_PI;
	PolyComplex pc(buf);
	size_t totalSamples = width * height;
	double totalErr = 0.0;		/* running total of delta squared */
	double currMax = 0.0;		/* squared! */
	
	for(size_t k=0; k<height; k++) {
		for(size_t j=0; j<width; j++) {
			/* arg := pi * (j^2/W + k^2/H) */
			double arg = ((double)(j * j)) / (double)width;
			arg += ((double)(k * k)) / (double)height;
			arg *= pi;
			double r = cos(arg);
			double i = sin(arg);

			double d = pc.real() - r;
			double deltaSquare = d * d;
			d = pc.imag() - i;
			deltaSquare += (d * d);
			if(deltaSquare > currMax) {
				currMax = deltaSquare;
			}
			totalErr += deltaSquare;

			++pc;
		}
		if((callback != NULL) && ((k & CALLBACK_INTERVAL_MASK_2D) == 0)) {
			/* let client know how we're doing */
			unsigned percent = ((k + 1) * 100) / totalSamples;
			callback(callbackArg, percent);
		}
	}
	if(callback != NULL) {
		callback(callbackArg, 100);
	}
	*mxe = sqrt(currMax);
	*rmse = sqrt(totalErr / totalSamples);
}
#pragma mark --- Analyze forward FFT of chirp signal ---

#define DEBUG_CHE	0

/* 
 * Return max error of the forward FFT of a chirp signal.
 * Used for 1-D and 2-D.
 *
 * The error is defined as CHE := max(k) |(|X[k]| / N^0.5) - 1|
 */
double fftMaxChirpError(
	const FFTComplex	*buf,
	size_t				numSamples,			/* # of complex samples */
	fftCallbackFcn		callback,			/* optional */
	void				*callbackArg)		/* optional */
{
	unsigned n;
	if(!fftIsPowerOfTwo(numSamples, &n)) {
		printf("***fftMaxChirpError: size not a power of two\n");
		return 0.0;
	}
	double oneOverN = 1.0 / (double)numSamples;
	
	double maxSqErr = 0.0;
	FFTFloat ampSqAtMax = 0.0;
	PolyComplex pc((FFTComplex *)buf);
	size_t dexAtMax = 0;
	
	for(size_t dex=0; dex<numSamples; dex++) {
		double r = pc.real();
		double i = pc.imag();
		++pc;
		double ampSq = (r*r + i*i);
		double err = fabs((ampSq * oneOverN) - 1.0);
		if(err > maxSqErr) {
			maxSqErr = err;
			ampSqAtMax = ampSq;
			dexAtMax = dex;
		}
		if((callback != NULL) && ((dex & CALLBACK_INTERVAL_MASK_1D) == 0)) {
			/* let client know how we're doing */
			unsigned percent = ((dex + 1) * 100) / numSamples;
			callback(callbackArg, percent);
		}
	}
	if(callback != NULL) {
		callback(callbackArg, 100);
	}
	
	double sqRtN = sqrt((double)numSamples);
	double ampl  = sqrt(ampSqAtMax);
	double drtn  = fabs((ampl / sqRtN) - 1.0);
	
	#if DEBUG_CHE
	
	pc.set((FFTComplex *)buf, dexAtMax);
	printf("max err at index %llu\n", (unsigned long long)dexAtMax);
	printf("real             %.8f\n", pc.real());
	printf("imag             %.8f\n", pc.imag());
	printf("ampl^2           %.8f\n", ampSqAtMax);
	printf("ampl             %.1f\n", ampl);
	printf("N                %llu\n", (unsigned long long)numSamples);
	printf("sqrt(N)          %.1f\n", sqRtN);
	printf("ampl / sqrt(N)   %.8f\n", ampl / sqRtN);
	printf("CHE              %.8f\n", drtn);
	#endif	/* DEBUG_CHE */

	return drtn;
}
