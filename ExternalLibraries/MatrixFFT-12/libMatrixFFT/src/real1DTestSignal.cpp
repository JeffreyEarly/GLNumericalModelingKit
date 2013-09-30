/*	File: real1DTestSignal.cpp 
	
	Description:
		real 1-D FFT test signal functions.
	
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
 * real1DTestSignal.cpp - real 1-D FFT test signal functions.
 *
 * This is an implementation of the "Real-signal test" algorithms described 
 * in section 7, "Error Analysis", of "Large-scale FFTs and convolutions on 
 * Apple hardware", available here:
 * http://images.apple.com/acg/pdf/FFTapps.pdf
 *
 * Created 7/21/2008. 
 * Copyright 2008, 2009 by Apple, Inc. 
 */

#include <stdlib.h>
#include <strings.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <libMatrixFFT/fftUtils.h>
#include <libMatrixFFT/real1DTestSignal.h>
#include <libMatrixFFT/MatrixFFT.h>
#include <libMatrixFFT/complexBufUtils.h>
#include "fftSinCos.h"
#include "PolyComplex.h"

#define TEST_SIGNAL_M	((double)0.05)
#define TEST_SIGNAL_2M	(2.0 * TEST_SIGNAL_M)

/* callbacks occur this often */
#define CALLBACK_INTERVAL_MASK		((size_t)0xffff)

/* ditto for 2-D-oriented loops */
#define CALLBACK_INTERVAL_2D_MASK	((size_t)0xff)

#pragma mark --- 1-D Real test signal ---

/* 
 * Generate 1-D real test signal.
 *
 * x[j] = e^(-Mj) - 2e^(-2Mj), j = 0...N-1
 */

/* core calculation, given j */
static void genTestSig(
	size_t j,
	FFTFloat *realp,
	FFTFloat *imagp)
{
	/* mMj = -Mj */
	double mMj = -(TEST_SIGNAL_M * (double)j);
	*realp = exp(mMj) - (2 * exp(2.0 * mMj));
	/* mMj = -M(j+1) */
	mMj -= TEST_SIGNAL_M;
	*imagp = exp(mMj) - (2 * exp(2.0 * mMj));
}

/* Public function to generate test signal */
int fftGenReal1DTestSignal(
	FFTComplex		*buf,
	size_t			numSamples,			/* number of real samples */
	fftCallbackFcn	callback,			/* optional */
	void			*callbackArg)		/* optional */
{	
	size_t complexSamples = numSamples >> 1;
	size_t j=0;
	PolyComplex pc(buf);
	
	/*
	 * pc.real contains even indices of real signal 
	 * pc.imag contains odd indices
	 */
	for(size_t dex=0; dex<complexSamples; dex++, j+=2) {
		genTestSig(j, pc.realP(), pc.imagP());
		++pc;
		if((callback != NULL) && ((dex & CALLBACK_INTERVAL_MASK) == 0)) {
			/* let client know how we're doing */
			unsigned percent = ((dex + 1) * 100) / complexSamples;
			callback(callbackArg, percent);
		}
	}
	if(callback != NULL) {
		callback(callbackArg, 100);
	}
	return 0;
}

/*
 * Analyze a 1-D real test signal in row order, comparing it against
 * the expected values. Return max delta (mxe) and RMSE between
 * actual and expected. 
 */
int fftAnalyzeReal1DTestSignal(
	FFTComplex		*buf,
	size_t			numSamples,			/* number of real samples */
	fftCallbackFcn	callback,			/* optional */
	void			*callbackArg,		/* optional */
	double			*mxe,				/* RETURNED */
	double			*rmse)				/* RETURNED */
{
	size_t complexSamples = numSamples >> 1;
	size_t j=0;
	PolyComplex pc(buf);
	double totalErr = 0.0;		/* running total of delta squared */
	double currMax = 0.0;		/* squared! */
	
	/*
	 * pc.real contains even indices of real signal 
	 * pc.imag contains odd indices
	 */
	for(size_t dex=0; dex<complexSamples; dex++, j+=2) {
		FFTFloat r, i;
		genTestSig(j, &r, &i);
		
		double d = pc.real() - r;
		double deltaSquare = d * d;
		if(deltaSquare > currMax) {
			currMax = deltaSquare;
		}
		totalErr += deltaSquare;
		
		d = pc.imag() - i;
		deltaSquare = d * d;
		if(deltaSquare > currMax) {
			currMax = deltaSquare;
		}
		totalErr += deltaSquare;
		++pc;

		if((callback != NULL) && ((dex & CALLBACK_INTERVAL_MASK) == 0)) {
			/* let client know how we're doing */
			unsigned percent = ((dex + 1) * 100) / complexSamples;
			callback(callbackArg, percent);
		}
	}
	
	*mxe  = sqrt(currMax);
	*rmse = sqrt(totalErr / numSamples);

	if(callback != NULL) {
		callback(callbackArg, 100);
	}
	return 0;
}

#pragma mark --- FFT output of 1-D test signal ---

/* 
 * Generate expected forward FFT of 1-D test signal.
 *
 * X[k] = F(M, k) - 2F(2M, k)
 *
 * Where
 *
 * F(L, k) = (1 - e^(-LN)) * (1 - (e^(-L) * e^(2*pi*i*k/N)))
 *           -----------------------------------------------
 *                1 - (2e^(-L) * cos(2*k*pi/N)) + e^(-2L)
 *
 * We scale the "expected" FFT here by a factor of 2 to match the output 
 * of both our own 1-D FFT and vDSP. 
 */
 
/* F(K,k) */
static void testSignalF(
	MatrixFFTPlan	mfftPlan,
	double			N,
	double			L,
	size_t			k,
	double			*r,		/* real, RETURNED */
	double			*i)		/* imaginary, RETURNED */
{
	/*
	 * numer := (e^(-L) * e^(2*pi*i*k/N))
	 *        = (e^(-L) * (cos(2*pi*k/N) + i * sin(2*pi*k/N))
	 *
	 * It would be nice if we could always do double-precision sin/cos here but 
	 * we really want to use the generated tables....
	 */ 
	FFTFloat f1Sin, f1Cos;
	fftSinCos(mfftPlan, k, &f1Cos, &f1Sin);
	double numerR, numerI;		/* numer */
	double mL = -L;
	double expML = exp(mL);		/* e^(-L) */
	numerR = expML * f1Cos;
	numerI = expML * f1Sin;
	
	/* 
	 * numer := 1 - numer 
	 *        = (1 - (e^(-L) * e^(2*pi*i*k/N)))
	 */
	numerR = 1.0 - numerR;
	numerI = -numerI;
	
	/* 
	 * numer := (1 - e^(-LN)) * numer 
	 *        = (1 - e^(-LN)) * (1 - (e^(-L) * e^(2*pi*i*k/N))) 
	 */
	double tmp = 1 - exp(mL * N);
	numerR *= tmp;
	numerI *= tmp;
	
	/*
	 * denom := 1 - (2e^(-L) * cos(2*k*pi)) + e^(-2L)
	 */
	double denom = 1.0 - (2.0 * expML * f1Cos) + exp(2.0 * mL);
	
	/* 
	 * result = numer / denom
	 */
	*r = numerR / denom;
	*i = numerI / denom;
}

/* X(k) = F(M, k) - 2F(2M, k) */
static void testSignalX(
	MatrixFFTPlan	mfftPlan,
	double			N,
	size_t			k,
	FFTFloat		*r,		/* real, RETURNED */
	FFTFloat		*i)		/* imaginary, RETURNED */
{
	double r1, r2, i1, i2;

	testSignalF(mfftPlan, N, TEST_SIGNAL_M, k, &r1, &i1);
	testSignalF(mfftPlan, N, TEST_SIGNAL_2M, k, &r2, &i2);
	
	/* scale factor injected here */
	*r = 2.0 * (r1 - (2.0 * r2));
	*i = 2.0 * (i1 - (2.0 * i2));
}

/* 
 * Produce output in column order with rows & columns swapped, matching our library's 
 * native output of 1-D real FFT 
 */
static int testSignal1D_column(
	MatrixFFTPlan	mfftPlan,
	FFTComplex		*buf,
	size_t			numSamples,			/* number of real samples */
	fftCallbackFcn	callback,			/* optional */
	void			*callbackArg)		/* optional */
{
	unsigned n;	
	if(!fftIsPowerOfTwo(numSamples, &n)) {
		printf("***testSignal1D_column: not a power of two\n");
		return -1;
	}
	
	/* rows and columns in real elements */
	size_t numRows;
	size_t numCols;
	mfftRectangle(mfftPlan, &numRows, &numCols);
	
	/* in complex */
	numCols >>= 1;
	

	/* 
	 * We're going thru the destination buffer in consecutive row order to 
	 * maximize memory performance. The pc pointers just increment on each op,
	 * unlike the value of k. 
	 */
	PolyComplex pc(buf, 1);
	double floatN = (double)numSamples;
	
	size_t k;
	for(size_t row=0; row<numRows; row++) {
		size_t startCol = 0;
		
		if(row == 0) {
			/* 
			 * caller already handled [0,0], the DC/Nyquist special case. 
			 * We start at column 1.
			 */
			startCol = 1;
			k = numRows;
		}
		else {
			k = row;
		}
		for(size_t col=startCol; col<numCols; col++, k+=numRows) {
			testSignalX(mfftPlan, floatN, k, pc.realP(), pc.imagP());
			++pc;
		}
		if((callback != NULL) && ((row & CALLBACK_INTERVAL_2D_MASK) == 0)) {
			/* let client know how we're doing */
			unsigned percent = ((row + 1) * 100) / numRows;
			callback(callbackArg, percent);
		}
	}
	if(callback != NULL) {
		callback(callbackArg, 100);
	}
	return 0;
}


/* public function */
int fftGenReal1DTestSignalFFT(
	MatrixFFTPlan	mfftPlan,
	FFTComplex		*buf,
	size_t			numSamples,			/* number of real samples */
	SignalFormat1D	format,				/* row, column, custom */
	fftCallbackFcn	callback,			/* optional */
	void			*callbackArg)		/* optional */
{	
	FFTFloat r, i;
	size_t numComplex = numSamples >> 1;
	double floatN = (FFTFloat)numSamples;
	
	PolyComplex pc(buf);
	
	/* DC portion - just the real at k=0 --> buf->real[0] */
	testSignalX(mfftPlan, floatN, 0, &r, &i);
	pc.real(r);
	
	/* Nyquist portion - just the real at k=N/2 --> buf->imag[0] */
	testSignalX(mfftPlan, floatN, numComplex, &r, &i);
	pc.imag(r);
	
	switch(format) {
		case ST1D_RowOrder:
			break;				/* we handle this below */
		case ST1D_Column:
			return testSignal1D_column(mfftPlan, buf, numSamples, callback, callbackArg);
	}
	
	/* produce output in normal row order */
	
	++pc;	
	for(size_t k=1; k<numComplex; k++) {
		testSignalX(mfftPlan, floatN, k, pc.realP(), pc.imagP());
		++pc;
		if((callback != NULL) && ((k & CALLBACK_INTERVAL_MASK) == 0)) {
			/* let client know how we're doing */
			unsigned percent = ((k + 1) * 100) / numComplex;
			callback(callbackArg, percent);
		}
	}
	if(callback != NULL) {
		callback(callbackArg, 100);
	}
	return 0;
}
	
#pragma mark ---- Analyze FFT of 1-D real test signal ----

/* signal in column order, matching our library's native output of 1-D real FFT */
static int analyzeTestSignal1D_FFTcolumn(
	MatrixFFTPlan	mfftPlan,
	FFTComplex		*buf,
	size_t			numSamples,			/* number of real samples */
	fftCallbackFcn	callback,			/* optional */
	void			*callbackArg,		/* optional */
	double			*mxe,				/* RETURNED */
	double			*rmse)				/* RETURNED */
{
	unsigned n;	
	if(!fftIsPowerOfTwo(numSamples, &n)) {
		printf("***testSignal1D_column: not a power of two\n");
		return -1;
	}
	
	/* rows and columns in real elements */
	size_t numRows;
	size_t numCols;
	mfftRectangle(mfftPlan, &numRows, &numCols);
	
	/* in complex */
	numCols >>= 1;
	
	/* 
	 * We're going thru the destination buffer in consecutive row order to 
	 * maximize memory performance. The pc pointers just increment on each op,
	 * unlike the value of k. 
	 */
	PolyComplex pc(buf);
	double floatN = (double)numSamples;

	double totalErr = 0.0;		/* running total of delta squared */
	double currMax = 0.0;		/* squared! */
	
	FFTFloat r, i, rDC, rNY;
	/* DC portion - just the real at k=0 - buf->real[0] */
	testSignalX(mfftPlan, floatN, 0, &rDC, &i);
	/* Nyquist portion - just the real at k=N/2 --> buf->imag[0] */
	testSignalX(mfftPlan, floatN, numSamples >> 1, &rNY, &i);

	// printf("rDC %.2f rNY %.2f pc.real %.2f pc.imag %.2f\n", rDC, rNY, pc.real(), pc.imag());

	double d = pc.real() - rDC;
	double deltaSquare = d * d;
	d = pc.imag() - rNY;
	deltaSquare += (d * d);
	if(deltaSquare > currMax) {
		currMax = deltaSquare;
	}
	totalErr += deltaSquare;
	
	++pc;
	
	for(size_t row=0; row<numRows; row++) {
		size_t k;
		size_t startCol = 0;
		
		if(row == 0) {
			/* 
			 * We already handled [0,0], the DC/Nyquist special case. 
			 * We start at column 1.
			 */
			startCol = 1;
			k = numRows;
		}
		else {
			k = row;
		}

		for(size_t col=startCol; col<numCols; col++, k+=numRows) {
			testSignalX(mfftPlan, floatN, k, &r, &i);
			
			//printf("r %.2f i %.2f pc.real %.2f pc.imag %.2f\n", r, i, pc.real(), pc.imag());

			d = pc.real() - r;
			deltaSquare = d * d;
			d = pc.imag() - i;
			deltaSquare += (d * d);
			if(deltaSquare > currMax) {
				currMax = deltaSquare;
			}
			totalErr += deltaSquare;

			++pc;
		}
		if((callback != NULL) && ((row & CALLBACK_INTERVAL_2D_MASK) == 0)) {
			/* let client know how we're doing */
			unsigned percent = ((row + 1) * 100) / numRows;
			callback(callbackArg, percent);
		}
	}
	
	*mxe  = sqrt(currMax);
	*rmse = sqrt(totalErr / numSamples);

	if(callback != NULL) {
		callback(callbackArg, 100);
	}
	return 0;
}


/* public function */
int fftAnalyzeReal1DTestSignalFFT(
	MatrixFFTPlan	mfftPlan,
	FFTComplex		*buf,
	size_t			numSamples,			/* number of real samples */
	SignalFormat1D	format,				/* row, column, custom */
	fftCallbackFcn	callback,			/* optional */
	void			*callbackArg,		/* optional */
	double			*mxe,				/* RETURNED */
	double			*rmse)				/* RETURNED */
{	
	size_t numComplex = numSamples >> 1;
	double floatN = (FFTFloat)numSamples;
	
	PolyComplex pc(buf);
	
	switch(format) {
		case ST1D_RowOrder:
			break;				/* we handle this below */
		case ST1D_Column:
			return analyzeTestSignal1D_FFTcolumn(mfftPlan, buf, numSamples, callback, callbackArg,
				mxe, rmse);
	}
	
	/* signal is in normal row order */

	double totalErr = 0.0;		/* running total of delta squared */
	double currMax = 0.0;		/* squared! */
	
	FFTFloat r, i, rDC, rNY;

	/* DC portion - just the real at k=0 - buf->real[0] */
	testSignalX(mfftPlan, floatN, 0, &rDC, &i);
	/* Nyquist portion - just the real at k=N/2 --> buf->imag[0] */
	testSignalX(mfftPlan, floatN, numComplex, &rNY, &i);

	double d = pc.real() - rDC;
	double deltaSquare = d * d;
	d = pc.imag() - rNY;
	deltaSquare += (d * d);
	if(deltaSquare > currMax) {
		currMax = deltaSquare;
	}
	totalErr += deltaSquare;
	
	++pc;
	
	for(size_t k=1; k<numComplex; k++) {
		testSignalX(mfftPlan, floatN, k, &r, &i);
		
		d = pc.real() - r;
		deltaSquare = d * d;
		d = pc.imag() - i;
		deltaSquare += (d * d);
		if(deltaSquare > currMax) {
			currMax = deltaSquare;
		}
		totalErr += deltaSquare;
		if((callback != NULL) && ((k & CALLBACK_INTERVAL_MASK) == 0)) {
			/* let client know how we're doing */
			unsigned percent = ((k + 1) * 100) / numComplex;
			callback(callbackArg, percent);
		}
		++pc;
	}
	
	*mxe  = sqrt(currMax);
	*rmse = sqrt(totalErr / numSamples);

	if(callback != NULL) {
		callback(callbackArg, 100);
	}
	return 0;
}
