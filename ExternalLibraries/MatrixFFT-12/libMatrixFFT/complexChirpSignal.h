/*	File: complexChirpSignal.h 
	
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
 
#ifndef	_COMPLEX_CHIRP_SIGNAL_H_
#define _COMPLEX_CHIRP_SIGNAL_H_

#include <libMatrixFFT/complexBuf.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Callback function for progress reporting. The percentDone value is 
 * between 0 and 100, with 100 indicating "done". The 'arg' value is 
 * whatever is was passed along with the callback pointer to 
 * fftGenReal1DTestSignal() or fftGenReal1DTestSignalFFT().
 */
typedef void (*fftCallbackFcn)(void *arg, unsigned percentDone);

#pragma mark --- 1-D chirp signal functions ---

/* 
 * Generate 1-dimension "chirp" signal.
 *
 * The signal is:
 *
 *		x[j] = e^(i pi j^2/n)
 *
 * Or
 *		x[j] = cos(pi j^2/n) + i sin(pi j^2/n)
 *
 * Where n = numSamples. 
 */
extern void fftGenChirp1D(
	FFTComplex		*buf,
	size_t			numSamples,			/* # of complex samples */
	fftCallbackFcn	callback,			/* optional */
	void			*callbackArg);		/* optional */

/* 
 * Analyze round-trip FFT of 1-D chirp signal, comparing result
 * against expected values. Return max delta (mxe) and RMSE between
 * actual and expected. 
 */
extern void fftAnalyzeChirp1D(
	FFTComplex		*buf,
	size_t			numSamples,			/* # of complex samples */
	fftCallbackFcn	callback,			/* optional */
	void			*callbackArg,		/* optional */
	double			*mxe,				/* RETURNED */
	double			*rmse);				/* RETURNED */


#pragma mark --- 2-D chirp signal functions ---

/* 
 * Generate the 2-dimension "chirp test" signal.
 *
 * The signal is:
 *
 *		x[k][j] = e ^ (i pi (j^2/W + k^2/H))
 *
 * Or
 *		x[k][j] = cos(pi * (j^2/W + k^2/H) + i sin(pi * (j^2/W + k^2/H))
 *
 * Where W = # columns and H = # rows. 
 */
extern void fftGenChirp2D(
	FFTComplex		*buf,
	size_t			width,				/* in complex samples */
	size_t			height,				/* # of rows */
	fftCallbackFcn	callback,			/* optional */
	void			*callbackArg);		/* optional */


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
	double			*rmse);				/* RETURNED */

#pragma mark --- Analyze forward FFT of chirp signal ---

/* 
 * Return max error of the forward FFT of a chirp signal.
 * Used for 1-D and 2-D.
 *
 * The error is defined as CHE := max(k) |(|X[k]| / N^0.5) - 1|
 */
extern double fftMaxChirpError(
	const FFTComplex	*buf,
	size_t				numSamples,			/* # of complex samples */
	fftCallbackFcn		callback,			/* optional */
	void				*callbackArg);		/* optional */


#ifdef __cplusplus
}
#endif

#endif	/* _COMPLEX_CHIRP_SIGNAL_H_ */
