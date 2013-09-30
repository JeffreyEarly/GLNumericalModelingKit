/*	File: NUFFT_Utils.h 
	
	Description:
		Common utility functions for NUFFT.
	
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
 * NUFFT_Utils.h - Common utility functions for NUFFT
 *			  
 * Created 04/30/2009. 
 * Copyright 2009 by Apple, Inc. 
 */
 
#ifndef	_NU_FFT_UTILS_H_
#define _NU_FFT_UTILS_H_

#ifndef __cplusplus
#error Hey! This module is C++ only!
#endif

#include <libMatrixFFT/NUFFT.h>
#include <stdint.h>
#include <sys/types.h>

#include "PolyComplex.h"

extern "C" {

#pragma mark --- Common routines implemented in NUFFT.cpp ---

/* Special case for this module: x^y = 1 if y==0, for any x including 0 */
extern FFTFloat fExpI(
    FFTFloat base,
    unsigned iexp);

/* Calculate (f * i)^ iexp, f is FFTFloat, iexp is unsigned integer */
extern void fiExpI(
    FFTFloat base,
    unsigned iexp,
	FFTFloat *rtnReal,		// RETURNED
	FFTFloat *rtnImag);		// RETURNED

extern FFTFloat factorial(unsigned i);

/* optimized real mod D when D is pwr of 2 */
extern FFTFloat realMod2D(
    FFTFloat realVal,
    unsigned log2D);

/* optimized ssize_t mod D when D is a power of 2 */
extern size_t ssizeMod2D(
    ssize_t iVal,
    unsigned log2D);

/* 1-D complex in-place FFT on a PolyComplex. */
extern MFFTReturn fftPoly(
	NUFFTPlan       nuFftPlan,
    PolyComplex     &pc,
    unsigned        log2NumElts,
    size_t          numElts);

/* Dump PolyComplex to stdout */
extern void dumpNuFftPC(
	const char *title,
	PolyComplex &pc,
	size_t n);

}

#endif  /* _NU_FFT_UTILS_H_ */

