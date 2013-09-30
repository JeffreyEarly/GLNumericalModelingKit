/*	File: fftPrecision.h
	
	Description:
		Precision-specific #defines.
	
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
 * Created 12/16/2008. 
 * Copyright 2008 by Apple, Inc. 
 */
 
#ifndef	_FFT_PRECISION_H_
#define _FFT_PRECISION_H_

#include <libMatrixFFT/MatrixFFTConfig.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef	FFT_DOUBLE_PREC
#error	You must #define FFT_DOUBLE_PREC.
#endif

#if		FFT_DOUBLE_PREC

/* double precision floating point */

/* basic types */
typedef double								FFTFloat;

/* redirected math.h calls */
#define FFTRound(f)							lround(f)
#define FFTSin(x)							sin(x)
#define FFTCos(x)							cos(x)
#define FFTAbs(x)							fabs(x)
#define FFTExp(x)							exp(x)
#define FFTCeil(x)                          ceil(x)
#define FFTLog2(x)                          log2(x)
#define FFTRoundl(x)                        lrint(x)
#define FFTFloor(x)                         floor(x)

#else	/* !FFT_DOUBLE_PREC */

/* single precision floating point */

/* basic types */
typedef float								FFTFloat;

/* redirected function calls */
#define FFTRound(f)							lroundf(f)
#define FFTSin(x)							sinf(x)
#define FFTCos(x)							cosf(x)
#define FFTAbs(x)							fabsf(x)
#define FFTExp(x)							expf(x)
#define FFTCeil(x)                          ceilf(x)
#define FFTLog2(x)                          log2f(x)
#define FFTRoundl(x)                        lrintf(x)
#define FFTFloor(x)                         floorf(x)

#endif	/* FFT_DOUBLE_PREC */

#ifdef __cplusplus
}
#endif

#endif	/* _FFT_PRECISION_H_ */

