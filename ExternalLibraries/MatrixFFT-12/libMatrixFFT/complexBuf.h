/*	File: complexBuf.h 
	
	Description:
		Native representation of a complex data buffer.
	
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
 * complexBuf.h - Native representation of a complex data buffer.
 *			  
 * Created 12/16/2008. 
 * Copyright 2008 by Apple, Inc. 
 */
 
#ifndef	_COMPLEX_BUF_H_
#define _COMPLEX_BUF_H_

#include <libMatrixFFT/MatrixFFTConfig.h>
#include <libMatrixFFT/fftPrecision.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef	FFT_SPLIT_COMPLEX
#error	You must #define FFT_SPLIT_COMPLEX.
#endif

#if		FFT_SPLIT_COMPLEX

/*
 * This is analagous to vDSP's DSPDoubleSplitComplex or DSPSplitComplex.
 */
typedef struct {
	FFTFloat	*real;
	FFTFloat	*imag;
} FFTComplex;

#else	/* interleaved */

typedef struct {
	FFTFloat	real;
	FFTFloat	imag;
} FFTComplex;

#endif	/* FFT_SPLIT_COMPLEX */


#if		FFT_SPLIT_COMPLEX

/* 
 * Copy pointers - NOT the actual data - between vDSPComplex and FFTComplex 
 */

#define vDSPComplexToFFT(src, dst)	{ (dst)->real = (src)->realp; (dst)->imag = (src)->imagp; }
#define fftComplexToVDSP(src, dst)	{ (dst)->realp = (src)->real; (dst)->imagp = (src)->imag; }

#endif	/* FFT_SPLIT_COMPLEX */

#ifdef __cplusplus
}
#endif

#endif	/* _COMPLEX_BUF_H_ */
