/*	File:  fftDebug.h 

	Description:
		Debugging support for MatrixFFT library.
	
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
 * Created 01/06/2009. 
 * Copyright 2009 by Apple, Inc. 
 */
 
#ifndef	_FFT_DEBUG_H_
#define _FFT_DEBUG_H_

#include <assert.h>
#include <stdio.h>
#include <libMatrixFFT/complexBufUtils.h>

#ifdef __cplusplus
extern "C" {
#endif

#pragma mark --- matrix and submatrix debugging display support ---

/* dump a submatrix to stdout */

#if		FFT_SPLIT_COMPLEX
extern void fftDumpSub(	
	const char			*title,
	const FFTFloat		*buf,
	size_t				rowSize);	/* in floats */
#else
extern void fftDumpSub(	
	const char			*title,
	const FFTComplex	*buf,
	size_t				rowSize);	/* in FFTComplex */
#endif

#define DUMP_MATRIX		0		/* dump matrix results of each op */
#define DUMP_SUB_SWAP	0		/* dump results of each submatrix transpose */
#define DUMP_TRANS		0		/* log transpose ops */

#if		DUMP_MATRIX
#define dumpMatrix(t, b, n)			fftDump1DComplex(t, b, n)
#define dumpMatrixRect(t, b, r, c)	fftDumpMatrixRect(t, b, r, c)
#define dumpRectReal(t, b, r, c)	fftDumpMatrixRect(t, b, r, c>>1)
#define dumpBuf(t, b, r, c)			fftDumpBuf(t, b, r, c)
#define dprintf(args...)			printf(args)
#else
#define dumpMatrix(t, b, n)			
#define dumpMatrixRect(t, b, r, c)	
#define dumpRectReal(t, b, r, c)
#define dumpBuf(t, b, r, c)		
#define dprintf(args...)
#endif

#if		DUMP_SUB_SWAP
#define dumpSub(t, b, r)			fftDumpSub(t, b, r)
#define dumpSubSwapPrint(x...)		printf(x)
#else
#define dumpSub(t, b, r)
#define dumpSubSwapPrint(x...)
#endif

#if		DUMP_TRANS
#define dtprintf(args...)			printf(args)
#else
#define dtprintf(args...)
#endif

#pragma mark --- misc. debugging macros ---

#ifdef	DEBUG
#define RFASSERT(s)			assert(s)
#define FFT_INLINE
#define dbprintf(args...)	printf(args);
#else
#define RFASSERT(s)
#define FFT_INLINE		inline
#define dbprintf(args...)
#endif

/* when true, log library-mallocd memory to stdout */
#define FFT_LOG_MALLOC		0
#if		FFT_LOG_MALLOC
#define mdprint(args...)	printf(args)
#else
#define mdprint(args...)
#endif	/* FFT_LOG_MALLOC */


#ifdef __cplusplus
}
#endif

#endif	/* _FFT_DEBUG_H_ */
