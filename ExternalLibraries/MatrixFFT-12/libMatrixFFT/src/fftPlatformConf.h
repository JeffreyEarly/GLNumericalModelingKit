/*	File:  fftPlatformConf.h 

	Description:
		Platform-dependent configuration support for MatrixFFT library.
	
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
 * fftPlatformConf.h - Platform-dependent configuration support for MatrixFFT library.
 *			  
 * Created 11/05/2009. 
 * Copyright 2009 by Apple, Inc. 
 */
 
#ifndef	_FFT_PLATFORM_CONF_H_
#define _FFT_PLATFORM_CONF_H_

#include <stdint.h>
#include <sys/types.h>      /* for size_t */
#include <libMatrixFFT/src/fftPriv.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Values for FFT_PLATFORM.
 * This is used is various places to select platform-dependent
 * configuration parameters. 
 * You can define this here, or via a compile time flag via -D.
 */
#define FFT_Xeon_8_32GB     0       /* 8-core Xeon, 32 GB RAM */
#define FFT_Xeon_4_4GB      1       /* 4-core Xeon, 4 GB RAM */
#define FFT_Nehalem_8_32GB  2       /* 8-core Nehalem, 32 GB RAM */
#define FFT_Other           3

#ifndef FFT_PLATFORM
#define FFT_PLATFORM        FFT_Nehalem_8_32GB
#endif  
    
/* 
 * Submatrix size in columns.
 * This is pretty confusing: but remember that column peel and unpeel ops
 * operate on submatrix-aligned data. In split complex config, peel and unpeel
 * ops operate on real and imaginary data completely separately, so a submatrix
 * size is in terms of floats. In interleave format, peel and unpeel operate
 * on combined real and imagtinary so a submatrix size is in terms of complex
 * elements. 
 */
#if		FFT_SPLIT_COMPLEX
#define FFT_SUBMATRIX_ATOM		FFT_FLOATS_PER_SUBMATRIX
#else
#define FFT_SUBMATRIX_ATOM		FFT_COMPLEX_PER_SUBMATRIX
#endif

#pragma mark --- Configuration accessors used at runtime ---

/***
 *** These are used internally to access platform-dependent configuration info.
 ***/
 
/*
 * Infer rectangle dimensions from 1-D size.
 */
extern void mfftMakeRectangle(
	unsigned    nIn,		// log2(size)
	uint32_t    optFlags,   // only MCF_HintTranspose used
	unsigned    *nRow,		// log2(numRows)
	unsigned    *nCol);     // log2(numCols)

/* 
 * Obtain the size of a column stripe, in complex elements, 
 * given FFT type, log2(numRows) and l2CachePerEngine.
 */
extern size_t mfftColumnStripeSize(
    unsigned    nRow,               // log2(numRows)
    bool        isReal,
    unsigned    numDims,
    uint32_t    l2CachePerEngine);
    
/*
 * Determine if it's best to fall back and use raw vDSP calls for
 * FFT of specified type and size. For 2-D ops, log2n is the log
 * of the total signal size (i.e. log2n = log2(numRows * numCols).
 * Returns true if vDSP calls should be used.
 */
extern bool mfftShouldUseVdsp(
    unsigned    log2n,               // log2(total signal size)
    bool        isReal,
    unsigned    numDims);

#pragma mark --- Configuration tweakers ---

/***
 *** Routines to modify the normal config parameters at runtime.
 *** These are used by the tools which create our config tables.
 ***/
 
/*
 * Specify signed offset to the normal derivation of 2-D coordinates from
 * a 1-D size. Note that this always only applies to the complex 1-D
 * rectangle. A 1-D real FFT has a 1-D complex subplan do this for the 
 * complex representation, so we only need one fo these.
 *
 * The specified offset applies to log2(numRows) when numRows is derived
 * from n[0] passed to mfftCreatePlan().
 */
extern void mfftSetRectangleOffset(
	int     rectOffset,
    bool    enableOffset);
	
/* 
 * Specify column stripe size in complex elements.
 * Specify 0 to revert to using built-in tables.
 */
extern void mfftSetColumnStripeSize(
    bool    isReal,
    size_t  colStripeSize);

/* 
 * Force vDSP or MatrixFFT in subsequent calls to mfftShouldUseVdsp.
 */
typedef enum {
    MF_Default = 0,         // use per config/signal
    MF_ForceVdsp,           // force raw vdSP
    MF_ForceMfft            // force MatrixFFT (raw vDSP disabled)
} MFFT_ForceVdsp;

extern void mfftSetForceVdsp(
    MFFT_ForceVdsp  forceVdsp);


#ifdef __cplusplus
}
#endif

#endif	/* _FFT_PLATFORM_CONF_H_ */
