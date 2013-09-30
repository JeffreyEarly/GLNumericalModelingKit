/*	File: NUFFT_Priv.h 
	
	Description:
		Private #defines and types for NUFFT.
	
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
 * NUFFT_Priv.h - Private #defines and types for NUFFT.
 *			  
 * Created 04/30/2009. 
 * Copyright 2009 by Apple, Inc. 
 */
 
#ifndef	_NU_FFT_PRIV_H_
#define _NU_FFT_PRIV_H_

#include <libMatrixFFT/MatrixFFT.h>
#include <stdint.h>
#include <sys/types.h>
#include "PolyComplex.h"

#ifdef __cplusplus
extern "C" {
#endif

#pragma mark --- NUFFTPlanStruct definition ---

/* Our version of an NUFFTPlan */
struct NUFFTPlanStruct
{
    size_t                  D;              // size of signal */
    unsigned                log2D;          // log2(D) 
    unsigned                Bref;           // 2b / lg(b)
    unsigned                Bopt;           // 2.7 * b / lg(b)
    uint32_t                implFlags;      // passed to nufftCreatePlan
    
    /* 
     * theta, mu: NF_Reference, NF_Optimized only
     * theta is well-aligned; freeTheta to be freed
     */
    FFTFloat                *theta;         // size = D
    void                    *freeTheta;
    size_t                  *mu;            // ditto
    
    /* 
     * NF_Reference only: 
     * F is a three-dimension array.
     * l.s. index = j,    0..D, size of an FFT
     * med  index = beta, 0..B
     * m.s. index = K,    0..NUFFT_STEPS
     */
    FFTComplex              *F;
    FFTComplex              *freeF;        // to free

    /* 
     * cached 1/beta! for beta 0..B.
     * NF_Reference only. 
     */
    double                  *betaFact;
    
    /* NF_Optimized only */
    FFTComplex              *s;
    FFTComplex              *freeS;         // to free
    
	#if		!FFT_SPLIT_COMPLEX
	/* 
	 * Buffer for conversion to/from vDSP format.
     * Vector-aligned, freed via vdspBufFree.
     * Only used for interleaved complex with small(ish) signals,
     * i.e. where we use vDSP FFTs.
	 */
	vDSPComplex				vdspBuf;
	vDSPComplex				vdspBufFree;
	#endif	/* FFT_SPLIT_COMPLEX */
	
	TP_ThreadPool			threadPool;

    /*
     * Either vdspSetup or mfftPlan is valid, depending on the 
     * size of the signal. Both are NULL if the plan is only 
     * used for NF_Discrete.
     */
	FFT_Setup				vdspSetup;
    MatrixFFTPlan           mfftPlan;
};

#pragma mark --- common typedefs and #defines ---

#define NF_IMPL_MASK    (NF_Discrete | NF_Reference | NF_Optimized)

/* NUFFT_STEPS must be a power of 2 to take advantage of realMod2D */
#define NUFFT_STEPS         8
#define NUFFT_LOG2_STEPS    3

/*
 * Threshhold for determining whether to perform the low-level FFTs using 
 * vDSP or MatrixFFT. Signals larger than or dqual to the thresshold use
 * MatrixFFT. The threshhold is expressed as log2(N), or the 'n' that's
 * passed to nufftCreatePlan().
 */
#if     FFT_DOUBLE_PREC
#define NUFFT_MFT_THRESH    14
#else   
/* single precision */
#define NUFFT_MFT_THRESH    16
#endif  /* FFT_DOUBLE_PREC */

#pragma mark --- debugging macros ---

/* dump lots of detailed info */
#define DUMP_NU_FFTS        0

#if		DUMP_NU_FFTS
#define dumpNuFft(title, vn, cnt)		fftDump1DDSPComplex(title, vn, cnt)
#define dumpNuFftFloat(title, f, cnt)   fftDump1DFloat(title, f, cnt)
#define dumpNuFftSize(title, f, cnt)    fftDump1DSize(title, f, cnt)
#define ndprintf(s...)					printf(s)
#else
#define dumpNuFft(title, vn, cnt)
#define dumpNuFftFloat(title, f, cnt)  
#define dumpNuFftSize(title, f, cnt)    
#define ndprintf(s...)
#endif  /* DUMP_NU_FFTS */

#pragma mark --- Reference routines implemented in NUFFT_ref.cpp ---

/*
 * Unoptimized NUFFT.
 */
extern MFFTReturn nuFftExecuteRef(
	NUFFTPlan       nuFftPlan,
	FFTComplex		*inBuf,	
    FFTFloat        *tau,  
    FFTComplex      *outBuf);

/* 
 * Brute force, literal nonuniform DFT. 
 */
extern MFFTReturn nuFftExecuteDft(
	NUFFTPlan       nuFftPlan,
	FFTComplex		*inBuf,	
    FFTFloat        *tau,   
    FFTComplex      *outBuf);

#ifdef __cplusplus
}
#endif

#endif  /* _NU_FFT_PRIV_H_ */
