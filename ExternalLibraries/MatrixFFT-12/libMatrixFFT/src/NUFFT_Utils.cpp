/*	File: NUFFT_Utils.cpp
	
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
 * NUFFT_Utils.cpp - Common utility functions for NUFFT
 *			  
 * Created 04/30/2009. 
 * Copyright 2009 by Apple, Inc. 
 */

#include "NUFFT_Utils.h"
#include "NUFFT_Priv.h"
#include <libMatrixFFT/vdspUtils.h>
#include "fftDebug.h"

/*
 * Various flavors of f^i.
 */
 
/* Calculate f^iexp, f is FFTFloat, iexp is unsigned integer */
static inline FFTFloat powerI(
    FFTFloat base,
    unsigned iexp)
{
    /* 
     * This module can skip this one, since it's handled by our (sole)
     * caller...normally you need this! */
    /*
    if(base == 0) {
        return 0.0;
    }
    */
	FFTFloat rtn = 0.0;
    if(iexp <= 4) {
        /* Unroll the switch here for max speed */
        switch(iexp) {
            /* not needed for this implementation, but normally you need it...
            case 0:
                rtn = 1.0;
                break;
             */
            case 1:
                rtn = base;
                break;
            /* remainders not strictly needed */
            case 2:
                rtn = base * base;
                break;
            case 3:
                rtn = base * base * base;
                break;
            case 4:
            default:
            {
                FFTFloat sq = base * base;
                rtn = sq * sq;
                break;
            }
        }
    }
    else {
        /* iexp > 4, recurse */
        if((iexp & 0x01) == 0) {
            /* even exponent */
            rtn = powerI(base*base, iexp>>1);
        }
        else {
            /* odd exponent */
            rtn = base * powerI(base*base, iexp>>1);
        }
    }
	return rtn;
}

/* Special case for this module: x^y = 1 if y==0, for any x including 0 */
FFTFloat fExpI(
    FFTFloat base,
    unsigned iexp)
{
	FFTFloat rtn = 0.0;
	if(iexp == 0) {
		rtn = 1.0;
	}
    else if(base == 0.0) {
        rtn = 0.0;
    }
	else {
		rtn = powerI(base, iexp);
	}

    return rtn;
}

/* Calculate (f * i)^ iexp, f is FFTFloat, iexp is unsigned integer */
void fiExpI(
    FFTFloat base,
    unsigned iexp,
	FFTFloat *rtnReal,		// RETURNED
	FFTFloat *rtnImag)		// RETURNED
{
	/*
 	 * First calculate the value without the i term.
	 */
	FFTFloat raw = fExpI(base, iexp);

	/* 
  	 * Now handle the imaginary portion.
     * Saving code space would require writing a zero to 
     * both output fields here, and updating one of them in 
     * the switch below...but for runtime optimization we
     * avoid that and write both fields exactly once. 
	 */	
	switch(iexp & 0x03) {   /* i.e., iexp mod 4 */
		case 0:
			/* i^4 = (-1)^2 = 1 */
			*rtnReal = raw;
            *rtnImag = 0.0;
			break;
		case 1:
			/* i */
			*rtnImag = raw;
            *rtnReal = 0.0;
            break;
		case 2:
			/* i^2 = -1 */
			*rtnReal = -raw;
            *rtnImag = 0.0;
			break;
		case 3:
		default:
			/* i^3 = -i */
			*rtnImag = -raw;
            *rtnReal = 0.0;
			break;
	}
}

FFTFloat factorial(unsigned i)
{
	FFTFloat rtn = 1.0;
	for(unsigned dex=2; dex<=i; dex++) {
		rtn *= (FFTFloat)dex;
	}
	return rtn;
}

/* optimized real mod D when D is pwr of 2 */
FFTFloat realMod2D(
    FFTFloat realVal,
    unsigned log2D)
{
    ssize_t iR = (ssize_t)FFTFloor(realVal); 
    iR >>= log2D;
    FFTFloat f = (FFTFloat)(iR << log2D);
    return realVal - f;
}

/* optimized ssize_t mod D when D is a power of 2 */
size_t ssizeMod2D(
    ssize_t iVal,
    unsigned log2D)
{
    ssize_t i = iVal >> log2D;
    ssize_t tmp = i << log2D;
    ssize_t rtn = iVal - tmp;
    RFASSERT(rtn >= 0);
    return (size_t)rtn;
}

/*
 * 1-D complex in-place FFT on a PolyComplex, three forms.
 */
MFFTReturn fftPoly(
	NUFFTPlan       nuFftPlan,
    PolyComplex     &pc,
    unsigned        log2NumElts,
    size_t          numElts)
{
    MFFTReturn mrtn = MR_Success;
    
    if(nuFftPlan->mfftPlan != NULL) {
        /* MatrixFFT */
        FFTComplex *fc = pc.fftCompl();
        mrtn = mfftExecute(nuFftPlan->mfftPlan, MEF_TransposeOutput, true, fc, fc);
        if(mrtn) {
            printf("***fftPoly: mfftExecute() error\n");
            mfftPrintErrInfo("mfftExecute", mrtn);
        }
    }
    else {
        /* vDSP - two forms */
        #if     FFT_SPLIT_COMPLEX
            /* vDSP native format */
            RFASSERT(nuFftPlan->vdspSetup != NULL);
            vDSPComplex vBuf = {pc.realP(), pc.imagP()};
            dumpNuFft("FFT In", &vBuf, numElts);
            FFTComplex1d(nuFftPlan->vdspSetup, &vBuf, log2NumElts, FFT_FORWARD);
            dumpNuFft("FFT Out", &vBuf, numElts);
        #else
            /* interleaved format - requires copy to/from temp buffer */
            RFASSERT(nuFftPlan->vdspSetup != NULL);
            RFASSERT(nuFftPlan->vdspBuf.realp != NULL);
            fftIntToVDSP(pc.fftCompl(), &nuFftPlan->vdspBuf, numElts);
            dumpNuFft("FFT In", &nuFftPlan->vdspBuf, numElts);
            FFTComplex1d(nuFftPlan->vdspSetup, &nuFftPlan->vdspBuf, log2NumElts, FFT_FORWARD);
            dumpNuFft("FFT Out", &nuFftPlan->vdspBuf, numElts);
            fftVDSPToInt(&nuFftPlan->vdspBuf, pc.fftCompl(), numElts);
        #endif  /* FFT_SPLIT_COMPLEX */
    }
    return mrtn;
}   

#if		FFT_SPLIT_COMPLEX

void dumpNuFftPC(
	const char *title,
	PolyComplex &pc,
	size_t n)
{
	FFTComplex fc = {pc.realP(), pc.imagP()};
	fftDump1DComplex(title, &fc, n);
}

#else	/* FFT_SPLIT_COMPLEX */

void dumpNuFftPC(
	const char *title,
	PolyComplex &pc,
	size_t n)
{
	fftDump1DComplex(title, pc.fftCompl(), n);
}

#endif	/* FFT_SPLIT_COMPLEX */


