/*	File: PolyComplex.h 
	
	Description:
		Polymorphic FFTComplex class.
	
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
 * Created 01/06/2009. 
 * Copyright 2009 by Apple, Inc. 
 */
 
#ifndef	_POLY_COMPLEX_H_
#define _POLY_COMPLEX_H_

#include <libMatrixFFT/complexBuf.h>
#include "fftIntel.h"

#ifndef FFT_SPLIT_COMPLEX
#error	You must #define FFT_SPLIT_COMPLEX.
#endif

/* 
 * The PolyComplex class can be used to write (some) code which manipulates complex
 * data in a manner which is independent of split/interleaved format. 
 * Although the underlying mechanisms differ for the two formats, the syntax of
 * PolyComplex is the same for both formats. 
 * PolyComplex can only be used when the real and imaginary streams operate
 * in parallel, side-by-side (i.e. you can't increment a real pointer without
 * incrementing the imaginary pointer in a split complex configuration). 
 * Offsets are conceptually in complex elements; in split complex format, adding
 * an offset results in adding the offset to both the real and imaginary
 * pointers; in interleaved format, an offset is in terms of an FFTComplex pointer. 
 */
 
#if		FFT_SPLIT_COMPLEX

class PolyComplex
{
public:
	PolyComplex(FFTComplex *fc)				: mFc(*fc) {}
	PolyComplex(FFTComplex *fc, size_t off) : mFc(*fc) { offset(off); }
	~PolyComplex() {}

    /* access underlying FFTComplex * */
    FFTComplex *fftCompl()                  { return &mFc; }

	/* set current pointers */
	void set(FFTComplex *fc)				{ mFc = *fc; }
	void set(FFTComplex *fc, size_t off)	{ mFc = *fc; offset(off); }
	
	/* offset current pointers by specified amount */
	void offset(size_t off)					{ mFc.real += off; mFc.imag += off; }
	void offset(int off)					{ mFc.real += off; mFc.imag += off; }

	/* simple offsets by +/- 1 - declared as prefix only to avoid copying object */
	PolyComplex &operator ++ ()				{ mFc.real++; mFc.imag++; return *this; }
	PolyComplex &operator -- ()				{ mFc.real--; mFc.imag--; return *this; }
	
	/* set current real/imaginary values */
	void real(FFTFloat r)					{ *mFc.real = r; }
	void imag(FFTFloat i)					{ *mFc.imag = i; }
	
	/* get current real/imaginary values */
	FFTFloat real()							{ return *mFc.real; }
	FFTFloat imag()							{ return *mFc.imag; }
	
	/* 
	 * Or if you prefer, access via ptrs...you can read or write, but don't 
	 * index like arrays.
	 */
	FFTFloat *realP()						{ return mFc.real; }
	FFTFloat *imagP()						{ return mFc.imag; }
	
	#if		FFT_INTEL
	
	/*** Intel SSE support ***/
	
	/* load two vectors, one real, one imaginary */
	void loadVect(FFTVector &r, FFTVector &i)
											{ r = FFTVectLoad(mFc.real);
											  i = FFTVectLoad(mFc.imag); }
	/* load two vectors, reverse order */
	void loadVectRev(FFTVector &r, FFTVector &i)
											{ r = FFTVectLoadRev(mFc.real);
											  i = FFTVectLoadRev(mFc.imag); }
	 /* store two vectors */
	 void storeVect(FFTVector r, FFTVector i)
											{ FFTVectStore(mFc.real, r);
											  FFTVectStore(mFc.imag, i); }
											  
	/* store two vectors, reverse order */
	void storeVectRev(FFTVector r, FFTVector i)
											{ FFTVectStoreRev(mFc.real, r);
											  FFTVectStoreRev(mFc.imag, i); }
	/* prefetch at both pointers */
	void prefetch()							{ fft_prefetch(mFc.real);
											  fft_prefetch(mFc.imag); }
	#endif	/* FFT_INTEL */
	
    /* fill contents with const data */
    void fill(FFTFloat val, size_t n)       { FFT_vfill(val, mFc.real, n);
                                              FFT_vfill(val, mFc.imag, n);
                                            }
                                            
    /* add a complex, expressed as a pair of FFTFloats */
    void add(FFTFloat r, FFTFloat i)        { *mFc.real += r; *mFc.imag += i; }
    
    /* multiply by a scalar value */
    void scalarMul(FFTFloat f)              { *mFc.real *= f; *mFc.imag *= f; }
    
private:
    FFTComplex mFc;
	
	/* not supported */
	PolyComplex(PolyComplex &src);
	void operator = (const PolyComplex &);
};

#else	/* !FFT_SPLIT_COMPLEX */

/* Interleaved complex format */

class PolyComplex
{
public:
	PolyComplex(FFTComplex *fc)				: mFc(fc) {}
	PolyComplex(FFTComplex *fc, size_t off) : mFc(fc) { offset(off); }
	~PolyComplex() {}
	
    /* access underlying FFTComplex * */
    FFTComplex *fftCompl()                  { return mFc; }
    
	/* set current pointer */
	void set(FFTComplex *fc)				{ mFc = fc; }
	void set(FFTComplex *fc, size_t off)	{ mFc = fc; offset(off); }
	
	/* offset current pointer by specified amount */
	void offset(size_t off)					{ mFc += off; }
	void offset(int off)					{ mFc += off; }
	
	/* simple offsets by +/- 1 - declared as prefix only to avoid copying object */
	PolyComplex &operator ++ ()				{ mFc++; return *this; }
	PolyComplex &operator -- ()				{ mFc--; return *this; }
	
	/* set current real/imaginary values */
	void real(FFTFloat r)					{ mFc->real = r; }
	void imag(FFTFloat i)					{ mFc->imag = i; }
	
	/* get current real/imaginary values */
	FFTFloat real()							{ return mFc->real; }
	FFTFloat imag()							{ return mFc->imag; }

	/* 
	 * Or if you prefer, access via ptrs...you can read or write, but don't 
	 * index like arrays. 
	 */
	FFTFloat *realP()						{ return &mFc->real; }
	FFTFloat *imagP()						{ return &mFc->imag; }
	
	#if		FFT_INTEL
	
	/*** Intel SSE support ***/
	
	/* load two vectors, one real, one imaginary */
	void loadVect(FFTVector &r, FFTVector &i)
		{ fftLoadComplex(mFc, i, r); }

	/* load two vectors, reverse order */
	void loadVectRev(FFTVector &r, FFTVector &i)
		{ fftLoadComplexRev(mFc, i, r); }

	 /* store two vectors */
	 void storeVect(FFTVector r, FFTVector i)
		{ fftStoreComplex(mFc, i, r); }
											  
	/* store two vectors, reverse order */
	void storeVectRev(FFTVector r, FFTVector i)
		{ fftStoreComplexRev(mFc, i, r); }

	/* prefetch */
	void prefetch()							{ fft_prefetch(mFc); }

	#endif	/* FFT_INTEL */
	
    /* fill contents with const data */
    void fill(FFTFloat val, size_t n)       { FFT_vfill(val, &mFc->real, n << 1); }
                                            
    /* add a complex, expressed as a pair of FFTFloats */
    void add(FFTFloat r, FFTFloat i)        { mFc->real += r; mFc->imag += i; }

    /* multiply by a scalar value */
    void scalarMul(FFTFloat f)              { mFc->real *= f; mFc->imag *= f; }

private:
	FFTComplex *mFc;
	
	/* not supported */
	PolyComplex(PolyComplex &src);
	void operator = (const PolyComplex &);
};


#endif	/* FFT_SPLIT_COMPLEX */

#endif	/* _POLY_COMPLEX_H_ */
