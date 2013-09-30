/*	File: fortranFFTWrapper.cpp 
	
	Description:
		Wrapper between Fortran and MatrixFFT.
	
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

#include <stdlib.h>
#include <strings.h>
#include <stdio.h>
#include <math.h>
#include <libMatrixFFT/MatrixFFT.h>
#include "fortranFFTWrapper.h"
#include <pthread.h>


/*
 * Note: this requires that the currently installed libMatrixFFT be configured
 * as double precision, interleaved complex. 
 */
#if     (FFT_DOUBLE_PREC == 0)
#error  fortranFFTWrapper requires double precision MatrixFFT.
#endif
#if     (FFT_SPLIT_COMPLEX == 0)
#error  fortranFFTWrapper requires Split complex format MatrixFFT.
#endif

/*
 * Global pthread_key_t and a lock to protect its creation. 
 */
static pthread_mutex_t planKeyLock = PTHREAD_MUTEX_INITIALIZER;
static pthread_key_t planKey = NULL;

/* 
 * Ensure this process has a valid planKey. Just create it once.
 */
static pthread_key_t getPlanKey()
{
    if(planKey != 0) {
        /* nonzero means we have a valid key */
        return planKey;
    }
    pthread_mutex_lock(&planKeyLock);
    if(planKey == 0) {
        /* avoid race condition with two threads creating this */
        pthread_key_create(&planKey, NULL);
    }
    pthread_mutex_unlock(&planKeyLock);
    return planKey;
}

/*
 * Each thread keeps one of these, kept on a per-thread basis via
 * pthread_setspecific().
 * Once an mfftPlan is created for a thread, that thread keeps on reusing
 * the same plan as long as the sizes, numDims, and isReal match. 
 * When they don't, the cached plan is freed and a new one replaces it.
 * The PerThreadPlanInfo itself is created (mallocd) at most once per
 * thread. 
 */
typedef struct {
    MatrixFFTPlan   mfftPlan;
    int             sizes[2];
    int             numDims;
    int             isReal;
} PerThreadPlanInfo;


/*
 * Obtain the per-thread MatrixFFTPlan for this thread, creating it
 * (along with its associated PerThreadPlanInfo) if necessary.
 * Returns NULL on any error.
 */
static MatrixFFTPlan perThreadPlan(
    int isReal,     
    int *sizes,      
    int numDims)    
{
    unsigned dims = (unsigned)numDims;
    if(dims > 2) {
        fprintf(stderr, "fft_wrapper_: Illegal numDims\n");
		
		// remove!!!
		printf("dims = %u\n", dims); 
        
		return NULL;
    }

    pthread_key_t theKey = getPlanKey();
    PerThreadPlanInfo *planInfo = (PerThreadPlanInfo *)pthread_getspecific(theKey);
    
    if(planInfo != NULL) {
        /* 
         * We have per-thread plan info for this thread.
         * See if it matches the caller's needs.
         */
        if((planInfo->isReal == isReal) && (planInfo->numDims == numDims)) {
            bool sizesMatch = true;
            for(int dim=0; dim<numDims; dim++) {
                if(planInfo->sizes[dim] != sizes[dim]) {
                    sizesMatch = false;
                    break;
                }
            }
            if(sizesMatch) {
                /* match found; use cached value */
                return planInfo->mfftPlan;
            }
        }
        
        /* 
         * We have a plan for this thread, but it doesn't match the caller's needs.
         * Free this cached plan and create a new one, keeping the PerThreadPlanInfo
         * intact.
         */
        mfftFreePlan(planInfo->mfftPlan);
    }
    else {
        /* 
         * First time through here for this thread.
         * Create new PerThreadPlanInfo from scratch.
         */
        planInfo = (PerThreadPlanInfo *)malloc(sizeof(PerThreadPlanInfo));
        if(planInfo == NULL) {
            fprintf(stderr, "fft_wrapper_: malloc failure\n");
            return NULL;
        }
        pthread_setspecific(theKey, planInfo);
    }
    
    /* 
     * One way or another, we have an empty PerThreadPlanInfo.
     * Cook up a new MatrixFFTPlan.
     */
    memset(planInfo, 0, sizeof(*planInfo));
    unsigned n[2];
    for(unsigned dim=0; dim<dims; dim++) {
        n[dim] = log2(sizes[dim]); 
    }
    bool isRealB = isReal ? true : false;
    
    MFFTReturn mrtn = mfftCreatePlan(dims, n, isRealB, 
        0,      // flags
        0,      // numThreads
        &planInfo->mfftPlan);
    if(mrtn) {
        mfftPrintErrInfo("mfftCreatePlan", mrtn);
        return NULL;
    }
    
    /*
     * OK, update planInfo to reflect the current plan.
     */
    for(unsigned dim=0; dim<dims; dim++) {
        planInfo->sizes[dim] = sizes[dim]; 
    }
    planInfo->numDims = numDims;
    planInfo->isReal = isReal;
    return planInfo->mfftPlan;
}

/* 
 * Public FFT wrapper function.
 */
void fft_wrapper_(     
    double real[], 
    double imag[], 
    int *forward,       /* 1 --> forward FFT; 0 --> inverse FFT */
    int *isReal,        /* 1 --> real-signal FFT; 0 --> complex FFT */
    int *sizes,         /* array of dimensions */
    int *numDims,       /* 1 or 2 dimensions */
    int *result)        /* RETURNED, actually an MFFTReturn */
{
	
	MFFTReturn mrtn;
	
	// Get per-thread plan.
    MatrixFFTPlan mfftPlan = perThreadPlan(*isReal, sizes, *numDims);
    if(mfftPlan == NULL) {
        *result = (int)MR_Internal;
        return;
    }
    
	// set up input, output
    FFTComplex fftSrc = {real, imag};
    FFTComplex fftDst = {real, imag};
    
    // Transpose output of forward FFT to   row order
    // Transpose input  of inverse FFT from row order
    // Normalize output of inverse FFT
	
    uint32_t optFlags;
    bool fwdFlag;
    if(*forward) {
        optFlags = MEF_TransposeOutput;
        fwdFlag = true;
    }
    else {
        optFlags = MEF_TransposeInput | MEF_NormOutput;
        fwdFlag = false;
    }
    
	// compute fft, in-place
	mrtn = mfftExecute(mfftPlan, optFlags, fwdFlag, &fftSrc, &fftDst);
	if(mrtn) {
        mfftPrintErrInfo("mfftCreatePlan", mrtn);
    }
    *result = (int)mrtn;
    return;
}

/* 
 * Release resources allocated by calls to fft_wrapper_() in the current thread.
 * Must be called from each thread which calls fft_wrapper_().
 * NOTE: we leak a pthread_key_t; we'd need one final release function
 * to be called once when all threads have called this function. 
 * I think that's overkill - it's just a long, per thread that uses it. 
 */
void fft_release_()
{
    pthread_key_t theKey = getPlanKey();
    PerThreadPlanInfo *planInfo = (PerThreadPlanInfo *)pthread_getspecific(theKey);
    
    if(planInfo == NULL) {
        /* No state cached for this thread; we're done */
        return;
    }
    
    if(planInfo->mfftPlan != NULL) {
        mfftFreePlan(planInfo->mfftPlan);
    }
    memset(planInfo, 0, sizeof(*planInfo));
    free(planInfo);
    pthread_setspecific(theKey, NULL);
}


