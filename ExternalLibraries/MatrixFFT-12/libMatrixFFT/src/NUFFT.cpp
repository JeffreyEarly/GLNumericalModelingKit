/*	File: NUFFT.cpp
	
	Description:
		Implementation of Nonuniform FFT. 
	
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
 * NUFFT.cpp - Implementation of Nonuniform FFT.
 *			  
 * Created 04/09/2009. 
 * Copyright 2009 by Apple, Inc. 
 */
 
#include <libMatrixFFT/NUFFT.h>
#include <libMatrixFFT/fftUtils.h>
#include <libMatrixFFT/complexBufUtils.h>
#include "NUFFT_Priv.h"
#include "NUFFT_Utils.h"
#include "fftDebug.h"
#include "PolyComplex.h"
#include <stdlib.h>
#include <CoreFoundation/CoreFoundation.h>

/* dump lots of detailed info - VERY verbose */
#define DUMP_NU_VERB        0

#if		DUMP_NU_VERB
#define nvprintf(s...)					printf(s)
#else
#define nvprintf(s...)
#endif  /* DUMP_NU_VERB */

#pragma mark --- Initial input signal modification ---

/*
 * Initialize input signal:
 *
 *   x[j] = x[j] * e ^ (-i * pi * theta[j])
 */
#if		!FFT_INTEL

static void nuFftInit(
    FFTComplex *x,
    FFTFloat *theta,
    size_t N)
{
    PolyComplex inPoly(x);
    FFTFloat negPi = -1.0 * M_PI;
    
    for(size_t j=0; j<N; j++) {
        /*
         * mul   = e ^ (-i * pi * theta[j])
         * expon = -pi * theta[j]
         * mul   = cos(expon) + (i * sin(expon))
         */
        FFTFloat expon = negPi * *theta++;
        FFTFloat mulR = FFTCos(expon);
        FFTFloat mulI = FFTSin(expon);

        /* x[j] *= mul */
        FFTFloat iR = inPoly.real();
        FFTFloat iI = inPoly.imag();
        FFTFloat xR = (iR * mulR) - (iI * mulI);
        FFTFloat xI = (iR * mulI) + (iI * mulR);
        inPoly.real(xR);
        inPoly.imag(xI);
        
        ++inPoly;
    }
}

#else   /* FFT_INTEL */

/*
 * Intel vectorized version of nuFftInit.
 */
static void nuFftInit(
    FFTComplex *x,
    FFTFloat *theta,
    size_t N)
{
    PolyComplex inPoly(x);
    FFTVector negPi = FFTVectSet1(-1.0 * M_PI);
    FFTVector *thetaV = (FFTVector *)theta;
    
    for(size_t j=0; j<N; j+=FFT_FLOATS_PER_VECTOR) {
        /*
         * mul   = e ^ (-i * pi * theta[j])
         * expon = -pi * theta[j]
         * mul   = cos(expon) + (i * sin(expon))
         */
        FFTVector expon = FFTVectMul(negPi, *thetaV);
        FFTVector mulR = fftVectCosine(expon);
        FFTVector mulI = fftVectSine(expon);

        /* 
         * x[j] *= mul 
         * first load x[j] info {iR,iJ}
         */
        FFTVector iR;// = inPoly.real();
        FFTVector iI;// = inPoly.imag();
        inPoly.loadVect(iR, iI);
        
        /* 
         * {xR,xI} := {iR,iI} * mul
         * xR = (iR * mulR) - (iI * mulI);
         * xI = (iR * mulI) + (iI * mulR);
         */
        FFTVector xR = FFTVectSub(FFTVectMul(iR, mulR), FFTVectMul(iI, mulI));
        FFTVector xI = FFTVectAdd(FFTVectMul(iR, mulI), FFTVectMul(iI, mulR));
        
        inPoly.storeVect(xR, xI);
        
        inPoly.offset(FFT_FLOATS_PER_VECTOR);
        thetaV++;
    }
}
 
#endif  /* FFT_INTEL */

#pragma mark --- Inner loop exponentiation ---

/* 
 * Perform a portion of crucial inner loop:
 *
 * outBuf[k] += ((-2 * pi * i * (k - N/2) / N) ^ m) * s[k] / m!
 *
 * Called out as a ThreadPool callout, or directly from nuFftLoop.
 */
#define NUFFT_DISABLE_INTEL     0

#if		(NUFFT_DISABLE_INTEL || !FFT_INTEL)

/* generic version */
static MFFTReturn nuFftLoopCallout(
	TP_TaskUnion *u)
{
    TP_NUFFT_Loop *tnl = &u->nuFftLoop;
    
    NUFFTPlan nuFftPlan = tnl->nuFftPlan;
    PolyComplex outPoly(tnl->outBuf, tnl->startK);
    PolyComplex s(nuFftPlan->s, tnl->startK);
    ssize_t endK = tnl->startK + tnl->kToProcess;
    ssize_t Nover2 = nuFftPlan->D >> 1;
    
    for(ssize_t k=tnl->startK; k<endK; k++) {
    
        /* expTerm := (-2 * pi * i * (k - N/2) / N) ^ m */
        FFTFloat base = tnl->negTwoPiOverN * (FFTFloat)(k - Nover2);
        FFTFloat expTermR;
        FFTFloat expTermI;
        fiExpI(base, tnl->m, &expTermR, &expTermI);
        
        /* addTerm := expTerm * s[k] */
        FFTFloat sI = s.imag();
        FFTFloat sR = s.real();
        FFTFloat addTermR = (sR * expTermR) - (sI * expTermI);
        FFTFloat addTermI = (sR * expTermI) + (sI * expTermR);
        
        nvprintf("--k=%lu base %.2f expTermR %.2f expTermI %.2f addTermR %.2f addTermI %.2f\n", 
            (unsigned long)k, base,
            expTermR, expTermI,
            addTermR, addTermI);
        
        /* addTerm := expTerm * s[k] / m! */
        addTermR *= tnl->oneOverMFact;
        addTermI *= tnl->oneOverMFact;
        
        nvprintf("--    addTermR %.2f addTermI %.2f\n", addTermR, addTermI);

        /* outBuf[k] += addTerm */
        outPoly.add(addTermR, addTermI);

        ++outPoly;
        ++s;
    }   /* for k */
    
    return MR_Success;
}

#else   /* FFT_INTEL */

/* 
 * Intel-optimized version of nuFftLoopCallout.
 * Currently only a win for single precision; a wash for double.
 */

static MFFTReturn nuFftLoopCallout(
	TP_TaskUnion *u)
{
    TP_NUFFT_Loop *tnl = &u->nuFftLoop;
    
    NUFFTPlan nuFftPlan = tnl->nuFftPlan;
    PolyComplex outPoly(tnl->outBuf, tnl->startK);
    PolyComplex s(nuFftPlan->s, tnl->startK);
    ssize_t endK = tnl->startK + tnl->kToProcess;
    ssize_t Nover2 = nuFftPlan->D >> 1;
    
    FFTVector neg2PiOverN  = FFTVectSet1(tnl->negTwoPiOverN);
    FFTVector oneOverMFact = FFTVectSet1(tnl->oneOverMFact);
    FFTVector fpv          = FFTVectSet1((FFTFloat)FFT_FLOATS_PER_VECTOR);
    FFTVectUnion kMinusNOver2;
    
    /*
     * Initial values for (k - N/2) vector; within the loop
     * this updates via += FFT_FLOATS_PER_VECTOR
     */
    for(ssize_t dex=0; dex<FFT_FLOATS_PER_VECTOR; dex++) {
        kMinusNOver2.f[dex] = (FFTFloat)(dex + (ssize_t)tnl->startK - Nover2);
    }
    
    for(ssize_t k=tnl->startK; k<endK; k+=FFT_FLOATS_PER_VECTOR) {
    
        /* 
         * base    :=  -2 * pi * (k - N/2) / N
         * expTerm := (-2 * pi * i * (k - N/2) / N) ^ m 
         *          =  base ^ m
         */
        FFTVectUnion base;
        base.v = FFTVectMul(neg2PiOverN, kMinusNOver2.v);
        FFTVectUnion expTermR;
        FFTVectUnion expTermI;
        for(unsigned dex=0; dex<FFT_FLOATS_PER_VECTOR; dex++) {
            fiExpI(base.f[dex], tnl->m, &expTermR.f[dex], &expTermI.f[dex]);
        }
        
        /* 
         * addTerm := expTerm * s[k] 
         *     real = (sR * expTermR) - (sI * expTermI
         *     imag = (sR * expTermI) + (sI * expTermR
         */
        FFTVector vR, vI;
        s.loadVect(vR, vI);
        FFTVector addTermR, addTermI;
        addTermR = FFTVectSub(FFTVectMul(vR, expTermR.v), FFTVectMul(vI, expTermI.v));
        addTermI = FFTVectAdd(FFTVectMul(vR, expTermI.v), FFTVectMul(vI, expTermR.v));
        
        /* addTerm := expTerm * s[k] / m! */
        addTermR = FFTVectMul(addTermR, oneOverMFact);
        addTermI = FFTVectMul(addTermI, oneOverMFact);
        
        /* outBuf[k] += addTerm */
        outPoly.loadVect(vR, vI);
        vR = FFTVectAdd(vR, addTermR);
        vI = FFTVectAdd(vI, addTermI);
        outPoly.storeVect(vR, vI);

        outPoly.offset(FFT_FLOATS_PER_VECTOR);
        s.offset(FFT_FLOATS_PER_VECTOR);
        kMinusNOver2.v = FFTVectAdd(kMinusNOver2.v, fpv);
    }   /* for k */
    
    return MR_Success;
}

#endif  /* FFT_INTEL */

/* min log2D to use threading */
#define NU_LOOP_THREAD_MIN      12

/* 
 * outBuf[k] += ((-2 * pi * i * (k - N/2) / N) ^ m) * s[k] / m!
 */
static void nuFftLoop(
    NUFFTPlan nuFftPlan,
    unsigned m,
    FFTFloat oneOverMFact,        /* optimization, caller calculates this */
    FFTFloat negTwoPiOverN,       /* ditto */
    FFTComplex *outBuf)
{
    unsigned numThreads = nuFftPlan->threadPool.numThreads;
    size_t N = nuFftPlan->D;

    if((nuFftPlan->log2D < NU_LOOP_THREAD_MIN) ||
       (numThreads <= 1)) {
        /* no threading */
        TP_TaskUnion tu;
        TP_NUFFT_Loop *tnl = &tu.nuFftLoop;
        
        tnl->nuFftPlan      = nuFftPlan;
        tnl->m              = m;
        tnl->oneOverMFact   = oneOverMFact;
        tnl->negTwoPiOverN  = negTwoPiOverN;
        tnl->outBuf         = outBuf;
        tnl->startK         = 0;
        tnl->kToProcess     = N;
        
        nuFftLoopCallout(&tu);
        return;
    }
    
    /* 
     * Divvy up the job amongst all configured threads 
     */
    size_t kPerThread = (N + numThreads - 1) / numThreads;
    size_t kRem = N;
    size_t currK = 0;
    
    for(unsigned dex=0; dex<numThreads; dex++) {
		TP_PerThread *pt   = &nuFftPlan->threadPool.perThread[dex];
		TP_Task *task      = &pt->task;
		TP_NUFFT_Loop *tnl = &task->u->nuFftLoop;

		task->op		   = TPO_App;
		task->threadFcn    = nuFftLoopCallout;
        tnl->nuFftPlan     = nuFftPlan;
        tnl->m             = m;
        tnl->oneOverMFact  = oneOverMFact;
        tnl->negTwoPiOverN = negTwoPiOverN;
        tnl->outBuf        = outBuf;
        tnl->startK        = currK;
        if(dex == (numThreads - 1)) {
            /* Might be more for this one, use all remaining */
            tnl->kToProcess = kRem;
        }
        else {
            tnl->kToProcess = kPerThread;
        }
        kRem  -= tnl->kToProcess;
        currK += tnl->kToProcess;
    }
    
	/* GO */
	MFFTReturn ourRtn = tpThreadDispatch(&nuFftPlan->threadPool, numThreads);
	if(ourRtn) {
		return;
	}
	tpThreadFinish(&nuFftPlan->threadPool, numThreads);
}

#pragma mark --- Public API ---

/* 
 * Create an NUFFTPlan for nonuniform FFT. 
 *
 * Returned NUFFTPlan must be freed via nufftFreePlan(). 
 */ 
MFFTReturn nufftCreatePlan(
	unsigned			dims,		
	unsigned			*n,			
    unsigned            bits,    
	uint32_t			optFlags,
	unsigned			numThreads,
	NUFFTPlan           *nuFftPlan)		/* RETURNED */
{
    if(dims != 1) {
        printf("***nufftCreatePlan: only one dimension supported\n");
        return MR_IllegalArg;
    }
    if((n == NULL) || (n[0] < 3) || (bits == 0)) {
        return MR_IllegalArg;
    }
    if((optFlags & NF_IMPL_MASK) == 0) {
        printf("***nufftCreatePlan: No implementation flag specified\n");
        return MR_IllegalArg;
    }
    
    bool implRef = optFlags & NF_Reference;
    bool implOpt = optFlags & NF_Optimized;
    
    NUFFTPlan rtnPlan = (NUFFTPlan)malloc(sizeof(*rtnPlan));
    if(rtnPlan == NULL) {
		return MR_Memory;
	}
	/* subsequent errors to errOut: */
    
    MFFTReturn ourRtn = MR_Success;
    memset(rtnPlan, 0, sizeof(*rtnPlan));

    rtnPlan->D = 1 << n[0];
    rtnPlan->log2D = n[0];
    rtnPlan->implFlags = optFlags;
    
    /* 
     * B        = ceil(2b/lg(b)) for reference, ceil(2.7b/lg(b)) for optimized
     * mu[j]    = round(tau[j])
     * theta[j] = tau[j] - mu[j]
     */
    size_t D = rtnPlan->D;
    if(implRef) {
        rtnPlan->Bref = (unsigned)FFTCeil(2.0 * bits / FFTLog2(bits));
    }
    if(implOpt) {
        rtnPlan->Bopt = (unsigned)FFTCeil(2.7 * (FFTFloat)bits / FFTLog2(bits));
    }
    
    if(implRef || implOpt) {
        rtnPlan->mu = (size_t *)malloc(D * sizeof(size_t));
        if(rtnPlan->mu == NULL) {
            ourRtn = MR_Memory;
            goto errOut;
        }
        
        rtnPlan->theta = (FFTFloat *)fftAllocAlign(D * sizeof(FFTFloat),
            &rtnPlan->freeTheta);
        if(rtnPlan->theta == NULL) {
            ourRtn = MR_Memory;
            goto errOut;
        }
    }
    if(implRef) {
        size_t fSize;
        
        rtnPlan->betaFact = (double *)malloc(rtnPlan->Bref * sizeof(double));
        if(rtnPlan->betaFact == NULL) {
            ourRtn = MR_Memory;
            goto errOut;
        }
        
        for(unsigned beta=0; beta<rtnPlan->Bref; beta++) {
            rtnPlan->betaFact[beta] = 1.0 / factorial(beta);
        }
    
        /* 
         * Three-D F declared as 1-D array
         *
         * m.s. index   0..NUFFT_STEPS
         * middle index 0..B
         * l.s. index   0..D
         */
        fSize = D * rtnPlan->Bref * NUFFT_STEPS;
        rtnPlan->F = fftAllocComplexArrayAlign(fSize, FFT_MEM_ALIGNMENT, 
            &rtnPlan->freeF);
        if(rtnPlan->F == NULL) {
            ourRtn = MR_Memory;
            goto errOut;
        }
    }   /* reference implementation only */
    
    if(implOpt) {
        rtnPlan->s = fftAllocComplexArrayAlign(D, FFT_MEM_ALIGNMENT, 
            &rtnPlan->freeS);
        if(rtnPlan->s == NULL) {
            ourRtn = MR_Memory;
            goto errOut;
        }
    }
    
    if(implOpt || implRef) {    
        /* DFT doesn't do any FFTs */
        
        if(rtnPlan->log2D >= NUFFT_MFT_THRESH) {
            /* 
             * Low-level FFTs via MatrixFFT. We call the private mfftCreatePlan1DComplex()
             * so we can disable ThreadPool creation there; the MFFT plan will use
             * our thread pool.
             * We pass in MCF_HintTranspose since we have to do a transpose; that hint
             * flag causes a square matrix config if possible (i.e. if n is even). 
             */
            ourRtn = mfftCreatePlan1DComplex(n, MCF_HintTranspose, numThreads,
                false,      // externalSineTable
                true,       // disableThreads
                false,      // expectColOut since we always cause a transpose
                &rtnPlan->mfftPlan);
            if(ourRtn) {
                printf("***Error creating MatrixFFT plan for log2(N) = %u\n", rtnPlan->log2D);
                goto errOut;
            }
        }
        else {
            /* 
             * Low-level FFTs via vDSP.
             */ 
            #if		!FFT_SPLIT_COMPLEX
            /* We need a vDSP buffer */
            if(fftAllocDSPComplexAlign(&rtnPlan->vdspBuf, D, VECTOR_SIZE, 
                    &rtnPlan->vdspBufFree)) {
                ourRtn = MR_Memory;
                goto errOut;
            }
            #endif

            rtnPlan->vdspSetup = FFTCreateSetup(rtnPlan->log2D);
            if(rtnPlan->vdspSetup == NULL) {
                printf("***Error creating vDSP setup for log2(N) = %u\n", rtnPlan->log2D);
                ourRtn = MR_Memory;
                goto errOut;
            }
        }
    }
    
    if(numThreads == 0) {
        numThreads = numCpuCores();
    }
    if(numThreads > 1) {
        ourRtn = tpThreadInit(&rtnPlan->threadPool, numThreads, 
            0,      // extraPerThread
            sizeof(TP_TaskUnion_U),
            TPO_None);
        if(ourRtn) {
            goto errOut;
        }
        if(rtnPlan->mfftPlan != NULL) {
            /*
             * Config the MatrixFFT plan to use our thread pool.
             */
            rtnPlan->mfftPlan->threadPoolP = &rtnPlan->threadPool;
        }
    }
    
errOut:
    if(ourRtn) {
        nufftFreePlan(rtnPlan);
    }
    else {
        *nuFftPlan = rtnPlan;
    }
    return MR_Success;
}

/*
 * Free an NUFFTPlan.
 */
void nufftFreePlan(
	NUFFTPlan nuFftPlan)
{
    if(nuFftPlan == NULL) {
        return;
    }
    COND_FREE(nuFftPlan->freeTheta);
    COND_FREE(nuFftPlan->mu);
    fftFreeComplexArrayAlign(nuFftPlan->F, nuFftPlan->freeF);
    fftFreeComplexArrayAlign(nuFftPlan->s, nuFftPlan->freeS);
    
	#if	!FFT_SPLIT_COMPLEX
	fftFreeDSPComplex(&nuFftPlan->vdspBufFree);
	#endif

	if(nuFftPlan->vdspSetup) {
		FFTFreeSetup(nuFftPlan->vdspSetup);
		nuFftPlan->vdspSetup = NULL;
	}
    if(nuFftPlan->mfftPlan) {
        mfftFreePlan(nuFftPlan->mfftPlan);
        nuFftPlan->mfftPlan = NULL;
    }
    COND_FREE(nuFftPlan->betaFact);
    
    tpThreadShutdown(&nuFftPlan->threadPool);
    
    free(nuFftPlan);
    return;
}


/* 
 * Execute 1-D Nonuniform FFT.
 */
MFFTReturn nuFftExecute(
	NUFFTPlan       nuFftPlan,
	uint32_t        optFlags,           /* NF_Optimized, etc. */
	FFTComplex		*inBuf,	
    FFTFloat        *tau,  
    FFTComplex      *outBuf) 
{
    if((nuFftPlan == NULL) || (inBuf == NULL) || (tau == NULL) || (outBuf == NULL)) {
        return MR_IllegalArg;
    }
    
    /* Dispatch to alternate implementations */
    switch(optFlags & NF_IMPL_MASK) {
        case NF_Discrete:
            return nuFftExecuteDft(nuFftPlan, inBuf, tau, outBuf);
        case NF_Reference:
            return nuFftExecuteRef(nuFftPlan, inBuf, tau, outBuf);
        case NF_Optimized:
            break;          /* break and execute, below */
        default:
            printf("nuFftExecute: Bad optFlags (must specify exactly one implementation)\n");
            return MR_IllegalArg;
    }
    
    if(!(nuFftPlan->implFlags & NF_Optimized)) {
        printf("***nuFftExecute(NF_Optimized): plan not created for this implementation\n");
        return MR_IllegalArg;
    }
    
    if(fftComplexEquiv(inBuf, outBuf)) {
        printf("***nuFftExecute: Can't perform NUFFT(NF_Optimized) in-place\n");
        return MR_IllegalArg;
    }
    
    /* 
     *  These should have been created in nufftCreatePlan since 
     * (implFlags & NF_Optimized) is true
     */
    RFASSERT(nuFftPlan->mu != NULL);
    RFASSERT(nuFftPlan->theta != NULL);
    RFASSERT(nuFftPlan->s != NULL);
    
    size_t   N = nuFftPlan->D;
    unsigned B = nuFftPlan->Bopt;
    
    /*
     * 1. Initialize:
     *
     *    outBuf   = zeroes
     *    mu[j]    = round(tau[j])
     *    theta[j] = tau[j] - mu[j]
     *    mu[j]    = mu[j] mod N
     */
    PolyComplex outPoly(outBuf);
    outPoly.fill(0.0, N);
    
    for(size_t j=0; j<N; j++) {
        ssize_t mu          = (ssize_t)FFTRoundl(tau[j]);
        nuFftPlan->theta[j] = tau[j] - (FFTFloat)mu;
        nuFftPlan->mu[j]    = ssizeMod2D(mu, nuFftPlan->log2D);
    }
    dumpNuFftFloat("theta", nuFftPlan->theta, N);
    dumpNuFftSize("mu", nuFftPlan->mu, N);
    
    /*
     * 2. Modify input signal:
     *
     *      x[j] = x[j] * e ^ (-i * pi * theta[j])
     */
    nuFftInit(inBuf, nuFftPlan->theta, N);
    
    /* 
     * constants and optimized variables used in the main loop
     */
    FFTFloat oneOverMFact = 1.0;      /* 1/m!, updated at top of loop */
    FFTFloat negTwoPiOverN = -2.0 * M_PI / (FFTFloat)N;

    /*
     * 3. Main loop: perform B standard FFTs.
     */
    for(unsigned m=0; m<B; m++) {
        if(m > 1) {
            oneOverMFact /= (FFTFloat)m;
        }
        
        /* s := zeroes */
        PolyComplex s(nuFftPlan->s);
        s.fill(0.0, N);
   
        /* 
         * s[mu[j]] += x[j] 
         *
         * This little loop takes an inordinate amount of CPU time, mainly because
         * the accesses to the s array are essentially random, so the cache miss behavior 
         * is terrible. Unfortunately we cannot multithread this loop due to the random 
         * updates to s.
         *
         * Doing the theta scaling here instead of at the end of the loop, even when 
         * using multiple threads and vector code there, is a win since we avoid a 
         * cache miss on every inBuf element - there would be one miss to read and 
         * write each element at the end of the loop, plus another miss here to read 
         * (assuming L2 cache overflow it by the end of the theta scale loop). 
         * Instead, we take 1 cache miss here to read x[j], and we scale it by theta, 
         * write it back, and use it to update s[] all at once (i.e. when x[j] is in 
         * L1 cache or even closer). 
         */
        PolyComplex inPoly(inBuf);
        FFTFloat *theta = nuFftPlan->theta;
        
        for(size_t j=0; j<N; j++) {
            if(m > 0) {
                inPoly.scalarMul(*theta++);
            }
            
            RFASSERT(nuFftPlan->mu[j] < N);
            s.set(nuFftPlan->s, nuFftPlan->mu[j]);
            s.add(inPoly.real(), inPoly.imag());            
            ++inPoly;
        }
        
        /* s := FFT(s) */
        ndprintf("--- Calculating FFT for m=%lu\n", (unsigned long)m);
        s.set(nuFftPlan->s);
        fftPoly(nuFftPlan, s, nuFftPlan->log2D, N); 
        
        /* 
         * outBuf[k] += ((-2 * pi * i * (k - N/2) / N) ^ m) * s[k] / m!
         */
        nuFftLoop(nuFftPlan, m, oneOverMFact, negTwoPiOverN, outBuf);
        
        #if		DUMP_NU_FFTS
        char title[100];
        sprintf(title, "outBuf after iteration for m = %lu", (unsigned long)m);
        outPoly.set(outBuf);
        dumpNuFftPC(title, outPoly, N);
        #endif
        
        /*
         * if (m < (B-1)) x[j] := x[j] * theta[j] 
         *
         * This is done at the top of the loop to minimize cache misses
         * on the x array.
         */
    }       /* for m */
        
    return MR_Success;
}

