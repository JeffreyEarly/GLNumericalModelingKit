/*	File: NUFFT_Ref.cpp
	
	Description:
		Unoptimized reference implementations of Nonuniform FFT and DFT. 
	
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
 * NUFFT_Ref.cpp - Unoptimized reference implementations of Nonuniform FFT and DFT.
 *			  
 * Created 04/30/2009. 
 * Copyright 2009 by Apple, Inc. 
 */
 
#include <libMatrixFFT/NUFFT.h>
#include <libMatrixFFT/fftUtils.h>
#include "NUFFT_Priv.h"
#include "NUFFT_Utils.h"
#include "fftDebug.h"
#include <stdlib.h>

#pragma mark --- Unoptimized Nonuniform FFT ---

MFFTReturn nuFftExecuteRef(
	NUFFTPlan       nuFftPlan,
	FFTComplex		*inBuf,	
    FFTFloat        *tau,  
    FFTComplex      *outBuf) 
{
    if(!(nuFftPlan->implFlags & NF_Reference)) {
        printf("***nuFftExecute(NF_Reference): plan not created for this implementation\n");
        return MR_IllegalArg;
    }

    /* Caller validated this stuff */
    RFASSERT(nuFftPlan != NULL);
    RFASSERT(inBuf != NULL);
    RFASSERT(tau != NULL);
    RFASSERT(outBuf != NULL);
    
    /* 
     * And these should have been created in nufftCreatePlan since 
     * (implFlags & NF_Reference) is true
     */
    RFASSERT(nuFftPlan->mu != NULL);
    RFASSERT(nuFftPlan->theta != NULL);
    RFASSERT(nuFftPlan->F != NULL);
    RFASSERT(nuFftPlan->betaFact != NULL);
    
    /* 
     * Note internal calculation is double precision regardless of 
     * configuration of the library
     */
    size_t   D = nuFftPlan->D;
    unsigned B = nuFftPlan->Bref;
    double dFloat = D;
    
    /*
     * 1. Initialize:
     *
     *    tauModD  = tau[j] % D
     *    mu[j]    = round(tauModD)
     *    theta[j] = tauModD - mu[j]
     *    mu[j]    = mu[j] mod D
     */
    for(size_t j=0; j<D; j++) {
        FFTFloat tauModD = realMod2D(tau[j], nuFftPlan->log2D);
        nuFftPlan->mu[j]    = (size_t)FFTRoundl(tauModD);
        nuFftPlan->theta[j] = tauModD - nuFftPlan->mu[j];
        if(nuFftPlan->mu[j] == D) {
            /* the only way mu[j] mod D would NOT be mu[j] */
            nuFftPlan->mu[j] = 0;
        }
    }
    
    /*
     * 2. Perform a total of 8B standard FFTs.
     */
    for(unsigned K=0; K<NUFFT_STEPS; K++) {
    
        /*
         * Partial exponent used below, const for this K
         */
        double exponPart = -2.0 * M_PI * K / NUFFT_STEPS;;

        for(unsigned beta=0; beta<B; beta++) {
        
            /* 
             * Working s array at this point is F[K][beta]
             */
            size_t fOff = (K * B * D) + (beta * D);
            PolyComplex s(nuFftPlan->F, fOff); 
            
            /* s := Zero signal of length D. */
            s.fill(0.0, D);
            
            /*
             * for j = 0..D
             *    u = mu[j] mod D
             *    s[u] += x[j] * (e ^ (-2 * pi * i * K * tau[j] / NUFFT_STEPS)) *
             *                theta[j]^beta
             */
            for(size_t j=0; j<D; j++) {
                size_t u = nuFftPlan->mu[j];
                PolyComplex su(nuFftPlan->F, fOff + u);
                
                /* 
                 * expTerm := e ^ (-2 * pi * i * K * tau[j] / NUFFT_STEPS) 
                 * expon    = -2 * pi * K * tau[j] / NUFFT_STEPS
                 * expTerm  = cos(expon) + (i * sin(expon))
                 *
                 * Since expTerm includes a 2*pi/8 term, we take tau[j]
                 * mod 8 here. 
                 */
                double tauMod8 = realMod2D(tau[j], NUFFT_LOG2_STEPS);
                double expon    = exponPart * tauMod8;
                double expTermR = cos(expon);
                double expTermI = sin(expon);
                
                /* xe := x[j] * expTerm */
                PolyComplex pe(inBuf, j);
                double xeR = (pe.real() * expTermR) - (pe.imag() * expTermI);
                double xeI = (pe.real() * expTermI) + (pe.imag() * expTermR);
                
                /* xe *= theta[j]^beta */
                double thetaToBeta = fExpI(nuFftPlan->theta[j], beta);
				// nprintf("theta[%lu] = %f  beta %u  theta^beta %f\n",
                //       (unsigned long)j, nuFftPlan->theta[j], beta, thetaToBeta);

                xeR *= thetaToBeta;
                xeI *= thetaToBeta;
                
                /* su[j] += x[j] * (e ^.... */
                su.add(xeR, xeI);
            }   /* for j */
            
            /*
             * F[K][beta] = s = FFT(s)
             */
			ndprintf("--- Calculating F[%u][%u] at %p\n", K, beta, s.realP());
            fftPoly(nuFftPlan, s, nuFftPlan->log2D, nuFftPlan->D); 
        }   /* for beta */
    }   /* for K */
    
   	size_t dOver8 = D >> 3;

    #if		DUMP_NU_FFTS
    for(unsigned K=0; K<NUFFT_STEPS; K++) {
        for(unsigned beta=0; beta<B; beta++) {
    		size_t thisF = (K * B * D) + (beta * D);
			char title[100];
			PolyComplex thisPc(nuFftPlan->F, thisF);
			sprintf(title, "F[%u][%u] at %p", K, beta, thisPc.realP());
			dumpNuFftPC(title, thisPc, dOver8);
		}
	}
	#endif

    /* 
     * Final step here, assemble result from F[][][].
     * The result is the concatentation of NUFFT_STEPS sets, each set D/8 elements. 
     */
	double neg2PiOverD = -2.0 * M_PI / dFloat;

    for(unsigned K=0; K<NUFFT_STEPS; K++) {
		PolyComplex currOut(outBuf, dOver8 * K);	
        for(size_t m=0; m<dOver8; m++) {
			/*
 			 * currOut = Sum[  f[k][bet][[m]] pow[-2 Pi I m/n, bet] 1/bet!, 
			 *    {bet, 0, B - 1}];
			 *
			 */
			double currR = 0.0;
			double currI = 0.0;
			double fBase = neg2PiOverD * (double)m;

			/* generate sum into currR, currI */
        	for(unsigned beta=0; beta<B; beta++) {
				/* powTerm = (-2 * pi * i * m / D) ^ beta */
				FFTFloat powTermR;
				FFTFloat powTermI;
				fiExpI(fBase, beta, &powTermR, &powTermI);

				/* powTerm /= beta! */
				double betaFact = nuFftPlan->betaFact[beta];
				powTermR *= betaFact;
				powTermI *= betaFact;

				/* get to F[k][beta][m] */
				size_t fOff = (K * B * D) + (beta * D) + m;
				PolyComplex currF(nuFftPlan->F, fOff);

				/* currTerm = F[k][beta][m] * powTerm */
				double ctR = (currF.real() * powTermR) - (currF.imag() * powTermI);
				double ctI = (currF.real() * powTermI) + (currF.imag() * powTermR);

				/* curr += ct */
				currR += ctR;
				currI += ctI;
			} /* for beta */

			/* Store current term */
			currOut.real(currR);
			currOut.imag(currI);
			++currOut;
            
		}	/* for m */
	} /* For K */
	
    return MR_Success;
}


#pragma mark --- Nonuniform DFT ---

/* ThreadPool callout. */
static MFFTReturn dftThreadCall(
	TP_TaskUnion *u)
{
    TP_NUDFT *dft = &u->nuDft;
    
    RFASSERT(dft->inBuf != NULL);
    RFASSERT(dft->outBuf != NULL);
    RFASSERT(dft->tau != NULL);
    
    PolyComplex X(dft->outBuf, dft->startElt);
    
    /* 
     * Note internal calculation is double precision regardless of 
     * configuration of the library
     */
    double constExpPart = -2.0 * M_PI / (double)dft->D;
    size_t endElt = dft->startElt + dft->eltsToProcess;
    
    for(size_t k=dft->startElt; k<endElt; k++) {
        double sumR = 0.0;
        double sumI = 0.0;
        PolyComplex x(dft->inBuf);
        double expPart = constExpPart * (FFTFloat)k;
        for(size_t j=0; j<dft->D; j++) {
            /*
             * expTerm := e ^ (-2 * pi * i * k * w[j] / D)
             * expon    = -2 * pi * k * w[j] / D
             * expTerm  = cos(expon) + (i * sin(expon))
             */
            double expon    = expPart * dft->tau[j];
            double expTermR = cos(expon);
            double expTermI = sin(expon);
            
            /* xe := x[j] * expTerm */
            double xeR = (x.real() * expTermR) - (x.imag() * expTermI);
            double xeI = (x.real() * expTermI) + (x.imag() * expTermR);
            ++x;
            
            /* Add in to running sum for X[k] */
            sumR += xeR;
            sumI += xeI;
        }
        X.real(sumR);
        X.imag(sumI);
        ++X;
    }
    return MR_Success;
}

#define DFT_THREAD_MIN  64      /* min D to use threading below */

/* 
 * Brute force, literal nonuniform DFT. 
 */
MFFTReturn nuFftExecuteDft(
	NUFFTPlan       nuFftPlan,
	FFTComplex		*inBuf,	
    FFTFloat        *tau,   
    FFTComplex      *outBuf)
{
    if(fftComplexEquiv(inBuf, outBuf)) {
        printf("***nuFftExecuteDft: Can't perform NUDFT in-place\n");
        return MR_IllegalArg;
    }
    
    /* Caller validated this stuff */
    RFASSERT(nuFftPlan != NULL);
    RFASSERT(inBuf != NULL);
    RFASSERT(tau != NULL);
    RFASSERT(outBuf != NULL);
    
    unsigned numThreads = nuFftPlan->threadPool.numThreads;
    if((nuFftPlan->D <= DFT_THREAD_MIN) ||
       (numThreads <= 1)) {
        /* Threading disabled */
        TP_TaskUnion tu;
        TP_NUDFT *dft       = &tu.nuDft;
        dft->inBuf          = inBuf;
        dft->tau            = tau;
        dft->outBuf         = outBuf;
        dft->D              = nuFftPlan->D;
        dft->startElt       = 0;
        dft->eltsToProcess  = dft->D;
        return dftThreadCall(&tu);
    }
    
    size_t eltsPerThread = nuFftPlan->D / numThreads;
    size_t eltsToGo = nuFftPlan->D;
    size_t currElt = 0;
    
    for(unsigned dex=0; dex<numThreads; dex++) {
		TP_PerThread *pt = &nuFftPlan->threadPool.perThread[dex];
		TP_Task *task    = &pt->task;
		TP_NUDFT *dft    = &task->u->nuDft;

		task->op		 = TPO_App;
		task->threadFcn  = dftThreadCall;
		dft->inBuf		 = inBuf;
		dft->tau         = tau;
        dft->D           = nuFftPlan->D;
        dft->outBuf      = outBuf;
        dft->startElt    = currElt;
        if(dex == (numThreads - 1)) {
            /* Might be more for this one, use all remaining */
            dft->eltsToProcess = eltsToGo;
        }
        else {
            dft->eltsToProcess = eltsPerThread;
        }
        eltsToGo -= dft->eltsToProcess;
        currElt  += dft->eltsToProcess;
    }
    
	/* GO */
	MFFTReturn ourRtn = tpThreadDispatch(&nuFftPlan->threadPool, numThreads);
	if(ourRtn) {
		return ourRtn;
	}
	ourRtn = tpThreadFinish(&nuFftPlan->threadPool, numThreads);
	return ourRtn;

}

