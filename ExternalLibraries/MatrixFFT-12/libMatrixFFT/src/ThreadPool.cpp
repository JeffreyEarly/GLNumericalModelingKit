/*	File: ThreadPool.cpp  
	
	Description:
		Persistent thread pool module.
	
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
 * ThreadPool.cpp - persistent thread pool module.
 *			  
 * Created Sep. 30 2008.
 * Copyright 2008 by Apple, Inc. 
 */
 
#include "ThreadPool.h"
#include <strings.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mach/thread_policy.h>
#include <mach/thread_act.h>
#include <mach/mach_init.h>

#ifdef	DEBUG
#include <assert.h>
#define TPASSERT(s)		assert(s)
#else
#define TPASSERT(s)
#endif

/*
 * Use Thread Affinity mechanism to place each thread on a different CPU.
 */
static void setThreadAffinity(
    TP_PerThread *pt)
{
    thread_affinity_policy ap;
    ap.affinity_tag = pt->threadNum + 1;    // non-null affinity tag
    
    int ret = thread_policy_set(mach_thread_self(),
          THREAD_AFFINITY_POLICY,
          (integer_t*) &ap,
          THREAD_AFFINITY_POLICY_COUNT);
    if(ret) {
        printf("setThreadAffinity: thread_policy_set returned %d\n", ret);
        /* oh well */
    }
}

/* 
 * Synchronization protocol for threads managed here:
 *
 * -- Each pthread has an associated TP_PerThread, a pointer to which is passed as
 *    the void * pthread argument.
 * -- A thread's TP_PerThread.task contains info about app-specific work to do.
 * -- When main thread has work for a thread to do, it does this:
 *     -- fill out TP_PerThread.task (done by clients of this module; remainder 
 *        done here)
 *     -- TP_PerThread.state = TPS_Working
 *     -- pthread_cond_broadcast(TP_PerThread.cond)
 *     -- possibly execute some other stuff in the main thread, then
 *     -- pthread_cond_wait(TP_PerThread.cond) until TP_PerThread.state = TPS_Idle
 *
 * Each of our managed pthreads does this:
 *     while(1) {
 *         pthread_cond_wait(TP_PerThread.cond) until TP_PerThread.state = TPS_Working;
 *         do the work (by calling TP_Task.threadFcn) or exit if op == TPO_Exit;
 *		   TP_PerThread.state = TPS_Idle;
 *         pthread_cond_broadcast(TP_PerThread.cond);
 *     }
 */
 
/* 
 * The worker thread.
 */
static void *tpThread(
	void *arg)
{
	TP_PerThread *pt = (TP_PerThread *)arg;
	TP_Task *task = &pt->task;
	
    if(pt->threadOptions & TPO_SeparateAffinity) {
        /*
         * Try to run each thread on a separate CPU core
         */
        setThreadAffinity(pt);
    }
    
	while(1) {
		tpThreadDebug("tpThread: thread %u waiting for work\n", pt->threadNum);
		
		/* First wait for work to do */
		if(pthread_mutex_lock(&pt->mutex)) {
			printf("***tpThread: Error acquiring lock; aborting.\n");
			/* What do we do now? */
			pthread_exit(NULL);
		}
		while(pt->state == TPS_Idle) {
			int rtn = pthread_cond_wait(&pt->cond, &pt->mutex);
			if(rtn) {
				printf("***tpThread: Error waiting on condition; error %d; aborting.\n", rtn);
				pthread_exit(NULL);
			}
		}
		pthread_mutex_unlock(&pt->mutex);
		
		/* do the work */
		switch(task->op) {
			case TPO_Exit:
				tpThreadDebug("tpThread: thread %u exiting\n", pt->threadNum);
				// hangs! --> pthread_exit(NULL);
				return NULL;
				/* NOT REACHED */
			default:
				tpThreadDebug("tpThread: thread %u dispatch\n", pt->threadNum);
				THR_TSTAMP(pt->startTime);
				task->status = task->threadFcn(task->u);
				THR_TSTAMP(pt->endTime);
				THR_ACCUM(pt->startTime, pt->endTime, pt->totalTime);
				tpThreadDebug("tpThread: thread %u return\n",	pt->threadNum);
				tpLogThreadSingle("=== thread %u: elapsed time %.2f ms\n", pt->threadNum,
					THR_ELAPS_MS(pt->startTime, pt->endTime));
				break;
		}
		tpThreadDebug("tpThread: thread %u finished work, notifying main thread\n",
			pt->threadNum);
		
		/* now let main thread know we're done */
		if(pthread_mutex_lock(&pt->mutex)) {
			printf("***Error acquiring lock; aborting.\n");
			/* What do we do now? */
			pthread_exit(NULL);
		}
		pt->state = TPS_Idle;
		if(pthread_cond_broadcast(&pt->cond)) {
			printf("***Error waking main thread; aborting.\n");
			pthread_exit(NULL);
		}
		if(pthread_mutex_unlock(&pt->mutex)) {
			printf("***Error acquiring server lock; aborting.\n");
			pthread_exit(NULL);
		}
	}
	/* NOT REACHED */
	return NULL;
}

/* 
 * Initialize a TP_ThreadPool.
 * Returns nonzero on error.
 */
MFFTReturn tpThreadInit(
	TP_ThreadPool	*threadPool,
	unsigned		numPthreads,
	unsigned		extraPerThreads,	/* optional addtional perThread data to allocate */
	unsigned		taskUnionBytes,		/* sizeof(TP_TaskUnion) */
    uint32_t        threadOptions)      /* TPO_xxx */
{
	threadPool->numThreads = numPthreads;
	threadPool->numPerThreads = 0;
	threadPool->perThread = NULL;
	if(numPthreads == 0) {
		return MR_Success;
	}

	threadPool->numThreads = numPthreads;
	threadPool->numPerThreads = numPthreads + extraPerThreads;
	size_t mallocSize = threadPool->numPerThreads * sizeof(TP_PerThread);
	threadPool->perThread = (TP_PerThread *)malloc(mallocSize);
	if(threadPool->perThread == NULL) {
		printf("***tpThreadInit: malloc failure\n");
		return MR_Memory;
	}
	/* subsequent errors to errOut: */
	memset(threadPool->perThread, 0, mallocSize);
	
	for(unsigned dex=0; dex<threadPool->numPerThreads; dex++) {
		TP_PerThread *pt = &threadPool->perThread[dex];
		pt->task.u = (TP_TaskUnion *)malloc(taskUnionBytes);
		if(pt->task.u == NULL) {
			printf("***tpThreadInit: malloc failure\n");
			return MR_Memory;
		}
		memset(pt->task.u, 0, taskUnionBytes);
		
		pt->state = TPS_Idle;
		pt->threadNum = dex;
		if(dex < threadPool->numThreads) {
			/* skip this for TP_PerThreads with no actual pthread associated with them */
			if(pthread_mutex_init(&pt->mutex, NULL)) {
				printf("***tpThreadInit: Error initializing mutex\n");
				goto errOut;
			}
			if(pthread_cond_init(&pt->cond, NULL)) {
				printf("***tpThreadInit: Error initializing pthreadCond\n");
				goto errOut;
			}
			if(pthread_create(&pt->thr, NULL, tpThread, pt)) {
				printf("***tpThreadInit: Error starting up server thread\n");
				goto errOut;
			}
		}
		pt->threadOptions = threadOptions;
        
		/*
		 * per-thread numSubmatrices coefficient. See comments at 
		 * fftSwapRowsColumnsThr().
		 * A bit of a hack, having this implementation-dependent code
		 * here. Maybe this should be a callback, or otherwise up
		 * to the caller. 
		 */
		if(dex < threadPool->numThreads) {
			float fract = (float)(threadPool->numThreads - dex) / (float)threadPool->numThreads;
			pt->squareCoeff = 1.0 - sqrtf(fract);
		}
	}
	return MR_Success;
	
errOut:
	free(threadPool->perThread);
	threadPool->perThread = NULL;
	return MR_Internal;
}

/* 
 * Clean up per-thread state.
 */
void tpThreadShutdown(
	TP_ThreadPool		*threadPool)
{
	if(threadPool->numThreads == 0) {
		return;
	}
	tpThreadDebug("tpThreadShutdown: shutting down %u threads\n",
		threadPool->numThreads);
		
	for(unsigned dex=0; dex<threadPool->numPerThreads; dex++) {
		TP_PerThread *pt = &threadPool->perThread[dex];
		TP_Task *task = &pt->task;
		
		if(dex < threadPool->numThreads) {
			/* ignore errors */
			pthread_mutex_lock(&pt->mutex);
			task->op = TPO_Exit;
			pt->state = TPS_Working;
			pthread_cond_broadcast(&pt->cond);
			pthread_mutex_unlock(&pt->mutex);
		}
		free(task->u);
		task->u = NULL;
	}
	
	/* 
	 * Wait for all threads to exit before we free up the resources 
	 * they're using.
	 */
	for(unsigned dex=0; dex<threadPool->numThreads; dex++) {
		TP_PerThread *pt = &threadPool->perThread[dex];
		void *thrRtn;
		pthread_join(pt->thr, &thrRtn);
		tpThreadDebug("tpThreadShutdown: thread %u exited\n", dex);
	}
	
	/* now free up the resources we created in tpThreadInit() */
	for(unsigned dex=0; dex<threadPool->numThreads; dex++) {
		TP_PerThread *pt = &threadPool->perThread[dex];
		pthread_mutex_destroy(&pt->mutex);
		pthread_cond_destroy(&pt->cond);
	}
	
	free(threadPool->perThread);
	threadPool->perThread = NULL;
}

/*
 * Dispatch tasks to specified number of threads. 
 * Caller has set up threadPool.perThread->task.u for each thread. 
 * Returns nonzero on error.
 */
MFFTReturn tpThreadDispatch(
	TP_ThreadPool		*threadPool,
	unsigned			numThreads)		/* might be different from threadPool->numThreads */
{
	tpThreadDebug("tpThreadDispatch: sending tasks to %u threads\n", numThreads);
	TPASSERT((threadPool->numThreads > 0) && 
			 (numThreads <= threadPool->numThreads));
	
	for(unsigned dex=0; dex<numThreads; dex++) {
		TP_PerThread *pt = &threadPool->perThread[dex];
		TP_Task *task = &pt->task;
		
		/* worker threads sets this on completion of task */
		task->status = MR_Internal;
		
		if(pthread_mutex_lock(&pt->mutex)) {
			printf("tpThreadDispatch: mutex error\n");
			return MR_Internal;
		}
		
		/* this thread better not be running! */
		TPASSERT(pt->state == TPS_Idle);
		if(pt->state != TPS_Idle) {
			printf("tpThreadDispatch: thread %u not idle\n", dex);
			return MR_Internal;
		}
		
		pt->state = TPS_Working;
		if(pthread_cond_broadcast(&pt->cond)) {
			printf("tpThreadDispatch: cond_broadcast error\n");
			return MR_Internal;
		}
		pthread_mutex_unlock(&pt->mutex);
	}
	return MR_Success;
}


/* 
 * Wait for numThreads to finish pending ops.
 * Returns first nonzero per-thread status found, else
 * returns MR_Success.
 */
MFFTReturn tpThreadFinish(
	TP_ThreadPool		*threadPool,
	unsigned			numThreads)		/* might be different from threadPool->numThreads */
{
	MFFTReturn ourRtn = MR_Success;
	THR_DECL(startTime);
	THR_DECL(endTime);

	THR_TSTAMP(startTime);
	for(unsigned dex=0; dex<numThreads; dex++) {
		TP_PerThread *pt = &threadPool->perThread[dex];
		TP_Task *task = &pt->task;

		tpThreadDebug("tpThreadFinish: waiting for thread %u\n", dex);
		if(pthread_mutex_lock(&pt->mutex)) {
			printf("***tpThreadFinish: Error acquiring lock\n");
			return MR_Internal;
		}
		while(pt->state == TPS_Working) {
			if(pthread_cond_wait(&pt->cond, &pt->mutex)) {
				printf("***tpThreadFinish: Error waiting on condition\n");
				return MR_Internal;
			}
		}
		pthread_mutex_unlock(&pt->mutex);
		tpThreadDebug("tpThreadFinish: thread %u FINISHED, status %d\n", dex, task->status);
		if((task->status != MR_Success) && (ourRtn == MR_Success)) {
			/* first error, report this as our return value */
			ourRtn = task->status;
		}
	}
	THR_TSTAMP(endTime);
	tpLogThreadSingle("=== thread M: time waiting for other threads %.2f ms\n", 
			THR_ELAPS_MS(startTime, endTime));
	return ourRtn;
}

#if		FFT_THREAD_TIME

/* reset all cumulative timers */
void tpThreadResetTimers(
	TP_ThreadPool		*threadPool,
	unsigned			numThreads)
{
	for(unsigned dex=0; dex<numThreads; dex++) {
		TP_PerThread *pt = &threadPool->perThread[dex];
		pt->totalTime = 0.0;
	}
}

/* log specified thread's cumulative time */
void tpThreadLogTime(
	TP_ThreadPool		*threadPool,
	unsigned			threadNum,
	const char			*title)
{
	tpLogTotalTime(title, threadPool->perThread[threadNum].totalTime);
}

#endif
