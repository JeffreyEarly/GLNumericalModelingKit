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
 * ThreadPool.h - persistent thread pool module. 
 *			  
 * Created Sep. 30 2008.
 * Copyright 2008 by Apple, Inc. 
 */
 
#ifndef	_THREAD_POOL_H_
#define _THREAD_POOL_H_

#include <pthread.h>

/* This #include only needed for MFFTReturn typedef */
#include <libMatrixFFT/MatrixFFT.h>

#ifdef __cplusplus
extern "C" {
#endif

#pragma mark --- debugging flags ---

/* dump ThreadPool-related info to stdout */
#define FFT_THREAD_DEBUG			0

/* enable per-thread timing */
#define FFT_THREAD_TIME				0

/* when true, show individual worker threads' single task times */
#define FFT_THREAD_TIME_SINGLE		0


#pragma mark --- Application-specific typedefs ---

/* 
 * Application-specific code defines this elsewhere; its specific definition doesn't
 * have to be visible here. 
 * Its size is passed to tpThreadInit() for app-specific per-thread allocation.
 * A pointer to the per-thread version of this is passed to the thread callout 
 * function to define work to be done. 
 */
union TP_TaskUnion_U;
typedef union TP_TaskUnion_U TP_TaskUnion;

#pragma mark --- ThreadPool typedefs ---

/* 
 * Each loop of a thread in this module performs one of these operations on
 * each loop - either perform some app-specific work, or exit.
 */
typedef enum {
    TPO_App,        /* app-specific work */
	TPO_Exit		/* pthread_exit() */
} TP_Op;

/*
 * tpThreadInit options. Also stored in TP_PerThread. 
 *
 * When TPO_SeparateAffinity is specified, the Thread Affinity mechanism
 * will be used to (attempt to) run each thread on a separate CPU core. 
 */
#define TPO_None                0x0000      /* no options / default */
#define TPO_SeparateAffinity    0x0001   

/* 
 * Client-provided callback, which is called from worker threads to
 * perform the actual work. 
 */
typedef MFFTReturn (*tpThreadFcn)(TP_TaskUnion *u);

/* 
 * Definition of one task; each thread performs one of these per loop.
 */
typedef struct {
	/* input, written by client, varies per task */
	TP_Op				op;
	TP_TaskUnion		*u;
	tpThreadFcn			threadFcn;
	
	/* 
	 * Output, owned by tpThread module, nonzero is error,
	 * reported via tpThreadFinish() 
	 */
	MFFTReturn			status;
} TP_Task;

/* 
 * Values for the semaphore used in synchronizing threads. 
 */
typedef enum {
	TPS_Idle = 0,		/* thread is idle, no work to do */
	TPS_Working			/* work to do */
} TP_State;

/*
 * Per-thread state, in TP_ThreadPool.perThread[].
 */
typedef struct {
	TP_State		state;
	pthread_mutex_t	mutex;
	pthread_cond_t	cond;
	pthread_t		thr;
	unsigned		threadNum;		
	float			squareCoeff;	/* for divvying up rows in a square operation */
	TP_Task			task;			/* per-task info */
    uint32_t        threadOptions;  /* TPO_xxx */
	#if		FFT_THREAD_TIME
	double			startTime;
	double			endTime;
	double			totalTime;
	#endif
} TP_PerThread;

/*
 * Public ThreadPool.
 */
typedef struct {
	unsigned		numThreads;		/* # of pthreads managed here */
	unsigned		numPerThreads;	/* size of perThread array, can be > numThreads */
	TP_PerThread	*perThread;		/* size numThreads+1 to provide main() a place to 
									 * store parameters for its work */
} TP_ThreadPool;

#pragma mark --- Public ThreadPool API ---

/*
 * The general scheme for using this module is:
 *
 * -- Call tpThreadInit() to set up specific number of pthreads, which 
 *    persist throughout the lifetime of a TP_ThreadPool object. 
 * -- to perform a threading op:
 *    -- initialize TP_Task.{op,u,threadFcn} for each
 *       element in TP_ThreadPool->perThread[] you want to use;
 *    -- call tpThreadDispatch() to start up numThreads threads 
 *       (numThreads might be less than TP_ThreadPool.numThreads);
 *    -- possibly do some work in the main thread;
 *    -- call tpThreadFinish() to wait for all threads to
 *       complete;
 */
 
/* 
 * Initialize a TP_ThreadPool.
 * The numPthreads argument here is the number of pthreads this object manages.
 * If you do some of your work in main(), and you want a total of n threads to 
 * do all the work, then call this with numPthreads == n-1.
 * We allocate additional TP_ThreadPool.perThread per the specified
 * extraPerThreads. If that's zero then the size of the perThread array is 
 * numPthreads. Note we NEVER allocate any perThread data if numPthreads is 0.
 * Returns nonzero on error.
 */
extern MFFTReturn tpThreadInit(
	TP_ThreadPool		*threadPool,
	unsigned			numPthreads,
	unsigned			extraPerThreads,	/* optional addtional perThread data to allocate */
	unsigned			taskUnionBytes,     /* sizeof(TP_TaskUnion) */
    uint32_t            threadOptions);     /* TPO_xxx, above */

/* 
 * Clean up a TP_ThreadPool. Shut down all threads and free memory.
 */
extern void tpThreadShutdown(
	TP_ThreadPool		*threadPool);
	
/*
 * Dispatch tasks to specified number of threads. 
 * Caller has set up threadPool->perThread->task.u for each thread. 
 * Returns nonzero on error.
 */
MFFTReturn tpThreadDispatch(
	TP_ThreadPool		*threadPool,
	unsigned			numThreads);	

/* 
 * Wait for numThreads to finish pending ops.
 * Returns first nonzero per-thread status found, else
 * returns MR_Success.
 */
MFFTReturn tpThreadFinish(
	TP_ThreadPool		*threadPool,
	unsigned			numThreads);		/* should be same as passed to tpThreadDispatch()  */

#pragma mark --- debugging and timing support ---

#if		FFT_THREAD_DEBUG
#define tpThreadDebug(s...)		printf(s)
#else
#define tpThreadDebug(s...)
#endif

#if		FFT_THREAD_TIME

#include <CoreFoundation/CoreFoundation.h>
#include <stdio.h>

#define tpLogThreadTime(s...)		printf(s)
#if		FFT_THREAD_TIME_SINGLE
#define tpLogThreadSingle(s...)		printf(s)
#else
#define tpLogThreadSingle(s...)
#endif	/* FFT_THREAD_TIME_SINGLE */
#define THR_DECL(var)				double var; var = 0.0
#define THR_TSTAMP(var)				var = CFAbsoluteTimeGetCurrent()
#define THR_ACCUM(start, end, cum)	cum += (end - start)
#define THR_RESET(cum)				cum = 0.0
#else

#define tpLogThreadTime(s...)
#define tpLogThreadSingle(s...)
#define THR_DECL(var)		
#define THR_TSTAMP(var)		
#define THR_ACCUM(start, end, cum)
#define THR_RESET(cum)			

#endif	/* FFT_THREAD_TIME */

#define THR_MS(t)					((t) * 1000.0)
#define THR_ELAPS_MS(start, end)	(THR_MS(end - start))


#if		FFT_THREAD_TIME

/* reset all cumulative timers */
extern void tpThreadResetTimers(
	TP_ThreadPool		*threadPool,
	unsigned			numThreads);		/* might be different from fftRowSetup->numThreads */

#define tpLogTotalTime(title, elaps)		printf("=== %s: %.2f ms\n", \
												title, THR_ELAPS_MS(0, elaps))

/* log specified thread's cumulative time */
extern void tpThreadLogTime(
	TP_ThreadPool		*threadPool,
	unsigned			threadNum,
	const char			*title);
	
/* 
 * add specified time to a thread's cumulative time.
 * Used when the main thread finished up work chiefly performed by a thread. 
 */
#define	tpAddThreadTime(threadPool, threadNum, elaps)	\
	(threadPool)->perThread[threadNum].totalTime += elaps
	
#else	/* !FFT_THREAD_TIME */
#define tpThreadResetTimers(threadPool, numThreads)
#define tpLogTotalTime(title, elaps)
#define tpThreadLogTime(threadPool, tn, title)
#define	tpAddThreadTime(threadPool, threadNum, elaps)
#endif	/* FFT_THREAD_TIME */

#ifdef __cplusplus
}
#endif

#endif	/* _THREAD_POOL_H_ */
