/*	File: fftUtils.cpp 
	
	Description:
		Common utility functions for FFT tests.
	
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
 * fftUtils.cpp - common utility functions for FFT tests.
 *
 * Created 7/17/2007. 
 * Copyright 2008 by Apple, Inc. 
 */

#include <stdlib.h>
#include <strings.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <libMatrixFFT/fftUtils.h>
#include <libMatrixFFT/devRandom.h>
#include "fftPriv.h"
#include "fileIo.h"
#include "fftDebug.h"
#include <CoreFoundation/CoreFoundation.h>
#include <mach/vm_param.h>
#include <sys/sysctl.h>

#pragma mark --- Test setup and logging ---

/*
 * Obtain CPU info. 
 */

#define SYSCTL_NAME		"machdep.cpu.brand_string"
#define SYSCTL_BUF_SIZE	1024

static void fftPrintCpuInfo()
{
	size_t len = SYSCTL_BUF_SIZE;
	char sysctlBuf[SYSCTL_BUF_SIZE];
	
	if(sysctlbyname(SYSCTL_NAME, sysctlBuf, &len, NULL, 0)) {
		printf("<unknown>\n");
		return;
	}
	
	/* 
	 * output of that sysctl looks like this:
	 *
	 * Genuine Intel(R) CPU            1500  @ 2.00GHz
	 */
	char *str = sysctlBuf;
	len = strlen(str);

	bool lastWasSpace = false;
	for(unsigned dex=0; dex<len; dex++) {
		char c = *str++;
		if(c == ' ') {
			/* compress multiple spaces into one */
			if(!lastWasSpace) {
				putchar(c);
			}
			lastWasSpace = true;
		}
		else {
			lastWasSpace = false;
			putchar(c);
		}
	}
	putchar('\n');
}

static const char *months[] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun", 
							   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };
							   
static void fftPrintCurrentTime()
{
	CFTimeZoneRef timeZone = CFTimeZoneCopyDefault();
	CFAbsoluteTime currTime = CFAbsoluteTimeGetCurrent();
	CFGregorianDate gd = CFAbsoluteTimeGetGregorianDate(currTime, timeZone);
	int is = (int)gd.second;
	printf("Test time   : %02d:%02d:%02d  %s %d %lu\n",
		gd.hour, gd.minute, is,
		months[gd.month - 1], gd.day, (unsigned long)gd.year);
}

/* print the contents of a CFString */
static void printCfStr(
	CFStringRef cfstr)
{
	CFDataRef strData = CFStringCreateExternalRepresentation(NULL, cfstr,
		kCFStringEncodingUTF8, true);
	if(strData == NULL) {
		printf("<<string decode error>>");
		return;
	}
	const char *cp = (const char *)CFDataGetBytePtr(strData);
	CFIndex len = CFDataGetLength(strData);
	for(CFIndex dex=0; dex<len; dex++) {
		putchar(*cp++);
	}
	CFRelease(strData);
}

/* This is declared in CFPriv.h, which we might not have access to */
extern "C" CFStringRef CFCopySystemVersionString(void);

static void printSystemVersion()
{
	CFStringRef versStr = CFCopySystemVersionString();
	if(versStr == NULL) {
		printf("<unknown>\n");
		return;
	}
	printCfStr(versStr);
	printf("\n");
	CFRelease(versStr);
}

void fftPrintTestBanner(
	const char *fftType,		/* e.g. "One-dimension real" */
	const char *impl,			/* e.g. "Accelerate", "FFTW" */
	bool doublePrec,
	const char *signalType,		/* e.g. "random", "chirp" */
	const char *options,		/* optional */
	unsigned loops,
	unsigned numThreads)
{
	printf("Library     : %s\n", impl);
	printf("FFT type    : %s\n", fftType);
	printf("Precision   : %s\n", doublePrec ? "Double" : "Single");
	printf("Signal type : %s\n", signalType);
	if(loops) {
		printf("Loops       : %u\n", loops);
	}
	if(numThreads > 1) {
		printf("Threads     : %u\n", numThreads);
	}
	if(options && (options[0] != '\0')) {
		printf("Options     : %s\n", options);
	}
	printf("System      : ");
	printSystemVersion();
	printf("CPU         : ");
	fftPrintCpuInfo();
	printf("Built as    : ");
	#if		defined(__x86_64__) || defined(__ppc64__)
	printf("64 bit\n");
	#else
	printf("32 bit\n");
	#endif
	printf("Complex     : %s\n", FFT_SPLIT_COMPLEX ? "Split" : "Interleaved");
	fftPrintCurrentTime();
}

/* 
 * Get current OS-measured "user time", or absolute (wall) time, in seconds, 
 * as a double.
 */
double fftGetTime(bool wallTime)		/* true --> CFAbsoluteTimeGetCurrent() 
										 * false --> user time via getrusage() */
{
	if(wallTime) {
		return CFAbsoluteTimeGetCurrent();
	}
	
	/* user time from kernel */
	struct rusage ru;
	if(getrusage(RUSAGE_SELF, &ru)) {
		perror("getrusage");
		return 0.0;
	}
	struct timeval *tvp = &ru.ru_utime;
	double rtn = (double)tvp->tv_sec;
	rtn += (double)tvp->tv_usec / 1000000.0;
	return rtn;
}

void appendOptStr(
	char *optStr,
	const char *newOpt)
{
	if(optStr[0] != '\0') {
		strcat(optStr, ", ");
		strcat(optStr, newOpt);
	}
	else {
		strcpy(optStr, newOpt);
	}
}

/* 
 * Obtain a compact string representation, in K or M as appropriate,
 * of specified input size_t. Caller must free() the result.
 */
#define MAX_OUT_LEN	128

char *fftStringRep(
	size_t inSize)
{
	char *outStr = (char *)malloc(MAX_OUT_LEN);
	if(outStr == NULL) {
		return strdup("Malloc error");
	}
	if(inSize < ((size_t)1 << 10)) {
		snprintf(outStr, MAX_OUT_LEN, "%llu", (unsigned long long)inSize);
		return outStr;
	}
	else if(inSize < ((size_t)1 << 20)) {
		size_t ks = inSize >> 10;
		snprintf(outStr, MAX_OUT_LEN, "%lluK", (unsigned long long)ks);
	}
	else if(inSize < ((size_t)1 << 30)) {
		size_t megs = inSize >> 20;
		snprintf(outStr, MAX_OUT_LEN, "%lluM", (unsigned long long)megs);
	}
    else {
		size_t gigs = inSize >> 30;
		snprintf(outStr, MAX_OUT_LEN, "%lluG", (unsigned long long)gigs);
    }
	return outStr;
}

/* 
 * Obtain compact string representation as power of 2 if possible, else
 * as an unsigned long long. Caller must free() the result. 
 */
char *fftStringRepPow2(
	size_t inSize)
{
	char *outStr = (char *)malloc(MAX_OUT_LEN);
	if(outStr == NULL) {
		return strdup("Malloc error");
	}
	unsigned log2n;
	if(fftIsPowerOfTwo(inSize, &log2n)) {
		sprintf(outStr, " 2^%u ", log2n);
	}
	else {
		sprintf(outStr, "%llu", (unsigned long long)inSize);
	}
	return outStr;
}

/* 
 * Parse an input string, possibly with a 'K', 'k', 'M', or 'm' suffix, 
 * or in 2^n form, into a size_t.
 */
size_t fftParseStringRep(
	const char *inStr)
{
	unsigned shiftCnt = 0;
	const char *endP = NULL;
	
	endP = strstr(inStr, "2^");
	if(endP) {
		size_t expo = strtol(endP + 2, NULL, 0);
		return (size_t)1 << expo;
	}
	
	endP = strchr(inStr, 'm');
	if(endP != NULL) {
		shiftCnt = 20;
	}
	if(endP == NULL) {
		endP = strchr(inStr, 'M');
		if(endP != NULL) {
			shiftCnt = 20;
		}
	}
	if(endP == NULL) {
		endP = strchr(inStr, 'k');
		if(endP != NULL) {
			shiftCnt = 10;
		}
	}
	if(endP == NULL) {
		endP = strchr(inStr, 'K');
		if(endP != NULL) {
			shiftCnt = 10;
		}
	}
	if(endP != NULL) {
		size_t len = endP - inStr;
		char num[len + 1];
		memmove(num, inStr, len);
		num[len] = '\0';
		size_t rtn = strtol(num, NULL, 0);
		return (rtn << shiftCnt);
	}
	else {
		return strtol(inStr, NULL, 0);
	}
}

/* 
 * Round up size to next power of 2.
 *
 * Returns 1, 2, 4, 8....
 * Returns 0, 1, 2, 3... in *log2NumSamples
 */
size_t fftRoundNumSamples(
	size_t numSamples,
	unsigned *log2NumSamples)		/* RETURNED */
{
	size_t pow2NumSamples = 1;		/* 1, 2, 4, 8... */
	*log2NumSamples = 0;			/* 0, 1, 2, 3... */
	
	while(numSamples > pow2NumSamples) {
		pow2NumSamples <<= 1;
		(*log2NumSamples)++;
	}
	return pow2NumSamples;
}

#pragma mark --- Misc. functions ---

/* 
 * Determine if a number is a power of two; if so, return true and  
 * (optionally) the exponent in exp. Else return false and return the
 * largest power of 2 less than num.
 */
bool fftIsPowerOfTwo(
	size_t num,
	unsigned *exp)		/* optionally RETURNED */
{
	unsigned _exp = 0;
	bool ourRtn = false;
	
	unsigned maxBits = sizeof(size_t) * 8;
	size_t pwrOfTwo = 1;
	
	for(_exp=0; _exp<maxBits; _exp++) {
		if(num == pwrOfTwo) {
			ourRtn = true;
			break;
		}
		if(num < pwrOfTwo) {
			/* Not power of two, terminate */
			if(exp != NULL) {
				*exp = _exp - 1;
			}
			return false;
		}
		pwrOfTwo <<= 1;
	}
	if(ourRtn && (exp != NULL)) {
		*exp = _exp;
	}
	return ourRtn;
}

/* 
 * Determine the number of active virtual CPU cores. 
 * On Xeon and prior Intel processors there is one virtual core
 * per physical core. Nehalem (X5570, W3520) has two
 * virtual cores per physical core. 
 */
unsigned numCpuCores()
{
	int     count;
	size_t  size = sizeof(count);

	if (sysctlbyname("hw.activecpu", &count, &size, NULL, 0)) {
		perror("sysctlbyname(hw.activecpu)");
		printf("***Error executing sysctlbyname(hw.activecpu)\n");
		return 1;
	}
	if(count <= 0) {
		return 1;
	}
	else {
		return (unsigned)count;
	}
}

/*
 * Determine L2+L3+... cache size per CPU core, i.e. total cache size - other than
 * L1 - per core.
 *
 * These sysctl ops are partially documented here:
 *
 * http://developer.apple.com/ReleaseNotes/Performance/RN-AffinityAPI/index.html
 *
 * hw.cacheconfig[0] reports the total number of logical processors.
 * hw.cacheconfig[1] reports the L1 sharing (the number of logical processors 
 *		sharing a level 1 cache).
 * hw.cacheconfig[2] reports the L2 sharing.
 * hw.cacheconfig[N] reports the L_N sharing; this array terminates with a zero.
 *      I.e. for a machine with no L3 cache, hw.cacheconfig[3] is 0.
 * 
 * hw.cachesize[0] reports the size of memory.
 * hw.cachesize[1] reports the size of the L1 data cache.
 * hw.cachesize[2] reports the size of the L2 cache.
 * hw.cachesize[N] reports the size of the L_N cache; this array is zero-terminated.
 *
 * The L2 cache size is not the total L2 cache on a processor die; it's the amount
 * of L2 cache shared by hw.cacheconfig[2] cores. So for example on a MacPro Xeon X5482,
 * a dual quad-core machine, each processor die has a total of 12 MB of L2 cache, 
 * but sysctl reports:
 *
 *   Processors per L2 cache : 2
 *   L2 size : 6291456
 *
 * That is, L2 is divided into 2 logical units, each of which is shared by 2
 * cores. We return 3145728 for that machine. 
 *
 * Returns 0 on error. 
 */
 
/*
 * As of 14 August 2009, on the only machine with an L3 cache - an 8-core Nehalem - 
 * it's not advantageous to use the L3 cache here. Future work might use it,
 * with a fair amount of tuning elsewhere. For now we just report the L2 cache 
 * size.
 */
#define FFT_CACHE_PER_CORE_USE_L3       0

uint64_t l2CachePerCore()
{
	size_t len = SYSCTL_BUF_SIZE;
	char configBuf[SYSCTL_BUF_SIZE];        /* hw.cacheconfig */
	char sizeBuf[SYSCTL_BUF_SIZE];          /* hw.cachesize */
    
	if(sysctlbyname("hw.cacheconfig", configBuf, &len, NULL, 0)) {
		dbprintf("sysctlbyname(hw.cacheconfig) error");
		return 0;
	}
    unsigned numEntries = len / sizeof(uint64_t);
	if(numEntries < 3) {
		dbprintf("sysctlbyname(hw.cacheconfig): short return");
		return 0;
	}
	uint64_t *procPerCache = (uint64_t *)configBuf;
	
	len = SYSCTL_BUF_SIZE;
	if(sysctlbyname("hw.cachesize", sizeBuf, &len, NULL, 0)) {
		dbprintf("sysctlbyname(hw.cachesize) error");
		return 0;
	}
    numEntries = len / sizeof(uint64_t);
	if(numEntries < 3) {
		dbprintf("sysctlbyname(hw.cachesize): short return");
		return 0;
	}
	uint64_t *cacheSize = (uint64_t *)sizeBuf;
	
    uint64_t cachePerCore = 0;
    
    unsigned maxEntries = numEntries;
    #if     !FFT_CACHE_PER_CORE_USE_L3
    maxEntries = 3;
    #endif
    
    /* add them up, one entry per cache level */
    for(unsigned dex=2; dex<maxEntries; dex++) {
        uint64_t procs = procPerCache[dex];
        uint64_t thisCacheSize = cacheSize[dex];
        if((procs == 0) || (thisCacheSize == 0)) {
            /* normal termination */
            break;
        }
        cachePerCore += (thisCacheSize / procs);
    }
    return cachePerCore;
}

/*
 * Flush all cache associated with specified memory range.
 * The first cut of this assumes that clflush() operates on all caches
 * for all CPUs int he system. Docs say:
 *  
 *   The invalidation is broadcast throughout the cache coherence 
 *   domain.
 *
 * But I can't find a good definition of a cache coherence domain - is
 * it all the cores in a system,. or just the cores that share one cache?
 */
void fftFlushCache(
    void *addr,
    size_t numBytes)
{
    #if     (FFT_INTEL && 1)

    unsigned char *cp = (unsigned char *)addr;
    for(size_t dex=0; dex<numBytes; dex+=CACHE_LINE_SIZE) {
        _mm_clflush(cp);
        cp += CACHE_LINE_SIZE;
    }
    _mm_mfence();
    
    #endif  /* FFT_INTEL */
}

/*
 * Determine number of iterations to measure, given signal size.
 */
#define ONEMEG  (1024 * 1024)

unsigned fftIterationsForSize(
    size_t numComplex)
{
    unsigned numIter = 4;
    if(numComplex > (8 * ONEMEG)) {
        numIter = 1;
    }
    else if(numComplex > (2 * ONEMEG)) {
        numIter = 2;
    }
    return numIter;
}

