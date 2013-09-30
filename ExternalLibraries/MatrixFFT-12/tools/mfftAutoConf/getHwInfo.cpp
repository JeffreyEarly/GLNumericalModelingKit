/*	File: getHwInfo.cpp
	
	Description:
		Obtain info about current CPU and memory.
	
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
 
#include "getHwInfo.h"
#include <sys/sysctl.h>
#include <strings.h>
#include <stdio.h>

#define CPU_BRAND_SYSCTL_NAME       "machdep.cpu.brand_string"


/*
 * Obtain CPU brand string as a NULL-terminated string, placed in caller's buffer
 */
static void cpuBrandString(
    char *cpuString)            /* at least CPU_STRING_LEN bytes */
{
	size_t len = CPU_STRING_LEN;
	char sysctlBuf[CPU_STRING_LEN];
	
	if(sysctlbyname(CPU_BRAND_SYSCTL_NAME, sysctlBuf, &len, NULL, 0)) {
		strcpy(cpuString, "<unknown>");
		return;
	}
	
	/* 
	 * output of that sysctl looks like this:
	 *
	 * Genuine Intel(R) CPU            1500  @ 2.00GHz
     *
     * Remove redundant spaces.
	 */
	char *inStr = sysctlBuf;
    char *outStr = cpuString;
	len = strlen(inStr);

	bool lastWasSpace = false;
	for(unsigned dex=0; dex<len; dex++) {
		char c = *inStr++;
		if(c == ' ') {
			/* compress multiple spaces into one */
			if(!lastWasSpace) {
				*outStr++ = c;
			}
			lastWasSpace = true;
		}
		else {
			lastWasSpace = false;
            *outStr++ = c;
		}
	}
    *outStr = '\0';
}

#define SYSCTL_BUF_SIZE 1024

int getHwInfo(HWInfo *hwInfo)
{
    memset(hwInfo, 0, sizeof(*hwInfo));
    
    /* 
     * Something like:
     * Genuine Intel(R) CPU 1500  @ 2.00GHz
     */
    cpuBrandString(hwInfo->cpuInfo);
    
    /*
     * hw.cacheconfig sysctl: 0-terminated array of 64-bit integers
     * val[0] = # virtual processors
     * val[n] = virtual processors per L[n] cache
     */
    char sysctlBuf[SYSCTL_BUF_SIZE];
	size_t len = SYSCTL_BUF_SIZE;
	if(sysctlbyname("hw.cacheconfig", sysctlBuf, &len, NULL, 0)) {
		perror("sysctlbyname");
		printf("sysctlbyname(%s) failed\n", "hw.cacheconfig");
		return -1;
	}
    
    uint64_t *up = (uint64_t *)sysctlBuf;
    unsigned numEntries = len / sizeof(uint64_t);
    if(numEntries < 2) {
        fprintf(stderr, "getHwInfo: short sysctl(hw.cacheconfig)");
        return -1;
    }
    
    /* skip CPU count */
    up++;
    numEntries--;
    if(numEntries > MAX_NUM_CACHES) {
        numEntries = MAX_NUM_CACHES;
    }
    for(unsigned dex=0; dex<numEntries; dex++) {
        hwInfo->cacheInfo[dex].processorsPerCache = *up++;
    }
    
    /*
     * hw.cachesize sysctl: 0-terminated array of 64-bit integers
     * val[0] = RAM size
     * val[n] = L[n] cache size
     */
	len = SYSCTL_BUF_SIZE;
	if(sysctlbyname("hw.cachesize", sysctlBuf, &len, NULL, 0)) {
		perror("sysctlbyname");
		printf("sysctlbyname(%s) failed\n", "hw.cacheconfig");
		return -1;
	}
    
    up = (uint64_t *)sysctlBuf;
    numEntries = len / sizeof(uint64_t);
    if(numEntries < 2) {
        fprintf(stderr, "getHwInfo: short sysctl(hw.cachesize)");
        return -1;
    }
    hwInfo->RAMSize = *up++;
    numEntries--;
    if(numEntries > MAX_NUM_CACHES) {
        numEntries = MAX_NUM_CACHES;
    }
    for(unsigned dex=1; dex<=numEntries; dex++) {
        hwInfo->cacheInfo[dex-1].cacheSize = *up++;
    }
    
    /*
     * hw.activecpu sysctl: returns one 32-bit integer indicating the number of 
     * virtual CPUs actually currently available.
     * Some processors, e.g. Nehalem, report multiple virtual CPUs per physical core.
     * Others, e.g. Xeon, have a 1-to-1 mapping between virtual and physical cores.
     */
	int32_t activeCpus;
	size_t  size = sizeof(activeCpus);

	if (sysctlbyname("hw.activecpu", &activeCpus, &size, NULL, 0)) {
		perror("sysctlbyname");
		printf("***Error executing sysctlbyname(hw.activecpu)\n");
		return 1;
	}
	if(activeCpus <= 0) {
        /* punt around this error */
		activeCpus = 1;
	}
    hwInfo->numVirtualCPUs = (uint32_t)activeCpus;
    
    /*
     * hw.physicalcpu - one 32-bit integer indicating the number of physical CPUs
     * currently configured
     */
	int32_t physCpus;
	size = sizeof(physCpus);

	if (sysctlbyname("hw.physicalcpu", &physCpus, &size, NULL, 0)) {
		perror("sysctlbyname");
		printf("***Error executing sysctlbyname(hw.physicalcpu)\n");
		return 1;
	}
	if(physCpus <= 0) {
        /* punt around this error */
		physCpus = 1;
	}
    hwInfo->numPhysicalCPUs = (uint32_t)physCpus;

    return 0;
}


