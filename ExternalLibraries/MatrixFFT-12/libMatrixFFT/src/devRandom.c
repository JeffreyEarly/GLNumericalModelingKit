/*	File: devRandom.c  
	
	Description:
		Simplified interface to /dev/random.
	
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
 *  devRandom.c - simplified interface to /dev/random
 *
 *  Created on 2/3/05.
 *	Copyright 2005 by Apple Computer. 
 */

#include <libMatrixFFT/devRandom.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <fcntl.h>
#include <limits.h>

/* 
 * Obtain random bytes. Returns nonzero on error. Caller allocates randBytes buf.
 * OS X happens to have an excellent source of random data in /dev/random.
 */
static int randFd = 0;
int getRandomBytes(unsigned char *randBytes, unsigned randBytesLen)
{
	int brtn = 0;
	ssize_t bytesRead;	
	
	if(randFd == 0) {
		randFd = open("/dev/random", O_RDONLY);
	}
	bytesRead = read(randFd, randBytes, (size_t)randBytesLen);
	if(bytesRead < 0) {
		perror("read(/dev/random)");
		brtn = -1;
	}
	else if(bytesRead != (ssize_t)randBytesLen) {
		fprintf(stderr, "getRandomBytes: short read\n");
		brtn = -1;
	}
	return brtn;
}

/* min <= return <= max */
unsigned genRandInt(unsigned min, unsigned max)
{
	unsigned i;
	if(min == max) {
		return min;
	}
	getRandomBytes((unsigned char *)&i, (unsigned)sizeof(unsigned));
	return (min + (i % (max - min + 1)));
}	

/* return random signed integer of native size */
int getRandomInt(void)
{
	int i;
	getRandomBytes((unsigned char *)&i, (unsigned)sizeof(int));
	return i;
}

/* return random double, -1.0 <= return <= 1.0 */
double getRandomDouble(void)
{
	double i = (double)getRandomInt();
	if(i > 0) {
		return i / (double)INT_MAX;
	}
	else return -(i / (double)INT_MIN);
}

