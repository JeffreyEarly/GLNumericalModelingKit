/*	File: fftwUtils.cpp
	
	Description:
		Common utility functions for FFTW tests.
	
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
 * fftwUtils.cpp - common utility functions for FFTW tests.
 *
 * Created 01/08/2009. 
 * Copyright 2009 by Apple, Inc. 
 */

#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <errno.h>
#include <stdlib.h>
#include "fftwUtils.h"

/* 
 * Locations of accumulate wisdom files.
 */
#define FFTW_WISDOM_FILE_F	"fftwWisdomFloat"
#define FFTW_WISDOM_FILE_D	"fftwWisdomDouble"
#define FFTW_WISDOM_DIR		"/tmp"
#
#ifdef	DEBUG
#define dprintf(args...)			printf(args)
#else
#define dprintf(args...)
#endif

static void fftwWisdomFilename(
	bool doublePrec,
	char *pathName)		// allocated by caller
{
	const char *wisdomDir = getenv("LOCAL_BUILD_DIR");
	if((wisdomDir == NULL) || (wisdomDir[0] != '/')) {
		wisdomDir = FFTW_WISDOM_DIR;
	}
	sprintf(pathName, "%s/%s", wisdomDir,
		doublePrec ? FFTW_WISDOM_FILE_D : FFTW_WISDOM_FILE_F);
}

/* 
 * Obtain and import wisdom.
 * Returns nonzero on error (which should generally not be a big deal; our wisdom files
 * are in /tmp/, so they get nuked on each boot).
 * Note test programs normally do *not* call this or fftwSaveWisdom if they are in 
 * ESTIMATE mode; that is used to show the performance degradation by failing to use
 * "good wisdom" - and we dfinitely don't want to save "bad wisdom" at the end of 
 * a test that's run in ESTIMATE mode.  
 */
int fftwGetWisdom(bool doublePrec)
{
	char fileName[1024];
	fftwWisdomFilename(doublePrec, fileName);
	
	FILE *inFile = fopen(fileName, "r");
	if(inFile == NULL) {
		dprintf("fftwGetWisdom(): no wisdom file found\n");
		return -1;
	}
	int ourRtn = 0;
	if(doublePrec) {
		ourRtn = fftw_import_wisdom_from_file(inFile);
	}
	else {
		ourRtn = fftwf_import_wisdom_from_file(inFile);
	}
	if(ourRtn) {
		/* success */
		ourRtn = 0;
	}
	else {
		ourRtn = 1;
		dprintf("fftwGetWisdom(): fftw*_import_wisdom_from_file error\n");
	}
	fclose(inFile);
	return ourRtn;
}

/*
 * Save current wisdom to file.
 */
int fftwSaveWisdom(bool doublePrec) 
{
	char fileName[1024];
	fftwWisdomFilename(doublePrec, fileName);
	
	FILE *outFile = fopen(fileName, "w");
	if(outFile == NULL) {
		printf("***fftwSaveWisdom: error creating wisdom file\n");
		return -1;
	}
	if(doublePrec) {	
		fftw_export_wisdom_to_file(outFile);
	}
	else {
		fftwf_export_wisdom_to_file(outFile);
	}
	fclose(outFile);
	return 0;
}
