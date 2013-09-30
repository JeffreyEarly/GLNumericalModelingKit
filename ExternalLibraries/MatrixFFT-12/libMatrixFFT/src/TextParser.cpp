/*	File: TextParser.cpp 
	
	Description:
		General purpose text parser class.
	
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
 * TextParser.cpp - general purpose text parser class
 *
 * Created 4 Jan 2006.
 * Copyright 2008 by Apple, Inc. 
 */

#include "TextParser.h"
#include "fileIo.h"
#include <strings.h>
#include <stdlib.h>
#include <stdexcept>
#include <ctype.h>

#define TP_DEBUG		0

using namespace std;

TextParser::TextParser(const char *_fileName, bool exitOnErr)
	: mBof(NULL), 
	  mEof(NULL), 
	  mCursor(NULL),
	  mLineLength(LINE_LENGTH_MAX),
	  mWeMallocd(false),
	  mFileName(NULL)
{
	unsigned dataLen;
	unsigned char *data;
	
	if(readFile(_fileName, &data, &dataLen)) {
		fprintf(stderr, "***TextParser: error reading %s\n", _fileName);
		if(exitOnErr) {
			exit(1);
		}
		throw invalid_argument("Bad file");
	}
	mBof = mCursor = (char *)data;
	mEof = mBof + dataLen;
	mWeMallocd = true;
	mFileName = strdup(_fileName);
}

TextParser::TextParser(const char *data, size_t dataLen)
	: mBof((char *)data), 
	  mEof((char *)(data + dataLen)),
	  mCursor((char *)data),
	  mLineLength(LINE_LENGTH_MAX),
	  mWeMallocd(false),
	  mFileName(NULL)
{

}

TextParser::~TextParser()
{
	if(mWeMallocd && (mBof != NULL)) {
		free(mBof);
	}
	if(mFileName) {
		free(mFileName);
	}
}
	
/*
 * Set cursor position, typically to either 0 or a value obtained
 * from getCursor().
 */
void TextParser::setCursor(size_t cursor)
{
	if(cursor == 0) {
		mCursor = mBof;
		return;
	}
	
	/* detect illegal or improper values */
	char *tentative = mBof + cursor;
	if(tentative > mEof) {
		printf("***TextParser::setCursor after EOF\n");
		throw invalid_argument("Bad cursor");
	}
	switch(tentative[-1]) {
		/* not at BOF, this should be newline */
		case '\n':
		case '\r':
			break;
		default:
			printf("***TextParser::setCursor not at beginning of line\n");
			/* that's enough; proceed */
			break;
	}
	mCursor = tentative;
}

/* 
 * Copy next line, up to but not including next newline or EOF, into
 * caller-supplied buffer. Caller's buffer is NULL-terminated. 
 * Cursor is left after trailing newline or EOF. 
 * Returns nonzero if cursor is at EOF. 
 */
int TextParser::getLine(char *lineBuf)
{
	if(mCursor == mEof) {
		/* 
		 * The only time we return 'eof' - subsequently we return a NULL-terminated
		 * string even if the file does not end in newline.
		 */
		return -1;
	}
	
	char *outp = lineBuf;
	
	while(mCursor < mEof) {
		char c = *mCursor++;
		int eol = false;
		switch(c) {
			/* 
			 * Note we only skip this line's newline terminator, not a sequence 
			 * of them, to allow caller to detect empty lines.
			 */
			case '\n':
				/* normal UNIX EOL */
				eol = true;
				break;
			case '\r':
				/* skip over Windows cr/lf */
				if((mCursor < mEof) && (*mCursor == '\n')) {
					mCursor++;
				}
				eol = true;
				break;
			default:
				break;
		}
		if(eol) {
			break;
		}
		*outp++ = c;
	}
	*outp = '\0';
	
	#if		TP_DEBUG
	printf("== found line '%s'\n", lineBuf);
	#endif	/* TP_DEBUG */
	
	return 0;
}

/*
 * Parse a line - as a NULL-terminated string - into whitespace-separated
 * tokens. Free the results with freeTokens(). Returns number of tokens.
 */
unsigned TextParser::parseLine(const char *lineBuf, const char **&tokens)
{
	unsigned numTokens = 0;
	char **tok = NULL;
	
	const char *cp = lineBuf;
	
	for(;;) {
		/* skip whitespace - cp is after last token or at beginning of line */
		for(;;) {
			char c = *cp;
			bool foundToken = false;
			switch(c) {
				case '\0':
				case ' ':
				case '\t':
					cp++;
					break;
				default:
					foundToken = true;
					break;
			}
			if(foundToken) {
				break;
			}
			if(c == '\0') {
				/* we're done */
				tokens = (const char **)tok;
				return numTokens;
			}
		}
		
		/* cp is at start of token */
		const char *end = cp + 1;
		for(end=cp+1; *end!='\0'; end++) {
			if(isspace(*end)) {
				break;
			}
		}
		
		/* end is one past end of token */
		unsigned tokenLen = end - cp;
		char *token = (char *)malloc(tokenLen + 1);
		memmove(token, cp, tokenLen);
		token[tokenLen] = '\0';
		numTokens++;
		tok = (char **)realloc(tok, (numTokens * sizeof(char *)));
		tok[numTokens-1] = token;
		
		#if		TP_DEBUG
		printf("== found token '%s'\n", token);
		#endif	/* TP_DEBUG */
		
		if(*end == '\0') {
			/* we're done */
			tokens = (const char **)tok;
			return numTokens;
		}
		
		/* prepare for next token search */
		cp = end;
	}
	/* not reached */
	return 0;
}

/*
 * Convenience function: getLine() followed by parseLine().
 * Returns nonzero if cursor is at EOF.
 */
int TextParser::getTokens(unsigned &numTokens, const char **&tokens)
{
	char lineBuf[mLineLength];
	
	int rtn = getLine(lineBuf);
	if(rtn) {
		tokens = NULL;
		numTokens = 0;
		return -1;
	}
	numTokens = parseLine(lineBuf, tokens);
	return 0;
}


/* dispose of tokens obtained from parseLine() or getTokens() */
void TextParser::freeTokens(unsigned numTokens, const char **tokens)
{
	unsigned dex;
	
	for(dex=0; dex<numTokens; dex++) {
		free((void *)(tokens[dex]));
	}
	free((void *)tokens);
}

/* 
 * Advance cursor past current line.
 * Returns nonzero if cursor at EOF.
 */
int TextParser::skipLine()
{
	char lineBuf[mLineLength];
	
	return getLine(lineBuf);
}

/*
 * Find line containing specified string. Returns the line in lineBuf
 * and returns true if found, else returns false. Cursor is positioned
 * after found line on success, else at EOF. 
 */
bool TextParser::findLine(
	const char *str,
	char *lineBuf)
{
	while(1) {
		if(getLine(lineBuf)) {
			/* EOF */
			return false;
		}
		if(strstr(lineBuf, str)) {
			return true;
		}	
	}
}

