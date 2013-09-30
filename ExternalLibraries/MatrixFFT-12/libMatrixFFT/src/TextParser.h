/*	File: TextParser.h  
	
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
 * TextParser.h - general purpose text parser class
 *
 * Created 4 Jan 2006.
 * Copyright 2008 by Apple, Inc. 
 */
 
#include <sys/param.h>

#define LINE_LENGTH_MAX		512			/* can be overridden via maxLineLength() */

/*
 * TextParser operates on a block of text obtained from a file or directly
 * from the caller. It processes text one line at a time, maintaining a cursor
 * which at any given time is at one of three places:
 *
 * -- beginning of file
 * -- after a newline (which could be /n, /r/n, or /r)
 * -- end of file
 *
 * Caller obtains the "next line" of text via getLine(), which copies the
 * line of text (sans newline) into a caller-supplied buffer. The cursor is
 * at the next character after the newline on return. 
 * 
 * The parseLine() function parses one line of text into separate mallocd tokens 
 * which are separated in the line of input text by spaces (' ' or '\t').
 *
 * Alternatively, getTokens() combines the previous two operations, returning
 * all the tokens in the "next line" of text and advancing the cursor to the 
 * next character after the newline. 
 */
class TextParser
{
private:
	/* not supported */
	TextParser(TextParser &src);
	void operator = (const TextParser &);

public:
	TextParser(const char *fileName, bool exitOnErr = true);
	TextParser(const char *data, size_t dataLen);
	~TextParser();
	
	void		maxLineLength(size_t lineLength)	{ mLineLength = lineLength; }
	size_t		maxLineLength()						{ return mLineLength; }
	
	/*
	 * Set cursor position, typically to either 0 or a value obtained
	 * from getCursor(). SetCursor() is the only function after our constructor
	 * that can throw. 
	 */
	void		setCursor(size_t cursor);
	size_t		getCursor()							{ return mCursor - mBof; }
	
	/* 
	 * Copy next line, up to but not including next newline or EOF, into
	 * caller-supplied buffer. Caller's buffer is NULL-terminated. 
	 * Cursor is left after trailing newline or at EOF. 
	 * Returns nonzero if cursor is at EOF on entry.  
	 * Multiple newlines are not coalesced; thus this can return an 
	 * empty (though properly NULL-terminated) string. 
	 */
	int			getLine(char *lineBuf);
	
	/*
	 * Parse a line - as a NULL-terminated string, typically obtained from
	 * getLine() - into tokens which are whitespace-separated in the input line. 
	 * Free the results with freeTokens(). Returns number of tokens.
	 */
	unsigned	parseLine(const char *lineBuf, const char **&tokens);
	
	/*
	 * Convenience function: getLine() followed by parseLine().
	 * Returns nonzero if cursor is at EOF.
	 */
	int			getTokens(unsigned &numTokens, const char **&tokens);
	
	/* dispose of tokens obtained from parseLine() or getTokens() */
	void		freeTokens(unsigned numTokens, const char **tokens);
	
	/* 
	 * Advance cursor past current line.
	 * Returns nonzero if cursor at EOF.
	 */
	int			skipLine();
	
	/* obtain filename used in constructor, if we have one */
	const char	*fileName()			{ return mFileName; }
	
	/*
	 * Find line containing specified string, starting at current cursor.
	 * Returns the line in lineBuf and returns true if found, else returns 
	 * false. Cursor is positioned after found line on success, else at EOF. 
	 */
	bool findLine(const char *str, char *lineBuf);
	
private:
	char		*mBof;				// beginning of data
	char		*mEof;				// one char past end of data
	char		*mCursor;
	unsigned	mLineLength;
	bool		mWeMallocd;			// we have to free mBof
	char		*mFileName;			// strdup'd, if we have it 
};
