//
// Copyright (c) 2005-2012, Matthijs van Leeuwen, Jilles Vreeken, Koen Smets
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
// 
// Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials // provided with the distribution.
// Neither the name of the Universiteit Utrecht, Universiteit Antwerpen, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
#ifndef __BINARYFICISCFILEASYNC_H
#define __BINARYFICISCFILEASYNC_H

#include "BinaryFicIscFile.h"

#define BUF_SIZE 4096


class BASS_API BinaryFicIscFileAsync : public BinaryFicIscFile {
public:
	BinaryFicIscFileAsync();
	virtual ~BinaryFicIscFileAsync();

protected:
	virtual void	OpenFile(const string &filename, const FileOpenType openType = FileReadable);

	virtual size_t	ReadLine(char* buffer, int size);

	virtual void	ReadBuf(char* buffer, uint32 size);

	char			GetChar();

private:
	void			RefillBuffer();

#if defined (_WINDOWS)
	HANDLE mHFile;
	OVERLAPPED mOverlapped;
	HANDLE mHEvent;
#endif
	char mFileBuf[2][BUF_SIZE];
	int mCurrentBuf;
	int mBytesRead;
	char* mPointer;
	int mCount;
	bool isOverlapped;
	bool eof;

};

#endif // __BINARYFICISCFILEASYNC_H
