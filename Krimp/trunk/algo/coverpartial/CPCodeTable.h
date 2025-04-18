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
#ifndef __CPCODETABLE_H
#define __CPCODETABLE_H

class CTISList;
class CoverSet;

#include "../CodeTable.h"

class CPCodeTable : public CodeTable {
public:
	CPCodeTable();
	CPCodeTable(const CPCodeTable& ct);
	virtual ~CPCodeTable();
	virtual CPCodeTable* Clone() { return new CPCodeTable(*this); }

	virtual string		GetShortName() { return "cp"; }
	static string		GetConfigBaseName() { return "coverpartial"; }

	virtual bool		NeedsDBOccurrences() { return true; }
	virtual void		UseThisStuff(Database *db, ItemSetType type, CTInitType ctinit, uint32 maxCTElemLength=0, uint32 toMinSup = 0);

	virtual void		CoverDB(CoverStats &stats);
	virtual void		CoverMaskedDB(CoverStats &stats, uint32 *rowList, uint32 numRows);
	virtual double		CalcTransactionCodeLength(ItemSet *transaction);

	virtual void		Add(ItemSet *is, uint64 candidateId = UINT64_MAX_VALUE);
	virtual void		Del(ItemSet *is, bool zap, bool keepList);

	virtual void		CommitAdd(bool updateLog=true);
	virtual void		CommitLastDel(bool zap,bool updateLog=true);
	virtual void		CommitAllDels(bool updateLog=true);
	virtual void		RollbackAdd();
	virtual void		RollbackLastDel();
	virtual void		RollbackAllDels();
	virtual void		RollbackCounts();

	virtual void		AddOneToEachUsageCount();

	virtual islist*		GetItemSetList();
	virtual islist*		GetPrePruneList();
	virtual islist*		GetPostPruneList();
	virtual islist*		GetSanitizePruneList();
	virtual void		UpdatePostPruneList(islist *pruneList, const islist::iterator &after);

	virtual void		SetAlphabetCount(uint32 item, uint32 count);
	virtual uint32		GetAlphabetCount(uint32 item);

	// Write state-change to log
	virtual void		UpdateCTLog();

protected:
	virtual void		CoverNoBinsDB(CoverStats &stats);
	virtual void		CoverBinnedDB(CoverStats &stats);

	bool				mIsBinnedDB;

	CoverSet			*mCoverSet;
	CoverSet			**mCDB;

	uint32				mCTLen;
	CTISList			**mCT;
	uint32				mMaxIdx;

	uint32				mNumDBRows;

	uint32*				mAlphabet;
	uint32*				mOldAlphabet;
	uint32				mABLen;
	uint32				mABUsedLen;
	uint32				mABUnusedCount;

	uint32*				mAllRows;

	DelInfo				mLastDel;
	dellist				mDeleted;

	uint32				mRollback;
	bool				mRollbackAlph;

	ChangeType			mChangeType;
	ItemSet*			mChangeItem;
	islist::iterator	mChangeIter;
	uint32				mChangeRow;

};

#endif // __CPCODETABLE_H
