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
#include "../clobal.h"

#include <stdlib.h>
#include <stdio.h>
#include "hash_map.h"

using namespace Afopt;

int CHash_Map::hashfunc(int nkey)
{
	return nkey%mnhashtable_size;
}

void CHash_Map::init(int size)
{
	int i;

	mnhashtable_size = size;
	mpentry_list = new HASH_NODE[mnhashtable_size];
	for(i=0;i<mnhashtable_size;i++)
	{
		mpentry_list[i].num_of_nodes = 0;
		mpentry_list[i].next = NULL;
	}

	mnmemory_size = mnhashtable_size*sizeof(HASH_NODE);
}

HASH_ENTRY* CHash_Map::find(int nkey)
{
	int nhash_value;
	HASH_ENTRY* phash_entry;

	nhash_value	= hashfunc(nkey);

	phash_entry = mpentry_list[nhash_value].next;
	while(phash_entry!=NULL)
	{
		if(phash_entry->key==nkey)
			break;
		phash_entry = phash_entry->next;
	}

	return phash_entry;
}

void CHash_Map::delete_entry(int nkey)
{
	int nhash_value;
	HASH_ENTRY* phash_entry, *pprev_entry;

	nhash_value	= hashfunc(nkey);

	pprev_entry = NULL;
	phash_entry = mpentry_list[nhash_value].next;
	while(phash_entry!=NULL)
	{
		if(phash_entry->key==nkey)
		{
			if(pprev_entry==NULL)
			{
				mpentry_list[nhash_value].next = phash_entry->next;
				delete phash_entry;
				phash_entry = mpentry_list[nhash_value].next;
			}
			else 
			{
				pprev_entry->next = phash_entry->next;
				delete phash_entry;
				phash_entry = pprev_entry->next;
			}
			
			mnmemory_size -= sizeof(HASH_ENTRY);
			break;
		}
		else 
		{
			pprev_entry = phash_entry;
			phash_entry = phash_entry->next;
		}
	}
}


HASH_ENTRY* CHash_Map::insert(int nkey, void* nvalue)
{
	int nhash_value;
	HASH_ENTRY *pnewentry;

	nhash_value = hashfunc(nkey);

	pnewentry = new HASH_ENTRY;
	pnewentry->key = nkey;
	pnewentry->pvalue = nvalue;
	pnewentry->next = mpentry_list[nhash_value].next;
	mpentry_list[nhash_value].next = pnewentry;
	mpentry_list[nhash_value].num_of_nodes++;

	mnmemory_size += sizeof(HASH_ENTRY);

	return pnewentry;
}

HASH_ENTRY* CHash_Map::insert(int nkey, int nvalue)
{
	int nhash_value;
	HASH_ENTRY *pnewentry;

	nhash_value = hashfunc(nkey);

	pnewentry = new HASH_ENTRY;
	pnewentry->key = nkey;
	pnewentry->value = nvalue;
	pnewentry->next = mpentry_list[nhash_value].next;
	mpentry_list[nhash_value].next = pnewentry;
	mpentry_list[nhash_value].num_of_nodes++;

	mnmemory_size += sizeof(HASH_ENTRY);

	return pnewentry;
}


void CHash_Map::destroy()
{
	int i;
	HASH_ENTRY* phash_entry, *pnext_entry;

	for(i=0;i<mnhashtable_size;i++)
	{
		phash_entry = mpentry_list[i].next;
		while(phash_entry!=NULL)
		{
			pnext_entry = phash_entry->next;
			delete phash_entry;
			mnmemory_size -= sizeof(HASH_ENTRY);
			phash_entry = pnext_entry;
		}
		mpentry_list[i].next = NULL;
	}
	delete []mpentry_list;
	mnmemory_size -= mnhashtable_size*sizeof(HASH_NODE);
	if(mnmemory_size!=0)
		printf("Error with the size of hash map memory\n");
	mnhashtable_size = 0;
}

/*
int CHash_Map::get_num_of_large(double threshold)
{
	int i, num_of_large_ones;
	HASH_ENTRY* phash_entry;

	num_of_large_ones = 0;
	for(i=0;i<mnhashtable_size;i++)
	{
		phash_entry = mpentry_list[i].next;
		while(phash_entry!=NULL)
		{
			if(phash_entry->value>=threshold)
				num_of_large_ones++;
			phash_entry = phash_entry->next;
		}
	}

	return num_of_large_ones; 
}

void CHash_Map::get_large_entries(double threshold, ITEM_COUNTER *plargeitems)
{
	int i, k;
	HASH_ENTRY* phash_entry;

	k = 0;
	for(i=0;i<mnhashtable_size;i++)
	{
		phash_entry = mpentry_list[i].next;
		while(phash_entry!=NULL)
		{
			if(phash_entry->value>=threshold)
			{
				plargeitems[k].item = phash_entry->key;
				plargeitems[k].support = phash_entry->value;
				k++;
			}
			phash_entry = phash_entry->next;
		}
	}
}

*/