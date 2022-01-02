// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Note that this file originally did not compile.  Me (holtgrew) change it
// it it compiled.  However, it might not make sense to fix it so it actually
// works.
// ==========================================================================

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#ifndef SEQAN_HEADER_MISC_MAP_H
#define SEQAN_HEADER_MISC_MAP_H

#include <algorithm>
#include "misc_base.h"


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

	//////////////////////////////////////////////////////////////////////////////
	//

	template <typename TPair>
	inline typename TPair::T1 & keyOf(TPair &pair) {
		return getValueI1(pair);
	}
	template <typename TPair>
	inline typename TPair::T1 const & keyOf(TPair const &pair) {
		return getValueI1(pair);
	}
	template <typename TPair>
	inline typename TPair::T2 & objectOf(TPair &pair) {
		return getValueI2(pair);
	}
	template <typename TPair>
	inline typename TPair::T2 const & objectOf(TPair const &pair) {
		return getValueI2(pair);
	}


	//////////////////////////////////////////////////////////////////////////////
	//

	template <typename TString, typename TLess = ::std::less<typename Value<TString>::Type::T1> >
	struct SequenceMap;


	template <typename TString>
	struct Value< SequenceMap<TString> > {
		typedef typename Value<TString>::Type Type;
	};
	template <typename TString>
	struct Size< SequenceMap<TString> > {
		typedef typename Size<TString>::Type Type;
	};

	template <typename TString>
	inline typename Size< SequenceMap<TString> >::Type length(SequenceMap<TString> const &set) {
		return length(set.string);
	}



	template <typename TKey, typename TObject, typename TSpec, typename TLess>
	struct SequenceMap< String< Pair<TKey, TObject>, TSpec >, TLess > {
		typedef Pair<TKey, TObject>						TValue;
		typedef String< Pair<TKey, TObject>, TSpec >	TSequence;
		typedef typename TSequence::size_type				TSize;

		TKey				maxKey;
		TLess				comp;
		Holder<TSequence>	string;

		SequenceMap():
			maxKey(0)	{}

		SequenceMap(TLess const &_comp):
			maxKey(0),
			comp(_comp) {}
	};


	//////////////////////////////////////////////////////////////////////////////
	//

//	template <typename TKey>
//	struct Map {
//		typedef ::std::map<TKey> Type;
//	};
	template <typename TKey, typename TObject>
	class Map< Pair<TKey, TObject> > {
//		typedef ::std::set< Pair<TKey, TObject>, SetLess< Pair<TKey, TObject> > > Type;
	};



	//////////////////////////////////////////////////////////////////////////////
	//

	template <typename TString, typename TLess>
	struct Iterator< SequenceMap<TString, TLess> > {
		typedef typename Iterator<TString>::Type	Type;
	};

	template <typename TString, typename TLess>
	struct Iterator< SequenceMap<TString, TLess> const > {
		typedef typename Iterator<TString const>::Type	Type;
	};


	//////////////////////////////////////////////////////////////////////////////
	//

	template <typename TString, typename TLess>
	inline void clear(SequenceMap<TString, TLess> &map) {
		clear(map.string);
		map.maxKey = 0;
	}

	template <typename TKey, typename TString, typename TLess>
	inline typename Iterator< SequenceMap<TString, TLess> >::Type find(TKey const &key, SequenceMap<TString, TLess> &map) {

		typedef typename Size< SequenceMap<TString, TLess> >::Type		TSize;
		typedef typename Iterator< SequenceMap<TString, TLess> >::Type	TIter;
		
		// accelerate binary search
		if (map.comp(map.maxKey, key))
			return end(map.string);

		TIter First_ = begin(map.string);
		TSize Count_ = end(map.string) - First_;

		for (; 0 < Count_; )
		{	// divide and conquer, find half that contains answer
			TSize Count2_ = Count_ / 2;
			TIter Mid_ = First_ + Count2_;

			if (map.comp((*Mid_).i1, key))
				First_ = ++Mid_, Count_ -= Count2_ + 1;
			else
				Count_ = Count2_;
		}
		return (First_);
	}

	template <typename TPair, typename TString, typename TLess>
	inline void insert(TPair &pair, SequenceMap<TString> &map)
	{
		typedef typename Iterator< SequenceMap<TString, TLess> >::Type	TIter;

		// accelerate binary search
		if (map.comp(map.maxKey, pair.i1)) {
// whatever
			map.maxKey = pair.i1;
			return append(map.string);
		}

		TIter iter = find(pair.i1, map);
// whatever
		insert(pair, map);
	}


	template <typename TKey, typename TString>
	inline void erase(TKey const &/*key*/, SequenceMap<TString> &/*map*/)
	{
            // TODO(holtgrew): Commented out the following, no idea how to make it work.
            /* 
		if (!max.comp(map.maxKey, key) && !max.comp(key, map.maxKey)) {
			// delete last element
			if (length(map.string))
				map.maxKey = *(map.string.end() - 1).i2;
			else
				map.maxKey = 0;
		}

		typedef typename Iterator< SequenceMap<TString, TLess> >::Type	TIter;
		TIter iter = find(pair.i1, map);
		delete iter;
            */
	}
}

#endif

