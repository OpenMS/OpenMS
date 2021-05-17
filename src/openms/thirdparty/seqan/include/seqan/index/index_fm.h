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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================

//SEQAN_NO_DDDOC:do not generate documentation for this file

#ifndef INDEX_FM_H
#define INDEX_FM_H

namespace seqan {

// ==========================================================================
// Forwards
// ==========================================================================

template <typename TChar, typename TSpec>
class PrefixSumTable;

// WT = WaveletTree
/**
.Tag.WT
..summary:Tag that specifies the @Spec.FMIndex@ to use a wavelet tree as the occurrence table.
..cat:Index
*/

template <typename TSpec = void>
class WT;

template <typename TOccSpec = WT<>, typename TSpec = void>
class FMIndex;

struct FinderFMIndex_;
typedef Tag<FinderFMIndex_> FinderFMIndex;

template <typename TText, typename TOccSpec, typename TSpec>
struct DefaultFinder<Index<TText, FMIndex<TOccSpec, TSpec> > >
{
	typedef FinderFMIndex Type;
};

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

struct FibrePrefixSumTable_;
struct FibreSA_;
struct FibreTempSA_;
struct FibreText_;
struct FibreLfTable_;
struct FibreSaLfTable_;
struct CompressText_;
struct Sentinel_;
struct Sentinels_;

typedef Tag<FibrePrefixSumTable_> const FibrePrefixSumTable;
typedef Tag<FibreSA_> const             FibreSA;
typedef Tag<FibreTempSA_> const         FibreTempSA;
typedef Tag<FibreText_> const           FibreText;
typedef Tag<FibreLfTable_> const        FibreLfTable;
typedef Tag<FibreSaLfTable_> const      FibreSaLfTable;
typedef Tag<CompressText_> const        CompressText;
typedef Tag<Sentinel_> const            Sentinel;
typedef Tag<Sentinels_> const           Sentinels;


// FM index fibres
/**
.Tag.FM Index Fibres
..summary:Tag to select a specific fibre of a @Spec.FMIndex@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a FM index.
..cat:Index

..tag.FibrePrefixSumTable:The prefix sum table of the index.
..tag.FibreSA:The compressed suffix array of the text.
..tag.FibreText:The original text of the index.
..tag.FibreLfTable:The lf table.
..tag.FibreSaLfTable:The lf table as well as the compressed suffix array.
...remarks:This tag can only be used with the functions @Function.indexRequire@ or @Function.indexSupplied@.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index_fm.h
*/

/**
.Tag.CompressText
..cat:Index
..summary:Tag to select a FM index variant that can be used such that it is 
not necessary to store the text after index construction. This index is very
space efficient.
*/

// ==========================================================================
// Metafunctions
// ==========================================================================

// ----------------------------------------------------------------------------
// Metafunction Fibre
// ----------------------------------------------------------------------------
    
/*
.Metafunction.Fibre:
..summary:Type of a specific FMIndex member (fibre).
..signature:Fibre<Index<TText, FMIndex<TSentinelRankDictionary, TSpec> >, TFibreSpec>::Type
..class:Spec.FMIndex
..cat:Index
..param.TText:The text type.
...type:Class.String
...type:Class.StringSet
..param.TSentinelRankDictionary:The type of the sentinel rank dictionary.
...type:Class.SentinelRankDictionary
..param.TSpec:Tag to specify a certain variant of the FM index.
...type:Tag.CompressText
...default;$void$
..param.TFibreSpec:Tag to specify the fibre.
...type:Tag.FM Index Fibres
..returns:Fibre type.
..remarks:Some containers, such as @Spec.FMIndex@, can be seen as a bundle consisting of various fibres. Because not 
every table is a fibre we did not call them tables, however, in many cases one can think of fibres as tables. The 
fibre interface was designed to unify the access to the members of the different fibres.
To get a reference or the type of a specific fibre use @Function.getFibre@ or @Metafunction.Fibre@.		
..include:seqan/index.h
*/

template <typename TText, typename TWaveletTreeSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<WT<TWaveletTreeSpec>, TSpec> >, FibreOccTable>
{
    typedef typename Value<TText>::Type TValue_;
	typedef SentinelRankDictionary<RankDictionary<WaveletTree<TValue_> >, Sentinel> Type;
};

template <typename TText, typename TStringSetSpec, typename TWaveletTreeSpec, typename TSpec>
struct Fibre<Index<StringSet<TText, TStringSetSpec>, FMIndex<WT<TWaveletTreeSpec>, TSpec > >, FibreOccTable>
{
    typedef typename Value<TText>::Type TValue_;
    typedef SentinelRankDictionary<RankDictionary<WaveletTree<TValue_> >, Sentinels> Type;
};

template <typename TText, typename TWaveletTreeSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<SBM<TWaveletTreeSpec>, TSpec> >, FibreOccTable>
{
    typedef typename Value<Index<TText, FMIndex<WT<TWaveletTreeSpec>, TSpec> > >::Type TValue_;
	typedef SentinelRankDictionary<RankDictionary<SequenceBitMask<TValue_> >, Sentinel> Type;
};

template <typename TText, typename TStringSetSpec, typename TWaveletTreeSpec, typename TSpec>
struct Fibre<Index<StringSet<TText, TStringSetSpec>, FMIndex<SBM<TWaveletTreeSpec>, TSpec > >, FibreOccTable>
{
    typedef typename Value<TText>::Type TValue_;
    typedef SentinelRankDictionary<RankDictionary<SequenceBitMask<TValue_> >, Sentinels> Type;
};

template <typename TText, typename TOccSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreLfTable>
{
    typedef typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreOccTable>::Type        TOccTable_;
	typedef typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibrePrefixSumTable>::Type  TPrefixSumTable_;
	typedef LfTable<TOccTable_, TPrefixSumTable_>   	                                        Type;
};

template <typename TText, typename TOccSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<TOccSpec, TSpec> > const, FibreLfTable>
{
    typedef typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> > const, FibreOccTable>::Type          TOccTable_;
	typedef typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> > const, FibrePrefixSumTable>::Type    TPrefixSumTable_;
	typedef LfTable<TOccTable_, TPrefixSumTable_>   	                                                Type;
};

template <typename TText, typename TOccSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibrePrefixSumTable>
{
    typedef typename Value<Index<TText, FMIndex<TOccSpec, TSpec> > >::Type  TChar_;
    typedef typename MakeUnsigned<TChar_>::Type                             TUChar_;
	typedef PrefixSumTable<TUChar_, void>                                   Type;
};

template <typename TText, typename TOccSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<TOccSpec, TSpec> > const, FibrePrefixSumTable>
{
    typedef typename Value<Index<TText, FMIndex<TOccSpec, TSpec> > const>::Type TChar_;
    typedef typename MakeUnsigned<TChar_>::Type                                 TUChar_;
	typedef PrefixSumTable<TUChar_, void>                                       Type;
};

template <typename TText, typename TOccSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreSA>
{
	typedef typename SAValue<Index<TText, FMIndex<TOccSpec, TSpec> > >::Type                TSAValue_;
	typedef SparseString<String<TSAValue_>, void>                                           TSparseString_;
	typedef typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreLfTable>::Type     TLfTable_;
	typedef CompressedSA<TSparseString_, TLfTable_, void>                                   Type;
};

template <typename TText, typename TOccSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreTempSA>
{
	typedef typename SAValue<Index<TText, FMIndex<TOccSpec, TSpec> > >::Type    TSAValue;
	typedef String<TSAValue, External<ExternalConfigLarge<> > >                 Type;
};

// ==========================================================================
// Classes
// ==========================================================================

// ----------------------------------------------------------------------------
// Class FmIndexInfo_ 
// ----------------------------------------------------------------------------

// Stores the information about an FM index file bundle and is written to the .fma file.

#ifdef PLATFORM_WINDOWS
    #pragma pack(push,1)
#endif

struct FmIndexInfo_
{
    // The compression factor.
    __uint32 compressionFactor;
    // The sizeof(TSAEntry) values for suffix array entries.
    __uint32 sizeOfSAEntry;
    // The length of the genome.
    __uint64 genomeLength;
}

#ifndef PLATFORM_WINDOWS
    __attribute__((packed))
#endif
    ;
#ifdef PLATFORM_WINDOWS
      #pragma pack(pop)
#endif

// ----------------------------------------------------------------------------
// Spec FMIndex 
// ----------------------------------------------------------------------------

/**
.Spec.FMIndex:
..summary:An index based on the Burrows-Wheeler transform.
..cat:Index
..general:Class.Index
..signature:Index<TText, FMIndex<TOccSpec, TSpec> >
..param.TText:The text type.
...type:Class.String
...type:Class.StringSet
..param.TOccSpec:Occurrence table specialisation. 
...type:Tag.WT
...type:Tag.SBM
...remarks:The tags are really shortcuts for the different @Class.SentinelRankDictionary@s
...default:Tag.WT
..param.TSpec:FM index specialisation.
...type:Tag.CompressText
...default:void
..include:seqan/index.h
*/

template <typename TText, typename TOccSpec, typename TSpec>
class Index<TText, FMIndex<TOccSpec, TSpec> >
{
    public:
    Holder<typename Fibre<Index, FibreText>::Type>  text;
	typename Fibre<Index, FibreLfTable>::Type       lfTable;
	typename Fibre<Index, FibreSA>::Type            compressedSA;
	typename Size<TText>::Type                      n;
    unsigned                                        compressionFactor; 

	Index() :
		n(0)
	{}

	Index(TText & text, unsigned compressionFactor = 10) :
	    text(text),
		n(_computeBwtLength(text)),
		compressionFactor(compressionFactor)
	{}

	inline bool operator==(const Index & b) const
    {
        return lfTable == b.lfTable &&
               compressedSA == b.compressedSA &&
               n == b.n &&
               compressionFactor == b.compressionFactor;
    }
};

// ==========================================================================
// Functions
// ==========================================================================
// ----------------------------------------------------------------------------
// Function clear
// ----------------------------------------------------------------------------

// Already documented
template <typename TText, typename TOccSpec, typename TSpec>
inline void clear(Index<TText, FMIndex<TOccSpec, TSpec> > & index)
{
    clear(getFibre(index, FibreLfTable()));
    clear(getFibre(index, FibreSA()));
}

// ----------------------------------------------------------------------------
// Helper function _computeBwtLength
// ----------------------------------------------------------------------------

// This function computes the length of the bwt string.
template <typename TText>
unsigned _computeBwtLength(TText const & text)
{
    return length(text) + 1;
}

// This function computes the length of the bwt string.
template <typename TText, typename TSetSpec>
unsigned _computeBwtLength(StringSet<TText, TSetSpec> const & text)
{
    return lengthSum(text) + countSequences(text);
}

// ----------------------------------------------------------------------------
// Helper function _createBwTable
// ----------------------------------------------------------------------------

// This function computes the BWT of a text. Note that the sentinel sign is substituted and its position stored.
// The function is tested implicitly in lfTableLfMapping() in test_index_fm.h
template <typename TBwt, typename TSentinelPosition, typename TText, typename TSA, typename TSentinelSub>
inline void _createBwTable(TBwt & bwt, TSentinelPosition & sentinelPos, TText const & text, TSA const & sa,
		                   TSentinelSub const sentinelSub)
{
	//typedefs
	typedef typename GetValue<TSA>::Type                    TSAValue;
    typedef typename Iterator<TSA const, Standard>::Type    TSAIter;
    typedef typename Iterator<TBwt, Standard>::Type         TBwtIter;

	//little helpers
    TSAIter saIt = begin(sa, Standard());
    TSAIter saItEnd = end(sa, Standard());
    TBwtIter bwtIt = begin(bwt, Standard());
    
	assignValue(bwtIt, back(text));
	for (; saIt != saItEnd; ++saIt)
	{
        ++bwtIt;
		TSAValue pos = getValue(saIt);
		if (pos != 0)
			assignValue(bwtIt, getValue(text, pos - 1));
		else
		{
			assignValue(bwtIt, sentinelSub);
			sentinelPos = bwtIt - begin(bwt, Standard());
		}
	}
}

// This function computes the BWT of a text. Note that the sentinel sign is substituted and its position stored.
// The function is tested implicitly in lfTableLfMapping() in test_index_fm.h
template <typename TBwt, typename TSentinelPosition, typename TText, typename TSetSpec, typename TSA, typename TSentinelSub>
inline void _createBwTable(TBwt & bwt, TSentinelPosition & sentinelPos, StringSet<TText, TSetSpec> const & text,
		                   TSA const & sa, TSentinelSub const sentinelSub)
{
    typedef typename Value<TSA>::Type                       TSAValue;
    typedef typename Size<TSA>::Type                        TSize;
    typedef typename Iterator<TSA const, Standard>::Type    TSAIter;
    typedef typename Iterator<TBwt, Standard>::Type         TBwtIter;

    // little Helpers
    TSize seqNum = countSequences(text);
    TSize totalLen = lengthSum(text);

    resize(sentinelPos, seqNum + totalLen, Exact());
    
    TSAIter saIt = begin(sa, Standard());
    TSAIter saItEnd = end(sa, Standard());
    TBwtIter bwtItBeg = begin(bwt, Standard());
    TBwtIter bwtIt = bwtItBeg;

    // fill the sentinel positions (they are all at the beginning of the bwt)
    for (TSize i = 1; i <= seqNum; ++i, ++bwtIt)
        assignValue(bwtIt, back(text[seqNum - i]));

    // compute the rest of the bwt
    for (; saIt != saItEnd; ++saIt, ++bwtIt)
    {
        TSAValue pos;    // = SA[i];
        posLocalize(pos, getValue(saIt), stringSetLimits(text));
        if (getSeqOffset(pos) != 0)
            assignValue(bwtIt, getValue(getValue(text, getSeqNo(pos)), getSeqOffset(pos) - 1));
        else
        {
            assignValue(bwtIt, sentinelSub);
            setBit(sentinelPos, bwtIt - bwtItBeg);
        }
    }

    // update the auxiliary rank support bit string information.
    _updateRanks(sentinelPos);
}

// ----------------------------------------------------------------------------
// Helper function _determineSentinelSubstitute
// ----------------------------------------------------------------------------

// This function determines the '$' substitute.
// The character with the smallest number of occurrences greater 0 is chosen.
template <typename TPrefixSumTable, typename TChar>
inline void _determineSentinelSubstitute(TPrefixSumTable const & pst, TChar & sub)
{
    typedef typename RemoveConst<TPrefixSumTable>::Type TNonConstPrefixSumTable;
	typedef typename Value<TNonConstPrefixSumTable>::Type TValue;

	TValue min = getPrefixSum(pst, length(pst) - 1);
	unsigned pos = length(pst) - 1;
	for (unsigned i = 0; i < length(pst) - 1; ++i)
	{
		unsigned diff = pst[i + 1] - pst[i];
		if (diff != 0 && diff < min)
		{
			min = diff;
			pos = i;
		}
	}
	sub = getCharacter(pst, pos);
}

// ----------------------------------------------------------------------------
// Function empty
// ----------------------------------------------------------------------------

// This function checks whether the index is empty. Its already documented.
template <typename TText, typename TIndexSpec, typename TFMISpeedEnhancement>
inline bool empty(Index<TText, FMIndex<TIndexSpec, TFMISpeedEnhancement> > const & index)
{
    return empty(getFibre(index, FibreLfTable()))
            && empty(getFibre(index, FibreSA()));
}

// ----------------------------------------------------------------------------
// Function _findFirstIndex
// ----------------------------------------------------------------------------

// This function is used by the finder interface. It initializes the range of the finder.
template <typename TText, typename TPattern, typename TIndexSpec, typename TFMISpeedEnhancement>
inline void
_findFirstIndex(Finder<Index<TText, FMIndex<TIndexSpec, TFMISpeedEnhancement> >, FinderFMIndex> & finder,
		        TPattern const & pattern, FinderFMIndex const &)
{
	typedef Index<TText, FMIndex<TIndexSpec, TFMISpeedEnhancement> >    TIndex;

	TIndex & index = haystack(finder);

	indexRequire(index, FibreSaLfTable());
    setContainer(finder.range.i1, getFibre(container(finder), FibreSA()));
    setContainer(finder.range.i2, getFibre(container(finder), FibreSA()));

	_range(index, pattern, finder.range);
}

// ----------------------------------------------------------------------------
// Function getFibre
// ----------------------------------------------------------------------------

/**
.Function.FMIndex#getFibre:
..summary:Returns a specific fibre of a fm index.
..signature:getFibre(index, fibreTag)
..class:Spec.FMIndex
..cat:Index
..param.index:The index holding the fibre.
...type:Spec.FMIndex
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.FM Index Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
..include:seqan/index.h
*/

template <typename TText, typename TIndexSpec, typename TFMISpeedEnhancement>
typename Fibre<Index<TText, FMIndex<TIndexSpec, TFMISpeedEnhancement> >, FibreLfTable >::Type &
getFibre(Index<TText, FMIndex<TIndexSpec, TFMISpeedEnhancement> > & index, FibreLfTable /*tag*/)
{
	return index.lfTable;
}

template <typename TText, typename TIndexSpec, typename TFMISpeedEnhancement>
typename Fibre<Index<TText, FMIndex<TIndexSpec, TFMISpeedEnhancement> >, FibreLfTable >::Type const &
getFibre(Index<TText, FMIndex<TIndexSpec, TFMISpeedEnhancement> > const & index, FibreLfTable /*tag*/)
{
	return index.lfTable;
}

template <typename TText, typename TIndexSpec, typename TFMISpeedEnhancement>
typename Fibre<Index<TText, FMIndex<TIndexSpec, TFMISpeedEnhancement> >, FibreSA >::Type &
getFibre(Index<TText, FMIndex<TIndexSpec, TFMISpeedEnhancement> > & index, FibreSA /*tag*/)
{
	return index.compressedSA;
}

template <typename TText, typename TIndexSpec, typename TFMISpeedEnhancement>
typename Fibre<Index<TText, FMIndex<TIndexSpec, TFMISpeedEnhancement> >, FibreSA >::Type const &
getFibre(Index<TText, FMIndex<TIndexSpec, TFMISpeedEnhancement> > const & index, FibreSA /*tag*/)
{
	return index.compressedSA;
}

// ----------------------------------------------------------------------------
// Function toSuffixPosition
// ----------------------------------------------------------------------------

/**
.Function.FMIndex#toSuffixPosition
..class:Spec.FMIndex
..summary:This function computes the position of a specified position in the suffix array (additionally containing 
entries for the sentinels. The returned position correspond to the suffix array of the original text without sentinels.
..signature:toSuffixPosition(fmIndex, pos, offset)
..param.fmIndex:The FM index.
...type:Spec.FMIndex
..param.pos:The position in the suffix array of the fm index (with sentinels).
...type:Concept.UnsignedIntegerConcept
..param.offset:The number of sequences in the original text.
...type:Concept.UnsignedIntegerConcept
..include:seqan/index.h
*/

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TPos, typename TSize>
inline typename SAValue<Index<TText, FMIndex<TOccSpec, TIndexSpec > > >::Type
toSuffixPosition(Index<TText, FMIndex<TOccSpec, TIndexSpec > > & index, TPos i, TSize offset)
{
    SEQAN_ASSERT_GEQ(suffixLength(i, index), offset);
    setSeqOffset(i, suffixLength(i, index) - offset);
    return i;
}

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TPos, typename TSize>
inline typename SAValue<Index<TText, FMIndex<TOccSpec, TIndexSpec > > const>::Type
toSuffixPosition(Index<TText, FMIndex<TOccSpec, TIndexSpec > > const & index, TPos i, TSize offset)
{
    SEQAN_ASSERT_GEQ(suffixLength(i, index), offset);
    setSeqOffset(i, suffixLength(i, index) - offset);
    return i;
}

// ----------------------------------------------------------------------------
// Helper function _getFrequencies
// ----------------------------------------------------------------------------

// This function determines the number of the different characters in the text.
template <typename TText, typename TSetSpec, typename TFreq>
inline void _getFrequencies(TFreq & freq,
						    StringSet<TText, TSetSpec> const & text)
{
	typedef typename Value<TText>::Type TChar;
	resize(freq, ValueSize<TChar>::VALUE, 0, Exact());

    typedef typename Size<TText>::Type TSize;
	for (TSize i = 0; i < length(text); ++i)
	    for (TSize j = 0; j < length(text[i]); ++j)
		    ++freq[getCharacterPosition(freq, text[i][j])];
}

// This function determines the number of the different characters in the text.
template <typename TText, typename TFreq>
inline void _getFrequencies(TFreq & freq,
		   	   	   	   	   TText const & text)
{
	typedef typename Value<TText>::Type TChar;
	resize(freq, ValueSize<TChar>::VALUE, 0, Exact());

    typedef typename Size<TText>::Type TSize;
	for (TSize i = 0; i < length(text); ++i)
		++freq[getCharacterPosition(freq, text[i])];
}

// ----------------------------------------------------------------------------
// Function _indexCreateSA
// ----------------------------------------------------------------------------

// This function computes the full and compressed suffix array. 
// Note, in contrast to indexCreate(index, FibreSA()) the full suffix array is also computed.
template <typename TText, typename TIndexSpec, typename TSpec, typename TSA>
inline bool _indexCreateSA(Index<TText, FMIndex<TIndexSpec, TSpec> > & index, TSA & fullSa, TText const & text)
{
	typedef Index<TText, FMIndex<TIndexSpec, TSpec> > TIndex;
	typedef typename Fibre<TIndex, FibreSA>::Type TCompressedSA;

    // TODO(singer): If there is a lfTable we do not need the Skew7
	// create the fulle sa

    // create the compressed sa
	TCompressedSA & compressedSA = getFibre(index, FibreSA());
    setLfTable(compressedSA, getFibre(index, FibreLfTable()));

    unsigned numSentinel = countSequences(text);
    createCompressedSa(compressedSA, fullSa, index.compressionFactor, numSentinel); 

	return true;
}

// ----------------------------------------------------------------------------
// Function _indexCreateLfTables
// ----------------------------------------------------------------------------

// This function creates all table of the lf table given a text and a suffix array.
template <typename TIndexSpec, typename TSpec, typename TText, typename TSA>
inline bool _indexCreateLfTables(Index<TText, FMIndex<TIndexSpec, TSpec> > & index, TText & text, TSA const & sa)
{
	typedef Index<TText, FMIndex<TIndexSpec, TSpec> >		        TIndex;
	typedef typename Fibre<TIndex, FibreLfTable>::Type              TLfTable;
	typedef typename Fibre<TLfTable, FibreOccTable>::Type           TOccTable;
	typedef typename Fibre<TOccTable, FibreSentinelPosition>::Type  TSentinelPosition;
	typedef typename Value<TIndex>::Type						    TAlphabet;

	createPrefixSumTable(index.lfTable.prefixSumTable, text);

	TAlphabet sentinelSub(0);
	_determineSentinelSubstitute(index.lfTable.prefixSumTable, sentinelSub);

	String<TAlphabet> bwt;
	resize(bwt, index.n, Exact());
	TSentinelPosition sentinelPos = _setDefaultSentinelPosition(length(bwt), TSentinelPosition());

	_createBwTable(bwt, sentinelPos, text, sa, sentinelSub);

    createSentinelRankDictionary(index.lfTable, bwt, sentinelSub, sentinelPos);

	_insertSentinel(index.lfTable.prefixSumTable, countSequences(text));

	return true;
}

// ----------------------------------------------------------------------------
// Function indexCreate
// ----------------------------------------------------------------------------

/**
.Function.FMIndex#indexCreate
..summary:Creates a specific @Metafunction.Fibre@.
..signature:indexCreate(index, fibreTag)
..param.index:The index to be created.
...type:Spec.FMIndex
..param.fibreTag:The fibre of the index to be computed.
...type:Tag.FM Index Fibres.tag.FibreSaLfTable
..remarks:If you call this function on the compressed text version of the FM index
you will get an error message: "Logic error. It is not possible to create this index without a text."
*/

// This function creates the index.
template <typename TText, typename TIndexSpec, typename TSpec>
inline bool _indexCreate(Index<TText, FMIndex<TIndexSpec, TSpec > > & index, TText & text)
{
	typedef Index<TText, FMIndex<TIndexSpec, TSpec> >   TIndex;
	typedef typename Fibre<TIndex, FibreTempSA>::Type   TTempSA;

    if (empty(text))
        return false;

    TTempSA tempSA;
    
	resize(tempSA, length(text), Exact());
	createSuffixArray(tempSA,
			text,
			Skew7());


	// create the compressed SA
	_indexCreateSA(index, tempSA, text);

	// create the lf table
	_indexCreateLfTables(index, text, tempSA);

	return true;
}


template <typename TText, typename TIndexSpec, typename TSpec>
inline bool indexCreate(Index<TText, FMIndex<TIndexSpec, TSpec> > & index, FibreSaLfTable const)
{
    return _indexCreate(index, getFibre(index, FibreText()));
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline bool indexCreate(Index<TText, FMIndex<TIndexSpec, TSpec> > & index)
{
    return indexCreate(index, FibreSaLfTable());
}

// ----------------------------------------------------------------------------
// Function indexSupplied
// ----------------------------------------------------------------------------

/**
.Function.FMIndex#indexSupplied:
..summary:Returns whether a specific @Metafunction.Fibre@ is present.
..param.fibreTag:
...type:Tag.FM Index Fibres
*/
template <typename TText, typename TIndexSpec, typename TSpec>
inline bool indexSupplied(Index<TText, FMIndex<TIndexSpec, TSpec > > & index, FibreSaLfTable const)
{
    return !(empty(getFibre(index, FibreSA())) || empty(getFibre(index, FibreLfTable())));
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline bool indexSupplied(Index<TText, FMIndex<TIndexSpec, TSpec > > const & index, FibreSaLfTable const)
{
    return !(empty(getFibre(index, FibreSA())) || empty(getFibre(index, FibreLfTable())));
}

// ----------------------------------------------------------------------------
// Function _range
// ----------------------------------------------------------------------------

// This function computes a range in the suffix array whose entries point to location
// in the text where the pattern occurs. 
template <typename TText, typename TOccSpec, typename TSpec, typename TPattern, typename TIter, typename TPairSpec>
inline void _range(const Index<TText, FMIndex<TOccSpec, TSpec> > & index, const TPattern & pattern, 
                   Pair<TIter, TPairSpec> & range)
{
    typedef Index<TText, FMIndex<TOccSpec, TSpec> >             TIndex;
    typedef typename Value<TIndex>::Type                        TAlphabet;
    typedef typename ValueSize<TAlphabet>::Type                 TAlphabetSize;
    typedef typename Size<TIndex>::Type                         TSize;
 	typedef typename Value<TPattern>::Type                      TChar;

    if (empty(pattern))
    {
	    setPosition(range.i1, countSequences(index));
	    setPosition(range.i2, index.n);
    }

	TSize i = length(pattern) - 1;
	TChar letter = pattern[i];

    // initilization
	TAlphabetSize letterPosition = getCharacterPosition(index.lfTable.prefixSumTable, letter);
	TSize sp = getPrefixSum(index.lfTable.prefixSumTable, letterPosition);
	TSize ep = getPrefixSum(index.lfTable.prefixSumTable, letterPosition + 1) - 1;

	// the search as proposed by Ferragina and Manzini
	while ((sp <= ep) && (i > 0))
	{
	    --i;
		letter = pattern[i];
		letterPosition = getCharacterPosition(index.lfTable.prefixSumTable, letter);
		TSize prefixSum = getPrefixSum(index.lfTable.prefixSumTable, letterPosition);
		sp = prefixSum + countOccurrences(index.lfTable.occTable, letter, sp - 1);
		ep = prefixSum + countOccurrences(index.lfTable.occTable, letter, ep) - 1;
	}

    setPosition(range.i1, sp);
	setPosition(range.i2, ep + 1);
}

// ----------------------------------------------------------------------------
// Function open
// ----------------------------------------------------------------------------

/**
.Function.Index#open
..class:Class.Index
..summary:This functions loads a dictionary from disk.
..signature:open(dictionary, fileName [, openMode])
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
*/

// This function can be used to open a previously saved index.
template <typename TText, typename TOccSpec, typename TSpec>
inline bool open(Index<TText, FMIndex<TOccSpec, TSpec> > & index, const char * fileName, int openMode)
{
    String<char> name;

//    typedef Index<TText, FMIndex<TOccSpec, TSpec> > TIndex;
//    typedef typename Fibre<TIndex, FibreSA>::Type TSAFibre;
//    typedef typename Value<TSAFibre>::Type TSAValue;

    String<FmIndexInfo_> infoString;

    name = fileName;    append(name, ".txt");
    if (!open(getFibre(index, FibreText()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".sa");
    if (!open(getFibre(index, FibreSA()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".lf");
    if (!open(getFibre(index, FibreLfTable()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".fma");
    if (!open(infoString, toCString(name), openMode)) return false;

    // TODO(weese): Why do we check the SA size here? Everywhere else we require the types being the same for reading and writing
    //              Moreover a correct size is not indicating a correct type
    // Check that the size of the SA entries is correct.
//    if (infoString[0].sizeOfSAEntry != 0 && infoString[0].sizeOfSAEntry != sizeof(TSAValue))
//        return false;

    index.compressionFactor = infoString[0].compressionFactor;
    index.n = infoString[0].genomeLength;
    getFibre(index, FibreSA()).lfTable = & getFibre(index, FibreLfTable());

    return true;
}

// This function can be used to open a previously saved index.
template <typename TText, typename TOccSpec, typename TSpec>
inline bool open(Index<TText, FMIndex<TOccSpec, TSpec> > & index, const char * fileName)
{
    return open(index, fileName, DefaultOpenMode<Index<TText, FMIndex<TOccSpec, TSpec> > >::VALUE);
}

// ----------------------------------------------------------------------------
// Function save
// ----------------------------------------------------------------------------

/**
.Function.Index#save
..class:Class.Index
..summary:This functions saves an index to disk.
..signature:save(index, fileName [, openMode])
..param.index:The index.
...type:Class.RankDictionary
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
*/

template <typename TText, typename TOccSpec, typename TSpec>
inline bool save(Index<TText, FMIndex<TOccSpec, TSpec> > const & index, const char * fileName, int openMode)
{
    String<char> name;

    typedef Index<TText, FMIndex<TOccSpec, TSpec> > TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type TSAFibre;
    typedef typename Value<TSAFibre>::Type TSAValue;

    String<FmIndexInfo_> infoString;
    FmIndexInfo_ info = { index.compressionFactor, sizeof(TSAValue), index.n };
    appendValue(infoString, info);

    name = fileName;    append(name, ".txt");
    if (!save(getFibre(index, FibreText()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".sa");
    if (!save(getFibre(index, FibreSA()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".lf");
    if (!save(getFibre(index, FibreLfTable()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".fma");
    if (!save(infoString, toCString(name), openMode)) return false;

    return true;
}

// This function can be used to save an index on disk.
template <typename TText, typename TOccSpec, typename TSpec>
inline bool save(Index<TText, FMIndex<TOccSpec, TSpec> > const & index, const char * fileName)
{
    return save(index, fileName, DefaultOpenMode<Index<TText, FMIndex<TOccSpec, TSpec> > >::VALUE);
}

}
#endif // INDEX_FM_H
