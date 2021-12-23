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

#ifndef INDEX_FM_RANKDICTIONARY_SBM
#define INDEX_FM_RANKDICTIONARY_SBM

namespace seqan {

// ==========================================================================
// Forwards
// ==========================================================================

template <typename TValue>
class SequenceBitMask;

template<typename TSpec> 
class RankDictionary;

// ==========================================================================
// Tags
// ==========================================================================

struct FibreBitStrings_;
typedef Tag<FibreBitStrings_> const FibreBitStrings;

/**
.Tag.SBM
..summary:Tag that specifies the @Spec.FMIndex@ to use a StringSet of rank support bis strings as the occurrence table.
..cat:Index
*/
template <typename TSpec = void>
class SBM;

// ==========================================================================
// Metafunctions
// ==========================================================================

/**
.Spec.SequenceBitMask Fibres
..cat:Index
..summary:Tag to select a specific fibre of a SequenceBitMask.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a SequenceBitMask.

..DISABLED.tag.FibreBitStrings:The string set containing a bit string for each node.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index.h
*/

// ----------------------------------------------------------------------------
// Metafunction Fibre
// ----------------------------------------------------------------------------

template <typename TValue>
struct Fibre<RankDictionary<SequenceBitMask<TValue> >, FibreBitStrings>
{
    typedef StringSet<RankSupportBitString<void> > Type;
};

template <typename TValue>
struct Fibre<RankDictionary<SequenceBitMask<TValue> > const, FibreBitStrings>
{
    typedef typename Fibre<RankDictionary<WaveletTree<TValue> >, FibreBitStrings>::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename TValue>
struct Size<RankDictionary<SequenceBitMask<TValue> > >
{
    typedef typename Size<String<TValue> >::Type Type;
};

template <typename TValue>
struct Size<RankDictionary<SequenceBitMask<TValue> > const> :
    public Size<RankDictionary<SequenceBitMask<TValue> > > {};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TValue>
struct Value<RankDictionary<SequenceBitMask<TValue> > >
{
    typedef TValue Type;
};

template <typename TValue>
struct Value<RankDictionary<SequenceBitMask<TValue> > const> :
    public Value<RankDictionary<SequenceBitMask<TValue> > > {};

// ==========================================================================
// Classes
// ==========================================================================

// ----------------------------------------------------------------------------
// Spec SequenceBitMask
// ----------------------------------------------------------------------------

/**
.Spec.SequenceBitMask:
..cat:Index
..summary:The string set bit string dictionary is a string set of rank support bit strings for constant time acces
of the rank of a specified character at a specified position.
..signature:SequenceBitMask<TValue>
..param.TValue:The value type of the .
..include:seqan/index.h
..remarks:This data structure is optimized for very small alphabets, such as @Spec.Dna@ or @Spec.Dna5@. Consider using a @Spec.WaveletTree@ if your alphabet size is larger.
*/
template <typename TValue>
class RankDictionary<SequenceBitMask<TValue> >
{
    typedef typename Fibre<RankDictionary<SequenceBitMask<TValue> >, FibreBitStrings>::Type    TBitStrings;

public:
    TBitStrings bitStrings;

    RankDictionary() {}

    template <typename TText>
    RankDictionary(TText const & text)
    {
        createRankDictionary(*this, text);
    }

    template <typename TText, typename TFreqTable>
    RankDictionary(TText const & text, TFreqTable const & freqTable)
    {
        createRankDictionary(*this, text, freqTable);
    }

    template <typename TText, typename TFreqTable, typename TPrefixSumTable>
    RankDictionary(TText const & text, TFreqTable const & freqTable, TPrefixSumTable const & prefixSumTable)
    {
        createRankDictionary(*this, text, freqTable, prefixSumTable);
    }

    RankDictionary & operator=(RankDictionary const & other)
    {
        bitStrings = other.bitStrings;
        return *this;
    }

    bool operator==(RankDictionary const & b) const
    {
        typedef typename Size<typename Fibre<RankDictionary<SequenceBitMask<TValue> >, FibreBitStrings>::Type>::Type TSize;
        
        if (length(bitStrings) != length(b.bitStrings))
            return false;
        
        for (TSize i = 0; i < length(bitStrings); ++i)
            if (!(bitStrings[i] == b.bitStrings[i]))
                return false;
        
        return true;
    }
};

// ==========================================================================
//Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function clear
// ----------------------------------------------------------------------------

template <typename TValue>
inline void clear(RankDictionary<SequenceBitMask<TValue> > & dictionary)
{
    clear(getFibre(dictionary, FibreBitStrings()));
}

// ----------------------------------------------------------------------------
// Function empty
// ----------------------------------------------------------------------------

template <typename TValue>
inline bool empty(RankDictionary<SequenceBitMask<TValue> > const & dictionary)
{
    return empty(getFibre(dictionary, FibreBitStrings()));
}

// ----------------------------------------------------------------------------
// Function getValue
// ----------------------------------------------------------------------------

template <typename TValue, typename TPos>
inline TValue
getValue(RankDictionary<SequenceBitMask<TValue> > const & dictionary, TPos pos)
{
    typedef typename Fibre<RankDictionary<SequenceBitMask<TValue> > const, FibreBitStrings>::Type TBitStrings;

    TBitStrings & bitStrings = getFibre(dictionary, FibreBitStrings());
    for (unsigned i = 0; i < ValueSize<TValue>::VALUE - 1; ++i)
        if (isBitSet(bitStrings[i], pos))
            return i;

    return maxValue<TValue>();
}

template <typename TValue, typename TPos>
inline TValue
getValue(RankDictionary<SequenceBitMask<TValue> > & dictionary, TPos pos)
{
    return getValue(const_cast<RankDictionary<SequenceBitMask<TValue> > const &>(dictionary), pos);
}

// ----------------------------------------------------------------------------
// Function getFibre
// ----------------------------------------------------------------------------

///.Function.RankDictionary#getFibre.param.fibreTag.type:Spec.SequenceBitMask Fibres
template <typename TValue>
inline typename Fibre<RankDictionary<SequenceBitMask<TValue> >, FibreBitStrings>::Type &
getFibre(RankDictionary<SequenceBitMask<TValue> >& dictionary, FibreBitStrings)
{
    return dictionary.bitStrings;
}

template <typename TValue>
inline typename Fibre<RankDictionary<SequenceBitMask<TValue> >, FibreBitStrings>::Type const &
getFibre(RankDictionary<SequenceBitMask<TValue> > const &dictionary, FibreBitStrings)
{
    return dictionary.bitStrings;
}

// ----------------------------------------------------------------------------
// Function countOccurrences
// ----------------------------------------------------------------------------

// This functions computes the number of occurrences of a specified character
// up to a specified position.
template <typename TValue, typename TCharIn, typename TPos>
inline typename Size<RankDictionary<SequenceBitMask<TValue> > const>::Type
countOccurrences(RankDictionary<SequenceBitMask<TValue> > const & dictionary,
                                 TCharIn const character, TPos const pos)
{
    return getRank(dictionary.bitStrings[ordValue(character)], pos);
}

template < typename TValue, typename TCharIn, typename TPos>
inline typename Size<RankDictionary<SequenceBitMask<TValue> > >::Type
countOccurrences(RankDictionary<SequenceBitMask<TValue> > & dictionary, TCharIn const character,
                                 TPos const pos)
{
    return countOccurrences(const_cast<RankDictionary<SequenceBitMask<TValue> > const &>(dictionary), character, pos);

    //return getRank(dictionary.bitStrings[ordValue(character)], pos);
}

// ----------------------------------------------------------------------------
// Function createRankDictionary
// ----------------------------------------------------------------------------

// TODO(singer): createRankDictionary
template <typename TValue, typename TText> 
inline void createRankDictionary(RankDictionary<SequenceBitMask<TValue> > & dictionary, TText const & text)
{ 
    typedef typename Fibre<RankDictionary<SequenceBitMask<TValue> >, FibreBitStrings>::Type TBitStrings;

    TBitStrings & bitStrings = getFibre(dictionary, FibreBitStrings());
    resize(bitStrings, ValueSize<TValue>::VALUE, Exact());

    for (unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i)
        resize(bitStrings[i], length(text), 0, Exact());
    
    for (unsigned i = 0; i < length(text); ++i)
        setBitTo(bitStrings[ordValue(text[i])], i, 1);

    for (unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i)
        _updateRanks(bitStrings[i]);
}

template <typename TValue, typename TSpec, typename TPrefixSumTable, typename TText> 
inline void createRankDictionary(LfTable<SentinelRankDictionary<RankDictionary<SequenceBitMask<TValue> >, TSpec >, TPrefixSumTable> & lfTable,
                              TText const & text)
{
    createRankDictionary(getFibre(getFibre(lfTable, FibreOccTable()), FibreRankDictionary()), text);
}

// ----------------------------------------------------------------------------
// Function open
// ----------------------------------------------------------------------------

template <typename TValue>
inline bool open(RankDictionary<SequenceBitMask<TValue> > & dictionary, const char * fileName, int openMode)
{
    String<char> name;
    name = fileName;    append(name, ".rd"); if (!open(getFibre(dictionary, FibreBitStrings()), toCString(name), openMode)) return false;
    return true;
}

template <typename TValue>
inline bool open(RankDictionary<SequenceBitMask<TValue> > & dictionary, const char * fileName)
{
    return open(dictionary, fileName, DefaultOpenMode<RankDictionary<SequenceBitMask<TValue> > >::VALUE);
}

// ----------------------------------------------------------------------------
// Function save
// ----------------------------------------------------------------------------

template <typename TValue>
inline bool save(RankDictionary<SequenceBitMask<TValue> > const & dictionary, const char * fileName, int openMode)
{
    String<char> name;
    name = fileName;    append(name, ".rd");   if (!save(getFibre(dictionary, FibreBitStrings()), toCString(name), openMode)) return false;
    
    return true;
}

template <typename TValue>
inline bool save(RankDictionary<SequenceBitMask<TValue> > const & dictionary, const char * fileName)
{
    return save(dictionary, fileName, DefaultOpenMode<RankDictionary<SequenceBitMask<TValue> > >::VALUE);
}
}
#endif  // INDEX_FM_RANKDICTIONARY_SBM
