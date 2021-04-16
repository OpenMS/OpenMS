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

// TODO(singer): sentinel

//SEQAN_NO_DDDOC:do not generate documentation for this file

#ifndef INDEX_FM_Sentinel_RANK_DICTIONARY_H_
#define INDEX_FM_Sentinel_RANK_DICTIONARY_H_

namespace seqan {

// ==========================================================================
//Forwards
// ==========================================================================

template <typename TValue>
class RankDictionary;

template<typename TRankDictionarySpec, typename TSpec> 
class SentinelRankDictionary;

// ==========================================================================
// Tags
// ==========================================================================

struct Sentinel_;
struct Sentinels_;
struct FibreRankDictionary_;
struct FibreSentinelPosition_;

typedef Tag<Sentinel_> const                Sentinel;
typedef Tag<Sentinels_> const               Sentinels;
typedef Tag<FibreRankDictionary_> const     FibreRankDictionary;
typedef Tag<FibreSentinelPosition_> const   FibreSentinelPosition;

// ----------------------------------------------------------------------------
// Spec SentinelRankDictionary Fibres
// ----------------------------------------------------------------------------

/**
.Tag.SentinelRankDictionary Fibres
..summary:Tag to select a specific fibre (e.g. table, object, ...) of a @Class.SentinelRankDictionary@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a SentinelRankDictionary.
..cat:SentinelRankDictionary

..tag.FibreRankDictionary:The rank dictionary.

..tag.FibreSentinelPosition:The bit string encoding the position of the sentinel sign.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index.h
*/

// ==========================================================================
// Metafunctions
// ==========================================================================

// ----------------------------------------------------------------------------
// Metafunction Fibre
// ----------------------------------------------------------------------------

// TODO(DOC)
///.Metafunction.Fibre.param.TSpec.type:Tag.SentinelRankDictionary Fibres

template <typename TValue, typename TSpec>
struct Fibre<SentinelRankDictionary<RankDictionary<WaveletTree<TValue> >, TSpec>, FibreRankDictionary>
{
    typedef RankDictionary<WaveletTree<TValue> > Type;
};

template <typename TValue, typename TSpec>
struct Fibre<SentinelRankDictionary<RankDictionary<SequenceBitMask<TValue> >, TSpec>, FibreRankDictionary>
{
    typedef RankDictionary<SequenceBitMask<TValue> > Type;
};

template <typename TRankDictionary, typename TSpec>
struct Fibre<SentinelRankDictionary<TRankDictionary, TSpec> const, FibreRankDictionary>
{
    typedef typename Fibre<SentinelRankDictionary<TRankDictionary, TSpec>, FibreRankDictionary>::Type const Type;
};

template <typename TRankDictionary>
struct Fibre<SentinelRankDictionary<TRankDictionary, Sentinel>, FibreSentinelPosition>
{
    typedef typename Size<TRankDictionary>::Type Type;
};

template <typename TRankDictionary>
struct Fibre<SentinelRankDictionary<TRankDictionary, Sentinel> const, FibreSentinelPosition>
{
    typedef typename Fibre<SentinelRankDictionary<TRankDictionary, Sentinel>, FibreSentinelPosition>::Type const Type;
};

template <typename TRankDictionary>
struct Fibre<SentinelRankDictionary<TRankDictionary, Sentinels>, FibreSentinelPosition>
{
    typedef RankSupportBitString<void> Type;
};

template <typename TRankDictionary>
struct Fibre<SentinelRankDictionary<TRankDictionary, Sentinels> const, FibreSentinelPosition>
{
    typedef typename Fibre<SentinelRankDictionary<TRankDictionary, Sentinels>, FibreSentinelPosition>::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TRankDictionary, typename TSpec>
struct Value<SentinelRankDictionary<TRankDictionary, TSpec> > :
    public Value<TRankDictionary> {};

template <typename TRankDictionary, typename TSpec>
struct Value<SentinelRankDictionary<TRankDictionary, TSpec> const> :    
    public Value<TRankDictionary const> {};

// ==========================================================================
// Classes
// ==========================================================================

// ----------------------------------------------------------------------------
// Class SentinelRankDictionary
// ----------------------------------------------------------------------------

/**
.Class.SentinelRankDictionary:
..cat:Index
..summary:A rank dictionary, additional storing sentinel character which are not accounted for in a rank querry.
..signature:SentinelRankDictionary<TRankDictionary, TSpec>
..param.TRankDictionary:The rank dictionary of a text.
...type:Class.RankDictionary
..param.TSpec:Specialisation
..include:seqan/index.h
*/

template <typename TLength, typename TTag>
TTag
_setDefaultSentinelPosition(TLength const _length, TTag const & /*tag*/)
{
    return _length;
}

template <typename TLength, typename TBitStringSpec>
RankSupportBitString<TBitStringSpec>
_setDefaultSentinelPosition(TLength const _length, RankSupportBitString<TBitStringSpec> const & /*tag*/)
{

    RankSupportBitString<TBitStringSpec> bitString;
    resize(bitString, _length, 0, Exact());
    return bitString;
}

template <typename TRankDictionary, typename TSpec>
class SentinelRankDictionary
{
public:
    typedef typename Value<SentinelRankDictionary>::Type TChar;
    typedef typename Fibre<SentinelRankDictionary<TRankDictionary, TSpec>, FibreSentinelPosition>::Type TSentinelPosition;

    TRankDictionary rankDictionary;
    TSentinelPosition sentinelPosition;
    TChar sentinelSubstitute;

    SentinelRankDictionary() :
        sentinelPosition(_setDefaultSentinelPosition(0u, TSentinelPosition())),
        sentinelSubstitute()
    {}
   
    // TODO(singer): Use concept sequence when available
    template <typename TValue, typename TStringSpec>
    SentinelRankDictionary(String<TValue, TStringSpec> const & text) :
        rankDictionary(text),
        sentinelPosition(_setDefaultSentinelPosition(length(text), TSentinelPosition())),
        sentinelSubstitute()
    {}
    
    template <typename THost, typename TSegmentSpec>
    SentinelRankDictionary(Segment<THost, TSegmentSpec> const & text) :
        rankDictionary(text),
        sentinelPosition(_setDefaultSentinelPosition(length(text), TSentinelPosition())),
        sentinelSubstitute()
    {}

    bool operator==(const SentinelRankDictionary & b) const
    {
        return rankDictionary == b.rankDictionary &&
               sentinelPosition == b.sentinelPosition &&
               sentinelSubstitute == b.sentinelSubstitute;
    }
};

//     // TODO(singer): Use concept sequence when available
//     template <typename TRankDictionary, typename TValue, typename TStringSpec>
//     SentinelRankDictionary<TRankDictionary, Sentinel>::SentinelRankDictionary(String<TValue, TStringSpec> const & text) :
//         rankDictionary(text),
//         sentinelPosition(length(text)),
//         sentinelSubstitute()
//     {}
//     
//     template <typename TRankDictionary, typename THost, typename TSegmentSpec>
//     SentinelRankDictionary::SentinelRankDictionary<TRankDictionary, Sentinel>(Segment<THost, TSegmentSpec> const & text) :
//         rankDictionary(text),
//         sentinelPosition(length(text)),
//         sentinelSubstitute()
//     {}
// 
//     template <typename TRankDictionary, typename TValue, typename TStringSpec>
//     SentinelRankDictionary::SentinelRankDictionary<TRankDictionary, Sentinels>(String<TValue, TStringSpec> const & text) :
//         rankDictionary(text),
//         sentinelPosition(),
//         sentinelSubstitute()
//     {}
//     
//     template <typename TRankDictionary, typename THost, typename TSegmentSpec>
//     SentinelRankDictionary::SentinelRankDictionary<TRankDictionary, Sentinels>(Segment<THost, TSegmentSpec> const & text) :
//         rankDictionary(text),
//         sentinelPosition(length(text)),
//         sentinelSubstitute()
//     {}

// ==========================================================================
//Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function clear
// ----------------------------------------------------------------------------

/**
.Function.SentinelRankDictionary#clear
..class:Class.SentinelRankDictionary
..summary:Clears the dictionary.
..signature:clear(dictionary)
..param.dictionary:The rank dictionary to be cleared.
...type:Class.SentinelRankDictionary
..include:seqan/index.h
*/

template <typename TRankDictionary>
inline void _clearSentinental(SentinelRankDictionary<TRankDictionary, Sentinel> &)
{}

template <typename TRankDictionary>
inline void _clearSentinental(SentinelRankDictionary<TRankDictionary, Sentinels> & dictionary)
{
    clear(getFibre(dictionary, FibreSentinelPosition()));
}

template <typename TRankDictionary, typename TSpec>
inline void clear(SentinelRankDictionary<TRankDictionary, TSpec> & dictionary)
{
    clear(getFibre(dictionary, FibreRankDictionary()));
    _clearSentinental(dictionary);
}

// ----------------------------------------------------------------------------
// Function sentinelPosition
// ----------------------------------------------------------------------------

/**
.Function.SentinelRankDictionary#sentinelPosition
..class:Class.SentinelRankDictionary
..summary:Returns whether a specified position is a sentinel position.
..signature:sentinelPosition(dictionary, pos)
..param.dictionary:The dictionary.
...type:Class.SentinelRankDictionary
..param.pos:The position.
...type:Concept.UnsignedIntegerConcept
..include:seqan/index.h
*/

//.Function.sentinelPosition.param.type:Class.RankDictionary
template <typename TRankDictionary, typename TPos>
inline bool sentinelPosition(SentinelRankDictionary<TRankDictionary, Sentinel> const & dictionary, TPos pos)
{
    return dictionary.sentinelPosition == pos;
}

template <typename TRankDictionary, typename TPos>
inline bool sentinelPosition(SentinelRankDictionary<TRankDictionary, Sentinels> const & dictionary, TPos pos)
{
    return isBitSet(getFibre(dictionary, FibreSentinelPosition()), pos);
}

// ----------------------------------------------------------------------------
// Function empty
// ----------------------------------------------------------------------------

/**
.Function.SentinelRankDictionary#empty
..class:Class.SentinelRankDictionary
..summary:Returns whether or not the dictionary is empty.
..signature:empty(dictionary)
..param.dictionary:The rank dictionary to be checked.
...type:Class.SentinelRankDictionary
..include:seqan/index.h
*/
template <typename TRankDictionary, typename TSpec>
inline bool empty(SentinelRankDictionary<TRankDictionary, TSpec> const & dictionary)
{
    return empty(getFibre(dictionary, FibreRankDictionary()));
}

// ----------------------------------------------------------------------------
// Function getValue
// ----------------------------------------------------------------------------

/**
.Function.SentinelRankDictionary#getValue
..class:Class.SentinelRankDictionary
..summary:Returns the character of a specified position.
..signature:getCharacter(dictionary, pos)
..param.dictionary:The rank dictionary.
...type:Class.RankDictionary
..param.pos:The position
..include:seqan/index.h
..example.code:
*/
template <typename TRankDictionary, typename TSpec, typename TPos>
inline typename Value<TRankDictionary>::Type
getValue(SentinelRankDictionary<TRankDictionary, TSpec > const & dictionary,
                 TPos pos)
{
    SEQAN_ASSERT_NEQ(sentinelPosition(dictionary, pos), true);
    return getValue(getFibre(dictionary, FibreRankDictionary()), pos);
}

template <typename TRankDictionary, typename TSpec, typename TPos>
inline typename Value<TRankDictionary>::Type
getValue(SentinelRankDictionary<TRankDictionary, TSpec > & dictionary,
                 TPos pos)
{
    SEQAN_ASSERT_NEQ(sentinelPosition(dictionary, pos), true);
    return getValue(getFibre(dictionary, FibreRankDictionary()), pos);
}

// ----------------------------------------------------------------------------
// Function getSentinelPosition
// ----------------------------------------------------------------------------

/*
.Function.SentinelRankDictionary#_getSentinelPosition
..class:Class.SentinelRankDictionary
..summary:Returns the position(s) of sentinel(s).
..signature:getSentinelPosition(dictionary)
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..remark:If there are several sentinels, a string containing the sentinel positions is returned.
 Should be used for debugging only.
..include:seqan/index.h
*/

template <typename TRankDictionary>
inline typename Size<SentinelRankDictionary<TRankDictionary, Sentinel> >::Type
_getSentinelPosition(SentinelRankDictionary<TRankDictionary, Sentinel> const & dictionary)
{
    return dictionary.sentinelPosition;
}

template <typename TRankDictionary>
inline String<typename Size<SentinelRankDictionary<TRankDictionary, Sentinel> >::Type > 
_getSentinelPosition(SentinelRankDictionary<TRankDictionary, Sentinels> const & dictionary)
{
    String<unsigned long> sentinelPositions;

    for (unsigned i = 0; i < length(dictionary.sentinelPosition); ++i)
        if (isBitSet(dictionary.sentinelPosition, i))
            appendValue(sentinelPositions, i);
    return sentinelPositions;
}

// ----------------------------------------------------------------------------
// Function getFibre
// ----------------------------------------------------------------------------

/**
.Function.SentinelRankDictionary#getFibre:
..class:Class.SentinelRankDictionary
..summary:Returns a specific fibre of a dictionary.
..signature:getFibre(dictionary, fibreTag)
..class:Class.RankDictionary
..cat:Index
..param.dictionary:The dictionary holding the fibre.
...type:Class.RankDictionary
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.SentinelRankDictionary Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
..include:seqan/index.h
*/
template <typename TRankDictionary, typename TSpec>
inline typename Fibre<SentinelRankDictionary<TRankDictionary, TSpec>, FibreRankDictionary>::Type &
getFibre(SentinelRankDictionary<TRankDictionary, TSpec>& dictionary, FibreRankDictionary)
{
    return dictionary.rankDictionary;
}

template <typename TRankDictionary, typename TSpec>
inline typename Fibre<SentinelRankDictionary<TRankDictionary, TSpec>, FibreRankDictionary>::Type const &
getFibre(SentinelRankDictionary<TRankDictionary, TSpec> const & dictionary, FibreRankDictionary)
{
    return dictionary.rankDictionary;
}

template <typename TRankDictionary, typename TSpec>
inline typename Fibre<SentinelRankDictionary<TRankDictionary, TSpec>, FibreSentinelPosition>::Type &
getFibre(SentinelRankDictionary<TRankDictionary, TSpec>& dictionary, FibreSentinelPosition)
{
    return dictionary.sentinelPosition;
}

template <typename TRankDictionary, typename TSpec>
inline typename Fibre<SentinelRankDictionary<TRankDictionary, TSpec>, FibreSentinelPosition>::Type const &
getFibre(SentinelRankDictionary<TRankDictionary, TSpec> const & dictionary, FibreSentinelPosition)
{
    return dictionary.sentinelPosition;
}

// ----------------------------------------------------------------------------
// Function countOccurrences
// ----------------------------------------------------------------------------

/**
.Function.SentinelRankDictionary#countOccurrences
..class:Class.SentinelRankDictionary
..summary:Returns the number of occurrences of a specified character from the start
to a specified position.
..signature:countOccurrences(dictionary, character, pos)
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..param.character:The character.
..param.pos:The position (which is included in the counting).
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";
RankDictionary<String<Dna5> > dictionary(genome);

std::cout << countOccurrences(dictionary, 'a', 3) << std::endl; // 1
std::cout << countOccurrences(dictionary, 'a', 4) << std::endl; // 2
*/

template <typename TRankDictionary, typename TChar, typename TPos>
inline unsigned countOccurrences(SentinelRankDictionary<TRankDictionary, Sentinel> const & dictionary,
                               TChar const character,
                               TPos const pos)
{
    unsigned occ = countOccurrences(getFibre(dictionary, FibreRankDictionary()), character, pos);
    if (ordEqual(getSentinelSubstitute(dictionary), character) && pos >= dictionary.sentinelPosition)
         --occ;
    return occ;
}

template <typename TRankDictionary, typename TChar, typename TPos>
inline unsigned countOccurrences(SentinelRankDictionary<TRankDictionary, Sentinels> const & dictionary,
                                TChar const character,
                                TPos const pos)
{
    unsigned occ = countOccurrences(getFibre(dictionary, FibreRankDictionary()), character, pos);
    if (ordEqual(getSentinelSubstitute(dictionary), character))
        return occ - getRank(getFibre(dictionary, FibreSentinelPosition()), pos);

    return occ;
}

// ----------------------------------------------------------------------------
// Function getSentinelSubstitute
// ----------------------------------------------------------------------------

/**
.Function.SentinelRankDictionary#getSentinelSubstitute
..class:Class.SentinelRankDictionary
..summary:Returns the character used to substitute the sentinel sign.
..signature:getSentinelSubstitute(dictionary)
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..include:seqan/index.h
*/
template <typename TRankDictionary, typename TSpec>
inline typename Value<TRankDictionary>::Type
getSentinelSubstitute(SentinelRankDictionary<TRankDictionary, TSpec> const & dictionary /*tag*/)
{
    return dictionary.sentinelSubstitute;
}

// ----------------------------------------------------------------------------
// Function setSentinelSubstitute
// ----------------------------------------------------------------------------

/**
.Function.SentinelRankDictionary#setSentinelSubstitute
..class:Class.SentinelRankDictionary
..summary:Sets the character used to substitute the sentinel sign.
..signature:setSentinelSubstitute(dictionary, character)
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..param.character:The sentinel substitute.
..include:seqan/index.h
*/

template <typename TRankDictionary, typename TSpec, typename TChar>
inline void setSentinelSubstitute(SentinelRankDictionary<TRankDictionary, TSpec> & dictionary,
                                TChar sentinelSubstitute)
{
    dictionary.sentinelSubstitute = sentinelSubstitute;
}

// ----------------------------------------------------------------------------
// Function setSentinelPosition
// ----------------------------------------------------------------------------

/**
.Function.SentinelRankDictionary#setSentinelPosition
..class:Class.SentinelRankDictionary
..summary:Sets the sentinel position..
..signature:setSentinelPosition(dictionary, pos)
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..param.pos:The sentinel position.
..include:seqan/index.h
*/

template <typename TRankDictionary, typename TSpec, typename TPos>
inline void setSentinelPosition(SentinelRankDictionary<TRankDictionary, TSpec> & dictionary,
                              TPos const & position)
{
    dictionary.sentinelPosition = position;
}

// ----------------------------------------------------------------------------
// Function sentinelRankDictionaryCreate
// ----------------------------------------------------------------------------

/**
.Function.SentinelRankDictionary#createSentinelRankDictionary
..class:Class.SentinelRankDictionary
..summary:This functions creates the dictionary structure.
..signature:void createSentinelRankDictionary(dictionary, text)
..param.dictionary:The dictionary.
...type:Class.RankDictionary.
..param.text:A text to be transfered into a dictionary.
...type:Class.String
..include:seqan/index.h
*/

template <typename TRankDictionary, typename TSpec, typename TText, typename TSentinelSub, typename TSentinelPos>
inline void createSentinelRankDictionary(SentinelRankDictionary<TRankDictionary, TSpec> & dictionary,
                              TText const & text,
                              TSentinelSub const & sentinelSub,
                              TSentinelPos const & sentinelPos)
{
    setSentinelSubstitute(dictionary, sentinelSub);
    setSentinelPosition(dictionary, sentinelPos);

    createRankDictionary(getFibre(dictionary, FibreRankDictionary()), text);
}

template <typename TRankDictionary, typename TSpec, typename TPrefixSumTable, typename TText, typename TDollarSub, typename TDollarPos>
inline void createSentinelRankDictionary(LfTable<SentinelRankDictionary<TRankDictionary, TSpec>, TPrefixSumTable> & lfTable,
                              TText const & text,
                              TDollarSub const & dollarSub,
                              TDollarPos const & dollarPos)
{
    setSentinelSubstitute(getFibre(lfTable, FibreOccTable()), dollarSub);
    setSentinelPosition(getFibre(lfTable, FibreOccTable()), dollarPos);

    createRankDictionary(lfTable, text);
}

// ----------------------------------------------------------------------------
// Function open
// ----------------------------------------------------------------------------

/**
.Function.SentinelRankDictionary#open
..class:Class.SentinelRankDictionary
..summary:This functions loads a dictionary from disk.
..signature:open(dictionary, fileName [, openMode])
..param.dictionary:The dictionary.
...type:Class.SentinelRankDictionary
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

template <typename TRankDictionary>
inline bool _openSentinelInformation(
    SentinelRankDictionary<TRankDictionary, Sentinel> & dictionary,
    const char * fileName,
    int openMode)
{
    String<char> name;

    typedef typename Fibre<SentinelRankDictionary<TRankDictionary, Sentinel>, FibreSentinelPosition>::Type TSentinelString;
    typedef typename Value<SentinelRankDictionary<TRankDictionary, Sentinel> >::Type TChar;

    String<Pair<TChar, TSentinelString, Pack> > sentinelValues;

    name = fileName;    append(name, ".dr");
    if (!open(sentinelValues, toCString(name), openMode) || empty(sentinelValues))
    {
        return false;
    }
    dictionary.sentinelSubstitute = sentinelValues[0].i1;
    dictionary.sentinelPosition = sentinelValues[0].i2;
    return true;
}

template <typename TRankDictionary>
inline bool _openSentinelInformation(
    SentinelRankDictionary<TRankDictionary, Sentinels> & dictionary,
    const char * fileName,
    int openMode)
{
    String<char> name;

    typedef typename Value<SentinelRankDictionary<TRankDictionary, Sentinels> >::Type TChar;
    String<TChar> sentinelSub;
    
    name = fileName;    append(name, ".drs"); if (!open(sentinelSub, toCString(name), openMode)) return false;
    name = fileName;    append(name, ".drp"); if (!open(dictionary.sentinelPosition, toCString(name), openMode)) return false;
    
    if (empty(sentinelSub))
        return false;

    dictionary.sentinelSubstitute = sentinelSub[0];
    
    return true;
}

template <typename TRankDictionary, typename TSpec>
inline bool open(
    SentinelRankDictionary<TRankDictionary, TSpec> & dictionary,
    const char * fileName,
    int openMode)
{
    if (!open(getFibre(dictionary, FibreRankDictionary()), fileName, openMode)) return false;
    if (!_openSentinelInformation(dictionary, fileName, openMode)) return false;
    return true;
}


    template <typename TRankDictionary, typename TSpec>
inline bool open(
    SentinelRankDictionary<TRankDictionary, TSpec> & dictionary,
    const char * fileName)
{
    return open(dictionary, fileName, DefaultOpenMode<SentinelRankDictionary<TRankDictionary, Sentinel> >::VALUE);
}

// ----------------------------------------------------------------------------
// Function save
// ----------------------------------------------------------------------------

/**
.Function.SentinelRankDictionary#save
..class:Class.SentinelRankDictionary
..summary:This functions saves a dictionary to disk.
..signature:open(dictionary, fileName [, openMode])
..param.dictionary:The dictionary.
...type:Class.SentinelRankDictionary
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

template <typename TRankDictionary>
inline bool _saveSentinelInformation(
    SentinelRankDictionary<TRankDictionary, Sentinel> const & dictionary,
    const char * fileName,
    int openMode)
{
    String<char> name;

    typedef typename Value<SentinelRankDictionary<TRankDictionary, Sentinel> >::Type TChar;
    typedef typename Fibre<SentinelRankDictionary<TRankDictionary, Sentinel>, FibreSentinelPosition>::Type TSentinelString;

    String<Pair<TChar, TSentinelString, Pack> > sentinelValues;
    appendValue(sentinelValues, Pair<TChar, TSentinelString>(dictionary.sentinelSubstitute, dictionary.sentinelPosition));
    
    name = fileName;    append(name, ".dr"); if (!save(sentinelValues, toCString(name), openMode)) return false;
    return true;
}

template <typename TRankDictionary>
inline bool _saveSentinelInformation(
    SentinelRankDictionary<TRankDictionary, Sentinels> const & dictionary,
    const char * fileName,
    int openMode)
{
    String<char> name;

    typedef typename Value<SentinelRankDictionary<TRankDictionary, Sentinels> >::Type TChar;
    String<TChar> sentinelSub;
    appendValue(sentinelSub, dictionary.sentinelSubstitute);

    name = fileName;    append(name, ".drs"); if(!save(sentinelSub, toCString(name), openMode)) return false;
    name = fileName;    append(name, ".drp"); if(!save(dictionary.sentinelPosition, toCString(name), openMode)) return false;
    return true;
}

template <typename TRankDictionary, typename TSpec>
inline bool save(
    SentinelRankDictionary<TRankDictionary, TSpec> const & dictionary,
    const char * fileName,
    int openMode)
{
    if (!save(getFibre(dictionary, FibreRankDictionary()), fileName, openMode)) return false;
    if (!_saveSentinelInformation(dictionary, fileName, openMode)) return false;
    
    return true;
}

template <typename TRankDictionary, typename TSpec>
inline bool save(
    SentinelRankDictionary<TRankDictionary, TSpec> const & dictionary,
    const char * fileName)
{
    return save(dictionary, fileName, DefaultOpenMode<SentinelRankDictionary<TRankDictionary, Sentinel> >::VALUE);
}

template <typename TRankDictionary, typename TSpec>
inline bool save(
    SentinelRankDictionary<TRankDictionary, TSpec> & dictionary,
    const char * fileName)
{
    return save(dictionary, fileName, DefaultOpenMode<SentinelRankDictionary<TRankDictionary, Sentinel> >::VALUE);
}

}
#endif  // INDEX_FM_Sentinel_RANK_DICTIONARY_H_

