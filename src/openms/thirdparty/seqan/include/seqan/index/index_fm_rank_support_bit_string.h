// ==========================================================================
//                 seqan - the library for sequence analysis
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

#ifndef INDEX_FM_RANK_SUPPORT_BIT_STRING_H_
#define INDEX_FM_RANK_SUPPORT_BIT_STRING_H_

#include <seqan/misc/misc_bit_twiddling.h>

namespace seqan {

//Rank Support Bit String
template <typename TSpec = void>
struct RankSupportBitString;

// ==========================================================================
// Tags
// ==========================================================================

// FM index fibres

/**
.Tag.RankSupportBitString Fibres
..summary:Tag to select a specific fibre (e.g. table, object, ...) of a @Class.RankSupportBitString@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a rank support bit string.
..cat:Index

..tag.FibreBits:The bit string.
..tag.FibreBlocks:The block string.
..tag.FibreSuperBlocks:The super block string.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index.h
*/

struct FibreBits_;
struct FibreBlocks_;
struct FibreSuperBlocks_;

typedef Tag<FibreBits_> const           FibreBits;         
typedef Tag<FibreBlocks_> const         FibreBlocks;      
typedef Tag<FibreSuperBlocks_> const    FibreSuperBlocks; 

typedef FibreBits           RankSupportBitStringBits;
typedef FibreBlocks         RankSupportBitStringBlocks;
typedef FibreSuperBlocks    RankSupportBitStringSuperBlocks;

// ==========================================================================
// Metafunctions
// ==========================================================================

template <typename TSpec>
struct DefaultOverflowImplicit<RankSupportBitString<TSpec> >
{
    typedef Generous Type;
};

//The limiting factor of the size is the underlying data type of the super block string
template <typename TSpec>
struct Size<RankSupportBitString<TSpec> >
{
    typedef unsigned long long Type;
};

template <typename TSpec>
struct Position<RankSupportBitString<TSpec> > :
	Size<RankSupportBitString<TSpec> > {};

template <typename TSpec>
struct Fibre<RankSupportBitString<TSpec>, FibreBits>
{
    typedef String<unsigned long long> Type;
};

template <typename TSpec>
struct Fibre<RankSupportBitString<TSpec>, FibreBlocks>
{
    typedef String<unsigned short> Type;
};

template <typename TSpec>
struct Fibre<RankSupportBitString<TSpec>, FibreSuperBlocks>
{
    typedef String<typename Size<RankSupportBitString<TSpec> >::Type> Type;
};

// ==========================================================================
// Classes
// ==========================================================================

/**
.Class.RankSupportBitString:
..summary:A bit string supporting rank queries in constant time.
..cat:Index
..signature:RankSupportBitString<TSpec>
..param.TSpec:Specialisation tag.
...default:void
..remarks:The constant rank query time is achieved by evaluating precomputed subsolutions. In order to do so, the bit string is divided into blocks of length l. A super block string stores for each block of l blocks the number of bits set from the beginning. In addition a block string stores the number of bits set in each block from the start of the last super block block. Therefore it is possible to compute the result of a rank query in constant time by adding information from the bit, block and super block string.
..include:seqan/index.h
*/

// Forward declaration because setBit is used in constructor
template <typename TSpec, typename TPos>
inline void setBitTo(RankSupportBitString<TSpec> & bitString, TPos pos, bool setBit);

template <typename TSpec>
struct RankSupportBitString
{
    typedef typename Fibre<RankSupportBitString, FibreBits>::Type         TBitString;
    typedef typename Fibre<RankSupportBitString, FibreBlocks>::Type       TBlockString;
    typedef typename Fibre<RankSupportBitString, FibreSuperBlocks>::Type  TSuperBlockString;

    TBitString                                  bits;
    TBlockString                                blocks;
    TSuperBlockString                           superBlocks;
    typename Size<RankSupportBitString>::Type   _length;

    RankSupportBitString() :
        _length(0)
    {}

    // TODO(singer): Use concept sequence when available
    template <typename TValue, typename TStringSpec>
    RankSupportBitString(String<TValue, TStringSpec> const & input) :
        _length(length(input))
    {
        typedef typename Iterator<String<TValue, TStringSpec> const>::Type TIter_;

        resize(*this, _length, Exact());
        typename Position<RankSupportBitString>::Type i = 0;
        for (TIter_ it = begin(input, Standard()); it != end(input, Standard()); ++it, ++i)
            setBitTo(*this, i, getValue(it));
        _updateRanks(*this);
    }
    
    template <typename THost, typename TStringSpec>
    RankSupportBitString(Segment<THost, TStringSpec> const & input) :
        _length(length(input))
    {
        typedef typename Iterator<Segment<THost, TStringSpec> const>::Type TIter_;

        resize(*this, _length, Exact());
        typename Position<RankSupportBitString>::Type i = 0;
        for (TIter_ it = begin(input, Standard()); it != end(input, Standard()); ++it, ++i)
            setBitTo(*this, i, getValue(it));
        _updateRanks(*this);
    }

    inline bool operator==(const RankSupportBitString & other) const
    {
        return _length == other._length &&
               bits == other.bits &&
               blocks == other.blocks &&
               superBlocks == other.superBlocks;
    }
};

// ==========================================================================
// Functions
// ==========================================================================

/**
.Function.RankSupportBitString#appendValue:
..signature:appendValue(target, value)
..cat:Index
..summary:Appends a value to a container.
..param.target:A container.
...type:Class.RankSupportBitString
..param.value:Value that is appended to $target$.
...type:Concept.UnsignedIntegerConcept
...type:nolink:bool
...remarks:If the value is different from 0 it is interpreted as 1.
..include:seqan/sequence.h
*/
template <typename TSpec>
inline void appendValue(RankSupportBitString<TSpec> & bitString, bool bit)
{
	typename Size<RankSupportBitString<TSpec> >::Type len = length(bitString);
    resize(bitString, len + 1);
    setBitTo(bitString, len, bit);
    _updateLastRank(bitString); 
}

// ==========================================================================
/**
.Function.clear
..param.object:
...type:Class.RankSupportBitString
*/
template <typename TSpec>
inline void clear(RankSupportBitString<TSpec> & bitString)
{
    clear(bitString.bits);
    clear(bitString.blocks);
    clear(bitString.superBlocks);
    bitString._length = 0;
}

// ==========================================================================
// This function returns the number of bits set in a block until a specified position
template <typename TValue>
inline unsigned _getRankInBlock(TValue const value, False)
{
    return popCount(value);
}

template <typename TValue>
inline unsigned _getRankInBlock(TValue const value, True)
{
    return popCount(value);
}

template <typename TValue>
inline unsigned _getRankInBlock(TValue const value)
{
    return _getRankInBlock(value, typename Eval<(BitsPerValue<TValue>::VALUE > 32)>::Type());
}

template <typename TSpec, typename TPos>
inline typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreSuperBlocks>::Type>::Type
_getRankInBlock(RankSupportBitString<TSpec> const & bitString, TPos const pos)
{
    //typedef RankSupportBitString<TSpec>                                     TRankSupportBitString;
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreBits>::Type    TFibreBits;
    typedef typename Value<TFibreBits>::Type                                TFibreBitsValue;

    TFibreBitsValue const mask = ((TFibreBitsValue)2u << _getPosInBlock(bitString, pos)) - 1;
    return _getRankInBlock(bitString.bits[_getBlockPos(bitString, pos)] & mask);
}

// ==========================================================================
// This function returns the number of bits set in a block until a specified position
template <typename TSpec, typename TPos>
inline typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreSuperBlocks>::Type>::Type
_getSuperBlockPos(RankSupportBitString<TSpec> const & /*bitString*/, TPos const pos)
{
    typedef RankSupportBitString<TSpec>                                 TRankSupportBitString;
    typedef typename Fibre<TRankSupportBitString, FibreBits>::Type      TFibreBits;
    typedef typename Fibre<TRankSupportBitString, FibreBlocks>::Type    TFibrelocks;
    typedef typename Value<TFibreBits>::Type                            TFibreBitsValue;
    typedef typename Value<TFibrelocks>::Type                     		TFibreBlocksValue;

    TFibreBlocksValue const _bitsPerValue = BitsPerValue<TFibreBitsValue>::VALUE;
    return pos / (_bitsPerValue * _bitsPerValue);
}

// ==========================================================================
/**
.Function.getRank
..summary:Returns the rank (the number of bits set from the start of the bit string) of a specified position.
..signature:getRank(bitString, pos)
..param.bitString:The bit string.
...type:Class.RankSupportBitString
..param.pos:Position of a bit.
..returns:Value type of the super block fibre (default unsigned long).
..include:seqan/index.h
..example.code:
*/
template <typename TSpec, typename TPos>
inline typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreSuperBlocks>::Type>::Type
getRank(RankSupportBitString<TSpec> const & bitString, TPos const pos)
{
    return bitString.superBlocks[_getSuperBlockPos(bitString, pos)]
         + bitString.blocks[_getBlockPos(bitString, pos)]
         + _getRankInBlock(bitString, pos);
}


/**
.Function.empty
..param.object:
...type:Class.RankSupportBitString
*/

// ==========================================================================
/**
.Function.isSetBit
..summary:Returns whether the bit with the given index is set to 1.
..signature:isSetBit(bitString, pos)
..param.bitString:The bit string.
...type:Class.RankSupportBitString
..param.pos:Position of the bit.
..returns:Returns whether a specified bit is set or not.
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";

RankSupportBitString<> bitString;
resize(bitString, length(genome));
...
mark all 'a's
...

for (unsigned i = 0; i < length(bitString); ++i)
    if(isSetBit(bitString, i))
        std::cout << "a found at: " << i << std::endl;
*/

template <typename TSpec, typename TPos>
inline bool isBitSet(RankSupportBitString<TSpec> const & bitString, TPos const pos)
{
    return isBitSet(bitString.bits[_getBlockPos(bitString, pos)], _getPosInBlock(bitString, pos));
}

// ==========================================================================
// This function returns the position in the block string of the corresponding block.
template <typename TSpec, typename TPos>
inline typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreSuperBlocks>::Type>::Type
_getBlockPos(RankSupportBitString<TSpec> const & /*bitString*/, TPos const pos)
{
    typedef RankSupportBitString<TSpec>                             TRankSupportBitString;
    typedef typename Fibre<TRankSupportBitString, FibreBits>::Type  TFibreBits;
    typedef typename Value<TFibreBits>::Type                        TFibreBitsValue;

    return pos / BitsPerValue<TFibreBitsValue>::VALUE;
}

// ==========================================================================
/**
.Function.RankSupportBitString#getFibre:
..summary:Returns a specific fibre of a container.
..signature:getFibre(container, fibreTag)
..class:Class.RankSupportBitString
..cat:Index
..param.container:The container holding the fibre.
...type:Class.RankSupportBitString
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.RankSupportBitString Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
..include:seqan/index.h
..example.code:
Index< String<char> > index_esa("tobeornottobe");

String<char> & text = getFibre(indexEsa, EsaText());
*/
template <typename TSpec>
inline typename Fibre<RankSupportBitString<TSpec>, FibreBits>::Type &
getFibre(RankSupportBitString<TSpec>&string, const FibreBits)
{
    return string.bits;
}

template <typename TSpec>
inline typename Fibre<RankSupportBitString<TSpec>, FibreBits>::Type const &
getFibre(RankSupportBitString<TSpec> const & string, const FibreBits)
{
    return string.bits;
}

template <typename TSpec>
inline typename Fibre<RankSupportBitString<TSpec>, FibreBlocks>::Type &
getFibre(RankSupportBitString<TSpec>&string, const FibreBlocks)
{
    return string.blocks;
}

template <typename TSpec>
inline typename Fibre<RankSupportBitString<TSpec>, FibreBlocks>::Type const &
getFibre(RankSupportBitString<TSpec> const & string, const FibreBlocks)
{
    return string.blocks;
}

template <typename TSpec>
inline typename Fibre<RankSupportBitString<TSpec>, FibreSuperBlocks>::Type &
getFibre(RankSupportBitString<TSpec>&string, const FibreSuperBlocks)
{
    return string.superBlocks;
}

template <typename TSpec>
inline typename Fibre<RankSupportBitString<TSpec>, FibreSuperBlocks>::Type const &
getFibre(const RankSupportBitString<TSpec>&string, const FibreSuperBlocks)
{
    return string.superBlocks;
}

// ==========================================================================
// This function returns the position of a specified bit within a block.
template <typename TSpec, typename TPos>
inline typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreBlocks>::Type>::Type
_getPosInBlock(RankSupportBitString<TSpec> const & /*bitString*/, TPos const pos)
{
    typedef RankSupportBitString<TSpec>                             TRankSupportBitString;
    typedef typename Fibre<TRankSupportBitString, FibreBits>::Type  TFibreBits;
    typedef typename Value<TFibreBits>::Type                        TFibreBitsValue;

    return pos % BitsPerValue<TFibreBitsValue>::VALUE;
}

// ==========================================================================
/**
.Function.length
..param.object:
...type:Class.RankSupportBitString
*/
template <typename TSpec>
inline typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreSuperBlocks>::Type>::Type
length(RankSupportBitString<TSpec> const & bitString)
{
    return bitString._length;
}


// ==========================================================================
/*
.Function._updateRanks
..summary:Adds the block and super block information to the bit string.
..signature:_updateRanks(bitString)
..param.bitString:The bit string to be completed.
...type:Class.RankSupportBitString
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";

RankSupportBitString<> bitString;
resize(bitString, length(genome));

for (unsigned i = 0; i < length(genome); ++i)
    if(genome[i] < Dna5('c'))
        setBitTo(bitString, 1);

_updateRanks(bitString);
*/

template <typename TSpec, typename TPos>
inline void _updateRanksImpl(RankSupportBitString<TSpec> & bitString, TPos pos)
{
    if (!empty(bitString))
    {
        typedef RankSupportBitString<TSpec>                                     TRankSupportBitString;
        typedef typename Fibre<RankSupportBitString<TSpec>, FibreBits>::Type    TFibreBits;
        typedef typename Fibre<TRankSupportBitString, FibreBlocks>::Type        TFibreBlocks;
        typedef typename Fibre<TRankSupportBitString, FibreSuperBlocks>::Type   TFibreSuperBlocks;
        typedef typename Value<TFibreBits>::Type                                TFibreBitsValue;
        typedef typename Value<TFibreBlocks>::Type                              TFibreBlocksValue;
        typedef typename Value<TFibreSuperBlocks>::Type                         TFibreSuperBlocksValue;

        TFibreSuperBlocksValue i = _getBlockPos(bitString, pos);

        TFibreBlocksValue _blockSum;
        if (i == 0)
            _blockSum = 0;
        else
            _blockSum = bitString.blocks[i - 1];

        TFibreSuperBlocksValue _sBlockSum;
        TFibreSuperBlocksValue superBlockPos = _getSuperBlockPos(bitString, pos);
        if (superBlockPos != 0)
        {
            _sBlockSum = bitString.superBlocks[superBlockPos];
        }
        else
            _sBlockSum = 0;
        
        if (i == 0)
            ++i;
        for (; i < length(bitString.bits); ++i)
        {
        	_blockSum += _getRankInBlock(bitString.bits[i - 1]);
            if ((i % BitsPerValue<TFibreBitsValue>::VALUE) == 0u)
            {
                _sBlockSum += _blockSum;
                bitString.superBlocks[++superBlockPos] = _sBlockSum;
                _blockSum = 0;
            }
            bitString.blocks[i] = _blockSum;
       }
    }
}

template <typename TSpec, typename TPos>
inline void _updateRanks(RankSupportBitString<TSpec> & bitString, TPos pos)
{
    _updateRanksImpl(bitString, pos);
}

template <typename TSpec>
inline void _updateRanks(RankSupportBitString<TSpec> & bitString)
{
    _updateRanksImpl(bitString, 0u);
}

// This function update the rank information of the last block.
// ==========================================================================
template <typename TSpec>
inline void _updateLastRank(RankSupportBitString<TSpec> & bitString)
{
    typedef RankSupportBitString<TSpec>                                     TRankSupportBitString;
    typedef typename Fibre<TRankSupportBitString, FibreSuperBlocks>::Type   TFibreSuperBlocks;
    typedef typename Value<TFibreSuperBlocks>::Type                         TFibreSuperBlocksValue;

    // It is only necessary to update the rank at the start of a new block
    TFibreSuperBlocksValue pos = length(bitString) - 1;
    if (_getPosInBlock(bitString, pos) != 0u) return;

    // Here we compute the new rank if the last bit does not reside in the 
    // first block and the last as well as the second last bit are in the 
    // same block.
    TFibreSuperBlocksValue superBlockPos = _getSuperBlockPos(bitString, pos);
    TFibreSuperBlocksValue blockPos = _getBlockPos(bitString, pos);
    if (blockPos > 0u && superBlockPos == _getSuperBlockPos(bitString, pos - 1))
        getFibre(bitString, FibreBlocks())[blockPos] = getFibre(bitString, FibreBlocks())[blockPos - 1]
                + _getRankInBlock(bitString.bits[blockPos - 1]);

    // Here we compute the new rank if the last and second last bit are in
    // different blocks.
    if ((superBlockPos > _getSuperBlockPos(bitString, pos - 1)) && superBlockPos > 0u)
        getFibre(bitString, FibreSuperBlocks())[superBlockPos] = getFibre(bitString, FibreSuperBlocks())[superBlockPos - 1]
                + getFibre(bitString, FibreBlocks())[blockPos - 1]
                + _getRankInBlock(bitString.bits[blockPos - 1]);
}


// ==========================================================================
template <typename TSpec, typename TSize, typename TExpand>
inline typename Size<RankSupportBitString<TSpec> >::Type
reserve(RankSupportBitString<TSpec> & bitString, TSize const size, Tag<TExpand> const tag)
{
    typedef RankSupportBitString<TSpec>                                     TRankSupportBitString;
    typedef typename Fibre<TRankSupportBitString, FibreBits>::Type          TFibreBits;
    typedef typename Fibre<TRankSupportBitString, FibreBlocks>::Type        TFibreBlocks;
    typedef typename Fibre<TRankSupportBitString, FibreSuperBlocks>::Type   TFibreSuperBlocks;
    typedef typename Value<TFibreBits>::Type                                TFibreBitsValue;
    typedef typename Value<TFibreBlocks>::Type                              TFibreBlocksValue;
    typedef typename Value<TFibreSuperBlocks>::Type                         TFibreSuperBlocksValue;

    TFibreBlocksValue _bitsPerBlock = BitsPerValue<TFibreBitsValue>::VALUE;
    TFibreSuperBlocksValue _numberOfBlocks = (size + _bitsPerBlock - 1) / _bitsPerBlock;

    resserve(bitString.blocks, _numberOfBlocks, tag);
    reserve(bitString.superBlocks, (_numberOfBlocks + _bitsPerBlock - 1) / _bitsPerBlock, tag);
    return reserve(bitString.bits, _numberOfBlocks, tag) * _bitsPerBlock;
}

// ==========================================================================
/**
.Function.resize
..param.object:
...type:Class.RankSupportBitString
*/
template <typename TSpec, typename TLength, typename TValue, typename TExpand>
inline typename Size<RankSupportBitString<TSpec> >::Type
resize(RankSupportBitString<TSpec> & bitString, TLength const _length, TValue const _value, Tag<TExpand> const tag)
{
    typedef RankSupportBitString<TSpec>                                     TRankSupportBitString;
    typedef typename Fibre<TRankSupportBitString, FibreBits>::Type          TFibreBits;
    typedef typename Fibre<TRankSupportBitString, FibreBlocks>::Type        TFibreBlocks;
    typedef typename Fibre<TRankSupportBitString, FibreSuperBlocks>::Type   TFibreSuperBlocks;
    typedef typename Value<TFibreBits>::Type                                TFibreBitsValue;
    typedef typename Value<TFibreBlocks>::Type                              TFibreBlocksValue;
    typedef typename Value<TFibreSuperBlocks>::Type                         TFibreSuperBlocksValue;

    TFibreBlocksValue _bitsPerBlock = BitsPerValue<TFibreBitsValue>::VALUE;
    TFibreSuperBlocksValue _numberOfBlocks = (_length + _bitsPerBlock - 1) / _bitsPerBlock;

    TLength currentLength = length(bitString);

    resize(bitString.bits, _numberOfBlocks, 0, tag);
    resize(bitString.blocks, _numberOfBlocks, 0, tag);
    resize(bitString.superBlocks, (_numberOfBlocks + _bitsPerBlock - 1) / _bitsPerBlock, 0, tag);

    if (_value != 0 && currentLength < _length)
    {
        for (unsigned i = currentLength; i < (typename MakeUnsigned<TLength const>::Type)_length; ++i)
            setBit(bitString, i);
        _updateRanks(bitString, currentLength);
    }

    return bitString._length = _length;
}

template <typename TSpec, typename TLength, typename TExpand>
inline typename Size<RankSupportBitString<TSpec> >::Type
resize(RankSupportBitString<TSpec> & bitString, TLength const _length, Tag<TExpand> const tag)
{
    return resize(bitString, _length, 0, tag);
}

// ==========================================================================
/**
.Function.setBitTo
..signature:setBitTo(bitString, pos, bit)
..param.bitString:The bit string.
...type:Class.RankSupportBitString
..param.pos:Position of the bit.
..param.bit:The value of the bit.
...remarks:Note that values different from 0 are interpreted as 1.
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";

RankSupportBitString<> bitString;
resize(bitString, length(genome));

for (unsigned i = 0; i < length(genome); ++i)
    if(genome[i] < Dna5('c'))
        setBitTo(bitString, 1);

updateRanks_(bitString);
*/
template <typename TSpec, typename TPos>
inline void setBit(RankSupportBitString<TSpec> & bitString, TPos pos)
{
    //typedef RankSupportBitString<TSpec>                                     TRankSupportBitString;
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreBits>::Type    TFibreBits;
    typedef typename Value<TFibreBits>::Type                                TFibreBitsValue;

    TFibreBitsValue const shiftValue = (TFibreBitsValue)1u << _getPosInBlock(bitString, pos);
    bitString.bits[_getBlockPos(bitString, pos)] |= shiftValue;
}

template <typename TSpec, typename TPos>
inline void clearBit(RankSupportBitString<TSpec> & bitString, TPos pos)
{
    //typedef RankSupportBitString<TSpec>                                     TRankSupportBitString;
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreBits>::Type    TFibreBits;
    typedef typename Value<TFibreBits>::Type                                TFibreBitsValue;

    TFibreBitsValue const shiftValue = (TFibreBitsValue)1u << _getPosInBlock(bitString, pos);
    bitString.bits[_getBlockPos(bitString, pos)] &= ~shiftValue;
}

template <typename TSpec, typename TPos>
inline void setBitTo(RankSupportBitString<TSpec> & bitString, TPos pos, bool value)
{
    if (value)
        setBit(bitString, pos);
    else
        clearBit(bitString, pos);
}

// ==========================================================================
/**
.Function.open
..signature:open(object, filename[ ,openMode])
..param.object:The object to be opened.
...type:Class.RankSupportBitString
*/

template <typename TSpec>
inline bool open(
    RankSupportBitString<TSpec> & string,
    const char * fileName,
    int openMode)
{
    typedef typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreSuperBlocks>::Type>::Type TValue;
    String<TValue> lengthString;
    resize(lengthString, 1, Exact());

    String<char> name;
    name = fileName;    append(name, ".bit");
    if (!open(getFibre(string, FibreBits()), toCString(name), openMode))
        return false;

    name = fileName;    append(name, ".bl");    open(getFibre(string, FibreBlocks()), toCString(name), openMode);
    name = fileName;    append(name, ".sbl");   open(getFibre(string, FibreSuperBlocks()), toCString(name), openMode);
    name = fileName;    append(name, ".len");    open(lengthString, toCString(name), openMode);
    string._length = lengthString[0];
    return true;
}

template <typename TSpec>
inline bool open(
    RankSupportBitString<TSpec> & string,
    const char * fileName)
{
    return open(string, fileName, OPEN_RDONLY);
}

template <typename TSpec, typename TSetSpec>
inline bool open(
    StringSet<RankSupportBitString<TSpec>, TSetSpec> & strings,
    const char * fileName)
{
    return open(strings, fileName, OPEN_RDONLY);
}

// ==========================================================================
/**
.Function.RankSupportBitString#save
..class:Class.RankSupportBitString
..summary:This functions saves a @Class.RankSupportBitString@ to disk.
..signature:save(bitString, fileName [, openMode])
..param.bitString:The bit string to be saved.
...type:Class.RightArrayBinaryTree
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
template <typename TSpec>
inline bool save(
    RankSupportBitString<TSpec> const & string,
    const char * fileName,
    int openMode)
{
    typedef typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreSuperBlocks>::Type>::Type TValue;
    String<TValue> lengthString;
    resize(lengthString, 1, Exact());
    lengthString[0] = length(string);
    String<char> name;
    name = fileName;    append(name, ".len");   save(lengthString, toCString(name), openMode);
    name = fileName;    append(name, ".bit");   save(getFibre(string, FibreBits()), toCString(name), openMode);
    name = fileName;    append(name, ".bl");    save(getFibre(string, FibreBlocks()), toCString(name), openMode);
    name = fileName;    append(name, ".sbl");   save(getFibre(string, FibreSuperBlocks()), toCString(name), openMode);
    return true;
}

template <typename TSpec>
inline bool save(
    RankSupportBitString<TSpec> const & string,
    const char * fileName)
{
    return save(string, fileName, DefaultOpenMode<RankSupportBitString<TSpec> >::VALUE);
}

template <typename TSpec, typename TSetSpec>
inline bool save(
    StringSet<RankSupportBitString<TSpec>, TSetSpec> const & strings,
    const char * fileName)
{
    return save(strings, fileName, DefaultOpenMode<RankSupportBitString<TSpec> >::VALUE);
}

// ==========================================================================
template <typename TValue>
inline void printBits(TValue entrie)
{
    unsigned bitsPerValue = BitsPerValue<TValue>::VALUE;
    TValue one = 1;
    std::cout << "entrie: " << entrie << std::endl;
    std::cout << bitsPerValue << std::endl;
    for (TValue i = 0; i < bitsPerValue; ++i)
    {
        std::cout << ((entrie >> i) & one);
    }
    std::cout << std::endl;
}

template <typename TValue, typename TSize>
inline std::ostream & printBits(std::ostream & stream, TValue entrie, TSize blockSize)
{
    unsigned bitsPerValue = BitsPerValue<TValue>::VALUE;
    bool temp;
    for (int i = bitsPerValue - 1; i >= 0; --i)
    {
        temp = (entrie >> i) & 1;
        stream << temp;
        if ((bitsPerValue - i) % blockSize == 0)
            stream << " ";
    }
    return stream;
}

template <typename TSpec>
inline std::ostream & operator<<(std::ostream & stream, const RankSupportBitString<TSpec> & rankSupportBitString)
{
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreBits>::Type                TBitString;
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreBlocks>::Type             TBlockString;
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreSuperBlocks>::Type        TSuperBlockString;

    typedef typename Value<TBitString>::Type          TBitStringValue;
    typedef typename Value<TBlockString>::Type       TBlockStringValue;
    typedef typename Value<TSuperBlockString>::Type  TSuperBlockStringValue;

    unsigned bitsPerBlock = BitsPerValue<typename Value<TBitString>::Type>::VALUE;

    TBitString const & bits = rankSupportBitString.bits;
    TBlockString const & blockString = rankSupportBitString.blocks;
    TSuperBlockString const & superBlocks = rankSupportBitString.superBlocks;

    stream << "  ";
    for (TBitStringValue i = 0; i < length(bits); i++)
    {
        printBits(stream, bits[i], bitsPerBlock);
    }
    stream << std::endl;

    for (TBlockStringValue i = 0; i < length(blockString); i++)
    {
        stream << blockString[i] << " ";
    }
    stream << std::endl;

    for (TSuperBlockStringValue i = 0; i < length(superBlocks); i++)
    {
        stream << superBlocks[i] << " ";
    }
    return stream;
}

}


#endif  // INDEX_FM_BRANK_SUPPORT_BIT_STRING_H_
