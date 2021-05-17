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

#ifndef INDEX_FM_RANK_SUPPORT_BIT_STRING_ITERATOR_H_
#define INDEX_FM_RANK_SUPPORT_BIT_STRING_ITERATOR_H_

namespace seqan {
// ==========================================================================
// Metafunctions
// ==========================================================================

template <typename TSpec>
struct Iterator<RankSupportBitString<TSpec> const, Standard>
{
    typedef Iter<RankSupportBitString<TSpec> const, PositionIterator> Type;
};

template <typename TSpec>
struct Iterator<RankSupportBitString<TSpec>, Standard>
{
    typedef Iter<RankSupportBitString<TSpec>, PositionIterator> Type;
};

template <typename TSpec>
struct Iterator<RankSupportBitString<TSpec>, Rooted>:
    Iterator<RankSupportBitString<TSpec>, Standard>{};

template <typename TSpec>
struct Iterator<RankSupportBitString<TSpec> const, Rooted>:
    Iterator<RankSupportBitString<TSpec> const, Standard>{};

// ==========================================================================
// Functions
// ==========================================================================
template <typename TSpec>
inline typename Iterator<RankSupportBitString<TSpec>, Standard>::Type
begin(RankSupportBitString<TSpec> & rsbs, Standard)
{
    return typename Iterator<RankSupportBitString<TSpec>, Standard>::Type(rsbs, 0u);
}

template <typename TSpec>
inline typename Iterator<RankSupportBitString<TSpec>, Standard>::Type
begin(RankSupportBitString<TSpec> const & rsbs, Standard)
{
    return typename Iterator<RankSupportBitString<TSpec> const, Standard>::Type(rsbs, 0u);
}

template <typename TSpec>
inline typename Iterator<RankSupportBitString<TSpec>, Rooted>::Type
begin(RankSupportBitString<TSpec> & rsbs, Rooted)
{
    return typename Iterator<RankSupportBitString<TSpec>, Rooted>::Type(rsbs, 0u);
}

template <typename TSpec>
inline typename Iterator<RankSupportBitString<TSpec>, Rooted>::Type
begin(RankSupportBitString<TSpec> const & rsbs, Rooted)
{
    return typename Iterator<RankSupportBitString<TSpec> const, Rooted>::Type(rsbs, 0u);
}

// ==========================================================================
template <typename TSpec>
inline typename Iterator<RankSupportBitString<TSpec>, Standard>::Type
end(RankSupportBitString<TSpec> & rsbs, Standard)
{
    return typename Iterator<RankSupportBitString<TSpec>, Standard >::Type(rsbs, length(rsbs));
}

template <typename TSpec>
inline typename Iterator<RankSupportBitString<TSpec> const, Standard>::Type
end(RankSupportBitString<TSpec> const & rsbs, Standard)
{
    return typename Iterator<RankSupportBitString<TSpec> const, Standard >::Type(rsbs, length(rsbs));
}

template <typename TSpec>
inline typename Iterator<RankSupportBitString<TSpec>, Rooted>::Type
end(RankSupportBitString<TSpec> & rsbs, Rooted)
{
    return typename Iterator<RankSupportBitString<TSpec>, Rooted >::Type(rsbs, length(rsbs));
}

template <typename TSpec>
inline typename Iterator<RankSupportBitString<TSpec> const, Rooted>::Type
end(RankSupportBitString<TSpec> const & rsbs, Rooted)
{
    return typename Iterator<RankSupportBitString<TSpec> const, Rooted >::Type(rsbs, length(rsbs));
}

// ==========================================================================
struct Bit_;
struct Rank_;

typedef Tag<Bit_>   const Bit;
typedef Tag<Rank_>  const Rank;

template <typename TSpec>
inline bool getValue(Iter<RankSupportBitString<TSpec>, PositionIterator> const & it)
{
    return isBitSet(value(it), position(it));
}

template <typename TSpec>
inline bool getValue(Iter<RankSupportBitString<TSpec>, PositionIterator> & it)
{
    return isBitSet(value(it), position(it));
}

template <typename TSpec>
inline bool getValue(Iter<RankSupportBitString<TSpec>, PositionIterator> const & it, Bit)
{
    return getValue(it);
}

template <typename TSpec>
inline bool getValue(Iter<RankSupportBitString<TSpec>, PositionIterator> & it, Bit)
{
    return getValue(it);
}

// ==========================================================================
template <typename TSpec>
inline typename Size<RankSupportBitString<TSpec> >::Type
getValue(Iter<RankSupportBitString<TSpec>, PositionIterator> const & it, Rank)
{
    return getRank(value(it), position(it));
}

template <typename TSpec>
inline typename Size<RankSupportBitString<TSpec> >::Type
getValue(Iter<RankSupportBitString<TSpec>, PositionIterator> & it, Rank)
{
    return getRank(value(it), position(it));
}

}
#endif // INDEX_FM_COMPRESSED_SA_ITERATOR_H_

