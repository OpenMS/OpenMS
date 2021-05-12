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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Specialization "Chained" for class Seed.
// ==========================================================================

#ifndef SEQAN_SEEDS_SEEDS_SEED_CHAINED_H_
#define SEQAN_SEEDS_SEEDS_SEED_CHAINED_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// ---------------------------------------------------------------------------
// Class Chained Seed
// ---------------------------------------------------------------------------

struct Chained_;
typedef Tag<Chained_> ChainedSeed;  // TODO(holtgrew): Chained already taken as template in file. Maybe prefer non-parameterized types for simpler names.

/**
.Spec.Chained Seed
..cat:Seed Handling
..summary:Describes a seed with start and end position2 and diagonal upper and lower bounds.
..description:Additionaly diagonal segments between start and end position2 are stored.
..general:Class.Seed
..signature:Seed<ChainedSeed, TConfig>
..param.TConfig:The configuration object used for the seed.
...default:DefaultSeedConfig.

.Memfunc.Chained Seed#Seed
..class:Spec.Chained Seed
..summary:Constructor
..signature: Seed<ChainedSeed, TConfig> ()
..signature: Seed<ChainedSeed, TConfig> (beginPosH, beginPosV, length)
..param.beginPosH: Begin position in database (horizontal).
..param.beginPosV: Begin position in query (vertical).
..param.length: Length of the seed.
..include:seqan/seeds.h
*/

template <typename TConfig>
class Seed<ChainedSeed, TConfig>
{
public:
    typedef typename TConfig::TPosition TPosition;
    typedef typename TConfig::TSize TSize;
    typedef typename TConfig::TDiagonal TDiagonal;
    typedef typename TConfig::TScoreValue TScoreValue;

    typedef SeedDiagonal<TPosition, TSize> TSeedDiagonal;

    std::list<TSeedDiagonal> _seedDiagonals;
    TDiagonal _lowerDiagonal;
    TDiagonal _upperDiagonal;

    TScoreValue _score;

    Seed() : _lowerDiagonal(0), _upperDiagonal(0), _score(0)
    {}

    Seed(TPosition beginPositionH, TPosition beginPositionV, TPosition seedLength) :
            _lowerDiagonal(beginPositionH - beginPositionV), _upperDiagonal(beginPositionH - beginPositionV), _score(0)
    {
        appendValue(_seedDiagonals, TSeedDiagonal(beginPositionH, beginPositionV, seedLength));
    }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

// ---------------------------------------------------------------------------
// Metafunction Value
// ---------------------------------------------------------------------------

/**
.Metafunction.Chained Seed#Value
..cat:Seed Handling
..class:Spec.Chained Seed
..summary:The seed diagonal type.
..signature:Value<TSeed>::Type
..param.TSeed:The seed to query for its diagonal object type.
...type:Spec.Chained Seed
..include:seqan/seeds.h
*/

template <typename TConfig>
struct Value<Seed<ChainedSeed, TConfig> >
{
    typedef Seed<ChainedSeed, TConfig> TSeed_;
    typedef typename Position<TSeed_>::Type TPosition_;
    typedef typename Size<TSeed_>::Type TSize_;

    typedef SeedDiagonal<TPosition_, TSize_> Type;
};

template <typename TConfig>
struct Value<Seed<ChainedSeed, TConfig> const>
{
    typedef Seed<ChainedSeed, TConfig> TSeed_;
    typedef typename Position<TSeed_>::Type TPosition_;
    typedef typename Size<TSeed_>::Type TSize_;

    typedef SeedDiagonal<TPosition_, TSize_> const Type;
};

// ---------------------------------------------------------------------------
// Metafunction Reference
// ---------------------------------------------------------------------------

/**
.Metafunction.Chained Seed#Reference
..cat:Seed Handling
..class:Spec.Chained Seed
..summary:The seed diagonal reference type.
..signature:Reference<TSeed>::Type
..param.TSeed:The seed to query for its seed diagonal reference type.
...type:Spec.Chained Seed
..include:seqan/seeds.h
*/

template <typename TConfig>
struct Reference<Seed<ChainedSeed, TConfig> >
{
    typedef Seed<ChainedSeed, TConfig> TSeed_;
    typedef typename Value<TSeed_>::Type TSeedDiagonal_;
    typedef TSeedDiagonal_ & Type;
};

template <typename TConfig>
struct Reference<Seed<ChainedSeed, TConfig> const>
{
    typedef Seed<ChainedSeed, TConfig> TSeed_;
    typedef typename Value<TSeed_>::Type TSeedDiagonal_;
    typedef TSeedDiagonal_ const & Type;
};

// ---------------------------------------------------------------------------
// Metafunction Iterator
// ---------------------------------------------------------------------------

/**
.Metafunction.Chained Seed#Iterator
..cat:Seed Handling
..class:Spec.Chained Seed
..summary:The seed diagonal iterator type.
..signature:Iterator<TSeed, Tag>::Type
..param.TSeed:The seed to query for its seed diagonal iterator type.
...type:Spec.Chained Seed
..param.Tag:The tag to select the iterator type with.
..include:seqan/seeds.h
*/

template <typename TConfig>
struct Iterator<Seed<ChainedSeed, TConfig>, Standard>
{
    typedef Seed<ChainedSeed, TConfig> TSeed_;
    typedef typename Value<TSeed_>::Type TSeedDiagonal_;
    typedef typename ::std::list<TSeedDiagonal_>::iterator Type;
};

template <typename TConfig>
struct Iterator<Seed<ChainedSeed, TConfig> const, Standard>
{
    typedef Seed<ChainedSeed, TConfig> TSeed_;
    typedef typename Value<TSeed_>::Type TSeedDiagonal_;
    typedef typename ::std::list<TSeedDiagonal_>::const_iterator Type;
};

// ===========================================================================
// Functions
// ===========================================================================

// ---------------------------------------------------------------------------
// Debug Function operator<<()
// ---------------------------------------------------------------------------

template <typename TStream, typename TConfig>
inline TStream &
operator<<(TStream & stream, Seed<ChainedSeed, TConfig> const & seed)
{
    typedef Seed<ChainedSeed, TConfig> const TSeed;
    typedef typename Iterator<TSeed>::Type TIterator;

    stream << "Seed<ChainedSeed, TConfig>([";
    for (TIterator it = begin(seed); it != end(seed); ++it) {
        if (it != begin(seed))
            stream << ", ";
        stream << *it;
    }
    stream << "])";
    return stream;
}

// ---------------------------------------------------------------------------
// Function operator==()
// ---------------------------------------------------------------------------

// TODO(holtgrew): Documentation comes through concept.

template <typename TConfig>
inline bool
operator==(Seed<ChainedSeed, TConfig> const & a, Seed<ChainedSeed, TConfig> const & b)
{
    return a._seedDiagonals == b._seedDiagonals &&
            a._upperDiagonal == b._upperDiagonal &&
            a._lowerDiagonal == b._lowerDiagonal;
}

// ---------------------------------------------------------------------------
// Function beginPositionH()
// ---------------------------------------------------------------------------

template <typename TConfig>
inline typename Position<Seed<ChainedSeed, TConfig> >::Type
beginPositionH(Seed<ChainedSeed, TConfig> const & seed)
{
	return front(seed._seedDiagonals).beginPositionH;
}

// ---------------------------------------------------------------------------
// Function endPositionH()
// ---------------------------------------------------------------------------

template <typename TConfig>
inline typename Position<Seed<ChainedSeed, TConfig> >::Type
endPositionH(Seed<ChainedSeed, TConfig> const & seed)
{
	return back(seed._seedDiagonals).beginPositionH + back(seed._seedDiagonals).length;
}

// ---------------------------------------------------------------------------
// Function beginPositionV()
// ---------------------------------------------------------------------------

template <typename TConfig>
inline typename Position<Seed<ChainedSeed, TConfig> >::Type
beginPositionV(Seed<ChainedSeed, TConfig> const & seed)
{
	return front(seed._seedDiagonals).beginPositionV;
}

// ---------------------------------------------------------------------------
// Function endPositionV()
// ---------------------------------------------------------------------------

template <typename TConfig>
inline typename Position<Seed<ChainedSeed, TConfig> >::Type
endPositionV(Seed<ChainedSeed, TConfig> const & seed)
{
	return back(seed._seedDiagonals).beginPositionV + back(seed._seedDiagonals).length;
}

// ---------------------------------------------------------------------------
// Function length()
// ---------------------------------------------------------------------------

/**
.Function.Chained Seed#length
..summary:Returns the number of diagonals in the chained seed.
..signature:TSize length(seed)
..class:Spec.Chained Seed
..param.seed:The seed to query.
...type:Spec.Chained Seed
..returns:The number of diagonals in the chained seed.
..include:seqan/seeds.h
*/

template <typename TConfig>
inline typename Size<Seed<ChainedSeed, TConfig> >::Type
length(Seed<ChainedSeed, TConfig> const & seed)
{
    return length(seed._seedDiagonals);
}

// ---------------------------------------------------------------------------
// Function appendDiagonal()
// ---------------------------------------------------------------------------

/**
.Function.appendDiagonal
..summary: Adds diagonal to the Chained Seed.
..cat:Seed Handling
..signature:appendDiag(seed, diagonal)
..class:Spec.Chained Seed
..param.seed: The seed to which the diagonal should be added.
...type:Spec.Chained Seed
..param.diag: The diagonal to add.
...type:Class.SeedDiagonal
...remarks: A diagonal consists of three values: 1: start in 1. sequence, 2: start in 2. sequence, 3: length of match
..include:seqan/seeds.h
*/

template <typename TConfig>
inline void
appendDiagonal(Seed<ChainedSeed, TConfig> & seed,
               typename Value<Seed<ChainedSeed, TConfig> >::Type const & diagonal)
{
    // TODO(holtgrew): Add empty().
    if (length(seed) > 0)
    {
        SEQAN_ASSERT_LEQ(back(seed._seedDiagonals).beginPositionH + back(seed._seedDiagonals).length, diagonal.beginPositionH);
        SEQAN_ASSERT_LEQ(back(seed._seedDiagonals).beginPositionV + back(seed._seedDiagonals).length, diagonal.beginPositionV);
    }

    appendValue(seed._seedDiagonals, diagonal);
}

// ---------------------------------------------------------------------------
// Function truncateDiagonals()
// ---------------------------------------------------------------------------

/**
.Function.truncateDiagonals
..summary:Removes diagonals from the given first one to the end of the seed's diagonals.
..cat:Seed Handling
..signature:truncateDiagonals(seed, first)
..class:Spec.Chained Seed
..param.seed: The seed to which the diagonal should be added.
...type:Spec.Chained Seed
..param.first: Iterator the first diagonal to remove.
..include:seqan/seeds.h
*/

template <typename TConfig>
inline void
truncateDiagonals(Seed<ChainedSeed, TConfig> & seed,
                  typename Iterator<Seed<ChainedSeed, TConfig> >::Type const & first)
{
     // TODO(holtgrew): Add erase() to std::list adaptors?
    seed._seedDiagonals.erase(first, seed._seedDiagonals.end());
}

// ---------------------------------------------------------------------------
// Function begin()
// ---------------------------------------------------------------------------

/**
.Function.Chained Seed#begin
..summary:Returns an iterator to the beginning of the seed digonals.
..class:Spec.Chained Seed
..signature:TIterator begin(seed, tag)
..include:seqan/seeds.h
*/

template <typename TConfig>
inline typename Iterator<Seed<ChainedSeed, TConfig> >::Type
begin(Seed<ChainedSeed, TConfig> & seed, Standard const &)
{
    return seed._seedDiagonals.begin();
}

template <typename TConfig>
inline typename Iterator<Seed<ChainedSeed, TConfig> const>::Type
begin(Seed<ChainedSeed, TConfig> const & seed, Standard const &)
{
    return seed._seedDiagonals.begin();
}

// ---------------------------------------------------------------------------
// Function end()
// ---------------------------------------------------------------------------

/**
.Function.Chained Seed#end
..summary:Returns an iterator to the end of the seed diagonals.
..class:Spec.Chained Seed
..signature:TIterator end(seed, tag)
..include:seqan/seeds.h
*/

template <typename TConfig>
inline typename Iterator<Seed<ChainedSeed, TConfig> >::Type
end(Seed<ChainedSeed, TConfig> & seed, Standard const &)
{
    return seed._seedDiagonals.end();
}

template <typename TConfig>
inline typename Iterator<Seed<ChainedSeed, TConfig> const>::Type
end(Seed<ChainedSeed, TConfig> const & seed, Standard const &)
{
    return seed._seedDiagonals.end();
}

// ---------------------------------------------------------------------------
// Function front()
// ---------------------------------------------------------------------------

/**
.Function.Chained Seed#front
..summary:Returns a reference to the first seed diagonal.
..class:Spec.Chained Seed
..signature:TReference front(seed)
..include:seqan/seeds.h
*/

template <typename TConfig>
inline typename Reference<Seed<ChainedSeed, TConfig> >::Type
front(Seed<ChainedSeed, TConfig> & seed)
{
    return front(seed._seedDiagonals);
}

template <typename TConfig>
inline typename Reference<Seed<ChainedSeed, TConfig> const>::Type
front(Seed<ChainedSeed, TConfig> const & seed)
{
    return front(seed._seedDiagonals);
}

// ---------------------------------------------------------------------------
// Function back()
// ---------------------------------------------------------------------------

/**
.Function.Chained Seed#back
..summary:Returns a reference to the last seed diagonal.
..class:Spec.Chained Seed
..signature:TReference back(seed)
..include:seqan/seeds.h
*/

template <typename TConfig>
inline typename Reference<Seed<ChainedSeed, TConfig> >::Type
back(Seed<ChainedSeed, TConfig> & seed)
{
    SEQAN_CHECKPOINT;
    return back(seed._seedDiagonals);
}

template <typename TConfig>
inline typename Reference<Seed<ChainedSeed, TConfig> const>::Type
back(Seed<ChainedSeed, TConfig> const & seed)
{
    SEQAN_CHECKPOINT;
    return back(seed._seedDiagonals);
}

// ---------------------------------------------------------------------------
// Debug Function
// ---------------------------------------------------------------------------

template <typename TStream, typename TConfig>
inline void
__write(TStream & stream, Seed<ChainedSeed, TConfig> const & seed, Tikz_ const &)
{
//IOREV _nodoc_ specialization not documented
    // Overall seed.
    stream << "\\draw[seed] (" << getBeginDim1(seed) << ", -" << getBeginDim0(seed) << ") -- (" << (getEndDim1(seed) - 1) << ", -" << (getEndDim0(seed) - 1) << ");" << std::endl;
    // Diagonals.
    typedef Seed<ChainedSeed, TConfig> TSeed;
    typedef typename Iterator<TSeed const, Standard>::Type TIterator;
    for (TIterator it = begin(seed); it != end(seed); ++it)
        stream << "\\draw[seed diagonal] (" << it->beginDim1 << ", -" << it->beginDim0 << ") -- (" << (it->beginDim1 + it->length - 1) << ", -" << (it->beginDim0 + it->length - 1) << ");" << std::endl;
}

}  // namespace seqan

#endif  // SEQAN_SEEDS_SEEDS_SEED_CHAINED_H_
