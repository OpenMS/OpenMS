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
// The class SeedSet.  ScoringScheme and related tags are defined in
// seeds_scoring_scheme.h
// ==========================================================================

#ifndef SEQAN_SEEDS_SEEDS_SEED_SET_BASE_H_
#define SEQAN_SEEDS_SEEDS_SEED_SET_BASE_H_

namespace seqan {

// ===========================================================================
// Forwards
// ===========================================================================

// TODO(holtgrew): Add DiagonalSorted as default specialization.

struct Unordered_;
typedef Tag<Unordered_> Unordered;

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// ---------------------------------------------------------------------------
// Class SeedSet
// ---------------------------------------------------------------------------

/**
.Class.SeedSet:
..summary:Handles a set of seeds with local chaining on adding seeds.
..cat:Seed Handling
..signature:SeedSet<TSeedSpec>
..param.TSeedSpec:Specialization of the underlying seed to use.
..include:seqan/seeds.h
*/

// ..param.TScored:Either UnScored or a seed set scoring scheme specification.
// ..param.TSpec:Specialization of the seed set.
// ..param.TSeedConfig:Configuration for the seeds.  Sensible defaults are chosen based on the other template parameters.

template <typename TSeedSpec, typename TSpec = Unordered>
class SeedSet;

// ===========================================================================
// Metafunctions
// ===========================================================================

///.Metafunction.Position.param.T.type:Class.SeedSet
///.Metafunction.Position.class:Class.SeedSet
///.Metafunction.Size.param.T.type:Class.SeedSet
///.Metafunction.Size.class:Class.SeedSet
///.Metafunction.Value.param.T.type:Class.SeedSet
///.Metafunction.Value.class:Class.SeedSet
///.Metafunction.GetValue.param.T.type:Class.SeedSet
///.Metafunction.GetValue.class:Class.SeedSet
///.Metafunction.Reference.param.T.type:Class.SeedSet
///.Metafunction.Reference.class:Class.SeedSet
///.Metafunction.Iterator.param.T.type:Class.SeedSet
///.Metafunction.Iterator.class:Class.SeedSet


// ===========================================================================
// Functions
// ===========================================================================

// Basic Container Functions

/**
.Function.begin.param.object.type:Class.SeedSet
..class:Class.SeedSet
.Function.end.param.object.type:Class.SeedSet
..class:Class.SeedSet
.Function.length.param.object.type:Class.SeedSet
..class:Class.SeedSet
.Function.front.param.container.type:Class.SeedSet
..class:Class.SeedSet
.Function.back.param.container.type:Class.SeedSet
..class:Class.SeedSet
*/
// TODO(holtgrew): dddoc {begin,end,length,front,back}All(T)

// SeedSet Functions

/**
.Function.SeedSet#addSeed:
..summary:Adds a seed to an existing @Class.SeedSet@ using different algorithms for local chaining.
..cat:Seed Handling
..signature:addSeed(set, seed, distance, bandwidth, score, seqH, seqV, tag)
..signature:addSeed(set, seed, distance, score, SimpleChain())
..signature:addSeed(set, seed, distance, Merge())
..signature:addSeed(set, seed, Single())
..class:Class.SeedSet
..param.set:The set to add the seed to.
...type:Class.SeedSet
..param.seed:The seed to be added.
...type:Class.Seed
..param.seqH: Database sequence (horizontal).
...type:Concept.SequenceConcept
...remarks:Only required for @Tag.Local Chaining|Chaos Chaining@.
..param.seqV: Query sequence (vertical).
...type:Concept.SequenceConcept
...remarks:Only required for @Tag.Local Chaining|Chaos Chaining@.
..param.score:The scoring scheme.
...type:Spec.Simple Score
...remarks:Note, only @Tag.Local Chaining|Chaos and SimpleChain@ require the score.
..param.distance:The maximal distance between the end point of the upper left and the begin point of the lower right @Class.Seed@ allowed for local chaining.
...type:Concept.IntegerConcept
...remarks: Note, only @Tag.Local Chaining| Chaos, SimpleChain and Merge@ require the distance information.
..param.bandwidth: The window size to search for a chainable @Class.Seed@.
...type:Concept.IntegerConcept
...remarks: Note, only @Tag.Local Chaining|Chaos@ requires the bandwidth information.
..param.tag: The algorithm that is used to add the new seed.
...type:Tag.Local Chaining
...remarks: Note that not every algorithm can be used with each type of @Class.Seed@. See special signatures above.
...remarks: The seed is copied and then added.
..returns:Boolean if successfully added.
...remarks:Always true for Tag Single.
..example.file:demos/seeds/seeds_add_seed.cpp
..example.text:The output is as follows:
..example.output:
Single Method:
Seed: Seed<Simple, TConfig>(4, 5, 8, 9, lower diag = -1, upper diag = -1)
Seed: Seed<Simple, TConfig>(10, 10, 15, 15, lower diag = 0, upper diag = 0)
Seed: Seed<Simple, TConfig>(14, 14, 18, 18, lower diag = 0, upper diag = 0)
Seed: Seed<Simple, TConfig>(21, 21, 24, 24, lower diag = 0, upper diag = 0)


Merge Method:
Seed: Seed<Simple, TConfig>(4, 5, 8, 9, lower diag = -1, upper diag = -1)
Seed: Seed<Simple, TConfig>(10, 10, 18, 18, lower diag = 0, upper diag = 0)
Seed: Seed<Simple, TConfig>(21, 21, 24, 24, lower diag = 0, upper diag = 0)


Chaos Method:
Seed: Seed<Simple, TConfig>(4, 5, 15, 15, lower diag = -1, upper diag = 0)
Seed: Seed<Simple, TConfig>(14, 14, 18, 18, lower diag = 0, upper diag = 0)
Seed: Seed<Simple, TConfig>(21, 21, 24, 24, lower diag = 0, upper diag = 0)
..include:seqan/seeds.h
*/

/**
.Intenral.Function.SeedSet#addSeeds:
..summary:Adds several seeds to an existing set. If a merging or chaining algorithm is used seeds are added if the merging or chaining fails.
..cat:Seed Handling
..signature:addSeed(set, container, tag)
..signature:addSeed(set, begin, end, tag)
..class:Class.SeedSet
..param.set:The set to which the new seed sould be added.
...type:Class.SeedSet
..param.container: Content is copied to set.
...type:Concept.ContainerConcept
..param.begin: Iterator pointing to the first value to add.
..param.end: Iterator pointing just behind the last value to add.
..param.tag: The algorithm that should be used to add the new @Class.Seed@.
...type:Tag.Local Chaining
...remarks: Note that not every algorithm can be used with each specialization of @Class.Seed@.
..include:seqan/seeds.h
*/

// ---------------------------------------------------------------------------
// Function minScore()
// ---------------------------------------------------------------------------

/**
.Function.SeedSet#minScore:
..summary:Returns the threshold to distinguish between high-scoring and low-scoring seeds.
..cat:Seed Handling
..signature:minScore(set)
..class:Class.SeedSet
..param.set:The set for which the threshold is set.
...type:Class.SeedSet
...remarks: If the score of a seed is higher than the given threshold, then it is virtually put
into a container storing the high-scoring seeds which can be iterated separately.
..see:Function.SeedSet#setMinScore
..include:seqan/seeds.h
*/

template <typename TSeedSpec, typename TSeedSetSpec>
typename SeedScore<typename Value<SeedSet<TSeedSpec, TSeedSetSpec> >::Type >::Type
minScore(SeedSet<TSeedSpec, TSeedSetSpec> const & seedSet)
{
    return seedSet._minScore;
}

// ---------------------------------------------------------------------------
// Function setMinScore()
// ---------------------------------------------------------------------------

/**
.Function.SeedSet#setMinScore:
..summary:Sets the threshold at which seeds are considered high-scoring.
..cat:Seed Handling
..signature:setMinScore(set, score)
..class:Class.SeedSet
..param.set:The seed set for which the threshold is to be set.
...type:Class.SeedSet
..param.score:The new threshold to set.
...remarks: If the score of a seed is higher than the given threshold, then it is virtually put
into a container storing the high-scoring seeds which can be iterated separately.
..see:Function.SeedSet#minScore
..include:seqan/seeds.h
*/

template <typename TSeedSpec, typename TSeedSetSpec, typename TScoreValue>
void setMinScore(SeedSet<TSeedSpec, TSeedSetSpec> & seedSet, TScoreValue val)
{
    seedSet._minScore = val;
}

// ---------------------------------------------------------------------------
// Function minSeedSize()
// ---------------------------------------------------------------------------

template <typename TSeedSpec, typename TSeedSetSpec>
typename Size<typename Value<SeedSet<TSeedSpec, TSeedSetSpec> >::Type >::Type
minSeedSize(SeedSet<TSeedSpec, TSeedSetSpec> const & seedSet)
{
    return seedSet._minSeedSize;
}

// ---------------------------------------------------------------------------
// Function setMinSeedSize()
// ---------------------------------------------------------------------------

template <typename TSeedSpec, typename TSeedSetSpec, typename TSize>
void setMinSeedSize(SeedSet<TSeedSpec, TSeedSetSpec> & seedSet, TSize siz)
{
    seedSet._minSeedSize = siz;
}

// ---------------------------------------------------------------------------
// Helper Function _qualityReached()
// ---------------------------------------------------------------------------

template <typename TSeedSpec, typename TSeedSetSpec, typename TSeedConfig>
inline bool _qualityReached(SeedSet<TSeedSpec, TSeedSetSpec> const & seedSet,
                            Seed<TSeedSpec, TSeedConfig> const & seed)
{
    return score(seed) >= minScore(seedSet) && seedSize(seed) >= minSeedSize(seedSet);
}

// ---------------------------------------------------------------------------
// Function clear()
// ---------------------------------------------------------------------------

///.Function.clear.param.object.type:Class.SeedSet
template <typename TSeedSpec, typename TSeedSetSpec>
inline void clear(SeedSet<TSeedSpec, TSeedSetSpec> & seedSet)
{
    seedSet._seeds.clear();
    seedSet._minScore = 0;
    seedSet._minSeedSize = 0;
}


// Debugging / TikZ Output

template <typename TStream, typename TQuerySequence, typename TDatabaseSequence, typename TSeedSetSpec, typename TSeedSpec>
inline void
__write(TStream & stream,
       TQuerySequence & sequence0,
       TDatabaseSequence & sequence1,
       SeedSet<TSeedSpec, TSeedSetSpec> const & seedSet,
       Tikz_ const &)
{
//IOREV _nodoc_ specialization not documented
    typedef SeedSet<TSeedSpec, TSeedSetSpec> TSeedSet;

    stream << "\\begin{tikzpicture}[" << std::endl
           << "    seed/.style={very thick}," << std::endl
           << "    seed diagonal/.style={red,<->}" << std::endl
           << "    ]" << std::endl;

    // Draw sequences.
    stream << "  \\draw";
    // Draw query / sequence 0;
    for (unsigned i = 0; i < length(sequence0); ++i)
        stream << std::endl << "    (0, -" << i << ") node {" << sequence0[i] << "}";
    stream << std::endl;
    // Draw database / sequence 1.
    for (unsigned i = 0; i < length(sequence1); ++i)
        stream << std::endl << "    (" << i << ", 0) node {" << sequence1[i] << "}";
    stream << ";" << std::endl;

    // Draw seeds.
    typedef typename Iterator<TSeedSet const, Standard>::Type TIterator;
    for (TIterator it = begin(seedSet); it != end(seedSet); ++it)
        _write(stream, value(it), Tikz_());
    stream << "\\end{tikzpicture}" << std::endl;
}

}  // namespace seqan

#endif  // SEQAN_SEEDS_SEEDS_SEED_SET_BASE_H_
