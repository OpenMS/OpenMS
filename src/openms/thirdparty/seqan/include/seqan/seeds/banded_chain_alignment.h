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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// The global interfaces to call to compute the banded chain alignment
// algorithm.
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_SEEDS_BANDED_CHAIN_ALIGNMENT_H_
#define CORE_INCLUDE_SEQAN_SEEDS_BANDED_CHAIN_ALIGNMENT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

/**
.Function.bandedChainAlignment
..summary:Computes the best global pairwise alignment between two sequences given a non-empty seed chain.
..cat:Alignments
..signature:bandedChainAlignment(align,          seedChain, scoringScheme1 [, scoringScheme2] [, alignConfig] [, k])
..signature:bandedChainAlignment(gapsH, gapsV,   seedChain, scoringScheme1 [, scoringScheme2] [, alignConfig] [, k])
..signature:bandedChainAlignment(frags, strings, seedChain, scoringScheme1 [, scoringScheme2] [, alignConfig] [, k])
..signature:bandedChainAlignment(alignmentGraph, seedChain, scoringScheme1 [, scoringScheme2] [, alignConfig] [, k])
..param.align:
An @Class.Align@ object that stores the alignment.
The number of rows must be 2 and the sequences must have already been set.
$row(align, 0)$ is the horizontal sequence in the alignment matrix, $row(align, 1)$ is the vertical sequence.
...type:Class.Align
..param.gapsH:Horizontal gapped sequence in alignment matrix.
...type:Class.Gaps
..param.gapsV:Vertical gapped sequence in alignment matrix.
...type:Class.Gaps
..param.frags:
String of @Class.Fragment@ objects.
The sequence with id $0$ is the horizontal one, the sequence with id $1$ is the vertical one.
..param.alignmentGraph:
@Spec.Alignment Graph@ object to store the alignment in.
...type:Spec.Alignment Graph
...remarks:The underlying @Class.StringSet@ must be an @Spec.Dependent|Dependent StringSet@.
..param.strings:A @Class.StringSet@ containing two sequences.
...type:Class.StringSet
..param.seedChain:
The container holding the @Class.Seed|seeds@. Note that the @Class.Seed|seeds@ have to be in montonic non-decreasing order
and the container has to implement a forward-iterator.
...see:Class.SeedSet
..param.scoringScheme1:The scoring scheme used for the alignment.
...remarks:If $scoringScheme2$ is specified, then $scoringScheme1$ is used for the regions around the seeds
and $scoringScheme2$ for the gap regions between two consecutive seeds.
...type:Class.Score
..param.scoringScheme2:The optional scoring scheme for the gap regions between two anchors.
...type:Class.Score
..param.alignConfig:The @Class.AlignConfig@ to use for the alignment.
...type:Class.AlignConfig
..param.k:Optional extension of the band around the seeds.
...type:nolink:$int$
...default:15
...remarks:At the moment only band extensions greater or equal $1$ are allowed.
..returns:An integer with the alignment score, as given by the @Metafunction.Value@ metafunction of the @Class.Score@ type.
If the seed chain is empty the metafunction @Metafunction.MinValue@ is used to return the minimal value of the selected score type and no alignment is computed.
..remarks:
There exist multiple overloads for this function with four configuration dimensions.
..remarks:
First, you can select whether begin and end gaps are free in either sequence using $alignConfig$.
..remarks:
Second, you can select the type of the target storing the alignment.
This can be either an @Class.Align@ object, two @Class.Gaps@ objects, a @Spec.Alignment Graph@, or a string of @Class.Fragment@ objects.
@Class.Align@ objects provide an interface to tabular alignments with the restriction of all rows having the same type.
Using two @Class.Gaps@ objects has the advantage that you can align sequences with different types, for example @Shortcut.DnaString@ and @Shortcut.Dna5String@.
@Spec.Alignment Graph|Alignment Graphs@ provide a graph-based representation of segment-based colinear alignments.
Using @Class.Fragment@ strings is useful for collecting many pairwise alignments, for example in the construction of @Spec.Alignment Graph|Alignment Graphs@ for multiple-sequence alignments (MSA).
..remarks:
Third, you can optionally give a second scoring scheme to fill the gaps between two consecutive seeds. Note that based on the specified scores
either an affine or linear gap cost function is used. This only depends on whether for one of the scoring schemes the scores for gap opening and gap extension differ or not.
If only one scoring scheme is defined the complete region is computed with the same scoring scheme.
..remarks:
Fourth, you can optinally select a proper band extension for the bands around the seeds. At the moment only band extensions of at least $1$ are allowed.
The default value is $15$ and is based on the default values for the LAGAN-algorithm described by Brudno et al., Genome Research, 2003.
..remarks:
The examples below show some common use cases.
..example.text:Banded chain alignment of two sequences using an @Class.Align@ object and using only one scoring scheme and no free end-gaps.
..example.code:
Dna5String seqH = "CGAATCCATCCCACACA";
Dna5String seqV = "GGCGATNNNCATGGCACA";

String<Seed<Simple> > seedChain;
appendValue(seedChain, Seed<Simple>(2, 0, 6, 5));
appendValue(seedChain, Seed<Simple>(9, 6, 12, 9));
appendValue(seedChain, Seed<Simple>(14, 11, 16, 17));

Align<Dna5String, ArrayGaps> alignment;
resize(rows(alignment), 2);
assignSource(row(alignment, 0), seqH);
assignSource(row(alignment, 1), seqV);

Score<int, Simple> scoringScheme(2, -1, -2);

int result = bandedChainAlignment(alignment, seedChain, scoringScheme, 2);

..example.text:Banded chain alignment of two sequences using two @Class.Gaps@ objects, an unordered seed set to hold the seeds,
two different scoring schemes for the gaps between the seeds and the seeds and free end-gaps.
..example.code:
DnaString seqH = "CGAATCCATCCCACACA";
Dna5String seqV = "GGCGATNNNCATGGCACA";

SeedSet<Simple, Unordered> seedChain;
addSeed(seedChain, Seed<Simple>(2, 0, 6, 5), Single());
addSeed(seedChain, Seed<Simple>(9, 6, 12, 9), Single());
addSeed(seedChain, Seed<Simple>(14, 11, 16, 17), Single());

Gaps<DnaString, ArrayGaps> gapsH(seqH);
Gaps<Dna5String, AnchorGaps<> > gapsV(seqV);

Score<int, Simple> scoringSchemeSeed(2, -1, -2);
Score<int, Simple> scoringSchemeGap(5, -3, -1, -5);
AlignConfig<true, true, true, true> alignConfig;

int result = globalAlignment(gapsH, gapsV, scoringSchemeSeed, scoringSchemeGap, alignConfig, 2);
..include:seqan/seeds.h
..wiki:Tutorial/Seed-and-Extend
..cite:Brudno M, Do CB, Cooper GM, et al.: LAGAN and Multi-LAGAN: Efficient Tools for Large-Scale Multiple Alignment of Genomic DNA. Genome Research 2003, 13: 721-731.
.
*/

// ----------------------------------------------------------------------------
// Function bandedChainAlignment()                                      [Align]
// ----------------------------------------------------------------------------

template <typename TSequence, typename TAlignSpec, typename TSeeds, typename TScoreValue, typename TScoreSpecAnchor,
          typename TScoreSpecGap, bool TFirstRow, bool TFirstColumn, bool TLastColumn, bool TLastRow, typename TACSpec>
inline TScoreValue
bandedChainAlignment(Align<TSequence, TAlignSpec> & align,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchor,
                     Score<TScoreValue, TScoreSpecGap> const & scoreSchemeGap,
                     AlignConfig<TFirstRow, TFirstColumn, TLastColumn, TLastRow, TACSpec> const & alignConfig,
                     unsigned bandExtension = 15)
{
    typedef typename Position<TSequence>::Type TPosition;
    typedef typename Size<TSequence>::Type TSize;
    typedef StringSet<String<TraceSegment_<TPosition, TSize> > > TTraceSegmentSet;

    TTraceSegmentSet traceSet;
    TScoreValue score =
        _setupAndRunBandedChainAlignment(traceSet, seedSet, source(row(align, 0)), source(row(align, 1)),
                                         scoreSchemeAnchor, scoreSchemeGap, alignConfig, bandExtension, GapsLeft());

    if (empty(traceSet))
        return score;

    _adaptTraceSegmentsTo(row(align,0), row(align,1), value(traceSet, 0));
    return score;
}

// With only one scoring scheme.
template <typename TSequence, typename TAlignSpec, typename TSeeds, typename TScoreValue, typename TScoreSpec,
          bool TFirstRow, bool TFirstColumn, bool TLastColumn, bool TLastRow, typename TACSpec>
inline TScoreValue
bandedChainAlignment(Align<TSequence, TAlignSpec> & align,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpec> const & scoreScheme,
                     AlignConfig<TFirstRow, TFirstColumn, TLastColumn, TLastRow, TACSpec> const & alignConfig,
                     unsigned bandExtension = 15)
{
    return bandedChainAlignment(align, seedSet, scoreScheme, scoreScheme, alignConfig, bandExtension);
}

// Without AlignConfig.
template <typename TSequence, typename TAlignSpec, typename TSeeds, typename TScoreValue, typename TScoreSpecAnchor,
          typename TScoreSpecGap>
inline TScoreValue
bandedChainAlignment(Align<TSequence, TAlignSpec> & align,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchor,
                     Score<TScoreValue, TScoreSpecGap> const & scoreSchemeGap,
                     unsigned bandExtension = 15)
{
    return bandedChainAlignment(align, seedSet, scoreSchemeAnchor, scoreSchemeGap, AlignConfig<>(), bandExtension);
}

// Without AlignConfig and with only one scoring scheme.
template <typename TSequence, typename TAlignSpec, typename TSeeds, typename TScoreValue, typename TScoreSpec>
inline TScoreValue
bandedChainAlignment(Align<TSequence, TAlignSpec> & align,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpec> const & scoreScheme,
                     unsigned bandExtension = 15)
{
    return bandedChainAlignment(align, seedSet, scoreScheme, scoreScheme, AlignConfig<>(), bandExtension);
}

// ----------------------------------------------------------------------------
// Function bandedChainAlignment()                                       [Gaps]
// ----------------------------------------------------------------------------

template <typename TSequenceH, typename TGapSpecH, typename TSequenceV, typename TGapSpecV, typename TSeeds,
          typename TScoreValue, typename TScoreSpecAnchor, typename TScoreSpecGap, bool TFirstRow, bool TFirstColumn,
          bool TLastColumn, bool TLastRow, typename TACSpec>
inline TScoreValue
bandedChainAlignment(Gaps<TSequenceH, TGapSpecH> & gapsHorizontal,
                     Gaps<TSequenceV, TGapSpecV> & gapsVertical,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchors,
                     Score<TScoreValue, TScoreSpecGap> const & scoreSchemeGaps,
                     AlignConfig<TFirstRow, TFirstColumn, TLastColumn, TLastRow, TACSpec> const & alignConfig,
                     unsigned bandExtension = 15)
{
    typedef typename Position<TSequenceH>::Type TPosition;
    typedef typename Size<TSequenceH>::Type TSize;
    typedef StringSet<String<TraceSegment_<TPosition, TSize> > > TTraceSegmentSet;

    TTraceSegmentSet traceSet;
    TScoreValue score =
        _setupAndRunBandedChainAlignment(traceSet, seedSet, source(gapsHorizontal), source(gapsVertical),
                                         scoreSchemeAnchors, scoreSchemeGaps, alignConfig, bandExtension, GapsLeft());

    if (empty(traceSet))
        return score;

    _adaptTraceSegmentsTo(gapsHorizontal, gapsVertical, value(traceSet, 0));
    return score;
}

// With only one scoring scheme.
template <typename TSequenceH, typename TGapSpecH, typename TSequenceV, typename TGapSpecV, typename TSeeds,
          typename TScoreValue, typename TScoreSpec, bool TFirstRow, bool TFirstColumn, bool TLastColumn, bool TLastRow,
          typename TACSpec>
inline TScoreValue
bandedChainAlignment(Gaps<TSequenceH, TGapSpecH> & gapsHorizontal,
                     Gaps<TSequenceV, TGapSpecV> & gapsVertical,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpec> const & scoreScheme,
                     AlignConfig<TFirstRow, TFirstColumn, TLastColumn, TLastRow, TACSpec> const & alignConfig,
                     unsigned bandExtension = 15)
{
    return bandedChainAlignment(gapsHorizontal, gapsVertical, seedSet, scoreScheme, scoreScheme, alignConfig,
                                bandExtension);
}

// Without AlignConfig.
template <typename TSequenceH, typename TGapSpecH, typename TSequenceV, typename TGapSpecV, typename TSeeds,
          typename TScoreValue, typename TScoreSpecAnchor, typename TScoreSpecGap>
inline TScoreValue
bandedChainAlignment(Gaps<TSequenceH, TGapSpecH> & gapsHorizontal,
                     Gaps<TSequenceV, TGapSpecV> & gapsVertical,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchor,
                     Score<TScoreValue, TScoreSpecGap> const & scoreSchemeGap,
                     unsigned bandExtension = 15)
{
    return bandedChainAlignment(gapsHorizontal, gapsVertical, seedSet, scoreSchemeAnchor, scoreSchemeGap,
                                AlignConfig<>(), bandExtension);
}

// Without AlignConfig and with only one scoring scheme.
template <typename TSequenceH, typename TGapSpecH, typename TSequenceV, typename TGapSpecV, typename TSeeds,
          typename TScoreValue, typename TScoreSpec>
inline TScoreValue
bandedChainAlignment(Gaps<TSequenceH, TGapSpecH> & gapsHorizontal,
                     Gaps<TSequenceV, TGapSpecV> & gapsVertical,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpec> const & scoreScheme,
                     unsigned bandExtension = 15)
{
    return bandedChainAlignment(gapsHorizontal, gapsVertical, seedSet, scoreScheme, scoreScheme, AlignConfig<>(),
                                bandExtension);
}

// ----------------------------------------------------------------------------
// Function bandedChainAlignment()                        [Graph<Alignment<> >]
// ----------------------------------------------------------------------------

template <typename TStringSet, typename TCargo, typename TGraphSpec, typename TSeeds, typename TScoreValue,
          typename TScoreSpecAnchor, typename TScoreSpecGap, bool TFirstRow, bool TFirstColumn, bool TLastColumn,
          bool TLastRow, typename TACSpec>
inline TScoreValue
bandedChainAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchors,
                     Score<TScoreValue, TScoreSpecGap> const & scoreSchemeGaps,
                     AlignConfig<TFirstRow, TFirstColumn, TLastColumn, TLastRow, TACSpec> const & alignConfig,
                     unsigned bandExtension = 15)
{
    typedef typename Position<TStringSet>::Type TPosition;
    typedef typename Size<TStringSet>::Type TSize;
    typedef StringSet<String<TraceSegment_<TPosition, TSize> > > TTraceSegmentSet;

    TTraceSegmentSet traceSet;
    TScoreValue score =
        _setupAndRunBandedChainAlignment(traceSet, seedSet, value(stringSet(alignmentGraph), 0),
                                         value(stringSet(alignmentGraph), 1), scoreSchemeAnchors, scoreSchemeGaps,
                                         alignConfig, bandExtension, GapsLeft());

    if (empty(traceSet))
        return score;

    _adaptTraceSegmentsTo(alignmentGraph, positionToId(stringSet(alignmentGraph), 0),
                          positionToId(stringSet(alignmentGraph), 1), value(traceSet, 0));
    return score;
}

// With only one scoring scheme.
template <typename TStringSet, typename TCargo, typename TGraphSpec, typename TSeeds, typename TScoreValue,
          typename TScoreSpec, bool TFirstRow, bool TFirstColumn, bool TLastColumn,
          bool TLastRow, typename TACSpec>
inline TScoreValue
bandedChainAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpec> const & scoreScheme,
                     AlignConfig<TFirstRow, TFirstColumn, TLastColumn, TLastRow, TACSpec> const & alignConfig,
                     unsigned bandExtension = 15)
{
    return bandedChainAlignment(alignmentGraph, seedSet, scoreScheme, scoreScheme, alignConfig, bandExtension);
}

// Without AlignConfig.
template <typename TStringSet, typename TCargo, typename TGraphSpec, typename TSeeds, typename TScoreValue,
          typename TScoreSpecAnchor, typename TScoreSpecGap>
inline TScoreValue
bandedChainAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchor,
                     Score<TScoreValue, TScoreSpecGap> const & scoreSchemeGap,
                     unsigned bandExtension = 15)
{
    return bandedChainAlignment(alignmentGraph, seedSet, scoreSchemeAnchor, scoreSchemeGap, AlignConfig<>(),
                                bandExtension);
}

// Without AlignConfig and with only one scoring scheme.
template <typename TStringSet, typename TCargo, typename TGraphSpec, typename TSeeds, typename TScoreValue,
          typename TScoreSpec>
inline TScoreValue
bandedChainAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpec> const & scoreScheme,
                     unsigned bandExtension = 15)
{
    return bandedChainAlignment(alignmentGraph, seedSet, scoreScheme, scoreScheme, AlignConfig<>(), bandExtension);
}

// ----------------------------------------------------------------------------
// Function bandedChainAlignment()                          [String<Fragments>]
// ----------------------------------------------------------------------------

template <typename TSize, typename TFragmentSpec, typename TStringSpec, typename TSequence, typename TStringSetSpec,
          typename TSeeds, typename TScoreValue, typename TScoreSpecAnchor, typename TScoreSpecGap,
          bool TFirstRow, bool TFirstColumn, bool TLastColumn, bool TLastRow, typename TACSpec>
inline TScoreValue
bandedChainAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                     StringSet<TSequence, TStringSetSpec> const & strings,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchors,
                     Score<TScoreValue, TScoreSpecGap> const & scoreSchemeGaps,
                     AlignConfig<TFirstRow, TFirstColumn, TLastColumn, TLastRow, TACSpec> const & alignConfig,
                     unsigned bandExtension = 15)
{
    typedef typename Position<TSequence>::Type TPosition;
    typedef StringSet<String<TraceSegment_<TPosition, TSize> > > TTraceSegmentSet;

    TTraceSegmentSet traceSet;
    TScoreValue score =
        _setupAndRunBandedChainAlignment(traceSet, seedSet, value(strings, 0), value(strings, 1), scoreSchemeAnchors,
                                         scoreSchemeGaps, alignConfig, bandExtension, GapsLeft());
    if (empty(traceSet))
        return score;

    _adaptTraceSegmentsTo(fragmentString, positionToId(strings, 0), positionToId(strings, 1), value(traceSet, 0));
    return score;
}

// With only one scoring scheme.
template <typename TSize, typename TFragmentSpec, typename TStringSpec, typename TSequence, typename TStringSetSpec,
          typename TSeeds, typename TScoreValue, typename TScoreSpec, bool TFirstRow, bool TFirstColumn,
          bool TLastColumn, bool TLastRow, typename TACSpec>
inline TScoreValue
bandedChainAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                     StringSet<TSequence, TStringSetSpec> const & strings,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpec> const & scoreScheme,
                     AlignConfig<TFirstRow, TFirstColumn, TLastColumn, TLastRow, TACSpec> const & alignConfig,
                     unsigned bandExtension = 15)
{
    return bandedChainAlignment(fragmentString, strings, seedSet, scoreScheme, scoreScheme, alignConfig, bandExtension);
}

// Without AlignConfig.
template <typename TSize, typename TFragmentSpec, typename TStringSpec, typename TSequence, typename TStringSetSpec,
          typename TSeeds, typename TScoreValue, typename TScoreSpecAnchor, typename TScoreSpecGap>
inline TScoreValue
bandedChainAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                     StringSet<TSequence, TStringSetSpec> const & strings,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchor,
                     Score<TScoreValue, TScoreSpecGap> const & scoreSchemeGap,
                     unsigned bandExtension = 15)
{
    return bandedChainAlignment(fragmentString, strings, seedSet, scoreSchemeAnchor, scoreSchemeGap, AlignConfig<>(),
                                bandExtension);
}

// Without AlignConfig and with only one scoring scheme.
template <typename TSize, typename TFragmentSpec, typename TStringSpec, typename TSequence, typename TStringSetSpec,
          typename TSeeds, typename TScoreValue, typename TScoreSpec>
inline TScoreValue
bandedChainAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                     StringSet<TSequence, TStringSetSpec> const & strings,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpec> const & scoreScheme,
                     unsigned bandExtension = 15)
{
    return bandedChainAlignment(fragmentString, strings, seedSet, scoreScheme, scoreScheme, AlignConfig<>(),
                                bandExtension);
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_SEEDS_BANDED_CHAIN_ALIGNMENT_H_
