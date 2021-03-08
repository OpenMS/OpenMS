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
//
// Based on the code by Carsten Kemena <carsten.kemena@crg.es>, debugged
// by Birte Kehr <birte.kehr@fu-berlin.de>.
// ==========================================================================
// Seed extension algorithms.
//
// The approach for gapped X-drop extension is based on the algorithm in
// Figure 2 from (Zhang et al., 2000).
//
//  Zhang Z, Schwartz S, Wagner L, Miller W.  A greedy algorithm for aligning
//  DNA sequences.  Journal of computational biologyA a journal of
//  computational molecular cell biology.  2000;7(1-2):203-14.
//  doi:10.1089/10665270050081478
// ==========================================================================

#ifndef SEQAN_SEEDS_SEEDS_EXTENSION_H_
#define SEQAN_SEEDS_SEEDS_EXTENSION_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// ---------------------------------------------------------------------------
// Tags Seed Extension
// ---------------------------------------------------------------------------

/**
.Tag.Seed Extension
..cat:Seed Handling
..summary:The algorithms used to extend a seed.
..see:Function.extendSeed
..tag.MatchExtend:Extends a seed until a mismatch occurs.
..tag.UngappedXDrop:Ungapped extension of a seed until score drops below a Value.
..tag.GappedXDrop:Gapped extension of a seed until score drops below a Value. Only @Spec.Simple Seed@s.
..include:seqan/seeds.h
*/
struct MatchExtend_;
typedef Tag<MatchExtend_> const MatchExtend;

struct UngappedXDrop_;
typedef Tag<UngappedXDrop_> const UnGappedXDrop;

struct GappedXDrop_;
typedef Tag<GappedXDrop_> const GappedXDrop;

// ---------------------------------------------------------------------------
// Enum ExtensionDirection
// ---------------------------------------------------------------------------

// TODO(holtgrew): Put into class?

/**
.Enum.Extension Direction
..cat:Seed Handling
..summary:The direction in which a seed should be extended.
..value.EXTEND_NONE:Perform no extension.
..value.EXTEND_LEFT:Extend the seed to the left.
..value.EXTEND_RIGHT:Extend the seed to the right.
..value.EXTEND_BOTH:Extend the seed in both directions.
..see:Function.extendSeed
..include:seqan/seeds.h
*/

enum ExtensionDirection
{
    EXTEND_NONE  = 0,
    EXTEND_LEFT  = 1,
    EXTEND_RIGHT = 2,
    EXTEND_BOTH  = 3
};

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

// ---------------------------------------------------------------------------
// Function extendSeed                                           [MatchExtend]
// ---------------------------------------------------------------------------

/**
.Function.extendSeed
..summary:Extends a seed.
..cat:Seed Handling
..signature:void extendSeed(seed, database, query, direction, MatchExtend);
..signature:void extendSeed(seed, database, query, direction, scoreMatrix, scoreDropOff, {UngappedXDrop, GappedXDrop});
..class:Class.Seed
..param.seed: The seed to extend.
...type:Class.Seed
..param.query: The query sequence (vertical).
...type:Class.String
..param.database: The database sequence (horizontal).
...type:Class.String
..param.direction: Defines the direction in which the seed should be extended.
...type:Enum.Extension Direction
..param.scoreDropOff: The score drop after which the extension should stop. The extension stops if this value is exceeded.
...remarks:Only used for the algorithms @Tag.Seed Extension.UngappedXDrop@ and @Tag.Seed Extension.GappedXDrop@
..param.scoreMatrix: The scoring scheme.
...type:Spec.Simple Score
...remarks:Only used for the algorithms @Tag.Seed Extension.UngappedXDrop@ and @Tag.Seed Extension.GappedXDrop@
..remarks:You can use the tags, @Tag.Seed Extension.MatchExtend@, @Tag.Seed Extension.UngappedXDrop@, and @Tag.Seed Extension.GappedXDrop@.
..remarks:Note that the diagonals updated in $seed$ do not necessarily reflect the diagonals for the optimal extension but the diagonals used in all traces of the extension.
However, they are guaranteed to include the optimal extension's trace.
..example.text:The documentation of the class @Class.Seed@ contains an example for seed extension.
..include:seqan/seeds.h
*/

// We need one specialization for each combination of the extension variants and seeds.  It is not worth to extract the
// common parts for simple and chained seeds.

template <typename TConfig, typename TDatabase, typename TQuery>
inline void 
extendSeed(Seed<Simple, TConfig> & seed, 
		   TDatabase const & database,
		   TQuery const & query,
		   ExtensionDirection direction,
		   MatchExtend const &)
{
    // For match extension of Simple Seeds, we can simply update the begin and end values in each dimension.
    typedef Seed<Simple, TConfig> TSeed;
    typedef typename Position<TSeed>::Type TPosition;
    typedef typename Size<TSeed>::Type TSize;
    
	// Extension to the left
	if (direction == EXTEND_LEFT || direction == EXTEND_BOTH)
    {
		TPosition posH = beginPositionH(seed);
		TPosition posV = beginPositionV(seed);
		while (posH >= 1 && posV >= 1 && database[posH - 1] == query[posV - 1])
        {
			--posH;
			--posV;
		}
		setBeginPositionH(seed, posH);
		setBeginPositionV(seed, posV);
	}

	// Extension to the right
	if (direction == EXTEND_RIGHT || direction == EXTEND_BOTH)
    {
		TSize lengthH = length(database);
		TSize lengthV = length(query);
		TPosition posH = endPositionH(seed);
		TPosition posV = endPositionV(seed);
		while (posH < lengthH && posV < lengthV && database[posH] == query[posV])
        {
			++posH;
			++posV;
		}
        setEndPositionH(seed, posH);
        setEndPositionV(seed, posV);
	}
}


template <typename TConfig, typename TDatabase, typename TQuery>
inline void 
extendSeed(Seed<ChainedSeed, TConfig> & seed, 
		   TDatabase const & database,
		   TQuery const & query,
		   ExtensionDirection direction,
		   MatchExtend const &)
{
    // For match extension of Chained Seeds, we extend the first and the last Seed Diagonal.
    SEQAN_ASSERT_GT(length(seed), 0u);

    typedef Seed<ChainedSeed, TConfig> TSeed;
    typedef typename Value<TSeed>::Type TSeedDiagonal;
    typedef typename Position<TSeedDiagonal>::Type TPosition;
    typedef typename Size<TSeedDiagonal>::Type TSize;
    
	// Extension to the left
	if (direction == EXTEND_LEFT || direction == EXTEND_BOTH)
    {
        TSeedDiagonal & diag = front(seed);
		TPosition posH = diag.beginPositionH;
		TPosition posV = diag.beginPositionV;
        TSize diagonalLength = diag.length;
		while (posH >= 1 && posV >= 1 && database[posH - 1] == query[posV - 1])
        {
			--posH;
			--posV;
            ++diagonalLength;
		}
        diag.beginPositionH = posH;
        diag.beginPositionV = posV;
        diag.length = diagonalLength;
	}

	// Extension to the right
	if (direction == EXTEND_RIGHT || direction == EXTEND_BOTH)
    {
		TSize lengthH = length(database);
		TSize lengthV = length(query);
        TSeedDiagonal & diag = back(seed);
		TPosition posH = diag.beginPositionH + diag.length;
		TPosition posV = diag.beginPositionV + diag.length;
        TSize diagonalLength = diag.length;
		while (posH < lengthH && posV < lengthV && database[posH] == query[posV])
        {
			++posH;
			++posV;
            ++diagonalLength;
		}
        diag.length = diagonalLength;
	}
}

// ---------------------------------------------------------------------------
// Function extendSeed                                         [UnGappedXDrop]
// ---------------------------------------------------------------------------

template <typename TConfig, typename TDatabase, typename TQuery, typename TScoreValue, typename TScoreSpec>
inline void 
extendSeed(Seed<Simple, TConfig> & seed,
           TDatabase const & database,
           TQuery const & query,
           ExtensionDirection direction,
           Score<TScoreValue, TScoreSpec> const & scoringScheme,
           TScoreValue scoreDropOff,
           UnGappedXDrop const &)
{
    // For ungapped X-drop extension of Simple Seeds, we can simply
    // update the begin and end values in each dimension.
    scoreDropOff = -scoreDropOff;

    typedef Seed<ChainedSeed, TConfig> TSeed;
    typedef typename Position<TSeed>::Type TPosition;
    typedef typename Size<TSeed>::Type TSize;
    
	// Extension to the left
	if (direction == EXTEND_LEFT || direction == EXTEND_BOTH)
    {
        TScoreValue tmpScore = 0;
		TPosition posH = beginPositionH(seed);
		TPosition posV = beginPositionV(seed);
        TPosition mismatchingSuffixLength = 0;
		while (posH >= 1 && posV >= 1 && tmpScore > scoreDropOff)
        {
            tmpScore += score(scoringScheme, sequenceEntryForScore(scoringScheme, database, posH),
                              sequenceEntryForScore(scoringScheme, query, posV));
            if (database[posH - 1] == query[posV - 1])
            {
                mismatchingSuffixLength = 0;
                if (tmpScore > static_cast<TScoreValue>(0))
                    tmpScore = 0;
            }
            else
            {
                mismatchingSuffixLength += 1;
            }
			--posH;
			--posV;
		}
		setBeginPositionH(seed, posH + mismatchingSuffixLength);
		setBeginPositionV(seed, posV + mismatchingSuffixLength);
	}

	// Extension to the right
	if (direction == EXTEND_RIGHT || direction == EXTEND_BOTH)
    {
        TScoreValue tmpScore = 0;
		TSize lengthH = length(database);
		TSize lengthV = length(query);
		TPosition posH = endPositionH(seed);
		TPosition posV = endPositionV(seed);
        TPosition mismatchingSuffixLength = 0;
		while (posH < lengthH && posV < lengthV && tmpScore > scoreDropOff)
        {
            tmpScore += score(scoringScheme, sequenceEntryForScore(scoringScheme, database, posH),
                              sequenceEntryForScore(scoringScheme, query, posV));
            if (database[posH] == query[posV])
            {
                mismatchingSuffixLength = 0;
                if (tmpScore > static_cast<TScoreValue>(0))
                    tmpScore = 0;
            }
            else
            {
                mismatchingSuffixLength += 1;
            }
            ++posH;
            ++posV;
		}
		setEndPositionH(seed, posH - mismatchingSuffixLength);
		setEndPositionV(seed, posV - mismatchingSuffixLength);
    }

    // TODO(holtgrew): Update score?!
}


template <typename TConfig, typename TDatabase, typename TQuery, typename TScoreValue, typename TScoreSpec>
inline void 
extendSeed(Seed<ChainedSeed, TConfig> & seed, 
		   TDatabase const & database,
		   TQuery const & query,
		   ExtensionDirection direction,
           Score<TScoreValue, TScoreSpec> const & scoringScheme,
           TScoreValue scoreDropOff,
		   UnGappedXDrop const &)
{
    // For ungapped X-drop extension of Chained Seeds, we extend the
    // first and the last Seed Diagonal.
    scoreDropOff = -scoreDropOff;

    typedef Seed<ChainedSeed, TConfig> TSeed;
    typedef typename Value<TSeed>::Type TSeedDiagonal;
    typedef typename Position<TSeedDiagonal>::Type TPosition;
    typedef typename Size<TSeedDiagonal>::Type TSize;
    
	// Extension to the left
	if (direction == EXTEND_LEFT || direction == EXTEND_BOTH)
    {
        TScoreValue tmpScore = 0;
        TPosition mismatchingSuffixLength = 0;
        TSeedDiagonal & diag = front(seed);
		TPosition posH = beginPositionH(seed);
		TPosition posV = beginPositionV(seed);
        TSize diagonalLength = diag.length;
		while (posH >= 1 && posV >= 1 && tmpScore > scoreDropOff)
        {
            tmpScore += score(scoringScheme, sequenceEntryForScore(scoringScheme, database, posH),
                              sequenceEntryForScore(scoringScheme, query, posV));
            if (database[posH - 1] == query[posV - 1])
            {
                mismatchingSuffixLength = 0;
                if (tmpScore > static_cast<TScoreValue>(0))
                    tmpScore = 0;
            }
            else
            {
                mismatchingSuffixLength += 1;
            }
			--posH;
			--posV;
            ++diagonalLength;
		}
        diag.beginPositionH = posH + mismatchingSuffixLength;
        diag.beginPositionV = posV + mismatchingSuffixLength;
        diag.length = diagonalLength - mismatchingSuffixLength;
	}

	// Extension to the right
	if (direction == EXTEND_RIGHT || direction == EXTEND_BOTH)
    {
        TScoreValue tmpScore = 0;
        TPosition mismatchingSuffixLength = 0;
		TSize lengthH = length(query);
		TSize lengthV = length(database);
        TSeedDiagonal & diag = back(seed);
		TPosition posH = diag.beginPositionH + diag.length;
		TPosition posV = diag.beginPositionV + diag.length;
        TSize diagonalLength = diag.length;
		while (posH < lengthH && posV < lengthV && tmpScore > scoreDropOff)
        {
            tmpScore += score(scoringScheme, sequenceEntryForScore(scoringScheme, database, posH),
                              sequenceEntryForScore(scoringScheme, query, posV));
            if (database[posH] == query[posV])
            {
                mismatchingSuffixLength = 0;
                if (tmpScore > static_cast<TScoreValue>(0))
                    tmpScore = 0;
            }
            else
            {
                mismatchingSuffixLength += 1;
            }
            ++posH;
            ++posV;
            ++diagonalLength;
		}
        diag.length = diagonalLength - mismatchingSuffixLength;
	}

    // TODO(holtgrew): Update score?!
}

// ---------------------------------------------------------------------------
// Function extendSeed                                           [GappedXDrop]
// ---------------------------------------------------------------------------

template<typename TAntiDiag, typename TDropOff, typename TScoreValue>
inline void
_initAntiDiags(TAntiDiag & ,
               TAntiDiag & antiDiag2,
               TAntiDiag & antiDiag3,
               TDropOff dropOff,
               TScoreValue gapCost,
               TScoreValue undefined)
{
	// antiDiagonals will be swaped in while loop BEFORE computation of antiDiag3 entries
	//  -> no initialization of antiDiag1 necessary

    resize(antiDiag2, 1);
    antiDiag2[0] = 0;

    resize(antiDiag3, 2);
    if (-gapCost > dropOff)
    {
        antiDiag3[0] = undefined;
        antiDiag3[1] = undefined;
    }
    else
    {
        antiDiag3[0] = gapCost;
        antiDiag3[1] = gapCost;
    }
}

template<typename TAntiDiag>
inline void
_swapAntiDiags(TAntiDiag & antiDiag1,
               TAntiDiag & antiDiag2,
			   TAntiDiag & antiDiag3)
{
    TAntiDiag temp;
    move(temp, antiDiag1);
    move(antiDiag1, antiDiag2);
    move(antiDiag2, antiDiag3);
    move(antiDiag3, temp);
}

template<typename TAntiDiag, typename TSize, typename TScoreValue>
inline TSize
_initAntiDiag3(TAntiDiag & antiDiag3,
			   TSize offset,
			   TSize maxCol,
			   TSize antiDiagNo,
			   TScoreValue minScore,
               TScoreValue gapCost,
               TScoreValue undefined)
{
	resize(antiDiag3, maxCol + 1 - offset);

    antiDiag3[0] = undefined;
	antiDiag3[maxCol - offset] = undefined;

	if ((int)antiDiagNo * gapCost > minScore)
    {
		if (offset == 0) // init first column
			antiDiag3[0] = antiDiagNo * gapCost;
		if (antiDiagNo - maxCol == 0) // init first row
			antiDiag3[maxCol - offset] = antiDiagNo * gapCost;
	}
	return offset;
}

template<typename TDiagonal, typename TSize>
inline void
_calcExtendedLowerDiag(TDiagonal & lowerDiag,
					   TSize minCol,
					   TSize antiDiagNo)
{
	TSize minRow = antiDiagNo - minCol;
	if ((TDiagonal)minCol - (TDiagonal)minRow < lowerDiag)
		lowerDiag = (TDiagonal)minCol - (TDiagonal)minRow;
}

template<typename TDiagonal, typename TSize>
inline void
_calcExtendedUpperDiag(TDiagonal & upperDiag,
					   TSize maxCol,
					   TSize antiDiagNo)
{
	TSize maxRow = antiDiagNo + 1 - maxCol;
	if ((TDiagonal)maxCol - 1 - (TDiagonal)maxRow > upperDiag)
		upperDiag = maxCol - 1 - maxRow;
}

template<typename TSeed, typename TSize, typename TDiagonal>
inline void
_updateExtendedSeed(TSeed & seed,
					ExtensionDirection direction,
					TSize cols,
					TSize rows,
					TDiagonal lowerDiag,
					TDiagonal upperDiag)
{
	if (direction == EXTEND_LEFT)
    {
		// Set lower and upper diagonals.
		TDiagonal beginDiag = beginDiagonal(seed);
		if (lowerDiagonal(seed) > beginDiag + lowerDiag)
			setLowerDiagonal(seed, beginDiag + lowerDiag);
		if (upperDiagonal(seed) < beginDiag + upperDiag)
			setUpperDiagonal(seed, beginDiag + upperDiag);

		// Set new start position of seed.
		setBeginPositionH(seed, beginPositionH(seed) - rows);
		setBeginPositionV(seed, beginPositionV(seed) - cols);
	} else {  // direction == EXTEND_RIGHT
		// Set new lower and upper diagonals.
		TDiagonal endDiag = endDiagonal(seed);
		if (upperDiagonal(seed) < endDiag - lowerDiag)
			setUpperDiagonal(seed, endDiag - lowerDiag);
		if (lowerDiagonal(seed) > endDiag - upperDiag)
			setLowerDiagonal(seed, endDiag - upperDiag);

		// Set new end position of seed.
        setEndPositionH(seed, endPositionH(seed) + rows);
        setEndPositionV(seed, endPositionV(seed) + cols);
	}
    SEQAN_ASSERT_GEQ(upperDiagonal(seed), lowerDiagonal(seed));
    SEQAN_ASSERT_GEQ(upperDiagonal(seed), beginDiagonal(seed));
    SEQAN_ASSERT_GEQ(upperDiagonal(seed), endDiagonal(seed));
    SEQAN_ASSERT_GEQ(beginDiagonal(seed), lowerDiagonal(seed));
    SEQAN_ASSERT_GEQ(endDiagonal(seed), lowerDiagonal(seed));
}

// Limit score;  In the general case we cannot do this so we simply perform a check on the score mismatch values.
template <typename TScoreValue, typename TScoreSpec, typename TAlphabet>
void
_extendSeedGappedXDropOneDirectionLimitScoreMismatch(Score<TScoreValue, TScoreSpec> & scoringScheme,
                                                     TScoreValue minErrScore,
                                                     TAlphabet * /*tag*/)
{
    // We cannot set a lower limit for the mismatch score since the score might be a scoring matrix such as Blosum62.
    // Instead, we perform a check on the matrix scores.
#if SEQAN_ENABLE_DEBUG
    {
        for (unsigned i = 0; i < valueSize<TAlphabet>(); ++i)
            for (unsigned j = 0; j <= i; ++j)
                SEQAN_ASSERT_GEQ_MSG(score(scoringScheme, TAlphabet(i), TAlphabet(j)), minErrScore,
                                     "Mismatch score too small!, i = %u, j = %u");
    }
#else
    (void)scoringScheme;
    (void)minErrScore;
#endif  // #if SEQAN_ENABLE_DEBUG
}

// In the case of a SimpleScore, however, we can set this.
template <typename TScoreValue, typename TAlphabet>
void
_extendSeedGappedXDropOneDirectionLimitScoreMismatch(Score<TScoreValue, Simple> & scoringScheme,
                                                     TScoreValue minErrScore,
                                                     TAlphabet * /*tag*/)
{
    setScoreMismatch(scoringScheme, std::max(scoreMismatch(scoringScheme), minErrScore));
}

template<typename TConfig, typename TQuerySegment, typename TDatabaseSegment, typename TScoreValue, typename TScoreSpec>
TScoreValue
_extendSeedGappedXDropOneDirection(
        Seed<Simple, TConfig> & seed,
        TQuerySegment const & querySeg,
        TDatabaseSegment const & databaseSeg,
        ExtensionDirection direction,
        Score<TScoreValue, TScoreSpec> scoringScheme,
		TScoreValue scoreDropOff)
{
    typedef typename Size<TQuerySegment>::Type TSize;
	typedef typename Seed<Simple,TConfig>::TDiagonal TDiagonal;

    TSize cols = length(querySeg)+1;
    TSize rows = length(databaseSeg)+1;
	if (rows == 1 || cols == 1)
        return 0;

    TScoreValue len = 2 * _max(cols, rows); // number of antidiagonals
    TScoreValue const minErrScore = minValue<TScoreValue>() / len; // minimal allowed error penalty
    setScoreGap(scoringScheme, _max(scoreGap(scoringScheme), minErrScore));
    typename Value<TQuerySegment>::Type * tag = 0;
    (void)tag;
    _extendSeedGappedXDropOneDirectionLimitScoreMismatch(scoringScheme, minErrScore, tag);

    TScoreValue gapCost = scoreGap(scoringScheme);
    TScoreValue undefined = minValue<TScoreValue>() - gapCost;

    // DP matrix is calculated by anti-diagonals
    String<TScoreValue> antiDiag1;	//smallest anti-diagonal
	String<TScoreValue> antiDiag2;
	String<TScoreValue> antiDiag3;	//current anti-diagonal

	// Indices on anti-diagonals include gap column/gap row:
	//   - decrease indices by 1 for position in query/database segment
