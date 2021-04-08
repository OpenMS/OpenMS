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

#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_DISTANCE_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_DISTANCE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Distance matrix calculation
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Distance Calculation:
..cat:Alignments
..summary:A tag to specify how to calculate distance matrices.
..include:seqan/graph_msa.h
*/

/**
.Tag.Distance Calculation.value.LibraryDistance:
	Using the library itself and heaviest common subsequence to determine a distance matrix
..include:seqan/graph_msa.h
*/
struct LibraryDistance_;
typedef Tag<LibraryDistance_> const LibraryDistance;


/**
.Tag.Distance Calculation.value.KmerDistance:
	Using a simple kmer count to determine a distance matrix
..include:seqan/graph_msa.h
*/
struct KmerDistance_;
typedef Tag<KmerDistance_> const KmerDistance;






//////////////////////////////////////////////////////////////////////////////
// LibraryDistance
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TMatrix>
inline void
getDistanceMatrix(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				  TMatrix& distanceMatrix,
				  LibraryDistance)
{
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	//typedef typename Id<TGraph>::Type TId;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	//typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	//typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Value<TMatrix>::Type TValue;

	// Initialization
	clear(distanceMatrix);
	TStringSet& str = stringSet(g);
	TSize nseq = length(str);
	resize(distanceMatrix, nseq * nseq, 0);

	// All pairwise alignments
	typedef String<String<TVertexDescriptor> > TSegmentString;
	TValue maxScore = 0;
	for(TSize i=0; i<nseq; ++i) {
		TSegmentString seq1;
		TSize len1 = length(str[i]);
		_buildLeafString(g, i, seq1);
		for(TSize j=i+1; j<nseq; ++j) {
			// Align the 2 strings
			TSegmentString seq2;
			TSize len2 = length(str[j]);
			_buildLeafString(g, j, seq2);
			TSegmentString alignSeq;
			TValue score = heaviestCommonSubsequence(g,seq1,seq2,alignSeq);
			
			// Normalize by distance
			if (len1 > len2) score /= len1;
			else score /= len2;
			if (score > maxScore) maxScore = score;
			
			// Remember the value
			distanceMatrix[i*nseq+j] = score;
		}
	}

	// Normalize values
	for(TSize i=0; i<nseq; ++i) 
		for(TSize j=i+1; j<nseq; ++j) 
			distanceMatrix[i*nseq+j] = SEQAN_DISTANCE_UNITY - ((distanceMatrix[i*nseq+j] * SEQAN_DISTANCE_UNITY) / maxScore );
}


//////////////////////////////////////////////////////////////////////////////
// KmerDistance
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


template<typename TStringSet, typename TCargo, typename TSpec, typename TMatrix, typename TSize, typename TAlphabet>
inline void
getDistanceMatrix(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				  TMatrix& distanceMatrix,
				  TSize ktup,
				  TAlphabet,
				  KmerDistance)
{
	//typedef typename Value<TMatrix>::Type TValue;
	typedef typename Iterator<TMatrix, Standard>::Type TMatrixIterator;

	getKmerSimilarityMatrix(stringSet(g), distanceMatrix, ktup, TAlphabet());

	// Similarity to distance conversion
	TMatrixIterator matIt = begin(distanceMatrix, Standard());
	TMatrixIterator endMatIt = end(distanceMatrix, Standard());
	for(;matIt != endMatIt;++matIt)
		*matIt = SEQAN_DISTANCE_UNITY - (*matIt);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TMatrix, typename TSize>
inline void 
getDistanceMatrix(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				  TMatrix& distanceMatrix,
				  TSize ktup,
				  KmerDistance)
{
	SEQAN_CHECKPOINT
	getDistanceMatrix(g, distanceMatrix, ktup, typename Value<typename Value<TStringSet>::Type>::Type(), KmerDistance() );
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TMatrix>
inline void 
getDistanceMatrix(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				  TMatrix& distanceMatrix,
				  KmerDistance)
{
	SEQAN_CHECKPOINT
	getDistanceMatrix(g, distanceMatrix, 3, KmerDistance() );
}


//////////////////////////////////////////////////////////////////////////////


/**
.Function.getDistanceMatrix
..class:Spec.Alignment Graph
..summary:Computes a pairwise distance matrix from an alignment graph.
..cat:Graph
..signature:
getDistanceMatrix(graph, mat [, tag])
getDistanceMatrix(graph, mat [, ktup] [, alphabet], KmerDistance)
..param.graph:An alignment graph containing the sequences and possible alignment edges.
...type:Spec.Alignment Graph
..param.mat:Out-parameter:Pairwise distance matrix.
...type:Class.String
..param.ktup:Length of k-mers.
...remarks:For KmerDistance the length of the k-mers.
..param.alphabet:Alphabet
...remarks:For KmerDistance the alphabet to use for k-mer counting (e.g., compressed alphabets).
..param.tag:Distance tag
...type:Tag.Distance Calculation
...remarks:Possible values are LibraryDistance or KmerDistance.
...default:KmerDistance
..returns:void
..include:seqan/graph_msa.h
*/
template<typename TStringSet, typename TCargo, typename TSpec, typename TMatrix>
inline void 
getDistanceMatrix(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				  TMatrix& distanceMatrix)
{
	SEQAN_CHECKPOINT
	getDistanceMatrix(g, distanceMatrix, KmerDistance() );
}



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
