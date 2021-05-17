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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_INDEX_ESA_BASE_H
#define SEQAN_HEADER_INDEX_ESA_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

	// dfs order
	struct Preorder_;
	struct Postorder_;

	template <typename TDfsOrder = Postorder_, typename THideEmptyEdges = True>
	struct VSTreeIteratorTraits {
		typedef TDfsOrder DfsOrder;
		typedef THideEmptyEdges HideEmptyEdges;
	};

/**
.Tag.DFS Order
..summary:Pre/postorder selection for depth-first search.
..cat:Index
..description:
These tags are given to @Function.goNext@ and trigger post-order or pre-order traversal of a suffix tree.
In case of $PreorderEmptyEdges$ and $PostorderEmptyEdges$, the empty edges are also traversed.
..tag.Preorder:Visit the node before its children.
..tag.Postorder:Visit the node after its children.
..tag.PreorderEmptyEdges:Visit the node before its children, visit empty edges.
..tag.PostorderEmptyEdges:Visit the node after its children, visit empty edges.
..include:seqan/index.h
*/
/*!
 * @defgroup DfsOrder DFS Order
 * 
 * @brief Pre/postorder selection for depth-first search.
 * 
 * These tags are given to @link goNext @endlink and trigger post-order or pre-
 * order traversal of a suffix tree. In case of <tt>PreorderEmptyEdges</tt> and
 * <tt>PostorderEmptyEdges</tt>, the empty edges are also traversed.
 * 
 * @tag DfsOrder#Preorder
 * 
 * @brief Visit the node before its children.
 * 
 * @tag DfsOrder#PostorderEmptyEdges
 * 
 * @brief Visit the node after its children, visit empty edges.
 * 
 * @tag DfsOrder#PreorderEmptyEdges
 * 
 * @brief Visit the node before its children, visit empty edges.
 * 
 * @tag DfsOrder#Postorder
 * 
 * @brief Visit the node after its children.
 */
	// predefined iterator traits
	struct Preorder:			VSTreeIteratorTraits<Preorder_,  True> {};
	struct Postorder:			VSTreeIteratorTraits<Postorder_, True> {};
	struct PreorderEmptyEdges:	VSTreeIteratorTraits<Preorder_,  False> {};	// also iterate over
	struct PostorderEmptyEdges:	VSTreeIteratorTraits<Postorder_, False> {};	// empty edges (with $-label)

	// traits for TopDown iterators (w/o ParentLinks) for which postorder/preorder is ignored
	struct HideEmptyEdges:		VSTreeIteratorTraits<Postorder_, True> {};
	struct EmptyEdges:			VSTreeIteratorTraits<Postorder_, False> {};	// empty edges (with $-label)
	
	// MultiMems are more specialized MaxRepeats
	template <typename TSpec = void>
	struct MaxRepeats_;	// base class
	struct MultiMems_;	// subclass tag



	// virtual suffix tree iterators
	template <typename TSpec = void>
	struct VSTree;

/**
.Tag.TopDown
..summary:Tag that specifies a @Spec.VSTree Iterator@ to traverse the virtual string tree from the root towards the leafs.
..cat:Index
..tag.Preorder:Pre-order traversal of the virtual string tree.
..tag.Postorder:Post-order traversal of the virtual string tree.
..tag.ParentLinks:A top down iterator with the possibility to go back up again.
..example
...text:The following example shows how a the @Tag.TopDown@ tag is used.
...file:demos/index/index_begin_atEnd_representative.cpp
...output:
A
AA
ATAA
TA
TAA
TATAA
--------------------------------
AA
ATAA
A
TAA
TATAA
TA

*/
		// top down traversal iterators
		template <typename TSpec = Preorder>
		struct TopDown {};					// starts in the suffix tree root and can go down and go right

			// allows an top-down iterator to go up
			template < typename TSpec = Preorder >
			struct ParentLinks {};			// .. can also go up

/**
.Tag.BottomUp
..summary:Tag that specifies a @Spec.VSTree Iterator@ to traverse the virtual string tree from the root towards the leafs.
..cat:Index
..tag.Postorder:Post-order traversal of the virtual string tree.
..example
...text:The following example shows how a the @Tag.BottomUp@ tag is used.
...file:demos/index/index_begin_atEnd_representative_bottomUp.cpp
...output:
AA
ATAA
A
TAA
TATAA
TA
*/
		// bottom up traversal iterators
		template <typename TSpec = Postorder>
		struct BottomUp {};					// starts in the first node of a depth-first-search and can go next

			struct	SuperMaxRepeats;					// maximal repeat and not part of a longer repeat
			struct	SuperMaxRepeatsFast;
			struct	Mums;								// Maximal Unique Match (unique in every sequence)

			typedef MaxRepeats_<void>		MaxRepeats;	// maximal repeat
			struct	MaxRepeatOccurrences;
			typedef MaxRepeats_<MultiMems_> MultiMems;	// Multiple Maximal Exact Match
			struct	MultiMemOccurences;					// i.e. maximal match over different sequences


/**
.Metafunction.GetVSTreeIteratorTraits:
..cat:Index
..summary:Default behaviour of @Function.goNext@ when no second parameter is given.
..signature:GetVSTreeIteratorTraits<TIterator>::Type
..class:Class.Index
..param.TIterator:A @Spec.VSTree Iterator@.
..returns:$Tag.Postorder$ by default and $Tag.Preorder$ if $TIterator$ is $VSTree<TopDown<ParentLinks<> > >$ or $VSTree<TopDown<ParentLinks<Preorder> > >$.
..include:seqan/index.h
*/
/*!
 * @mfn Index#GetVSTreeIteratorTraits
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Default behaviour of @link goNext @endlink when no second parameter is
 *        given.
 * 
 * @signature GetVSTreeIteratorTraits<TIterator>::Type
 * 
 * @tparam TIterator A @link VSTreeIterator @endlink.
 * 
 * @return TReturn @link DfsOrder#Postorder @endlink by default and @link DfsOrde#Preorder @endlink
 *                 if <tt>TIterator</tt> is <tt>VSTree<TopDown<ParentLinks<> >
 *                 ></tt> or <tt>VSTree<TopDown<ParentLinks<Preorder> > ></tt>.
 */

	template <typename TIterator>
	struct GetVSTreeIteratorTraits:
		DeepestSpec<TIterator> {};

//////////////////////////////////////////////////////////////////////////////

	template <typename TSize>
	struct VertexEsa {
		Pair<TSize> range;			// current SA interval of hits (unique node identifier)
		TSize		parentRight;	// right boundary of parent node's range (allows to go right)

		VertexEsa() : range(0, 0), parentRight(0) {}

		VertexEsa(MinimalCtor):
			range(0,0),
			parentRight(0) {}

		VertexEsa(TSize otherRangeLeft, TSize otherRangeRight, TSize otherParentRight):
			range(Pair<TSize>(otherRangeLeft, otherRangeRight)),
			parentRight(otherParentRight) {}

		VertexEsa(Pair<TSize> const &otherRange, TSize otherParentRight):
			range(otherRange),
			parentRight(otherParentRight) {}

		VertexEsa(VertexEsa const &other):
			range(other.range),
			parentRight(other.parentRight) {}
	};
	
	template <typename TSize>
	inline bool operator==(VertexEsa<TSize> const &a, VertexEsa<TSize> const &b)
	{
		return a.range == b.range;
	}

	template <typename TSize>
	inline bool operator!=(VertexEsa<TSize> const &a, VertexEsa<TSize> const &b)
	{
		return a.range != b.range;
	}

//////////////////////////////////////////////////////////////////////////////
///.Metafunction.VertexDescriptor.param.T.type:Spec.IndexEsa
///.Metafunction.VertexDescriptor.class:Spec.IndexEsa

	template < typename TText, typename TSpec >
	struct VertexDescriptor< Index<TText, IndexEsa<TSpec> > > {
		typedef typename Size< Index<TText, IndexEsa<TSpec> > >::Type TSize;
		typedef VertexEsa<TSize> Type;
	};


//////////////////////////////////////////////////////////////////////////////
// needful forward declarations

	struct ArrayGaps_;
    typedef Tag<ArrayGaps_> ArrayGaps;

	template <typename TSource, typename TSpec>
	class Align;


//////////////////////////////////////////////////////////////////////////////
// ESA fibres

/**
.Tag.ESA Index Fibres
..summary:Tag to select a specific fibre (e.g. table, object, ...) of an @Spec.IndexEsa.ESA@ index.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of an Enhanced Suffix Array based @Spec.IndexEsa.Index@.
..cat:Index

..tag.EsaText:The original text the index should be based on.

..tag.EsaRawText:The raw text the index is really based on.
...remarks:$EsaText$ and $EsaRawText$ fibres are equal by default.
They differ if the index text is a set of strings. Then, raw text is the concatenation of all strings in this set.

..tag.EsaSA:The suffix array.
...remarks:The suffix array contains the indices of all suffices of $EsaRawText$ in lexicographical order.
...remarks:@Metafunction.Fibre@ returns a @Class.String@ over the alphabet of the @Metafunction.SAValue@ of $TIndex$.

..tag.EsaLcp:The lcp table.
...remarks:The lcp table contains the lcp-value of two adjacent suffices in the suffix array $EsaSA$.
...remarks:@Metafunction.Fibre@ returns a @Class.String@ over the alphabet of a size type.

..tag.EsaChildtab:The child table.
...remarks:The child table contains structural information of the suffix tree (see Abhouelda et al.).
...remarks:@Metafunction.Fibre@ returns a @Class.String@ over the alphabet of a size type.

..tag.EsaBwt:The Burrows-Wheeler table.
...remarks:The Burrows-Wheeler table contains the Burrows-Wheeler transformation of $EsaRawText$.
The entries are the characters left of the corresponding suffix in the suffix array $EsaSA$.
...remarks:@Metafunction.Fibre@ returns the same type for $EsaRawText$ and for $EsaBwt$.

..see:Metafunction.Fibre
..see:Function.getFibre
..see:Spec.IndexEsa
..include:seqan/index.h
*/

/*!
 * @defgroup IndexEsaFibres Index Esa Fibres
 * 
 * @brief Tag to select a specific fibre (e.g. table, object, ...) of an @link
 *        IndexEsa @endlink index.
 * 
 * @section Remarks
 * 
 * These tags can be used to get @link Fibre.Fibres @endlink of an Enhanced
 * Suffix Array based @link IndexEsa.Index @endlink.
 * 
 * @see Fibre
 * @see getFibre
 * @see IndexEsa
 * 
 * @tag IndexEsaFibres#EsaSA
 * 
 * @headerfile seqan/index.h
 *
 * @brief The suffix array.
 * 
 * @section Remarks
 * 
 * The suffix array contains the indices of all suffices of <tt>EsaRawText</tt>
 * in lexicographical order.
 * 
 * @link Fibre @endlink returns a @link String @endlink over the alphabet of the
 * @link SAValue @endlink of <tt>TIndex</tt>.
 * 
 * @tag IndexEsaFibres#EsaChildtab
 * 
 * @headerfile seqan/index.h
 *
 * @brief The child table.
 * 
 * @section Remarks
 * 
 * The child table contains structural information of the suffix tree (see
 * Abhouelda et al.).
 * 
 * @link Fibre @endlink returns a @link String @endlink over the alphabet of a
 * size type.
 * 
 * @tag IndexEsaFibres#EsaRawText
 * 
 * @headerfile seqan/index.h
 *
 * @brief The raw text the index is really based on.
 * 
 * @section Remarks
 * 
 * <tt>EsaText</tt> and <tt>EsaRawText</tt> fibres are equal by default. They
 * differ if the index text is a set of strings. Then, raw text is the
 * concatenation of all strings in this set.
 * 
 * @tag IndexEsaFibres#EsaText
 * 
 * @headerfile seqan/index.h
 *
 * @brief The original text the index should be based on.
 * 
 * @tag IndexEsaFibres#EsaBwt
 * 
 * @headerfile seqan/index.h
 *
 * @brief The Burrows-Wheeler table.
 * 
 * @section Remarks
 * 
 * The Burrows-Wheeler table contains the Burrows-Wheeler transformation of
 * <tt>EsaRawText</tt>. The entries are the characters left of the corresponding
 * suffix in the suffix array <tt>EsaSA</tt>.
 * 
 * @link Fibre @endlink returns the same type for <tt>EsaRawText</tt> and for
 * <tt>EsaBwt</tt>.
 * 
 * @tag IndexEsaFibres#EsaLcp
 * 
 * @headerfile seqan/index.h
 *
 * @brief The lcp table.
 * 
 * @section Remarks
 * 
 * The lcp table contains the lcp-value of two adjacent suffices in the suffix
 * array <tt>EsaSA</tt>.
 * 
 * @link Fibre @endlink returns a @link String @endlink over the alphabet of a
 * size type.
 */

///.Metafunction.Fibre.param.TSpec.type:Tag.ESA Index Fibres

	typedef FibreText		EsaText;
	typedef FibreRawText	EsaRawText;
	typedef FibreSA         EsaSA;
	typedef FibreRawSA		EsaRawSA;
	typedef FibreSae		EsaSae;
	typedef FibreLcp		EsaLcp;
	typedef FibreLcpe		EsaLcpe;
	typedef FibreChildtab	EsaChildtab;
	typedef FibreBwt		EsaBwt;


//////////////////////////////////////////////////////////////////////////////
// ESA index

/**
.Spec.IndexEsa:
..summary:The enhanced suffix array index is very fast index, requiring more memory than other indices.
In addition to the suffix array an lcp (longest common prefix) table and a child table (containing structural 
information of the suffix tree) are provided.
..cat:Index
..general:Class.Index
..signature:Index<TText, IndexEsa<> >
..param.TText:The text type.
...type:Class.String
...type:Class.StringSet
..remarks:The main fibres (see @Class.Index@ and @Metafunction.Fibre@) of this index are a suffix array
(see @Tag.ESA Index Fibres.EsaSA@), a lcp table (see @Tag.ESA Index Fibres.EsaLcp@) and a child table (see @Tag.ESA Index Fibres.EsaChildtab@).
..remarks:This index can be accessed as a Suffix Tree using the @Spec.VSTree Iterator@ classes.  ..include:seqan/index.h
*/
/*!
 * @class IndexEsa
 * 
 * @extends Index
 * 
 * @headerfile seqan/index.h
 * 
 * @brief An index based on an enhanced suffix array.
 * 
 * @signature Index<TText, IndexEsa<> >
 * 
 * @tparam TText The text type. Types: String
 * 
 * @section Remarks
 * 
 * The fibres (see @link Index @endlink and @link Fibre @endlink) of this index
 * are a suffix array (see @link ESA Index Fibres.EsaSA @endlink), a lcp table
 * (see @link ESA Index Fibres.EsaLcp @endlink), etc.
 * 
 * This index can be accessed as a Suffix Tree using the @link VSTree Iterator
 * @endlink classes.
 * 
 * @see ESA Index Fibres
 */

/*
	already defined in index_base.h

	template <typename TSpec = void>
	struct IndexEsa;
*/

	template < typename TText, typename TSpec >
	class Index<TText, IndexEsa<TSpec> > {
	public:
		Holder<typename Fibre<Index, EsaText>::Type>	text;
		typename Fibre<Index, EsaSA>::Type				sa;			// suffix array 
		typename Fibre<Index, EsaLcp>::Type			lcp;		// longest-common-prefix table
		typename Fibre<Index, EsaLcpe>::Type			lcpe;		// extended lcp table
		typename Fibre<Index, EsaChildtab>::Type		childtab;	// child table (tree topology)
		typename Fibre<Index, EsaBwt>::Type			bwt;		// burrows-wheeler table
		typename Cargo<Index>::Type						cargo;		// user-defined cargo

		Index() {}

		Index(Index &other):
			text(other.text),
			sa(other.sa),
			lcp(other.lcp),
			lcpe(other.lcpe),
			childtab(other.childtab),
			bwt(other.bwt),
			cargo(other.cargo) {}

		Index(Index const &other):
			text(other.text),
			sa(other.sa),
			lcp(other.lcp),
			lcpe(other.lcpe),
			childtab(other.childtab),
			bwt(other.bwt),
			cargo(other.cargo) {}

		template <typename TText_>
		Index(TText_ &_text):
			text(_text) {}

		template <typename TText_>
		Index(TText_ const &_text):
			text(_text) {}
	};

//////////////////////////////////////////////////////////////////////////////

	template < typename TText, typename TSpec >
	void _indexRequireTopDownIteration(Index<TText, IndexEsa<TSpec> > &index) 
	{
		indexRequire(index, EsaSA());
		indexRequire(index, EsaLcp());
		indexRequire(index, EsaChildtab());
	}

	template < typename TText, typename TSpec >
	void _indexRequireBottomUpIteration(Index<TText, IndexEsa<TSpec> > &index) 
	{
		indexRequire(index, EsaSA());
		indexRequire(index, EsaLcp());
	}

//////////////////////////////////////////////////////////////////////////////
///.Function.clear.param.object.type:Class.Index

	template <typename TText, typename TSpec>
	inline void clear(Index<TText, IndexEsa<TSpec> > &index) {
		clear(getFibre(index, EsaSA()));
		clear(getFibre(index, EsaLcp()));
		clear(getFibre(index, EsaLcpe()));
		clear(getFibre(index, EsaChildtab()));
		clear(getFibre(index, EsaBwt()));
	}

// ----------------------------------------------------------------------------
// Function open
// ----------------------------------------------------------------------------

	template < typename TObject, typename TSpec >
	inline bool open(
		Index< TObject, IndexEsa<TSpec> > &index, 
		const char *fileName,
		int openMode)
	{
		String<char> name;

		name = fileName;	append(name, ".txt");
		if ((!open(getFibre(index, EsaText()), toCString(name), openMode)) && 
			(!open(getFibre(index, EsaText()), fileName, openMode))) return false;

		name = fileName;	append(name, ".sa");
        if (!open(getFibre(index, EsaSA()), toCString(name), openMode)) return false;

		name = fileName;	append(name, ".lcp");
        if (!open(getFibre(index, EsaLcp()), toCString(name), openMode)) return false;

		name = fileName;	append(name, ".child");
        if (!open(getFibre(index, EsaChildtab()), toCString(name), openMode)) return false;

		name = fileName;	append(name, ".bwt");
        if (!open(getFibre(index, EsaBwt()), toCString(name), openMode)) return false;

		return true;
	}
	template < typename TObject, typename TSpec >
	inline bool open(
		Index< TObject, IndexEsa<TSpec> > &index, 
		const char *fileName) 
	{
		return open(index, fileName, DefaultOpenMode<Index< TObject, IndexEsa<TSpec> > >::VALUE);
	}

// ----------------------------------------------------------------------------
// Function save
// ----------------------------------------------------------------------------

	template < typename TObject, typename TSpec >
	inline bool save(
		Index< TObject, IndexEsa<TSpec> > &index, 
		const char *fileName,
		int openMode)
	{
		String<char> name;

		name = fileName;	append(name, ".txt");	
		if ((!save(getFibre(index, EsaText()), toCString(name), openMode)) && 
			(!save(getFibre(index, EsaText()), fileName, openMode))) return false;

		name = fileName;	append(name, ".sa");
        if (!save(getFibre(index, EsaSA()), toCString(name), openMode)) return false;

		name = fileName;	append(name, ".lcp");
        if (!save(getFibre(index, EsaLcp()), toCString(name), openMode)) return false;

		name = fileName;	append(name, ".child");
        if (!save(getFibre(index, EsaChildtab()), toCString(name), openMode)) return false;

		name = fileName;	append(name, ".bwt");
        if (!save(getFibre(index, EsaBwt()), toCString(name), openMode)) return false;

		return true;
	}
	template < typename TObject, typename TSpec >
	inline bool save(
		Index< TObject, IndexEsa<TSpec> > &index, 
		const char *fileName)
	{
		return save(index, fileName, DefaultOpenMode<Index< TObject, IndexEsa<TSpec> > >::VALUE);
	}

}

#endif
