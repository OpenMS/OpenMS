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

#ifndef SEQAN_HEADER_INDEX_BASE_H
#define SEQAN_HEADER_INDEX_BASE_H

//#define SEQAN_TEST_INDEX

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// needful forward declarations

	// suffix array construction specs
	struct Skew3;
	struct Skew7;
	struct LarssonSadakane;
	struct ManberMyers;
	struct SAQSort;
	struct QGramAlg;

	// lcp table construction algorithms
	struct Kasai;
	struct KasaiOriginal;	// original, but more space-consuming algorithm

	// enhanced suffix array construction algorithms
	struct Childtab;
	struct Bwt;

/*
.Tag.Index Find Algorithm
..summary:Tag to specify the index search algorithm.
..remarks:These tag can be used to specify the @Function.find@ algorithm 
for @Class.Index@ based substring searches.
..cat:Index

..tag.EsaFindMlr:Binary search with mlr-heuristic.
...remarks:Exact string matching using a suffix array binary search with the mlr-heuristic.

..tag.EsaFindLcpe:Binary search using lcp values.
...remarks:Exact string matching using a suffix array binary search and a lcp-interval tree.

..tag.FinderSTree:Suffix tree search.
...remarks:Exact string matching using a suffix tree.

..see:Class.Finder
..see:Spec.IndexEsa
..see:Spec.IndexQGram
..include:seqan/index.h
*/

/*!
 * @defgroup IndexFindAlgorithm Index Find Algorithm
 * 
 * @brief Tag to specify the index search algorithm.
 * 
 * @section Remarks
 * 
 * These tags can be used to specify the @link find @endlink algorithm for @link
 * Index @endlink based substring searches.
 * 
 * @see Finder
 * 
 * @tag IndexFindAlgorithm#FinderSTree
 * 
 * @brief Suffix tree search.
 * 
 * @section Remarks
 * 
 * Exact string matching using a suffix tree.
 * 
 * @tag IndexFindAlgorithm#PizzaChiliFinder
 * 
 * @brief Finds an occurrence in a @link Pizza & Chili Index @endlink index.
 * 
 * @section Remarks
 * 
 * The actual algorithm used for searching depends on the @link Pizza & Chili
 * Index Tags @endlink used.
 * 
 * @tag IndexFindAlgorithm#QGramFindLookup
 * 
 * @brief q-gram search. Finds q-grams in a @link IndexQGram @endlink index
 *        using the hash table.
 * 
 * @tag IndexFindAlgorithm#EsaFindLcpe
 * 
 * @brief Binary search using lcp values.
 * 
 * @section Remarks
 * 
 * Exact string matching using a suffix array binary search and a lcp-interval
 * tree.
 * 
 * @tag IndexFindAlgorithm#EsaFindMlr
 * 
 * @brief Binary search with mlr-heuristic.
 * 
 * @section Remarks
 * 
 * Exact string matching using a suffix array binary search with the mlr-
 * heuristic.
 */
	// finder tags
    struct FinderMlr_;     // simple Suffix Array finder with mlr-heuristic
    struct FinderLcpe_;    // Suffix Array finder using an enhanced LCP-Table
    struct FinderSTree_;    // Suffix Array finder using an enhanced LCP-Table

    typedef Tag<FinderMlr_> const EsaFindMlr;
    typedef Tag<FinderLcpe_> const EsaFindLcpe;
    typedef Tag<FinderSTree_> const FinderSTree;

	template <typename TSpec = void>
	struct IndexEsa {};


// ----------------------------------------------------------------------------
// Metafunction DefaultIndexSpec
// ----------------------------------------------------------------------------
/**
.Metafunction.DefaultIndexSpec:
..cat:Index
..summary:Default @Class.Index@ specialization type.
..signature:DefaultIndexSpec<TText>::Type
..class:Class.Index
..param.TText:The given text type.
..returns:Can be @Spec.IndexEsa@ or $IndexQGram$, etc.
..remarks:Currently @Spec.IndexEsa@ is default if $TText$ is a @Class.String@.
..include:seqan/index.h
*/
/*!
 * @mfn Index#DefaultIndexSpec
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Default @link Index @endlink specialization type.
 * 
 * @signature DefaultIndexSpec<TText>::Type
 * 
 * @tparam TText The given text type.
 * 
 * @return TReturn Currently the return type is @link IndexEsa @endlink.
 * 
 * @section Remarks
 */
    template < typename TObject >
    struct DefaultIndexSpec {
        typedef IndexEsa<> Type;
    };

// ----------------------------------------------------------------------------
// Metafunction DefaultIndexStringSpec
// ----------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////
/**
.Metafunction.DefaultIndexStringSpec:
..cat:Index
..summary:Default @Class.String@ specialization type of the @Metafunction.Fibre@ of an @Class.Index@. 
..signature:DefaultIndexStringSpec<TIndex>::Type
..class:Class.Index
..param.TIndex:An @Class.Index@ Type.
..returns:If the underlying text is a @Class.String@ or a set of Strings (see @Class.StringSet@) the String's spec. type is returned.
..remarks:Most of the @Class.Index@ fibres are strings. The @Class.String@ specialization type is chosen by this meta-function.
..include:seqan/index.h
*/


//TODO(singer): Does not belong here but to the text concept
/*!
 * @mfn Index#DefaultIndexStringSpec
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Default @link String @endlink specialization type of the @link Fibre
 *        @endlink of an @link Index @endlink.
 * 
 * @signature DefaultIndexStringSpec<TIndex>::Type
 * 
 * @tparam TIndex An @link Index @endlink Type.
 * 
 * @return TReturn If the underlying text is a @link String @endlink or a set of
 *                 Strings (see @link StringSet @endlink) the String's spec.
 *                 type is returned.
 * 
 * @section Remarks
 * 
 * Most of the @link Index @endlink fibres are strings. The @link String
 * @endlink specialization type is chosen by this meta-function.
 */	
    
    // default which should actually never been used
    template < typename TIndex >
    struct DefaultIndexStringSpec {
        typedef Alloc<> Type;
    };

    template < typename TValue, typename TSpec >
    struct DefaultIndexStringSpec< String<TValue, External<TSpec> > > {
        typedef External<TSpec> Type;
    };

	template < typename TString, typename TSpec >
	struct DefaultIndexStringSpec< StringSet<TString, TSpec> >:
		DefaultIndexStringSpec<TString> {};

//////////////////////////////////////////////////////////////////////////////
/**
.Class.Index:
..summary:Contains preprocessing data of a fixed text. In combination with a @Class.Finder@ or a @Spec.VSTree Iterator@ it allows fast dictionary look-up and advanced computations.
..cat:Index
..signature:Index<TText[, TSpec]>
..param.TText:The text type.
...type:Class.String
...type:Class.StringSet
...metafunction:Metafunction.Host
..param.TSpec:The index type.
...type:Spec.IndexEsa
...type:Spec.IndexWotd
...type:Spec.IndexQGram
...type:Spec.FMIndex
...default:The result of @Metafunction.DefaultIndexSpec@: @Spec.IndexEsa@
...metafunction:Metafunction.Spec
..remarks:Indices allow fast dictionary look-ups and other advanced computations because they contain pre-computated
information about the underlying text. These information are stored in so called fibres (see @Metafunction.Fibre@).
..remarks:In order to search for a pattern one can use a @Class.Finder@ or a @Spec.VSTree Iterator@. The @Class.Finder@ 
is especially useful when searching for exact matches while the @Spec.VSTree Iterator@ allows to iterate an index as if
traversing a tree/trie.
..remarks:These fibres are created on demand depending on the requirements of an algorithm.
..include:seqan/index.h
..example
...text:The following code shows how to search for exact matches between the reference "tobeornottobe" and the
pattern "to" with the means of a Finder.
...file:demos/index/index_finder.cpp
...output:Hit at position: 9
Hit at position: 0
...text:This code shows how an index can be used with iterators to achieve a pre-order tree like traversal
in DFS of the text "tobeornottobe". In order to do so a Top-Down History iterator is used.
...file:demos/index/index_iterator.cpp
...output:

be
beornottobe
e
eornottobe
nottobe
o
obe
obeornottobe
ornottobe
ottobe
rnottobe
t
tobe
tobeornottobe
ttobe
...text:
Note that you can also use specialised iterators such as:
...code:Iterator<TIndex, TopDown<ParentLinks<PreOrder> > >::Type 
...text:or
...code:Iterator<TIndex, TopDown<ParentLinks<PostOrder> > >::Type
...text:You can achieve a post-order traversal like this:
...snippet:demos/index/index_iterator_short.cpp|iteration
...output:beornottobe
be
eornottobe
e
nottobe
obeornottobe
obe
ornottobe
ottobe
o
rnottobe
tobeornottobe
tobe
ttobe
t

*/

///.Function.setHaystack.param.haystack.type:Class.Index

template < 
        typename TObject, 
        typename TSpec = typename DefaultIndexSpec<TObject>::Type > 
	class Index;

/*!
 * @class Index
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Indices are data structures which contain preprocessing data of a
 *        fixed text and in doing so allow fast dictionary look-up and advanced
 *        computations. There are various indices implemented in SeqAn, which
 *        can be used for different scenarious.
 * 
 * @signature Index<TText[, TSpec]>
 * 
 * @tparam TSpec The index type. Default: The result of @link DefaultIndexSpec @endlink
 * @tparam TText The text type. Types: String, StringSet
 * 
 * @section Remarks
 * 
 * An index contains various arrays or objects, also called fibres (see @link Index#Fibre @endlink).
 * 
 * These fibres are created on demand depending on the requirements of an algorithm.
 *
 * The following are common index fibres:
 *
 * <table border="1">
 * <tr>
 *   <td>FibreSA</td>
 *   <td>Suffix array fibre.  A string of @link TextConcept#SAValue SAValue<TText>::Type @endlink.</td>
 * </tr>
 * </table>
 *
 * The list of fibres is available in @link IndexFibres @endlink
 *
 * @see IndexFibres
 */ 

	template <typename TObject, typename TSpec>
	struct Host< Index<TObject, TSpec> > {
		typedef TObject Type;
	};

	template <typename TObject, typename TSpec>
	struct Spec< Index<TObject, TSpec> > {
		typedef TSpec Type;
	};
/**
.Metafunction.Fibre:
..summary:Type of a specific container member (fibre).
..signature:Fibre<TContainer, TSpec>::Type
..class:Class.Index
..cat:Index
..param.TContainer:The container type.
...type:Class.Index
..param.TSpec:Tag to specify the fibre.
..returns:Fibre type.
..remarks:Some containers, such as @Class.Index@, can be seen as a bundle consisting of various fibres. Because not every table is a fibre we did not call them tables, however, in many cases one can think of fibres as tables. The fibre interface was designed to unify the access to the members of the different fibres.
To get a reference or the type of a specific fibre use @Function.getFibre@ or @Metafunction.Fibre@.		
..remarks:A @Metafunction.Fibre@ does not need to be a real container. It can also be a view (see @Tag.ESA Index Fibres.EsaRawText@).
..include:seqan/index.h
*/

/*!
 * @mfn Index#Fibre
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Type of a specific container member (fibre).
 * 
 * @signature Fibre<TContainer, TSpec>::Type
 * 
 * @tparam TSpec Tag to specify the fibre. Types: @link IndexEsaFibres @endlink, @link WaveletTreeFibres @endlink, 
 *               @link RightArrayBinaryTreeFibres @endlink, @link SentinelRankDictionaryFibres @endlink
 * @tparam TContainer The container type. Types: Index, RankDictionary,
 *                    RightArrayBinaryTree
 * 
 * @return Type Fibre type.
 * 
 * @section Remarks
 * 
 * Some containers, such as @link Index @endlink, can be seen as a bundle
 * consisting of various fibres. Because not every table is a fibre we did not
 * call them tables, however, in many cases one can think of fibres as tables.
 * The fibre interface was designed to unify the access to the members of the
 * different fibres. To get a reference or the type of a specific fibre use
 * @link getFibre @endlink or @link Fibre @endlink.
 * 
 * A @link Fibre @endlink does not need to be a real container. It can also be a
 * view (see @link ESA Index Fibres.EsaRawText @endlink).
 * 
 * @see Index#getFibre
 */

	// meta function to get the type of a bundle fibre
	template < typename TIndex, typename TSpec >
	struct Fibre {
		typedef String< typename Size<TIndex>::Type > Type;
	};

	template < typename TIndex, typename TSpec >
	struct Fibre<TIndex const, TSpec> {
		typedef typename Fibre<TIndex, TSpec>::Type const Type;
	};

	struct FibreRecord {
		unsigned	id;
		void*		ptr;
		bool		owner;
	};

	// less function to search in sorted list for fibre id
	struct FibreLess: public ::binary_function<FibreRecord, unsigned, bool>
	{	// functor for operator>
		inline bool operator()(FibreRecord const & _Left, unsigned const Right_) const
		{	// apply operator> to operands
			return (_Left.id < Right_);
		}
	};

//////////////////////////////////////////////////////////////////////////////
/**
.Metafunction.DefaultIndexCreator:
..cat:Index
..summary:Default algorithm to create a demanded and not yet existing @Metafunction.Fibre@.
..signature:DefaultIndexCreator<TIndex, TFibre>::Type
..class:Class.Index
..param.TIndex:An @Class.Index@ Type.
..param.TFibre:A tag specifying the fibre (e.g. @Tag.ESA Index Fibres.EsaSA@).
..returns:A tag specifying the default algorithm to create the fibre with.
..include:seqan/index.h
*/

/*!
 * @mfn Index#DefaultIndexCreator
 * 
 * @headerfile seqan/index.h
 * 
 * @deprecated advanced
 * 
 * @brief Default algorithm to create a demanded and not yet existing @link
 *        Fibre @endlink.
 *
 * @signature DefaultIndexCreator<TIndex, TFibre>::Type
 * 
 * @tparam TIndex An @link Index @endlink Type.
 * @tparam TFibre A tag specifying the fibre (e.g. @link ESA Index Fibres.EsaSA
 *                @endlink).
 * 
 * @return TReturn A tag specifying the default algorithm to create the fibre
 *                 with.    // standard algorithm for indices creation
 */

    template < typename TIndex, typename TFibre >
	struct DefaultIndexCreator {
		typedef Default Type;
	};


//////////////////////////////////////////////////////////////////////////////

	template < 
		typename TSA,
		typename TText,
        typename TAlgSpec >
    struct SACreatorRandomAccess_
    {
        typedef typename AllowsFastRandomAccess<TSA>::Type   TRandomSA;
        typedef typename AllowsFastRandomAccess<TText>::Type TRandomText;
        typedef typename And<TRandomText,TRandomSA>::Type Type;
    };

	template < 
        typename TLCP,
		typename TText,
		typename TSA,
        typename TAlgSpec >
    struct LcpCreatorRandomAccess_
    {
        typedef typename AllowsFastRandomAccess<TText>::Type TRandomText;
        typedef typename AllowsFastRandomAccess<TLCP>::Type  TRandomLCP;
        typedef typename AllowsFastRandomAccess<TSA>::Type   TRandomSA;
        typedef typename And<TRandomLCP, typename And<TRandomText,TRandomSA>::Type>::Type Type;
    };


//////////////////////////////////////////////////////////////////////////////
/*
	.Class.Bundle:
	..summary:General purpose container of various members.
	..signature:Bundle<TValue, TSize>
	..param.TValue:The value type, that is the type of the items/characters stored in the string.
	...remarks:Use @Metafunction.Value@ to get the value type for a given class.
	..param.TSpec:The specializing type.
	...default:$Alloc<>$, see @Spec.Alloc String@.
*/
/*
	template < typename TSpec = void >
	truct Bundle {
		typedef ::std::vector<FibreRecord>	TFibreRecords;
		TFibreRecords						fibres;
	};

	template < typename TBundleSpec, typename TFibreSpec >
	inline FibreRecord& getRecord(Bundle<TBundleSpec> &bundle, TFibreSpec const) {
		unsigned id = (unsigned)ClassIdentifier_<TFibreSpec>::getID();

		typename Bundle<TBundleSpec>::TFibreRecords::iterator first = lower_bound(bundle.fibres.begin(), bundle.fibres.end(), id, FibreLess());
		if (!first->id != id) {
			FibreRecord rec;
			rec.id = id;
			rec.ptr = NULL;
			rec.owner = true;
			bundle.fibres.insert(first, rec);
		} else
			return *first;
	}

	template < typename TBundleSpec, typename TFibreSpec >
	inline typename Fibre<Bundle<TBundleSpec>, TFibreSpec>::Type & getFibre(Bundle<TBundleSpec> &bundle, TFibreSpec const) {
		typedef typename Fibre<Bundle<TBundleSpec>, TFibreSpec>::Type Type;
		unsigned id = (unsigned)ClassIdentifier_<TFibreSpec>::getID();

		FibreRecord &rec = getRecord(bundle, TFibreSpec());
		if (!rec.ptr)
			rec.ptr = new Type();
		return *reinterpret_cast<Type*>(rec.ptr);
	}

	template < typename TBundleSpec, typename TFibreSpec >
	inline typename Fibre<Bundle<TBundleSpec>, TFibreSpec>::Type const & getFibre(Bundle<TBundleSpec> const &bundle, TFibreSpec const) {
		typedef typename Fibre<Bundle<TBundleSpec>, TFibreSpec>::Type Type;
		unsigned id = (unsigned)ClassIdentifier_<TFibreSpec>::getID();

		FibreRecord &rec = getRecord(bundle, TFibreSpec());
		return *reinterpret_cast<Type*>(rec.ptr);
	}
*/

//////////////////////////////////////////////////////////////////////////////
// various fibre specs for enhanced suffix arrays

	struct FibreText_;		// Original text. Can be a String or a StringSet
	struct FibreRawText_;	// Concatenation of the strings above
	struct FibreSA_;		// suffix array (of raw text with virtual $-delimiters) with Pair entries
	struct FibreRawSA_;	// suffix array with integer entries
	struct FibreSae_;		// suffix array reordered in a b-tree
	struct FibreLcp_;		// lcp table of raw text
	struct FibreLcpe_;		// lcp interval tree
	struct FibreChildtab_;	// childtab (Kurtz et al.) of raw text
	struct FibreBwt_;		// burrows wheeler table of raw text

	typedef Tag<FibreText_> const		FibreText;
	typedef Tag<FibreRawText_> const	FibreRawText;
	typedef Tag<FibreSA_> const         FibreSA;
	typedef Tag<FibreRawSA_> const		FibreRawSA;
	typedef Tag<FibreSae_> const		FibreSae;
	typedef Tag<FibreLcp_> const		FibreLcp;
	typedef Tag<FibreLcpe_> const		FibreLcpe;
	typedef Tag<FibreChildtab_> const	FibreChildtab;
	typedef Tag<FibreBwt_> const		FibreBwt;

//////////////////////////////////////////////////////////////////////////////
/**
.Metafunction.SAValue:
..cat:Index
..summary:The default alphabet type of a suffix array, i.e. the type to store a position of a string or string set.
..signature:SAValue<TObject>::Type
..class:Class.Index
..param.TObject:A string, string set, or index type.
...type:Class.String
...type:Class.StringSet
...type:Class.Index
..returns:A type to store a position. 
...text:If $TObject$ is a @Class.String@, it is a single integer value. By default this is the @Metafunction.Size@ type of $TObject$.
...text:If $TObject$ is a @Class.StringSet@, it could be a single integer too (called global position, see @Spec.ConcatDirect@) or a @Class.Pair@ (called local position, see @Spec.Owner@).
Currently SeqAn defaults to a local position for @Class.StringSet@ classes (index_base.h).
..example.code:template < typename TString, typename TSpec >
struct SAValue< StringSet<TString, TSpec> > {
	typedef Pair<
		typename Size< StringSet<TString, TSpec> >::Type,
		typename SAValue<TString>::Type,
		Pack
	> Type;
};
..remarks.note:SAValue is the return type of various function, e.g. @Function.position@ for the @Class.Index@ @Class.Finder@ class, @Function.getOccurrence@, @Function.getOccurrences@ etc.
You should always use the type of this meta-function to store the return values.
If you want to write algorithms for both variants (local and global positions) you 
should use the functions @Function.posLocalize@, @Function.posGlobalize@, @Function.getSeqNo@ and @Function.getSeqOffset@.
..remarks.note:If $TObject$ is an @Class.Index@, @Metafunction.Position@ returns the same value as $SAValue$. You can change the position type of an index by overloading $SAValue$, not @Metafunction.Position@.
..include:seqan/index.h
*/

/*!
 * @mfn Index#SAValue
 * 
 * @headerfile seqan/index.h
 * 
 * @brief The default alphabet type of a suffix array, i.e. the type to store a
 *        position of a string or string set.
 * 
 * @signature SAValue<TObject>::Type
 * 
 * @tparam TObject A string, string set, or index type. Types: String,
 *                 StringSet, Index
 * 
 * @return TReturn A type to store a position.If <tt>TObject</tt> is a @link
 *                 String @endlink, it is a single integer value. By default
 *                 this is the @link Size @endlink type of <tt>TObject</tt>.If
 *                 <tt>TObject</tt> is a @link StringSet @endlink, it could be a
 *                 single integer too (called global position, see @link
 *                 ConcatDirect @endlink) or a @link Pair @endlink (called local
 *                 position, see @link Owner @endlink). Currently SeqAn defaults
 *                 to a local position for @link StringSet @endlink classes
 *                 (index_base.h).
 * 
 * @section Remarks
 * 
 * type=note:SAValue is the return type of various function, e.g. @link position
 * @endlink for the @link Index @endlink @link Finder @endlink class, @link
 * getOccurrence @endlink, @link getOccurrences @endlink etc. You should always
 * use the type of this meta-function to store the return values. If you want to
 * write algorithms for both variants (local and global positions) you should
 * use the functions @link posLocalize @endlink, @link posGlobalize @endlink,
 * @link getSeqNo @endlink and @link getSeqOffset @endlink.
 * 
 * type=note:If <tt>TObject</tt> is an @link Index @endlink, @link Position
 * @endlink returns the same value as <tt>SAValue</tt>. You can change the
 * position type of an index by overloading <tt>SAValue</tt>, not @link Position
 * @endlink.
 * 
 * @section Examples
 * 
 * @code{.cpp}
 * template < typename TString, typename TSpec >
 * struct SAValue< StringSet<TString, TSpec> > {
 * 	typedef Pair<
 * 		typename Size< StringSet<TString, TSpec> >::Type,
 * 		typename SAValue<TString>::Type,
 * 		Pack
 * 	> Type;
 * };
 * @endcode
 * @see orderOccurrences	
 */ 
    
    template <typename TObject>
	struct SAValue:
		Size<TObject> {};
	
	template <typename TObject>
	struct SAValue<TObject const>:
		SAValue<TObject> {};
	
	// to speed up sequence number computation
	// we use a pair of seqNo and localPosition
	template < typename TString, typename TSpec >
	struct SAValue< StringSet<TString, TSpec> > {
		typedef Pair<
			typename Size< StringSet<TString, TSpec> >::Type,
			typename SAValue<TString>::Type,
			Pack
		> Type;
	};

/*
	template < typename TString, typename TSpec >
	struct SAValue< StringSet<TString, TSpec> > {
		typedef Pair<
			typename Size< StringSet<TString, TSpec> >::Type,
			typename SAValue<TString>::Type,
			BitPacked<2,30>						    // max. 4 sequences 
		> Type;										// max. 2^30 characters each
	};
*/
	template < typename TText, typename TSpec >
	struct SAValue< Index<TText, TSpec> >:
		SAValue<TText> {};

    template < typename TObject, typename TSpec >
	struct DefaultIndexStringSpec< Index<TObject, TSpec> >:
		DefaultIndexStringSpec<TObject> {};

//////////////////////////////////////////////////////////////////////////////
// value and size type of an index

	template < typename TText, typename TSpec >
    struct Value< Index<TText, TSpec> > {
		typedef typename Value<
			typename Fibre< Index<TText, TSpec>, FibreRawText>::Type 
		>::Type Type;
    };

	template < typename TText, typename TSpec >
    struct Size< Index<TText, TSpec> > {
		typedef typename Size<
			typename Fibre< Index<TText, TSpec>, FibreRawText>::Type 
		>::Type Type;
    };

	template < typename TText, typename TSpec >
	struct Position< Index<TText, TSpec> >:
		SAValue< Index<TText, TSpec> > {};

//////////////////////////////////////////////////////////////////////////////
// infix of an index

	template < typename TText, typename TSpec >
    struct Infix< Index<TText, TSpec> >:
		public Infix<TText> {};

	template < typename TText, typename TSpec >
    struct Infix< Index<TText, TSpec> const >:
		public Infix<TText> {};

//////////////////////////////////////////////////////////////////////////////
// default table type

	template < typename TObject, typename TSpec, typename TFibre >
	struct Fibre< Index<TObject, TSpec>, Tag<TFibre> const > {
		typedef String< 
			typename Size< Index<TObject, TSpec> >::Type,
			typename DefaultIndexStringSpec< Index<TObject, TSpec> >::Type 
		> Type;
	};

//////////////////////////////////////////////////////////////////////////////
// original text

	template < typename TText, typename TSpec >
	struct Fibre< Index<TText, TSpec>, FibreText> {
		typedef TText Type;
	};

//////////////////////////////////////////////////////////////////////////////
// concatenated text

	template < typename TText, typename TSpec >
	struct Fibre< Index<TText, TSpec>, FibreRawText> {
		typedef typename Concatenator<TText>::Type Type;
	};

//////////////////////////////////////////////////////////////////////////////
// suffix array type

	template < typename TText, typename TSpec >
	struct Fibre< Index<TText, TSpec>, FibreSA> {
		typedef String<
			typename SAValue< Index<TText, TSpec> >::Type,
			typename DefaultIndexStringSpec< Index<TText, TSpec> >::Type 
		> Type;
	};

//////////////////////////////////////////////////////////////////////////////
// globalize functor

	template <typename InType, typename TLimitsString, typename Result = typename Value<TLimitsString>::Type>
	struct FunctorGlobalize : public ::unary_function<InType,Result>
	{
		TLimitsString const *limits;
		FunctorGlobalize() {}
		FunctorGlobalize(TLimitsString const &_limits) : limits(&_limits) {}

		inline Result operator()(InType const &x) const
		{
			return posGlobalize(x, *limits);
		}
    };

	template <typename InType, typename Result>
	struct FunctorGlobalize<InType, Nothing, Result> : public ::unary_function<InType,InType>
	{
		FunctorGlobalize() {}
		FunctorGlobalize(Nothing const &) {}

        inline InType operator()(InType const &x) const
        {
			return x;
		}
    };

//////////////////////////////////////////////////////////////////////////////
// raw suffix array contains integer offsets relative to raw text
/*
	template < typename TString, typename TSpec >
	struct Fibre< Index<TString, TSpec>, FibreRawSA>:
		public Fibre< Index<TString, TSpec> const, FibreSA> {};

	template < typename TString, typename TSSetSpec, typename TSpec >
	struct Fibre< Index<StringSet<TString, TSSetSpec>, TSpec>, FibreRawSA> 
	{
		typedef Index< StringSet<TString, TSSetSpec>, TSpec> TIndex;
		typedef ModifiedString<
			typename Fibre<TIndex, FibreSA>::Type,
			ModView< FunctorGlobalize< 
				typename Value< typename Fibre<TIndex, FibreSA>::Type >::Type,
				typename StringSetLimits<TString>::Type >
			>
		> Type;
	};
*/
	template < typename TText, typename TSpec >
	struct Fibre< Index<TText, TSpec>, FibreRawSA> 
	{
		typedef Index<TText, TSpec> TIndex;
		typedef ModifiedString<
			typename Fibre<TIndex, FibreSA>::Type,
			ModView< FunctorGlobalize< 
				typename Value< typename Fibre<TIndex, FibreSA>::Type >::Type,
				typename StringSetLimits<TText>::Type >
			>
		> Type;
	};

//////////////////////////////////////////////////////////////////////////////
// default burrows-wheeler table

	template < typename TText, typename TSpec >
	struct Fibre< Index<TText, TSpec>, FibreBwt> {
		typedef String <
			typename Value< Index<TText, TSpec> >::Type,
			typename DefaultIndexStringSpec< Index<TText, TSpec> >::Type
		> Type;
	};


//////////////////////////////////////////////////////////////////////////////
// default fibre creators

	template < typename TText, typename TSpec >
	struct DefaultIndexCreator<Index<TText, TSpec>, FibreSA> {
        typedef Skew7 Type;							// standard suffix array creator is skew7
    };

	template < typename TText, typename TSpec >
	struct DefaultIndexCreator<Index<TText, TSpec>, FibreLcp> {
        typedef Kasai Type;
    };

	template < typename TText, typename TSpec >
	struct DefaultIndexCreator<Index<TText, TSpec>, FibreBwt> {
        typedef Bwt Type;
    };

	template < typename TText, typename TSpec >
	struct DefaultIndexCreator<Index<TText, TSpec>, FibreChildtab> {
        typedef Childtab Type;
    };


	template <typename TText, typename TSpec>
	inline Holder<TText> & _dataHost(Index<TText, TSpec> &index) {
		return index.text;
	}
	template <typename TText, typename TSpec>
	inline Holder<TText> const & _dataHost(Index<TText, TSpec> const &index) {
		return index.text;
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.getFibre:
..summary:Returns a specific fibre of a container.
..signature:getFibre(index, fibreTag)
..class:Class.Index
..cat:Index
..param.index:The index holding the fibre.
...type:Class.Index
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.ESA Index Fibres
...type:Tag.QGram Index Fibres
...type:Tag.WOTD Index Fibres
...type:Tag.FM Index Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
..include:seqan/index.h
..example
...text:The following code shows a simple example how the function @Function.getFibre@ is used.
...file:demos/index/index_begin_range_goDown_representative_repLength.cpp
...output:The string ISSI occurs 2 times in MISSISSIPPI and has 4 characters.
The string ISSI occurs 2 times in MISSISSIPPI and has 4 characters.
*/

/*!
 * @fn Index#getFibre
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns a specific fibre of a container.
 * 
 * @signature getFibre(container, fibreTag)
 * 
 * @param fibreTag A tag that identifies the @link Index#Fibre @endlink. Types: @link IndexEsaFibres Index Esa Fibres @endlink, 
 * @link FMIndexFibres FM Index Fibres @endlink, @link IndexSaFibres Index SA Fibres @endlink, @link IndexWotdFibres Index Wotd Fibres @endlink, @link IndexDfiFibres Index Dfi Fibres @endlink and @link IndexQGramFibres Index QGram Fibres @endlink.
 *
 * @param container The container holding the fibre. Types: @link Index Index @endlink
 * 
 * @return TReturn A reference to the @link Fibre @endlink object.
 * 
 * @section Examples
 * 
 * @code{.cpp}
 * Index< String<char> > index_esa("tobeornottobe");
 *  
 * String<char> & text = getFibre(indexEsa, EsaText());
 * @endcode
 *
 * @see Index#Fibre
 */

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreText>::Type & 
	getFibre(Index<TText, TSpec> &index, FibreText) {
		return value(index.text);
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreText>::Type & 
	getFibre(Index<TText, TSpec> const &index, FibreText) {
		return value(index.text);
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreRawText>::Type & 
	getFibre(Index<TText, TSpec> &index, FibreRawText) {
		return concat(value(index.text));
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreRawText>::Type & 
	getFibre(Index<TText, TSpec> const &index, FibreRawText) {
		return concat(value(index.text));
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreSA>::Type & 
	getFibre(Index<TText, TSpec> &index, FibreSA) {
		return index.sa;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreSA>::Type & 
	getFibre(Index<TText, TSpec> const &index, FibreSA) {
		return index.sa;
	}

//////////////////////////////////////////////////////////////////////////////
/*
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreRawSA>::Type const & 
	getFibre(Index<TText, TSpec> &index, FibreRawSA) {
		return indexSA(index);
	}

	template <typename TString, typename TSSetSpec, typename TSpec>
	inline typename Fibre<Index<StringSet<TString, TSSetSpec>, TSpec>, FibreRawSA>::Type
	getFibre(Index<StringSet<TString, TSSetSpec>, TSpec> &index, FibreRawSA) 
	{
		typedef Index< StringSet<TString, TSSetSpec>, TSpec> TIndex;
		
		typedef FunctorGlobalize<
			typename Value< typename Fibre<TIndex, FibreSA>::Type >::Type,
			typename StringSetLimits<StringSet<TString, TSSetSpec> >::Type
		> TFunctor;
		
		typedef ModifiedString<
			typename Fibre<Index<StringSet<TString, TSSetSpec>, TSpec>, FibreSA>::Type,
			ModView< TFunctor >
		> ModString;

		return ModString(indexSA(index), TFunctor(stringSetLimits(indexText(index))));
	}
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreRawSA>::Type
	getFibre(Index<TText, TSpec> &index, FibreRawSA) 
	{
		typedef Index<TText, TSpec> TIndex;
		
		typedef FunctorGlobalize<
			typename Value< typename Fibre<TIndex, FibreSA>::Type >::Type,
			typename StringSetLimits<TText>::Type
		> TFunctor;
		
		typedef ModifiedString<
			typename Fibre<Index<TText, TSpec>, FibreSA>::Type,
			ModView< TFunctor >
		> ModString;

		return ModString(indexSA(index), TFunctor(stringSetLimits(indexText(index))));
	}
//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreLcp>::Type & 
	getFibre(Index<TText, TSpec> &index, FibreLcp) {
		return index.lcp;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreLcp>::Type & 
	getFibre(Index<TText, TSpec> const &index, FibreLcp) {
		return index.lcp;
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreLcpe>::Type & 
	getFibre(Index<TText, TSpec> &index, FibreLcpe) {
		return index.lcpe;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreLcpe>::Type & 
	getFibre(Index<TText, TSpec> const &index, FibreLcpe) {
		return index.lcpe;
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreChildtab>::Type & 
	getFibre(Index<TText, TSpec> &index, FibreChildtab) {
		return index.childtab;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreChildtab>::Type & 
	getFibre(Index<TText, TSpec> const &index, FibreChildtab) {
		return index.childtab;
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreBwt>::Type & 
	getFibre(Index<TText, TSpec> &index, FibreBwt) {
		return index.bwt;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreBwt>::Type & 
	getFibre(Index<TText, TSpec> const &index, FibreBwt) {
		return index.bwt;
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.Index#length
..cat:Index
..summary:The number of characters in the underlying text of the index is returned.
..signature:length(index)
..class:Class.Index
..param.index:The index to return the number of characters of.
...type:Class.Index
..returns:The number of characters in the raw underlying text of the index is returned.
...metafunction:Metafunction.Size
..remarks:If the underlying text is a @Class.StringSet@ then the sum of all characters of the sequneces in the string
set is returned.
..include:seqan/index.h
..example
...text:The following code shows how @Function.length@ can be used on an index in order to determine the number of characters in the underlying text.
...file:demos/index/index_length_countSequences.cpp
...output:Hit at position: < 1 , 2 >
Hit at position: < 0 , 0 >
 */

/*!
 * @fn Index#length
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the number of characters in the raw underlying text of the
 *        index.
 * 
 * @signature length(index)
 * 
 * @param index An index of a text. Types: @link Index @endlink
 * 
 * @return TSize Returns the number of characters in the raw underlying text of the
 *        index with TSize being the result of the @link Size @endlink metafunction
 *        of @link Index @endlink.
 */

	template <typename TText, typename TSpec>
	inline typename Size<Index<TText, TSpec> >::Type 
	length(Index<TText, TSpec> const &index) {
		return length(indexRawText(index));
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.countSequences
..cat:Index
..summary:Return the number of sequences in an index' underlying text.
..signature:countSequences(index)
..class:Class.Index
..param.index:The index to return the number of sequences of.
...type:Class.Index
..returns:The number of sequences in the index' underlying text.
...metafunction:Metafunction.Size
..include:seqan/index.h
..example
...text:The following code shows how @Function.countSequences@ can be used on an index in order to determine the number of sequences in the underlying text (which can be a @Class.StringSet@).
...file:demos/index/index_length_countSequences.cpp
...output:Hit at position: < 1 , 2 >
Hit at position: < 0 , 0 >

 */

/*!
 * @fn Index#countSequences
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Return the number of sequences in an index' underlying text.
 * 
 * @signature countSequences(index)
 * 
 * @param index The index to return the number of sequences of. Types: @link Index @endlink
 * 
 * @return TReturn The number of sequences in the index' underlying text.
 *                 Metafunctions: Metafunction.Size
 */

	template <typename TText, typename TSpec>
	inline typename Size<TText>::Type 
	countSequences(Index<TText, TSpec> const &index) {
		return countSequences(indexText(index));
	}

//////////////////////////////////////////////////////////////////////////////
// TODO(singer): Since this is a public function it should be documented
	template <typename TText, typename TSpec>
	struct GetSequenceByNo< Index<TText, TSpec> >
	{
		typedef typename GetSequenceByNo<TText>::Type Type;
	};

	template <typename TText, typename TSpec>
	struct GetSequenceByNo< Index<TText, TSpec> const>
	{
		typedef typename GetSequenceByNo<TText const>::Type Type;
	};

//////////////////////////////////////////////////////////////////////////////
// TODO(singer): Since this is a public function it should be documented
	template <typename TSeqNo, typename TText, typename TSpec>
	inline typename GetSequenceByNo< Index<TText, TSpec> >::Type
	getSequenceByNo(TSeqNo seqNo, Index<TText, TSpec> &index)
	{
		return getSequenceByNo(seqNo, indexText(index));
	}

	template <typename TSeqNo, typename TText, typename TSpec>
	inline typename GetSequenceByNo< Index<TText, TSpec> const>::Type
	getSequenceByNo(TSeqNo seqNo, Index<TText, TSpec> const &index)
	{
		return getSequenceByNo(seqNo, indexText(index));
	}

//////////////////////////////////////////////////////////////////////////////
// TODO(singer): Since this is a public function it should be documented
	template <typename TSeqNo, typename TText, typename TSpec>
	inline typename Size<Index<TText, TSpec> >::Type 
	sequenceLength(TSeqNo seqNo, Index<TText, TSpec> const &index) {
		return sequenceLength(seqNo, indexText(index));
	}

//////////////////////////////////////////////////////////////////////////////
// TODO(singer): Since this is a public function it should be documented
	template <typename TPos, typename TText, typename TSpec>
	inline typename Size<Index<TText, TSpec> >::Type 
	suffixLength(TPos pos, Index<TText, TSpec> const &index)
    {
		return length(indexText(index)) - pos;
	}

	template <typename TPos, typename TString, typename TSSetSpec, typename TSpec>
	inline typename Size<Index<StringSet<TString, TSSetSpec>, TSpec> >::Type 
	suffixLength(TPos pos, Index<StringSet<TString, TSSetSpec>, TSpec> const &index)
    {
        typename StringSetLimits<StringSet<TString, TSSetSpec> >::Type const &limits = stringSetLimits(index);
		return sequenceLength(getSeqNo(pos, limits), index) - getSeqOffset(pos, limits);
	}


//////////////////////////////////////////////////////////////////////////////
/**
.Function.textAt:
..summary:Shortcut for $value(indexText(..), ..)$.
..cat:Index
..signature:textAt(position, index)
..class:Class.Index
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference or proxy to the value.
..include:seqan/index.h
..example
...text:The following code shows how the BWT of a text can be computed.
...file:demos/index/index_textAt_indexText_saAt_indexRequire.cpp
...output:BWT	Suffices
P	PI
S	SIPPI
S	SISSIPPI
M	MISSISSIPPI
I	I
P	PPI
I	IPPI
S	SSIPPI
S	SSISSIPPI
I	ISSIPPI
I	ISSISSIPPI*/

/*!
 * @fn Index#textAt
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>value(indexText(..), ..)</tt>.
 * 
 * @signature textAt(position, index)
 * 
 * @param index The @link Index @endlink object. Types: @link Index @endlink
 *
 * @param position A position in the array on which the value should be
 *                 accessed.
 * 
 * @return TReturn A reference or proxy to the value.
 *
 * @section Note The result of this function when used on an Index<TText, FMIndex<TOccSpec, Compress> > is not defined.
 */

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, FibreRawText>::Type>::Type 
	textAt(TPos i, TIndex &index) {
		return value(getFibre(index, FibreRawText()), i);
	}
	template <typename TPos, typename TString, typename TSSetSpec, typename TSpec>
	inline typename Reference<typename Fibre< Index< StringSet<TString, TSSetSpec>, TSpec>, FibreRawText>::Type>::Type 
	textAt(TPos i, Index< StringSet<TString, TSSetSpec>, TSpec> &index) {
		return value(getFibre(index, FibreRawText()), posGlobalize(i, stringSetLimits(index)));
	}
	template <typename TPos, typename TString, typename TSSetSpec, typename TSpec>
	inline typename Reference<typename Fibre< Index< StringSet<TString, TSSetSpec>, TSpec> const, FibreRawText>::Type>::Type 
	textAt(TPos i, Index< StringSet<TString, TSSetSpec>, TSpec> const &index) {
		return value(getFibre(index, FibreRawText()), posGlobalize(i, stringSetLimits(index)));
	}
	template <typename TPos, typename TString, typename TSpec>
	inline typename Reference<typename Fibre< Index< StringSet<TString, Owner<Default> >, TSpec>, FibreRawText>::Type>::Type 
	textAt(TPos i, Index< StringSet<TString, Owner<Default> >, TSpec> &index) {
		Pair <
			typename Size< StringSet<TString, Owner<Default> > >::Type,
			typename Size< TString >::Type > locPos;
		posLocalize(locPos, i, stringSetLimits(index));
		return value(value(getFibre(index, FibreText()), getValueI1(locPos)), getValueI2(locPos));
	}
	template <typename TPos, typename TString, typename TSpec>
	inline typename Reference<typename Fibre< Index< StringSet<TString, Owner<Default> >, TSpec> const, FibreRawText>::Type>::Type 
	textAt(TPos i, Index< StringSet<TString, Owner<Default> >, TSpec> const &index) {
		Pair <
        typename Size< StringSet<TString, Owner<Default> > >::Type,
        typename Size< TString >::Type > locPos;
		posLocalize(locPos, i, stringSetLimits(index));
		return value(value(getFibre(index, FibreText()), getValueI1(locPos)), getValueI2(locPos));
	}

//////////////////////////////////////////////////////////////////////////////
// infix

	template <typename TText, typename TSpec, typename TPosBegin, typename TPosEnd>
	inline typename Infix<TText>::Type
	infix(Index<TText, TSpec> &index, TPosBegin pos_begin, TPosEnd pos_end)
	{
		return infix(indexText(index), pos_begin, pos_end);
	}

	template <typename TText, typename TSpec, typename TPosBegin, typename TPosEnd>
	inline typename Infix<TText>::Type
	infix(Index<TText, TSpec> const &index, TPosBegin pos_begin, TPosEnd pos_end)
	{
		return infix(indexText(index), pos_begin, pos_end);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.rawtextAt:
..summary:Shortcut for $value(indexRawText(..), ..)$.
..cat:Index
..signature:rawtextAt(position, index)
..class:Class.Index
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference or proxy to the value.
..include:seqan/index.h
*/
/*!
 * @fn Index#rawtextAt
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>value(indexRawText(..), ..)</tt>.
 * 
 * @signature rawtextAt(position, index)
 * 
 * @param index The @link Index @endlink object. Types: @link Index @endlink
 *
 * @param position A position in the array on which the value should be
 *                 accessed.
 * 
 * @return TReturn A reference or proxy to the value.
 *
 * @section Note The result of this function when used on an Index<TText, FMIndex<TOccSpec, Compress> > is not defined.
 */

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, FibreRawText>::Type>::Type rawtextAt(TPos i, TIndex &index) {
		return value(getFibre(index, FibreRawText()), i);
	}
	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex const, FibreRawText>::Type>::Type rawtextAt(TPos i, TIndex const &index) {
		return value(getFibre(index, FibreRawText()), i);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.saAt:
..summary:Shortcut for $value(indexSA(..), ..)$.
..cat:Index
..signature:saAt(position, index)
..class:Class.Index
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference or proxy to the value.
..include:seqan/index.h
..example
...text:The following code shows how the BWT of a text can be computed.
...file:demos/index/index_textAt_indexText_saAt_indexRequire.cpp
...output:BWT	Suffices
P	PI
S	SIPPI
S	SISSIPPI
M	MISSISSIPPI
I	I
P	PPI
I	IPPI
S	SSIPPI
S	SSISSIPPI
I	ISSIPPI
I	ISSISSIPPI*/
/*!
 * @fn IndexEsa#saAt
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>value(indexSA(..), ..)</tt>.
 *
 * @deprecated advanced
 * 
 * @signature saAt(position, index)
 * 
 * @param index The @link Index @endlink object holding the fibre. Types:
 *              @link IndexEsa @endlink
 * @param position A position in the array on which the value should be
 *                 accessed.
 * 
 * @return TReturn A reference or proxy to the value.
 */

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, FibreSA>::Type>::Type saAt(TPos i, TIndex &index) {
		return value(getFibre(index, FibreSA()), i);
	}
	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex const, FibreSA>::Type>::Type saAt(TPos i, TIndex const &index) {
		return value(getFibre(index, FibreSA()), i);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.rawsaAt:
..summary:Shortcut for $value(indexRawSA(..), ..)$.
..cat:Index
..signature:saAt(position, index)
..class:Class.Index
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference or proxy to the value.
..include:seqan/index.h
*/

/*!
 * @fn IndexEsa#rawsaAt
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>value(indexRawSA(..), ..)</tt>.
 * 
 * @deprecated. advanced
 *
 * @signature rawsaAt(position, index)
 * 
 * @param index The @link Index @endlink object holding the fibre. Types:
 *              @link IndexEsa @endlink
 * @param position A position in the array on which the value should be
 *                 accessed.
 */

	template <typename TPos, typename TIndex>
	inline typename Value<typename Fibre<TIndex const, FibreRawSA>::Type>::Type rawsaAt(TPos i, TIndex const &index) {
		return posGlobalize(saAt(i, index), stringSetLimits(indexText(index)));
	}


//////////////////////////////////////////////////////////////////////////////
/**
.Function.IndexEsa#lcpAt:
..summary:Shortcut for $value(indexLcp(..), ..)$.
..cat:Index
..signature:lcpAt(position, index)
..class:Spec.IndexEsa
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Spec.IndexEsa@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference or proxy to the value.
..include:seqan/index.h
*/
/*!
 * @fn IndexEsa#lcpAt
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>value(indexLcp(..), ..)</tt>.
 * 
 * @signature lcpAt(position, index)
 * 
 * @param index The @link Index @endlink object holding the fibre. Types:
 *              @link IndexEsa @endlink
 * @param position A position in the array on which the value should be
 *                 accessed.
 * 
 * @return TReturn A reference or proxy to the value.
 */

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, FibreLcp>::Type>::Type lcpAt(TPos i, TIndex &index) {
		return value(getFibre(index, FibreLcp()), i);
	}
	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex const, FibreLcp>::Type>::Type lcpAt(TPos i, TIndex const &index) {
		return value(getFibre(index, FibreLcp()), i);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.IndexEsa#lcpeAt:
..summary:Shortcut for $value(indexLcpe(..), ..)$.
..cat:Index
..signature:lcpeAt(position, index)
..class:Spec.IndexEsa
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Spec.IndexEsa@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference or proxy to the value.
..include:seqan/index.h
*/
/*!
 * @fn IndexEsa#lcpeAt
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>value(indexLcpe(..), ..)</tt>.
 * 
 * @signature lcpeAt(position, index)
 * 
 * @param index The @link Index @endlink object holding the fibre. Types:
 *              @link IndexEsa @endlink
 * @param position A position in the array on which the value should be
 *                 accessed.
 * 
 * @return TReturn A reference or proxy to the value.
 */

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, FibreLcpe>::Type>::Type lcpeAt(TPos i, TIndex &index) {
		return value(getFibre(index, FibreLcpe()), i);
	}
	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex const, FibreLcpe>::Type>::Type lcpeAt(TPos i, TIndex const &index) {
		return value(getFibre(index, FibreLcpe()), i);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.IndexEsa#childAt:
..summary:Shortcut for $value(indexChildtab(..), ..)$.
..cat:Index
..signature:childAt(position, index)
..class:Spec.IndexEsa
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Spec.IndexEsa@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference or proxy to the value.
..include:seqan/index.h
*/

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, FibreChildtab>::Type>::Type childAt(TPos i, TIndex &index) {
		return value(getFibre(index, FibreChildtab()), i);
	}
	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex const, FibreChildtab>::Type>::Type childAt(TPos i, TIndex const &index) {
		return value(getFibre(index, FibreChildtab()), i);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.IndexEsa#bwtAt:
..summary:Shortcut for $value(indexBwt(..), ..)$.
..cat:Index
..signature:bwtAt(position, index)
..class:Spec.IndexEsa
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Spec.IndexEsa@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference or proxy to the value.
..include:seqan/index.h
*/
/*!
 * @fn IndexEsa#childAt
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>value(indexChildtab(..), ..)</tt>.
 * 
 * @signature childAt(position, index)
 * 
 * @param index The @link Index @endlink object holding the fibre. Types:
 *              @link IndexEsa @endlink
 * @param position A position in the array on which the value should be
 *                 accessed.
 * 
 * @return TReturn A reference or proxy to the value.
 */

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, FibreBwt>::Type>::Type bwtAt(TPos i, TIndex &index) {
		return value(getFibre(index, FibreBwt()), i);
	}
	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex const, FibreBwt>::Type>::Type bwtAt(TPos i, TIndex const &index) {
		return value(getFibre(index, FibreBwt()), i);
	}

//////////////////////////////////////////////////////////////////////////////

    template <typename TIndex, typename TPos, typename TSize>
	inline typename SAValue<TIndex>::Type toSuffixPosition(TIndex &, TPos i, TSize) {
        return i;
	}
    template <typename TIndex, typename TPos, typename TSize>
	inline typename SAValue<TIndex const>::Type toSuffixPosition(TIndex const &, TPos i, TSize) {
        return i;
	}

//////////////////////////////////////////////////////////////////////////////
// interface for infinity/invalid values

	template <typename TValue>
	inline void _setSizeInval(TValue &v) {
		v = MaxValue<TValue>::VALUE;
	}

	template <typename TValue>
	inline bool _isSizeInval(TValue const &v) {
//IOREV _notio_
		return v == MaxValue<TValue>::VALUE;
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexText:
..summary:Shortcut for $getFibre(index, FibreText())$.
..cat:Index
..signature:indexText(index)
..class:Class.Index
..param.index:The @Class.Index@ object holding the fibre.
...type:Class.Index
..returns:A reference to the text fibre (original text).
..include:seqan/index.h
..example
...text:The following code shows how the BWT of a text can be computed.
...file:demos/index/index_textAt_indexText_saAt_indexRequire.cpp
...output:BWT	Suffices
P	PI
S	SIPPI
S	SISSIPPI
M	MISSISSIPPI
I	I
P	PPI
I	IPPI
S	SSIPPI
S	SSISSIPPI
I	ISSIPPI
I	ISSISSIPPI*/
/*!
 * @fn Index#indexText
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>getFibre(.., EsaText)</tt>.
 * 
 * @signature indexText(index)
 * 
 * @param index The @link Index @endlink object holding the fibre. Types:
 *              @link Index @endlink
 * 
 * @return TReturn A reference to the text of the index.
 *
 * @section Note The result of this function when used on an Index<TText, FMIndex<TOccSpec, Compress> > is not defined.
 *
 */

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreText>::Type & indexText(Index<TText, TSpec> &index) { return getFibre(index, FibreText()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreText>::Type & indexText(Index<TText, TSpec> const &index) { return getFibre(index, FibreText()); }

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename StringSetLimits<TText const>::Type
	stringSetLimits(Index<TText, TSpec> &) { 
		return Nothing(); 
	}

	template <typename TText, typename TSpec>
	inline typename StringSetLimits<TText const>::Type
	stringSetLimits(Index<TText, TSpec> const &) { 
		return Nothing(); 
	}

	template <typename TString, typename TSSetSpec, typename TSpec>
	inline typename StringSetLimits< StringSet<TString, TSSetSpec> const >::Type & 
	stringSetLimits(Index<StringSet<TString, TSSetSpec>, TSpec> &index) {
		return stringSetLimits(indexText(index)); 
	}

	template <typename TString, typename TSSetSpec, typename TSpec>
	inline typename StringSetLimits< StringSet<TString, TSSetSpec> const >::Type & 
	stringSetLimits(Index<StringSet<TString, TSSetSpec>, TSpec> const &index) {
		return stringSetLimits(indexText(index)); 
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexRawText:
..summary:Shortcut for $getFibre(.., EsaRawText)$.
..cat:Index
..signature:indexRawText(index)
..class:Class.Index
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference to the @Tag.ESA Index Fibres.EsaRawText@ fibre (concatenated input text).
..include:seqan/index.h
*/
//TODO(singer) The RawText Fibre exist for more then the Esa index
/*!
 * @fn IndexEsa#indexRawText
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>$getFibre(.., EsaRawText)</tt>.
 * 
 * @signature rawtextAt(position, index)
 * 
 * @param index The @link Index @endlink object. Types: @link Index @endlink
 *
 * @param position A position in the array on which the value should be
 *                 accessed.
 * 
 * @return TReturn A reference or proxy to the value.
 */

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreRawText>::Type & indexRawText(Index<TText, TSpec> &index) { return getFibre(index, FibreRawText()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreRawText>::Type & indexRawText(Index<TText, TSpec> const &index) { return getFibre(index, FibreRawText()); }

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexSA:
..summary:Shortcut for $getFibre(.., FibreSA)$.
..cat:Index
..signature:indexSA(index)
..class:Class.Index
..param.index:The @Class.Index@ object holding the fibre.
...type:Class.Index
..returns:A reference to the suffix array fibre.
..include:seqan/index.h
..example
...text:The following code shows how the BWT of a text can be computed.
...file:demos/index/index_textAt_indexText_saAt_indexRequire.cpp
...output:BWT	Suffices
P	PI
S	SIPPI
S	SISSIPPI
M	MISSISSIPPI
I	I
P	PPI
I	IPPI
S	SSIPPI
S	SSISSIPPI
I	ISSIPPI
I	ISSISSIPPI*/
//TODO(singer) The function in not only defined for the esa index
/*!
 * @fn IndexEsa#indexSA
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>getFibre(.., EsaSA)</tt>.
 * 
 * @signature indexSA(index)
 * 
 * @param index The @link Index @endlink object holding the fibre. Types:
 *              @link IndexEsa @endlink
 * 
 * @return TReturn A reference to the @link ESA Index Fibres.EsaSA @endlink
 *                 fibre (suffix array).
 */

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreSA>::Type & indexSA(Index<TText, TSpec> &index) { return getFibre(index, FibreSA()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreSA>::Type & indexSA(Index<TText, TSpec> const &index) { return getFibre(index, FibreSA()); }

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexRawSA:
..summary:Shortcut for $getFibre(.., EsaRawSA)$.
..cat:Index
..signature:indexRawSA(index)
..class:Class.Index
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference to the @Tag.ESA Index Fibres.EsaRawSA@ fibre (suffix array).
..include:seqan/index.h
*/
/*!
 * @fn Index#indexRawSA
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>getFibre(.., EsaRawSA)</tt>.
 * 
 * @signature indexRawSA(index)
 * 
 * @param index The @link Index @endlink object holding the fibre. Types:
 *              @link IndexEsa @endlink
 * 
 * @return TReturn A reference to the @link ESA Index Fibres.EsaRawSA @endlink
 *                 fibre (suffix array).
 */

/*
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreRawSA>::Type const & indexRawSA(Index<TText, TSpec> &index) { return getFibre(index, FibreRawSA()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreRawSA>::Type const & indexRawSA(Index<TText, TSpec> const &index) { return getFibre(index, FibreRawSA()); }
*/
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreRawSA>::Type indexRawSA(Index<TText, TSpec> &index) { return getFibre(index, FibreRawSA()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreRawSA>::Type indexRawSA(Index<TText, TSpec> const &index) { return getFibre(index, FibreRawSA()); }

//////////////////////////////////////////////////////////////////////////////
/**
.Function.IndexEsa#indexLcp:
..summary:Shortcut for $getFibre(.., EsaLcp)$.
..cat:Index
..signature:indexLcp(index)
..class:Spec.IndexEsa
..param.index:The @Spec.IndexEsa@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference to the @Tag.ESA Index Fibres.EsaLcp@ fibre (lcp table).
..include:seqan/index.h
*/
/*!
 * @fn IndexEsa#indexLcp
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>getFibre(.., EsaLcp)</tt>.
 *
 * deprecated. advanced
 * 
 * @signature indexLcp(index)
 * 
 * @param index The @link Index @endlink object holding the fibre. Types:
 *              @link IndexEsa @endlink
 * 
 * @return TReturn A reference to the @link ESA Index Fibres.EsaLcp @endlink
 *                 fibre (lcp table).
 */

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreLcp>::Type & indexLcp(Index<TText, TSpec> &index) { return getFibre(index, FibreLcp()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreLcp>::Type & indexLcp(Index<TText, TSpec> const &index) { return getFibre(index, FibreLcp()); }

//////////////////////////////////////////////////////////////////////////////
/**
.Function.IndexEsa#indexLcpe:
..summary:Shortcut for $getFibre(.., EsaLcpe)$.
..cat:Index
..signature:indexLcpe(index)
..class:Spec.IndexEsa
..param.index:The @Spec.IndexEsa@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference to the @Tag.ESA Index Fibres.EsaLcpe@ fibre (enhanced lcp table).
..include:seqan/index.h
*/
/*!
 * @fn IndexEsa#indexLcpe
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>getFibre(.., EsaLcpe)</tt>.
 * 
 * @signature indexLcpe(index)
 * 
 * @param index The @link Index @endlink object holding the fibre. Types:
 *              @link IndexEsa @endlink
 * 
 * @return TReturn A reference to the @link ESA Index Fibres.EsaLcpe @endlink
 *                 fibre (enhanced lcp table).
 */

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreLcpe>::Type & indexLcpe(Index<TText, TSpec> &index) { return getFibre(index, FibreLcpe()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreLcpe>::Type & indexLcpe(Index<TText, TSpec> const &index) { return getFibre(index, FibreLcpe()); }

//////////////////////////////////////////////////////////////////////////////
/**
.Function.IndexEsa#indexBwt:
..summary:Shortcut for $getFibre(.., EsaBwt)$.
..cat:Index
..signature:indexBwt(index)
..class:Spec.IndexEsa
..param.index:The @Spec.IndexEsa@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference to the @Tag.ESA Index Fibres.EsaBwt@ fibre (Burrows-Wheeler table).
..include:seqan/index.h
*/
/*
 * @fn IndexEsa#indexBwt
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>getFibre(.., EsaBwt)</tt>.
 * 
 * @signature indexBwt(index)
 * 
 * @param index The @link IndexEsa @endlink object holding the fibre. Types:
 *              @link IndexEsa @endlink
 * 
 * @return TReturn A reference to the @link ESA Index Fibres.EsaBwt @endlink
 *                 fibre (Burrows-Wheeler table).
 */

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreBwt>::Type & indexBwt(Index<TText, TSpec> &index) { return getFibre(index, FibreBwt()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreBwt>::Type & indexBwt(Index<TText, TSpec> const &index) { return getFibre(index, FibreBwt()); }

//////////////////////////////////////////////////////////////////////////////
/**
.Function.IndexEsa#indexChildtab:
..summary:Shortcut for $getFibre(.., EsaChildtab)$.
..cat:Index
..signature:indexChildtab(index)
..class:Spec.IndexEsa
..param.index:The @Spec.IndexEsa@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference to the @Tag.ESA Index Fibres.EsaChildtab@ fibre (child table).
..include:seqan/index.h
*/
/*!
 * @fn IndexEsa#indexChildtab
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Shortcut for <tt>getFibre(.., EsaChildtab)</tt>.
 * 
 * @signature indexChildtab(index)
 * 
 * @param index The @link Index @endlink object holding the fibre. Types:
 *              @link IndexEsa @endlink
 * 
 * @return TReturn A reference to the @link ESA Index Fibres.EsaChildtab
 *                 @endlink fibre (child table).
 */

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreChildtab>::Type & indexChildtab(Index<TText, TSpec> &index) { return getFibre(index, FibreChildtab()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreChildtab>::Type & indexChildtab(Index<TText, TSpec> const &index) { return getFibre(index, FibreChildtab()); }


// ----------------------------------------------------------------------------
// Function open
// ----------------------------------------------------------------------------
/**
.Function.Index#open
..class:Class.Index
..summary:This functions opens an index from disk.
..signature:open(index, fileName [, mode])
..param.index:The index to be opened.
...type:Class.Index
..param.fileName:C-style character string containing the file name.
..param.mode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
..example
...text:The following code shows how the function @Function.open@ is used with indices.
...file:demos/index/index_open_save.cpp
...output:1
1
*/
/*!
 * @fn Index#open
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions opens an index from disk.
 * 
 * @signature open(index, fileName [, mode])
 * 
 * @param mode The combination of flags defining how the file should be
 *                 opened.To open a file read-only, write-only or to read and
 *                 write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                 <tt>OPEN_RDWR</tt>.To create or overwrite a file add
 *                 <tt>OPEN_CREATE</tt>.To append a file if existing add
 *                 <tt>OPEN_APPEND</tt>.To circumvent problems, files are always
 *                 opened in binary mode. Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                 OPEN_APPEND</tt>
 *
 * @param index The index to be opened. Types: @link Index @endlink
 *
 * @param fileName C-style character string containing the file name.
 * 
 * @return TReturn A <tt>bool</tt> which is <tt>true</tt> on success.
 */

// ----------------------------------------------------------------------------
// Function save
// ----------------------------------------------------------------------------

/*! 
 * @fn Index#save
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions saves an index to disk.
 * 
 * @signature save(index, fileName [, mode])
 * 
 * @param mode The combination of flags defining how the file should be
 *                 opened.To open a file read-only, write-only or to read and
 *                 write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                 <tt>OPEN_RDWR</tt>.To create or overwrite a file add
 *                 <tt>OPEN_CREATE</tt>.To append a file if existing add
 *                 <tt>OPEN_APPEND</tt>.To circumvent problems, files are always
 *                 opened in binary mode. Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                 OPEN_APPEND</tt>
 *
 * @param index The index to be saved to disk. Types: @link Index @endlink
 *
 * @param fileName C-style character string containing the file name.
 * 
 * @return TReturn A <tt>bool</tt> which is <tt>true</tt> on success.
 */
/**
.Function.Index#save
..class:Class.Index
..summary:This functions saves an index to disk.
..signature:save(index, fileName [, mode])
..param.index:The index to be saved to disk.
...type:Class.Index
..param.fileName:C-style character string containing the file name.
..param.mode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
..example
...text:The following code shows how the function @Function.open@ is used with indices.
...file:demos/index/index_open_save.cpp
...output:1
1
*/

}

#endif

