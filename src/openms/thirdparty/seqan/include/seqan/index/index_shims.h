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

#ifndef SEQAN_HEADER_INDEX_SHIMS_H
#define SEQAN_HEADER_INDEX_SHIMS_H

namespace SEQAN_NAMESPACE_MAIN
{

	//////////////////////////////////////////////////////////////////////////////
	// Suffix Array creation wrappers
	//////////////////////////////////////////////////////////////////////////////

	// build suffix array with an external pipeling algorithm (skew3, skew7, ...)
	template < 
		typename TSA, 
		typename TObject, 
		typename TAlgSpec >
	void _createSuffixArrayPipelining(
		TSA &suffixArray,
		TObject const &text,
		TAlgSpec const)
	{
	SEQAN_CHECKPOINT
        // signed characters behave different than unsigned when compared
        // to get the same index with signed or unsigned chars we simply cast them to unsigned
        // before feeding them into the pipeline
        typedef typename MakeUnsigned_< typename Value<TObject>::Type >::Type TUValue;

        // specialization
		typedef Pipe< TObject, Source<> >				src_t;
        typedef Pipe< src_t, Caster<TUValue> >          unsigner_t;
		typedef Pipe< unsigner_t, TAlgSpec >	        creator_t;

		// instantiation and processing
		src_t		src(text);
        unsigner_t  unsigner(src);
		creator_t	creator(unsigner);

		suffixArray << creator;
		#ifdef SEQAN_TEST_INDEX
			isSuffixArray(suffixArray, text);
		#endif
	}


	// build suffix array (external) for mutliple sequences
	template < 
		typename TSA, 
		typename TString, 
		typename TSpec,
		typename TAlgSpec >
	void _createSuffixArrayPipelining(
		TSA &suffixArray,
		StringSet<TString, TSpec> const &stringSet,
		TAlgSpec const)
	{
	SEQAN_CHECKPOINT
        // signed characters behave different than unsigned when compared
        // to get the same index with signed or unsigned chars we simply cast them to unsigned
        // before feeding them into the pipeline
		typedef typename Concatenator<StringSet<TString, TSpec> >::Type			TConcat;
        typedef typename MakeUnsigned_< typename Value<TConcat>::Type >::Type	TUValue;
		typedef Multi<
			TAlgSpec, 
			typename Value<TSA>::Type, 
			typename StringSetLimits<StringSet<TString, TSpec> >::Type >        MultiConstrSpec;

        // specialization
		typedef Pipe< TConcat, Source<> >				src_t;
        typedef Pipe< src_t, Caster<TUValue> >          unsigner_t;
		typedef Pipe< unsigner_t, MultiConstrSpec >	    creator_t;

		// instantiation and processing
		src_t		src(concat(stringSet));
        unsigner_t  unsigner(src);
		creator_t	creator(unsigner, stringSetLimits(stringSet));

		suffixArray << creator;
		#ifdef SEQAN_TEST_INDEX
			//isSuffixArray(suffixArray, stringSet);
		#endif
	}


/**
.Function.createSuffixArray:
..summary:Creates a suffix array from a given text.
..cat:Index
..signature:createSuffixArray(suffixArray, text[, algo_tag])
..param.suffixArray:The resulting suffix array.
..param.text:A given text.
..param.algo_tag:A tag that identifies the algorithm which is used for creation.
..remarks:This function should not be called directly. Please use @Function.indexCreate@ or @Function.indexRequire@.
The size of $suffixArray$ must be at least $length(text)$ before calling this function.
..include:seqan/index.h
*/
/*!
 * @fn createSuffixArray
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Creates a suffix array from a given text.
 * 
 * @signature createSuffixArray(suffixArray, text[, algo_tag])
 * 
 * @param text A given text. Types: @link SequenceConcept @endlink
 * @param algo_tag A tag that identifies the algorithm which is used for
 *                 creation.
 * @param suffixArray The resulting suffix array.
 * 
 * @section Remarks
 * 
 * This function should not be called directly. Please use @link indexCreate
 * @endlink or @link indexRequire @endlink. The size of <tt>suffixArray</tt>
 * must be at least <tt>length(text)</tt> before calling this function.
 * 
 * Demo: Demo.Suffix Array
 */
    template < typename TSA,
               typename TText,
			   typename TAlgSpec >
    inline void _createSuffixArrayRandomAccess(
		TSA &sa,
		TText const &s,
		TAlgSpec const &alg)
	{
	SEQAN_CHECKPOINT
		// -> call internal memory algorithm with an extended interface (+ alphabet size, max_depth)
		if (BitsPerValue< typename Value<TText>::Type >::VALUE > 16)
			createSuffixArray(sa, s, alg, length(s), 0);
		else
			createSuffixArray(sa, s, alg, ValueSize< typename Value<TText>::Type >::VALUE, 0);
	}

	template < 
		typename TSA, 
		typename TText, 
		typename TAlgSpec >
	inline void _createSuffixArrayWrapper(
		TSA &sa,
		TText const &s,
		TAlgSpec const &alg,
        True)
	{
	SEQAN_CHECKPOINT
        _createSuffixArrayRandomAccess(sa, s, alg);
	}

    // always use external Skew7 for multiple strings
	template < 
		typename TSA, 
		typename TSequence, 
		typename TSetSpec,
		typename TAlgSpec >
	inline void _createSuffixArrayWrapper(
		TSA &sa,
		StringSet<TSequence, TSetSpec> const &s,
		TAlgSpec const &,
        True)
	{
	SEQAN_CHECKPOINT
        _createSuffixArrayPipelining(sa, s, Skew7());
	}

	template < 
		typename TSA, 
		typename TText,
		typename TAlgSpec >
	inline void _createSuffixArrayWrapper(
		TSA &sa,
		TText const &s,
		TAlgSpec const &alg,
        False)
	{
	SEQAN_CHECKPOINT
        _createSuffixArrayPipelining(sa, s, alg);
	}

	template < 
		typename TSA,
		typename TText,
        typename TAlgSpec >
	inline void createSuffixArray(
		TSA &sa,
		TText const &s,
		TAlgSpec const &alg)
	{
	SEQAN_CHECKPOINT
        _createSuffixArrayWrapper(sa, s, alg, typename SACreatorRandomAccess_<TSA, TText, TAlgSpec>::Type());
    }



//____________________________________________________________________________



	//////////////////////////////////////////////////////////////////////////////
	// LCP Table creation wrappers
	//////////////////////////////////////////////////////////////////////////////

	template < 
        typename TLCPTable,
		typename TObject, 
        typename TSA,
		typename TAlgSpec >
	void _createLCPTablePipelining(
		TLCPTable &LCP,
		TObject const &text,
		TSA const &suffixArray,
		TAlgSpec const)
	{
	SEQAN_CHECKPOINT
		// specialization
		typedef Pipe< TObject, Source<> >							srcText_t;
		typedef Pipe< TSA, Source<> >   							srcSA_t;
	    typedef Pipe< Bundle2< srcText_t, srcSA_t >, TAlgSpec >		creator_t;

		// instantiation and processing
		srcText_t	srcText(text);
		srcSA_t		srcSA(suffixArray);
		creator_t	creator(bundle2(srcText, srcSA));

		LCP << creator;
		#ifdef SEQAN_TEST_INDEX
			isLCPTable(LCP, suffixArray, text);
		#endif
	}


	// build lcp table (external) for mutliple sequences
	template < 
        typename TLCPTable,
		typename TString,
		typename TSpec,
        typename TSA,
		typename TAlgSpec >
	void _createLCPTablePipelining(
		TLCPTable &LCP,
		StringSet<TString, TSpec> const &stringSet,
		TSA const &suffixArray,
		TAlgSpec const)
	{
	SEQAN_CHECKPOINT
		typedef typename Concatenator<StringSet<TString, TSpec> >::Type TConcat;
		typedef Multi<
			TAlgSpec, 
			typename Value<TSA>::Type, 
			typename StringSetLimits<StringSet<TString, TSpec> >::Type > MultiConstrSpec;

		// specialization
		typedef Pipe< TConcat, Source<> >								srcText_t;
		typedef Pipe< TSA, Source<> >   								srcSA_t;
	    typedef Pipe< Bundle2< srcText_t, srcSA_t >, MultiConstrSpec >	creator_t;

		// instantiation and processing
		srcText_t	srcText(concat(stringSet));
		srcSA_t		srcSA(suffixArray);
		creator_t	creator(bundle2(srcText, srcSA), stringSetLimits(stringSet));

		LCP << creator;
		#ifdef SEQAN_TEST_INDEX
			isLCPTable(LCP, suffixArray, text);
		#endif
	}


/**
.Function.createLcpTable:
..summary:Creates a lcp table from a given text and suffix array.
..cat:Index
..signature:createLcpTable(lcp, text, suffixArray[, algo_tag])
..param.lcp:The resulting lcp table.
..param.text:A given text.
..param.suffixArray:The suffix array of $text$.
..param.algo_tag:A tag that identifies the algorithm which is used for creation.
..remarks:This function should not be called directly. Please use @Function.indexCreate@ or @Function.indexRequire@.
The size of $lcp$ must be at least $length(text)$ before calling this function.
..include:seqan/index.h
*/
/*!
 * @fn createLcpTable
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Creates a lcp table from a given text and suffix array.
 * 
 * @signature createLcpTable(lcp, text, suffixArray[, algo_tag])
 * 
 * @param text A given text. Types: @link SequenceConcept @endlink
 * @param algo_tag A tag that identifies the algorithm which is used for
 *                 creation.
 * @param suffixArray The suffix array of <tt>text</tt>.
 * @param lcp The resulting lcp table.
 * 
 * @section Remarks
 * 
 * This function should not be called directly. Please use @link indexCreate
 * @endlink or @link indexRequire @endlink. The size of <tt>lcp</tt> must be at
 * least <tt>length(text)</tt> before calling this function.
 */

	template < 
        typename TLCP,
		typename TText,
		typename TSA,
        typename TAlgSpec >
	inline void _createLCPTableWrapper(
		TLCP &lcp,
		TText const &s,
		TSA const &sa,
		TAlgSpec const &alg,
        True)
	{
	SEQAN_CHECKPOINT
        _createLCPTableRandomAccess(lcp, s, sa, alg);
	}

	template < 
        typename TLCP,
		typename TText,
		typename TSA,
        typename TAlgSpec >
	inline void _createLCPTableWrapper(
		TLCP &lcp,
		TText const &s,
		TSA const &sa,
		TAlgSpec const &alg,
        False)
	{
	SEQAN_CHECKPOINT
        _createLCPTablePipelining(lcp, s, sa, alg);
	}

	template < 
        typename TLCP,
		typename TText,
		typename TSA,
        typename TAlgSpec >
	inline void createLcpTable(
		TLCP &lcp,
		TText const &s,
		TSA const &sa,
		TAlgSpec const &alg)
	{
	SEQAN_CHECKPOINT
        _createLCPTableWrapper(lcp, s, sa, alg, typename LcpCreatorRandomAccess_<TLCP, TText, TSA, TAlgSpec>::Type());
    }


//____________________________________________________________________________



	//////////////////////////////////////////////////////////////////////////////
	// Enhanced LCP Table creation wrappers
	//////////////////////////////////////////////////////////////////////////////

    // build enhanced LCP table with an external pipelining algorithm (ext kasai, ...)
	// and a dynamic programming tree construction alg
	// (in contrast to the LCP table the enhanced LCP table contains the lcp-values 
	// of suffix intervals used in the binary search)
	template < 
		typename TValue, 
        typename TSpec,
		typename TObject, 
        typename TSA,
        typename TLCP,
		typename TAlgSpec >
	void createLcpeTableExt(
		String< TValue, TSpec > &LCPE,
		TObject const &text,
		TSA const &suffixArray,
        TLCP const & /*LCP*/,
		TAlgSpec const)
	{
	SEQAN_CHECKPOINT
		typedef typename Concatenator<TObject>::Type				TConcat;

		// specialization
		typedef Pipe< TConcat, Source<> >						    srcText_t;
		typedef Pipe< TSA, Source<> >					    	    srcSA_t;
	    typedef Pipe< Bundle2< srcText_t, srcSA_t >, TAlgSpec >		creator_t;

		// instantiation and processing
		srcText_t	srcText(concat(text));
		srcSA_t		srcSA(suffixArray);
		creator_t	creator(bundle2(srcText, srcSA));

		#ifdef SEQAN_TEST_INDEX
			isLCPTable(creator, suffixArray, text);
		#endif
		createLcpBinTree(LCPE, creator);
	}

    // build enhanced LCP table with an lcp algorithm
	// and a dynamic programming tree construction alg
    template <
        typename TValue,
        typename TSpec,
        typename TText,
        typename TSA,
        typename TLCP,
		typename TAlgSpec >
    void createLcpeTable(
		String< TValue, TSpec > &LCPE,
		TText const &s,
		TSA const &,
        TLCP const &LCP,
		TAlgSpec const)
	{
	SEQAN_CHECKPOINT
        //TSA LCP;
        //resize(LCP, length(s), Exact());
		// we use LCPE[n-lcpSize..n-1] as a temporary buffer instead of allocating one
        typename Size<TText>::Type lcpSize = length(s) > 1? length(s) - 1: 0;
		typename Suffix<String< TValue, TSpec > >::Type LCPcopy = suffix(LCPE, length(LCPE) - lcpSize);
        LCPcopy = prefix(LCP, lcpSize);
        createLcpBinTree(LCPE, LCP);
    }

    template <
        typename TValue,
        typename TConfig,
        typename TText,
        typename TSA,
        typename TLCP,
		typename TAlgSpec >
    void createLcpeTable(
		String< TValue, External<TConfig> > &LCPE,
		TText const &s,
		TSA const &SA,
        TLCP const &LCP,
		TAlgSpec const alg)
	{
	SEQAN_CHECKPOINT
        createLcpeTableExt(LCPE, s, SA, LCP, alg);
    }

    template <
        typename TValue,
        typename TSpec,
        typename TText,
        typename TSA,
        typename TLCP>
    inline void createLcpeTable(
		String< TValue, TSpec > &LCPE,
		TText &s,
		TSA &SA,
        TLCP &LCP)
	{
	SEQAN_CHECKPOINT
		createLcpeTable(LCPE, s, SA, LCP, Kasai());
    }


//____________________________________________________________________________



	//////////////////////////////////////////////////////////////////////////////
	// Burrows-Wheeler-Table creation wrappers
	//////////////////////////////////////////////////////////////////////////////

	template < typename TBWT,
               typename TText,
			   typename TSA >
    void createBWTableExt(
		TBWT &bwt,
		TText const &s,
		TSA const &SA)
	{
	SEQAN_CHECKPOINT
		// specialization
		typedef Pipe< TText, Source<> >						srcText_t;
		typedef Pipe< TSA, Source<> >   					srcSA_t;
	    typedef Pipe< Bundle2< srcText_t, srcSA_t >, Bwt >	creator_t;

		// instantiation and processing
		srcText_t	srcText(s);
		srcSA_t		srcSA(SA);
		creator_t	creator(bundle2(srcText, srcSA));

		bwt << creator;
	}

/**
.Function.createBWTable:
..summary:Creates a Burrows-Wheeler table from a given text and suffix array.
..cat:Index
..signature:createBWTable(bwt, text, suffixArray[, algo_tag])
..param.bwt:The resulting Burrows-Wheeler table.
..param.text:A given text.
..param.suffixArray:The suffix array of $text$.
..param.algo_tag:A tag that identifies the algorithm which is used for creation.
..remarks:This function should not be called directly. Please use @Function.indexCreate@ or @Function.indexRequire@.
The size of $bwt$ must be at least $length(text)$ before calling this function.
..include:seqan/index.h
*/
/*!
 * @fn createBWTable
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Creates a Burrows-Wheeler table from a given text and suffix array.
 * 
 * @signature createBWTable(bwt, text, suffixArray[, algo_tag])
 * 
 * @param bwt The resulting Burrows-Wheeler table.
 * @param algo_tag A tag that identifies the algorithm which is used for
 *                 creation.
 * @param suffixArray The suffix array of <tt>text</tt>.
 * @param text A given text. Types: @link SequenceConcept @endlink
 * 
 * @section Remarks
 * 
 * This function should not be called directly. Please use @link indexCreate
 * @endlink or @link indexRequire @endlink. The size of <tt>bwt</tt> must be at
 * least <tt>length(text)</tt> before calling this function.
 */
	// default
	template < typename TBWT, typename TText, typename TSA, typename TTextRandom_ >
    inline void _createBWTableWrapper(TBWT &bwt, TText const &s, TSA const &sa,		TTextRandom_ const)
	{
	SEQAN_CHECKPOINT
		createBWTableExt(bwt, concat(s), sa);
	}

	// text supports fast random access
	template < typename TBWT, typename TText, typename TSA >
    inline void _createBWTableWrapper(TBWT &bwt, TText const &s, TSA const &sa,		True const)
	{
	SEQAN_CHECKPOINT
		createBWTableInt(bwt, concat(s), sa);
	}

	template < typename TBWT, typename TText, typename TSA >
    inline void createBWTable(TBWT &bwt, TText const &s, TSA const &sa)
	{
	SEQAN_CHECKPOINT
		_createBWTableWrapper(bwt, s, sa, typename AllowsFastRandomAccess<TText>::Type());
	}


//////////////////////////////////////////////////////////////////////////////

	template <typename TOccValue>
	struct SAValueLess_:
		public ::std::less<TOccValue> {};

	template <typename T1, typename T2, typename TPack>
	struct SAValueLess_< Pair<T1,T2,TPack> >:
		public ::binary_function< Pair<T1,T2,TPack>, Pair<T1,T2,TPack>, bool> 
	{
		inline bool operator()(Pair<T1,T2,TPack> const &a, Pair<T1,T2,TPack> const &b) const {
			return	(getValueI1(a) < getValueI1(b)) ||
					((getValueI1(a) == getValueI1(b)) && (getValueI2(a) < getValueI2(b)));
		}
	};

//TODO(singer): Make this internal
/**
.Function.orderOccurrences:
..summary:Sorts a string of occurrences.
..cat:Index
..signature:orderOccurrences(occstring)
..param.occstring:String of occurrences.
...remarks:Contains suffix array values returned by @Function.getOccurrences@.
..remarks:The occurrences are sorted by increasing positions.
..see:Function.getOccurrences
..see:Metafunction.SAValue
..include:seqan/index.h
*/
/*!
 * @fn orderOccurrences
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Sorts a string of occurrences.
 * 
 * @signature orderOccurrences(occString)
 * 
 * @param occString String of occurrences. Contains suffix array values returned
 *                  by @link getOccurrences @endlink.
 * 
 * @section Remarks
 * 
 * The occurrences are sorted by increasing positions.
 * 
 * Demo: Demo.Mummy
 * 
 * Demo: Demo.Maximal Repeats
 * 
 * Demo: Demo.Maximal Unique Matches
 * 
 * @see getOccurrences
 * @see SAValue
 */
	template <typename TValue, typename TSpec>
	inline void orderOccurrences(String<TValue, TSpec> &occString)
	{
	SEQAN_CHECKPOINT
		::std::sort(begin(occString, Standard()), end(occString, Standard()), SAValueLess_<TValue>());
	}


//////////////////////////////////////////////////////////////////////////////
// fibre creators

/**
.Function.indexCreate:
..summary:Creates a specific @Metafunction.Fibre@.
..cat:Index
..signature:indexCreate(index, fibreTag[, algoTag])
..class:Class.Index
..param.index:The @Class.Index@ object holding the fibre.
...type:Class.Index
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.ESA Index Fibres
...type:Tag.QGram Index Fibres
...type:Tag.WOTD Index Fibres
...type:Tag.FM Index Fibres
..param.algoTag:A tag that identifies the algorithm which is used to create the fibre.
...default:The result of @Metafunction.DefaultIndexCreator@.
..returns:A $bool$ which is $true$ on a successful creation.
..remarks:$indexCreate$ calls the fibre corresponding $createXXX(..)$ function (e.g. @Function.createSuffixArray@).
..include:seqan/index.h
*/
/*!
 * @fn Index#indexCreate
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Creates a specific @link Fibre @endlink.
 * 
 * @signature indexCreate(index, fibreTag[, algoTag])
 * 
 * @param fibreTag A tag that identifies the @link Fibre @endlink
 * @param algoTag A tag that identifies the algorithm which is used to create
 *                 the fibre. Default: The result of @link DefaultIndexCreator
 *                 @endlink.
 * @param index The @link Index @endlink object holding the fibre. Types: @link Index @endlink
 * 
 * @return TReturn A <tt>bool</tt> which is <tt>true</tt> on a successful
 *                 creation.
 * 
 * @section Remarks
 * 
 * <tt>indexCreate</tt> calls the fibre corresponding <tt>createXXX(..)</tt>
 * function (e.g. @link createSuffixArray @endlink).
 */
	template <typename TText, typename TSpec, typename TSpecAlg>
	inline bool indexCreate(Index<TText, TSpec> &index, FibreSA, TSpecAlg const alg) {
	SEQAN_CHECKPOINT
		resize(indexSA(index), length(indexRawText(index)), Exact());
		createSuffixArray(indexSA(index), indexText(index), alg);
		return true;
	}

	template <typename TText, typename TSpec, typename TSpecAlg>
	inline bool indexCreate(Index<TText, TSpec> &index, FibreLcp, TSpecAlg const alg) {
	SEQAN_CHECKPOINT
		resize(indexLcp(index), length(indexRawText(index)), Exact());
		createLcpTable(indexLcp(index), indexText(index), indexSA(index), alg);
		return true;
	}

	template <typename TText, typename TSpec, typename TSpecAlg>
	inline bool indexCreate(Index<TText, TSpec> &index, FibreLcpe, TSpecAlg const alg) {
	SEQAN_CHECKPOINT
	//TODO: separate LCP from LCPE (for now LCPE = LCP + extra)
		resize(indexLcpe(index), sizeofLcpe(lengthSum(index)), Exact());
		createLcpeTable(indexLcpe(index), indexRawText(index), indexSA(index), indexLcp(index), alg);
		return true;
//		return false;
	}

	template <typename TText, typename TSpec>
	inline bool indexCreate(Index<TText, TSpec> &index, FibreBwt, Bwt const) {
	SEQAN_CHECKPOINT
		resize(indexBwt(index), length(indexRawText(index)), Exact());
		createBWTable(indexBwt(index), indexText(index), indexRawSA(index));
		return true;
	}

	template <typename TText, typename TSpec>
	inline bool indexCreate(Index<TText, TSpec> &index, FibreChildtab, Childtab const) {
	SEQAN_CHECKPOINT
		resize(indexChildtab(index), length(indexRawText(index)), Exact());
		createChildtab(indexChildtab(index), indexLcp(index));
		return true;
	}

	template <typename TText, typename TSpec, typename TFibre>
	inline bool indexCreate(Index<TText, TSpec> &index, Tag<TFibre> const fibre) {
	SEQAN_CHECKPOINT
		return indexCreate(index, fibre, typename DefaultIndexCreator<Index<TText, TSpec>, Tag<TFibre> const>::Type());
	}


//////////////////////////////////////////////////////////////////////////////
// automatic fibre creation

/**
.Function.indexSupplied:
..summary:Returns whether a specific @Metafunction.Fibre@ is present.
..cat:Index
..signature:indexSupplied(index, fibreTag)
..class:Class.Index
..param.index:The @Class.Index@ object holding the fibre.
...type:Class.Index
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.ESA Index Fibres
..returns:A $bool$ which is $true$, iff the fibre is present.
..include:seqan/index.h
*/
/*!
 * @fn Index#indexSupplied
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns whether a specific @link Fibre @endlink is present.
 * 
 * @signature indexSupplied(index, fibreTag)
 * 
 * @param index The @link Index @endlink object holding the fibre. Types: @link Index @endlink
 * @param fibreTag A tag that identifies the @link Fibre @endlink. 
 *                 Index Fibres
 * 
 * @return TReturn A <tt>bool</tt> which is <tt>true</tt>, iff the fibre is
 *                 present.
 */
	template <typename TText, typename TSpec, typename TFibre>
	inline bool indexSupplied(Index<TText, TSpec> &index, Tag<TFibre> const fibre) {
	SEQAN_CHECKPOINT
		return !empty(getFibre(index, fibre));
	}


//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexRequire:
..summary:On-demand creation of a specific @Metafunction.Fibre@.
..cat:Index
..signature:indexRequire(index, fibreTag)
..class:Class.Index
..param.index:The @Class.Index@ object holding the fibre.
...type:Class.Index
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.ESA Index Fibres
...type:Tag.QGram Index Fibres
...type:Tag.WOTD Index Fibres
...type:Tag.FM Index Fibres
..returns:A $bool$ which is $true$ on a successful creation.
..remarks:If the fibre already exists (@Function.indexSupplied@ is true) then $indexRequire$ does nothing.
If the fibre doesn't exist then @Function.indexCreate@ is called to create it.
..include:seqan/index.h
..example
...text:The following code shows how the BWT of an text can be computed.
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
I	ISSISSIPPI
*/
/*!
 * @fn Index#indexRequire
 * 
 * @headerfile seqan/index.h
 * 
 * @brief On-demand creation of a specific @link Fibre @endlink.
 * 
 * @signature indexRequire(index, fibre_tag)
 * 
 * @param index The @link Index @endlink object holding the fibre. Types: @link Index @endlink
 * @param fibre_tag A tag that identifies the @link Fibre @endlink
 * 
 * @return TReturn A <tt>bool</tt> which is <tt>true</tt> on a successful
 *                 creation.
 * 
 * @section Remarks
 * 
 * If the fibre already exists (@link indexSupplied @endlink is true) then
 * <tt>indexRequire</tt> does nothing. If the fibre doesn't exist then @link
 * indexCreate @endlink is called to create it.
 */

	template <typename TText, typename TSpec, typename TFibre>
	inline bool indexRequire(Index<TText, TSpec> &index, Tag<TFibre> const fibre) {
	SEQAN_CHECKPOINT
		if (indexSupplied(index, fibre)) return true;				// if the table doesn't exist,
		if (!indexSolveDependencies(index, fibre)) return false;	// fulfill requirements
		return indexCreate(index, fibre);							// and create table
	}


//////////////////////////////////////////////////////////////////////////////
// index cargo interface

	template <typename TText, typename TSpec>
	inline typename Reference< typename Cargo<Index<TText, TSpec> >::Type >::Type
	cargo(Index<TText, TSpec> & me)
	{
	SEQAN_CHECKPOINT
		return me.cargo;
	}

	template <typename TText, typename TSpec>
	inline typename Reference< typename Cargo<Index<TText, TSpec> const>::Type >::Type
	cargo(Index<TText, TSpec> const & me)
	{
	SEQAN_CHECKPOINT
		return me.cargo;
	}

//////////////////////////////////////////////////////////////////////////////
// solve dependencies

	template <typename TText, typename TSpec, typename TFibre>
	inline bool indexSolveDependencies(Index<TText, TSpec> &, Tag<TFibre> const) {
	SEQAN_CHECKPOINT
		return true;
	}

	template <typename TText, typename TSpec>
	inline bool indexSolveDependencies(Index<TText, TSpec> &index, FibreLcp) {
	SEQAN_CHECKPOINT
		return indexRequire(index, FibreSA());
	}

	template <typename TText, typename TSpec>
	inline bool indexSolveDependencies(Index<TText, TSpec> &index, FibreLcpe) {
	SEQAN_CHECKPOINT
		return indexRequire(index, FibreLcp());
	}

	template <typename TText, typename TSpec>
	inline bool indexSolveDependencies(Index<TText, TSpec> &index, FibreChildtab) {
	SEQAN_CHECKPOINT
		return indexRequire(index, FibreLcp());
	}

	template <typename TText, typename TSpec>
	inline bool indexSolveDependencies(Index<TText, TSpec> &index, FibreBwt) {
	SEQAN_CHECKPOINT
		return indexRequire(index, FibreSA());
	}


//////////////////////////////////////////////////////////////////////////////
// open

	template < typename TValue, typename TSpec >
	inline bool open(String<TValue, TSpec> &string, const char *fileName, int openMode) {
	SEQAN_CHECKPOINT
		String<TValue, External< ExternalConfigLarge<> > > extString;
		if (!open(extString, fileName, openMode & ~OPEN_CREATE)) return false;
		assign(string, extString, Exact());
		return true;
	}
	template < typename TValue, typename TSpec >
	inline bool open(String<TValue, TSpec> &string, const char *fileName) {
	SEQAN_CHECKPOINT
		return open(string, fileName, OPEN_RDONLY);
	}

	template < typename THost, typename TSpec >
	inline bool open(Segment<THost, TSpec> &string, const char *fileName, int openMode) {
	SEQAN_CHECKPOINT
		String<typename Value<THost>::Type, External< ExternalConfigLarge<> > > extString;
		if (!open(extString, fileName, openMode & ~OPEN_CREATE)) return false;
		assign(string, extString, Exact());
		return true;
	}
	template < typename THost, typename TSpec >
	inline bool open(Segment<THost, TSpec> &string, const char *fileName) {
	SEQAN_CHECKPOINT
		return open(string, fileName, OPEN_RDONLY);
	}

	// ATTENTION:
	// This implementation of open doesn't work with external memory StringSets (External<>, MMap<>)
	// If you need a persistent external StringSet you have to use a Owner<ConcatDirect<> > StringSet.
	template < typename TString, typename TSSSpec >
	inline bool open(StringSet<TString, TSSSpec> &multi, const char *fileName, int openMode) {
	SEQAN_CHECKPOINT
		char id[12]; // 2^32 has 10 decimal digits + 1 (0x00)
		unsigned i = 0;
		clear(multi);
		CharString name;
		while (true)
		{
			sprintf(id, ".%u", i);
			name = fileName;
			append(name, id);
			{
				resize(multi, i + 1);
				if (!open(multi[i], toCString(name), (openMode & ~OPEN_CREATE) | OPEN_QUIET))
				{
					resize(multi, i);
					break;
				}
			}
			++i;
		}
		return i > 0;
	}

	template < typename TValue, typename TSpec, typename TSSSpec >
        inline bool open(StringSet<String<TValue, TSpec>, Dependent<TSSSpec> > &, const char *, int) {
        SEQAN_CHECKPOINT
        // Do nothing for dependent string sets
        return true;
    }

	template < typename TString, typename TSSSpec >
	inline bool open(StringSet<TString, Owner<ConcatDirect<TSSSpec> > > &multi, const char *fileName, int openMode) {
	SEQAN_CHECKPOINT
		CharString name;
		name = fileName;
		append(name, ".concat");
		if (!open(multi.concat, toCString(name), openMode | OPEN_QUIET)) return false;
		name = fileName;
		append(name, ".limits");
		if (!open(multi.limits, toCString(name), openMode | OPEN_QUIET) && !empty(multi.concat))
		{
			clear(multi);
			return false;
		}
		// limits file was just created
		if (empty(multi.limits))
			appendValue(multi.limits, 0);
		return true;
	}

	template < typename TValue, typename TSpec, typename TSSSpec>
	inline bool open(StringSet<String<TValue, TSpec>, TSSSpec> &multi, const char *fileName) {
	SEQAN_CHECKPOINT
		return open(multi, fileName, OPEN_RDONLY);
	}


//////////////////////////////////////////////////////////////////////////////
// save

	template < typename TValue, typename TSpec >
	inline bool save(String<TValue, TSpec> const &string, const char *fileName, int openMode) {
	SEQAN_CHECKPOINT
//
//		if (length(string) == 0) return true;
		String<TValue, External< ExternalConfigLarge<> > > extString;
		if (!open(extString, fileName, openMode)) return false;
		assign(extString, string, Exact());
		return true;
	}
	template < typename TValue, typename TSpec >
	inline bool save(String<TValue, TSpec> const &string, const char *fileName) {
	SEQAN_CHECKPOINT
		return save(string, fileName, OPEN_WRONLY | OPEN_CREATE);
	}

	template < typename THost, typename TSpec >
	inline bool save(Segment<THost, TSpec> const &string, const char *fileName, int openMode) {
	SEQAN_CHECKPOINT
		if (length(string) == 0) return true;
		String<typename Value<THost>::Type, External< ExternalConfigLarge<> > > extString;
		if (!open(extString, fileName, openMode)) return false;
		assign(extString, string, Exact());
		return true;
	}
	template < typename THost, typename TSpec >
	inline bool save(Segment<THost, TSpec> const &string, const char *fileName) {
	SEQAN_CHECKPOINT
		return save(string, fileName, OPEN_WRONLY | OPEN_CREATE);
	}

	template < typename TString, typename TSSSpec>
	inline bool save(StringSet<TString, TSSSpec> const &multi, const char *fileName, int openMode) {
	SEQAN_CHECKPOINT
		if (length(multi) == 0) return true;
		char id[12]; // 2^32 has 10 decimal digits + 2 ('.' and 0x00)
		CharString name;
		for(unsigned i = 0; i < length(multi); ++i)
		{
			sprintf(id, ".%u", i);
			name = fileName;
			append(name, &(id[0]));
			if (!save(multi[i], toCString(name), openMode))
				return false;
		}
		return true;
	}

	template < typename TValue, typename TSpec, typename TSSSpec >
    inline bool save(StringSet<String<TValue, TSpec>, Dependent<TSSSpec> > const &, const char *, int) {
        SEQAN_CHECKPOINT
        // Do nothing for dependent string sets
        return true;
    }

	template < typename TString, typename TSSSpec >
	inline bool save(StringSet<TString, Owner<ConcatDirect<TSSSpec> > > const &multi, const char *fileName, int openMode) {
	SEQAN_CHECKPOINT
		CharString name;
		name = fileName;
		append(name, ".concat");
		if (!save(multi.concat, toCString(name), openMode)) return false;
		name = fileName;
		append(name, ".limits");
		if (!save(multi.limits, toCString(name), openMode)) return false;
		return true;
	}
	template < typename TValue, typename TSpec, typename TSSSpec>
	inline bool save(StringSet<String<TValue, TSpec>, TSSSpec> const &multi, const char *fileName) {
	SEQAN_CHECKPOINT
		return save(multi, fileName, OPEN_WRONLY | OPEN_CREATE);
	}

}

#endif
