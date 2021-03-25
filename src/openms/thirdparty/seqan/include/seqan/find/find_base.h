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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Definition of the class Finder and supporting tags and metafunctions.
// ==========================================================================

#ifndef SEQAN_HEADER_FIND_BASE_H
#define SEQAN_HEADER_FIND_BASE_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Tags

/**
.Tag.FindInfix:
..cat:Searching
..summary:Find needle as a substring of haystack (infix search).
..see:Tag.FindPrefix
..see:Spec.Myers
..see:Spec.DPSearch
..include:seqan/find.h
 */
struct FindInfix;

	
/**
.Tag.FindPrefix:
..cat:Searching
..summary:Find needle as a prefix of the haystack (prefix serach)
..see:Tag.FindInfix
..see:Spec.Myers
..see:Spec.DPSearch
..include:seqan/find.h
*/
struct FindPrefix;

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.DefaultFinder:
..cat:Searching
..summary:Default @Class.Finder@ specialization type.
..signature:DefaultFinder<THaystack>::Type
..param.THaystack:The given haystack type.
..returns:Is $void$ by default and @Tag.Index Find Algorithm.EsaFindMlr@ if $THaystack$ is an @Class.Index@.
..include:seqan/find.h
*/
template < typename TObject >
struct DefaultFinder 
{
	typedef void Type;
};

/**
.Metafunction.DefaultPattern:
..cat:Searching
..summary:Default @Class.Pattern@ specialization type.
..signature:DefaultPattern<TNeedle>::Type
..param.TNeedle:The given needle type.
..returns:Is $void$ by default.
..include:seqan/find.h
*/
template < typename TObject >
struct DefaultPattern 
{
	typedef void Type;
};

/**
.Metafunction.Haystack:
..summary:Returns the haystack type of a @Class.Finder@ type.
..cat:Searching
..signature:Haystack<TFinder>::Type
..class:Class.Finder
..param.TFinder:A @Class.Finder@ type.
...type:Class.Finder
..returns:The haystack type of $TFinder$, i.e. $THaystack$ for $Finder<THaystack, TSpec>$.
This is an alias to function @Function.host@ of the pattern function.
..see:Function.host
..include:seqan/find.h
*/

template <typename TFinder>
struct Haystack 
{
	typedef typename Container<TFinder>::Type Type;
};

/**
.Metafunction.Needle:
..summary:Returns the needle type of a @Class.Pattern@ type.
..cat:Searching
..signature:Needle<TPattern>::Type
..class:Class.Pattern
..param.TPattern:A @Class.Pattern@ type.
...type:Class.Pattern
..returns:The needle type of $TPattern$, i.e. $TNeedle$ for $Pattern<TNeedle, TSpec>$.
..include:seqan/find.h
*/

template <typename TPattern>
struct Needle 
{
	typedef typename Host<TPattern>::Type Type;
};

template <typename THost, typename TSpec>
struct Needle<Segment<THost, TSpec> > 
{
	typedef Segment<THost, TSpec> Type;
};

template <typename THost, typename TSpec>
struct Needle<Segment<THost, TSpec> const> 
{
	typedef Segment<THost, TSpec> const Type;
};


//////////////////////////////////////////////////////////////////////////////

/**
.Function.find:
..summary:Search for a @Class.Pattern@ in a @Class.Finder@ object.
..cat:Searching
..signature:find(finder, pattern)
..signature:find(finder, pattern, k)
..class:Class.Finder
..param.finder:The @Class.Finder@ object to search through.
...remarks:For online-algorithm $patterns$, finder can also be an arbitrary @Concept.RootedIteratorConcept|Rooted Iterator@.
...type:Class.Finder
...type:Concept.RootedIteratorConcept
..param.pattern:The @Class.Pattern@ object to search for.
...remarks:For index $finders$, pattern can also be a Sequence.
...type:Class.Pattern
..param.k:Desired minimal score (for approximate matching).
...remarks:$k$ has to be a number <= 0.
...remarks:Differences are deletions, insertions and substitutions.
..returns:$boolean$ that indicates whether an occurrence of $pattern$ was found or not.
..remarks:Repeated calls of this function iterate through all occurrences of $pattern$.
..include:seqan/find.h
*/

/**
.Class.Finder:
..summary:Holds the haystack and a current search context.
..cat:Searching
..signature:Finder<THaystack[, TSpec]>
..param.THaystack:The haystack type.
...type:Class.String
...type:Class.Index
..param.TSpec:The index-algorithm to search with (Optional).
...default:The result of @Metafunction.DefaultFinder@
...remarks:Leave empty for online pattern matching (see @Class.Pattern@).
...remarks:If $THaystack$ is an @Class.Index@, then $TSpec$ specifies the index search algorithm.
..remarks:$position(finder)$ returns the position of the current hit in the haystack.
If $THaystack$ is a set of strings or an index of a set of strings, then $position(finder)$ returns a @Class.Pair@ $(hayNo, pos)$,
in which $hayNo$ is the haystack index and $pos$ the local position of the hit.
..remarks:To reset the finder object and use it on another text or different text position, use $clear(finder)$
Note that $clear(finder)$ doesn't move the text iterator. To start the search from the beginning or somewhere else in the text, use
@Function.goBegin@ or @Function.setPosition@.
..example
...text:The following example shows how one can search online for a pattern in a haystack. Note that it is neccessary to reset the finder befor searching for another pattern.
...file:demos/find/finder_online.cpp
...output:Hit at position: 4
Hit at position: 12
Hit at position: 22
Hit at position: 43
Hit at position: 8
Hit at position: 16
Hit at position: 17
Hit at position: 29
Hit at position: 35
Hit at position: 38
...text:In contrast to the example above the code below shows how one can use a Finder with an index as base. Again, note that it is neccessary to reset the finder befor searching for another pattern.
...file:demos/find/finder_index.cpp
...output:Hit at position: 12
Hit at position: 4
Hit at position: 22
Hit at position: 43
Hit at position: 38
Hit at position: 8
Hit at position: 35
Hit at position: 29
Hit at position: 17
Hit at position: 16
..include:seqan/find.h
*/

/*!
 * @class Finder
 * 
 * @headerfile seqan/find.h
 * 
 * @brief Holds the haystack and a current search context.
 * 
 * @signature Finder<THaystack[, TSpec]>
 * 
 * @tparam TSpec The index-algorithm to search with (Optional).Leave empty for
 *               online pattern matching (see @link Pattern @endlink).If
 *               <tt>THaystack</tt> is an @link Index @endlink, then
 *               <tt>TSpec</tt> specifies the index search algorithm. Types:
 *               Pigeonhole, Swift, Backtracking Default: The result of @link
 *               DefaultFinder @endlink
 * @tparam THaystack The haystack type. Types: String, Index
 * 
 * @section Remarks
 * 
 * <tt>position(finder)</tt> returns the position of the current hit in the
 * haystack. If <tt>THaystack</tt> is a set of strings or an index of a set of
 * strings, then <tt>position(finder)</tt> returns a @link Pair @endlink
 * <tt>(hayNo, pos)</tt>, in which <tt>hayNo</tt> is the haystack index and
 * <tt>pos</tt> the local position of the hit.
 * 
 * To reset the finder object and use it on another text or different text
 * position, use <tt>clear(finder)</tt> Note that <tt>clear(finder)</tt> doesn't
 * move the text iterator. To start the search from the beginning or somewhere
 * else in the text, use @link goBegin @endlink or @link setPosition @endlink.
 * 
 * @section Examples
 * 
 * The following example shows how to restart a search from the beginning of a
 * text.
 * 
 * @code{.cpp}
 * CharString hstck = "I spy with my little eye something that is yellow";
 * Finder<CharString> finder(hstck);
 *  
 * Pattern<CharString, Horspool> p1("y");
 * findAll(finder, p1);
 *  
 * goBegin(finder);    // move Finder to the beginning of the text
 * clear(finder);      // reset Finder
 *  
 * Pattern<CharString, Horspool> p2("t");
 * findAll(finder, p2);
 * @endcode
 * Demo: Demo.Index Finder StringSet
 * 
 * Demo: Demo.Index Finder
 */

///.Function.clear.param.object.type:Class.Finder
///.Function.clear.class:Class.Finder
///.Function.position.param.iterator.type:Class.Finder
///.Function.position.class:Class.Finder

template < typename THaystack, typename TSpec = typename DefaultFinder<THaystack>::Type >
class Finder
{
public:
	typedef typename Iterator<THaystack, Rooted>::Type TIterator;
	typedef typename Position<THaystack>::Type TPosition;
	typedef typename Size<THaystack>::Type TSize;

	TIterator data_iterator;
	TPosition data_endPos; //note: we need this since iterator could point to begin or end (depending on pattern type)
	TSize data_length;
	bool _needReinit;					// if true, the Pattern needs to be reinitialized
	bool _beginFind_called;					// if false, then findBegin was not yet called for this match position (see findBegin default implementation)

/**
.Memfunc.Finder#Finder:
..class:Class.Finder
..summary:Constructor
..signature:Finder()
 */
	Finder()
		: data_endPos(0)
		, data_length(0)
		, _needReinit(true)
		, _beginFind_called(false)
	{}

	/**
.Memfunc.Finder#Finder:
..signature:Finder(haystack)
..param.haystack:The haystack to work on, $THaystack$.
	 */
	Finder(THaystack & haystack)
		: data_iterator(begin(haystack, Rooted()))
		, data_endPos(0)
		, data_length(0)
		, _needReinit(true) 
		, _beginFind_called(false)
	{}

	/**
.Memfunc.Finder#Finder:
..signature:Finder(iter)
..param.iter:The iterator to work on, either const or non-const.
	 */
	Finder(TIterator &iter)
		: data_iterator(iter)
		, data_endPos(0)
		, data_length(0)
		, _needReinit(true) 
		, _beginFind_called(false)
	{}

	Finder(TIterator const &iter)
		: data_iterator(iter)
		, data_endPos(0)
		, data_length(0)
		, _needReinit(true) 
		, _beginFind_called(false)
	{}
	
	/**
.Memfunc.Finder#Finder:
..signature:Finder(orig)
..param.orig:Finder object to copy from (copy constructor).
...type:Class.Finder
	 */
	Finder(Finder const &orig)
		: data_iterator(orig.data_iterator)
		, data_endPos(orig.data_endPos)
		, data_length(orig.data_length)
		, _needReinit(orig._needReinit) 
		, _beginFind_called(orig._beginFind_called)
	{}

	~Finder() {}

//____________________________________________________________________________

	Finder const &
	operator = (Finder const & other)
	{
		data_iterator = other.data_iterator;
		data_endPos = other.data_endPos;
		data_length = other.data_length;
		_needReinit = other._needReinit;
		_beginFind_called = other._beginFind_called;
		return *this;
	}

//____________________________________________________________________________

	inline typename Reference<TIterator>::Type 
	operator* () 
	{
SEQAN_CHECKPOINT
		return value(hostIterator(*this));
	}

	inline typename Reference<TIterator const>::Type 
	operator* () const
	{
SEQAN_CHECKPOINT
		return value(hostIterator(*this));
	}

//____________________________________________________________________________

	operator TIterator () const
	{
SEQAN_CHECKPOINT
		return data_iterator;
	}

//____________________________________________________________________________

};

//////////////////////////////////////////////////////////////////////////////

template <typename T>
inline void
_setFinderEnd(T &) {}

template <typename T, typename TPosition>
inline void
_setFinderEnd(T &, TPosition) {}


template <typename THaystack, typename TSpec>
inline void
_setFinderEnd(Finder<THaystack, TSpec> & me)
{//shortcut: move end position to iterator position +1
	me._beginFind_called = false;
	me.data_endPos = position(me)+1;
}
template <typename THaystack, typename TSpec, typename TPosition>
inline void
_setFinderEnd(Finder<THaystack, TSpec> & me,
			  TPosition end_pos)
{
	me._beginFind_called = false;
	me.data_endPos = end_pos;
}

//____________________________________________________________________________

template <typename T, typename TSize>
inline void
_setFinderLength(T &, TSize) {}

template <typename THaystack, typename TSpec, typename TSize>
inline void
_setFinderLength(Finder<THaystack, TSpec> & me,
				 TSize _length)
{
	me.data_length = _length;
}

//////////////////////////////////////////////////////////////////////////////
///.Function.beginPosition.param.object.type:Class.Finder
///.Function.beginPosition.class:Class.Finder

template <typename THaystack, typename TSpec>
inline typename Position<THaystack>::Type 
beginPosition(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	return me.data_endPos - me.data_length;
}

template <typename THaystack, typename TSpec>
inline typename Position<THaystack const>::Type
beginPosition(Finder<THaystack, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return me.data_endPos - me.data_length;
}

//////////////////////////////////////////////////////////////////////////////
///.Function.begin.param.object.type:Class.Finder
///.Function.begin.class:Class.Finder

template <typename THaystack, typename TSpec, typename TTag>
inline typename Iterator<THaystack, Tag<TTag> const>::Type
begin(Finder<THaystack, TSpec> & me,
	  Tag<TTag> const tag)
{
SEQAN_CHECKPOINT
	return iter(haystack(me), beginPosition(me), tag);
}

template <typename THaystack, typename TSpec, typename TTag>
inline typename Iterator<THaystack const, Tag<TTag> const>::Type
begin(Finder<THaystack, TSpec> const & me,
	  Tag<TTag> const tag)
{
SEQAN_CHECKPOINT
	return iter(haystack(me), beginPosition(me), tag);
}

//////////////////////////////////////////////////////////////////////////////
///.Function.endPosition.param.object.type:Class.Finder
///.Function.endPosition.class:Class.Finder

template <typename THaystack, typename TSpec>
inline typename Position<THaystack>::Type 
endPosition(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	return me.data_endPos;
}

template <typename THaystack, typename TSpec>
inline typename Position<THaystack const>::Type
endPosition(Finder<THaystack, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return me.data_endPos;
}

//////////////////////////////////////////////////////////////////////////////
///.Function.end.param.object.type:Class.Finder
///.Function.end.class:Class.Finder

template <typename THaystack, typename TSpec, typename TTag>
inline typename Iterator<THaystack, Tag<TTag> const>::Type
end(Finder<THaystack, TSpec> & me,
	Tag<TTag> const tag)
{
SEQAN_CHECKPOINT
	return iter(haystack(me), endPosition(me), tag);
}

template <typename THaystack, typename TSpec, typename TTag>
inline typename Iterator<THaystack const, Tag<TTag> const>::Type
end(Finder<THaystack, TSpec> const & me,
	Tag<TTag> const tag)
{
SEQAN_CHECKPOINT
	return iter(haystack(me), endPosition(me), tag);
}

//////////////////////////////////////////////////////////////////////////////
///.Function.length.param.object.type:Class.Finder
///.Function.length.class:Class.Finder

template <typename THaystack, typename TSpec>
inline typename Size<THaystack>::Type 
length(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	return me.data_length;
}
template <typename THaystack, typename TSpec>
inline typename Size<THaystack const>::Type 
length(Finder<THaystack, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return me.data_length;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.Finder#infix:
..summary:Returns the segment of the last found match in haystack.
..signature:Infix infix(finder)
..param.finder:An online finder.
...type:Class.Finder
..returns:An @Metafunction.Infix.infix@ of the @Metafunction.Haystack.haystack@.
...metafunction:Metafunction.Infix
..remarks:This function works only correct if the begin position of the match was already found,
see @Function.findBegin@.
..include:seqan/find.h
*/
template <typename THaystack, typename TSpec>
inline typename Infix<THaystack>::Type
infix(Finder<THaystack, TSpec> & me)
{
	return infix(haystack(me), beginPosition(me), endPosition(me));
}

template <typename THaystack, typename TSpec>
inline typename Infix<THaystack const>::Type
infix(Finder<THaystack, TSpec> const & me)
{
	return infix(haystack(me), beginPosition(me), endPosition(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename THaystack, typename TSpec>
inline typename Parameter_<THaystack>::Type 
host(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	return container(hostIterator(me));
}

template <typename THaystack, typename TSpec>
inline typename Parameter_<THaystack>::Type 
host(Finder<THaystack, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return container(hostIterator(me));
}

template <typename THaystack, typename TSpec>
inline typename Parameter_<THaystack>::Type 
container(Finder<THaystack, TSpec> & me)
{
	return container(hostIterator(me));
}

template <typename THaystack, typename TSpec>
inline typename Parameter_<THaystack>::Type 
container(Finder<THaystack, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return container(hostIterator(me));
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline void
setHost(Finder<THaystack, TSpec> & me, 
		typename Parameter_<THaystack>::Type container_)
{
SEQAN_CHECKPOINT
	setContainer(hostIterator(me), container_);
	goBegin(me);
}

template <typename THaystack, typename TSpec>
inline void
setContainer(Finder<THaystack, TSpec> & me, 
			 typename Parameter_<THaystack>::Type container_)
{
SEQAN_CHECKPOINT
	setContainer(hostIterator(me), container_);
	goBegin(me);
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline typename Iterator<THaystack, Rooted>::Type &
hostIterator(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	return me.data_iterator;
}

template <typename THaystack, typename TSpec>
inline typename Iterator<THaystack, Rooted>::Type const &
hostIterator(Finder<THaystack, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return me.data_iterator;
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline bool
empty(Finder<THaystack, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return me._needReinit;
}

template <typename THaystack, typename TSpec>
inline void
clear(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	me._needReinit = true;
}

//____________________________________________________________________________

template <typename T>
inline void
_finderSetNonEmpty(T & me)
{
SEQAN_CHECKPOINT
	goBegin(me);
}


template <typename THaystack, typename TSpec>
inline void
_finderSetNonEmpty(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	me._needReinit = false;
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline bool
atBegin(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	return (!empty(me) && atBegin(hostIterator(me)));
}

template <typename THaystack, typename TSpec>
inline bool
atEnd(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	return (!empty(me) && atEnd(hostIterator(me)));
}

//____________________________________________________________________________

/**
.Function.goBegin
..signature:goBegin(finder)
..class:Class.Finder
...param.finder:Finder object to go to beginning in.
..include:seqan/find.h
 */
template <typename THaystack, typename TSpec>
inline void
goBegin(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	//_finderSetNonEmpty(me);
	goBegin(hostIterator(me));
}

/**
.Function.goEnd
..signature:goEnd(finder)
..class:Class.Finder
...param.finder:Finder object to go to end in.
..include:seqan/find.h
 */
template <typename THaystack, typename TSpec>
inline void
goEnd(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	//_finderSetNonEmpty(me);
	goEnd(hostIterator(me));
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, TSpec> >::Type
position(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	if (empty(me)) return 0;
	return position(hostIterator(me));
}

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, TSpec> >::Type
position(Finder<THaystack, TSpec> const & me)
{
SEQAN_CHECKPOINT
	if (empty(me)) return 0;
	return position(hostIterator(me));
}

//____________________________________________________________________________
/**
.Function.setPosition:
..cat:Searching
..summary:Sets the position of a finder.
..signature:setPosition(finder, pos)
..class:Class.Finder
..param.finder:A finder.
...type:Class.Finder
..param.pos:A position.
...metafunction:Metafunction.Position
..see:Function.position
..include:seqan/find.h
*/

template <typename THaystack, typename TSpec, typename TPosition>
inline void 
setPosition(Finder<THaystack, TSpec> & me, TPosition pos_)
{
SEQAN_CHECKPOINT
	setPosition(hostIterator(me), pos_);
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline Finder<THaystack, TSpec> &
operator--(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	--hostIterator(me);
	return me;
}

template <typename THaystack, typename TSpec>
inline Finder<THaystack, TSpec> &
operator++(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
/*			if (beforeBegin()) {
		goBegin(hostIterator(me));
	} else*/
		++hostIterator(me);
	return me;
}

//////////////////////////////////////////////////////////////////////////////
// operator +
//////////////////////////////////////////////////////////////////////////////

template <typename THaystack, typename TSpec, typename TIntegral>
inline Finder<THaystack, TSpec> const
operator + (Finder<THaystack, TSpec> const & left, TIntegral right)
{
SEQAN_CHECKPOINT
	return Finder<THaystack, TSpec>(hostIterator(left) + right);
}

//////////////////////////////////////////////////////////////////////////////
// operator +=
//////////////////////////////////////////////////////////////////////////////

template <typename THaystack, typename TSpec, typename TIntegral>
inline Finder<THaystack, TSpec> &
operator += (Finder<THaystack, TSpec> & left,
				TIntegral right)
{
SEQAN_CHECKPOINT
	hostIterator(left) += right;
	return left;
}

//////////////////////////////////////////////////////////////////////////////
// operator -
//////////////////////////////////////////////////////////////////////////////

template <typename THaystack, typename TSpec, typename TIntegral>
inline Finder<THaystack, TSpec> const
operator - (Finder<THaystack, TSpec> const & left, TIntegral right)
{
SEQAN_CHECKPOINT
	return Finder<THaystack, TSpec>(hostIterator(left) - right);
}

template <typename THaystack, typename TSpec, typename TIntegral>
inline typename Difference<Finder<THaystack, TSpec> const>::Type
operator - (Finder<THaystack, TSpec> const & left, Finder<THaystack, TSpec> const & right)
{
SEQAN_CHECKPOINT
	return hostIterator(left) - hostIterator(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator -=
//////////////////////////////////////////////////////////////////////////////

template <typename THaystack, typename TSpec, typename TIntegral>
inline Finder<THaystack, TSpec> &
operator -= (Finder<THaystack, TSpec> & left,
				TIntegral right)
{
SEQAN_CHECKPOINT
	hostIterator(left) -= right;
	return left;
}

//____________________________________________________________________________


/**
.Function.setHaystack:
..summary:Sets the haystack of a @Class.Finder@ object.
..cat:Searching
..signature:setHaystack(finder, haystack)
..class:Class.Finder
..param.finder:The @Class.Finder@ object to search with.
...type:Class.Finder
..param.haystack:The haystack object the finder searches through.
...type:Class.String
..include:seqan/find.h
*/

template < typename THaystack, typename TSpec >
inline void
setHaystack(Finder<THaystack, TSpec> &obj, THaystack const &hstk) {
	setHost(obj, hstk);
}

/**
.Function.haystack:
..summary:Returns the haystack of a @Class.Finder@ object.
..cat:Searching
..signature:haystack(finder)
..class:Class.Finder
..param.finder:The @Class.Finder@ object to search through.
...type:Class.Finder
..returns:The haystack object.
..remarks:The result type is @Metafunction.Haystack@$<TFinder>::Type$ for finder of type $TFinder$.
..include:seqan/find.h
*/

template < typename TObject >
inline typename Parameter_<typename Haystack<TObject>::Type>::Type
haystack(TObject &obj) {
	return container(obj);
}

template < typename TObject >
inline typename Parameter_<typename Haystack<TObject const>::Type>::Type
haystack(TObject const &obj) {
	return container(obj);
}



//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Container.param.T.type:Class.Finder
///.Metafunction.Container.class:Class.Finder
template <typename THaystack, typename TSpec>
struct Container< Finder<THaystack, TSpec> > {
	typedef THaystack Type;
};

template <typename THaystack, typename TSpec>
struct Container< Finder<THaystack, TSpec> const> {
	typedef THaystack const Type;
};


///.Metafunction.Host.param.T.type:Class.Finder
///.Metafunction.Host.class:Class.Finder
template <typename THaystack, typename TSpec>
struct Host< Finder<THaystack, TSpec> > {
	typedef THaystack Type;
};

template <typename THaystack, typename TSpec>
struct Host< Finder<THaystack, TSpec> const> {
	typedef THaystack const Type;
};


///.Metafunction.Value.param.T.type:Class.Finder
///.Metafunction.Value.class:Class.Finder
template <typename THaystack, typename TSpec>
struct Value< Finder<THaystack, TSpec> > {
	typedef typename Value<THaystack>::Type Type;
};

template <typename THaystack, typename TSpec>
struct Position< Finder<THaystack, TSpec> >:
	Position<THaystack> {};


///.Metafunction.Difference.param.T.type:Class.Finder
///.Metafunction.Difference.class:Class.Finder
template <typename THaystack, typename TSpec>
struct Difference< Finder<THaystack, TSpec> > {
	typedef typename Difference<THaystack>::Type Type;
};

template <typename THaystack, typename TSpec>
struct Size< Finder<THaystack, TSpec> > {
	typedef typename Size<THaystack>::Type Type;
};


///.Metafunction.Iterator.param.T.type:Class.Finder
///.Metafunction.Iterator.class:Class.Finder
template <typename THaystack, typename TSpec, typename TIteratorSpec>
struct Iterator< Finder<THaystack, TSpec>, TIteratorSpec >
{
	typedef typename Iterator<THaystack>::Type Type;
};

template <typename THaystack, typename TSpec, typename TIteratorSpec>
struct Iterator< Finder<THaystack, TSpec> const, TIteratorSpec >
{
	typedef typename Iterator<THaystack>::Type Type;
};


// TODO(holtgrew): Document DefaultGetIterator at main location, first.
// .Metafunction.DefaultGetIterator.param.T.type:Class.Finder
template <typename THaystack, typename TSpec>
struct DefaultGetIteratorSpec< Finder<THaystack, TSpec> >:
	DefaultGetIteratorSpec< THaystack >
{
};

template <typename THaystack, typename TSpec>
struct DefaultGetIteratorSpec< Finder<THaystack, TSpec> const>:
	DefaultGetIteratorSpec< THaystack const>
{
};

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
