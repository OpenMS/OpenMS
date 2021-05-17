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

 //TODO(weese:) sync new documentation to old dddoc updates

#ifndef SEQAN_HEADER_SHAPE_BASE_H
#define SEQAN_HEADER_SHAPE_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

	template <unsigned q>
	struct UngappedShape {};
	typedef UngappedShape<0> SimpleShape;

	template <typename TSpec>
	struct GappedShape {};
	typedef GappedShape<Default> GenericShape;


/**
.Class.Shape
..cat:Index
..summary:Stores hash value and shape for an ungapped or gapped q-gram.
..signature:Shape<TValue, TSpec>
..param.TValue:The @Metafunction.Value@ type of the string the shape is applied to (e.g. $Dna$).
..param.TSpec:The specializing type.
...default:@Spec.SimpleShape@, for ungapped q-grams.
..remarks:The @Metafunction.ValueSize@ of Shape is the ValueSize of TValue which is the alphabet size.
..remarks:To get the span or the weight of a shape call @Function.length@ or @Function.weight@.
..example
...text:The following code shows how one can use a gapped shape to search for the pattern "ACxA" in a reference. First
we assign a form to the shape and then compute the corresponding hash value. The hash value of a string and a Shape 
object is unique, such that one can retrieve the string from a shape if the hash value is known.
...file:demos/index/shape.cpp
...output:The hash is: 4
Hit at position: 0
Hit at position: 14
Hit at position: 17
.Memfunc.Shape#Shape:
..class:Class.Shape
..summary:Constructor
..signature:Shape<TValue, TSpec> ()
..signature:Shape<TValue, TSpec> (shape)
..param.shape:Other Shape object. (copy constructor)
..include:seqan/index.h
*/
/*!
 * @class Shape
 * 
 * @brief Stores hash value and shape for an ungapped or gapped q-gram.
 * 
 * @signature Shape<TValue, TSpec>
 * 
 * @tparam TSpec The specializing type. Default: @link SimpleShape @endlink, for
 *               ungapped q-grams.
 * @tparam TValue The @link Value @endlink type of the string the shape is
 *                applied to (e.g. <tt>Dna</tt>).
 * 
 * @section Remarks
 * 
 * The @link ValueSize @endlink of Shape is the ValueSize of TValue which is the
 * alphabet size.
 * 
 * To get the span or the weight of a shape call @link length @endlink or @link
 * weight @endlink.
 */
	template <typename TValue = Dna, typename TSpec = SimpleShape>
	class Shape;

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.Shape
///.Metafunction.Value.class:Class.Shape
	template <typename TValue, typename TSpec>
	struct Value<Shape<TValue,TSpec> >
	{
		typedef __uint64 Type;
	};

///.Metafunction.Size.param.T.type:Class.Shape
///.Metafunction.Size.class:Class.Shape
	template <typename TValue, typename TSpec>
	struct Size<Shape<TValue,TSpec> >
	{
		typedef unsigned long Type;
	};

///.Metafunction.LENGTH.param.T.type:Class.Shape
///.Metafunction.LENGTH.class:Class.Shape
    template <typename TValue, unsigned q>
	struct LENGTH< Shape<TValue, UngappedShape<q> > >
	{
		enum { VALUE = q };
	};

///.Metafunction.WEIGHT.param.T.type:Class.Shape
///.Metafunction.WEIGHT.class:Class.Shape
    template <typename TValue, unsigned q>
	struct WEIGHT< Shape<TValue, UngappedShape<q> > >
	{
		enum { VALUE = q };
	};

///.Metafunction.ValueSize.param.T.type:Class.Shape
///.Metafunction.ValueSize.class:Class.Shape
	template <typename TValue, typename TSpec>
	struct ValueSize< Shape<TValue, TSpec> > 
	{
		typedef typename Value<Shape<TValue, TSpec> >::Type THashValue;
		static const THashValue VALUE = Power<
						ValueSize<TValue>::VALUE, 
						WEIGHT< Shape<TValue, TSpec> >::VALUE >::VALUE;
	};

///.Metafunction.Host.param.T.type:Class.Shape
///.Metafunction.Host.class:Class.Shape

	template <typename TValue, typename TSpec>
	struct Host<Shape<TValue,TSpec> >
	{
		typedef TValue Type;
	};


//////////////////////////////////////////////////////////////////////////////

/**
.Spec.SimpleShape:
..cat:Index
..summary:A variable length ungapped shape (also called q-gram or k-mer).
..general:Class.Shape
..signature:Shape<TValue, SimpleShape>
..param.TValue:The @Metafunction.Value@ type of the string the shape is applied to (e.g. $Dna$).
..remarks:A SimpleShape must be resized first to a valid length. To do so, call @Function.resize@.
..see:Spec.UngappedShape
..include:seqan/index.h
*/
/*!
 * @class SimpleShape
 * 
 * @extends Shape
 * 
 * @headerfile seqan/index.h
 * 
 * @brief A variable length ungapped shape (also called q-gram or k-mer).
 * 
 * @signature Shape<TValue, SimpleShape>
 * 
 * @tparam TValue The @link Value @endlink type of the string the shape is
 *                applied to (e.g. <tt>Dna</tt>).
 * 
 * @section Remarks
 * 
 * A SimpleShape must be resized first to a valid length. To do so, call @link
 * resize @endlink.
 * 
 * @see UngappedShape
 */

	//////////////////////////////////////////////////////////////////////////////
	// ungapped shape with variable length
	//////////////////////////////////////////////////////////////////////////////

	template <typename TValue>
	class Shape<TValue, SimpleShape>
	{
	public:
//____________________________________________________________________________

		unsigned					span;
		typename Value<Shape>::Type	hValue;
		typename Value<Shape>::Type	XValue;
		typename Value<Shape>::Type	leftFactor;
		typename Value<Shape>::Type	leftFactor2;
		TValue						leftChar;
//____________________________________________________________________________
		
/**
.Memfunc.SimpleShape#Shape:
..class:Spec.SimpleShape
..summary:Constructor
..signature:Shape<TValue, SimpleShape> ()
..signature:Shape<TValue, SimpleShape> (shape)
..signature:Shape<TValue, SimpleShape> (q)
..param.shape:Other Shape object. (copy constructor)
..param.q:Length of the ungapped q-gram.
*/
/*!
 * @fn SimpleShape::Shape
 * 
 * @brief Constructor
 * 
 * @signature Shape<TValue, SimpleShape> ()
 * @signature Shape<TValue, SimpleShape> (shape)
 * @signature Shape<TValue, SimpleShape> (q)
 * 
 * @param q Length of the ungapped q-gram.
 * @param shape Other Shape object. (copy constructor)
 */ 
		Shape():
			span(0),
			hValue(0),
			XValue(0),
			leftFactor(0),
			leftFactor2(0),
            leftChar(0) {}
		
		Shape(unsigned _span):
			hValue(0),
			XValue(0),
			leftFactor(0),
			leftFactor2(0),
			leftChar(0)
		{
			resize(*this, _span);
		}

		template <unsigned q>
		Shape(Shape<TValue, UngappedShape<q> > const &other)
		{
			*this = other;
		}	

//____________________________________________________________________________

		template <unsigned q>
		inline Shape &
		operator=(Shape<TValue, UngappedShape<q> > const &other)
		{
			span = other.span;
			hValue = other.hValue;
			XValue = other.XValue;
			leftFactor = other.leftFactor;
			leftFactor2 = other.leftFactor2;
			leftChar = other.leftChar;
			return *this;
		}
	};

	//////////////////////////////////////////////////////////////////////////////
	// ungapped shape with fixed length q
	//////////////////////////////////////////////////////////////////////////////

/**
.Spec.UngappedShape:
..cat:Index
..summary:A fixed length ungapped shape (also called q-gram or k-mer).
..general:Class.Shape
..signature:Shape<TValue, UngappedShape<q> >
..param.TValue:The @Metafunction.Value@ type of the sequence the shape is applied to (e.g. $Dna$).
..param.q:The length of the shape.
..include:seqan/index.h
*/
/*!
 * @class UngappedShape
 * 
 * @extends Shape
 * 
 * @headerfile seqan/index.h
 * 
 * @brief A fixed length ungapped shape (also called q-gram or k-mer).
 * 
 * @signature Shape<TValue, UngappedShape<q> >
 * 
 * @tparam q The length of the shape.
 * @tparam TValue The @link Value @endlink type of the sequence the shape is
 *                applied to (e.g. <tt>Dna</tt>).
 * 
 * @see SimpleShape
 */

	template <typename TValue, unsigned q>
	class Shape<TValue, UngappedShape<q> >
	{
	public:
		typedef typename Value<Shape>::Type THashValue;
//____________________________________________________________________________

		static const unsigned span = q;
		static const THashValue leftFactor = Power<ValueSize<TValue>::VALUE, q - 1>::VALUE;
		static const THashValue leftFactor2 = (Power<ValueSize<TValue>::VALUE, q>::VALUE - 1) / (ValueSize<TValue>::VALUE - 1);
		// Sigma^(q-1) + Sigma^(q-2) + ... + Sigma + 1

		THashValue	hValue;		// current hash value
		THashValue	XValue;		// Sum_{i=0..q-1} (x_i + 1)
		TValue		leftChar;	// leftmost character
//____________________________________________________________________________
		Shape():
			hValue(0),
			XValue(0),
            leftChar(0) {}
	};



//////////////////////////////////////////////////////////////////////////////

///.Function.value.param.object.type:Class.Shape
///.Function.value.class:Class.Shape
	template <typename TValue, typename TSpec>
	inline typename Value< Shape<TValue, TSpec> >::Type
	value(Shape<TValue, TSpec> &me)
	{
		return me.hValue;
	}

	template <typename TValue, typename TSpec>
	inline typename Value< Shape<TValue, TSpec> >::Type
	value(Shape<TValue, TSpec> const &me)
	{
		return me.hValue;
	}


//____________________________________________________________________________

///.Function.length.param.object.type:Class.Shape
///.Function.length.class:Class.Shape
	template <typename TValue, typename TSpec>
	inline typename Size< Shape<TValue, TSpec> >::Type
	length(Shape<TValue, TSpec> const &me)
	{
	SEQAN_CHECKPOINT
		return me.span;
	}

//____________________________________________________________________________

/**.Function.weight:
..cat:Index
..summary:Number of relevant positions in a shape.
..signature:weight(shape)
..class:Class.Shape
..param.shape:Shape object for which the number of relevant positions is determined.
...type:Class.Shape
..returns:Number of relevant positions.
..remarks.text:For ungapped shapes the return value is the result of the @Function.length@ function.
For gapped shapes this is the number of '1's.
*/
/*!
 * @fn Shape#weight
 * 
 * @brief Number of relevant positions in a shape.
 * 
 * @signature weight(shape)
 * 
 * @param shape Shape object for which the number of relevant positions is
 *              determined. Types: @link Shape @endlink
 * 
 * @return TReturn Number of relevant positions.
 * 
 * @section Remarks
 * 
 * For ungapped shapes the return value is the result of the @link length
 * @endlink function. For gapped shapes this is the number of '1's.
 */
	template <typename TValue, typename TSpec>
	inline typename Size< Shape<TValue, TSpec> >::Type
	weight(Shape<TValue, TSpec> const &me)
	{
	SEQAN_CHECKPOINT
		return length(me);
	}

//____________________________________________________________________________

///.Function.resize.param.object.type:Spec.SimpleShape
///.Function.resize.class:Spec.SimpleShape
	template <typename TValue, typename TSize>
	inline typename Size< Shape<TValue, SimpleShape> >::Type
	resize(Shape<TValue, SimpleShape> & me, TSize new_length)
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, SimpleShape> >::Type	THValue;
		me.leftFactor = _intPow((THValue)ValueSize<TValue>::VALUE, new_length - 1);
		me.leftFactor2 = (_intPow((THValue)ValueSize<TValue>::VALUE, new_length) - 1) / (ValueSize<TValue>::VALUE - 1);
		return me.span = new_length;
	}

//____________________________________________________________________________

/**
.Function.hash:
..cat:Index
..summary:Computes a (lower) hash value for a shape applied to a sequence.
..description:The hash value (a.k.a. code) of a q-gram is the lexicographical rank of this q-gram in the set of all possible q-grams.
For example, the hash value of the @Spec.Dna@ 3-gram AAG is 2 as there are only two 3-grams (AAA and AAC) having a smaller lexicographical rank.
If @Function.hash@ is called with a gapped shape, the q-gram is the text subsequence of no-gap shape positions relative to the text iterator, e.g. a shape 1101 at the beginning of text ACGT corresponds to the 3-gram ACT.

..signature:hash(shape, it)
..signature:hash(shape, it, charsLeft)
..class:Class.Shape
..param.shape:Shape to be used for hashing.
...type:Class.Shape
..param.it:Sequence iterator pointing to the first character of the shape.
..param.charsLeft:The distance of $it$ to the string end. 
If $charsLeft$ is smaller than the shape's span, the hash value corresponds to the lexicographically smallest shape beginning with $charsLeft$ characters. The hash value of such a truncated shape corresponds to the shape applied to a text padded with the smallest alphabet characters.
..returns:Hash value of the shape.
..see:Function.hashNext
..see:Function.hashUpper
..see:Function.hash2
..example:
...text:Code example that computes hash values of 4-grams with different shapes starting at the beginning of a text.
...file:demos/index/shape_hash.cpp
...text:The resulting hexadecimal hash values of the three 4-mers GATT, GATC and GATA are:
...output:
0x8f
0x8d
0x8c
..include:seqan/index.h
*/
/*!
 * @fn Shape#hash
 * 
 * @brief Computes a (lower) hash value for a shape applied to a sequence.
 * 
 * @signature hash(shape, it)
 * @signature hash(shape, it, charsLeft)
 * 
 * @param charsLeft The distance of <tt>it</tt> to the string end. If
 *                  <tt>charsLeft</tt> is smaller than the shape's span, the
 *                  hash value corresponds to the smallest shape beginning with
 *                  <tt>charsLeft</tt> characters.
 * @param shape Shape to be used for hashing. Types: @link Shape @endlink
 * @param it Sequence iterator pointing to the first character of the shape.
 * 
 * @return TReturn Hash value of the shape.
 * 
 * @see hashNext
 * @see hashUpper
 * @see hash2
 * @see unhash
 */
	template <typename TValue, typename TIter>
	inline typename Value< Shape<TValue, SimpleShape> >::Type
	hash(Shape<TValue, SimpleShape> &me, TIter it)
	{
		//typedef typename Value< Shape<TValue, SimpleShape> >::Type	THValue;
		typedef typename Size< Shape<TValue, SimpleShape> >::Type	TSize;

		SEQAN_ASSERT_GT((unsigned)me.span, 0u);

		me.hValue = ordValue(me.leftChar = *it);
		for(TSize i = 1; i < me.span; ++i) {
			++it;
			me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
		}
		return me.hValue;
	}

/**
.Function.hashInit:
..cat:Index
..summary:Preprocessing step of a pure @Function.hashNext@ loop.
..description:Overlapping q-grams can efficiently be hashed by calling @Function.hash@ on the first text position and @Function.hashNext@ on succeeding, adjacent positions. One drawback of this scenario is that for-loops cannot start with the first position directly and become more complicated. As a remedy, @Function.hashInit@ was introduced which initializes the @Class.Shape@ to be used with @Function.hashNext@ on the first position directly.
..signature:hashInit(shape, it)
..class:Class.Shape
..param.shape:Shape to be used for hashing.
...type:Class.Shape
..param.it:Sequence iterator pointing to the first text character of the shape.
..see:Function.hashNext
..example:
...text:Two hash loop examples. The first loop uses @Function.hash@/@Function.hashNext@ while the second use @Function.hashInit@/@Function.hashNext@ and can process all hashes within the loop.
...file:demos/index/shape_hash_init.cpp
...text:The two loops produce the same hash values:
...output:
0	0	1	4	17	4	18	11	47	63	62	56	
0	0	1	4	17	4	18	11	47	63	62	56
..include:seqan/index.h
*/
	template <typename TValue, typename TIter>
	inline void
	hashInit(Shape<TValue, SimpleShape> &me, TIter it)
	{
		//typedef typename Value< Shape<TValue, SimpleShape> >::Type	THValue;
		typedef typename Size< Shape<TValue, SimpleShape> >::Type	TSize;

		SEQAN_ASSERT_GT((unsigned)me.span, 0u);

        me.leftChar = 0;
		me.hValue = ordValue(*it);
		for(TSize i = 2; i < me.span; ++i) {
			++it;
			me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
		}
	}

//____________________________________________________________________________
// fixed ungapped shapes

	// loop unrolling ...
	template <typename THValue, typename TValue, typename TIter>
	inline THValue
	_hashFixedShape(THValue hash, TIter &, TValue const, UngappedShape<1> const) {
		return hash;
	}

	template <typename THValue, typename TValue, typename TIter, unsigned q>
	inline THValue
	_hashFixedShape(THValue hash, TIter &it, TValue const, UngappedShape<q> const) {
		++it;
		return _hashFixedShape(
			hash * ValueSize<TValue>::VALUE + ordValue((TValue)*it),
			it, TValue(), UngappedShape<q - 1>());
	}

	// ... for fixed ungapped shapes
	template <typename TValue, unsigned q, typename TIter>
	inline typename Value< Shape<TValue, UngappedShape<q> > >::Type
	hash(Shape<TValue, UngappedShape<q> > &me, TIter it)
	{
		//typedef typename Value< Shape<TValue, UngappedShape<q> > >::Type	THValue;
		//typedef typename Size< Shape<TValue, UngappedShape<q> > >::Type     TSize;

		me.hValue = ordValue(me.leftChar = *it);
		return me.hValue = _hashFixedShape(me.hValue, it, TValue(), UngappedShape<q>());
	}

	template <typename TValue, unsigned q, typename TIter>
	inline typename Value< Shape<TValue, UngappedShape<q> > >::Type
	hashInit(Shape<TValue, UngappedShape<q> > &me, TIter it)
	{
		//typedef typename Value< Shape<TValue, UngappedShape<q> > >::Type	THValue;
		//typedef typename Size< Shape<TValue, UngappedShape<q> > >::Type	TSize;

        me.leftChar = 0;
		me.hValue = ordValue(*it);

        if (q > 1)
            me.hValue = _hashFixedShape(me.hValue, it, TValue(), UngappedShape<q-1>());

		return me.hValue;
    }

	template <typename TValue, typename TSpec, typename TIter, typename TSize>
	inline typename Value< Shape<TValue, TSpec> >::Type
	hash(Shape<TValue, TSpec> &me, TIter it, TSize charsLeft)
	{
		//typedef typename Value< Shape<TValue, TSpec> >::Type	THValue;

		SEQAN_ASSERT_GT((unsigned)me.span, 0u);

		TSize iEnd = me.span;
		if (iEnd > charsLeft) iEnd = charsLeft;

		TSize i = 0;
		if (iEnd > 0) {
			me.hValue = ordValue(me.leftChar = *it);
			for(i = 1; i < iEnd; ++i) {
				++it;
				me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
			}
		} else
			return me.hValue = 0;

		// fill shape with zeros
		for(; i < (TSize)me.span; ++i)
			me.hValue *= ValueSize<TValue>::VALUE;
		return me.hValue;
	}

//____________________________________________________________________________
// Tuple -> fixed ungapped shapes

	template <typename THValue, typename TValue, typename TTValue, unsigned SIZE, typename TPack>
	inline THValue
	_hashTuple2FixedShape(
		THValue const, 
		Tuple<TTValue, SIZE, TPack> const &tuple,
		TValue const,
		UngappedShape<1> const) 
	{
		return ordValue(tuple[0]);
	}

	template <typename THValue, typename TValue, typename TTValue, unsigned SIZE, typename TPack, unsigned q>
	inline THValue
	_hashTuple2FixedShape(
		THValue const, 
		Tuple<TTValue, SIZE, TPack> const &tuple,
		TValue const,
		UngappedShape<q> const) 
	{
		return _hashTuple2FixedShape(THValue(), tuple, TValue(), UngappedShape<q - 1>()) 
			* ValueSize<TValue>::VALUE + ordValue(tuple[q-1]);
	}

	// ... for fixed ungapped shapes
	template <
		typename TValue,
		typename TTValue, 
		unsigned SIZE, 
		unsigned q>
	typename Value< Shape<TValue, UngappedShape<q> > >::Type
	hash(
		Shape<TValue, UngappedShape<q> > &me, 
		Tuple<TTValue, SIZE, BitPacked<> > /*const &*/tuple)
	{
	SEQAN_CHECKPOINT
		if (ValueSize<TValue>::VALUE == (1 << BitsPerValue<TTValue>::VALUE))
			if (q == SIZE)
				return tuple.i;
			else
				return tuple >> (q - SIZE);
		else
			return me.hValue = _hashTuple2FixedShape(me.hValue, tuple, TValue(), UngappedShape<q>());
	}

	template <
		typename TValue,
		typename TTValue, 
		unsigned SIZE, 
		typename TPack, 
		unsigned q>
	typename Value< Shape<TValue, UngappedShape<q> > >::Type
	hash(
		Shape<TValue, UngappedShape<q> > &me, 
		Tuple<TTValue, SIZE, TPack> /*const &*/tuple)
	{
	SEQAN_CHECKPOINT
		return me.hValue = _hashTuple2FixedShape(me.hValue, tuple, TValue(), UngappedShape<q>());
	}

//____________________________________________________________________________

/**.Function.hashUpper:
..cat:Index
..summary:Computes an upper hash value for a shape applied to a sequence.
..signature:hashUpper(shape, it, charsLeft)
..class:Class.Shape
..param.shape:Shape to be used for hashing.
...type:Class.Shape
..param.it:Sequence iterator pointing to the first character of the shape.
..param.charsLeft:The distance of $it$ to the string end. 
..returns:Upper hash value of the shape.
The hash value corresponds to the maximal @Function.hash@ value of a shape beginning with $min(charsLeft,length(shape))$ characters + 1.
..remarks:This function in conjunction with @Function.hash@ is useful to search a q-gram index for p-grams with p<q.
*/
/*!
 * @fn Shape#hashUpper
 * 
 * @brief Computes an upper hash value for a shape applied to a sequence.
 * 
 * @signature hashUpper(shape, it, charsLeft)
 * 
 * @param charsLeft The distance of <tt>it</tt> to the string end.
 * @param shape Shape to be used for hashing. Types: @link Shape @endlink
 * @param it Sequence iterator pointing to the first character of the shape.
 * 
 * @return TReturn Upper hash value of the shape. The hash value corresponds to
 *                 the maximal @link hash @endlink value of a shape beginning
 *                 with <tt>min(charsLeft,length(shape))</tt> characters + 1.
 * 
 * @section Remarks
 * 
 * This function in conjunction with @link hash @endlink is useful to search a
 * q-gram index for p-grams with p<q.
 * 
 * @see hash
 */

	template <typename TValue, typename TSpec, typename TIter, typename TSize>
	inline typename Value< Shape<TValue, TSpec> >::Type
	hashUpper(Shape<TValue, TSpec> &me, TIter it, TSize charsLeft)
	{
		//typedef typename Value< Shape<TValue, TSpec> >::Type	THValue;

		SEQAN_ASSERT_GT((unsigned)me.span, 0u);

		TSize iEnd = me.span;
		if (iEnd > charsLeft) iEnd = charsLeft;

		TSize i = 0;
		if (iEnd > 0) {
			me.hValue = ordValue(me.leftChar = *it);
			for(i = 1; i < iEnd; ++i) {
				++it;
				me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
			}
			++me.hValue;
		} else
			me.hValue = 1;

		// fill shape with zeros
		for(; i < (TSize)me.span; ++i)
			me.hValue *= ValueSize<TValue>::VALUE;
		return me.hValue;
	}

//____________________________________________________________________________

/**
.Function.hashNext:
..cat:Index
..summary:Computes the hash value for the adjacent shape.
..description:Overlapping q-grams can efficiently be hashed by calling @Function.hash@ on the first text position and @Function.hashNext@ on succeeding, adjacent positions. Alternatively, @Function.hashInit@ can be used to replace the first @Function.hash@ call by another @Function.hashNext@ call and hence ease loops iterating all overlapping q-grams.
..signature:hashNext(shape, it)
..class:Class.Shape
..param.shape:Shape to be used for hashing.
...type:Class.Shape
..param.it:Sequence iterator pointing to the first character of the adjacent shape.
..returns:Hash value of the q-gram.
..remarks:@Function.hash@ has to be called before.
..example:Hash loop example that outputs the hash values of all overlapping 3-grams using @Function.hash@ and @Function.hashNext@.
...file:demos/index/shape_hash_next.cpp
...text:The overlapping hash values are:
...output:
0	0	1	4	17	4	18	11	47	63	62	56	
..include:seqan/index.h
*/
/*!
 * @fn Shape#hashNext
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Computes the hash value for the adjacent shape.
 * 
 * @signature hashNext(shape, it)
 * 
 * @param shape Shape to be used for hashing. Types: @link Shape @endlink
 * @param it Sequence iterator pointing to the first character of the adjacent
 *           shape.
 * 
 * @return TReturn Hash value of the q-gram.
 * 
 * @section Remarks
 * 
 * @link hash @endlink has to be called before.
 * 
 * @see hash
 */
	template <typename TValue, typename TSpec, typename TIter>
	inline typename Value< Shape<TValue, TSpec> >::Type
	hashNext(Shape<TValue, TSpec> &me, TIter const &it)
	{
	SEQAN_CHECKPOINT
		// remove first, shift left, and add next character
		typedef typename Value< Shape<TValue, TSpec> >::Type	THValue;
		typedef typename Size< Shape<TValue, TSpec> >::Type		TSize;

		SEQAN_ASSERT_GT((unsigned)me.span, 0u);

		me.hValue = 
			(me.hValue - ordValue(me.leftChar) * (THValue)me.leftFactor) * ValueSize<TValue>::VALUE
			+ ordValue((TValue)*(it + ((TSize)me.span - 1)));
		me.leftChar = *it;
		return me.hValue;
	}

//____________________________________________________________________________

/**.Function.hash2:
..cat:Index
..summary:Computes an unique hash value of a shape applied to a sequence, even if the sequence is shorter than the shape span
..signature:hash2(shape, it, charsLeft)
..class:Class.Shape
..param.shape:Shape to be used for hashing.
...type:Class.Shape
..param.it:Sequence iterator pointing to the first character of the shape.
..param.charsLeft:The distance of $it$ to the string end. 
..returns:Hash value of the shape.
..see:Function.hash2Next
..see:Function.hash2Upper
*/
/*!
 * @fn Shape#hash2
 * 
 * @brief Computes an unique hash value of a shape applied to a sequence, even
 *        if the sequence is shorter than the shape span
 * 
 * @signature hash2(shape, it, charsLeft)
 * 
 * @param charsLeft The distance of <tt>it</tt> to the string end.
 * @param shape Shape to be used for hashing. Types: @link Shape @endlink
 * @param it Sequence iterator pointing to the first character of the shape.
 * 
 * @return TReturn Hash value of the shape.
 * 
 * @see hash2Next
 * @see hash2Upper
 * @see unhash
 * @see hash
 */

	template <typename TValue, typename TSpec, typename TIter, typename TSize>
	inline typename Value< Shape<TValue, TSpec> >::Type
	hash2(Shape<TValue, TSpec> &me, TIter it, TSize charsLeft)
	{
		//typedef typename Value< Shape<TValue, TSpec> >::Type	THValue;

		SEQAN_ASSERT_GT((unsigned)me.span, 0u);

		TSize iEnd = me.span;
		if (iEnd > charsLeft) iEnd = charsLeft;

		TSize i = 0;
		if (iEnd > 0) {
			me.hValue = me.XValue = ordValue(me.leftChar = *it);
			for(i = 1; i < iEnd; ++i) {
				++it;
				// update sum of x_i
				me.XValue += ordValue((TValue)*it);
				// shift hash
				me.hValue = me.hValue * ValueSize<TValue>::VALUE + me.XValue;
			}
		} else
			return me.hValue = me.XValue = 0;

		// fill shape with zeros
		for(; i < (TSize)me.span; ++i)
			me.hValue = me.hValue * ValueSize<TValue>::VALUE + me.XValue;
		return me.hValue += iEnd;
	}

/**.Function.hash2Upper:
..cat:Index
..summary:Computes an upper unique hash value of a shape applied to a sequence, even if the sequence is shorter than the shape span.
..signature:hash2Upper(shape, it, charsLeft)
..class:Class.Shape
..param.shape:Shape to be used for hashing.
...type:Class.Shape
..param.it:Sequence iterator pointing to the first character of the shape.
..param.charsLeft:The distance of $it$ to the string end. 
..returns:Upper hash value of the shape.
The hash value corresponds to the maximal @Function.hash2@ value of a shape beginning with the $min(charsLeft,length(shape))$ characters + 1
..remarks:This function in conjunction with @Function.hash2@ is useful to search a q-gram index for p-grams with p<q.
*/
/*!
 * @fn Shape#hash2Upper
 * 
 * @brief Computes an upper unique hash value of a shape applied to a sequence,
 *        even if the sequence is shorter than the shape span.
 * 
 * @signature hash2Upper(shape, it, charsLeft)
 * 
 * @param charsLeft The distance of <tt>it</tt> to the string end.
 * @param shape Shape to be used for hashing. Types: @link Shape @endlink
 * @param it Sequence iterator pointing to the first character of the shape.
 * 
 * @return TReturn Upper hash value of the shape. The hash value corresponds to
 *                 the maximal @link hash2 @endlink value of a shape beginning
 *                 with the <tt>min(charsLeft,length(shape))</tt> characters + 1
 * 
 * @section Remarks
 * 
 * This function in conjunction with @link hash2 @endlink is useful to search a
 * q-gram index for p-grams with p<q.
 * 
 * @see hash2
 */
	template <typename TValue, typename TSpec, typename TIter, typename TSize>
	inline typename Value< Shape<TValue, TSpec> >::Type
	hash2Upper(Shape<TValue, TSpec> &me, TIter it, TSize charsLeft)
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, TSpec> >::Type	THValue;

		SEQAN_ASSERT_GT((unsigned)me.span, 0u);

		TSize iEnd = me.span;
		if (iEnd > charsLeft) iEnd = charsLeft;

		THValue hValue, XValue;
		TSize i = 0;
		if (iEnd > 0) {
			hValue = XValue = ordValue((TValue)*it);
			for(i = 1; i < iEnd; ++i) {
				++it;
				// update sum of x_i
				XValue += ordValue((TValue)*it);
				// shift hash
				hValue = hValue * ValueSize<TValue>::VALUE + XValue;
			}
		} else
			hValue = XValue = 0;

		if (charsLeft <= me.span) {
			++XValue;
			++hValue;
		}

		// fill shape with zeros
		for(; i < (TSize)me.span; ++i)
			hValue = hValue * ValueSize<TValue>::VALUE + XValue;
		return hValue += iEnd;
	}

//____________________________________________________________________________

/**
.Function.hash2Next:
..cat:Index
..summary:Computes a unique hash value for the adjacent shape, even if it is shorter than q.
..signature:hash2Next(shape, it)
..class:Class.Shape
..param.shape:Shape to be used for hashing the q-gram.
...type:Class.Shape
..param.it:Sequence iterator pointing to the first character of the adjacent shape.
..returns:Hash value of the shape.
..remarks:@Function.hash@ has to be called before with $shape$ on the left adjacent q-gram.
..include:seqan/index.h
*/
/*!
 * @fn Shape#hash2Next
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Computes a unique hash value for the adjacent shape, even if it is
 *        shorter than q.
 * 
 * @signature hash2Next(shape, it)
 * 
 * @param shape Shape to be used for hashing the q-gram. Types: @link Shape @endlink
 * @param it Sequence iterator pointing to the first character of the adjacent
 *           shape.
 * 
 * @return TReturn Hash value of the shape.
 * 
 * @section Remarks
 * 
 * @link hash @endlink has to be called before with <tt>shape</tt> on the left
 * adjacent q-gram.
 * 
 * @see hash2
 */

	template <typename TValue, typename TSpec, typename TIter, typename TSize>
	inline typename Value< Shape<TValue, TSpec> >::Type
	hash2Next(Shape<TValue, TSpec> &me, TIter &it, TSize charsLeft)
	{
	SEQAN_CHECKPOINT
		// remove first, shift left, and add next character
		typedef typename Value< Shape<TValue, TSpec> >::Type	THValue;

		SEQAN_ASSERT_GT((unsigned)me.span, 0u);

		if (charsLeft >= me.span) {
			// update sum of x_i
			me.XValue = me.XValue + ordValue((TValue)*(it + me.span - 1)) - ordValue(me.leftChar);
			// shift hash
			me.hValue = (me.hValue - ordValue(me.leftChar) * (THValue)me.leftFactor2) * ValueSize<TValue>::VALUE + me.XValue
						- me.span * (ValueSize<TValue>::VALUE - 1);
		} else {
			// update sum of x_i
			me.XValue -= ordValue(me.leftChar);
			// shift hash
			me.hValue = (me.hValue - ordValue(me.leftChar) * (THValue)me.leftFactor2) * ValueSize<TValue>::VALUE + me.XValue
				        - charsLeft * (ValueSize<TValue>::VALUE - 1) - ValueSize<TValue>::VALUE;
		}

		me.leftChar = *it;
		return me.hValue;
	}

/**
.Function.unhash:
..cat:Index
..summary:Inverse of the @Function.hash@ function; for ungapped shapes.
..signature:unhash(result, hash, q)
..class:Class.Shape
..param.result:String to write the result to.
...type:Class.String
..param.hash:The hash value previously computed with @Function.hash@.
...type:nolink:A number.
..param.q:The $q$-gram length.
...type:nolink:$unsigned$
..remarks:
..see:Function.hash
..see:Function.hash2
..include:seqan/index.h
*/
/*!
 * @fn Shape#unhash
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Inverse of the @link hash @endlink function; for ungapped shapes.
 * 
 * @signature unhash(result, hash, q)
 * 
 * @param q The <tt>q</tt>-gram length. Types: nolink:<tt>unsigned</tt>
 * @param hash The hash value previously computed with @link hash @endlink.
 *             Types:
 * @param result String to write the result to. Types: @link String @endlink
 * 
 * @section Remarks
 * 
 * @see hash
 * @see hash2
 */
	template <typename TString, typename THash>
	inline void unhash(TString &result, THash hash, unsigned q)
	{
	SEQAN_CHECKPOINT
		typedef typename Value<TString>::Type	TValue;

		resize(result, q);
		for (unsigned i = q; i > 0; ) 
		{
			result[--i] = (TValue)(hash % ValueSize<TValue>::VALUE);
			hash /= ValueSize<TValue>::VALUE;
		}
	}

//____________________________________________________________________________

/**.Function.stringToShape:
..cat:Index
..summary:Takes a shape given as a string of '1' (relevant position) and '0' 
(irrelevant position) and converts it into a Shape object.
..signature:stringToShape(shape, bitmap)
..class:Class.Shape
..param.shape:Shape object that is manipulated.
...type:Spec.SimpleShape
..param.bitmap:A character string of '1' and '0' representing relevant and irrelevant positions (blanks) respectively.
...remarks:This string must begin with a '1'. Trailing '0's are ignored.
...remarks:If $shape$ is a @Spec.SimpleShape@ at most one contiguous sequences of $1$s is allowed.
...type:Class.String
..see:Function.shapeToString
..see:Function.reverse
*/
/*!
 * @fn Shape#stringToShape
 * 
 * @brief Takes a shape given as a string of '1' (relevant position) and '0'
 *        (irrelevant position) and converts it into a Shape object.
 * 
 * @signature stringToShape(shape, bitmap)
 * 
 * @param shape Shape object that is manipulated. Types: SimpleShape,
 *              GenericShape, OneGappedShape
 * @param bitmap A character string of '1' and '0' representing relevant and
 *               irrelevant positions (blanks) respectively.This string must
 *               begin with a '1'. Trailing '0's are ignored.If <tt>shape</tt>
 *               is a @link SimpleShape @endlink at most one contiguous
 *               sequences of <tt>1</tt>s is allowed.If <tt>shape</tt> is a
 *               @link OneGappedShape @endlink at most two contiguous sequences
 *               of '1's are allowed. Types: @link String @endlink
 * 
 * @see shapeToString
 * @see reverse
 */
	template <typename TValue, typename TShapeString>
	inline bool
	stringToShape(
		Shape<TValue, SimpleShape> &me, 
		TShapeString const &bitmap)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TShapeString const>::Type		TIter;
		typedef typename Size<TShapeString const>::Type			TSize;

		TIter it = begin(bitmap, Standard());
		TIter itEnd = end(bitmap, Standard());

		TSize ones = 0;
		for(; it != itEnd && *it == '0' ; ++it) ;
		for(; it != itEnd && *it == '1' ; ++it)	++ones;
		for(; it != itEnd && *it == '0' ; ++it) ;

		resize(me, ones);

		return it == itEnd;
	}

//____________________________________________________________________________

/**.Function.shapeToString:
..cat:Index
..class:Class.Shape
..summary:Converts a given shape into a sequence of '1' (relevant position) and '0' 
(irrelevant position).
..signature:shapeToString(bitmap, shape)
..class:Class.Shape
..param.bitmap:The resulting sequence object.
...type:Class.String
..param.shape:Shape object.
...type:Class.Shape
..see:Function.stringToShape
*/
/*!
 * @fn Shape#shapeToString
 * 
 * @brief Converts a given shape into a sequence of '1' (relevant position) and
 *        '0' (irrelevant position).
 * 
 * @signature shapeToString(bitmap, shape)
 * 
 * @param shape Shape object. Types: @link Shape @endlink
 * @param bitmap The resulting sequence object. Types: @link String @endlink
 * 
 * @see stringToShape
 */

	template <typename TShapeString, typename TValue, unsigned q>
	inline void
	shapeToString(
		TShapeString &bitmap,
		Shape<TValue, UngappedShape<q> > const &me)
	{
	SEQAN_CHECKPOINT

		clear(bitmap);
		resize(bitmap, length(me), '1');
	}

//____________________________________________________________________________
	
///.Function.reverse.param.object.type:Spec.SimpleShape
///.Function.reverse.class:Spec.SimpleShape

	template <typename TValue, typename TSpec>
	inline void
	reverse(Shape<TValue, TSpec> &)
	{
	}
	
}	// namespace seqan

#endif
