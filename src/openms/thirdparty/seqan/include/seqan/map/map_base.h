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

#ifndef SEQAN_HEADER_MAP_BASE_H
#define SEQAN_HEADER_MAP_BASE_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

// Manual Forwards
template <typename TKey, typename TCargo, typename TCompare, typename TAlloc, typename TKey2>
inline typename Cargo< ::std::map<TKey,TCargo, TCompare, TAlloc> >::Type &
cargo(::std::map<TKey,TCargo, TCompare, TAlloc> & me, TKey2 const & _key);

//////////////////////////////////////////////////////////////////////////////
//insertion tags

template <typename TSpec = Default>
struct Skiplist;


/**
.Class.Map:
..cat:Map
..summary:Set/dictionary container.
..signature:Map<TValue, TSpec >
..param.TValue:Type of values that are stored in the map.
...metafunction:Metafunction.Value
...remarks:Use a @Class.Pair.Pair<Key, Cargo>@ to implement a dictionary mapping from $Key$ to $Cargo$.
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:@Spec.Skiplist@
..include:seqan/map.h
*/
template <typename TElement, typename TSpec = Skiplist<> >
class Map;


//////////////////////////////////////////////////////////////////////////////
// In SeqAn sets and maps store elements as pairs of (key,cargo) 
// the elements of sets without objects are the keys.
//////////////////////////////////////////////////////////////////////////////

/*moved to basic_aggregates.h
template <typename TKey, typename TObject, typename TSpec>
struct Key< Pair<TKey, TObject, TSpec> > 
{
	typedef TKey Type;
};

template <typename TKey, typename TCargo, typename TSpec>
struct Cargo< Pair<TKey, TCargo, TSpec> > 
{
	typedef TCargo Type;
};
*/

//////////////////////////////////////////////////////////////////////////////
// Type for mapValue function that implements [] for map types 

template <typename TMap, typename TCargo>
struct MapValueImpl_
{
	typedef TCargo & Type;
};
template <typename TMap>
struct MapValueImpl_<TMap, Nothing>
{
	typedef bool Type;
};

/**
.Metafunction.MapValue:
..cat:Map
..summary:Type of the map value type.
..signature:MapValue<T>::Type
..class:Class.Map
..param.T:A map type.
...type:Class.Map
..returns.param.Type:The type of the value of T.
..see:Metafunction.Cols
..include:seqan/map.h
 */
template <typename TMap>
struct MapValue :
	MapValueImpl_< TMap, typename Cargo<TMap>::Type >
{
};



template <typename TCargo>
struct ImplMapValue_
{
	template <typename TMap, typename TKey2>
	static inline TCargo &
	mapValue_(TMap & me,
		TKey2 const & _key)
	{
		return cargo(me, _key);
	}
};

template <>
struct ImplMapValue_<Nothing>
{
	template <typename TMap, typename TKey2>
	static inline bool
	mapValue_(TMap & me,
		TKey2 const & _key)
	{
		return hasKey(me, _key);
	}
};

/**
.Function.mapValue:
..cat:Map
..summary:Subscript operator $[ ]$ of maps. 
..signature:MapValue mapValue(map, key)
..class:Class.Map
..param.map:A map.
...type:Class.Map
..param.key:A key.
...metafunction:Metafunction.Key
..returns:If $map$ is a set: The same as @Function.Map#hasKey@.
...text:If $map$ is a dictionary: The same as @Function.Map#value@.
...metafunction:Metafunction.MapValue
..remarks: Usually, @Function.Map#value@ implements the subscript operator $[ ]$, but for maps, 
this operator is implemented in $mapValue$. 
The semantic of this operator depends on the kind of map: If the map has a @Metafunction.Cargo.cargo@, 
than $mapValue(map, key)$ returns the cargo of the (first) value in the map of the given key.
If the map has no @Metafunction.Cargo.cargo@, than the function returns a $true$, if $key$ is in $map$, or $false$ otherwise.
...note:There is no way to create a set of @Class.Pair@, since it is always interpreted as a key/value pair.
If you need a key type that holds two members, define your own key type.
..remarks:You may overload @Metafunction.Key@ and @Metafunction.Cargo@ for your own value type in order to define, what part of your value type is used as key and what as cargo.
..see:Function.Map#value
..see:Function.Map#cargo
..see:Function.Map#hasKey
..include:seqan/map.h
*/

template <typename TMap, typename TKey>
inline typename MapValue<TMap>::Type
mapValue(TMap & me,
		 TKey const & _key)
{
	typedef typename Cargo<TMap>::Type TCargo;
	return ImplMapValue_<TCargo>::mapValue_(me, _key);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TElement>
inline TElement & 
key(TElement & element) 
{
	return element;
}
template <typename TElement>
inline TElement const & 
key(TElement const & element) 
{
	return element;
}

template <typename TKey, typename TObject, typename TSpec>
inline TKey & 
key(Pair<TKey, TObject, TSpec> & element) 
{
	return element.i1;
}
template <typename TKey, typename TObject, typename TSpec>
inline TKey const &
key(Pair<TKey, TObject, TSpec> const & element) 
{
	return element.i1;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TElement, typename TSource>
inline void
setKey(TElement & element,
	   TSource const & source) 
{
	element = source;
}
template <typename TKey, typename TObject, typename TSpec, typename TSource>
inline void 
setKey(Pair<TKey, TObject, TSpec> & element,
	   TSource const & source) 
{
	element.i1 = source;
}

//////////////////////////////////////////////////////////////////////////////
//no default cargo function

template <typename TKey, typename TObject, typename TSpec>
inline TObject & 
cargo(Pair<TKey, TObject, TSpec> & element) 
{
	return element.i2;
}
template <typename TKey, typename TObject, typename TSpec>
inline TObject const &
cargo(Pair<TKey, TObject, TSpec> const & element) 
{
	return element.i2;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TObject, typename TSpec, typename TSource>
inline void 
setCargo(Pair<TKey, TObject, TSpec> & element,
	   TSource const & source) 
{
	element.i2 = source;
}

//////////////////////////////////////////////////////////////////////////////

}

#endif

