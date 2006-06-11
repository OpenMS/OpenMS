// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Id: SortedIteratorWrapper.h,v 1.5 2006/06/01 07:39:54 elange Exp $
// $Author: elange $
// $Maintainer:  $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_SORTEDITERATORWRAPPER_H
#define OPENMS_DATASTRUCTURES_SORTEDITERATORWRAPPER_H

#include <vector>
// for std::vector
#include <functional>
// for std::less
#include <algorithm>
// for std::sort

namespace OpenMS
{
  /**
    * @brief An iterator for sorted iteration.
    *
    * This class is a wrapper around a pair of ForwardIterators
    * to iterate over a range of elements in sorted order,
    * according to the compare functor @p Cmp. On construction,
    * pointers to the elements in the speecified range are copied
    * into a vector which gets sorted after that.
    *
    * Objects of this class are Forward Iterators, meaning that
    * it is not possible to iterate backwards. For efficiency
    * reasons, the post-increment operator is NOT provided.
    *
    * Note also that this wrapper class is only useful if you
    * want to iterate over ALL objects in the specified range.
    *
    *
    * @param It The original iterator type.
    * @param Cmp The compare functor. Note that the functor does
    *            NOT introduce an additional level of pointer
    *            indirection when comparing elements.
    */
  template<class It, class Cmp = std::less<typename It::value_type> >
  class SortedIteratorWrapper
  {
  private:
    template<class Functor>
    class Deref
    {
    public:
      typedef typename Functor::Param Param;

      inline Deref() {}
      inline Deref(Functor f) : func(f) {}

      inline bool operator()(Param* a, Param* b)
      {
        return func(*a, *b);
      }

    private:
      Functor func;
    };

  public:
    /**
      * This type is actually not used in this class, but
      * required by some STL algorithms.
      * @sa difference_type, iterator_category
      */
    typedef std::size_t size_type;

    /**
      * This type is actually not used in this class, but
      * required by some STL algorithms.
      * @sa size_type, iterator_category
      */
    typedef std::ptrdiff_t difference_type;

    /**
      * This type is actually not used in this class, but
      * required by some STL algorithms. It hints algorithms
      * what they can expect from this iterator class.
      * @sa size_type, difference_type
      */
    typedef std::forward_iterator_tag iterator_category;

    /**
      * This is the standard interface to the Type template
      * parameter. It is used by most STL algorithms.
      * @sa reference, pointer, ValueType
      */
    typedef typename std::iterator_traits<It>::value_type value_type;

    /**
      * This is the standard interface to the Ref template
      * parameter, denoting the iterator's reference type.
      * It is used by most STL algorithms.
      * @sa value_type, pointer, Reference
      */
    typedef typename std::iterator_traits<It>::reference reference;

    /**
      * This is the standard interface to the Type template
      * parameter, denoting the iterator's pointer type.
      * It is used by most STL algorithms.
      * @sa reference, pointer, Pointer
      */
    typedef typename std::iterator_traits<It>::pointer pointer;

    /**
      * This type is provided, because the OpenMS coding
      * convention states that types start with a capital
      * letter, in contrast to the STL types.
      * @sa value_type
      */
    typedef value_type ValueType;

    /**
      * This type is provided, because the OpenMS coding
      * convention states that types start with a capital
      * letter, in contrast to the STL types.
      * @sa reference
      */
    typedef reference Reference;

    /**
      * This type is provided, because the OpenMS coding
      * convention states that types start with a capital
      * letter, in contrast to the STL types.
      * @sa pointer
      */
    typedef pointer Pointer;

    /**
      * Default constructor. Constructs an iterator pointing
      * to nowhere. This is used by the QuadTree::end() function.
      * You should not use this constructor. Use the QuadTree
      * iterator interface instead.
      */
    inline SortedIteratorWrapper();

    /**
      * Constructor. Creates an iterator which can be used to
      * iterate over the elements in the range [@p start, @p end)
      * in sorted order.
      * @param start An iterator pointing to the first element.
      * @param end An iterator pointing to the last element.
      */
    inline SortedIteratorWrapper(const It& start, const It& end);

    inline SortedIteratorWrapper(const SortedIteratorWrapper& it);

    inline SortedIteratorWrapper& operator=(const SortedIteratorWrapper& it);

    /**
      * Returns a reference to the data pointed to by the iterator.
      * Note that for ConstIterators, this is a const reference.
      * This function's usage is the same as for STL iterators.
      * @return A reference to the data pointed to by the iterator.
      */
    inline Reference operator*() const;

    /**
      * Returns a pointer to the data pointed to by the iterator.
      * Note that for ConstIterators, this is a const pointer.
      * This function's usage is the same as for STL iterators.
      * @return A pointer to the data pointed to by the iterator.
      */
    inline Pointer operator->() const;

    /**
      * Moves the iterator to next item stored in the sorted vector.
      */
    inline SortedIteratorWrapper& operator++();

    /**
      * Checks this iterator and the @p it iterator for equality.
      * Two iterators are considered equal, if both point to the
      * end of their range.
      */
    inline bool operator==(const SortedIteratorWrapper& it);

    /**
      * Checks this iterator and the @p it iterator for inequality.
      * Two iterators are considered inequal, if they are not equal,
      * according to operator==().
      */
    inline bool operator!=(const SortedIteratorWrapper& it);

  private:
    std::vector<Pointer> sorted_;
    typename std::vector<Pointer>::iterator current_;
  };
}

// OpenMS::SortedIteratorWrapper
template<class It, class Cmp>
inline OpenMS::SortedIteratorWrapper<It, Cmp>::SortedIteratorWrapper()
  : current_(sorted_.end())
{

}

template<class It, class Cmp>
inline OpenMS::SortedIteratorWrapper<It, Cmp>::SortedIteratorWrapper(const It& start, const It& end)
{
  // write pointers to all objects into the vector
  for (It i = start; i != end; ++i)
    sorted_.push_back(&*i);
  // sort according to functor Cmp
  std::sort(sorted_.begin(), sorted_.end(), Deref<Cmp>());
  current_ = sorted_.begin();
}

template<class It, class Cmp>
inline OpenMS::SortedIteratorWrapper<It, Cmp>::SortedIteratorWrapper(const OpenMS::SortedIteratorWrapper<It, Cmp>& it)
  : sorted_(it.sorted_), current_(sorted_.begin())
{

}

template<class It, class Cmp>
inline OpenMS::SortedIteratorWrapper<It, Cmp>&
  OpenMS::SortedIteratorWrapper<It, Cmp>::operator=(const OpenMS::SortedIteratorWrapper<It, Cmp>& it)
{
  sorted_ = it.sorted_;
  current_ = sorted_.begin();
}

template<class It, class Cmp>
inline typename OpenMS::SortedIteratorWrapper<It, Cmp>::Reference
  OpenMS::SortedIteratorWrapper<It, Cmp>::operator*() const
{
  return **current_;
}

template<class It, class Cmp>
inline typename OpenMS::SortedIteratorWrapper<It, Cmp>::Pointer
  OpenMS::SortedIteratorWrapper<It, Cmp>::operator->() const
{
  return &(operator*());
}

template<class It, class Cmp>
inline OpenMS::SortedIteratorWrapper<It, Cmp>&
  OpenMS::SortedIteratorWrapper<It, Cmp>::operator++()
{
  current_++;
  return *this;
}

template<class It, class Cmp>
inline bool OpenMS::SortedIteratorWrapper<It, Cmp>::operator==(const SortedIteratorWrapper& it)
{
  // equality is only true, if both iterators point
  // to sorted_.end()
  return current_ == sorted_.end() && it.current_ == it.sorted_.end();
}

template<class It, class Cmp>
inline bool OpenMS::SortedIteratorWrapper<It, Cmp>::operator!=(const SortedIteratorWrapper& it)
{
  return !operator==(it);
}

#endif
