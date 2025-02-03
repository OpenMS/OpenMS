// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>

#include <algorithm>
#include <typeinfo>
#include <vector>

namespace OpenMS
{

  /**
    @brief This vector holds pointer to the elements of another container.

    If you for example want to sort the elements of a constant container, you would have to copy the whole container.
    @n To avoid copy actions this class only holds pointer to the constant elements of a container.
    @n You can insert new elements, but it is not possible to change existing ones.

    The following code demonstrates the use of this class:
    @code
    FeatureMap<> map;
    map.resize(10); //...fill map with data

    //Create pointer vector to the map
    ConstRefVector<FeatureMap<> > ref_vector(map);
    //Sort the pointer vector without changing the original data
    ref_vector.sortByIntensity();
    @endcode

    @improvement Check if we can omit the iterator template arguments (Clemens)

    @ingroup Datastructures
  */
  template <typename ContainerT>
  class ConstRefVector
  {
public:

    ///ConstIterator for the ConstRefVector
    template <class ValueT>
    class ConstRefVectorConstIterator
    {
      friend class ConstRefVector;

public:
      typedef ValueT ValueType;
      typedef ValueType value_type;
      typedef typename std::vector<ValueType*>::difference_type difference_type;
      typedef const value_type& reference;
      typedef const value_type* pointer;
      typedef std::random_access_iterator_tag iterator_category;

      ConstRefVectorConstIterator() = default;

      ConstRefVectorConstIterator(const ConstRefVectorConstIterator&) = default;

      ConstRefVectorConstIterator& operator=(const ConstRefVectorConstIterator&) = default;

      ~ConstRefVectorConstIterator() = default;

      ConstRefVectorConstIterator(const typename std::vector<ValueType*>* vec, unsigned int position)
      {
        vector_ = (typename std::vector<ValueType*>*)vec;
        position_ = position;
      }

      ConstRefVectorConstIterator(typename std::vector<ValueType*>* vec, unsigned int position)
      {
        vector_ = vec;
        position_ = position;
      }

      bool operator<(const ConstRefVectorConstIterator& it) const
      {
        return position_ < it.position_;
      }

      bool operator>(const ConstRefVectorConstIterator& it) const
      {
        return position_ > it.position_;
      }

      bool operator<=(const ConstRefVectorConstIterator& it) const
      {
        return position_ < it.position_ || position_ == it.position_;
      }

      bool operator>=(const ConstRefVectorConstIterator& it) const
      {
        return position_ > it.position_ || position_ == it.position_;
      }

      bool operator==(const ConstRefVectorConstIterator& it) const
      {
        return position_ == it.position_ && vector_ == it.vector_;
      }

      bool operator!=(const ConstRefVectorConstIterator& it) const
      {
        return position_ != it.position_ || vector_ != it.vector_;
      }

      ConstRefVectorConstIterator& operator++()
      {
        position_ += 1;
        return *this;
      }

      ConstRefVectorConstIterator operator++(int)
      {
        ConstRefVectorConstIterator tmp(*this);
        ++(*this);
        return tmp;
      }

      ConstRefVectorConstIterator& operator--()
      {
        position_ -= 1;
        return *this;
      }

      ConstRefVectorConstIterator operator--(int)
      {
        ConstRefVectorConstIterator tmp(*this);
        --(*this);
        return tmp;
      }

      ConstRefVectorConstIterator operator-(difference_type n) const
      {
        ConstRefVectorConstIterator tmp(*this);
        tmp.position_ -= n;
        return tmp;
      }

      ConstRefVectorConstIterator operator+(difference_type n) const
      {
        ConstRefVectorConstIterator tmp(*this);
        tmp.position_ += n;
        return tmp;
      }

      ConstRefVectorConstIterator& operator+=(difference_type n)
      {
        position_ += n;
        return *this;
      }

      ConstRefVectorConstIterator& operator-=(difference_type n)
      {
        position_ -= n;
        return *this;
      }

      friend difference_type operator-(const ConstRefVectorConstIterator& i1, const ConstRefVectorConstIterator& i2)
      {
        return i1.position_ - i2.position_;
      }

      friend ConstRefVectorConstIterator operator+(difference_type n, const ConstRefVectorConstIterator& i)
      {
        ConstRefVectorConstIterator tmp(i);
        tmp.position_ += n;
        return tmp;
      }

      reference operator*() const
      {
        return *((*vector_)[position_]);
      }

      pointer operator->() const
      {
        return (*vector_)[position_];
      }

protected:

      typename std::vector<ValueType*>* vector_;
      unsigned int position_;
    };


    /// Mutable iterator for the ConstRefVector
    template <class ValueT>
    class ConstRefVectorIterator :
      public ConstRefVectorConstIterator<ValueT>
    {
      friend class ConstRefVector;

public:

      typedef ValueT ValueType;
      typedef typename ConstRefVectorConstIterator<ValueType>::value_type& reference;
      typedef typename ConstRefVectorConstIterator<ValueType>::value_type* pointer;

      using ConstRefVectorConstIterator<ValueType>::vector_;
      using ConstRefVectorConstIterator<ValueType>::position_;

      ConstRefVectorIterator() = default;

      ConstRefVectorIterator(const ConstRefVectorIterator&) = default;
      
      ConstRefVectorIterator(typename std::vector<ValueType*>* vec, unsigned int position) :
        ConstRefVectorConstIterator<ValueType>(vec, position)
      {
      }

      ~ConstRefVectorIterator() = default;

      ConstRefVectorIterator& operator=(const ConstRefVectorIterator& rhs) = default;

      reference operator*() const
      {
        return *((*vector_)[position_]);
      }

      pointer operator->() const
      {
        return (*vector_)[position_];
      }

      ConstRefVectorIterator& operator++()
      {
        ConstRefVectorConstIterator<ValueType>::position_ += 1;
        return *this;
      }

      ConstRefVectorIterator operator++(int)
      {
        ConstRefVectorIterator tmp(*this);
        ++(*this);
        return tmp;
      }

      ConstRefVectorIterator& operator--()
      {
        ConstRefVectorConstIterator<ValueType>::position_ -= 1;
        return *this;
      }

      ConstRefVectorIterator operator--(int)
      {
        ConstRefVectorIterator tmp(*this);
        --(*this);
        return tmp;
      }

      ConstRefVectorIterator operator-(typename ConstRefVectorIterator::difference_type n) const
      {
        ConstRefVectorIterator tmp(*this);
        tmp.position_ -= n;
        return tmp;
      }

      ConstRefVectorIterator operator+(typename ConstRefVectorIterator::difference_type n) const
      {
        ConstRefVectorIterator tmp(*this);
        tmp.position_ += n;
        return tmp;
      }

      friend ConstRefVectorIterator operator+(typename ConstRefVectorIterator::difference_type n, const ConstRefVectorIterator& i)
      {
        ConstRefVectorIterator tmp(i);
        tmp.position_ += n;
        return tmp;
      }

      ConstRefVectorIterator& operator+=(typename ConstRefVectorIterator::difference_type n)
      {
        ConstRefVectorConstIterator<ValueType>::position_ += n;
        return *this;
      }

      ConstRefVectorIterator& operator-=(typename ConstRefVectorIterator::difference_type n)
      {
        ConstRefVectorConstIterator<ValueType>::position_ -= n;
        return *this;
      }

      friend void swap(ConstRefVectorIterator& i1, ConstRefVectorIterator& i2)
      {
        unsigned int tmp = i1.position_;
        i1.position_ = i2.position_;
        i2.position_ = tmp;
      }

    };


    ///Type definitions
    //@{
    /// Wrapped container type
    typedef ContainerT ContainerType;
    ///
    typedef typename ContainerType::value_type ValueType;
    typedef ConstRefVectorIterator<const ValueType> Iterator;
    typedef ConstRefVectorConstIterator<const ValueType> ConstIterator;
    typedef std::reverse_iterator<Iterator> ReverseIterator;
    typedef std::reverse_iterator<ConstIterator> ConstReverseIterator;
    //@}

    ///STL-compliance type definitions
    //@{
    typedef ValueType value_type;
    typedef typename ContainerType::size_type size_type;
    typedef typename ContainerType::difference_type difference_type;
    typedef typename ContainerType::reference reference;
    typedef typename ContainerType::const_reference const_reference;
    typedef typename ContainerType::pointer pointer;
    typedef Iterator iterator;
    typedef ConstIterator const_iterator;
    typedef ReverseIterator reverse_iterator;
    typedef ConstReverseIterator const_reverse_iterator;
    //@}


    /// See std::vector documentation.
    void push_back(const ValueType& x)
    {
      const ValueType* element = &x;
      vector_.push_back(element);
    }

    /// See std::vector documentation.
    void pop_back()
    {
      vector_.pop_back();
    }

    /// See std::vector documentation.
    size_type size() const
    {
      return vector_.size();
    }

    /// See std::vector documentation.
    size_type capacity() const
    {
      return std::max(vector_.size(), capacity_);
    }

    /// See std::vector documentation.
    void reserve(size_type n)
    {
      size_type cap = capacity();

      if (n > cap)
      {
        vector_.reserve(n);
        capacity_ = n;
      }
    }

    /// See std::vector documentation.
    size_type max_size() const
    {
      return vector_.max_size();
    }

    /// See std::vector documentation.
    Iterator begin()
    {
      return Iterator((std::vector<const ValueType*>*) & vector_, (unsigned int)0);
    }

    /// See std::vector documentation.
    Iterator end()
    {
      return Iterator((std::vector<const ValueType*>*) & vector_, (unsigned int)(vector_.size()));
    }

    /// See std::vector documentation.
    ConstIterator begin() const
    {
      return ConstIterator((const std::vector<const ValueType*>*) & vector_, (unsigned int)0);
    }

    /// See std::vector documentation.
    ConstIterator end() const
    {
      return ConstIterator((const std::vector<const ValueType*>*) & vector_, (unsigned int)(vector_.size()));
    }

    /// See std::vector documentation.
    ReverseIterator rbegin()
    {
      return ReverseIterator(end());
    }

    /// See std::vector documentation.
    ReverseIterator rend()
    {
      return ReverseIterator(begin());
    }

    /// See std::vector documentation.
    ConstReverseIterator rbegin() const
    {
      return ConstReverseIterator(end());
    }

    /// See std::vector documentation.
    ConstReverseIterator rend() const
    {
      return ConstReverseIterator(begin());
    }

    /// See std::vector documentation.
    void resize(size_type new_size)
    {
      vector_.resize(new_size);
      capacity_ = vector_.capacity();
    }

    /// See std::vector documentation.
    void resize(size_type new_size, const ValueType& t)
    {
      vector_.resize(new_size, &t);
      capacity_ = vector_.capacity();
    }

    /// See std::vector documentation.
    const_reference front() const
    {
      return *(begin());
    }

    /// See std::vector documentation.
    const_reference back() const
    {
      return *(end() - 1);
    }

    /// See std::vector documentation.
    void clear()
    {
      vector_.clear();
    }

    /// See std::vector documentation.
    bool empty() const
    {
      return vector_.empty();
    }

    /// See std::vector documentation.
    const_reference operator[](size_type n) const
    {
      return *(vector_[n]);
    }

    /// See std::vector documentation.
    bool operator==(const ConstRefVector& array) const
    {
      if (base_container_ptr_ != array.base_container_ptr_)
      {
        return false;
      }
      if (size() != array.size())
      {
        return false;
      }
      for (Size i = 0; i < size(); i++)
      {
        if (typeid(*(vector_[i])) != typeid(*(array.vector_[i])))
        {
          return false;
        }
        if (vector_[i]->operator!=(* array.vector_[i]))
        {
          return false;
        }
      }
      return true;
    }

    /// See std::vector documentation.
    bool operator!=(const ConstRefVector& array) const
    {
      return !(operator==(array));
    }

    /// Comparison of container sizes
    bool operator<(const ConstRefVector& array) const
    {
      return size() < array.size();
    }

    /// Comparison of container sizes
    bool operator>(const ConstRefVector& array) const
    {
      return size() > array.size();
    }

    /// Comparison of container sizes
    bool operator<=(const ConstRefVector& array) const
    {
      return operator<(array) || operator==(array);
    }

    /// Comparison of container sizes
    bool operator>=(const ConstRefVector& array) const
    {
      return operator>(array) || operator==(array);
    }

    /// See std::vector documentation.
    void swap(ConstRefVector& array)
    {
      vector_.swap(array.vector_);
    }

    /// See std::vector documentation.
    friend void swap(ConstRefVector& a1, ConstRefVector& a2)
    {
      a1.vector_.swap(a2.vector_);
    }

    /// See std::vector documentation.
    Iterator insert(Iterator pos, const ValueType& element)
    {
      const ValueType* pointer = &element;
      vector_.insert(vector_.begin() + pos.position_, pointer);
      return pos;
    }

    /// See std::vector documentation.
    void insert(Iterator pos, size_type n, const ValueType& element)
    {
      const ValueType* pointer;
      std::vector<const ValueType*> tmp;
      for (size_type i = 0; i < n; i++)
      {
        pointer = &element;
        tmp.push_back(pointer);
      }
      vector_.insert(vector_.begin() + pos.position_, tmp.begin(), tmp.end());
    }

    /// See std::vector documentation.
    template <class InputIterator>
    void insert(Iterator pos, InputIterator f, InputIterator l)
    {
      const ValueType* pointer;
      std::vector<const ValueType*> tmp;
      for (InputIterator it = f; it != l; ++it)
      {
        pointer = &(*it);
        tmp.push_back(pointer);
      }
      vector_.insert(vector_.begin() + pos.position_, tmp.begin(), tmp.end());
    }

    /// See std::vector documentation.
    Iterator erase(Iterator pos)
    {
      vector_.erase(vector_.begin() + pos.position_);
      return pos;
    }

    /// See std::vector documentation.
    Iterator erase(Iterator first, Iterator last)
    {
      vector_.erase(vector_.begin() + first.position_, vector_.begin() + last.position_);
      return first;
    }

    ///@name Constructors and Destructor
    //@{
    /// See std::vector documentation.
    ConstRefVector() :
      capacity_(0),
      base_container_ptr_(nullptr)
    {
    }

    /// See std::vector documentation.
    ConstRefVector(size_type n) :
      capacity_(0),
      base_container_ptr_(0)
    {
      vector_ = std::vector<const ValueType*>(n);
    }

    /// See std::vector documentation.
    ConstRefVector(size_type n, const ValueType& element) :
      capacity_(0),
      base_container_ptr_(0)
    {
      vector_ = std::vector<const ValueType*>(n, &element);
    }

    /// See std::vector documentation.
    ConstRefVector(const ConstRefVector& p) :
      capacity_(0),
      base_container_ptr_(p.base_container_ptr_)
    {
      const ValueType* element;
      for (ConstIterator it = p.begin(); it != p.end(); ++it)
      {
        element = &(*it);
        vector_.push_back(element);
      }
    }

    /// See std::vector documentation.
    template <class InputIterator>
    ConstRefVector(InputIterator f, InputIterator l) :
      capacity_(0),
      base_container_ptr_(nullptr)
    {
      const ValueType* pointer;
      for (InputIterator it = f; it != l; ++it)
      {
        pointer = &(*it);
        vector_.push_back(pointer);
      }
    }

    /// See std::vector documentation.
    ConstRefVector(ContainerType& p) :
      capacity_(0),
      base_container_ptr_(&p)
    {
      const ValueType* element;
      for (typename ContainerType::iterator it = p.begin(); it != p.end(); ++it)
      {
        element = &(*it);
        vector_.push_back(element);
      }
    }

    /// See std::vector documentation.
    ~ConstRefVector()
    {
    }

    //@}

    /// See std::vector documentation.
    ConstRefVector& operator=(const ConstRefVector& rhs)
    {
      if (this == &rhs) return *this;

      base_container_ptr_ = rhs.base_container_ptr_;
      clear();
      reserve(rhs.size());
      const ValueType* element;
      for (ConstIterator it = rhs.begin(); it != rhs.end(); ++it)
      {
        element = &(*it);
        vector_.push_back(element);
      }

      return *this;
    }

    /// See std::vector documentation.
    template <class InputIterator>
    void assign(InputIterator f, InputIterator l)
    {
      clear();
      insert(end(), f, l);
    }

    /// See std::vector documentation.
    void assign(size_type n, const ValueType& x)
    {
      clear();
      insert(end(), n, x);
    }

    /**
      @brief Sorting.

      These simplified sorting methods are supported for convenience.

      @note To use these methods, the ValueType must provide an interface similar to Peak1D or Peak2D.
    */
    //@{

    /// Sorts the elements by intensity
    void sortByIntensity(bool reverse = false)
    {
      if (reverse)
      {
        auto pointerCmp = [](auto& left, auto& right){typename ValueType::IntensityLess cmp; return cmp(*left, *right);};
        std::sort(vector_.begin(), vector_.end(), [&](auto &left, auto &right) {return pointerCmp(right, left);});
      }
      else
      {
        std::sort(vector_.begin(), vector_.end(), [](auto& left, auto& right){typename ValueType::IntensityLess cmp; return cmp(*left, *right);}); 
      }
    }

    /// Lexicographically sorts the elements by their position.
    void sortByPosition()
    {
      std::sort(vector_.begin(), vector_.end(), [](auto& left, auto& right){typename ValueType::PositionLess cmp; return cmp(*left, *right);}); 
    }

    //@}

    /**
      @name Generic sorting function templates.
      Any element comparator can be
      given as template argument.  You can also give the comparator as an
      argument to the function template (this is useful if the comparator is
      not default constructed, but keep in mind that STL copies comparators
      a lot).

      <p> Thus your can e.g. write <code>elements.sortByComparator <
      DFeature<1>::IntensityLess > ()</code>, if elements have type
      <code>ConstRefVector < 1, DFeature <1> ></code>.
    */
    //@{
    template <typename ComparatorType>
    void sortByComparator(ComparatorType const& comparator = ComparatorType())
    {
      std::sort(vector_.begin(), vector_.end(), [&](auto& left, auto& right){return comparator(*left, *right);});
    }

    //@}

    //----------------------------------------------------------------------

protected:

    ///the internal vector of ValueType pointers
    std::vector<const ValueType*> vector_;
    ///the current capacity
    size_type capacity_;
    /// Pointer to the base container
    const ContainerType* base_container_ptr_;
  };

} // namespace OpenMS

