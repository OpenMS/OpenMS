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
// $Id: DFeatureMap.h,v 1.19 2006/06/01 15:26:18 groepl Exp $
// $Author: groepl $
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_DPEAKREFERENCEARRAY_H
#define OPENMS_KERNEL_DPEAKREFERENCEARRAY_H

#include <OpenMS/config.h>
#include <OpenMS/KERNEL/DFeatureMap.h>
#include <OpenMS/KERNEL/KernelTraits.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>

#include <vector>

namespace OpenMS
{

  /**
    @brief A container for (composite) features.

    A map is a container holding D-dimensional features,
    which in turn represent chemical entities (peptides, proteins, etc.) found
    in a D-dimensional experiment.
    To avoid copy actions this class only holds pointer to the const elements of a container.
    You can insert new elements, but it is not possible to change existing ones. The class does not
    Maps are implemented as vectors of features and have basically the same interface
    as an STL vector has (model of Random Access Container and Back Insertion Sequence).
    Maps are typically created from peak data of 2D runs through the FeatureFinder.

    @ingroup Kernel, Serialization
  */
  template <Size D, typename TraitsT = KernelTraits, typename MapT = DFeatureMap<D, DFeature<D, TraitsT> > >
  class DPeakReferenceArray
  {
  public:

    ///ConstIterator for the DPeakArray
    template <class IteratorPeakT>
    class DPeakReferenceArrayConstIterator
    {

      friend class DPeakReferenceArray;

    public:
      /** @name Type definitions */
      //@{
      typedef IteratorPeakT IteratorPeakType;
      typedef IteratorPeakType value_type;
      typedef typename std::vector<IteratorPeakType*>::difference_type difference_type;
      typedef const value_type& reference;
      typedef const value_type* pointer;
      typedef std::random_access_iterator_tag iterator_category;
      //@}

      DPeakReferenceArrayConstIterator()
      {}

      DPeakReferenceArrayConstIterator(const typename std::vector<IteratorPeakType*>* vec , unsigned int position)
      {
        vector_ = (typename std::vector<IteratorPeakType*>*)vec;
        position_ = position;
      }

      DPeakReferenceArrayConstIterator(typename std::vector<IteratorPeakType*>* vec , unsigned int position)
      {
        vector_ = vec;
        position_ = position;
      }

      DPeakReferenceArrayConstIterator(const DPeakReferenceArrayConstIterator& it)
      {
        vector_=it.vector_;
        position_=it.position_;
      }

      ~DPeakReferenceArrayConstIterator()
      {}

      DPeakReferenceArrayConstIterator& operator = (const DPeakReferenceArrayConstIterator& rhs)
      {
        if (this==&rhs) return *this;

        vector_=rhs.vector_;
        position_=rhs.position_;

        return *this;
      }

      bool operator < (const DPeakReferenceArrayConstIterator& it) const
      {
        return position_ < it.position_;
      }

      bool operator > (const DPeakReferenceArrayConstIterator& it) const
      {
        return position_ > it.position_;
      }

      bool operator <= (const DPeakReferenceArrayConstIterator& it) const
      {
        return (position_ < it.position_ || position_ == it.position_);
      }

      bool operator >= (const DPeakReferenceArrayConstIterator& it) const
      {
        return (position_ > it.position_ || position_ == it.position_);
      }

      bool operator == (const DPeakReferenceArrayConstIterator& it) const
      {
        return position_ == it.position_ && vector_ == it.vector_;
      }

      bool operator != (const DPeakReferenceArrayConstIterator& it) const
      {
        return position_ != it.position_ || vector_ != it.vector_;
      }

      DPeakReferenceArrayConstIterator& operator ++ ()
      {
        position_ += 1;
        return *this;
      }

      DPeakReferenceArrayConstIterator operator ++ (int)
      {
        DPeakReferenceArrayConstIterator tmp(*this);
        ++(*this);
        return tmp;
      }

      DPeakReferenceArrayConstIterator& operator -- ()
      {
        position_-=1;
        return *this;
      }

      DPeakReferenceArrayConstIterator operator -- (int)
      {
        DPeakReferenceArrayConstIterator tmp(*this);
        --(*this);
        return tmp;
      }

      DPeakReferenceArrayConstIterator operator - (difference_type n) const
      {
        DPeakReferenceArrayConstIterator tmp(*this);
        tmp.position_ -= n;
        return tmp;
      }

      DPeakReferenceArrayConstIterator operator + (difference_type n) const
      {
        DPeakReferenceArrayConstIterator tmp(*this);
        tmp.position_ += n;
        return tmp;
      }

      DPeakReferenceArrayConstIterator& operator += (difference_type n)
      {
        position_ += n;
        return *this;
      }

      DPeakReferenceArrayConstIterator& operator -= (difference_type n)
      {
        position_ -= n;
        return *this;
      }

      friend difference_type operator - ( const DPeakReferenceArrayConstIterator& i1, const DPeakReferenceArrayConstIterator& i2 )
      {
        return (i1.position_ - i2.position_);
      }

      friend DPeakReferenceArrayConstIterator operator + ( difference_type n, const DPeakReferenceArrayConstIterator& i )
      {
        DPeakReferenceArrayConstIterator tmp(i);
        tmp.position_ += n;
        return tmp;
      }

      reference operator * ()
      {
        return *((*vector_)[position_]);
      }

      pointer operator -> ()
      {
        return (*vector_)[position_];
      }

      pointer operator -> () const
      {
        return (*vector_)[position_];
      }

      reference operator [] (difference_type n)
      {
        return *((*this)+n);
      }

    protected:

      typename std::vector<IteratorPeakType*>* vector_;
      unsigned int position_;
    };


    /// Mutable iterator for the DPeakArray
    template <class IteratorPeakT>
  class DPeakReferenceArrayIterator : public DPeakReferenceArrayConstIterator<IteratorPeakT>
    {
      friend class DPeakReferenceArray;

    public:
      typedef IteratorPeakT IteratorPeakType;
      typedef typename DPeakReferenceArrayConstIterator<IteratorPeakType>::value_type& reference;
      typedef typename DPeakReferenceArrayConstIterator<IteratorPeakType>::value_type* pointer;

      using DPeakReferenceArrayConstIterator<IteratorPeakType>::vector_;
      using DPeakReferenceArrayConstIterator<IteratorPeakType>::position_;


      DPeakReferenceArrayIterator()
      {}

      DPeakReferenceArrayIterator(typename std::vector<IteratorPeakType*>* vec, unsigned int position): DPeakReferenceArrayConstIterator<IteratorPeakType>(vec,position)
      {}

      DPeakReferenceArrayIterator(const DPeakReferenceArrayIterator<IteratorPeakType>& it): DPeakReferenceArrayConstIterator<IteratorPeakType>(it)
      {}

      ~DPeakReferenceArrayIterator()
      {}

      reference operator * ()
      {
        return *((*vector_)[position_]);
      }

      pointer operator -> ()
      {
        return (*vector_)[position_];
      }

      const pointer operator -> () const
      {
        return (*vector_)[position_];
      }

      typename DPeakReferenceArrayIterator::reference operator [] (typename DPeakReferenceArrayIterator::difference_type n)
      {
        return *((*this)+n);
      }

      DPeakReferenceArrayIterator& operator ++ ()
      {
        DPeakReferenceArrayConstIterator<IteratorPeakType>::position_+=1;
        return *this;
      }

      DPeakReferenceArrayIterator operator ++ (int)
      {
        DPeakReferenceArrayIterator tmp(*this);
        ++(*this);
        return tmp;
      }

      DPeakReferenceArrayIterator& operator -- ()
      {
        DPeakReferenceArrayConstIterator<IteratorPeakType>::position_-=1;
        return *this;
      }

      DPeakReferenceArrayIterator operator -- (int)
      {
        DPeakReferenceArrayIterator tmp(*this);
        --(*this);
        return tmp;
      }

      DPeakReferenceArrayIterator operator - (typename DPeakReferenceArrayIterator::difference_type n) const
      {
        DPeakReferenceArrayIterator tmp(*this);
        tmp.position_ -= n;
        return tmp;
      }

      DPeakReferenceArrayIterator operator + (typename DPeakReferenceArrayIterator::difference_type n) const
      {
        DPeakReferenceArrayIterator tmp(*this);
        tmp.position_ += n;
        return tmp;
      }

      friend DPeakReferenceArrayIterator operator + (typename DPeakReferenceArrayIterator::difference_type n, const DPeakReferenceArrayIterator& i )
      {
        DPeakReferenceArrayIterator tmp(i);
        tmp.position_ += n;
        return tmp;
      }

      DPeakReferenceArrayIterator& operator += (typename DPeakReferenceArrayIterator::difference_type n)
      {
        DPeakReferenceArrayConstIterator<IteratorPeakType>::position_ += n;
        return *this;
      }

      DPeakReferenceArrayIterator& operator -= (typename DPeakReferenceArrayIterator::difference_type n)
      {
        DPeakReferenceArrayConstIterator<IteratorPeakType>::position_ -= n;
        return *this;
      }

      friend void swap(DPeakReferenceArrayIterator& i1, DPeakReferenceArrayIterator& i2)
      {
        unsigned int tmp = i1.position_;
        i1.position_ = i2.position_;
        i2.position_ = tmp;
      }

    protected:

    };

    /**
       @name Type definitions
    */
    //@{
    typedef TraitsT TraitsType;
    typedef MapT BaseMapType;
    typedef typename BaseMapType::value_type PeakType;
    typedef DPeakReferenceArrayIterator<PeakType> Iterator;
    typedef DPeakReferenceArrayConstIterator<PeakType> ConstIterator;
    typedef std::reverse_iterator<Iterator> ReverseIterator;
    typedef std::reverse_iterator<ConstIterator> ConstReverseIterator;
    typedef typename std::vector<PeakType>::value_type value_type;
    typedef typename std::vector<PeakType>::size_type size_type;
    typedef typename std::vector<PeakType>::difference_type difference_type;
    typedef typename std::vector<PeakType>::reference reference;
    typedef typename std::vector<PeakType>::const_reference const_reference;
    typedef typename std::vector<PeakType>::pointer pointer;
    typedef Iterator iterator;
    typedef ConstIterator const_iterator;
    typedef ReverseIterator reverse_iterator;
    typedef ConstReverseIterator const_reverse_iterator;
    //@}


    /// See std::vector documentation.
    void push_back(const PeakType& x)
    {
      PeakType* feature = &x;
      vector_.push_back(feature);
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
      return std::max(vector_.size(),capacity_);
    }

    /// See std::vector documentation.
    void reserve(size_type n)
    {
      size_type cap = capacity();
      if (n>cap)
      {
        vector_.reserve(n);
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
      return Iterator((std::vector<PeakType*>*)&vector_,(unsigned int)0);
    }

    /// See std::vector documentation.
    Iterator end()
    {
      return Iterator((std::vector<PeakType*>*)&vector_,(unsigned int)(vector_.size()));
    }

    /// See std::vector documentation.
    ConstIterator begin() const
    {
      return ConstIterator((const std::vector<PeakType*>*)&vector_,(unsigned int)0);
    }

    /// See std::vector documentation.
    ConstIterator end() const
    {
      return ConstIterator((const std::vector< PeakType*>*)&vector_,(unsigned int)(vector_.size()));
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
    void resize(size_type new_size, const PeakType& t=PeakType())
    {
      size_type old_size = vector_.size();
      if (new_size<old_size)
      {
        vector_.resize(new_size);
      }
      else if (new_size>old_size)
      {
        vector_.resize(new_size);
      }
    }

    /// See std::vector documentation.
    reference front()
    {
      return *(begin());
    }

    /// See std::vector documentation.
    const_reference front() const
    {
      return *(begin());
    }

    /// See std::vector documentation.
    reference back()
    {
      return *(end()-1);
    }

    /// See std::vector documentation.
    const_reference back() const
    {
      return *(end()-1);
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
    reference operator [](size_type n)
    {
      return *(vector_[n]);
    }

    /// See std::vector documentation.
    const_reference operator [](size_type n) const
    {
      return *(vector_[n]);
    }

    /// See std::vector documentation.
    bool operator == (const DPeakReferenceArray& array) const
    {
      if (base_container_ptr_ != array.base_container_ptr_)
      {
        return false;
      }
      if (size()!=array.size())
      {
        return false;
      }
      for (unsigned int i=0;i<size();i++)
      {
        if (typeid(*(vector_[i]))!=typeid(*(array.vector_[i])))
        {
          return false;
        }
        if ( vector_[i]->operator!=(*array.vector_[i]) )
        {
          return false;
        }
      }
      return true;
    }

    /// See std::vector documentation.
    bool operator !=(const DPeakReferenceArray& array) const
    {
      return !(operator==(array));
    }

    /// Comparison of container sizes
    bool operator < (const DPeakReferenceArray& array) const
    {
      return size() < array.size();
    }

    /// Comparison of container sizes
    bool operator > (const DPeakReferenceArray& array) const
    {
      return size() > array.size();
    }

    /// Comparison of container sizes
    bool operator <= (const DPeakReferenceArray& array) const
    {
      return operator<(array) || operator==(array);
    }

    /// Comparison of container sizes
    bool operator >= (const DPeakReferenceArray& array) const
    {
      return operator>(array) || operator==(array);
    }

    /// See std::vector documentation.
    void swap(DPeakReferenceArray& array)
    {
      vector_.swap(array.vector_);
    }

    /// See std::vector documentation.
    friend void swap(DPeakReferenceArray& a1, DPeakReferenceArray& a2)
    {
      a1.vector_.swap(a2.vector_);
    }

    /// See std::vector documentation.
    Iterator insert(Iterator pos, const PeakType& feature)
    {
      PeakType* pointer = &feature;
      vector_.insert(vector_.begin()+pos.position_,pointer);
      return pos;
    }

    /// See std::vector documentation.
    void insert(Iterator pos, size_type n, const PeakType& feature)
    {
      PeakType* pointer;
      std::vector<PeakType*> tmp;
      for (size_type i=0;i<n;i++)
      {
        pointer = &feature;
        tmp.push_back(pointer);
      }
      vector_.insert(vector_.begin()+pos.position_,tmp.begin(),tmp.end());
    }

    /// See std::vector documentation.
    template <class InputIterator>
    void insert(Iterator pos, InputIterator f, InputIterator l)
    {
      PeakType* pointer;
      std::vector<PeakType*> tmp;
      for (InputIterator it=f;it!=l;++it)
      {
        pointer = (&it->clone());
        tmp.push_back(pointer);
      }
      vector_.insert(vector_.begin()+pos.position_,tmp.begin(),tmp.end());
    }

    /// See std::vector documentation.
    Iterator erase(Iterator pos)
    {
      vector_.erase(vector_.begin()+pos.position_);
      return pos;
    }

    /// See std::vector documentation.
    Iterator erase(Iterator first,Iterator last)
    {
      vector_.erase(vector_.begin()+first.position_,vector_.begin()+last.position_);
      return first;
    }

    /** @name Accesssor methods
     */
    //@{
    /// Set pointer to the base container
    void setBaseContainerPointer(const BaseMapType& base_map) { base_container_ptr_ = &base_map; }
    /// Set pointer to the base container
    void setBaseContainerPointer(const BaseMapType* base_map_pointer) { base_container_ptr_ = base_map_pointer; }
    /// Get grid
    BaseMapType* getBaseContainerPointer() { return *base_container_ptr_; }
    /// Get grid (non-mutable)
    BaseMapType const* getBaseContainerPointer() const { return *base_container_ptr_; }
    //@}


    /** @name Constructors and Destructor
      */
    //@{
    /// See std::vector documentation.
    DPeakReferenceArray()
        : capacity_(0),
    base_container_ptr_(0) {}

    /// See std::vector documentation.
    DPeakReferenceArray(size_type n)
        : capacity_(0),
        base_container_ptr_(0)
    {
      vector_=std::vector<PeakType*>(n);
    }
    /// See std::vector documentation.
    DPeakReferenceArray(size_type n, const PeakType& feature)
        : capacity_(0),
        base_container_ptr_(0)
    {
      vector_=std::vector<PeakType*>(n);
    }
    /// See std::vector documentation.
    DPeakReferenceArray(const DPeakReferenceArray& p)
        :  capacity_(0),
        base_container_ptr_(p.base_container_ptr_)
    {
      PeakType* feature;
      for (ConstIterator it=p.begin(); it!=p.end();++it)
      {
        feature = const_cast<PeakType*>(&(*it));
        vector_.push_back(feature);
      }
    }
    /// See std::vector documentation.
    template <class InputIterator>
    DPeakReferenceArray(InputIterator f, InputIterator l)
        : capacity_(0),
        base_container_ptr_(0)
    {
      PeakType* feature;
      for (InputIterator it=f;it!=l;++it)
      {
        feature = const_cast<PeakType*>(&(*it));
        vector_.push_back(feature);
      }
    }
    ///
    DPeakReferenceArray(const BaseMapType& p)
        :  capacity_(0),
        base_container_ptr_(&p)
    {
      PeakType* feature;
      for (typename BaseMapType::iterator it = p.begin(); it != p.end(); ++it)
      {
        feature = const_cast<PeakType*>(&(*it));
        vector_.push_back(feature);
      }
    }

    /// See std::vector documentation.
    ~DPeakReferenceArray()
    {}
    //@}

    /// See std::vector documentation.
    DPeakReferenceArray& operator = (const DPeakReferenceArray& rhs)
    {
      if (this==&rhs) return *this;

      base_container_ptr_ = rhs.base_container_ptr_;
      clear();
      reserve(rhs.size());
      PeakType* feature;
      for (ConstIterator it=rhs.begin(); it!=rhs.end();++it)
      {
        feature = const_cast<PeakType*>(&(*it));
        vector_.push_back(feature);
      }

      return *this;
    }

    /// See std::vector documentation.
    template <class InputIterator>
    void assign(InputIterator f , InputIterator l)
    {
      clear();
      insert(end(),f,l);
    }

    /// See std::vector documentation.
    void assign(size_type n , const PeakType& x)
    {
      clear();
      insert(end(),n,x);
    }

    /**
      @brief Sorting.

      These simplified sorting methods are supported in addition to the standard sorting methods of std::vector.
    */
    //@{

    /// Sorts the features by intensity
    void sortByIntensity()
    {
      std::sort(vector_.begin(), vector_.end(), PointerComparator < typename PeakType::IntensityLess > () );
    }

    /// Lexicographically sorts the features by their position.
    void sortByPosition()
    {
      std::sort(vector_.begin(), vector_.end(), PointerComparator<typename PeakType::PositionLess>() );
    }

    /**
      @brief Sorts the features by one dimension of their position.

      It is only sorted according to dimentsion @p i .
    */
    void sortByNthPosition(UnsignedInt i) throw (Exception::NotImplemented);

    //@}

    /**
      @name Generic sorting function templates.
      Any feature comparator can be
      given as template argument.  You can also give the comparator as an
      argument to the function template (this is useful if the comparator is
      not default constructed, but keep in mind that STL copies comparators
      a lot).

      <p> Thus your can e.g. write <code>features.sortByComparator <
      DFeature<1>::IntensityLess > ()</code>, if features have type
      <code>DPeakReferenceArray < 1, DFeature <1> ></code>.
    */
    //@{
    template < typename ComparatorType >
    void sortByComparator ( ComparatorType const & comparator )
    {
      std::sort(vector_.begin(), vector_.end(), PointerComparator < ComparatorType > ( comparator ) );
    }
    template < typename ComparatorType >
    void sortByComparator ()
    {
      std::sort(vector_.begin(), vector_.end(), PointerComparator< ComparatorType >() );
    }
    //@}

    //----------------------------------------------------------------------

  protected:

    typedef std::vector<PeakType*> InternalPointerVector;

    /**
      @brief Do not use this.

      ...unless you know what you are doing!  It is used by PointerPriorityQueue
      which has a special constructor for vector arguments, which is faster than
      using begin() and end() because it does not rely on push_back().
    */
    InternalPointerVector const & internalPointerVector() const
    {
      return vector_;
    }

    ///the internal vector of PeakType pointers
    std::vector<PeakType*> vector_;
    ///the current capacity
    size_type capacity_;
    /// Pointer to the base container
    BaseMapType* base_container_ptr_;

    ///@name Serialization
    //@{
  private:
    /// Serialization interface
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /* version */ )
    {
      ar & boost::serialization::make_nvp("vector",vector_);
      ar & boost::serialization::make_nvp("capacity",capacity_);
    }
    //@}

    /// Serialization
    friend class boost::serialization::access;

  };

  ///Print the contents to a stream.
  template <Size D, typename Feature>
  std::ostream& operator << (std::ostream& os, const DPeakReferenceArray<D, Feature>& array)
  {
    os << "-- DFEATUREARRAY BEGIN --"<<std::endl;
    for (typename DPeakReferenceArray<D, Feature>::const_iterator it = array.begin(); it!=array.end(); ++it)
    {
      os << *it << std::endl;
    }
    os << "-- DFEATUREARRAY END --"<<std::endl;
    return os;
  }

  //---------------------------------------------------------------
  //  Implementation of the inline / template functions
  //---------------------------------------------------------------

  template <Size D, typename TraitsT, typename PeakT >
  void DPeakReferenceArray<D,TraitsT,PeakT>::sortByNthPosition(UnsignedInt i) throw (Exception::NotImplemented)
  {
    OPENMS_PRECONDITION(i < Index(D), "illegal dimension")
    if (i==0)
    {
      std::sort(vector_.begin(), vector_.end(), PointerComparator< typename PeakType::template NthPositionLess<0> >() );
    }
    else if (i==1)
    {
      std::sort(vector_.begin(), vector_.end(), PointerComparator< typename PeakType::template NthPositionLess<1> >() );
    }
    else if (i==2)
    {
      std::sort(vector_.begin(), vector_.end(), PointerComparator< typename PeakType::template NthPositionLess<2> >() );
    }
    else
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,__FUNCTION__);
    }
  }


} // namespace OpenMS

#endif // OPENMS_KERNEL_DPeakReferenceArray_H

