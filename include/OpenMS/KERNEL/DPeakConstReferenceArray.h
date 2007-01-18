// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_DPEAKCONSTREFERENCEARRAY_H
#define OPENMS_KERNEL_DPEAKCONSTREFERENCEARRAY_H

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/FORMAT/Serialization.h>
#include <OpenMS/CONCEPT/Macros.h>

#include <vector>

namespace OpenMS
{

  /**
    @brief This container holds pointer to the elements of another container.
    
    If you for example want to sort the elements of a constant container, you have to copy the whole container.
    To avoid copy actions this class only holds pointer to the constant elements of a container. 
    It behaves like a DPeakArray. You can insert new elements, but it is not possible to change existing ones.
    (E.g. generating a DPeakConstReferenceArray pointer_array of a DFeatureMap feature_map is done by:
    pointer_array(feature_map.begin(),feature_map.end()))
    
    @ingroup Kernel, Serialization
  */
  template <typename MapT>
  class DPeakConstReferenceArray
  {
  public:

    ///ConstIterator for the DPeakArray
    template <class IteratorPeakT>
    class DPeakConstReferenceArrayConstIterator
    {

      friend class DPeakConstReferenceArray;

    public:
      typedef IteratorPeakT IteratorPeakType;
      typedef IteratorPeakType value_type;
      typedef typename std::vector<IteratorPeakType*>::difference_type difference_type;
      typedef const value_type& reference;
      typedef const value_type* pointer;
      typedef std::random_access_iterator_tag iterator_category;

      DPeakConstReferenceArrayConstIterator()
      {}

      DPeakConstReferenceArrayConstIterator(const typename std::vector<IteratorPeakType*>* vec , unsigned int position)
      {
        vector_ = (typename std::vector<IteratorPeakType*>*)vec;
        position_ = position;
      }

      DPeakConstReferenceArrayConstIterator(typename std::vector<IteratorPeakType*>* vec , unsigned int position)
      {
        vector_ = vec;
        position_ = position;
      }

      DPeakConstReferenceArrayConstIterator(const DPeakConstReferenceArrayConstIterator& it)
      {
        vector_=it.vector_;
        position_=it.position_;
      }

      ~DPeakConstReferenceArrayConstIterator()
      {}

      DPeakConstReferenceArrayConstIterator& operator = (const DPeakConstReferenceArrayConstIterator& rhs)
      {
        if (this==&rhs) return *this;

        vector_=rhs.vector_;
        position_=rhs.position_;

        return *this;
      }

      bool operator < (const DPeakConstReferenceArrayConstIterator& it) const
      {
        return position_ < it.position_;
      }

      bool operator > (const DPeakConstReferenceArrayConstIterator& it) const
      {
        return position_ > it.position_;
      }

      bool operator <= (const DPeakConstReferenceArrayConstIterator& it) const
      {
        return (position_ < it.position_ || position_ == it.position_);
      }

      bool operator >= (const DPeakConstReferenceArrayConstIterator& it) const
      {
        return (position_ > it.position_ || position_ == it.position_);
      }

      bool operator == (const DPeakConstReferenceArrayConstIterator& it) const
      {
        return position_ == it.position_ && vector_ == it.vector_;
      }

      bool operator != (const DPeakConstReferenceArrayConstIterator& it) const
      {
        return position_ != it.position_ || vector_ != it.vector_;
      }

      DPeakConstReferenceArrayConstIterator& operator ++ ()
      {
        position_ += 1;
        return *this;
      }

      DPeakConstReferenceArrayConstIterator operator ++ (int)
      {
        DPeakConstReferenceArrayConstIterator tmp(*this);
        ++(*this);
        return tmp;
      }

      DPeakConstReferenceArrayConstIterator& operator -- ()
      {
        position_-=1;
        return *this;
      }

      DPeakConstReferenceArrayConstIterator operator -- (int)
      {
        DPeakConstReferenceArrayConstIterator tmp(*this);
        --(*this);
        return tmp;
      }

      DPeakConstReferenceArrayConstIterator operator - (difference_type n) const
      {
        DPeakConstReferenceArrayConstIterator tmp(*this);
        tmp.position_ -= n;
        return tmp;
      }

      DPeakConstReferenceArrayConstIterator operator + (difference_type n) const
      {
        DPeakConstReferenceArrayConstIterator tmp(*this);
        tmp.position_ += n;
        return tmp;
      }

      DPeakConstReferenceArrayConstIterator& operator += (difference_type n)
      {
        position_ += n;
        return *this;
      }

      DPeakConstReferenceArrayConstIterator& operator -= (difference_type n)
      {
        position_ -= n;
        return *this;
      }

      friend difference_type operator - ( const DPeakConstReferenceArrayConstIterator& i1, const DPeakConstReferenceArrayConstIterator& i2 )
      {
        return (i1.position_ - i2.position_);
      }

      friend DPeakConstReferenceArrayConstIterator operator + ( difference_type n, const DPeakConstReferenceArrayConstIterator& i )
      {
        DPeakConstReferenceArrayConstIterator tmp(i);
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

      //      reference operator [] (difference_type n)
      //      {
      //        return *((*this)+n);
      //      }

    protected:

      typename std::vector<IteratorPeakType*>* vector_;
      unsigned int position_;
    };


    /// Mutable iterator for the DPeakArray
    template <class IteratorPeakT>
  class DPeakConstReferenceArrayIterator : public DPeakConstReferenceArrayConstIterator<IteratorPeakT>
    {
      friend class DPeakConstReferenceArray;

    public:
      typedef IteratorPeakT IteratorPeakType;
      typedef typename DPeakConstReferenceArrayConstIterator<IteratorPeakType>::value_type& reference;
      typedef typename DPeakConstReferenceArrayConstIterator<IteratorPeakType>::value_type* pointer;

      using DPeakConstReferenceArrayConstIterator<IteratorPeakType>::vector_;
      using DPeakConstReferenceArrayConstIterator<IteratorPeakType>::position_;


      DPeakConstReferenceArrayIterator()
      {}

      DPeakConstReferenceArrayIterator(typename std::vector<IteratorPeakType*>* vec, unsigned int position): DPeakConstReferenceArrayConstIterator<IteratorPeakType>(vec,position)
      {}

      DPeakConstReferenceArrayIterator(const DPeakConstReferenceArrayIterator<IteratorPeakType>& it): DPeakConstReferenceArrayConstIterator<IteratorPeakType>(it)
      {}

      ~DPeakConstReferenceArrayIterator()
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

      //      typename DPeakConstReferenceArrayIterator::reference operator [] (typename DPeakConstReferenceArrayIterator::difference_type n)
      //      {
      //        return *((*this)+n);
      //      }

      DPeakConstReferenceArrayIterator& operator ++ ()
      {
        DPeakConstReferenceArrayConstIterator<IteratorPeakType>::position_+=1;
        return *this;
      }

      DPeakConstReferenceArrayIterator operator ++ (int)
      {
        DPeakConstReferenceArrayIterator tmp(*this);
        ++(*this);
        return tmp;
      }

      DPeakConstReferenceArrayIterator& operator -- ()
      {
        DPeakConstReferenceArrayConstIterator<IteratorPeakType>::position_-=1;
        return *this;
      }

      DPeakConstReferenceArrayIterator operator -- (int)
      {
        DPeakConstReferenceArrayIterator tmp(*this);
        --(*this);
        return tmp;
      }

      DPeakConstReferenceArrayIterator operator - (typename DPeakConstReferenceArrayIterator::difference_type n) const
      {
        DPeakConstReferenceArrayIterator tmp(*this);
        tmp.position_ -= n;
        return tmp;
      }

      DPeakConstReferenceArrayIterator operator + (typename DPeakConstReferenceArrayIterator::difference_type n) const
      {
        DPeakConstReferenceArrayIterator tmp(*this);
        tmp.position_ += n;
        return tmp;
      }

      friend DPeakConstReferenceArrayIterator operator + (typename DPeakConstReferenceArrayIterator::difference_type n, const DPeakConstReferenceArrayIterator& i )
      {
        DPeakConstReferenceArrayIterator tmp(i);
        tmp.position_ += n;
        return tmp;
      }

      DPeakConstReferenceArrayIterator& operator += (typename DPeakConstReferenceArrayIterator::difference_type n)
      {
        DPeakConstReferenceArrayConstIterator<IteratorPeakType>::position_ += n;
        return *this;
      }

      DPeakConstReferenceArrayIterator& operator -= (typename DPeakConstReferenceArrayIterator::difference_type n)
      {
        DPeakConstReferenceArrayConstIterator<IteratorPeakType>::position_ -= n;
        return *this;
      }

      friend void swap(DPeakConstReferenceArrayIterator& i1, DPeakConstReferenceArrayIterator& i2)
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
    typedef MapT BaseMapType;
    typedef typename BaseMapType::value_type PeakType;
    enum { DIMENSION=PeakType::DIMENSION };
    typedef typename PeakType::TraitsType TraitsType;
    typedef DPeakConstReferenceArrayIterator<const PeakType> Iterator;
    typedef DPeakConstReferenceArrayConstIterator<const PeakType> ConstIterator;
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
      const PeakType* element = &x;
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
      return std::max(vector_.size(),capacity_);
    }

    /// See std::vector documentation.
    void reserve(size_type n)
    {
      size_type cap = capacity();

      if (n>cap)
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
      return Iterator((std::vector<const PeakType*>*)&vector_,(unsigned int)0);
    }

    /// See std::vector documentation.
    Iterator end()
    {
      return Iterator((std::vector<const PeakType*>*)&vector_,(unsigned int)(vector_.size()));
    }

    /// See std::vector documentation.
    ConstIterator begin() const
    {
      return ConstIterator((const std::vector<const PeakType*>*)&vector_,(unsigned int)0);
    }

    /// See std::vector documentation.
    ConstIterator end() const
    {
      return ConstIterator((const std::vector< const PeakType*>*)&vector_,(unsigned int)(vector_.size()));
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
    void resize(size_type new_size, const PeakType& t)
    {
      vector_.resize(new_size,&t); 
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
    //     reference operator [](size_type n)
    //     {
    //       return *(vector_[n]);
    //     }

    /// See std::vector documentation.
    const_reference operator [](size_type n) const
    {
      return *(vector_[n]);
    }

    /// See std::vector documentation.
    bool operator == (const DPeakConstReferenceArray& array) const
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
    bool operator !=(const DPeakConstReferenceArray& array) const
    {
      return !(operator==(array));
    }

    /// Comparison of container sizes
    bool operator < (const DPeakConstReferenceArray& array) const
    {
      return size() < array.size();
    }

    /// Comparison of container sizes
    bool operator > (const DPeakConstReferenceArray& array) const
    {
      return size() > array.size();
    }

    /// Comparison of container sizes
    bool operator <= (const DPeakConstReferenceArray& array) const
    {
      return operator<(array) || operator==(array);
    }

    /// Comparison of container sizes
    bool operator >= (const DPeakConstReferenceArray& array) const
    {
      return operator>(array) || operator==(array);
    }

    /// See std::vector documentation.
    void swap(DPeakConstReferenceArray& array)
    {
      vector_.swap(array.vector_);
    }

    /// See std::vector documentation.
    friend void swap(DPeakConstReferenceArray& a1, DPeakConstReferenceArray& a2)
    {
      a1.vector_.swap(a2.vector_);
    }

    /// See std::vector documentation.
    Iterator insert(Iterator pos, const PeakType& element)
    {
      const PeakType* pointer = &element;
      vector_.insert(vector_.begin()+pos.position_,pointer);
      return pos;
    }

    /// See std::vector documentation.
    void insert(Iterator pos, size_type n, const PeakType& element)
    {
      const PeakType* pointer;
      std::vector<const PeakType*> tmp;
      for (size_type i=0;i<n;i++)
      {
        pointer = &element;
        tmp.push_back(pointer);
      }
      vector_.insert(vector_.begin()+pos.position_,tmp.begin(),tmp.end());
    }

    /// See std::vector documentation.
    template <class InputIterator>
    void insert(Iterator pos, InputIterator f, InputIterator l)
    {
      const PeakType* pointer;
      std::vector<const PeakType*> tmp;
      for (InputIterator it=f;it!=l;++it)
      {
        pointer = &(*it);
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
    DPeakConstReferenceArray()
        : capacity_(0),
    base_container_ptr_(0) {}

    /// See std::vector documentation.
    DPeakConstReferenceArray(size_type n)
        : capacity_(0),
        base_container_ptr_(0)
    {
      vector_=std::vector<const PeakType*>(n);
    }
    /// See std::vector documentation.
    DPeakConstReferenceArray(size_type n, const PeakType& element)
        : capacity_(0),
        base_container_ptr_(0)
    {
      vector_=std::vector<const PeakType*>(n,&element);
    }
    /// See std::vector documentation.
    DPeakConstReferenceArray(const DPeakConstReferenceArray& p)
        :  capacity_(0),
        base_container_ptr_(p.base_container_ptr_)
    {
      const PeakType* element;
      for (ConstIterator it=p.begin(); it!=p.end();++it)
      {
        element = &(*it);
        vector_.push_back(element);
      }
    }
    /// See std::vector documentation.
    template <class InputIterator>
    DPeakConstReferenceArray(InputIterator f, InputIterator l)
        : capacity_(0),
        base_container_ptr_(0)
    {
      const PeakType* pointer;
      for (InputIterator it=f;it!=l;++it)
      {
        pointer = &(*it);
        vector_.push_back(pointer);
      }
    }
    ///
    DPeakConstReferenceArray(BaseMapType& p)
        :  capacity_(0),
        base_container_ptr_(&p)
    {
      const PeakType* element;
      for (typename BaseMapType::iterator it = p.begin(); it != p.end(); ++it)
      {
        element = &(*it);
        vector_.push_back(element);
      }
    }

    /// See std::vector documentation.
    ~DPeakConstReferenceArray()
    {}
    //@}

    /// See std::vector documentation.
    DPeakConstReferenceArray& operator = (const DPeakConstReferenceArray& rhs)
    {
      if (this==&rhs) return *this;

      base_container_ptr_ = rhs.base_container_ptr_;
      clear();
      reserve(rhs.size());
      const PeakType* element;
      for (ConstIterator it=rhs.begin(); it!=rhs.end();++it)
      {
        element = &(*it);
        vector_.push_back(element);
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

    /// Sorts the elements by intensity
    void sortByIntensity()
    {
      std::sort(vector_.begin(), vector_.end(), PointerComparator < typename PeakType::IntensityLess > () );
    }

    /// Lexicographically sorts the elements by their position.
    void sortByPosition()
    {
      std::sort(vector_.begin(), vector_.end(), PointerComparator<typename PeakType::PositionLess>() );
    }

    /**
      @brief Sorts the elements by one dimension of their position.

      It is only sorted according to dimentsion @p i .
    */
    void sortByNthPosition(UnsignedInt i) throw (Exception::NotImplemented);

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
      <code>DPeakConstReferenceArray < 1, DFeature <1> ></code>.
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

    typedef std::vector<const PeakType*> InternalPointerVector;

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
    std::vector<const PeakType*> vector_;
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
  template <typename MapT>
  std::ostream& operator << (std::ostream& os, const DPeakConstReferenceArray<MapT>& array)
  {
    os << "-- DFEATUREARRAY BEGIN --"<<std::endl;
    for (typename DPeakConstReferenceArray<MapT>::const_iterator it = array.begin();
         it!=array.end();
         ++it
        )
    {
      os << *it << std::endl;
    }
    os << "-- DFEATUREARRAY END --"<<std::endl;
    return os;
  }

  //---------------------------------------------------------------
  //  Implementation of the inline / template functions
  //---------------------------------------------------------------

  template <typename MapT>
  void DPeakConstReferenceArray<MapT>::sortByNthPosition(UnsignedInt i)
  throw (Exception::NotImplemented)
  {
    OPENMS_PRECONDITION(i < Index(PeakType::DIMENSION), "illegal dimension")
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

#endif // OPENMS_KERNEL_DPEAKCONSTREFERENCEARRAY_H

