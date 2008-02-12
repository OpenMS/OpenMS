// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <OpenMS/DATASTRUCTURES/SparseVector.h>

#include <stdexcept>
#include <cassert>
#include <cmath>
#include <sstream>
#include <OpenMS/CONCEPT/Exception.h>

#include <iostream>

using namespace std;

namespace OpenMS
{

  SparseVector::DoubleProxy::DoubleProxy(SparseVector& vec,UInt index)
    : vec_(vec), index_(index)
  {
  }

  // if there is a entry in the map from SparseVector, return that
  // if not it is a zero, so return 0
  SparseVector::DoubleProxy::operator double() const
  {
    double value = 0;
    map<UInt,double>::const_iterator cmit = vec_.values_.find(index_);
    if ( cmit != vec_.values_.end() )
    {
      value = cmit->second;
    }
    return value;
  }
 
  // only save zeros
  SparseVector::DoubleProxy& SparseVector::DoubleProxy::operator = (const SparseVector::DoubleProxy& rhs)
  {
		if (this != &rhs)
		{
    	double value = 0;
    	map<UInt, double>::const_iterator cmit = rhs.vec_.values_.find(rhs.index_);
    	if (cmit != rhs.vec_.values_.end())
    	{
      	value = cmit->second;
    	}
    	if (fabs(value) > 1e-8)
    	{
      	vec_.values_[index_] = value;
    	}
    	else 
			{
				if (vec_.values_.find(index_) != vec_.values_.end())
    		{
      		vec_.values_[index_] = value;
    		}
			}
		}
    return *this;
  }

  // only save zeros
  SparseVector::DoubleProxy& SparseVector::DoubleProxy::operator=(double val)
  {
    if (fabs(val) > 1e-8)
    {
     	vec_.values_[index_] = val;
    }
    else 
		{
			if (vec_.values_.find(index_) != vec_.values_.end())
    	{
     		vec_.values_[index_] = val;
    	}
		}
    return *this;
  }
 
  SparseVector::SparseVector()
    :values_(),size_(0)
  {
  }
  
  SparseVector::SparseVector(int size)
    :values_(),size_(size)
  {
  }

  SparseVector::SparseVector(const SparseVector& source)
    :values_(source.values_),size_(source.size_)
  {
  }

  SparseVector& SparseVector::operator = (const SparseVector& source)
  {
		if (this != &source)
		{
    	values_ = source.values_;
    	size_ = source.size_;
		}
    return *this;
  }
  
  SparseVector::~SparseVector()
  {
  }

  UInt SparseVector::nonzero_size() const
  {
    return values_.size();
  }

  UInt SparseVector::size() const
  {
    return size_;
  }

  void SparseVector::push_back(double value)
  {
    operator[](size_++) = value;
  }

  const SparseVector::DoubleProxy SparseVector::operator[](UInt pos) const 
  {
    assert(pos < size_);
    return SparseVector::DoubleProxy(const_cast<SparseVector&>(*this),pos);
  }

  SparseVector::DoubleProxy SparseVector::operator[](UInt pos)
  {
    assert(pos < size_);
    return SparseVector::DoubleProxy(*this,pos);
  }
  
  double SparseVector::at(UInt pos) const 
  {
    if (pos >= size_)
    {
      stringstream ss;
      ss << "tried to access " << pos << " while size was " << size_;
      throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__,"out_of_range in SparseVector",ss.str().c_str());
    }
    else 
    {
      return operator[](pos);
    }
  }

  void SparseVector::clear()
  {
    values_.clear();
    size_ = 0;
  }

  void SparseVector::resize(UInt newsize)
  {
    // if the vector is to be smaller
    // delete all invalid entries
    if (newsize < size_)
    {
      for (map<UInt,double>::iterator mit = values_.begin(); mit != values_.end();)
      {
        if (mit->first >= newsize)
        {
          UInt nextvalue = (++mit)->first;
          values_.erase(--mit);
          mit = values_.find(nextvalue);
          
        }
        else ++mit;
      }
    }
    size_ = newsize;
  }
 
  SparseVector::iterator SparseVector::begin()
  {
    return SparseVectorIterator(*this,0);;
  }

  SparseVector::iterator SparseVector::end()
  {
    return SparseVectorIterator(*this,this->size());
  }
  
  SparseVector::SparseVectorConstIterator& SparseVector::SparseVectorConstIterator::hop()
  {
    //debug
    assert(valit_ != vector_.values_.end() );
    ++valit_;
    if ( valit_ == vector_.values_.end() ) position_ = vector_.size_;
    else position_ = valit_->first;
    return *this;
  }
  
  UInt SparseVector::SparseVectorIterator::position() const
  {
    return position_;
  }
  
  UInt SparseVector::SparseVectorConstIterator::position() const
  {
    return position_;
  }
  
  SparseVector::SparseVectorIterator& SparseVector::SparseVectorIterator::hop()
  {
    //debug
    assert(valit_ != vector_.values_.end() );
    if (position_ != valit_->first )
    {
      position_ = vector_.values_.upper_bound(position_)->first;
    }
    else
    {
      ++valit_;
    }
    position_ = valit_->first;
    if ( valit_ == vector_.values_.end() ) position_ = vector_.size_;
    return *this;
  }
  
  SparseVector::const_iterator SparseVector::begin() const
  {
    return SparseVectorConstIterator(*this,0);
  }

  SparseVector::const_iterator SparseVector::end() const
  {
    return SparseVectorConstIterator(*this,this->size());
  }

  bool SparseVector::SparseVectorIterator::operator!=(const SparseVector::SparseVectorIterator& other)
  {
    return (position_ != other.position_ || &vector_ != &other.vector_);
  }

  SparseVector::SparseVectorIterator::SparseVectorIterator(SparseVector& vector, int position)
    :position_(position),vector_(vector),valit_(vector.values_.begin())
  {
  }
  
  SparseVector::SparseVectorIterator::SparseVectorIterator(const SparseVector::SparseVectorIterator& source)
    :position_(source.position_), vector_(source.vector_),valit_(source.valit_)
  {
  }

  SparseVector::SparseVectorIterator::~SparseVectorIterator()
  {
  }
  
  SparseVector::SparseVectorIterator& SparseVector::SparseVectorIterator::operator++()
  {
    ++position_;
    return *this;
  }

  SparseVector::SparseVectorIterator SparseVector::SparseVectorIterator::operator++(int)
  {
    SparseVector::SparseVectorIterator tmp(*this);
    ++position_;
    return tmp;
  }

  SparseVector::DoubleProxy SparseVector::SparseVectorIterator::operator*()
  {
    assert(position_ < vector_.size_);
    return SparseVector::DoubleProxy(this->vector_,position_);
  }
  
  bool SparseVector::SparseVectorConstIterator::operator!=(const SparseVector::SparseVectorConstIterator& other)
  {
    return (position_ != other.position_ || &vector_ != &other.vector_);
  }

  SparseVector::SparseVectorConstIterator::SparseVectorConstIterator(const SparseVector::SparseVectorIterator& source)
    :position_(source.position_), vector_(source.vector_),valit_(source.valit_)
  { 
  }
  
  SparseVector::SparseVectorConstIterator::SparseVectorConstIterator(const SparseVector& vector, int position)
    :position_(position),vector_(vector),valit_(vector.values_.begin())
  {
  }
  
  SparseVector::SparseVectorConstIterator::SparseVectorConstIterator(const SparseVector::SparseVectorConstIterator& source)
    :position_(source.position_), vector_(source.vector_),valit_(source.valit_)
  {
  }

  SparseVector::SparseVectorConstIterator::~SparseVectorConstIterator()
  {
  }
  
  SparseVector::SparseVectorConstIterator&SparseVector::SparseVectorConstIterator::operator++()
  {
    assert(position_ <= vector_.size_);
    ++position_;
    return *this;
  }

  SparseVector::SparseVectorConstIterator SparseVector::SparseVectorConstIterator::operator++(int)
  {
   SparseVector::SparseVectorConstIterator tmp(*this);
    ++position_;
    assert(position_ <= vector_.size_);
    return tmp;
  }

  double SparseVector::SparseVectorConstIterator::operator*()
  {
    assert(position_ < vector_.size_);
    if (position_ == valit_->first)
    {
      return valit_->second;
    }
    else 
    {
      map<UInt,double>::const_iterator cmit = vector_.values_.find(position_);
      if ( cmit != vector_.values_.end() )
      {
        return cmit->second;
      }
      else return 0;
    }
  }
 
}
