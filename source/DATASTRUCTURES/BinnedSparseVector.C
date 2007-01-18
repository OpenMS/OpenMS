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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <OpenMS/DATASTRUCTURES/BinnedSparseVector.h>

#include <stdexcept>
#include <cassert>
#include <cmath>
#include <sstream>
#include <OpenMS/CONCEPT/Exception.h>

#include <iostream>

using namespace std;

namespace OpenMS
{

  BinnedSparseVector::DoubleProxy::DoubleProxy(BinnedSparseVector& vec,uint index)
    : vec_(vec), index_(index)
  {
  }

  // if there is a entry in the map from BinnedSparseVector, return that
  // if not it is a zero, so return 0
  BinnedSparseVector::DoubleProxy::operator double() const
  {
    double value = 0;
    map<uint,double>::const_iterator cmit = vec_.values_.find(index_);
    if ( cmit != vec_.values_.end() )
    {
      value = cmit->second;
    }
    return value;
  }
 
  // only save zeros
  BinnedSparseVector::DoubleProxy& BinnedSparseVector::DoubleProxy::operator = (const BinnedSparseVector::DoubleProxy& rhs)
  {
		if (this != &rhs)
		{
    	double value = 0;
    	map<uint, double>::const_iterator cmit = rhs.vec_.values_.find(rhs.index_);
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
  BinnedSparseVector::DoubleProxy& BinnedSparseVector::DoubleProxy::operator=(double val)
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
 
  BinnedSparseVector::BinnedSparseVector()
    :values_(),size_(0)
  {
  }
  
  BinnedSparseVector::BinnedSparseVector(int size)
    :values_(),size_(size)
  {
  }

  BinnedSparseVector::BinnedSparseVector(const BinnedSparseVector& source)
    :values_(source.values_),size_(source.size_)
  {
  }

  BinnedSparseVector& BinnedSparseVector::operator = (const BinnedSparseVector& source)
  {
		if (this != &source)
		{
    	values_ = source.values_;
    	size_ = source.size_;
		}
    return *this;
  }
  
  BinnedSparseVector::~BinnedSparseVector()
  {
  }

  uint BinnedSparseVector::nonzero_size() const
  {
    return values_.size();
  }

  uint BinnedSparseVector::size() const
  {
    return size_;
  }

  void BinnedSparseVector::push_back(double value)
  {
    operator[](size_++) = value;
  }

  const BinnedSparseVector::DoubleProxy BinnedSparseVector::operator[](uint pos) const 
  {
    assert(pos < size_);
    return BinnedSparseVector::DoubleProxy(const_cast<BinnedSparseVector&>(*this),pos);
  }

  BinnedSparseVector::DoubleProxy BinnedSparseVector::operator[](uint pos)
  {
    assert(pos < size_);
    return BinnedSparseVector::DoubleProxy(*this,pos);
  }
  
  double BinnedSparseVector::at(uint pos) const 
  {
    if (pos >= size_)
    {
      stringstream ss;
      ss << "tried to access " << pos << " while size was " << size_;
      throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__,"out_of_range in BinnedSparseVector",ss.str().c_str());
    }
    else 
    {
      return operator[](pos);
    }
  }

  void BinnedSparseVector::clear()
  {
    values_.clear();
    size_ = 0;
  }

  void BinnedSparseVector::resize(uint newsize)
  {
    // if the vector is to be smaller
    // delete all invalid entries
    if (newsize < size_)
    {
      for (map<uint,double>::iterator mit = values_.begin(); mit != values_.end();)
      {
        if (mit->first >= newsize)
        {
          uint nextvalue = (++mit)->first;
          values_.erase(--mit);
          mit = values_.find(nextvalue);
          
        }
        else ++mit;
      }
    }
    size_ = newsize;
  }
 
  BinnedSparseVector::iterator BinnedSparseVector::begin()
  {
    return BinnedSparseVectorIterator(*this,0);;
  }

  BinnedSparseVector::iterator BinnedSparseVector::end()
  {
    return BinnedSparseVectorIterator(*this,this->size());
  }
  
  BinnedSparseVector::BinnedSparseVectorConstIterator& BinnedSparseVector::BinnedSparseVectorConstIterator::hop()
  {
    //debug
    assert(valit_ != vector_.values_.end() );
    ++valit_;
    if ( valit_ == vector_.values_.end() ) position_ = vector_.size_;
    else position_ = valit_->first;
    return *this;
  }
  
  uint BinnedSparseVector::BinnedSparseVectorIterator::position() const
  {
    return position_;
  }
  
  uint BinnedSparseVector::BinnedSparseVectorConstIterator::position() const
  {
    return position_;
  }
  
  BinnedSparseVector::BinnedSparseVectorIterator& BinnedSparseVector::BinnedSparseVectorIterator::hop()
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
  
  BinnedSparseVector::const_iterator BinnedSparseVector::begin() const
  {
    return BinnedSparseVectorConstIterator(*this,0);
  }

  BinnedSparseVector::const_iterator BinnedSparseVector::end() const
  {
    return BinnedSparseVectorConstIterator(*this,this->size());
  }

  bool BinnedSparseVector::BinnedSparseVectorIterator::operator!=(const BinnedSparseVector::BinnedSparseVectorIterator& other)
  {
    return (position_ != other.position_ || &vector_ != &other.vector_);
  }

  BinnedSparseVector::BinnedSparseVectorIterator::BinnedSparseVectorIterator(BinnedSparseVector& vector, int position)
    :position_(position),vector_(vector),valit_(vector.values_.begin())
  {
  }
  
  BinnedSparseVector::BinnedSparseVectorIterator::BinnedSparseVectorIterator(const BinnedSparseVector::BinnedSparseVectorIterator& source)
    :position_(source.position_), vector_(source.vector_),valit_(source.valit_)
  {
  }

  BinnedSparseVector::BinnedSparseVectorIterator::~BinnedSparseVectorIterator()
  {
  }
  
  BinnedSparseVector::BinnedSparseVectorIterator& BinnedSparseVector::BinnedSparseVectorIterator::operator++()
  {
    ++position_;
    return *this;
  }

  BinnedSparseVector::BinnedSparseVectorIterator BinnedSparseVector::BinnedSparseVectorIterator::operator++(int)
  {
    BinnedSparseVector::BinnedSparseVectorIterator tmp(*this);
    ++position_;
    return tmp;
  }

  BinnedSparseVector::DoubleProxy BinnedSparseVector::BinnedSparseVectorIterator::operator*()
  {
    assert(position_ < vector_.size_);
    return BinnedSparseVector::DoubleProxy(this->vector_,position_);
  }
  
  bool BinnedSparseVector::BinnedSparseVectorConstIterator::operator!=(const BinnedSparseVector::BinnedSparseVectorConstIterator& other)
  {
    return (position_ != other.position_ || &vector_ != &other.vector_);
  }

  BinnedSparseVector::BinnedSparseVectorConstIterator::BinnedSparseVectorConstIterator(const BinnedSparseVector::BinnedSparseVectorIterator& source)
    :position_(source.position_), vector_(source.vector_),valit_(source.valit_)
  { 
  }
  
  BinnedSparseVector::BinnedSparseVectorConstIterator::BinnedSparseVectorConstIterator(const BinnedSparseVector& vector, int position)
    :position_(position),vector_(vector),valit_(vector.values_.begin())
  {
  }
  
  BinnedSparseVector::BinnedSparseVectorConstIterator::BinnedSparseVectorConstIterator(const BinnedSparseVector::BinnedSparseVectorConstIterator& source)
    :position_(source.position_), vector_(source.vector_),valit_(source.valit_)
  {
  }

  BinnedSparseVector::BinnedSparseVectorConstIterator::~BinnedSparseVectorConstIterator()
  {
  }
  
  BinnedSparseVector::BinnedSparseVectorConstIterator&BinnedSparseVector::BinnedSparseVectorConstIterator::operator++()
  {
    assert(position_ <= vector_.size_);
    ++position_;
    return *this;
  }

  BinnedSparseVector::BinnedSparseVectorConstIterator BinnedSparseVector::BinnedSparseVectorConstIterator::operator++(int)
  {
   BinnedSparseVector::BinnedSparseVectorConstIterator tmp(*this);
    ++position_;
    assert(position_ <= vector_.size_);
    return tmp;
  }

  double BinnedSparseVector::BinnedSparseVectorConstIterator::operator*()
  {
    assert(position_ < vector_.size_);
    if (position_ == valit_->first)
    {
      return valit_->second;
    }
    else 
    {
      map<uint,double>::const_iterator cmit = vector_.values_.find(position_);
      if ( cmit != vector_.values_.end() )
      {
        return cmit->second;
      }
      else return 0;
    }
  }
 
}
