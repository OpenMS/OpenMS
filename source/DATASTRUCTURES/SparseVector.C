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
#include <OpenMS/DATASTRUCTURES/SparseVector.h>
#include <iostream>

#include <cmath>

using namespace std;

namespace OpenMS
{

  SparseVector::DoubleProxy::DoubleProxy(SparseVector& vec, uint index)
    : vec_(vec), index_(index)
  {
  }

  SparseVector::DoubleProxy::operator SparseVector::valuetype() const
  {
    if (index_ >= vec_.middle_ && index_ < vec_.middle_ + vec_.rightarray_.size())
    {
      return vec_.rightarray_.at(index_ - vec_.middle_);
    }
    else if (index_ >= vec_.firstentry_ && index_ < vec_.middle_)
    {
      return vec_.leftarray_.at(index_ - vec_.firstentry_);
    }
    else
    {
      return 0;
    }
  }

  SparseVector::DoubleProxy& SparseVector::DoubleProxy::operator = (const SparseVector::DoubleProxy& rhs)
  {
		if (this != &rhs)
		{
    	SparseVector::valuetype value = 0;
    	if (rhs.index_ >= rhs.vec_.middle_ && rhs.index_ < rhs.vec_.middle_ + rhs.vec_.rightarray_.size())
    	{
      	value = rhs.vec_.rightarray_.at(rhs.index_ - rhs.vec_.middle_);
    	}
    	else 
			{
				if (rhs.index_ >= rhs.vec_.firstentry_ && rhs.index_ < rhs.vec_.middle_)
    		{
      		value = rhs.vec_.leftarray_.at(rhs.index_ - rhs.vec_.firstentry_);
    		}
			}
    	if(fabs(value) > 1e-8)
    	{
      	vec_.insert(index_,value);
    	}
		}
    return *this;
  }

  SparseVector::DoubleProxy& SparseVector::DoubleProxy::operator = (SparseVector::valuetype val)
  {
   	if (fabs(val) > 1e-8)
   	{
     	vec_.insert(index_,val);
   	}

    return *this;
  }
  
  SparseVector::SparseVector()
    :firstentry_(0),middle_(0)
  {
  }

  /**
  \param offset beginning position (nonzero elements)
  \param startsize expected size
  */
  SparseVector::SparseVector(uint offset, uint startsize) 
  {
    middle_ = firstentry_ = offset;
    rightarray_ = vector<SparseVector::valuetype>(startsize);
    leftarray_ = vector<SparseVector::valuetype>(0);
  }

  SparseVector::SparseVector(const SparseVector& source)
    :firstentry_(source.firstentry_), middle_(source.middle_), rightarray_(source.rightarray_), leftarray_(source.leftarray_)
  {
  }
 
  SparseVector::~SparseVector()
  {
  }
  
  SparseVector& SparseVector::operator = (const SparseVector& source)
  {
		if (this != &source)
		{
	    middle_ = source.middle_;
  	  firstentry_ = source.firstentry_;
    	rightarray_ = source.rightarray_;
    	leftarray_ = source.leftarray_;
		}
    return *this;
  }
  
  //until operator[] cant detect write access it has to be explicit
  uint SparseVector::insert(uint pos,SparseVector::valuetype value)
  {
    if (pos > middle_+rightarray_.size() && fabs(value) < 1e-8 ) return 0;
    else if (pos < middle_ -leftarray_.size() && fabs(value) < 1e-8 ) return 0;
    if (rightarray_.size() == 0 && leftarray_.size() == 0) // first insert should go in middle
    {
      middle_ = firstentry_= pos;
      rightarray_ = vector<SparseVector::valuetype>(10); // bigger?
    }
    if ( pos >= middle_)  //rightarray
    {
      if (pos >= middle_ + rightarray_.size())
      {
        growback(pos - middle_ - rightarray_.size() +1);
      }
      rightarray_.at(pos - middle_) = value;
    }
    else //leftarray
    {
      if (pos < firstentry_) 
      {
        growfront(firstentry_ - pos);
      }
      leftarray_.at(pos - firstentry_ ) = value;
    }
    return 0;
  }

  const SparseVector::DoubleProxy SparseVector::operator[](uint pos) const
  {
    return SparseVector::DoubleProxy(const_cast<SparseVector&>(*this),pos);
  }

  SparseVector::DoubleProxy SparseVector::operator[](uint pos)
  {
    return SparseVector::DoubleProxy(*this,pos);
  }

  void SparseVector::growfront(uint add)
  {
    firstentry_ -= add;
    leftarray_.resize(leftarray_.size() + add);
  }

  void SparseVector::growback(uint add)
  {
    rightarray_.resize(rightarray_.size() + add);
  }

}
