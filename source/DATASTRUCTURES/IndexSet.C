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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/IndexSet.h>
#include <algorithm>
#include <iostream>

using namespace std;

namespace OpenMS
{
	const UnsignedInt IndexSet::DOTS = std::numeric_limits<UnsignedInt>::max();
	const UnsignedInt IndexSet::END = std::numeric_limits<UnsignedInt>::max();

	IndexSet::IndexSet() : set_(), is_sorted_(true), size_(0)
	{}

	IndexSet::~IndexSet(){}  

	IndexSet::IndexSet(const IndexSet& source) : 
		set_(source.set_),
		is_sorted_(source.is_sorted_), 
		size_(source.size_)
	{}
    
	IndexSet& IndexSet::operator = (const IndexSet& source)
	{
		set_ = source.set_;
		is_sorted_ = source.is_sorted_;
		size_ = source.size_;
		return *this;
	}
	
	IndexSet::IndexSet(UnsignedInt index): set_(), is_sorted_(true), size_(0)
	{
		add(index);
	}
	
	IndexSet::IndexSet(UnsignedInt index_from, UnsignedInt index_to): set_(), is_sorted_(true), size_(0)
	{
		add(index_from,index_to);
	}

	bool IndexSet::operator == (const IndexSet& rhs) const
	{
		if (!is_sorted_) throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The set is not sorted! It has to be in order to access it!");
		return rhs.set_ == set_;
	}
	
	bool IndexSet::operator != (const IndexSet& rhs) const
	{
		return !(operator == (rhs));		
	}

	bool IndexSet::isEmpty()
	{
		return set_.empty();	
	}


	void IndexSet::clear()
	{
		is_sorted_ = true;
		size_ = 0;
		set_.clear();	
	}

  IndexSet& IndexSet::add(UnsignedInt index)
	{
		if (is_sorted_ && set_.size()>0 && index < set_[set_.size()-1]) 
			is_sorted_ = false;

		++size_;
		set_.push_back(index);
		return *this;
	}


  IndexSet& IndexSet::add(UnsignedInt index_from, UnsignedInt index_to)
	{
		if (index_from > index_to)
			add(index_to, index_from);
		else
		{
			if (is_sorted_ && set_.size()>0 && index_from < set_[set_.size()-1]) 
				is_sorted_ = false;

			size_ += index_to-index_from+1;
			set_.push_back(index_from);			
			if (index_from != index_to) // skip if: add(n,n)
			{
				if (index_from+1 < index_to) set_.push_back(DOTS);      // skip if: add(n,n+1)
				set_.push_back(index_to);
			}
		}
		return *this;
	}

  IndexSet& IndexSet::remove(UnsignedInt index)
	{
		sort_(index,index);
		return *this;
	}

  IndexSet& IndexSet::remove(UnsignedInt index_from, UnsignedInt index_to)
	{
		if (index_from > index_to)
			remove(index_to, index_from);
		else
			sort_(index_from,index_to);
		return *this;
	}

	IndexSet::ConstIterator IndexSet::begin() const throw(Exception::Precondition)
	{
		if (!is_sorted_) throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The set is not sorted! It has to be in order to access it!");
		if (set_.size()==0) //empty IndexSet
			return ConstIterator(END,END,this);
		else
			return ConstIterator(set_[0],0,this);
	}

	IndexSet::ConstIterator IndexSet::begin(UnsignedInt index) const throw(Exception::Precondition)
	{
		if (!is_sorted_) throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The set is not sorted! It has to be in order to access it!");
		ConstIterator it(set_[0],0,this);
		while (it!=end() && *it<index) it++;
		return it;
	}

	IndexSet::ConstIterator IndexSet::end() const throw(Exception::Precondition)
	{
		if (!is_sorted_) throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The set is not sorted! It has to be in order to access it!");

		// index_ = "one bigger than last Index" for comperators
		// pos_ = "set size" for decrement
		return ConstIterator(END,END,this);
	}

	void IndexSet::sort()
	{
		if (!is_sorted_) sort_(); 
	}

	void IndexSet::sort_(UnsignedInt skip_from, UnsignedInt skip_to)
	{
		is_sorted_ = true;
		vector<UnsignedInt> tmp;
		for (const_iterator it=begin(); it!=end(); it++)	
			if ( skip_from == END || *it < skip_from || *it > skip_to) //user defined skip
				tmp.push_back(*it);

		std::sort(tmp.begin(),tmp.end());
		vector<UnsignedInt>::const_iterator vec_it=tmp.begin();
		UnsignedInt actual;
		UnsignedInt start=0;
		set_.clear();
		size_ = 0;

		while(vec_it!=tmp.end())  // compress
		{
			actual = *vec_it;
			++vec_it;

			if ( vec_it==tmp.end() || actual+1 != (*vec_it) ) // not consecutive indices
			{
				if (vec_it!=tmp.end() && actual == (*vec_it))  //skip replicates
					continue;

				if (start==0)  // no index range started
					add(actual);
				else           // finish index range
				{
					add(start,actual);
					start = 0;
				}
			}
			else if (start==0)
				start = actual;  // start new index range if not started yet

		}
		is_sorted_ = true;
	}

	std::ostream& operator << (std::ostream& os, const IndexSet& set)
 	{
		UnsignedInt size = set.set_.size();
		for (UnsignedInt i=0; i<size; ++i)
		{
			if (set.set_[i] == IndexSet::DOTS)
				os << "..";
			else
			{
				os << set.set_[i];
				if (i+1<size && set.set_[i+1] != IndexSet::DOTS)	os << ",";
			}
		}
		return os;
	}

	Size IndexSet::size() const
	{
		return size_;
	}

}
