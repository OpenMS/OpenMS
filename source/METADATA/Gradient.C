// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/Gradient.h>

#include <algorithm>

using namespace std;

namespace OpenMS 
{
  Gradient::Gradient(): 
 //   //PersistentObject(),
    eluents_(),
    times_(),
    percentages_()
  {
    
  }
  
  Gradient::Gradient(const Gradient& source):
    //PersistentObject(source),
    eluents_(source.eluents_),
    times_(source.times_),
    percentages_(source.percentages_)
  {
    
  }
   
  Gradient::~Gradient()
  {
    
  }
  
  Gradient& Gradient::operator = (const Gradient& source)
  {
    if (source == *this) return *this;

    eluents_ = source.eluents_;
    times_ = source.times_;
    percentages_ = source.percentages_;
    //PersistentObject::operator = (source);

    return *this;
  }

  bool Gradient::operator == (const Gradient& rhs) const
  {
    return ( eluents_ == rhs.eluents_ ) &&
    			 ( times_ == rhs.times_) &&
    			 ( percentages_ == rhs.percentages_)
    			 ;
  }

  bool Gradient::operator != (const Gradient& rhs) const
  {
    return !(operator == (rhs));
  }

	void Gradient::addEluent(const String& eluent)
	{
		//check if the eluent is valid
		std::vector<String>::iterator elu_it = find(eluents_.begin(), eluents_.end(),eluent);
		if (elu_it != eluents_.end())
		{
			throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, "A eluent with this name already exists!",eluent);
		}
		
		eluents_.push_back(eluent);
		// add zero values to percentages
		percentages_.push_back( vector < UInt >( times_.size(), 0 ) );
	}

	void Gradient::clearEluents()
	{
		eluents_.clear();
	}

	const std::vector <String>& Gradient::getEluents() const
	{
		return eluents_;
	}

	void Gradient::addTimepoint(Int timepoint)
	{
		if ((times_.size()>0) && (timepoint <= times_[times_.size()-1])) 
		{
			throw Exception::OutOfRange(__FILE__,__LINE__,__PRETTY_FUNCTION__);
		}
		times_.push_back(timepoint);
		
		// add zero values to percentages
		for (Size i = 0 ;  i < eluents_.size(); ++i)
		{
			percentages_[i].push_back(0);
		}
	}

	void Gradient::clearTimepoints()
	{
		times_.clear();
	}

	const std::vector<Int>& Gradient::getTimepoints() const
	{
		return times_;
	}

	void Gradient::setPercentage(const String& eluent, Int timepoint, UInt percentage)
	{
		// (1) validity check
		
		//check if the eluent is valid
		std::vector<String>::iterator elu_it = find(eluents_.begin(), eluents_.end(),eluent);
		if (elu_it == eluents_.end())
		{
			throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, "The given eluent does not exist in the list of eluents!",eluent);
		}

		//check if the timepoint is valid
		std::vector<Int>::iterator time_it = find(times_.begin(),times_.end(),timepoint);
		if (time_it == times_.end())
		{
			throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, "The given timepoint does not exist in the list of timepoints!",String(timepoint));
		}
		
		// percentage is valid?
		if (percentage>100)
		{
			throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, "The percentage should be between 0 and 100!",String(percentage));
		}		
		
		// (2) Look up indices
		
		UInt elu_index(0),time_index(0);
		//look up eluents index
		for (std::vector<String>::iterator it = eluents_.begin(); it != eluents_.end(); ++it)
		{
			if (*it == eluent) break;
			++elu_index;
		}
		//look up timepoint index
		for (std::vector<Int>::iterator it = times_.begin(); it != times_.end(); ++it)
		{
			if (*it == timepoint) break;
			++time_index;
		}
		
		// (3) set percentage
		percentages_[elu_index][time_index] = percentage;
	}

	const std::vector< std::vector< UInt > >& Gradient::getPercentages() const
	{
		return percentages_;
	}

	UInt Gradient::getPercentage(const String& eluent, Int timepoint) const
	{
		// (1) validity check
		
		//check if the eluent is valid
		std::vector<String>::const_iterator elu_it = find(eluents_.begin(), eluents_.end(),eluent);
		if (elu_it == eluents_.end())
		{
			throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, "The given eluent does not exist in the list of eluents!",eluent);
		}

		//check if the timepoint is valid
		std::vector<Int>::const_iterator time_it = find(times_.begin(),times_.end(),timepoint);
		if (time_it == times_.end())
		{
			throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, "The given timepoint does not exist in the list of timepoints!",String(timepoint));
		}
		
		// (2) Look up indices
		
		UInt elu_index(0),time_index(0);
		//look up eluents index
		for (std::vector<String>::const_iterator it = eluents_.begin(); it != eluents_.end(); ++it)
		{
			if (*it == eluent) break;
			++elu_index;
		}
		//look up timepoint index
		for (std::vector<Int>::const_iterator it = times_.begin(); it != times_.end(); ++it)
		{
			if (*it == timepoint) break;
			++time_index;
		}
		
		// (3) return percentage
		return percentages_[elu_index][time_index];			
	}
	
	void Gradient::clearPercentages()
	{
		percentages_.clear();
		// fill all precentages with 0
		percentages_.insert(percentages_.begin(), eluents_.size(), vector < UInt>( times_.size(), 0 ) );
	}
	
	bool Gradient::isValid() const
	{
		for (Size j = 0 ;  j < times_.size(); ++j)
		{
			Int sum =0;
			for (Size i = 0 ;  i < eluents_.size(); ++i)
			{
				sum += percentages_[i][j];
			}
			if (sum !=100)
			{
				return false;
			}
		}
		return true;
	}

} // namespace OpenMS

