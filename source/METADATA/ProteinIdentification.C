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
// $Id: ProteinIdentification.C,v 1.2 2006/06/10 06:40:18 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <sstream>
#include <algorithm>

using namespace std;

namespace OpenMS {

  ProteinIdentification::ProteinIdentification()
    : date_(),
    	protein_hits_(), 
    	protein_significance_threshold_(0)
  {
  }

  ProteinIdentification::ProteinIdentification(const ProteinIdentification& source) : 
  	date_(source.date_), 
  	protein_hits_(source.protein_hits_), 
  	protein_significance_threshold_(source.protein_significance_threshold_)
  {
  }
  
  ProteinIdentification::~ProteinIdentification() 
  {
  	protein_hits_.clear();
  }

  const std::vector<ProteinHit>& ProteinIdentification::getProteinHits() const 
  { 
  	return protein_hits_;
  }

  float ProteinIdentification::getProteinSignificanceThreshold() const
  { 
  	return protein_significance_threshold_;
  }

	void ProteinIdentification::setProteinSignificanceThreshold(float value) 
	{ 
		protein_significance_threshold_ = value;
	}

  DateTime& ProteinIdentification::getDateTime()
  {
    return date_;
  }

  const DateTime& ProteinIdentification::getDateTime() const
  {
    return date_;
  }

  void ProteinIdentification::setDateTime(const DateTime& date)
  {
    date_ = date;
  }

  void ProteinIdentification::clear()
  {
   date_.clear();
   protein_significance_threshold_ = 0;
   protein_hits_.clear();
  }
  
  ProteinIdentification& ProteinIdentification::operator=(const ProteinIdentification& source) 
  {
  	if (this == &source)
  	{
  		return *this;		
  	}
    //PersistentObject::operator=(source);
    date_ = source.date_;
    protein_hits_ = source.protein_hits_;
    protein_significance_threshold_ = source.protein_significance_threshold_;
    return *this;  
  }

	bool ProteinIdentification::operator == (const ProteinIdentification& rhs) const
	{
		return date_ == rhs.getDateTime() 
						&& protein_hits_ == rhs.getProteinHits()
						&& protein_significance_threshold_ == rhs.getProteinSignificanceThreshold();						 
	}
		
	bool ProteinIdentification::operator != (const ProteinIdentification& rhs) const
	{
		return date_ != rhs.getDateTime() 
						|| protein_hits_ != rhs.getProteinHits()
						|| protein_significance_threshold_ != rhs.getProteinSignificanceThreshold();						 
	}
	
  void ProteinIdentification::insertProteinHit(const ProteinHit& input)
  {
    if ( !protein_hits_.size() )
    {  
      protein_hits_.push_back(input);
    }
    else if ( protein_hits_.begin()->getScoreType() == input.getScoreType() )
    {
      protein_hits_.push_back(input);
    }
    else
    {
      stringstream ss;
      ss << protein_hits_.begin()->getScoreType() << " != " <<  input.getScoreType();
      throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Incompatible ProteinHit.score_type",ss.str().c_str());
    }
  }

  void ProteinIdentification::assignRanks()
  {
    sort();
  }
    
  void ProteinIdentification::sort()
  {		
		std::sort(protein_hits_.begin(), protein_hits_.end(), RankLess());
  }
  
	bool ProteinIdentification::empty() const
	{
		DateTime temp_date;
		
		return (protein_significance_threshold_ == 0
						&& protein_hits_.size() == 0
						&& date_ == temp_date);		
	}

}// namespace OpenMS
