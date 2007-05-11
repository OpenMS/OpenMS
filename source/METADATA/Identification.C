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
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <sstream>
#include <algorithm>

using namespace std;

namespace OpenMS {

  Identification::Identification()
    : MetaInfoInterface(),
			higher_score_better_(true),
			protein_significance_threshold_(0.0)
  {
  }

  Identification::Identification(const Identification& source) 
		: MetaInfoInterface(source), 
			id_(source.id_),
			search_engine_(source.search_engine_),
			search_engine_version_(source.search_engine_version_),
			search_parameters_(source.search_parameters_),
			date_(source.date_),
			protein_score_type_(source.protein_score_type_),
			higher_score_better_(source.higher_score_better_),
			protein_hits_(source.protein_hits_),
	  	protein_significance_threshold_(source.protein_significance_threshold_) 
  {
  }
  
  Identification::~Identification() 
  {
  }
 
	void Identification::setDateTime(const DateTime& date)
	{
		date_ = date;
	}

	const DateTime& Identification::getDateTime() const
	{
		return date_;
	}

	void Identification::setHits(const vector<ProteinHit>& protein_hits)
	{
		protein_hits_ = protein_hits;
	}

	const vector<ProteinHit>& Identification::getHits() const
	{
		return protein_hits_;
	}

	// retrival of the peptide significance threshold value
  Real Identification::getSignificanceThreshold() const 
  { 
  	return protein_significance_threshold_;
  }

	// setting of the peptide significance threshold value
	void Identification::setSignificanceThreshold(Real value) 
	{ 
		protein_significance_threshold_ = value;
	}
 
	void Identification::setScoreType(const String& type)
	{
		protein_score_type_ = type;
	}

	const String& Identification::getScoreType() const
	{
		return protein_score_type_;
	}

	void Identification::insertHit(const ProteinHit& protein_hit)
	{
		protein_hits_.push_back(protein_hit);
	}

  Identification& Identification::operator=(const Identification& source) 
  {
  	if (this == &source)
  	{
  		return *this;		
  	}
    MetaInfoInterface::operator=(source);
		id_ = source.id_;
		search_engine_ = source.search_engine_;
		search_engine_version_ = source.search_engine_version_;
		search_parameters_ = source.search_parameters_;
		date_ = source.date_;
    protein_hits_ = source.protein_hits_;
		protein_score_type_ = source.protein_score_type_;
    protein_significance_threshold_ = source.protein_significance_threshold_;
		higher_score_better_ = source.higher_score_better_;
    return *this;  
  }

	// Equality operator
	bool Identification::operator == (const Identification& rhs) const
	{
		return 	MetaInfoInterface::operator==(rhs) &&
						id_ == rhs.id_ &&
						search_engine_ == rhs.search_engine_ &&
				    search_engine_version_ == rhs.search_engine_version_ &&
						search_parameters_ == rhs.search_parameters_ &&
						date_ == rhs.date_ &&
						protein_hits_ == rhs.protein_hits_ &&
						protein_score_type_ == rhs.protein_score_type_ &&
						protein_significance_threshold_ == rhs.protein_significance_threshold_ &&
						higher_score_better_ == rhs.higher_score_better_;

	}
		
	// Inequality operator
	bool Identification::operator != (const Identification& rhs) const
	{
		return !operator==(rhs);						 
	}
	
  void Identification::assignRanks()
  {
    UInt rank = 1;
    sort();
    for ( vector<ProteinHit>::iterator lit = protein_hits_.begin(); lit != protein_hits_.end(); ++lit )
    {
      lit->setRank(rank++);
    }
  }
    
  void Identification::sort()
  {
	
 		if (higher_score_better_)
  	{
			std::sort(protein_hits_.begin(), protein_hits_.end(), PeptideHit::ScoreMore());
  	}
  	else
  	{
  		std::sort(protein_hits_.begin(), protein_hits_.end(), PeptideHit::ScoreLess());
  	}
  }

	bool Identification::isHigherScoreBetter() const
	{
		return higher_score_better_;
	} 
	   
	void Identification::setHigherScoreBetter(bool value)
	{
		higher_score_better_ = value;
	} 

	const String& Identification::getIdentifier() const
	{
		return id_;
	} 
	
	void Identification::setIdentifier(const String& id)
	{
		id_ = id;
	}

	void Identification::setSearchEngine(const String& search_engine)
	{
		search_engine_ = search_engine;
	}
	
	const String& Identification::getSearchEngine() const
	{
		return search_engine_;
	}
	
	void Identification::setSearchEngineVersion(const String& search_engine_version)
	{
		search_engine_version_ = search_engine_version;
	}
	
	const String& Identification::getSearchEngineVersion() const
	{
		return search_engine_version_;
	}
	
	void Identification::setSearchParameters(const SearchParameters& search_parameters)
	{
		search_parameters_ = search_parameters;
	}
	
	const Identification::SearchParameters& Identification::getSearchParameters() const
	{
		return search_parameters_;
	}

}// namespace OpenMS
