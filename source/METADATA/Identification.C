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
#include <sstream>
#include <algorithm>

using namespace std;

namespace OpenMS {

  Identification::Identification()
    : ProteinIdentification(),
    	peptide_hits_(), 
    	peptide_significance_threshold_(0) 
  {
  }

  Identification::Identification(const Identification& source) : 
  	ProteinIdentification(source), 
  	peptide_hits_(source.peptide_hits_), 
  	peptide_significance_threshold_(source.peptide_significance_threshold_) 
  {
  }
  
  Identification::~Identification() 
  {
  	peptide_hits_.clear();
  	protein_hits_.clear();
  }
  
  // read access to peptide hits
  const std::vector<PeptideHit>& Identification::getPeptideHits() const 
 	{ 
 		return peptide_hits_;
 	}

  // mutable access to peptide hits
  std::vector<PeptideHit>& Identification::getPeptideHits() 
 	{ 
 		return peptide_hits_;
 	}

	// retrival of the peptide significance threshold value
  float Identification::getPeptideSignificanceThreshold() const 
  { 
  	return peptide_significance_threshold_;
  }

	// setting of the peptide significance threshold value
	void Identification::setPeptideSignificanceThreshold(float value) 
	{ 
		peptide_significance_threshold_ = value;
	}

  void Identification::clear()
  {
   date_.clear();
   peptide_significance_threshold_ = 0;
   protein_significance_threshold_ = 0;
   peptide_hits_.clear();
   protein_hits_.clear();
  }
  
  Identification& Identification::operator=(const Identification& source) 
  {
  	if (this == &source)
  	{
  		return *this;		
  	}
    //PersistentObject::operator=(source);
    date_ = source.date_;
    peptide_hits_ = source.peptide_hits_;
    protein_hits_ = source.protein_hits_;
    peptide_significance_threshold_ = source.peptide_significance_threshold_;
    protein_significance_threshold_ = source.protein_significance_threshold_;
    return *this;  
  }

	// Equality operator
	bool Identification::operator == (const Identification& rhs) const
	{
		return date_ == rhs.getDateTime() 
						&& peptide_hits_ == rhs.getPeptideHits()
						&& protein_hits_ == rhs.getProteinHits()
						&& peptide_significance_threshold_ == rhs.getPeptideSignificanceThreshold()
						&& protein_significance_threshold_ == rhs.getProteinSignificanceThreshold();						 
	}
		
	// Inequality operator
	bool Identification::operator != (const Identification& rhs) const
	{
		return !(*this == rhs);						 
	}
	
	// inserts a peptide hit if the peptide hit has the same score type as the existing hits
  void Identification::insertPeptideHit(const PeptideHit& input)
  {
    if ( !peptide_hits_.size() )
    {  
      peptide_hits_.push_back(input);
    }
    else if ( peptide_hits_.begin()->getScoreType() == input.getScoreType() )
    {
      peptide_hits_.push_back(input);
    }
    else
    {
      stringstream ss;
      ss << peptide_hits_.begin()->getScoreType() << " != " <<  input.getScoreType();
      throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Incompatible PeptideHit.score_type",ss.str().c_str());
    }
  }

	// Stores peptide and protein hits and clears all previous hit information
  void Identification::setPeptideAndProteinHits(const std::vector<PeptideHit>& peptide_hits, const std::vector<ProteinHit>& protein_hits)
  {
  	peptide_hits_.clear();
  	protein_hits_.clear();
  	peptide_hits_ = peptide_hits;
  	protein_hits_ = protein_hits;
  }
  
  void Identification::assignRanks()
  {
    UnsignedInt rank = 1;
    sort();
    for ( vector<PeptideHit>::iterator lit = peptide_hits_.begin(); lit != peptide_hits_.end(); ++lit )
    {
      lit->setRank(rank++);
    }
    rank = 1;
    for ( vector<ProteinHit>::iterator lit = protein_hits_.begin(); lit != protein_hits_.end(); ++lit )
    {
      lit->setRank(rank++);
    }
  }
    
  void Identification::sort()
  {
  	ProteinIdentification::sort();		
		std::sort(peptide_hits_.begin(), peptide_hits_.end(), ScoreMore());
  }
  
	bool Identification::empty() const
	{
		DateTime temp_date;
		
		return (peptide_significance_threshold_ == 0
						&& protein_significance_threshold_ == 0
						&& protein_hits_.size() == 0
						&& peptide_hits_.size() == 0
						&& date_ == temp_date);		
	}
	
  vector<PeptideHit>* Identification::getReferencingHits(String date_time, String accession) const
  {
  	vector<PeptideHit>* found_hits = new vector<PeptideHit>();
  	
  	// for every peptide hit
		for(UnsignedInt i = 0; i < peptide_hits_.size(); i++)
		{
			const vector< pair<String, String> >& references = peptide_hits_[i].getProteinIndices();
			//for every reference of the peptide hit
			for(UnsignedInt j = 0; j < references.size(); j++)
			{
				if (references[j].first == date_time && references[j].second == accession)
				{
					found_hits->push_back(peptide_hits_[i]);
				}
			}
		}
		return found_hits;
  	
  }

	vector<PeptideHit>* Identification::getNonReferencingHits(const multimap< String, ProteinHit >& protein_hits) const
  {
  	vector<PeptideHit>* 								found_hits;
  	vector<PeptideHit>* 								all_found_hits 		= new vector<PeptideHit>();
  	String 															temp_key 					= "";
		vector<ProteinHit> 									temp_hits;

		for(multimap< String, ProteinHit >::const_iterator hits_it = protein_hits.begin();
				hits_it != protein_hits.end();
				hits_it++)
		{
			if (temp_key == "")
			{
				temp_key = hits_it->first;
			}
			if (temp_key == hits_it->first)
			{ 
				temp_hits.push_back(hits_it->second);
			}
			else
			{
				found_hits = getNonReferencingHits(temp_hits.begin(), 
																					 temp_hits.end(), 
																					 temp_key);
	  		all_found_hits->insert(all_found_hits->end(), found_hits->begin(), found_hits->end());
				delete found_hits;  							
	
				temp_hits.clear();
				temp_hits.push_back(hits_it->second);
				temp_key = hits_it->first;																					 				
			}
		}
		found_hits = getNonReferencingHits(temp_hits.begin(), 
																			 temp_hits.end(), 
																			 temp_key);
		all_found_hits->insert(all_found_hits->end(), found_hits->begin(), found_hits->end());
		delete found_hits;  							

  	return all_found_hits;
	}

  void Identification::insertProteinHit(const ProteinHit& input)
  {
    if ( protein_hits_.size() == 0)
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


}// namespace OpenMS
