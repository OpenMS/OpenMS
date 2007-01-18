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

#include <OpenMS/METADATA/PeptideHit.h>

#include <algorithm>

using namespace std;

namespace OpenMS 
{
	// default constructor
  PeptideHit::PeptideHit()
			:	score_(0), 
				score_type_(""), 
				rank_(0), 
				charge_(0), 
				sequence_("")
  {
  }
  
	// values constructor
  PeptideHit::PeptideHit(double score, 
  											 std::string score_type, 
  											 uint rank, 
												 SignedInt charge,
  											 String sequence)
    	: score_(score), 
    		score_type_(score_type), 
    		rank_(rank),
				charge_(charge),
    		sequence_(sequence)
  {
  	sequence_.trim();
  }
  
	// copy constructor
  PeptideHit::PeptideHit(const PeptideHit& source)
    	: //PersistentObject(source),
				score_(source.score_), 
				score_type_(source.score_type_), 
				rank_(source.rank_),
				charge_(source.charge_),
				sequence_(source.sequence_),
				corresponding_protein_indices_(source.corresponding_protein_indices_)
  {
  }
  
	// destructor
  PeptideHit::~PeptideHit()
  {  	
  	score_type_.erase();
  	sequence_.erase();
  }
  
  void PeptideHit::clear() {
    //clearId();
    score_ = 0;
    score_type_.erase();
    sequence_.erase();
    rank_ = 0;
		charge_ = 0;
    corresponding_protein_indices_.clear();
  }
   
  PeptideHit& PeptideHit::operator= (const PeptideHit& source)
  {
	 	if (this == &source)
  	{
  		return *this;
  	}  			
    //PersistentObject::operator= (source);
    score_ = source.score_;
		charge_ = source.charge_;
    score_type_ = source.score_type_;
		rank_  = source.rank_;
    sequence_ = source.sequence_;
    corresponding_protein_indices_ = source.corresponding_protein_indices_;
    return *this;
  }
	
	bool PeptideHit::operator == (const PeptideHit& rhs) const	
	{
		return score_ == rhs.score_ 
			&& score_type_ == rhs.score_type_ 
			&& rank_ == rhs.rank_ 
			&& charge_ == rhs.charge_
			&& sequence_ == rhs.sequence_
			&& corresponding_protein_indices_ == rhs.corresponding_protein_indices_;
	}

	bool PeptideHit::operator != (const PeptideHit& rhs) const	
	{
		return !(*this == rhs);
	}
	
	void PeptideHit::addProteinIndex(const pair<String, String>& index) 
	{ 
		bool found = false;
		
		for(vector< pair<String, String> >::iterator it = corresponding_protein_indices_.begin(); 
				it != corresponding_protein_indices_.end();
				it++)
		{
			if (it->first == index.first && it->second == index.second)
			{
				found = true;
			}
		}
		if (!found)
		{
			corresponding_protein_indices_.push_back(index);		
		}
	}

	void PeptideHit::addProteinIndex(const DateTime& date, const String& accession) 
	{ 
		String date_time = "";
		
		date.get(date_time);
		addProteinIndex(make_pair(date_time, accession));		
	}

  // returns the score of the peptide hit 
  float PeptideHit::getScore() const 
  {
  	return score_;
  }
  
  // returns the type of the score
  const std::string& PeptideHit::getScoreType() const 
  {
  	return score_type_;
  }
	
	// returns the rank of the peptide hit
  UnsignedInt PeptideHit::getRank() const 
  {
  	return rank_;
  }
	
	// returns the peptide sequence without trailing or following spaces
	String PeptideHit::getSequence() const 
	{
		return sequence_;
	}

	SignedInt PeptideHit::getCharge() const
	{
		return charge_;
	}
	
	void PeptideHit::setSequence(const String& sequence) 
	{
		sequence_ = sequence; 
		sequence_.trim();
	}

	void PeptideHit::setCharge(SignedInt charge)
	{
		charge_ = charge;
	}
	
	// returns the corresponding protein indices
	const vector< pair<String, String> >& PeptideHit::getProteinIndices() const 
	{
		return corresponding_protein_indices_;
	}

	// returns the corresponding protein indices
	vector< pair<String, String> >& PeptideHit::getProteinIndices()
	{
		return corresponding_protein_indices_;
	}

  // sets the score of the peptide hit 
  void PeptideHit::setScore(const double& score) 
  {
  	score_ = score;
  }
  
  // sets the type of the score
  void PeptideHit::setScoreType(const std::string& score_type) 
  {
  	score_type_ = score_type;
  }

	// sets the rank
  void PeptideHit::setRank(uint newrank) 
  {
  	rank_ = newrank;
  }
  
  void PeptideHit::setProteinIndices(const std::vector< std::pair<String, String> >& indices)
  {
  	corresponding_protein_indices_ = indices;
  }
  
} // namespace OpenMS
