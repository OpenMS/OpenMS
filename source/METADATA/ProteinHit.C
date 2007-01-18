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

#include <OpenMS/METADATA/ProteinHit.h>

using namespace std;

namespace OpenMS 
{
	// default constructor
  ProteinHit::ProteinHit()
			:	//PersistentObject(),
				score_(0), 
				score_type_(""), 
				rank_(0), 
				accession_(""), 
				accession_type_(""), 
				sequence_("")
  {
  }
  
	// values constructor
  ProteinHit::ProteinHit(DoubleReal score, 
  											 std::string score_type, 
  											 UnsignedInt rank, 
  											 String accession, 
  											 std::string accession_type, 
  											 String sequence)
    	: //PersistentObject(),
    		score_(score), 
    		score_type_(score_type), 
    		rank_(rank), 
    		accession_(accession), 
    		accession_type_(accession_type),
    		sequence_(sequence)
  {
  	accession_.trim();
  	sequence_.trim();
  }
  
	// copy constructor
  ProteinHit::ProteinHit(const ProteinHit& source)
    	: //PersistentObject(source),
				score_(source.score_), 
				score_type_(source.score_type_), 
				rank_(source.rank_), 
				accession_(source.accession_), 
				accession_type_(source.accession_type_),
				sequence_(source.sequence_)
  {
  }
  
	// destructor
  ProteinHit::~ProteinHit()
  {  	
  	score_type_.erase();
  	sequence_.erase();
  	accession_.erase();
  	accession_type_.erase();
  }
  
  // clears all data of the protein hit
  void ProteinHit::clear() {
    //clearId();
    score_ = 0;
    score_type_.erase();
    sequence_.erase();
    accession_.erase();
    accession_type_.erase();
    rank_ = 0;
  }
   
  // assignment operator
  ProteinHit& ProteinHit::operator= (const ProteinHit& source)
  {
  	if (this == &source)
  	{
  		return *this;
  	}  			
    //PersistentObject::operator= (source);
    score_ = source.score_;
    score_type_ = source.score_type_;
		rank_  = source.rank_;
    sequence_ = source.sequence_;
    accession_ = source.accession_;
    accession_type_ = source.accession_type_;
    return *this;
  }
	
	// equality operator
	bool ProteinHit::operator == (const ProteinHit& rhs) const	
	{
		return score_ == rhs.score_ 
			&& score_type_ == rhs.score_type_ 
			&& rank_ == rhs.rank_ 
			&& accession_ == rhs.accession_ 
			&& accession_type_ == rhs.accession_type_ 
			&& sequence_ == rhs.sequence_;
	}

	// inequality operator
	bool ProteinHit::operator != (const ProteinHit& rhs) const	
	{
		return score_ != rhs.score_ 
			|| score_type_ != rhs.score_type_ 
			|| rank_ != rhs.rank_ 
			|| accession_ != rhs.accession_ 
			|| accession_type_ != rhs.accession_type_ 
			|| sequence_ != rhs.sequence_;
	}
	
  // returns the score of the protein hit 
  Real ProteinHit::getScore() const 
  {
  	return score_;
  }
  
  // returns the type of the score
  const std::string& ProteinHit::getScoreType() const 
  {
  	return score_type_;
  }
	
	// returns the rank of the protein hit
  UnsignedInt ProteinHit::getRank() const 
  {
  	return rank_;
  }
	
	// returns the protein sequence
	const String& ProteinHit::getSequence() const 
	{
		return sequence_;
	}
	
	// returns the accession of the protein
	const String& ProteinHit::getAccession() const 
	{
		return accession_;
	}
	
	// returns the type of the accession string of the protein
	const std::string& ProteinHit::getAccessionType() const 
	{
		return accession_type_;
	}    	
	
  // sets the score of the protein hit 
  void ProteinHit::setScore(const DoubleReal& score) 
  {
  	score_ = score;
  }

  // sets the type of the score
  void ProteinHit::setScoreType(const std::string& score_type) 
  {
  	score_type_ = score_type;
  }
  
	// sets the rank
  void ProteinHit::setRank(UnsignedInt newrank) 
  {
  	rank_ = newrank;
  }
	
	// sets the protein sequence
	void ProteinHit::setSequence(const String& sequence) 
	{
		sequence_ = sequence;
		sequence_.trim();
	}
	
	// sets the accession of the protein
	void ProteinHit::setAccession(const String& accession) 
	{
		accession_ = accession;
		accession_.trim();
	}
	
	// sets the type of the accession string of the protein
	void ProteinHit::setAccessionType(const std::string& accession_type) 
	{
		accession_type_ = accession_type;
	}    	

//	void ProteinHit::deserialize_(const string& attribute_name, 
//																const string& value)
//  {
//		//PersistentObject::deserialize_(attribute_name, value);
//		//todo attribute setzen
//  }

} // namespace OpenMS
