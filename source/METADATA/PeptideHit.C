// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/PeptideHit.h>

#include <algorithm>

using namespace std;

namespace OpenMS 
{
	// default constructor
  PeptideHit::PeptideHit()
		:	MetaInfoInterface(),
			score_(0), 
			rank_(0), 
			charge_(0),
			aa_before_(' '),
			aa_after_(' ')
  {
  }
  
	// values constructor
  PeptideHit::PeptideHit(DoubleReal score, UInt rank, Int charge, const AASequence& sequence)
    	: MetaInfoInterface(),
    		score_(score), 
    		rank_(rank),
				charge_(charge),
    		sequence_(sequence),
				aa_before_(' '),
				aa_after_(' ')
  {
  }
  
	// copy constructor
  PeptideHit::PeptideHit(const PeptideHit& source)
		:	MetaInfoInterface(source),
			score_(source.score_), 
			rank_(source.rank_),
			charge_(source.charge_),
			sequence_(source.sequence_),
			aa_before_(source.aa_before_),
			aa_after_(source.aa_after_),
			corresponding_protein_accessions_(source.corresponding_protein_accessions_)
  {
  }
  
	// destructor
  PeptideHit::~PeptideHit()
  {  	
  }
   
  PeptideHit& PeptideHit::operator= (const PeptideHit& source)
  {
	 	if (this == &source)
  	{
  		return *this;
  	}
  	
  	MetaInfoInterface::operator=(source);
    score_ = source.score_;
		charge_ = source.charge_;
		rank_  = source.rank_;
    sequence_ = source.sequence_;
    aa_before_ = source.aa_before_;
    aa_after_ = source.aa_after_;
    corresponding_protein_accessions_ = source.corresponding_protein_accessions_;

    return *this;
  }
	
	bool PeptideHit::operator == (const PeptideHit& rhs) const	
	{
		return MetaInfoInterface::operator==(rhs)
			&& score_ == rhs.score_ 
			&& rank_ == rhs.rank_ 
			&& charge_ == rhs.charge_
			&& sequence_ == rhs.sequence_
			&& aa_before_ == rhs.aa_before_
			&& aa_after_ == rhs.aa_after_
			&& corresponding_protein_accessions_ == rhs.corresponding_protein_accessions_;
	}

	bool PeptideHit::operator != (const PeptideHit& rhs) const	
	{
		return !operator==(rhs);
	}

	void PeptideHit::addProteinAccession(const String& accession) 
	{
		if (find(corresponding_protein_accessions_.begin(), corresponding_protein_accessions_.end(), accession) == corresponding_protein_accessions_.end())
		{
			corresponding_protein_accessions_.push_back(accession);
		}
	}

  // returns the score of the peptide hit 
  DoubleReal PeptideHit::getScore() const 
  {
  	return score_;
  }
  
	// returns the rank of the peptide hit
  UInt PeptideHit::getRank() const 
  {
  	return rank_;
  }
	
	// returns the peptide sequence without trailing or following spaces
	const AASequence& PeptideHit::getSequence() const 
	{
		return sequence_;
	}

	Int PeptideHit::getCharge() const
	{
		return charge_;
	}
	
	void PeptideHit::setSequence(const AASequence& sequence) 
	{
		sequence_ = sequence; 
	}

	void PeptideHit::setCharge(Int charge)
	{
		charge_ = charge;
	}
	
	// returns the corresponding protein accessions
	const vector<String>& PeptideHit::getProteinAccessions() const 
	{
		return corresponding_protein_accessions_;
	}

	void PeptideHit::setProteinAccessions(const vector<String>& accessions)
	{
		corresponding_protein_accessions_ = accessions;
	}
	
  // sets the score of the peptide hit 
  void PeptideHit::setScore(DoubleReal score) 
  {
  	score_ = score;
  }
  
	// sets the rank
  void PeptideHit::setRank(UInt newrank) 
  {
  	rank_ = newrank;
  }

	void PeptideHit::setAABefore(char acid)
	{
		aa_before_ = acid;
	}
	
	char PeptideHit::getAABefore() const
	{
		return aa_before_;
	}
	
	void PeptideHit::setAAAfter(char acid)
	{
		aa_after_ = acid;
	}
	
	char PeptideHit::getAAAfter() const
	{
		return aa_after_;
	}


} // namespace OpenMS
