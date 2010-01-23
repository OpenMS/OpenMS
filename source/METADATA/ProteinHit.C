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

#include <OpenMS/METADATA/ProteinHit.h>

using namespace std;

namespace OpenMS 
{
	// default constructor
  ProteinHit::ProteinHit()
		:	MetaInfoInterface(),
			score_(0), 
			rank_(0), 
			accession_(""), 
			sequence_(""),
			coverage_(0)
  {
  }
  
	// values constructor
  ProteinHit::ProteinHit(DoubleReal score, UInt rank, String accession, String sequence)
    :	MetaInfoInterface(),
  		score_(score), 
  		rank_(rank), 
  		accession_(accession.trim()), 
  		sequence_(sequence.trim()),
  		coverage_(0)
  {
  }
  
	// copy constructor
  ProteinHit::ProteinHit(const ProteinHit& source)
		:	MetaInfoInterface(source),
			score_(source.score_), 
			rank_(source.rank_), 
			accession_(source.accession_), 
			sequence_(source.sequence_),
			coverage_(source.coverage_)
  {
  }
  
	// destructor
  ProteinHit::~ProteinHit()
  {  	
  }
   
  // assignment operator
  ProteinHit& ProteinHit::operator= (const ProteinHit& source)
  {
  	if (this == &source)
  	{
  		return *this;
  	}
  	
    MetaInfoInterface::operator=(source);
    score_ = source.score_;
		rank_  = source.rank_;
    sequence_ = source.sequence_;
    accession_ = source.accession_;
    coverage_ = source.coverage_;
    
    return *this;
  }
	
	// equality operator
	bool ProteinHit::operator == (const ProteinHit& rhs) const	
	{
		return MetaInfoInterface::operator==(rhs)
			&& score_ == rhs.score_ 
			&& rank_ == rhs.rank_
			&& accession_ == rhs.accession_ 
			&& sequence_ == rhs.sequence_
			&& coverage_ == rhs.coverage_;
	}

	// inequality operator
	bool ProteinHit::operator != (const ProteinHit& rhs) const	
	{
		return ! operator==(rhs);
	}
	
  // returns the score of the protein hit 
  Real ProteinHit::getScore() const 
  {
  	return score_;
  }
  
	// returns the rank of the protein hit
  UInt ProteinHit::getRank() const 
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
	
	// returns the coverage (in percent) of the protein hit based upon matched peptides
	DoubleReal ProteinHit::getCoverage() const
	{
		return coverage_;
	}
	
  // sets the score of the protein hit 
  void ProteinHit::setScore(const DoubleReal score) 
  {
  	score_ = score;
  }

	// sets the rank
  void ProteinHit::setRank(UInt newrank) 
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

	// sets the coverage (in percent) of the protein hit based upon matched peptides
	void ProteinHit::setCoverage(const DoubleReal coverage)
	{
		coverage_ = coverage;
	}

} // namespace OpenMS
