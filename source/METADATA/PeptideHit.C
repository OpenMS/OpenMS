// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
