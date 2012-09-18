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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/IdentificationHit.h>

using namespace std;

namespace OpenMS {

  IdentificationHit::IdentificationHit()
    : MetaInfoInterface(),
    	id_(),
			charge_(0),
			calculated_mass_to_charge_(0.0),
			name_(""),
			pass_threshold_(true),
			rank_(0)
  {
  }

  IdentificationHit::IdentificationHit(const IdentificationHit& rhs)
  	: MetaInfoInterface(rhs),
  		id_(rhs.id_),
			charge_(rhs.charge_),
			calculated_mass_to_charge_(rhs.calculated_mass_to_charge_),
			experimental_mass_to_charge_(rhs.experimental_mass_to_charge_),
			name_(rhs.name_),
			pass_threshold_(rhs.pass_threshold_),
			rank_(rhs.rank_)
  {
  }

  IdentificationHit::~IdentificationHit()
  {
  }

  IdentificationHit& IdentificationHit::operator=(const IdentificationHit& rhs)
  {
  	if (this == &rhs)
  	{
  		return *this;
  	}

    MetaInfoInterface::operator=(rhs);
    id_ = rhs.id_;
		charge_ = rhs.charge_;
		calculated_mass_to_charge_ = rhs.calculated_mass_to_charge_;
		experimental_mass_to_charge_ = rhs.experimental_mass_to_charge_;
		name_ = rhs.name_;
		pass_threshold_ = rhs.pass_threshold_;
		rank_ = rhs.rank_;

    return *this;
  }

	// Equality operator
	bool IdentificationHit::operator == (const IdentificationHit& rhs) const
	{
		return MetaInfoInterface::operator==(rhs)
				&& id_ == rhs.id_
				&& charge_ == rhs.charge_
				&& calculated_mass_to_charge_ == rhs.calculated_mass_to_charge_
				&& experimental_mass_to_charge_ == rhs.experimental_mass_to_charge_
				&& name_ == rhs.name_
				&& pass_threshold_ == rhs.pass_threshold_
				&& rank_ == rhs.rank_
				;
	}

	// Inequality operator
	bool IdentificationHit::operator != (const IdentificationHit& rhs) const
	{
		return !(*this == rhs);
	}


  void IdentificationHit::setId(const String& id)
	{
		id_ = id;
	}

  const String& IdentificationHit::getId() const
	{
		return id_;
	}

  void IdentificationHit::setCharge(Int charge)
	{
		charge_ = charge;
	}

  Int IdentificationHit::getCharge() const
	{
		return charge_;
	}

	void IdentificationHit::setCalculatedMassToCharge(DoubleReal mz)
	{
		calculated_mass_to_charge_ = mz;
	}

  DoubleReal IdentificationHit::getCalculatedMassToCharge() const
	{
		return calculated_mass_to_charge_;
	}

	void IdentificationHit::setExperimentalMassToCharge(DoubleReal mz)
	{
		experimental_mass_to_charge_ = mz;
	}

	DoubleReal IdentificationHit::getExperimentalMassToCharge() const
	{
		return experimental_mass_to_charge_;
	}

  void IdentificationHit::setName(const String& name)
	{
		name_ = name;
	}

  const String& IdentificationHit::getName() const
	{
		return name_;
	}

  void IdentificationHit::setPassThreshold(bool pass)
	{
		pass_threshold_ = pass;
	}

  bool IdentificationHit::getPassThreshold() const
	{
		return pass_threshold_;
	}

  void IdentificationHit::setRank(Int rank)
	{
		rank_ = rank;
	}
  
	Int IdentificationHit::getRank() const
	{
		return rank_;
	}


}// namespace OpenMS
