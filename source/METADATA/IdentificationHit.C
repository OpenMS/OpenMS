// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
