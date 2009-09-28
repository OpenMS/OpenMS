// -*- Mode: C++; tab-width: 2; -*-
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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>

#include <algorithm>

using namespace std;


namespace OpenMS
{
	ReactionMonitoringTransition::ReactionMonitoringTransition()
		:	MetaInfoInterface(),
			CVTermList(),
			precursor_mz_(numeric_limits<DoubleReal>::max()),
			precursor_charge_(numeric_limits<Int>::max()),
			product_mz_(numeric_limits<DoubleReal>::max()),
			product_charge_(numeric_limits<Int>::max())
	{
	}

  ReactionMonitoringTransition::ReactionMonitoringTransition(const ReactionMonitoringTransition& rhs)
		:	MetaInfoInterface(rhs),
			CVTermList(rhs),
			name_(rhs.name_),
			precursor_mz_(rhs.precursor_mz_),
      precursor_charge_(rhs.precursor_charge_),
      product_mz_(rhs.product_mz_),
      product_charge_(rhs.product_charge_),
      interpretation_list_(rhs.interpretation_list_),
			peptide_ref_(rhs.peptide_ref_),
			compound_ref_(rhs.compound_ref_),
			configurations_(rhs.configurations_)
	{
	}

	ReactionMonitoringTransition::~ReactionMonitoringTransition()
	{
	}

	ReactionMonitoringTransition& ReactionMonitoringTransition::operator = (const ReactionMonitoringTransition& rhs)
	{
		if (&rhs != this)
		{
			MetaInfoInterface::operator = (rhs);
			CVTermList::operator = (rhs);
			name_ = rhs.name_;
			precursor_mz_ = rhs.precursor_mz_;
			precursor_charge_ = rhs.precursor_charge_;
			product_mz_ = rhs.product_mz_;
			product_charge_ = rhs.product_charge_;
			interpretation_list_ = rhs.interpretation_list_;
			peptide_ref_ = rhs.peptide_ref_;
			compound_ref_ = rhs.compound_ref_;
			configurations_ = rhs.configurations_;
		}
		return *this;
	}

	bool ReactionMonitoringTransition::operator == (const ReactionMonitoringTransition& rhs) const
	{
		return  MetaInfoInterface::operator == (rhs) &&
			      CVTermList::operator == (rhs) &&
			      name_ == rhs.name_ &&
			      precursor_mz_ == rhs.precursor_mz_ &&
			      precursor_charge_ == rhs.precursor_charge_ &&
			      product_mz_ == rhs.product_mz_ &&
			      product_charge_ == rhs.product_charge_ &&
			      interpretation_list_ == rhs.interpretation_list_ &&
			      peptide_ref_ == rhs.peptide_ref_ &&
			      compound_ref_ == rhs.compound_ref_ &&
			      configurations_ == rhs.configurations_;
	}

	bool ReactionMonitoringTransition::operator != (const ReactionMonitoringTransition& rhs) const
	{
		return !(*this == rhs);
	}

	void ReactionMonitoringTransition::setName(const String& name)
	{
		name_ = name;
	}

	const String& ReactionMonitoringTransition::getName() const
	{
		return name_;
	}

	void ReactionMonitoringTransition::setPeptideRef(const String& peptide_ref)
	{
		peptide_ref_ = peptide_ref;
	}

	const String& ReactionMonitoringTransition::getPeptideRef() const
	{
		return peptide_ref_;
	}

	void ReactionMonitoringTransition::setCompoundRef(const String& compound_ref)
	{
		compound_ref_ = compound_ref;
	}

	const String& ReactionMonitoringTransition::getCompoundRef() const
	{
		return compound_ref_;
	}

	void ReactionMonitoringTransition::setPrecursorMZ(DoubleReal mz)
	{
		precursor_mz_ = mz;
	}

	DoubleReal ReactionMonitoringTransition::getPrecursorMZ() const
	{
		return precursor_mz_;
	}

	void ReactionMonitoringTransition::setPrecursorCharge(Int charge)
	{
		precursor_charge_ = charge;
	}

	Int ReactionMonitoringTransition::getPrecursorCharge() const
	{
		return precursor_charge_;
	}

  void ReactionMonitoringTransition::setProductMZ(DoubleReal mz)
  {
    product_mz_ = mz;
  }

  DoubleReal ReactionMonitoringTransition::getProductMZ() const
  {
    return product_mz_;
  }

  void ReactionMonitoringTransition::setProductCharge(Int charge)
  {
    product_charge_ = charge;
  }

  Int ReactionMonitoringTransition::getProductCharge() const
  {
    return product_charge_;
  }

	void ReactionMonitoringTransition::setInterpretations(const vector<TransitionInterpretation>& interpretations)
	{
		interpretation_list_ = interpretations;
	}

	const vector<TransitionInterpretation>& ReactionMonitoringTransition::getInterpretations() const
	{
		return interpretation_list_;
	}

	void ReactionMonitoringTransition::addInterpretation(const TransitionInterpretation& interpretation)
	{
		interpretation_list_.push_back(interpretation);
	}

	void ReactionMonitoringTransition::setConfigurations(const vector<Configuration>& configurations)
	{
		configurations_ = configurations;
	}

	const vector<ReactionMonitoringTransition::Configuration>& ReactionMonitoringTransition::getConfigurations() const
	{
		return configurations_;
	}

	void ReactionMonitoringTransition::addConfiguration(const Configuration& configuration)
	{
		configurations_.push_back(configuration);
	}

  void ReactionMonitoringTransition::updateMembers_()
	{
	}
}


