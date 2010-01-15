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
		:	CVTermList(),
			precursor_mz_(numeric_limits<DoubleReal>::max()),
			product_mz_(numeric_limits<DoubleReal>::max())
	{
	}

  ReactionMonitoringTransition::ReactionMonitoringTransition(const ReactionMonitoringTransition& rhs)
		: CVTermList(rhs),
			name_(rhs.name_),
			precursor_mz_(rhs.precursor_mz_),
      precursor_cv_terms_(rhs.precursor_cv_terms_),
      product_mz_(rhs.product_mz_),
      product_cv_terms_(rhs.product_cv_terms_),
      interpretation_list_(rhs.interpretation_list_),
			peptide_ref_(rhs.peptide_ref_),
			compound_ref_(rhs.compound_ref_),
			configurations_(rhs.configurations_),
			prediction_(rhs.prediction_)
	{
	}

	ReactionMonitoringTransition::~ReactionMonitoringTransition()
	{
	}

	ReactionMonitoringTransition& ReactionMonitoringTransition::operator = (const ReactionMonitoringTransition& rhs)
	{
		if (&rhs != this)
		{
			CVTermList::operator = (rhs);
			name_ = rhs.name_;
			precursor_mz_ = rhs.precursor_mz_;
			precursor_cv_terms_ = rhs.precursor_cv_terms_;
			product_mz_ = rhs.product_mz_;
			product_cv_terms_ = rhs.product_cv_terms_;
			interpretation_list_ = rhs.interpretation_list_;
			peptide_ref_ = rhs.peptide_ref_;
			compound_ref_ = rhs.compound_ref_;
			configurations_ = rhs.configurations_;
			prediction_ = rhs.prediction_;
		}
		return *this;
	}

	bool ReactionMonitoringTransition::operator == (const ReactionMonitoringTransition& rhs) const
	{
		return  CVTermList::operator == (rhs) &&
			      name_ == rhs.name_ &&
			      precursor_mz_ == rhs.precursor_mz_ &&
			      precursor_cv_terms_ == rhs.precursor_cv_terms_ &&
			      product_mz_ == rhs.product_mz_ &&
			      product_cv_terms_ == rhs.product_cv_terms_ &&
			      interpretation_list_ == rhs.interpretation_list_ &&
			      peptide_ref_ == rhs.peptide_ref_ &&
			      compound_ref_ == rhs.compound_ref_ &&
			      configurations_ == rhs.configurations_ &&
						prediction_ == rhs.prediction_;
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

	void ReactionMonitoringTransition::setPrecursorCVTermList(const CVTermList& list)
	{
		precursor_cv_terms_ = list;
	}

	void ReactionMonitoringTransition::addPrecursorCVTerm(const CVTerm& cv_term)
	{
		precursor_cv_terms_.addCVTerm(cv_term);
	}

	const CVTermList& ReactionMonitoringTransition::getPrecursorCVTermList() const
	{
		return precursor_cv_terms_;
	}


  void ReactionMonitoringTransition::setProductMZ(DoubleReal mz)
  {
    product_mz_ = mz;
  }

  DoubleReal ReactionMonitoringTransition::getProductMZ() const
  {
    return product_mz_;
  }

  void ReactionMonitoringTransition::setProductCVTermList(const CVTermList& list)
  {
    product_cv_terms_ = list;
  }

	void ReactionMonitoringTransition::addProductCVTerm(const CVTerm& cv_term)
	{
		product_cv_terms_.addCVTerm(cv_term);
	}

  const CVTermList& ReactionMonitoringTransition::getProductCVTermList() const
  {
    return product_cv_terms_;
  }

	void ReactionMonitoringTransition::setInterpretations(const vector<CVTermList>& interpretations)
	{
		interpretation_list_ = interpretations;
	}

	const vector<CVTermList>& ReactionMonitoringTransition::getInterpretations() const
	{
		return interpretation_list_;
	}

	void ReactionMonitoringTransition::addInterpretation(const CVTermList& interpretation)
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

	void ReactionMonitoringTransition::setPrediction(const CVTermList& prediction)
	{
		prediction_ = prediction;
	}

	const CVTermList& ReactionMonitoringTransition::getPrediction() const
	{
		return prediction_;
	}

	void ReactionMonitoringTransition::addPredictionTerm(const CVTerm& term)
	{
		prediction_.addCVTerm(term);
	}
	

  void ReactionMonitoringTransition::updateMembers_()
	{
	}
}


