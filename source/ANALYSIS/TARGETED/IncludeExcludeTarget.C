// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TARGETED/IncludeExcludeTarget.h>

#include <algorithm>

using namespace std;


namespace OpenMS
{
	IncludeExcludeTarget::IncludeExcludeTarget()
		:	CVTermList(),
			precursor_mz_(numeric_limits<DoubleReal>::max()),
			product_mz_(numeric_limits<DoubleReal>::max())
	{
	}

  IncludeExcludeTarget::IncludeExcludeTarget(const IncludeExcludeTarget& rhs)
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

	IncludeExcludeTarget::~IncludeExcludeTarget()
	{
	}

	IncludeExcludeTarget& IncludeExcludeTarget::operator = (const IncludeExcludeTarget& rhs)
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

	bool IncludeExcludeTarget::operator == (const IncludeExcludeTarget& rhs) const
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

	bool IncludeExcludeTarget::operator != (const IncludeExcludeTarget& rhs) const
	{
		return !(*this == rhs);
	}

	void IncludeExcludeTarget::setName(const String& name)
	{
		name_ = name;
	}

	const String& IncludeExcludeTarget::getName() const
	{
		return name_;
	}

	void IncludeExcludeTarget::setPeptideRef(const String& peptide_ref)
	{
		peptide_ref_ = peptide_ref;
	}

	const String& IncludeExcludeTarget::getPeptideRef() const
	{
		return peptide_ref_;
	}

	void IncludeExcludeTarget::setCompoundRef(const String& compound_ref)
	{
		compound_ref_ = compound_ref;
	}

	const String& IncludeExcludeTarget::getCompoundRef() const
	{
		return compound_ref_;
	}

	void IncludeExcludeTarget::setPrecursorMZ(DoubleReal mz)
	{
		precursor_mz_ = mz;
	}

	DoubleReal IncludeExcludeTarget::getPrecursorMZ() const
	{
		return precursor_mz_;
	}

	void IncludeExcludeTarget::setPrecursorCVTermList(const CVTermList& list)
	{
		precursor_cv_terms_ = list;
	}

	void IncludeExcludeTarget::addPrecursorCVTerm(const CVTerm& cv_term)
	{
		precursor_cv_terms_.addCVTerm(cv_term);
	}

	const CVTermList& IncludeExcludeTarget::getPrecursorCVTermList() const
	{
		return precursor_cv_terms_;
	}


  void IncludeExcludeTarget::setProductMZ(DoubleReal mz)
  {
    product_mz_ = mz;
  }

  DoubleReal IncludeExcludeTarget::getProductMZ() const
  {
    return product_mz_;
  }

  void IncludeExcludeTarget::setProductCVTermList(const CVTermList& list)
  {
    product_cv_terms_ = list;
  }

	void IncludeExcludeTarget::addProductCVTerm(const CVTerm& cv_term)
	{
		product_cv_terms_.addCVTerm(cv_term);
	}

  const CVTermList& IncludeExcludeTarget::getProductCVTermList() const
  {
    return product_cv_terms_;
  }

	void IncludeExcludeTarget::setInterpretations(const vector<CVTermList>& interpretations)
	{
		interpretation_list_ = interpretations;
	}

	const vector<CVTermList>& IncludeExcludeTarget::getInterpretations() const
	{
		return interpretation_list_;
	}

	void IncludeExcludeTarget::addInterpretation(const CVTermList& interpretation)
	{
		interpretation_list_.push_back(interpretation);
	}

	void IncludeExcludeTarget::setConfigurations(const vector<Configuration>& configurations)
	{
		configurations_ = configurations;
	}

	const vector<IncludeExcludeTarget::Configuration>& IncludeExcludeTarget::getConfigurations() const
	{
		return configurations_;
	}

	void IncludeExcludeTarget::addConfiguration(const Configuration& configuration)
	{
		configurations_.push_back(configuration);
	}

	void IncludeExcludeTarget::setPrediction(const CVTermList& prediction)
	{
		prediction_ = prediction;
	}

	const CVTermList& IncludeExcludeTarget::getPrediction() const
	{
		return prediction_;
	}

	void IncludeExcludeTarget::addPredictionTerm(const CVTerm& term)
	{
		prediction_.addCVTerm(term);
	}
	

  void IncludeExcludeTarget::updateMembers_()
	{
	}
}


