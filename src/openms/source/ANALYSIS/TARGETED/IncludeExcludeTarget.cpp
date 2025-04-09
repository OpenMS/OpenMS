// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TARGETED/IncludeExcludeTarget.h>

#include <utility>

namespace OpenMS
{
  IncludeExcludeTarget::IncludeExcludeTarget() :
    CVTermList(),
    precursor_mz_(std::numeric_limits<double>::max()),
    product_mz_(std::numeric_limits<double>::max())
  {
  }

  IncludeExcludeTarget::IncludeExcludeTarget(const IncludeExcludeTarget & rhs) = default;

  IncludeExcludeTarget::~IncludeExcludeTarget() = default;

  IncludeExcludeTarget & IncludeExcludeTarget::operator=(const IncludeExcludeTarget & rhs)
  {
    if (&rhs != this)
    {
      CVTermList::operator=(rhs);
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
      rts_ = rhs.rts_;
    }
    return *this;
  }

  bool IncludeExcludeTarget::operator==(const IncludeExcludeTarget & rhs) const
  {
    return CVTermList::operator==(rhs) &&
           name_ == rhs.name_ &&
           precursor_mz_ == rhs.precursor_mz_ &&
           precursor_cv_terms_ == rhs.precursor_cv_terms_ &&
           product_mz_ == rhs.product_mz_ &&
           product_cv_terms_ == rhs.product_cv_terms_ &&
           interpretation_list_ == rhs.interpretation_list_ &&
           peptide_ref_ == rhs.peptide_ref_ &&
           compound_ref_ == rhs.compound_ref_ &&
           configurations_ == rhs.configurations_ &&
           prediction_ == rhs.prediction_ &&
           rts_ == rhs.rts_;
  }

  bool IncludeExcludeTarget::operator!=(const IncludeExcludeTarget & rhs) const
  {
    return !(*this == rhs);
  }

  void IncludeExcludeTarget::setName(const String & name)
  {
    name_ = name;
  }

  const String & IncludeExcludeTarget::getName() const
  {
    return name_;
  }

  void IncludeExcludeTarget::setPeptideRef(const String & peptide_ref)
  {
    peptide_ref_ = peptide_ref;
  }

  const String & IncludeExcludeTarget::getPeptideRef() const
  {
    return peptide_ref_;
  }

  void IncludeExcludeTarget::setCompoundRef(const String & compound_ref)
  {
    compound_ref_ = compound_ref;
  }

  const String & IncludeExcludeTarget::getCompoundRef() const
  {
    return compound_ref_;
  }

  void IncludeExcludeTarget::setPrecursorMZ(double mz)
  {
    precursor_mz_ = mz;
  }

  double IncludeExcludeTarget::getPrecursorMZ() const
  {
    return precursor_mz_;
  }

  void IncludeExcludeTarget::setPrecursorCVTermList(const CVTermList & list)
  {
    precursor_cv_terms_ = list;
  }

  void IncludeExcludeTarget::addPrecursorCVTerm(const CVTerm & cv_term)
  {
    precursor_cv_terms_.addCVTerm(cv_term);
  }

  const CVTermList & IncludeExcludeTarget::getPrecursorCVTermList() const
  {
    return precursor_cv_terms_;
  }

  void IncludeExcludeTarget::setProductMZ(double mz)
  {
    product_mz_ = mz;
  }

  double IncludeExcludeTarget::getProductMZ() const
  {
    return product_mz_;
  }

  void IncludeExcludeTarget::setProductCVTermList(const CVTermList & list)
  {
    product_cv_terms_ = list;
  }

  void IncludeExcludeTarget::addProductCVTerm(const CVTerm & cv_term)
  {
    product_cv_terms_.addCVTerm(cv_term);
  }

  const CVTermList & IncludeExcludeTarget::getProductCVTermList() const
  {
    return product_cv_terms_;
  }

  void IncludeExcludeTarget::setInterpretations(const std::vector<CVTermList> & interpretations)
  {
    interpretation_list_ = interpretations;
  }

  const std::vector<CVTermList> & IncludeExcludeTarget::getInterpretations() const
  {
    return interpretation_list_;
  }

  void IncludeExcludeTarget::addInterpretation(const CVTermList & interpretation)
  {
    interpretation_list_.push_back(interpretation);
  }

  void IncludeExcludeTarget::setConfigurations(const std::vector<Configuration> & configurations)
  {
    configurations_ = configurations;
  }

  const std::vector<IncludeExcludeTarget::Configuration> & IncludeExcludeTarget::getConfigurations() const
  {
    return configurations_;
  }

  void IncludeExcludeTarget::addConfiguration(const Configuration & configuration)
  {
    configurations_.push_back(configuration);
  }

  void IncludeExcludeTarget::setPrediction(const CVTermList & prediction)
  {
    prediction_ = prediction;
  }

  const CVTermList & IncludeExcludeTarget::getPrediction() const
  {
    return prediction_;
  }

  void IncludeExcludeTarget::addPredictionTerm(const CVTerm & term)
  {
    prediction_.addCVTerm(term);
  }

  void IncludeExcludeTarget::updateMembers_()
  {
  }

  void IncludeExcludeTarget::setRetentionTime(IncludeExcludeTarget::RetentionTime rt)
  {
    rts_ = std::move(rt);
  }

  const IncludeExcludeTarget::RetentionTime & IncludeExcludeTarget::getRetentionTime() const
  {
    return rts_;
  }

}
