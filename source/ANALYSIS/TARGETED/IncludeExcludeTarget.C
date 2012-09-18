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
// $Maintainer: Hannes Roest $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TARGETED/IncludeExcludeTarget.h>

#include <algorithm>

namespace OpenMS
{
  IncludeExcludeTarget::IncludeExcludeTarget() :
    CVTermList(),
    precursor_mz_(std::numeric_limits<DoubleReal>::max()),
    product_mz_(std::numeric_limits<DoubleReal>::max())
  {
  }

  IncludeExcludeTarget::IncludeExcludeTarget(const IncludeExcludeTarget & rhs) :
    CVTermList(rhs),
    name_(rhs.name_),
    precursor_mz_(rhs.precursor_mz_),
    precursor_cv_terms_(rhs.precursor_cv_terms_),
    product_mz_(rhs.product_mz_),
    product_cv_terms_(rhs.product_cv_terms_),
    interpretation_list_(rhs.interpretation_list_),
    peptide_ref_(rhs.peptide_ref_),
    compound_ref_(rhs.compound_ref_),
    configurations_(rhs.configurations_),
    prediction_(rhs.prediction_),
    rts(rhs.rts)
  {
  }

  IncludeExcludeTarget::~IncludeExcludeTarget()
  {
  }

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
      rts = rhs.rts;
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
           rts == rhs.rts;
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

  void IncludeExcludeTarget::setPrecursorMZ(DoubleReal mz)
  {
    precursor_mz_ = mz;
  }

  DoubleReal IncludeExcludeTarget::getPrecursorMZ() const
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

  void IncludeExcludeTarget::setProductMZ(DoubleReal mz)
  {
    product_mz_ = mz;
  }

  DoubleReal IncludeExcludeTarget::getProductMZ() const
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
    rts = rt;
  }

  const IncludeExcludeTarget::RetentionTime & IncludeExcludeTarget::getRetentionTime() const
  {
    return rts;
  }

}
