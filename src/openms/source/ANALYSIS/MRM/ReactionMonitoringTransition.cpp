// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>

#include <OpenMS/CONCEPT/Helpers.h>
#include <algorithm>

namespace OpenMS
{

  static const unsigned char DETECTING_TRANSITION_LOC = 0;
  static const unsigned char IDENTIFYING_TRANSITION_LOC = 1;
  static const unsigned char QUANTIFYING_TRANSITION_LOC = 2;

  ReactionMonitoringTransition::ReactionMonitoringTransition() :
    CVTermList(),
    precursor_mz_(0.0),
    library_intensity_(-101),
    decoy_type_(UNKNOWN),
    precursor_cv_terms_(nullptr),
    prediction_(nullptr)
  {
    // Default is: true, false, true
    // NOTE: do not change that, the same default is implicitely assumed in TraMLHandler
    transition_flags_[DETECTING_TRANSITION_LOC] = true;
    transition_flags_[IDENTIFYING_TRANSITION_LOC] = false;
    transition_flags_[QUANTIFYING_TRANSITION_LOC] = true;
  }

  ReactionMonitoringTransition::ReactionMonitoringTransition(const ReactionMonitoringTransition & rhs) :
    CVTermList(rhs),
    name_(rhs.name_),
    peptide_ref_(rhs.peptide_ref_),
    compound_ref_(rhs.compound_ref_),
    precursor_mz_(rhs.precursor_mz_),
    library_intensity_(rhs.library_intensity_),
    decoy_type_(rhs.decoy_type_),
    precursor_cv_terms_(nullptr),
    product_(rhs.product_),
    intermediate_products_(rhs.intermediate_products_),
    rts(rhs.rts),
    prediction_(nullptr),
    transition_flags_(rhs.transition_flags_)
  {
    // We copy the internal object (not just the ptr)
    if (rhs.precursor_cv_terms_ != nullptr)
    {
      precursor_cv_terms_ = new CVTermList(*rhs.precursor_cv_terms_);
    }
    if (rhs.prediction_ != nullptr)
    {
      prediction_ = new Prediction(*rhs.prediction_);
    }
  }

  ReactionMonitoringTransition::~ReactionMonitoringTransition()
  {
    delete precursor_cv_terms_;
    delete prediction_;
  }

  ReactionMonitoringTransition & ReactionMonitoringTransition::operator=(const ReactionMonitoringTransition & rhs)
  {
    if (&rhs != this)
    {
      CVTermList::operator=(rhs);
      name_ = rhs.name_;
      peptide_ref_ = rhs.peptide_ref_;
      compound_ref_ = rhs.compound_ref_;
      precursor_mz_ = rhs.precursor_mz_;
      intermediate_products_ = rhs.intermediate_products_;
      product_ = rhs.product_;
      rts = rhs.rts;
      library_intensity_ = rhs.library_intensity_;
      decoy_type_ = rhs.decoy_type_;
      transition_flags_ = rhs.transition_flags_;

      // We copy the internal object (not just the ptr)
      delete precursor_cv_terms_;
      precursor_cv_terms_ = nullptr;
      if (rhs.precursor_cv_terms_ != nullptr)
      {
        precursor_cv_terms_ = new CVTermList(*rhs.precursor_cv_terms_);
      }

      // We copy the internal object (not just the ptr)
      delete prediction_;
      prediction_ = nullptr;
      if (rhs.prediction_ != nullptr)
      {
        prediction_ = new Prediction(*rhs.prediction_);
      }
    }
    return *this;
  }

  bool ReactionMonitoringTransition::operator==(const ReactionMonitoringTransition & rhs) const
  {
    return CVTermList::operator==(rhs) &&
           name_ == rhs.name_ &&
           peptide_ref_ == rhs.peptide_ref_ &&
           compound_ref_ == rhs.compound_ref_ &&
           precursor_mz_ == rhs.precursor_mz_ &&
           OpenMS::Helpers::cmpPtrSafe< CVTermList* >(precursor_cv_terms_, rhs.precursor_cv_terms_) &&
           product_ == rhs.product_ &&
           intermediate_products_ == rhs.intermediate_products_ &&
           rts == rhs.rts &&
           OpenMS::Helpers::cmpPtrSafe< Prediction* >(prediction_, rhs.prediction_) &&
           library_intensity_ == rhs.library_intensity_ &&
           decoy_type_ == rhs.decoy_type_ &&
           transition_flags_ == rhs.transition_flags_;
  }

  bool ReactionMonitoringTransition::operator!=(const ReactionMonitoringTransition & rhs) const
  {
    return !(*this == rhs);
  }

  void ReactionMonitoringTransition::setName(const String & name)
  {
    name_ = name;
  }

  const String & ReactionMonitoringTransition::getName() const
  {
    return name_;
  }

  void ReactionMonitoringTransition::setNativeID(const String & name)
  {
    name_ = name;
  }

  const String & ReactionMonitoringTransition::getNativeID() const
  {
    return name_;
  }

  void ReactionMonitoringTransition::setPeptideRef(const String & peptide_ref)
  {
    peptide_ref_ = peptide_ref;
  }

  const String & ReactionMonitoringTransition::getPeptideRef() const
  {
    return peptide_ref_;
  }

  void ReactionMonitoringTransition::setCompoundRef(const String & compound_ref)
  {
    compound_ref_ = compound_ref;
  }

  const String & ReactionMonitoringTransition::getCompoundRef() const
  {
    return compound_ref_;
  }

  void ReactionMonitoringTransition::setPrecursorMZ(double mz)
  {
    precursor_mz_ = mz;
  }

  double ReactionMonitoringTransition::getPrecursorMZ() const
  {
    return precursor_mz_;
  }

  bool ReactionMonitoringTransition::hasPrecursorCVTerms() const
  {
    return (precursor_cv_terms_ != nullptr);
  }

  void ReactionMonitoringTransition::setPrecursorCVTermList(const CVTermList & list)
  {
    delete precursor_cv_terms_;
    precursor_cv_terms_ = new CVTermList(list);
  }

  void ReactionMonitoringTransition::addPrecursorCVTerm(const CVTerm & cv_term)
  {
    if (!precursor_cv_terms_)
    {
      precursor_cv_terms_ = new CVTermList();
    }
    precursor_cv_terms_->addCVTerm(cv_term);
  }

  const CVTermList & ReactionMonitoringTransition::getPrecursorCVTermList() const
  {
    OPENMS_PRECONDITION(hasPrecursorCVTerms(), "ReactionMonitoringTransition has no PrecursorCVTerms, check first with hasPrecursorCVTerms()")
    return *precursor_cv_terms_;
  }

  void ReactionMonitoringTransition::setProductMZ(double mz)
  {
    product_.setMZ(mz);
  }

  double ReactionMonitoringTransition::getProductMZ() const
  {
    return product_.getMZ();
  }

  int ReactionMonitoringTransition::getProductChargeState() const
  { 
    return product_.getChargeState();
  }

  bool ReactionMonitoringTransition::isProductChargeStateSet() const
  { 
    return product_.hasCharge();
  }

  void ReactionMonitoringTransition::addProductCVTerm(const CVTerm & cv_term)
  {
    product_.addCVTerm(cv_term);
  }

  const std::vector<ReactionMonitoringTransition::Product> & ReactionMonitoringTransition::getIntermediateProducts() const
  {
    return intermediate_products_;
  }

  void ReactionMonitoringTransition::addIntermediateProduct(ReactionMonitoringTransition::Product product)
  {
    intermediate_products_.push_back(product);
  }

  void ReactionMonitoringTransition::setIntermediateProducts(const std::vector<ReactionMonitoringTransition::Product> & intermediate_products)
  {
    intermediate_products_ = intermediate_products;
  }

  void ReactionMonitoringTransition::setProduct(ReactionMonitoringTransition::Product product)
  {
    product_ = product;
  }

  const ReactionMonitoringTransition::Product & ReactionMonitoringTransition::getProduct() const
  {
    return product_;
  }

  void ReactionMonitoringTransition::setRetentionTime(ReactionMonitoringTransition::RetentionTime rt)
  {
    rts = rt;
  }

  const ReactionMonitoringTransition::RetentionTime & ReactionMonitoringTransition::getRetentionTime() const
  {
    return rts;
  }

  bool ReactionMonitoringTransition::hasPrediction() const
  {
    return (prediction_ != nullptr);
  }

  void ReactionMonitoringTransition::setPrediction(const Prediction & prediction)
  {
    delete prediction_;
    prediction_ = new Prediction(prediction);
  }

  const ReactionMonitoringTransition::Prediction & ReactionMonitoringTransition::getPrediction() const
  {
    OPENMS_PRECONDITION(hasPrecursorCVTerms(), "ReactionMonitoringTransition has no Prediction object, check first with hasPrediction()")
    return *prediction_;
  }

  void ReactionMonitoringTransition::addPredictionTerm(const CVTerm & term)
  {
    if (!prediction_)
    {
      prediction_ = new Prediction();
    }
    prediction_->addCVTerm(term);
  }

  void ReactionMonitoringTransition::updateMembers_()
  {
  }

  ReactionMonitoringTransition::DecoyTransitionType ReactionMonitoringTransition::getDecoyTransitionType() const
  {
    return decoy_type_;
  }

  void ReactionMonitoringTransition::setDecoyTransitionType(const DecoyTransitionType & d)
  {
    decoy_type_ = d;
  }

  double ReactionMonitoringTransition::getLibraryIntensity() const
  {
    return library_intensity_;
  }

  void ReactionMonitoringTransition::setLibraryIntensity(const double intensity)
  {
    library_intensity_ = intensity;
  }

  bool ReactionMonitoringTransition::isDetectingTransition() const
  {
    return transition_flags_[DETECTING_TRANSITION_LOC];
  }

  void ReactionMonitoringTransition::setDetectingTransition(bool val)
  {
    transition_flags_[DETECTING_TRANSITION_LOC] = val;
  }

  bool ReactionMonitoringTransition::isIdentifyingTransition() const
  {
    return transition_flags_[IDENTIFYING_TRANSITION_LOC];
  }

  void ReactionMonitoringTransition::setIdentifyingTransition(bool val)
  {
    transition_flags_[IDENTIFYING_TRANSITION_LOC] = val;
  }

  bool ReactionMonitoringTransition::isQuantifyingTransition() const
  {
    return transition_flags_[QUANTIFYING_TRANSITION_LOC];
  }

  void ReactionMonitoringTransition::setQuantifyingTransition(bool val)
  {
    transition_flags_[QUANTIFYING_TRANSITION_LOC] = val;
  }

} // namespace OpenMS
