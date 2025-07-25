// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>

#include <OpenMS/CONCEPT/Helpers.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <utility>

namespace OpenMS
{

  static const unsigned char DETECTING_TRANSITION_LOC = 0;
  static const unsigned char IDENTIFYING_TRANSITION_LOC = 1;
  static const unsigned char QUANTIFYING_TRANSITION_LOC = 2;

  ReactionMonitoringTransition::ReactionMonitoringTransition() :
    CVTermList(),
    trans_type_(OpenSwath::TransType::UNKNOWN),
    library_intensity_(-101),
    decoy_type_(UNKNOWN),
    precursor_mz_(0.0),
    precursor_cv_terms_(nullptr),
    prediction_(nullptr)
  {
    // Default is: true, false, true
    // NOTE: do not change that, the same default is implicitly assumed in TraMLHandler
    transition_flags_[DETECTING_TRANSITION_LOC] = true;
    transition_flags_[IDENTIFYING_TRANSITION_LOC] = false;
    transition_flags_[QUANTIFYING_TRANSITION_LOC] = true;
  }

  ReactionMonitoringTransition::ReactionMonitoringTransition(const ReactionMonitoringTransition & rhs) :
    CVTermList(rhs),
    name_(rhs.name_),
    transition_ref_(rhs.transition_ref_),
    trans_type_(rhs.trans_type_),
    library_intensity_(rhs.library_intensity_),
    decoy_type_(rhs.decoy_type_),
    precursor_mz_(rhs.precursor_mz_),
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

  ReactionMonitoringTransition::ReactionMonitoringTransition(ReactionMonitoringTransition && rhs) noexcept :
    CVTermList(std::move(rhs)),
    name_(std::move(rhs.name_)),
    transition_ref_(std::move(rhs.transition_ref_)),
    trans_type_(std::move(rhs.trans_type_)),
    library_intensity_(std::move(rhs.library_intensity_)),
    decoy_type_(std::move(rhs.decoy_type_)),
    precursor_mz_(std::move(rhs.precursor_mz_)),
    precursor_cv_terms_(std::move(rhs.precursor_cv_terms_)),
    product_(std::move(rhs.product_)),
    intermediate_products_(std::move(rhs.intermediate_products_)),
    rts(std::move(rhs.rts)),
    prediction_(std::move(rhs.prediction_)),
    transition_flags_(std::move(rhs.transition_flags_))
  {
    rhs.precursor_cv_terms_ = nullptr;
    rhs.prediction_ = nullptr;
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
      transition_ref_ = rhs.transition_ref_;
      trans_type_ = rhs.trans_type_;
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

  ReactionMonitoringTransition & ReactionMonitoringTransition::operator=(ReactionMonitoringTransition && rhs) noexcept
  {
    if (&rhs != this)
    {
      CVTermList::operator=(std::move(rhs));
      name_ = std::move(rhs.name_);
      transition_ref_ = std::move(rhs.transition_ref_);
      trans_type_ = std::move(rhs.trans_type_);
      precursor_mz_ = std::move(rhs.precursor_mz_);
      intermediate_products_ = std::move(rhs.intermediate_products_);
      product_ = std::move(rhs.product_);
      rts = std::move(rhs.rts);
      library_intensity_ = std::move(rhs.library_intensity_);
      decoy_type_ = std::move(rhs.decoy_type_);
      transition_flags_ = std::move(rhs.transition_flags_);

      // Move the ptr-based objects to the current objects and delete them in the rhs
      delete precursor_cv_terms_;
      precursor_cv_terms_ = rhs.precursor_cv_terms_;
      rhs.precursor_cv_terms_ = nullptr;

      delete prediction_;
      prediction_ = rhs.prediction_;
      rhs.prediction_ = nullptr;
    }
    return *this;
  }

  bool ReactionMonitoringTransition::operator==(const ReactionMonitoringTransition & rhs) const
  {
    return CVTermList::operator==(rhs) &&
           name_ == rhs.name_ &&
          transition_ref_ == rhs.transition_ref_ &&
           trans_type_ == rhs.trans_type_ &&
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

  void ReactionMonitoringTransition::setTransGroupRef(const String & transition_ref, OpenSwath::TransType type)
  {
    trans_type_ = type;
    transition_ref_ = transition_ref;
  }

  const String & ReactionMonitoringTransition::getTransGroupRef() const
  {
    return transition_ref_;
  }

  void ReactionMonitoringTransition::setTransGroupType(OpenSwath::TransType type)
  {
    trans_type_ = type;
  }

  OpenSwath::TransType ReactionMonitoringTransition::getTransGroupType() const
  {
    return trans_type_;
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

  void ReactionMonitoringTransition::addIntermediateProduct(const ReactionMonitoringTransition::Product& product)
  {
    intermediate_products_.push_back(product);
  }

  void ReactionMonitoringTransition::setIntermediateProducts(const std::vector<ReactionMonitoringTransition::Product> & intermediate_products)
  {
    intermediate_products_ = intermediate_products;
  }

  void ReactionMonitoringTransition::setProduct(ReactionMonitoringTransition::Product product)
  {
    product_ = std::move(product);
  }

  const ReactionMonitoringTransition::Product & ReactionMonitoringTransition::getProduct() const
  {
    return product_;
  }

  void ReactionMonitoringTransition::setRetentionTime(ReactionMonitoringTransition::RetentionTime rt)
  {
    rts = std::move(rt);
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
