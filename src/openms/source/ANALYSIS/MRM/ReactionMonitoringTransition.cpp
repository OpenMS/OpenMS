// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>

#include <algorithm>

namespace OpenMS
{

  ReactionMonitoringTransition::ReactionMonitoringTransition() :
    CVTermList(),
    precursor_mz_(0.0),
    decoy_type_(UNKNOWN),
    library_intensity_(-101)
  {
  }

  ReactionMonitoringTransition::ReactionMonitoringTransition(const ReactionMonitoringTransition & rhs) :
    CVTermList(rhs),
    name_(rhs.name_),
    peptide_ref_(rhs.peptide_ref_),
    compound_ref_(rhs.compound_ref_),
    precursor_mz_(rhs.precursor_mz_),
    precursor_cv_terms_(rhs.precursor_cv_terms_),
    product_(rhs.product_),
    intermediate_products_(rhs.intermediate_products_),
    rts(rhs.rts),
    prediction_(rhs.prediction_),
    decoy_type_(rhs.decoy_type_),
    library_intensity_(rhs.library_intensity_)
  {
  }

  ReactionMonitoringTransition::~ReactionMonitoringTransition()
  {
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
      precursor_cv_terms_ = rhs.precursor_cv_terms_;
      intermediate_products_ = rhs.intermediate_products_;
      product_ = rhs.product_;
      rts = rhs.rts;
      prediction_ = rhs.prediction_;
      decoy_type_ = rhs.decoy_type_;
      library_intensity_ = rhs.library_intensity_;
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
           precursor_cv_terms_ == rhs.precursor_cv_terms_ &&
           product_ == rhs.product_ &&
           intermediate_products_ == rhs.intermediate_products_ &&
           rts == rhs.rts &&
           prediction_ == rhs.prediction_ &&
           decoy_type_ == rhs.decoy_type_ &&
           library_intensity_ == rhs.library_intensity_;
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

  void ReactionMonitoringTransition::setPrecursorCVTermList(const CVTermList & list)
  {
    precursor_cv_terms_ = list;
  }

  void ReactionMonitoringTransition::addPrecursorCVTerm(const CVTerm & cv_term)
  {
    precursor_cv_terms_.addCVTerm(cv_term);
  }

  const CVTermList & ReactionMonitoringTransition::getPrecursorCVTermList() const
  {
    return precursor_cv_terms_;
  }

  void ReactionMonitoringTransition::setProductMZ(double mz)
  {
    CVTerm product_mz;
    std::vector<CVTerm> product_cvterms;
    product_mz.setCVIdentifierRef("MS");
    product_mz.setAccession("MS:1000827");
    product_mz.setName("isolation window target m/z");
    product_mz.setValue(mz);
    product_cvterms.push_back(product_mz);

    Map<String, std::vector<CVTerm> >  cvtermlist = product_.getCVTerms();
    cvtermlist[product_mz.getAccession()] = product_cvterms;

    product_.replaceCVTerms(cvtermlist);
  }

  double ReactionMonitoringTransition::getProductMZ() const
  {
    try
    {
      if (!product_.getCVTerms().has("MS:1000827"))
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Product mz has never been set");
      }
      return product_.getCVTerms()["MS:1000827"][0].getValue().toString().toDouble();
    }
    catch (char * /*str*/)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Product mz has never been set");
    }
  }

  int ReactionMonitoringTransition::getProductChargeState() const
  { 
    try
    { 
      if (!product_.getChargeState() > 0)
      { 
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Product charge has never been set");
      }
      return product_.getChargeState();
    }
    catch (char * /*str*/)
    { 
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Product charge has never been set");
    }
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

  void ReactionMonitoringTransition::setPrediction(const Prediction & prediction)
  {
    prediction_ = prediction;
  }

  const ReactionMonitoringTransition::Prediction & ReactionMonitoringTransition::getPrediction() const
  {
    return prediction_;
  }

  void ReactionMonitoringTransition::addPredictionTerm(const CVTerm & term)
  {
    prediction_.addCVTerm(term);
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
    bool detecting = true;
    if (this->metaValueExists("detecting_transition"))
    {
      if (!this->getMetaValue("detecting_transition").toBool())
      {
        detecting = false;
      }
      else if (this->getMetaValue("detecting_transition").toBool())
      {
        detecting = true;
      }
    }

    return detecting;
  }

  bool ReactionMonitoringTransition::isIdentifyingTransition() const
  {
    bool identifying = false;
    if (this->metaValueExists("identifying_transition"))
    {
      if (!this->getMetaValue("identifying_transition").toBool())
      {
        identifying = false;
      }
      else if (this->getMetaValue("identifying_transition").toBool())
      {
        identifying = true;
      }
    }

    return identifying;
  }

  bool ReactionMonitoringTransition::isQuantifyingTransition() const
  {
    bool quantifying = true;
    if (this->metaValueExists("quantifying_transition"))
    {
      if (!this->getMetaValue("quantifying_transition").toBool())
      {
        quantifying = false;
      }
      else if (this->getMetaValue("quantifying_transition").toBool())
      {
        quantifying = true;
      }
    }

    return quantifying;
  }

} // namespace OpenMS
