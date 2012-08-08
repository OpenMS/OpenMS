// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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

namespace OpenMS
{

  ReactionMonitoringTransition::ReactionMonitoringTransition() :
    CVTermList(),
    precursor_mz_(std::numeric_limits<DoubleReal>::max()),
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
    prediction_(rhs.prediction_),
    rts(rhs.rts),
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
      prediction_ = rhs.prediction_;
      product_ = rhs.product_;
      rts = rhs.rts;
      decoy_type_ = rhs.decoy_type_;
      library_intensity_ = rhs.library_intensity_;
    }
    return *this;
  }

  bool ReactionMonitoringTransition::operator==(const ReactionMonitoringTransition & rhs) const
  {
    return CVTermList::operator == (rhs) &&
           name_ == rhs.name_ &&
           peptide_ref_ == rhs.peptide_ref_ &&
           compound_ref_ == rhs.compound_ref_ &&
           precursor_mz_ == rhs.precursor_mz_ &&
           precursor_cv_terms_ == rhs.precursor_cv_terms_ &&
           product_ == rhs.product_ &&
           intermediate_products_ == rhs.intermediate_products_ &&
           prediction_ == rhs.prediction_ &&
           rts == rhs.rts &&
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

  void ReactionMonitoringTransition::setPrecursorMZ(DoubleReal mz)
  {
    precursor_mz_ = mz;
  }

  DoubleReal ReactionMonitoringTransition::getPrecursorMZ() const
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

  void ReactionMonitoringTransition::setProductMZ(DoubleReal mz)
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

  DoubleReal ReactionMonitoringTransition::getProductMZ() const
  {
    try
    {
      return product_.getCVTerms()["MS:1000827"][0].getValue().toString().toDouble();
    }
    catch (char * /*str*/)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Product mz has never been set");
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

  DoubleReal ReactionMonitoringTransition::getLibraryIntensity() const
  {
    return library_intensity_;
  }

  void ReactionMonitoringTransition::setLibraryIntensity(const DoubleReal intensity)
  {
    library_intensity_ = intensity;
  }

} // namespace OpenMS
