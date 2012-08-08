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

#ifndef OPENMS_ANALYSIS_MRM_REACTIONMONITORINGTRANSITION_H
#define OPENMS_ANALYSIS_MRM_REACTIONMONITORINGTRANSITION_H

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/CVTermList.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief This class stores a SRM/MRM transition

    The default values for precursor and product m/z values are
    set to numeric_limits<DoubleReal>::max(). Default values for
    precursor an product charge is set to numeric_limits<Int>::max().
  */
  class OPENMS_DLLAPI ReactionMonitoringTransition :
    public CVTermList
  {

public:

    typedef TargetedExperimentHelper::Configuration Configuration;
    typedef TargetedExperimentHelper::RetentionTime RetentionTime;
    typedef TargetedExperimentHelper::TraMLProduct Product;
    typedef TargetedExperimentHelper::Prediction Prediction;

    enum DecoyTransitionType
    {
      UNKNOWN,
      TARGET,
      DECOY
    };

    /** @name Constructors and destructors
    */
    //@{
    /// default constructor
    ReactionMonitoringTransition();

    /// copy constructor
    ReactionMonitoringTransition(const ReactionMonitoringTransition & rhs);

    /// destructor
    virtual ~ReactionMonitoringTransition();
    //@}

    /// assignment operator
    ReactionMonitoringTransition & operator=(const ReactionMonitoringTransition & rhs);

    /** @name Accessors
    */
    //@{
    void setName(const String & name);

    const String & getName() const;

    void setNativeID(const String & name);

    const String & getNativeID() const;

    void setPeptideRef(const String & peptide_ref);

    const String & getPeptideRef() const;

    void setCompoundRef(const String & compound_ref);

    const String & getCompoundRef() const;

    /// sets the precursor mz (Q1 value)
    void setPrecursorMZ(DoubleReal mz);

    DoubleReal getPrecursorMZ() const;

    void setPrecursorCVTermList(const CVTermList & list);

    void addPrecursorCVTerm(const CVTerm & cv_term);

    const CVTermList & getPrecursorCVTermList() const;

    void setProductMZ(DoubleReal mz);

    DoubleReal getProductMZ() const;

    //void setProductCVTermList(const CVTermList& list);

    void addProductCVTerm(const CVTerm & cv_term);

    //const CVTermList& getProductCVTermList() const;

    const std::vector<Product> & getIntermediateProducts() const;

    void addIntermediateProduct(Product product);

    void setIntermediateProducts(const std::vector<Product> & products);

    void setProduct(Product product);

    const Product & getProduct() const;

    void setRetentionTime(RetentionTime rt);

    const RetentionTime & getRetentionTime() const;

    void setPrediction(const Prediction & prediction);

    void addPredictionTerm(const CVTerm & prediction);

    const Prediction & getPrediction() const;

    DecoyTransitionType getDecoyTransitionType() const;

    void setDecoyTransitionType(const DecoyTransitionType & d);

    DoubleReal getLibraryIntensity() const;

    void setLibraryIntensity(DoubleReal intensity);

    //@}

    /** @name Predicates
    */
    //@{
    /// equality operator
    bool operator==(const ReactionMonitoringTransition & rhs) const;

    /// inequality operator
    bool operator!=(const ReactionMonitoringTransition & rhs) const;
    //@}

    /**  @name  Comparator classes.
       These classes implement binary predicates that can be used
       to compare two transitions with respect to their product ions.
   */
    //@{
    /// Comparator by Product ion MZ
    struct ProductMZLess :
      std::binary_function<ReactionMonitoringTransition, ReactionMonitoringTransition, bool>
    {
      inline bool operator()(ReactionMonitoringTransition const & left, ReactionMonitoringTransition const & right) const
      {
        return left.getProductMZ() < right.getProductMZ();
      }

    };
    //@}

protected:

    void updateMembers_();

    //@{ 
    /// Attributes:

    String name_; // id, required attribute

    // attributes to a peptide / compound (optional)
    String peptide_ref_;
    String compound_ref_;
    //@}

    //@{ 
    /// Subelements:

    // Precursor
    // Product
    // IntermediateProduct
    // RetentionTime
    // Prediction
    // cvparam / userParam

    // A transition has exactly one precursor and it must supply the CV Term 1000827 (isolation window target m/z
    DoubleReal precursor_mz_;
    CVTermList precursor_cv_terms_;

    Product product_;

    std::vector<Product> intermediate_products_;

    RetentionTime rts;

    Prediction prediction_;
    //@}

    /// specific properties of a transition (e.g. specific CV terms)
    DecoyTransitionType decoy_type_;

    DoubleReal library_intensity_;

  };
}

#endif // OPENMS_ANALYSIS_MRM_REACTIONMONITORINGTRANSITION_H
