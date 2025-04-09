// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/CVTermList.h>

#include <vector>
#include <bitset>

namespace OpenMS
{

  /**
    @brief This class stores a SRM/MRM transition

    This class is capable of representing a \<Transition\> tag in a TraML
    document completely and contains all associated information.

    The default values for precursor m/z is 0.0 which indicates that it is
    uninitialized.
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
      UNKNOWN, ///< Unknown type
      TARGET, ///< Target transition
      DECOY ///< Decoy transition 
    };

    /** @name Constructors and destructors
    */
    //@{
    /// default constructor
    ReactionMonitoringTransition();

    /// copy constructor
    ReactionMonitoringTransition(const ReactionMonitoringTransition & rhs);

    /// Move constructor
    ReactionMonitoringTransition(ReactionMonitoringTransition&&) noexcept;

    /// destructor
    ~ReactionMonitoringTransition() override;
    //@}

    /// assignment operator
    ReactionMonitoringTransition & operator=(const ReactionMonitoringTransition & rhs);

    /// move assignment operator
    ReactionMonitoringTransition & operator=(ReactionMonitoringTransition && rhs) noexcept ;

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
    void setPrecursorMZ(double mz);

    /// get the precursor mz (Q1 value)
    double getPrecursorMZ() const;

    /// Returns true if precursor CV Terms exist (means it is safe to call getPrecursorCVTermList)
    bool hasPrecursorCVTerms() const;

    void setPrecursorCVTermList(const CVTermList & list);

    void addPrecursorCVTerm(const CVTerm & cv_term);

    /** @brief Obtain the list of CV Terms for the precursor
     *
     * @note You first need to check whether they exist using hasPrecursorCVTerms() 
    */
    const CVTermList & getPrecursorCVTermList() const;

    void setProductMZ(double mz);

    double getProductMZ() const;

    int getProductChargeState() const;

    bool isProductChargeStateSet() const;

    void addProductCVTerm(const CVTerm & cv_term);

    const std::vector<Product> & getIntermediateProducts() const;

    void addIntermediateProduct(const Product& product);

    void setIntermediateProducts(const std::vector<Product> & products);

    void setProduct(Product product);

    const Product & getProduct() const;

    /// Returns true if a Prediction object exists (means it is safe to call getPrediction)
    bool hasPrediction() const;

    void setPrediction(const Prediction & prediction);

    void addPredictionTerm(const CVTerm & prediction);

    /** @brief Obtain the Prediction object 
     *
     * @note You first need to check whether the object is accessible using hasPrediction() 
    */
    const Prediction & getPrediction() const;

    /// Returns the type of transition (target or decoy)
    DecoyTransitionType getDecoyTransitionType() const;

    /// Sets the type of transition (target or decoy)
    void setDecoyTransitionType(const DecoyTransitionType & d);

    /// Returns the library intensity (ion count or normalized ion count from a spectral library)
    double getLibraryIntensity() const;

    /// Sets the library intensity (ion count or normalized ion count from a spectral library)
    void setLibraryIntensity(double intensity);

    /** @name Retention time
     * 
     * @brief Retention time set on transition level
     *
     * @note Retention time can be set on analyte level (peptide / compound
     * level) or transition level -- usually the retention time set here on
     * transition level is ignored, therefore only use this if you have a good
     * reason.
     *
    */
    //@{
    /** @brief Sets the retention time on transition level
     *
     * Stores a retention time on the level of an individual transition
     *
     * @note This is very likely not what you want, to set the retention time
     * of an analyte use TargetedExperimentHelper::PeptideCompound::getRetentionTime()
    */
    void setRetentionTime(RetentionTime rt);

    /** @brief Gets the retention time on transition level
     *
     * Stores a retention time on the level of an individual transition
     *
     * @note This is very likely not what you want, to set the retention time
     * of an analyte use TargetedExperimentHelper::PeptideCompound::getRetentionTime()
    */
    const RetentionTime & getRetentionTime() const;

    //@}

    /** @name Detecting transitions
     * 
     * @brief Transitions used for detection of analyte
     *
     * Detecting transitions represent the set of transitions of an assay that
     * should be used to enable candidate peak group detection. Ideally they
     * were observed and validated in previous experiments and have an
     * associated library intensity.
     *
     * For an overview, see Schubert OT et al., Building high-quality assay libraries for
     * targeted analysis of SWATH MS data., Nat Protoc. 2015 Mar;10(3):426-41.
     * doi: 10.1038/nprot.2015.015. Epub 2015 Feb 12. PMID: 25675208 
     *
    */
    //@{
    bool isDetectingTransition() const;

    void setDetectingTransition(bool val);
    //@}

    /** @name Identifying transitions
     * 
     * @brief Transitions used for differentiation between different analyte isomers (peptidoforms).
     *
     * Identifying transitions represent the set of transitions of an assay
     * that should be used for the independent identification of a candidate
     * peak group (often in relation to another analyte isomer or peptidoform). 
     * These transitions could for example be site-determining transitions in a 
     * phospho-proteomic DIA experiment. These transitions will be scored independently 
     * of the detecting transitions.
     *
    */
    //@{
    bool isIdentifyingTransition() const;

    void setIdentifyingTransition(bool val);
    //@}

    /** @name Quantifying transitions
     * 
     * @brief Transitions used for quantification
     *
     * Quantifying transitions represent the set of transitions of an assay
     * that should be used for the quantification of the peptide. This includes
     * exclusion of e.g. interfered transitions (example: light/heavy peptide
     * pairs isolated in the same swath window), that should not be used for
     * quantification of the peptide.
     *
    */
    //@{
    bool isQuantifyingTransition() const;

    void setQuantifyingTransition(bool val);
    //@}

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
    struct ProductMZLess
    {
      inline bool operator()(ReactionMonitoringTransition const & left, ReactionMonitoringTransition const & right) const
      {
        return left.getProductMZ() < right.getProductMZ();
      }
    };
    //@}

    /**  @name  Comparator classes.
    These classes implement binary predicates that can be used
    to compare two transitions with respect to their name.
    */
    //@{
    /// Comparator by name
    struct NameLess
    {
      inline bool operator()(ReactionMonitoringTransition const & left, ReactionMonitoringTransition const & right) const
      {
        return left.getName() < right.getName();
      }
    };
    //@}

protected:

    void updateMembers_();

    /** @name Attributes
    */
    //@{

    String name_; ///< id, required attribute

    // attributes to a peptide / compound (optional)
    String peptide_ref_; ///< Reference to a specific peptide
    String compound_ref_; ///< Reference to a specific compound

    /// Intensity of the product (q3) ion (stored in CV Term 1001226 inside the \<Transition\> tag)
    double library_intensity_;

    /// specific properties of a transition (e.g. specific CV terms)
    DecoyTransitionType decoy_type_;
    //@}

    /** @name Subelements
    */
    //@{

    // Precursor
    // Product
    // IntermediateProduct
    // RetentionTime
    // Prediction
    // cvparam / userParam

    /// A transition has exactly one precursor and it must supply the CV Term 1000827 (isolation window target m/z)
    double precursor_mz_;

    /// (Other) CV Terms of the Precursor (Q1) of the transition or target
    CVTermList* precursor_cv_terms_;

    /// Product (Q3) of the transition
    Product product_;

    /// Intermediate product ion information of the transition when using MS3 or above (optional)
    std::vector<Product> intermediate_products_;

    /// Information about predicted or calibrated retention time (optional)
    RetentionTime rts;

    /// Information about a prediction for a suitable transition using some software (optional)
    Prediction* prediction_;

    /// A set of flags to store information about the transition at hand
    /// (currently: detecting, identifying and quantifying)
    std::bitset<3> transition_flags_;
    //@}
  };
}

