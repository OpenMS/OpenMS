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

#ifndef OPENMS_ANALYSIS_MRM_REACTIONMONITORINGTRANSITION_H
#define OPENMS_ANALYSIS_MRM_REACTIONMONITORINGTRANSITION_H

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/CVTermList.h>

#include <vector>
#include <bitset>

namespace OpenMS
{
  /**
    @brief This class stores a SRM/MRM transition

    This class is capable of representing a <Transition> tag in a TraML
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
    ~ReactionMonitoringTransition() override;
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
    void setPrecursorMZ(double mz);

    /// get the precursor mz (Q1 value)
    double getPrecursorMZ() const;

    /// Returns true if precursor CV Terms exist (means it is safe to call getPrecursorCVTermList)
    bool hasPrecursorCVTerms() const;

    void setPrecursorCVTermList(const CVTermList & list);

    void addPrecursorCVTerm(const CVTerm & cv_term);

    /* @brief Obtain the list of CV Terms for the precursor
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

    void addIntermediateProduct(Product product);

    void setIntermediateProducts(const std::vector<Product> & products);

    void setProduct(Product product);

    const Product & getProduct() const;

    void setRetentionTime(RetentionTime rt);

    const RetentionTime & getRetentionTime() const;

    /// Returns true if a Prediction object exists (means it is safe to call getPrediction)
    bool hasPrediction() const;

    void setPrediction(const Prediction & prediction);

    void addPredictionTerm(const CVTerm & prediction);

    /* @brief Obtain the Prediction object 
     *
     * @note You first need to check whether the object is accessible using hasPrediction() 
    */
    const Prediction & getPrediction() const;

    DecoyTransitionType getDecoyTransitionType() const;

    void setDecoyTransitionType(const DecoyTransitionType & d);

    double getLibraryIntensity() const;

    void setLibraryIntensity(double intensity);

    /** @brief Detecting transitions
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
    bool isDetectingTransition() const;

    void setDetectingTransition(bool val);

    /** @brief Identifying transitions
     * 
     * Identifying transitions represent the set of transitions of an assay
     * that should be used for the independent identification of a candidate
     * peak group. These transitions will be scored independently of the
     * detecting transitions.
     *
    */
    bool isIdentifyingTransition() const;

    void setIdentifyingTransition(bool val);

    /** @brief Quantifying transitions
     * 
     * Quantifying transitions represent the set of transitions of an assay
     * that should be used for the quantification of the peptide. This includes
     * exclusion of e.g. interfered transitions (example: light/heavy peptide
     * pairs isolated in the same swath window), that should not be used for
     * quantification of the peptide.
     *
    */
    bool isQuantifyingTransition() const;

    void setQuantifyingTransition(bool val);

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

    /// A transition has exactly one precursor and it must supply the CV Term 1000827 (isolation window target m/z)
    double precursor_mz_;

    /// Intensity of the product (q3) ion (stored in CV Term 1001226 inside the <Transition> tag)
    double library_intensity_;

    /// specific properties of a transition (e.g. specific CV terms)
    DecoyTransitionType decoy_type_;

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
    // (currently: detecting, identifying and quantifying)
    std::bitset<3> transition_flags_;
    //@}
  };
}

#endif // OPENMS_ANALYSIS_MRM_REACTIONMONITORINGTRANSITION_H
