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
    void setPrecursorMZ(double mz);

    double getPrecursorMZ() const;

    void setPrecursorCVTermList(const CVTermList & list);

    void addPrecursorCVTerm(const CVTerm & cv_term);

    const CVTermList & getPrecursorCVTermList() const;

    void setProductMZ(double mz);

    double getProductMZ() const;

    int getProductChargeState() const;

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

    double getLibraryIntensity() const;

    void setLibraryIntensity(double intensity);

    /** @brief Detecting transitions
     * 
     * Detecting transitions represent the set of transitions of an assay that should be used
     * to enable candidate peak group detection. Ideally they were observed and validated in 
     * previous experiments and have an associated library intensity.
     * For an overview, see Schubert OT et al., Building high-quality assay libraries for
     * targeted analysis of SWATH MS data., Nat Protoc. 2015 Mar;10(3):426-41.
     * doi: 10.1038/nprot.2015.015. Epub 2015 Feb 12. PMID: 25675208 
     *
    */
    bool isDetectingTransition() const;

    /** @brief Identifying transitions
     * 
     * Identifying transitions represent the set of transitions of an assay that should be used
     * for the independent identification of a candidate peak group. These transitions will
     * be scored independently of the detecting transitions.
     *
    */
    bool isIdentifyingTransition() const;

    /** @brief Quantifying transitions
     * 
     * Quantifying transitions represent the set of transitions of an assay that should be used
     * for the quantification of the peptide. This includes exclusion of e.g. interfered
     * transitions (example: light/heavy peptide pairs isolated in the same swath), that
     * should not be used for quantification of the peptide.
     *
    */
    bool isQuantifyingTransition() const;

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

    // A transition has exactly one precursor and it must supply the CV Term 1000827 (isolation window target m/z)
    double precursor_mz_;

    CVTermList precursor_cv_terms_;

    Product product_;

    std::vector<Product> intermediate_products_;

    RetentionTime rts;

    Prediction prediction_;
    //@}

    /// specific properties of a transition (e.g. specific CV terms)
    DecoyTransitionType decoy_type_;

    double library_intensity_;
  };
}

#endif // OPENMS_ANALYSIS_MRM_REACTIONMONITORINGTRANSITION_H
