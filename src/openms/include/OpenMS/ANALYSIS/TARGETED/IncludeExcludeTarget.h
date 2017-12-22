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
// $Maintainer: Hannes Roest $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_TARGETED_INCLUDEEXCLUDETARGET_H
#define OPENMS_ANALYSIS_TARGETED_INCLUDEEXCLUDETARGET_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/CVTermList.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief This class stores a SRM/MRM transition

    The default values for precursor and product m/z values are
    set to numeric_limits<double>::max(). Default values for
    precursor an product charge is set to numeric_limits<Int>::max().
  */
  class OPENMS_DLLAPI IncludeExcludeTarget :
    public CVTermList
  {

public:

    typedef TargetedExperimentHelper::Configuration Configuration;
    typedef TargetedExperimentHelper::RetentionTime RetentionTime;

    /** @name Constructors and destructors
    */
    //@{
    /// default constructor
    IncludeExcludeTarget();

    /// copy constructor
    IncludeExcludeTarget(const IncludeExcludeTarget & rhs);

    /// destructor
    ~IncludeExcludeTarget() override;
    //@}

    /// assignment operator
    IncludeExcludeTarget & operator=(const IncludeExcludeTarget & rhs);

    /** @name Accessors
    */
    //@{
    void setName(const String & name);

    const String & getName() const;

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

    void setProductCVTermList(const CVTermList & list);

    void addProductCVTerm(const CVTerm & cv_term);

    const CVTermList & getProductCVTermList() const;

    void setInterpretations(const std::vector<CVTermList> & interpretations);

    const std::vector<CVTermList> & getInterpretations() const;

    void addInterpretation(const CVTermList & interpretation);

    void setConfigurations(const std::vector<Configuration> & configuration);

    const std::vector<Configuration> & getConfigurations() const;

    void addConfiguration(const Configuration & configuration);

    void setPrediction(const CVTermList & prediction);

    void addPredictionTerm(const CVTerm & prediction);

    const CVTermList & getPrediction() const;

    void setRetentionTime(RetentionTime rt);

    const RetentionTime & getRetentionTime() const;
    //@}

    /** @name Predicates
    */
    //@{
    /// equality operator
    bool operator==(const IncludeExcludeTarget & rhs) const;

    /// inequality operator
    bool operator!=(const IncludeExcludeTarget & rhs) const;
    //@}

protected:

    void updateMembers_();

    String name_;

    double precursor_mz_;

    CVTermList precursor_cv_terms_;

    double product_mz_;

    CVTermList product_cv_terms_;

    std::vector<CVTermList> interpretation_list_;

    String peptide_ref_;

    String compound_ref_;

    std::vector<Configuration> configurations_;

    CVTermList prediction_;

    RetentionTime rts_;

  };
}

#endif // OPENMS_ANALYSIS_TARGETED_INCLUDEEXCLUDETARGET_H
