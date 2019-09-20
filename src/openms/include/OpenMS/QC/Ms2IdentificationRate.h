// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Chris Bielow$
// $Authors: Patricia Scheil, Swenja Wagner$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include "OpenMS/QC/QCBase.h"
#include <string>
#include <vector>

namespace OpenMS
{
  class FeatureMap;
  class MSExperiment;
  /**
   @brief This class is a metric for the QualityControl-ToppTool.
   
   This class computes the MS2 Identification Rate (as #identified PSMs divided by total number of MS2 scans) given a FeatureMap and an MSExperiment.
   Only pep-ids with FDR metavalue 'target_decoy' equal to 'target' are counted, unless force_fdr flag is set (assumes all pep-ids are target peptides)

   */
  class OPENMS_DLLAPI Ms2IdentificationRate : public QCBase
  {
  public:
    /// Structure for storing results
    struct IdentificationRateData
    {
      UInt64 num_peptide_identification;
      UInt64 num_ms2_spectra;
      double identification_rate;
    };

  private:
    /// name of the metric
    const String name_ = "Ms2IdentificationRate";
    
    /// container that stores results
    std::vector<IdentificationRateData> rate_result_;

  public:
    /// Default constructor
    Ms2IdentificationRate() = default;

    /// Destructor
    virtual ~Ms2IdentificationRate() = default;

    /**
     * @brief computes Ms2 Identification Rate
     *
     * stores results as a struct in a vector
     * Only pep-ids with FDR metavalue annotation as 'target' are counted, unless force_fdr flag is set (assumes all pep-ids are target peptides)
     *
     * @param feature_map Input FeatureMap with target/decoy annotation
     * @param exp MSExperiment for counting number of MS2 spectra
     * @param force_fdr Count all(!) PepIDs towards number of identified MS2 spectra (ignore target/decoy information if any)
     * @exception Exception::MissingInformation is thrown if the FeatureXML is empty
     * @exception Exception::MissingInformation is thrown if the mzML is empty
     * @exception Exception::MissingInformation is thrown if the experiment doesn't contain MS2 spectra
     * @exception Exception::Precondition is thrown if there are more identifications than MS2 spectra
     */
    void compute(const FeatureMap& feature_map, const MSExperiment& exp, bool force_fdr = false);

    /// returns the name of the metric
    const String& getName() const override;
    
    /// returns results
    const std::vector<IdentificationRateData>& getResults() const;

    /**
     * @brief Returns the input data requirements of the compute(...) function
     * @return Status for RAWMZML and POSTFDRFEAT
     */
    QCBase::Status requires() const override;
  };

} // namespace OpenMS
