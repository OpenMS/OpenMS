// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow$
// $Authors: Patricia Scheil, Swenja Wagner$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/QC/QCBase.h>

#include <OpenMS/CONCEPT/Types.h>
#include <string>
#include <vector>

namespace OpenMS
{
  class FeatureMap;
  class MSExperiment;
  class MzTabMetaData;
  class PeptideIdentification;

  /**
   @brief This class is a metric for the QualityControl-ToppTool.

   This class computes the MS2 Identification Rate (as \#identified PSMs divided by total number of MS2 scans) given a FeatureMap and an MSExperiment.
   Only pep-ids with FDR metavalue 'target_decoy' equal to 'target' are counted, unless assume_all_target flag is set (assumes all pep-ids are target peptides)

   */
  class OPENMS_DLLAPI Ms2IdentificationRate : public QCBase
  {
  public:
    /// Structure for storing results
    struct IdentificationRateData {
      Size num_peptide_identification = 0;
      Size num_ms2_spectra = 0;
      double identification_rate = 0.;
    };

  private:
    /// name of the metric
    const String name_ = "Ms2IdentificationRate";

    /// container that stores results
    std::vector<IdentificationRateData> rate_result_;

    /// returns number of all ms2 spectra in an MSExperiment
    Size getMS2Count_(const MSExperiment& exp);

    /*
     * @brief Checks pepID for target/decoy
     *
     * Only checks the first (!) hit, all other hits are ignored
     * Is static so that it can be used with MapUtilities::applyFunctionOnPeptideIDs() without creating a new object for each ID
     *
     * @param id             pepID to be checked
     * @param all_targets    always returns true (if the hits aren't empty)
     * @return               true/false
     * @throws               MissingInformation if target/decoy annotation is missing
     */
    static bool isTargetPeptide_(const PeptideIdentification& id, bool all_targets);

    /*
     * @brief Calculates id-rate and writes the result into a IdentificationRateData object which is appended to rate_result_
     *
     * @param ms2_spectra_count  number of found ms2 spectra
     * @param pep_ids_count      number of found (target) peptide identifications
     */
    void writeResults_(Size pep_ids_count, Size ms2_spectra_count);

  public:
    /// Default constructor
    Ms2IdentificationRate() = default;

    /// Destructor
    virtual ~Ms2IdentificationRate() = default;

    /**
     * @brief computes Ms2 Identification Rate with FeatureMap
     *
     * stores results as a struct in a vector
     * Only pep-ids with target/decoy annotation as 'target' are counted, unless force_index flag is set (assumes all pep-ids are target peptides)
     *
     * @param feature_map       Input FeatureMap with target/decoy annotation
     * @param exp               MSExperiment for counting number of MS2 spectra
     * @param assume_all_target Count all(!) PepIDs towards number of identified MS2 spectra (ignore target/decoy information if any)
     * @exception               MissingInformation is thrown if the mzML is empty
     * @exception               MissingInformation is thrown if the experiment doesn't contain MS2 spectra
     * @exception               Precondition is thrown if there are more identifications than MS2 spectra
     */
    void compute(const FeatureMap& feature_map, const MSExperiment& exp, bool assume_all_target = false);

    /**
     * @brief computes Ms2 Identification Rate with PeptideIdentifications
     *
     * stores results as a struct in a vector
     * Only pep-ids with target/decoy annotation as 'target' are counted, unless force_index flag is set (assumes all pep-ids are target peptides)
     *
     * @param pep_ids           Input PeptideIdentifications with target/decoy annotation
     * @param exp               MSExperiment for counting number of MS2 spectra
     * @param assume_all_target Count all(!) PepIDs towards number of identified MS2 spectra (ignore target/decoy information if any)
     * @exception               MissingInformation is thrown if the mzML is empty
     * @exception               MissingInformation is thrown if the experiment doesn't contain MS2 spectra
     * @exception               Precondition is thrown if there are more identifications than MS2 spectra
     */
    void compute(const std::vector<PeptideIdentification>& pep_ids, const MSExperiment& exp, bool assume_all_target = false);

    /// returns the name of the metric
    const String& getName() const override;

    /// returns results
    const std::vector<IdentificationRateData>& getResults() const;

    /**
     * @brief Returns the input data requirements of the compute(...) function
     * @return Status for RAWMZML and POSTFDRFEAT
     */
    QCBase::Status requirements() const override;

    void addMetaDataMetricsToMzTab(MzTabMetaData& meta) const;
  };

} // namespace OpenMS
