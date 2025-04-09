// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow$
// $Authors: Patricia Scheil, Swenja Wagner$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/QC/QCBase.h>
#include <vector>

namespace OpenMS
{
  class FeatureMap;
  class MSExperiment;
  class PeptideIdentification;
  class WindowMower;

  class OPENMS_DLLAPI FragmentMassError : public QCBase
  {
  public:
    /// Default constructor
    FragmentMassError() = default;

    /// Destructor
    virtual ~FragmentMassError() = default;

    /**
     * @brief Structure for storing results: average and variance of all FragmentMassErrors in ppm
     */
    struct Statistics {
      double average_ppm = 0;
      double variance_ppm = 0;
    };

    /**
     * @brief computes FragmentMassError (FME) in ppm and Dalton (only of the first PeptideHit of each PepID)
     *
     * Stores average FME over all spectra (one for each PeptideIdentification) and its variance in ppm as a struct in a vector.
     * Each FME (in ppm) is stored at the first PeptideHit of the corresponding PeptideIdentification as metavalue Constants::UserParam::FRAGMENT_ERROR_PPM_METAVALUE_USERPARAM
     * and contains the FME for each peak in the corresponding spectrum.
     * Same is done for the FME in Da - as metavalue Constants::UserParam::FRAGMENT_ERROR_DA_METAVALUE_USERPARAM.
     * For both tolerance units the variance of FMEs over the spectrum is also stored as a metavalue with the extension "_variance" to the metavalue name.
     * Note: Variance will not be written if 1 or less FMEs were calculated.
     * Note: If the metavalues already exist, they will be overwritten.
     *
     * @param fmap Input FeatureMap for annotation and data for theoretical spectra
     * @param exp Input MSExperiment for MS2 spectra; spectra should be sorted (ascending RT)
     * @param map_to_spectrum Map to find index of spectrum given by meta value at PepID
     * @param tolerance_unit Tolerance in ppm or Dalton (if auto was chosen, the unit and value will taken from FeatureMap metadata)
     * @param tolerance Search window for matching peaks; distance has to be lower than tolerance value (Will be overwritten if tolerance_unit AUTO is chosen)
     * @throws Exceptions::MissingInformation If fragment mass tolerance is missing in metadata of FeatureMap
     * @throws Exception::InvalidParameter PeptideID is missing meta value 'spectrum_reference'
     * @throws Exception::IllegalArgument Spectrum for a PepID has ms-level of 1
     * @throws Exception::MissingInformation If no fragmentation method given in a MS2 precursor
     * @throws Exception::InvalidParameter If the fragmentation method is not ECD, ETD, CID or HCD
     */
    void compute(FeatureMap& fmap, const MSExperiment& exp, const QCBase::SpectraMap& map_to_spectrum, ToleranceUnit tolerance_unit = ToleranceUnit::AUTO, double tolerance = 20);

    /**
     * @brief computes FragmentMassError (FME) in ppm and Dalton (only of the first PeptideHit of each PepID)
     *
     * Stores average FME over all spectra and its variance in ppm as a struct in a vector.
     * Each FME (in ppm) is stored at the first PeptideHit of the corresponding PeptideIdentification as metavalue Constants::UserParam::FRAGMENT_ERROR_PPM_METAVALUE_USERPARAM
     * and contains the FME for each peak in the corresponding spectrum.
     * Same is done for the FME in Da - as metavalue Constants::UserParam::FRAGMENT_ERROR_DA_METAVALUE_USERPARAM.
     * For both tolerance units the variance of FMEs over the spectrum is also stored as a metavalue with the extension "_variance" to the metavalue name.
     * Note: Variance will not be written if 1 or less FMEs were calculated.
     * Note: If the metavalues already exist, they will be overwritten.
     *
     * @param pep_ids Input vector of peptide identifications for annotation and data for theoretical spectra
     * @param search_params Input search parameters (corresponding to ID search that generated @p pep_ids) for finding fragment mass tolerance and unit automatically
     * @param exp Input MSExperiment for MS2 spectra; spectra should be sorted (ascending RT)
     * @param map_to_spectrum Map to find index of spectrum given by meta value at PepID
     * @param tolerance_unit Tolerance in ppm or Dalton (if auto was chosen, the unit and value will taken from FeatureMap metadata)
     * @param tolerance Search window for matching peaks; distance has to be lower than tolerance value (Will be overwritten if tolerance_unit AUTO is chosen)
     * @throws Exceptions::MissingInformation If fragment mass tolerance is missing in @p search_params
     * @throws Exception::InvalidParameter PeptideID is missing meta value 'spectrum_reference'
     * @throws Exception::IllegalArgument Spectrum for a PepID has ms-level of 1
     * @throws Exception::MissingInformation If no fragmentation method given in a MS2 precursor
     * @throws Exception::InvalidParameter If the fragmentation method is not ECD, ETD, CID or HCD
     */
    void compute(std::vector<PeptideIdentification>& pep_ids, const ProteinIdentification::SearchParameters& search_params, const MSExperiment& exp, const QCBase::SpectraMap& map_to_spectrum,
                 ToleranceUnit tolerance_unit = ToleranceUnit::AUTO, double tolerance = 20);

    /// returns the name of the metric
    const String& getName() const override;

    /// returns results
    const std::vector<Statistics>& getResults() const;


    /**
     * @brief Returns the input data requirements of the compute(...) function
     * @return Status for RAWMZML and POSTFDRFEAT
     */
    QCBase::Status requirements() const override;

  private:
    /// container that stores results
    std::vector<Statistics> results_;

    static void calculateFME_(PeptideIdentification& pep_id, const MSExperiment& exp, const QCBase::SpectraMap& map_to_spectrum, bool& print_warning, double tolerance,
                              FragmentMassError::ToleranceUnit tolerance_unit, double& accumulator_ppm, UInt32& counter_ppm, WindowMower& window_mower_filter);

    static void calculateVariance_(FragmentMassError::Statistics& result, const PeptideIdentification& pep_id, const UInt num_ppm);
  };

} // namespace OpenMS
