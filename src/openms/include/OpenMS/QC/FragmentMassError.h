// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
    struct Statistics
    {
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
     * @param search_params Input search parameters (corresponding to ID search that generated @param pep_ids) for finding fragment mass tolerance and unit automaticly
     * @param exp Input MSExperiment for MS2 spectra; spectra should be sorted (ascending RT)
     * @param map_to_spectrum Map to find index of spectrum given by meta value at PepID
     * @param tolerance_unit Tolerance in ppm or Dalton (if auto was chosen, the unit and value will taken from FeatureMap metadata)
     * @param tolerance Search window for matching peaks; distance has to be lower than tolerance value (Will be overwritten if tolerance_unit AUTO is chosen)
     * @throws Exceptions::MissingInformation If fragment mass tolerance is missing in @search_params
     * @throws Exception::InvalidParameter PeptideID is missing meta value 'spectrum_reference'
     * @throws Exception::IllegalArgument Spectrum for a PepID has ms-level of 1
     * @throws Exception::MissingInformation If no fragmentation method given in a MS2 precursor
     * @throws Exception::InvalidParameter If the fragmentation method is not ECD, ETD, CID or HCD
     */
    void compute(std::vector<PeptideIdentification>& pep_ids, const ProteinIdentification::SearchParameters& search_params, const MSExperiment& exp, const QCBase::SpectraMap& map_to_spectrum, ToleranceUnit tolerance_unit = ToleranceUnit::AUTO, double tolerance = 20);

    /// returns the name of the metric
    const String& getName() const override;
    
    /// returns results
    const std::vector<Statistics>& getResults() const;


    /**
    * @brief Returns the input data requirements of the compute(...) function
    * @return Status for RAWMZML and POSTFDRFEAT
    */
    QCBase::Status requires() const override;

  private:
    /// container that stores results
    std::vector<Statistics> results_;

    static void calculateFME_(PeptideIdentification& pep_id, const MSExperiment& exp, const QCBase::SpectraMap& map_to_spectrum, bool& print_warning, double tolerance, FragmentMassError::ToleranceUnit tolerance_unit, double& accumulator_ppm, UInt32& counter_ppm, WindowMower& window_mower_filter);

    static void calculateVariance_(FragmentMassError::Statistics& result, const PeptideIdentification& pep_id, const UInt num_ppm);
  };

} //namespace OpenMS
