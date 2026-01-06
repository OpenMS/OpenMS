// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/config.h>
#include <iomanip>

namespace OpenMS
{
  /**
    @brief FLASHDeconv Spectrum level output *.tsv, *.msalign (for TopPIC) file formats
     @ingroup FileIO
**/

  class OPENMS_DLLAPI FLASHDeconvSpectrumFile
  {
  public:
    /**
            @brief write the header in the tsv output file (spectrum level)
            @param[out] os output stream to the output file
            @param[in] ms_level ms level of the spectrum
            @param[out] detail if set true, detailed information of the mass (e.g., peak list for the mass) is written
            @param[out] report_decoy if set true, decoy and qvalue information will be written.
       */
    static void writeDeconvolvedMassesHeader(std::ostream& os,
                                             uint ms_level,
                                             bool detail,
                                             bool report_decoy);
    /**
          @brief write the deconvolved masses in the output file (spectrum level)
          @param[in] dspec deconvolved spectrum to write
          @param[out] fs file stream to the output file
          @param[out] file_name the output file name that the deconvolved masses will be written.
          @param[in] avg averagine information to calculate monoisotopic and average mass difference within this function. In PeakGroup (peaks of DeconvolvedSpectrum) only monoisotopic mass is recorded. To write both monoisotopic and average masses, their mass difference should be calculated using this averagine information.
          @param[in] decoy_avg averagine for noise decoy
          @param[in] tol mass tolerance
          @param[out] write_detail if this is set, more detailed information on each mass will be written in the output file.
          @param[out] record_decoy if set true, decoy and qvalue information will be written.
          @param[out] noise_decoy_weight noise decoy weight. Determines how often the noise decoy masses will be written
          Default MS1 headers are:
            FileName, ScanNum, TargetDecoyType, RetentionTime, MassCountInSpec, AverageMass, MonoisotopicMass,
            SumIntensity, MinCharge, MaxCharge,
            PeakCount, IsotopeCosine, ChargeScore, MassSNR, ChargeSNR, RepresentativeCharge, RepresentativeMzStart, RepresentativeMzEnd, setQscore, PerChargeIntensity, PerIsotopeIntensity

      Default MS2 headers include MS1 headers plus:
        PrecursorScanNum, PrecursorMz, PrecursorIntensity, PrecursorCharge, PrecursorSNR, PrecursorMonoisotopicMass, PrecursorQscore

      Detailed MS1 and MS2 headers include all corresponding headers above plus:
        PeakMZs, PeakIntensities, PeakCharges, PeakMasses, PeakIsotopeIndices, PeakPPMErrors
    */
    static void writeDeconvolvedMasses(const DeconvolvedSpectrum& dspec,
                                       std::ostream& os,
                                       const String& file_name,
                                       const FLASHHelperClasses::PrecalculatedAveragine& avg,
                                       const FLASHHelperClasses::PrecalculatedAveragine& decoy_avg,
                                       double tol,
                                       bool write_detail,
                                       bool record_decoy, double noise_decoy_weight);

    /**
     *
     * @param[in] map
     * @param[in] deconvolved_spectra
     * @param[in] deconvolved_mzML_file
     * @param[in] annotated_mzML_file
     * @param[in] mzml_charge
     * @param[in] tols
     */
    static void writeMzML(const MSExperiment& map,
                          std::vector<DeconvolvedSpectrum>& deconvolved_spectra,
                          const String& deconvolved_mzML_file,
                          const String& annotated_mzML_file,
                          int mzml_charge, DoubleList tols);

    /**
     * write isobaric quantification results
     * @param[out] os output stream
     * @param[in] deconvolved_spectra deconvolved spectra to write
     */
    static void writeIsobaricQuantification(std::ostream& os, std::vector<DeconvolvedSpectrum>& deconvolved_spectra);

    /// write topFD header
    static void writeTopFDHeader(std::ostream& os, const Param& param);

    /**
      @brief write the deconvolved masses TopFD output (*.msalign)
      @param[in] dspec deconvolved spectrum to write
      @param[out] os output stream to the output file
      @param[in] filename mzml file name
      @param[in] qval_threshold qvalue threshold to filter out high qvalue precursors.
      @param[in] min_ms_level min ms level of the dataset
      @param[in] randomize_precursor_mass if set, a random number between -100 to 100 is added to precursor mass
      @param[in] randomize_fragment_mass if set, a random number between -100 to 100 is added to fragment mass
    */
    static void writeTopFD(const DeconvolvedSpectrum& dspec, std::ostream& os, const String& filename,
                           double qval_threshold = 1.0,
                           uint min_ms_level = 1,
                           bool randomize_precursor_mass = false,
                           bool randomize_fragment_mass = false);

  private:

    /// number of minimum peak count in topFD msalign file
    static const int topFD_min_peak_count_ = 3;
    /// number of maximum peak count in topFD msalign file
    static const int topFD_max_peak_count_ = 500;
  };
}// namespace OpenMS
