// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>

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
      @param fs file stream to the output file
      @param ms_level ms level of the spectrum
      @param detail if set true, detailed information of the mass (e.g., peak list for the mass) is written
      @param dummy if set true, dummy and qvalue information will be written.
       */
    static void writeDeconvolvedMassesHeader(std::fstream& fs,
                                             uint ms_level,
                                             bool detail,
                                             bool dummy);

    /**
      @brief write the deconvolved masses in the output file (spectrum level)
      @param dspec deconvolved spectrum to write
      @param target_spec target spectrum only used for dummy spectrum output
      @param fs file stream to the output file
      @param file_name the output file name that the deconvolved masses will be written.
      @param avg averagine information to calculate monoisotopic and average mass difference within this function. In PeakGroup (peaks of DeconvolvedSpectrum) only monoisotopic mass is recorded. To write both monoisotopic and average masses, their mass difference should be calculated using this averagine information.
      @param tol mass tolerance
      @param write_detail if this is set, more detailed information on each mass will be written in the output file.
      @param record_dummy if set true, dummy and qvalue information will be written.
      Default MS1 headers are:
        FileName, ScanNum, TargetDummyType, RetentionTime, MassCountInSpec, AverageMass, MonoisotopicMass,
        SumIntensity, MinCharge, MaxCharge,
        PeakCount, IsotopeCosine, ChargeScore, MassSNR, ChargeSNR, RepresentativeCharge, RepresentativeMzStart, RepresentativeMzEnd, Qscore, PerChargeIntensity, PerIsotopeIntensity

      Default MS2 headers include MS1 headers plus:
        PrecursorScanNum, PrecursorMz, PrecursorIntensity, PrecursorCharge, PrecursorSNR, PrecursorMonoisotopicMass, PrecursorQscore

      Detailed MS1 and MS2 headers include all corresponding headers above plus:
        PeakMZs, PeakIntensities, PeakCharges, PeakMasses, PeakIsotopeIndices, PeakPPMErrors
    */
    static void writeDeconvolvedMasses(DeconvolvedSpectrum& dspec,
                                       DeconvolvedSpectrum& target_spec,
                                       std::fstream& fs,
                                       const String& file_name,
                                       const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg,
                                       double tol,
                                       bool write_detail,
                                       bool record_dummy);

    /**
      @brief write the deconvolved masses TopFD output (*.msalign)
      @param dspec deconvolved spectrum to write
      @param fs file stream to the output file
      @param snr_threshold SNR threshold to filter out low SNR precursors. Even if a PeakGroup has a high deconvolution quality, it should be still discarded for identification when its precursor SNR (SNR within the isolation window) is too low.
      @param min_ms_level min ms level of the dataset
      @param randomize_precursor_mass if set, a random number between -100 to 100 is added to precursor mass
      @param randomize_fragment_mass if set, a random number between -100 to 100 is added to fragment mass
    */
    //      @param avg averagine information to calculate monoisotopic and average mass difference
    static void writeTopFD(DeconvolvedSpectrum& dspec, std::fstream& fs,
                           double snr_threshold = 1.0,
                           const uint min_ms_level = 1,
                           bool randomize_precursor_mass = false,
                           bool randomize_fragment_mass = false);

  private:

    /// number of minimum peak count in topFD msalign file
    static const int topFD_min_peak_count_ = 3;
    /// number of maximum peak count in topFD msalign file
    static const int topFD_max_peak_count_ = 500;

  };
}// namespace OpenMS