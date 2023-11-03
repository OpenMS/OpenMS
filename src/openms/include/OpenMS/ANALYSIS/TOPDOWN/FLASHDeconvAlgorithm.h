// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <boost/dynamic_bitset.hpp>
#include <iostream>

namespace OpenMS
{
  /**
  @brief FLASHDeconv algorithm: ultrafast mass deconvolution algorithm for top down mass spectrometry dataset
  From MSSpectrum, this class outputs DeconvolvedSpectrum.
  Deconvolution takes three steps:
   i) decharging and select candidate masses - speed up via binning
   ii) collecting isotopes from the candidate masses and deisotope - peak groups are defined here
   iii) scoring and filter out low scoring masses (i.e., peak groups)
  @ingroup Topdown
*/

  class OPENMS_DLLAPI FLASHDeconvAlgorithm : public DefaultParamHandler, public ProgressLogger
  {
  public:
    /// default constructor
    FLASHDeconvAlgorithm();

    /// copy constructor

    /// move constructor
    FLASHDeconvAlgorithm(FLASHDeconvAlgorithm&& other) = default;

    /// assignment operator
    FLASHDeconvAlgorithm& operator=(const FLASHDeconvAlgorithm& fd) = default;

    /// move assignment operator
    FLASHDeconvAlgorithm& operator=(FLASHDeconvAlgorithm&& fd) = default;

    /// destructor
    ~FLASHDeconvAlgorithm() = default;

    void run(MSExperiment& map, std::vector<DeconvolvedSpectrum>& deconvolved_spectra, std::vector<FLASHDeconvHelperStructs::MassFeature>& deconvolved_feature);

    /// get calculated averagine. Call after calculateAveragine is called.
    const FLASHDeconvHelperStructs::PrecalculatedAveragine& getAveragine();

  protected:
    void updateMembers_() override;

  private:
    SpectralDeconvolution sd_, sd_charge_decoy_, sd_noise_decoy_, sd_isotope_decoy_;

    int merge_spec_ = 0;

    int forced_ms_level_ = 0;
    UInt max_ms_level_ = 0;

    UInt current_max_ms_level_ = 0;
    UInt current_min_ms_level_ = 0;

    int preceding_MS1_count_ = 0;

    String ida_log_file_;

    DoubleList tols_, min_cos_;
    bool use_RNA_averagine_ = false;
    bool report_decoy_ = false;

    /// default precursor isolation window size.
    double isolation_window_size_;

    std::map<int, std::vector<std::vector<float>>> precursor_map_for_ida_;
    std::map<int, PeakGroup> ms2scan_to_precursor_peak_group_map_; // MS2 scan number, peak group

    int scanMap_(MSExperiment& map);

    void mergeSpectra_(MSExperiment& map, uint ms_level) const;

    int runFD_(MSExperiment& map, std::vector<DeconvolvedSpectrum>& deconvolved_spectra);

    void runFeatureFinding_(std::vector<DeconvolvedSpectrum>& deconvolved_spectra, std::vector<FLASHDeconvHelperStructs::MassFeature>& deconvolved_features);

    void updatePrecursorQScores_(std::vector<DeconvolvedSpectrum>& deconvolved_spectra);

    void registerPrecursorFromIda_(std::map<int, PeakGroup>& precursor_peak_group_map, int scan_number, double start_mz, double end_mz);

    static int getScanNumber_(const MSExperiment& map, Size index);
    /**
    @brief register the precursor peak as well as the precursor peak group (or mass) if possible for MSn (n>1) spectrum.
    Given a precursor peak (found in the original MS n-1 Spectrum) the masses containing the precursor peak are searched.
    If multiple masses are detected, the one with the best setQscore is selected. For the selected mass, its corresponding peak group (along with precursor peak) is registered.
    If no such mass exists, only the precursor peak is registered.
    @param survey_scans the candidate precursor spectra - the user may allow search of previous N survey scans.
    @param precursor_map_for_real_time_acquisition this contains the deconvolved mass information from FLASHIda runs.
    */
    void findPrecursorPeakGroups_(std::map<int, PeakGroup>& precursor_peak_group_map, const MSExperiment& map, const std::vector<DeconvolvedSpectrum>& deconvolved_spectra, uint ms_level);

    static void filterLowPeaks_(MSExperiment& map, Size count);
  };
} // namespace OpenMS