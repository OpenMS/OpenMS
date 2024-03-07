// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
    FLASHDeconvAlgorithm(const FLASHDeconvAlgorithm&) = default;

    /// move constructor
    FLASHDeconvAlgorithm(FLASHDeconvAlgorithm&& other) = default;

    /// assignment operator
    FLASHDeconvAlgorithm& operator=(const FLASHDeconvAlgorithm& fd) = default;

    /// move assignment operator
    FLASHDeconvAlgorithm& operator=(FLASHDeconvAlgorithm&& fd) = default;

    /// destructor
    ~FLASHDeconvAlgorithm() = default;

    /**
     * @brief Run FLASHDeconv algorithm for @p map and store @p deconvolved_spectra and @p deconvolved_feature
     * @param map the dataset
     * @param deconvolved_spectra the deconvolved spectra will be stored in here
     * @param deconvolved_feature tje deconvolved features wll be strored in here
     */
    void run(MSExperiment& map, std::vector<DeconvolvedSpectrum>& deconvolved_spectra, std::vector<FLASHDeconvHelperStructs::MassFeature>& deconvolved_feature);

    /// get calculated averagine. Call after calculateAveragine is called.
    const FLASHDeconvHelperStructs::PrecalculatedAveragine& getAveragine();

    /// get noise decoy weight
    double getNoiseDecoyWeight()
    {
      return noise_decoy_weight_;
    }

  protected:
    void updateMembers_() override;

  private:
    /// SpectralDeconvolution  instances for spectral deconvolution for target and decoys
    SpectralDeconvolution sd_, sd_charge_decoy_, sd_noise_decoy_, sd_isotope_decoy_;

    /// to merge spectra.
    int merge_spec_ = 0;

    /// forced MS level
    int forced_ms_level_ = 0;

    /// maximum MS level, which is 4.
    UInt max_ms_level_ = 4;

    /// current maximum MS level - i.e., for MS2, this is the precursor charge
    UInt current_max_ms_level_ = 0;

    /// currment minimum MS level.
    UInt current_min_ms_level_ = 0;

    /// the number of preceding full scans from which MS2 precursor mass will be searched.
    int preceding_MS1_count_ = 0;

    /// FLASHIda log file name
    String ida_log_file_;

    /// mass tolerances, and minimum cosine scores per MS level
    DoubleList tols_, min_cos_;

    /// to use RNA averagine model
    bool use_RNA_averagine_ = false;

    /// should decoy deconvolution be done?
    bool report_decoy_ = false;

    /// default precursor isolation window size.
    double isolation_window_size_;

    /// noise decoy weight determined with qvalue calcualtion.
    double noise_decoy_weight_ = 1;
    /// FLASHIda parsing information is stored here: MS1 scan - information
    std::map<int, std::vector<std::vector<float>>> precursor_map_for_ida_;
    /// a map from native ID to precursor peak group
    std::map<String, PeakGroup> native_id_precursor_peak_group_map_;

    /// read dataset to update ms level information
    void updateMSLevels_(MSExperiment& map);

    /// merge spectra
    void mergeSpectra_(MSExperiment& map, uint ms_level);

    /// run spectral deconvolution
    void runSpectralDeconvolution_(MSExperiment& map, std::vector<DeconvolvedSpectrum>& deconvolved_spectra);

    /// run feature finding to get deconvolved features
    void runFeatureFinding_(std::vector<DeconvolvedSpectrum>& deconvolved_spectra, std::vector<FLASHDeconvHelperStructs::MassFeature>& deconvolved_features);

    /// with found deconvolved features, update QScores for masses that are contained in features.
    static void updatePrecursorQScores_(std::vector<DeconvolvedSpectrum>& deconvolved_spectra, int ms_level);

    /// find precursor peak groups from FLASHIda log file
    void findPrecursorPeakGroupsFormIdaLog_(const MSExperiment& map, Size index, double start_mz, double end_mz);

    /// register the precursor peak group (or mass) if possible for MSn (n>1) spectrum.
    void findPrecursorPeakGroupsForMSnSpectra_(const MSExperiment& map, const std::vector<DeconvolvedSpectrum>& deconvolved_spectra, uint ms_level);

    /// get scan number of the spectrum in @p index -th in @p map
    static int getScanNumber_(const MSExperiment& map, Size index);

    /// filter low intensity peaks
    static void filterLowPeaks_(MSExperiment& map);
  };
} // namespace OpenMS
