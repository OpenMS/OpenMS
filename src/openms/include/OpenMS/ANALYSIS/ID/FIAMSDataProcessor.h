// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Svetlana Kutuzova, Douglas McCloskey $
// $Authors: Svetlana Kutuzova, Douglas McCloskey $
// --------------------------------------------------------------------------
 
#pragma once
#include <OpenMS/PROCESSING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/FORMAT/MzTabM.h>
#include <OpenMS/FORMAT/MzTabFile.h>

namespace OpenMS
{
  /**
    @brief Data processing for FIA-MS data
    
    Flow injection analysis (FIA) omits the separation step by removal of the column.
    It allows for much faster processing time with the cost of ambiguities in the data 
    interpretation. The compounds are identified through the accurate mass search.

    Flow injection analysis class implements the basic FIA-MS data processing steps such as: 
    acquiring the data for the certain time interval, summing along the time axis, smoothing 
    the peaks, peak picking and accurate mass search. The batch runs are to be managed with 
    the FIAMSSchedule class that takes a simple csv file as an input.

    The workflow is inspired by the data processing from Fuhrer et al https://pubs.acs.org/doi/10.1021/ac201267k
    though it is not the exact implementation.
  */
  class OPENMS_DLLAPI FIAMSDataProcessor  :
    public DefaultParamHandler
  {
public:
    /// Constructor
    FIAMSDataProcessor();

    /// Default destructor
    ~FIAMSDataProcessor() override = default;

    /// Copy constructor
    FIAMSDataProcessor(const FIAMSDataProcessor& cp) = default;

    /// Assignment
    FIAMSDataProcessor& operator=(const FIAMSDataProcessor& fdp) = default;

    /**
      @brief Run the full analysis for the experiment for the given time interval

      The workflow steps are:
      - the time axis of the experiment is cut to the interval from 0 to n_seconds
      - the spectra are summed into one along the time axis with the bin size determined by mz and instrument resolution
      - data is smoothed by applying the Savitzky-Golay filter
      - peaks are picked
      - the accurate mass search for all the picked peaks is performed

      The intermediate summed spectra and picked peaks can be saved to the filesystem. 
      Also, the results of the accurate mass search and the signal-to-noise information 
      of the resulting spectrum is saved.

      @param experiment  Input MSExperiment
      @param n_seconds Input number of seconds
      @param load_cached_spectrum Load the cached picked spectrum if exists
      @param[out] output Output of the accurate mass search results
      @return a boolean indicating if the picked spectrum was loaded from the cached file
    */
    bool run(const MSExperiment& experiment, const float n_seconds, OpenMS::MzTab& output, const bool load_cached_spectrum = true);

    /**
      @brief Cut the time axis of the experiment from 0 to @p n_seconds

      @param experiment  Input MSExperiment
      @param n_seconds Input number of seconds
      @param output   [out] Spectra with retention time less than @p n_seconds
    */
    void cutForTime(const MSExperiment& experiment, const float n_seconds, std::vector<MSSpectrum>& output);

    /**
      @brief Sum the spectra with different retention times into one.

      The bin size for summing the intensities is defined as mz / (resolution*4) 
      for all the mzs taken with the @p bin_step defined in the parameters.
      Uses `SpectrumAddition::addUpSpectra` function with the sliding bin size parameter. 

      @param input  Input vector of spectra
      @return a spectrum
    */
    MSSpectrum mergeAlongTime(const std::vector<OpenMS::MSSpectrum>& input);

    /**
      @brief Pick peaks from the summed spectrum

      @param input  Input vector of spectra
      @return a spectrum with picked peaks
    */
    MSSpectrum extractPeaks(const MSSpectrum& input);

    /**
      @brief Convert a spectrum to a feature map with the corresponding polarity

      Applies `SavitzkyGolayFilter` and `PeakPickerHiRes`

      @param input  Input a picked spectrum
      @return a feature map with the peaks converted to features and polarity from the parameters
    */
    FeatureMap convertToFeatureMap(const MSSpectrum& input);

    /**
      @brief Estimate noise for each peak

      Uses `SignalToNoiseEstimatorMedianRapid`

      @param input  Input a picked spectrum
      @return a spectrum object storing logSN information
    */
    MSSpectrum trackNoise(const MSSpectrum& input);

    /**
      @brief Perform accurate mass search

      Uses `AccurateMassSearchEngine`

      @param input  Input a feature map
      @param output  [out] mzTab file with the accurate mass search results
    */
    void runAccurateMassSearch(FeatureMap& input, OpenMS::MzTab& output);

    /**
      @brief Get mass-to-charge ratios to base the summing the spectra along the time axis upon
    */
    const std::vector<float>& getMZs();

    /**
      @brief Get the sliding bin sizes for summing the spectra along the time axis
    */
    const std::vector<float>& getBinSizes();

protected:
    void updateMembers_() override;

private:
    /**
      @brief Store the spectrum to the given filepath
    */
    void storeSpectrum_(const MSSpectrum& input, const String& filename);

    std::vector<float> mzs_; 
    std::vector<float> bin_sizes_;
    SavitzkyGolayFilter sgfilter_;
    PeakPickerHiRes picker_;
  };
} // namespace OpenMS
