// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/ChromatogramPeak.h>

#include <OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
#include <OpenMS/PROCESSING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/PROCESSING/SMOOTHING/GaussFilter.h>

#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>

#ifdef WITH_CRAWDAD
#include <CrawdadWrapper.h>
#endif

namespace OpenMS
{

  /**
  @brief The PeakPickerChromatogram finds peaks a single chromatogram.

  @htmlinclude OpenMS_PeakPickerChromatogram.parameters

  It uses the PeakPickerHiRes internally to find interesting seed candidates.
  These candidates are then expanded and a right/left border of the peak is
  searched.
  Additionally, overlapping peaks can be removed.

  The picked chromatogram's FloatDataArrays report, per peak (indexed via the
  FLOATINDICES enum): "FWHM", "IntegratedIntensity", "leftWidth", "rightWidth",
  and "SN" (the apex signal-to-noise ratio, as estimated by the internal
  SignalToNoiseEstimatorMedian). The "SN" array is only populated with real
  values if the 'report_sn' parameter is enabled or 'signal_to_noise' > 0;
  otherwise it is filled with -1 for every peak.

  */

  class OPENMS_DLLAPI PeakPickerChromatogram :
    public DefaultParamHandler
  {

public:

    //@{
    /// Constructor
    PeakPickerChromatogram();

    /// Destructor
    ~PeakPickerChromatogram() override {}
    //@}

	/// indices into FloatDataArrays of resulting picked chromatograms
	enum FLOATINDICES { IDX_FWHM = 0, IDX_ABUNDANCE = 1, IDX_LEFTBORDER = 2, IDX_RIGHTBORDER = 3, IDX_SN = 4, SIZE_OF_FLOATINDICES };

    /**
      @brief Finds peaks in a single chromatogram and annotates left/right borders

      It uses a modified algorithm of the PeakPickerHiRes

      This function will return a picked chromatogram
    */
    void pickChromatogram(const MSChromatogram& chromatogram, MSChromatogram& picked_chrom);
    

    /**
      @brief Finds peaks in a single chromatogram and annotates left/right borders

      It uses a modified algorithm of the PeakPickerHiRes

      This function will return a picked chromatogram and a smoothed chromatogram
    */
    void pickChromatogram(const MSChromatogram& chromatogram, MSChromatogram& picked_chrom, MSChromatogram& smoothed_chrom);

protected:

    void pickChromatogramCrawdad_(const MSChromatogram& chromatogram, MSChromatogram& picked_chrom);

    void pickChromatogram_(const MSChromatogram& chromatogram, MSChromatogram& picked_chrom);

    /**
      @brief Compute peak area (peak integration)
    */
    void integratePeaks_(const MSChromatogram& chromatogram);

    /**
      @brief Helper function to find the closest peak in a chromatogram to "target_rt" 

      The search will start from the index current_peak, so the function is
      assuming the closest peak is to the right of current_peak.

      It will return the index of the closest peak in the chromatogram.
    */
    Size findClosestPeak_(const MSChromatogram& chromatogram, double target_rt, Size current_peak = 0);

    /**
      @brief Helper function to remove overlapping peaks in a single Chromatogram

    */
    void removeOverlappingPeaks_(const MSChromatogram& chromatogram, MSChromatogram& picked_chrom);


    /// Synchronize members with param class
    void updateMembers_() override;

    /// Assignment operator is protected for algorithm
    PeakPickerChromatogram& operator=(const PeakPickerChromatogram& rhs);

    // Members
    /// Frame length for the SGolay smoothing
    UInt sgolay_frame_length_;
    /// Polynomial order for the SGolay smoothing
    UInt sgolay_polynomial_order_;
    /// Width of the Gaussian smoothing
    double gauss_width_;
    /// Whether to use Gaussian smoothing
    bool use_gauss_;
    /// Whether to resolve overlapping peaks 
    bool remove_overlapping_;

    /// Forced peak with
    double peak_width_;
    /// Signal to noise threshold
    double signal_to_noise_;

    /// Signal to noise window length
    double sn_win_len_;
    /// Signal to noise bin count
    UInt sn_bin_count_;
    /// Whether to write out log messages of the SN estimator
    bool write_sn_log_messages_;
    /// Whether to report the apex signal-to-noise ratio of each picked peak (forces
    /// initialization of the S/N estimator even if 'signal_to_noise' is 0)
    bool report_sn_;
    /// Peak picker method
    std::string method_;

    /// Temporary vector to hold the integrated intensities
    std::vector<double> integrated_intensities_;
    /// Temporary vector to hold the peak left widths
    std::vector<int> left_width_;
    /// Temporary vector to hold the peak right widths
    std::vector<int> right_width_;
    /// Temporary vector to hold the apex signal-to-noise ratios
    std::vector<double> apex_sn_;

    PeakPickerHiRes pp_;
    SavitzkyGolayFilter sgolay_;
    GaussFilter gauss_;
    SignalToNoiseEstimatorMedian<MSChromatogram > snt_;
  };
}


