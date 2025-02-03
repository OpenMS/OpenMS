// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Author: Erhan Kenar $
// $Maintainer: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#define DEBUG_PEAK_PICKING
#undef DEBUG_PEAK_PICKING
//#undef DEBUG_DECONV
namespace OpenMS
{
  class MSChromatogram;
  class OnDiscMSExperiment;

  /**
    @brief This class implements a fast peak-picking algorithm best suited for
    high resolution MS data (FT-ICR-MS, Orbitrap). In high resolution data, the
    signals of ions with similar mass-to-charge ratios (m/z) exhibit little or
    no overlapping and therefore allow for a clear separation. Furthermore, ion
    signals tend to show well-defined peak shapes with narrow peak width.

    This peak-picking algorithm detects ion signals in profile data and
    reconstructs the corresponding peak shape by cubic spline interpolation.
    Signal detection depends on the signal-to-noise ratio which is adjustable
    by the user (see parameter signal_to_noise). A picked peak's m/z and
    intensity value is given by the maximum of the underlying peak spline.

    So far, this peak picker was mainly tested on high resolution data. With
    appropriate preprocessing steps (e.g. noise reduction and baseline
    subtraction), it might be also applied to low resolution data.

    This implementation performs peak picking in a single dimension (m/z);
    two-dimensional data such as ion mobility separated data needs additional
    pre-processing. The current implementation treats these data as
    one-dimensional data, performs peak picking in the m/z dimension and
    reports the intensity weighted ion mobility of the picked peaks (which will
    produce correct results if the data has been binned previously but
    incorrect results if fully 2D data is provided as input).

    @htmlinclude OpenMS_PeakPickerHiRes.parameters

    @note The peaks must be sorted according to ascending m/z!

    @ingroup PeakPicking
  */
  class OPENMS_DLLAPI PeakPickerHiRes :
    public DefaultParamHandler,
    public ProgressLogger
  {
public:
    /// Constructor
    PeakPickerHiRes();

    /// Destructor
    ~PeakPickerHiRes() override;

    /// structure for peak boundaries
    struct PeakBoundary
    {
        double mz_min;
        double mz_max;
    };

    /**
      @brief Applies the peak-picking algorithm to a single spectrum
      (MSSpectrum). The resulting picked peaks are written to the output
      spectrum.
     
      @param input  input spectrum in profile mode
      @param output  output spectrum with picked peaks
     */
    void pick(const MSSpectrum& input, MSSpectrum& output) const;

     /**
      @brief Applies the peak-picking algorithm to a single chromatogram
      (MSChromatogram). The resulting picked peaks are written to the output chromatogram.
     
      @param input  input chromatogram in profile mode
      @param output  output chromatogram with picked peaks
     */
    void pick(const MSChromatogram& input, MSChromatogram& output) const;

    /**
      @brief Applies the peak-picking algorithm to a single spectrum
      (MSSpectrum). The resulting picked peaks are written to the output
      spectrum. Peak boundaries are written to a separate structure.
     
      @param input  input spectrum in profile mode
      @param output  output spectrum with picked peaks
      @param boundaries  boundaries of the picked peaks
      @param check_spacings  check spacing constraints? (yes for spectra, no for chromatograms)
     */
    void pick(const MSSpectrum& input, MSSpectrum& output, std::vector<PeakBoundary>& boundaries, bool check_spacings = true) const;

    /**
      @brief Applies the peak-picking algorithm to a single chromatogram
      (MSChromatogram). The resulting picked peaks are written to the output chromatogram.
     
      @param input  input chromatogram in profile mode
      @param output  output chromatogram with picked peaks
      @param boundaries  boundaries of the picked peaks
      @param check_spacings  check spacing constraints? (yes for spectra, no for chromatograms)
     */
    void pick(const MSChromatogram& input, MSChromatogram& output, std::vector<PeakBoundary>& boundaries, bool check_spacings = false) const;

    /**
      @brief Applies the peak-picking algorithm to a map (MSExperiment). This
      method picks peaks for each scan in the map consecutively. The resulting
      picked peaks are written to the output map.
     
      @param input  input map in profile mode
      @param output  output map with picked peaks
      @param check_spectrum_type  if set, checks spectrum type and throws an exception if a centroided spectrum is passed 
     */
    void pickExperiment(const PeakMap& input, PeakMap& output, const bool check_spectrum_type = true) const;

    /**
      @brief Applies the peak-picking algorithm to a map (MSExperiment). This
      method picks peaks for each scan in the map consecutively. The resulting
      picked peaks are written to the output map.
     
      @param input  input map in profile mode
      @param output  output map with picked peaks
      @param boundaries_spec  boundaries of the picked peaks in spectra
      @param boundaries_chrom  boundaries of the picked peaks in chromatograms
      @param check_spectrum_type  if set, checks spectrum type and throws an exception if a centroided spectrum is passed 
     */
    void pickExperiment(const PeakMap& input,
                        PeakMap& output,
                        std::vector<std::vector<PeakBoundary> >& boundaries_spec,
                        std::vector<std::vector<PeakBoundary> >& boundaries_chrom,
                        const bool check_spectrum_type = true) const;

    /**
      @brief Applies the peak-picking algorithm to a map (MSExperiment). This
      method picks peaks for each scan in the map consecutively. The resulting
      picked peaks are written to the output map.

      Currently we have to give up const-correctness but we know that everything on disc is constant
    */
    void pickExperiment(/* const */ OnDiscMSExperiment& input, PeakMap& output, const bool check_spectrum_type = true) const;

protected:

    template <typename ContainerType>
    void pick_(const ContainerType& input, ContainerType& output, std::vector<PeakBoundary>& boundaries, bool check_spacings = true, int im_index = -1) const;

    // signal-to-noise parameter
    double signal_to_noise_;

    // maximal spacing difference defining a large gap
    double spacing_difference_gap_;
    
    // maximal spacing difference defining a missing data point
    double spacing_difference_;

    // maximum number of missing points
    unsigned missing_;

    // MS levels to which peak picking is applied
    std::vector<Int> ms_levels_;

    /// add floatDataArray 'FWHM'/'FWHM_ppm' to spectra with peak FWHM
    bool report_FWHM_;

    /// unit of 'FWHM' float data array (can be absolute or ppm).
    bool report_FWHM_as_ppm_;

    // docu in base class
    void updateMembers_() override;

  }; // end PeakPickerHiRes

} // namespace OpenMS

