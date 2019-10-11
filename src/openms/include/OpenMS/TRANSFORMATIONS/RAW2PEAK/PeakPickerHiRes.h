// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
     * @brief Applies the peak-picking algorithm to a single spectrum
     * (MSSpectrum). The resulting picked peaks are written to the output
     * spectrum.
     *
     * @param input  input spectrum in profile mode
     * @param output  output spectrum with picked peaks
     */
    void pick(const MSSpectrum& input, MSSpectrum& output) const;

     /**
     * @brief Applies the peak-picking algorithm to a single chromatogram
     * (MSChromatogram). The resulting picked peaks are written to the output chromatogram.
     *
     * @param input  input chromatogram in profile mode
     * @param output  output chromatogram with picked peaks
     */
    void pick(const MSChromatogram& input, MSChromatogram& output) const;

    /**
     * @brief Applies the peak-picking algorithm to a single spectrum
     * (MSSpectrum). The resulting picked peaks are written to the output
     * spectrum. Peak boundaries are written to a separate structure.
     *
     * @param input  input spectrum in profile mode
     * @param output  output spectrum with picked peaks
     * @param boundaries  boundaries of the picked peaks
     * @param check_spacings  check spacing constraints? (yes for spectra, no for chromatograms)
     */
    void pick(const MSSpectrum& input, MSSpectrum& output, std::vector<PeakBoundary>& boundaries, bool check_spacings = true) const;

    /**
     * @brief Applies the peak-picking algorithm to a single chromatogram
     * (MSChromatogram). The resulting picked peaks are written to the output chromatogram.
     *
     * @param input  input chromatogram in profile mode
     * @param output  output chromatogram with picked peaks
     * @param boundaries  boundaries of the picked peaks
     */
    void pick(const MSChromatogram& input, MSChromatogram& output, std::vector<PeakBoundary>& boundaries, bool check_spacings = false) const;

    /**
     * @brief Applies the peak-picking algorithm to a map (MSExperiment). This
     * method picks peaks for each scan in the map consecutively. The resulting
     * picked peaks are written to the output map.
     *
     * @param input  input map in profile mode
     * @param output  output map with picked peaks
     * @param check_spectrum_type  if set, checks spectrum type and throws an exception if a centroided spectrum is passed 
     */
    void pickExperiment(const PeakMap& input, PeakMap& output, const bool check_spectrum_type = true) const;

    /**
     * @brief Applies the peak-picking algorithm to a map (MSExperiment). This
     * method picks peaks for each scan in the map consecutively. The resulting
     * picked peaks are written to the output map.
     *
     * @param input  input map in profile mode
     * @param output  output map with picked peaks
     * @param boundaries_spec  boundaries of the picked peaks in spectra
     * @param boundaries_chrom  boundaries of the picked peaks in chromatograms
     * @param check_spectrum_type  if set, checks spectrum type and throws an exception if a centroided spectrum is passed 
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
    void pick_(const ContainerType& input, ContainerType& output, std::vector<PeakBoundary>& boundaries, bool check_spacings = true) const;

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

