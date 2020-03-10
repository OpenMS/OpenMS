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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <vector>

namespace OpenMS
{

  /**
    @brief This class implements a fast peak-picking algorithm best suited for
    high resolution MS data (FT-ICR-MS, Orbitrap). In high resolution data, the
    signals of ions with similar mass-to-charge ratios (m/z) exhibit little or
    no overlapping and therefore allow for a clear separation. Furthermore, ion
    signals tend to show well-defined peak shapes with narrow peak width.

    This peak-picking algorithm detects ion signals in raw data and
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
  class OPENMS_DLLAPI PeakPickerMaxima 
  {
public:

    /// Constructor
    PeakPickerMaxima(double signal_to_noise, double spacing_difference = 1.5,
        double spacing_difference_gap = 4.0, double sn_window_length = 200,
        unsigned missing = 2);

    /// Destructor
    virtual ~PeakPickerMaxima() {}

    /**
      @brief The PeakCandidate describes the output of the peak picker

      It contains the m/z and intensity value of the peak candidate.

      It also contains the original index in the m/z axis where the peak was
      found as well as an estimate of its right and left boundary. 
    */
    struct OPENMS_DLLAPI PeakCandidate 
    {
      /// index of the peak apex (relative to the input data) 
      int pos;
      /// index of the left boundary (relative to the input data) 
      int left_boundary;
      /// index of the right boundary (relative to the input data) 
      int right_boundary;
      /// m/z value of the peak apex
      double mz_max;
      /// intensity value of the peak apex
      double int_max;
    };

    /**
      @brief Will find local maxima in raw data

      @param mz_array The array containing m/z values
      @param int_array The array containing intensity values
      @param pc The resulting array containing the peak candidates
      @param check_spacings check spacing constraints? (recommended settings: yes for spectra, no for chromatograms)

      @note This function will directly report peak apices with right and left
      boundaries but will not use any fitting to estimate the true m/z and
      intensity of the peak. Note that the mz_max and int_max fields will be
      empty in the result (set to -1).

    */
    void findMaxima(const std::vector<double>& mz_array,
        const std::vector<double>& int_array,
        std::vector<PeakCandidate>& pc,
        bool check_spacings = true);

    /**
      @brief Will pick peaks in a spectrum

      @param mz_array The array containing m/z values
      @param int_array The array containing intensity values
      @param pc The resulting array containing the peak candidates
      @param check_spacings check spacing constraints? (recommended settings: yes for spectra, no for chromatograms)

      @note This function will first find maxima in the intensity domain and
      then use a spline function to estimate the best m/z and intensity for
      each peak candidate.
    */
    void pick(std::vector<double>& mz_array, 
        std::vector<double>& int_array, 
        std::vector<PeakCandidate>& pc,
        bool check_spacings = true);

protected:

    // signal-to-noise parameter
    double signal_to_noise_;

    // signal-to-noise window length 
    double sn_window_length_;

    // maximal spacing difference defining a missing data point
    double spacing_difference_;

    // maximal spacing difference defining a large gap
    double spacing_difference_gap_;

    // maximum number of missing points
    unsigned missing_;

  }; // end PeakPickerMaxima

} // namespace OpenMS

