// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange, Alexandra Zerck $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKERCWT_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKERCWT_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePeakDeconvolution.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransform.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransformNumIntegration.h>

//#define DEBUG_PEAK_PICKING
#undef DEBUG_PEAK_PICKING
//#define DEBUG_DECONV
namespace OpenMS
{
  class PeakShape;

  /**
         @brief This class implements a peak picking algorithm using wavelet techniques

         The algorithm is described in detail in Lange et al. (2006) Proc. PSB-06.

         This peak picking algorithm uses the continuous wavelet transform of a profile data signal to detect mass peaks.
         Afterwards a given asymmetric peak function is fitted to the profile data and important peak parameters (e.g. fwhm)
         are extracted.
         In an optional step these parameters can be optimized using a non-linear optimization method.

         The peak parameters are stored in the meta data arrays of the spectra (see MSSpectrum) in this order:
         - rValue
         - area
         - fwhm
         - leftWidth
         - rightWidth
         - peakShape
         - SignalToNoise
         .

         @note The peaks must be sorted according to ascending m/z!

         @htmlinclude OpenMS_PeakPickerCWT.parameters

         @ingroup PeakPicking
  */
  class OPENMS_DLLAPI PeakPickerCWT :
    public DefaultParamHandler,
    public ProgressLogger
  {
public:
    /// Profile data iterator type
    typedef MSSpectrum::iterator PeakIterator;
    /// Const profile data iterator type
    typedef MSSpectrum::const_iterator ConstPeakIterator;

    /// Constructor
    PeakPickerCWT();

    /// Destructor
    ~PeakPickerCWT() override;

    /**
                @brief Applies the peak picking algorithm to a single spectrum.

                Picks the peaks in the input spectrum and writes the resulting peaks to the output container.
    */
    void pick(const MSSpectrum & input, MSSpectrum & output) const;

    /**
                @brief Picks the peaks in an MSExperiment.

                Picks the peaks successive in every scan in the spectrum range. The detected peaks are stored in the output MSExperiment.

        @throws Exception::UnableToFit() if peak width cannot be determined (if estimation is set to auto)
    */
    void pickExperiment(const PeakMap & input, PeakMap & output);

    /**
         @brief Estimates average peak width that can then be used for peak picking.

         The spectra with the highest TICs are used to estimate an average peak width that
         can be used as the peak_width parameter for picking the complete data set.
         Typically, the number of peaks increases with decreasing peak width until a plateau
         is reached. The beginning of this plateau is our estimate for the peak width.
         This estimate is averaged over several spectra.

    */
    double estimatePeakWidth(const PeakMap & input);

protected:

    /// Threshold for the peak height in the MS 1 level
    float peak_bound_;

    /// Threshold for the peak height in the MS 2 level
    float peak_bound_ms2_level_;

    /// Signal to noise threshold
    float signal_to_noise_;

    /// The minimal full width at half maximum
    float fwhm_bound_;

    /// The search radius for the determination of a peak's maximum position
    UInt radius_;

    /// The dilation of the wavelet
    float scale_;

    /// The threshold for correlation
    float peak_corr_bound_;

    /// The threshold for the noise level (TODO: Use the information of the signal to noise estimator)
    float noise_level_;

    /// Switch for the optimization of peak parameters
    bool optimization_;

    /// Switch for the deconvolution of peak parameters
    bool deconvolution_;

    /// Switch for the 2D optimization of peak parameters
    bool two_d_optimization_;


    void updateMembers_() override;

    /**
      @brief Class for the internal peak representation

      A regular Data-Object which contains some additional useful information
      for analyzing peaks and their properties
      The left and right iterators delimit a range in the profile data which represents a profile peak.
      They define the profile peak endpoints. @p max points to the profile data point in [left, right] with the highest intensity, the
      maximum of the profile peak.

    */
    struct OPENMS_DLLAPI PeakArea_
    {
      typedef MSSpectrum::iterator PeakIterator;
      PeakIterator left;  ///< iterator to the leftmost valid point
      PeakIterator max;   ///< iterator to the maximum position
      PeakIterator right; ///< iterator to the rightmost valid point (inclusive)
      DPosition<1> centroid_position; ///< The estimated centroid position in m/z
    };

    /// Computes the peak's left and right area
    void getPeakArea_(const PeakArea_ & area, double & area_left, double & area_right) const;

    /// Returns the best fitting peakshape
    PeakShape fitPeakShape_(const PeakArea_ & area) const;

    /**
                @brief Returns the squared Pearson coefficient.

                Computes the correlation of the peak and the original data given by the peak endpoints area.left and area.right.
                If the value is near 1, the fitted peakshape and the profile data are expected to be very similar.
    */
    double correlate_(const PeakShape & peak, const PeakArea_ & area, Int direction = 0) const;


    /**
                @brief Finds the next maximum position in the wavelet transform wt.

                If the maximum is greater than peak_bound_cwt we search for the corresponding maximum in the profile data interval [first,last)
                given a predefined search radius radius. Only peaks with intensities greater than peak_bound_
                are relevant. If no peak is detected the method return false.
                For direction=1, the method runs from first to last given direction=-1 it runs the other way around.
    */
    bool getMaxPosition_(const PeakIterator first, const PeakIterator last, const ContinuousWaveletTransform & wt, 
                         PeakArea_ & area, const Int distance_from_scan_border, 
                         const double peak_bound_cwt, const double peak_bound_ms2_level_cwt, const Int direction = 1) const;


    /**
                @brief Determines a peaks's endpoints.

                The algorithm does the following:
                - let x_m be the position of the maximum in the data and let (x_l, x_r) be
                the left and right neighbours
                -	(1) starting from x_l', walk left until one of the following happens
                - the new point is lower than the original bound => we found our left endpoint
                - the new point is larger than the last, but the point left from
                the new point is smaller. In that case, we either ran into another
                peak, or we encounter some noise. Therefore we now look in the cwt
                at the position corresponding to this value. If the cwt here is
                monotonous, we consider the point as noise and continue further to the
                left. Otherwise, we probably found the beginning of a new peak and
                therefore stop here.
                .
                -	(2) analogous procedure to the right of x_r
                .
    */
    bool getPeakEndPoints_(PeakIterator first, PeakIterator last, PeakArea_ & area, Int distance_from_scan_border,
                           Int & peak_left_index, Int & peak_right_index, ContinuousWaveletTransformNumIntegration & wt) const;


    /**
                @brief Estimates a peak's centroid position.

                Computes the centroid position of the peak using all profile data points which are greater than
                'centroid_percentage' (user-param) of the most intensive profile data point.
    */
    void getPeakCentroid_(PeakArea_ & area) const;

    /// Computes the value of a theoretical Lorentz peak at position x
    inline double lorentz_(const double height, const double lambda, const double pos, const double x) const
    {
      const double x2 = lambda * (x - pos);
      return height / (1 + x2*x2);
    }

    /**
                @brief Computes the threshold for the peak height in the wavelet transform and initializes the wavelet transform.

                Given the threshold for the peak height a corresponding value peak_bound_cwt can be computed
                for the continuous wavelet transform.
                Therefore we compute a theoretical Lorentzian peakshape with height=peak_bound_ and a width which
                is similar to the width of the wavelet. Taking the maximum in the wavelet transform of the
                Lorentzian peak we have a peak bound in the wavelet transform.
    */
    void initializeWT_(ContinuousWaveletTransformNumIntegration& wt, const double peak_bound_in, double& peak_bound_ms_cwt) const;

    /** @name Methods needed for separation of overlapping peaks
     */
    //@{

    /**
            @brief Separates overlapping peaks.

            It determines the number of peaks lying underneath the initial peak using the cwt with different scales.
            Then a nonlinear optimization procedure is applied to optimize the peak parameters.
    */
    bool deconvolutePeak_(PeakShape & shape, std::vector<PeakShape> & peak_shapes, double peak_bound_cwt) const;

    /// Determines the number of peaks in the given mass range using the cwt
    Int getNumberOfPeaks_(ConstPeakIterator first, ConstPeakIterator last, std::vector<double> & peak_values,
                          Int direction, double resolution, ContinuousWaveletTransformNumIntegration & wt, double peak_bound_cwt) const;

    /// Estimate the charge state of the peaks
    Int determineChargeState_(std::vector<double> & peak_values) const;

    /// Add a peak
    void addPeak_(std::vector<PeakShape> & peaks_DC, PeakArea_ & area, double left_width, double right_width, OptimizePeakDeconvolution::Data & data) const;
    //@}
  };  // end PeakPickerCWT

} // namespace OpenMS

#endif
