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
// $Author: Erhan Kenar $
// $Maintainer: Timo Sachsenberg $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKERHIRES_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKERHIRES_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>

#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>


#define DEBUG_PEAK_PICKING
#undef DEBUG_PEAK_PICKING
//#undef DEBUG_DECONV
namespace OpenMS
{
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
    void pick(const MSSpectrum& input, MSSpectrum& output) const
    {
      std::vector<PeakBoundary> boundaries;
      pick(input, output, boundaries);
    }

     /**
     * @brief Applies the peak-picking algorithm to a single chromatogram
     * (MSChromatogram). The resulting picked peaks are written to the output chromatogram.
     *
     * @param input  input chromatogram in profile mode
     * @param output  output chromatogram with picked peaks
     */
    void pick(const MSChromatogram& input, MSChromatogram& output) const
    {
      std::vector<PeakBoundary> boundaries;
      pick(input, output, boundaries);
    }

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
    void pick(const MSSpectrum& input, MSSpectrum& output, std::vector<PeakBoundary>& boundaries, bool check_spacings = true) const
    {
      // copy meta data of the input spectrum
      output.clear(true);
      output.SpectrumSettings::operator=(input);
      output.MetaInfoInterface::operator=(input);
      output.setRT(input.getRT());
      output.setMSLevel(input.getMSLevel());
      output.setName(input.getName());
      output.setType(SpectrumSettings::CENTROID);
      if (report_FWHM_)
      {
        output.getFloatDataArrays().resize(1);
        output.getFloatDataArrays()[0].setName( report_FWHM_as_ppm_ ? "FWHM_ppm" : "FWHM");
      }
      
      // don't pick a spectrum with less than 5 data points
      if (input.size() < 5) return;

      // if both spacing constraints are disabled, don't check spacings at all:
      if ((spacing_difference_ == std::numeric_limits<double>::infinity()) &&
          (spacing_difference_gap_ == std::numeric_limits<double>::infinity()))
      {
        check_spacings = false;
      }

      // signal-to-noise estimation
      SignalToNoiseEstimatorMedian<MSSpectrum > snt;
      snt.setParameters(param_.copy("SignalToNoise:", true));

      if (signal_to_noise_ > 0.0)
      {
        snt.init(input);
      }

      // find local maxima in profile data
      for (Size i = 2; i < input.size() - 2; ++i)
      {
        double central_peak_mz = input[i].getMZ(), central_peak_int = input[i].getIntensity();
        double left_neighbor_mz = input[i - 1].getMZ(), left_neighbor_int = input[i - 1].getIntensity();
        double right_neighbor_mz = input[i + 1].getMZ(), right_neighbor_int = input[i + 1].getIntensity();

        // do not interpolate when the left or right support is a zero-data-point
        if (std::fabs(left_neighbor_int) < std::numeric_limits<double>::epsilon()) continue;
        if (std::fabs(right_neighbor_int) < std::numeric_limits<double>::epsilon()) continue;

        // MZ spacing sanity checks
        double left_to_central = 0.0, central_to_right = 0.0, min_spacing = 0.0;
        if (check_spacings)
        {
          left_to_central = central_peak_mz - left_neighbor_mz;
          central_to_right = right_neighbor_mz - central_peak_mz;
          min_spacing = (left_to_central < central_to_right) ? left_to_central : central_to_right;
        }

        double act_snt = 0.0, act_snt_l1 = 0.0, act_snt_r1 = 0.0;
        if (signal_to_noise_ > 0.0)
        {
          act_snt = snt.getSignalToNoise(input[i]);
          act_snt_l1 = snt.getSignalToNoise(input[i - 1]);
          act_snt_r1 = snt.getSignalToNoise(input[i + 1]);
        }

        // look for peak cores meeting MZ and intensity/SNT criteria
        if ((central_peak_int > left_neighbor_int) && 
            (central_peak_int > right_neighbor_int) && 
            (act_snt >= signal_to_noise_) && 
            (act_snt_l1 >= signal_to_noise_) && 
            (act_snt_r1 >= signal_to_noise_) &&
            (!check_spacings || 
             ((left_to_central < spacing_difference_ * min_spacing) && 
              (central_to_right < spacing_difference_ * min_spacing))))
        {
          // special case: if a peak core is surrounded by more intense
          // satellite peaks (indicates oscillation rather than
          // real peaks) -> remove

          double act_snt_l2 = 0.0, act_snt_r2 = 0.0;

          if (signal_to_noise_ > 0.0)
          {
            act_snt_l2 = snt.getSignalToNoise(input[i - 2]);
            act_snt_r2 = snt.getSignalToNoise(input[i + 2]);
          }

          // checking signal-to-noise?
          if ((i > 1) &&
              (i + 2 < input.size()) &&
              (left_neighbor_int < input[i - 2].getIntensity()) &&
              (right_neighbor_int < input[i + 2].getIntensity()) &&
              (act_snt_l2 >= signal_to_noise_) &&
              (act_snt_r2 >= signal_to_noise_) &&
              (!check_spacings ||
               ((left_neighbor_mz - input[i - 2].getMZ() < spacing_difference_ * min_spacing) && 
                (input[i + 2].getMZ() - right_neighbor_mz < spacing_difference_ * min_spacing))))
          {
            ++i;
            continue;
          }

          std::map<double, double> peak_raw_data;

          peak_raw_data[central_peak_mz] = central_peak_int;
          peak_raw_data[left_neighbor_mz] = left_neighbor_int;
          peak_raw_data[right_neighbor_mz] = right_neighbor_int;

          // peak core found, now extend it
          // to the left
          Size k = 2;

          bool previous_zero_left(false); // no need to extend peak if previous intensity was zero
          Size missing_left(0);
          Size left_boundary(i - 1); // index of the left boundary for the spline interpolation
          
          while ((k <= i) && // prevent underflow
                 (i - k + 1 > 0) && 
                 !previous_zero_left && 
                 (missing_left <= missing_) && 
                 (input[i - k].getIntensity() <= peak_raw_data.begin()->second) &&
                 (!check_spacings || 
                  (peak_raw_data.begin()->first - input[i - k].getMZ() < spacing_difference_gap_ * min_spacing)))
          {
            double act_snt_lk = 0.0;

            if (signal_to_noise_ > 0.0)
            {
              act_snt_lk = snt.getSignalToNoise(input[i - k]);
            }

            if ((act_snt_lk >= signal_to_noise_) && 
                (!check_spacings ||
                 (peak_raw_data.begin()->first - input[i - k].getMZ() < spacing_difference_ * min_spacing)))
            {
              peak_raw_data[input[i - k].getMZ()] = input[i - k].getIntensity();
            }
            else
            {
              ++missing_left;
              if (missing_left <= missing_)
              {
                peak_raw_data[input[i - k].getMZ()] = input[i - k].getIntensity();
              }
            }

            previous_zero_left = (input[i - k].getIntensity() == 0);
            left_boundary = i - k;
            ++k;
          }

          // to the right
          k = 2;

          bool previous_zero_right(false); // no need to extend peak if previous intensity was zero
          Size missing_right(0);
          Size right_boundary(i+1); // index of the right boundary for the spline interpolation

          while ((i + k < input.size()) && 
                 !previous_zero_right && 
                 (missing_right <= missing_) && 
                 (input[i + k].getIntensity() <= peak_raw_data.rbegin()->second) &&
                 (!check_spacings ||
                  (input[i + k].getMZ() - peak_raw_data.rbegin()->first < spacing_difference_gap_ * min_spacing)))
          {
            double act_snt_rk = 0.0;

            if (signal_to_noise_ > 0.0)
            {
              act_snt_rk = snt.getSignalToNoise(input[i + k]);
            }

            if ((act_snt_rk >= signal_to_noise_) && 
                (!check_spacings ||
                 (input[i + k].getMZ() - peak_raw_data.rbegin()->first < spacing_difference_ * min_spacing)))
            {
              peak_raw_data[input[i + k].getMZ()] = input[i + k].getIntensity();
            }
            else
            {
              ++missing_right;
              if (missing_right <= missing_)
              {
                peak_raw_data[input[i + k].getMZ()] = input[i + k].getIntensity();
              }
            }

            previous_zero_right = (input[i + k].getIntensity() == 0);
            right_boundary = i + k;
            ++k;
          }

          // skip if the minimal number of 3 points for fitting is not reached
          if (peak_raw_data.size() < 3) continue;

          CubicSpline2d peak_spline (peak_raw_data);

          // calculate maximum by evaluating the spline's 1st derivative
          // (bisection method)
          double max_peak_mz = central_peak_mz;
          double max_peak_int = central_peak_int;
          double threshold = 0.000001;
          double lefthand = left_neighbor_mz;
          double righthand = right_neighbor_mz;

          bool lefthand_sign = 1;
          double eps = std::numeric_limits<double>::epsilon();

          // bisection
          do
          {
            double mid = (lefthand + righthand) / 2.0;
            double midpoint_deriv_val = peak_spline.derivatives(mid, 1);

            // if deriv nearly zero then maximum already found
            if (!(std::fabs(midpoint_deriv_val) > eps))
            {
              break;
            }

            bool midpoint_sign = (midpoint_deriv_val < 0.0) ? 0 : 1;

            if (lefthand_sign ^ midpoint_sign)
            {
              righthand = mid;
            }
            else
            {
              lefthand = mid;
            }
          }
          while (righthand - lefthand > threshold);

          max_peak_mz = (lefthand + righthand) / 2;
          max_peak_int = peak_spline.eval(max_peak_mz);

          //
          // compute FWHM
          //
          if (report_FWHM_)
          {
            double fwhm_int = max_peak_int / 2.0;
            threshold = 0.01 * fwhm_int;
            double mz_mid, int_mid; 
            // left:
            double mz_left = peak_raw_data.begin()->first;
            double mz_center = max_peak_mz;
            if (peak_spline.eval(mz_left) > fwhm_int)
            { // the spline ends before half max is reached -- take the leftmost point (probably an underestimation)
              mz_mid = mz_left;
            } else
            {
              do 
              {
                mz_mid = mz_left / 2 + mz_center / 2;
                int_mid = peak_spline.eval(mz_mid);
                if (int_mid < fwhm_int)
                {
                  mz_left = mz_mid;
                }
                else
                {
                  mz_center = mz_mid;
                }
              } while(fabs(int_mid - fwhm_int) > threshold);
            }
            const double fwhm_left_mz = mz_mid;

            // right ...
            double mz_right = peak_raw_data.rbegin()->first;
            mz_center = max_peak_mz;
            if (peak_spline.eval(mz_right) > fwhm_int)
            { // the spline ends before half max is reached -- take the rightmost point (probably an underestimation)
              mz_mid = mz_right;
            } else
              {
              do 
              {
                mz_mid = mz_right / 2 + mz_center / 2;
                int_mid = peak_spline.eval(mz_mid);
                if (int_mid < fwhm_int)
                {
                  mz_right = mz_mid;
                }
                else
                {
                  mz_center = mz_mid;
                }

              } while(fabs(int_mid - fwhm_int) > threshold);
            }
            const double fwhm_right_mz = mz_mid;
            const double fwhm_absolute = fwhm_right_mz - fwhm_left_mz;
            output.getFloatDataArrays()[0].push_back( report_FWHM_as_ppm_ ? fwhm_absolute / max_peak_mz  * 1e6 : fwhm_absolute);
          } // FWHM

          // save picked peak into output spectrum
          Peak1D peak;
          PeakBoundary peak_boundary;
          peak.setMZ(max_peak_mz);
          peak.setIntensity(max_peak_int);
          peak_boundary.mz_min = input[left_boundary].getMZ();
          peak_boundary.mz_max = input[right_boundary].getMZ();
          output.push_back(peak);
          
          boundaries.push_back(peak_boundary);

          // jump over profile data points that have been considered already
          i = i + k - 1;
        }
      }

      return;
    }


    /**
     * @brief Applies the peak-picking algorithm to a single chromatogram
     * (MSChromatogram). The resulting picked peaks are written to the output chromatogram.
     *
     * @param input  input chromatogram in profile mode
     * @param output  output chromatogram with picked peaks
     * @param boundaries  boundaries of the picked peaks
     */
    void pick(const MSChromatogram& input, MSChromatogram& output, std::vector<PeakBoundary>& boundaries) const
    {
      // copy meta data of the input chromatogram
      output.clear(true);
      output.ChromatogramSettings::operator=(input);
      output.MetaInfoInterface::operator=(input);
      output.setName(input.getName());

      MSSpectrum input_spectrum;
      MSSpectrum output_spectrum;
      for (MSChromatogram::const_iterator it = input.begin(); it != input.end(); ++it)
      {
        Peak1D p;
        p.setMZ(it->getRT());
        p.setIntensity(it->getIntensity());
        input_spectrum.push_back(p);
      }

      pick(input_spectrum, output_spectrum, boundaries, false); // no spacing checks!

      for (MSSpectrum::const_iterator it = output_spectrum.begin(); it != output_spectrum.end(); ++it)
      {
        ChromatogramPeak p;
        p.setRT(it->getMZ());
        p.setIntensity(it->getIntensity());
        output.push_back(p);
      }

      // copy float data arrays (for FWHM)
      output.getFloatDataArrays().resize(output_spectrum.getFloatDataArrays().size());
      for (Size i = 0; i < output_spectrum.getFloatDataArrays().size(); ++i)
      {
        output.getFloatDataArrays()[i].insert(output.getFloatDataArrays()[i].begin(), output_spectrum.getFloatDataArrays()[i].begin(), output_spectrum.getFloatDataArrays()[i].end());
        output.getFloatDataArrays()[i].setName(output_spectrum.getFloatDataArrays()[i].getName());
      }
    }

    /**
     * @brief Applies the peak-picking algorithm to a map (MSExperiment). This
     * method picks peaks for each scan in the map consecutively. The resulting
     * picked peaks are written to the output map.
     *
     * @param input  input map in profile mode
     * @param output  output map with picked peaks
     * @param check_spectrum_type  if set, checks spectrum type and throws an exception if a centroided spectrum is passed 
     */
    void pickExperiment(const PeakMap& input, PeakMap& output, const bool check_spectrum_type = true) const
    {
        std::vector<std::vector<PeakBoundary> > boundaries_spec;
        std::vector<std::vector<PeakBoundary> > boundaries_chrom;
        pickExperiment(input, output, boundaries_spec, boundaries_chrom, check_spectrum_type);
    }

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
    void pickExperiment(const PeakMap& input, PeakMap& output, std::vector<std::vector<PeakBoundary> >& boundaries_spec, std::vector<std::vector<PeakBoundary> >& boundaries_chrom, const bool check_spectrum_type = true) const
    {
      // make sure that output is clear
      output.clear(true);

      // copy experimental settings
      static_cast<ExperimentalSettings &>(output) = input;

      // resize output with respect to input
      output.resize(input.size());

      Size progress = 0;
      startProgress(0, input.size() + input.getChromatograms().size(), "picking peaks");

      if (input.getNrSpectra() > 0)
      {
        for (Size scan_idx = 0; scan_idx != input.size(); ++scan_idx)
        {
          if (!ListUtils::contains(ms_levels_, input[scan_idx].getMSLevel()))
          {
            output[scan_idx] = input[scan_idx];
          }
          else
          {
            std::vector<PeakBoundary> boundaries_s; // peak boundaries of a single spectrum

            // determine type of spectral data (profile or centroided)
            SpectrumSettings::SpectrumType spectrum_type = input[scan_idx].getType();

            if (spectrum_type == SpectrumSettings::CENTROID && check_spectrum_type)
            {
              throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, __FUNCTION__, "Error: Centroided data provided but profile spectra expected.");
            }

            pick(input[scan_idx], output[scan_idx], boundaries_s);
            boundaries_spec.push_back(boundaries_s);
          }
          setProgress(++progress);
        }
      }


      for (Size i = 0; i < input.getChromatograms().size(); ++i)
      {
        MSChromatogram chromatogram;
        std::vector<PeakBoundary> boundaries_c; // peak boundaries of a single chromatogram
        pick(input.getChromatograms()[i], chromatogram, boundaries_c);
        output.addChromatogram(chromatogram);
        boundaries_chrom.push_back(boundaries_c);
        setProgress(++progress);
      }
      endProgress();

      return;
    }

    /**
      @brief Applies the peak-picking algorithm to a map (MSExperiment). This
      method picks peaks for each scan in the map consecutively. The resulting
      picked peaks are written to the output map.

      Currently we have to give up const-correctness but we know that everything on disc is constant
    */
    void pickExperiment(/* const */ OnDiscPeakMap& input, PeakMap& output, const bool check_spectrum_type = true) const
    {
      // make sure that output is clear
      output.clear(true);

      // copy experimental settings
      static_cast<ExperimentalSettings &>(output) = *input.getExperimentalSettings();

      Size progress = 0;
      startProgress(0, input.size() + input.getNrChromatograms(), "picking peaks");

      if (input.getNrSpectra() > 0)
      {

        // resize output with respect to input
        output.resize(input.size());

        for (Size scan_idx = 0; scan_idx != input.size(); ++scan_idx)
        {
          if (!ListUtils::contains(ms_levels_, input[scan_idx].getMSLevel()))
          {
            output[scan_idx] = input[scan_idx];
          }
          else
          {
            MSSpectrum s = input[scan_idx];
            s.sortByPosition();

            // determine type of spectral data (profile or centroided)
            SpectrumSettings::SpectrumType spectrum_type = s.getType();

            if (spectrum_type == SpectrumSettings::CENTROID && check_spectrum_type)
            {
              throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, __FUNCTION__, "Error: Centroided data provided but profile spectra expected.");
            }

            pick(s, output[scan_idx]);
          }
          setProgress(++progress);
        }
      }

      for (Size i = 0; i < input.getNrChromatograms(); ++i)
      {
        MSChromatogram chromatogram;
        pick(input.getChromatogram(i), chromatogram);
        output.addChromatogram(chromatogram);
        setProgress(++progress);
      }
      endProgress();

      return;
    }

protected:
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

#endif
