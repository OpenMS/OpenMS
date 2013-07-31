// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Erhan Kenar $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKERHIRES_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKERHIRES_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

#include <map>


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
  class OPENMS_DLLAPI PeakPickerHiRes :
    public DefaultParamHandler,
    public ProgressLogger
  {
public:
    /// Constructor
    PeakPickerHiRes();

    /// Destructor
    virtual ~PeakPickerHiRes();

    /**
      @brief Applies the peak-picking algorithm to a single spectrum
      (MSSpectrum). The resulting picked peaks are written to the output
      spectrum.
    */
    template <typename PeakType>
    void pick(const MSSpectrum<PeakType> & input, MSSpectrum<PeakType> & output) const
    {
      // copy meta data of the input spectrum
      output.clear(true);
      output.SpectrumSettings::operator=(input);
      output.MetaInfoInterface::operator=(input);
      output.setRT(input.getRT());
      output.setMSLevel(input.getMSLevel());
      output.setName(input.getName());
      output.setType(SpectrumSettings::PEAKS);

      // don't pick a spectrum with less than 5 data points
      if (input.size() < 5) return;

      // signal-to-noise estimation
      SignalToNoiseEstimatorMedian<MSSpectrum<PeakType> > snt;

      if (signal_to_noise_ > 0.0)
      {
        snt.init(input);
      }

      // find local maxima in raw data
      for (Size i = 2; i < input.size() - 2; ++i)
      {
        double central_peak_mz = input[i].getMZ(), central_peak_int = input[i].getIntensity();
        double left_neighbor_mz = input[i - 1].getMZ(), left_neighbor_int = input[i - 1].getIntensity();
        double right_neighbor_mz = input[i + 1].getMZ(), right_neighbor_int = input[i + 1].getIntensity();

        // MZ spacing sanity checks
        double left_to_central = std::fabs(central_peak_mz - left_neighbor_mz);
        double central_to_right = std::fabs(right_neighbor_mz - central_peak_mz);
        double min_spacing = (left_to_central < central_to_right) ? left_to_central : central_to_right;

        double act_snt = 0.0, act_snt_l1 = 0.0, act_snt_r1 = 0.0;

        if (signal_to_noise_ > 0.0)
        {
          act_snt = snt.getSignalToNoise(input[i]);
          act_snt_l1 = snt.getSignalToNoise(input[i - 1]);
          act_snt_r1 = snt.getSignalToNoise(input[i + 1]);
        }

        // look for peak cores meeting MZ and intensity/SNT criteria
        if (act_snt >= signal_to_noise_
           && left_to_central < 1.5 * min_spacing
           && central_peak_int > left_neighbor_int
           && act_snt_l1 >= signal_to_noise_
           && central_to_right < 1.5 * min_spacing
           && central_peak_int > right_neighbor_int
           && act_snt_r1 >= signal_to_noise_)
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

          if ((i > 1
              && std::fabs(left_neighbor_mz - input[i - 2].getMZ()) < 1.5 * min_spacing
              && left_neighbor_int < input[i - 2].getIntensity()
              && act_snt_l2 >= signal_to_noise_)
             &&
              ((i + 2) < input.size()
              && std::fabs(input[i + 2].getMZ() - right_neighbor_mz) < 1.5 * min_spacing
              && right_neighbor_int < input[i + 2].getIntensity()
              && act_snt_r2 >= signal_to_noise_)
              )
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

          Size missing_left(0);
          Size missing_right(0);

          while ((i - k + 1) > 0
                && (missing_left < 2)
                && input[i - k].getIntensity() <= peak_raw_data.begin()->second)
          {

            double act_snt_lk = 0.0;

            if (signal_to_noise_ > 0.0)
            {
              act_snt_lk = snt.getSignalToNoise(input[i - k]);
            }


            if (act_snt_lk >= signal_to_noise_ && std::fabs(input[i - k].getMZ() - peak_raw_data.begin()->first) < 1.5 * min_spacing)
            {
              peak_raw_data[input[i - k].getMZ()] = input[i - k].getIntensity();
            }
            else
            {
              peak_raw_data[input[i - k].getMZ()] = input[i - k].getIntensity();
              ++missing_left;
            }

            ++k;

          }

          // to the right
          k = 2;
          while ((i + k) < input.size()
                && (missing_right < 2)
                && input[i + k].getIntensity() <= peak_raw_data.rbegin()->second)
          {

            double act_snt_rk = 0.0;

            if (signal_to_noise_ > 0.0)
            {
              act_snt_rk = snt.getSignalToNoise(input[i + k]);
            }

            if (act_snt_rk >= signal_to_noise_ && std::fabs(input[i + k].getMZ() - peak_raw_data.rbegin()->first) < 1.5 * min_spacing)
            {
              peak_raw_data[input[i + k].getMZ()] = input[i + k].getIntensity();
            }
            else
            {
              peak_raw_data[input[i + k].getMZ()] = input[i + k].getIntensity();
              ++missing_right;
            }

            ++k;
          }


          // output all raw data points selected for one peak
          // TODO: #ifdef DEBUG_ ...
          // for (std::map<double, double>::const_iterator map_it = peak_raw_data.begin(); map_it != peak_raw_data.end(); ++map_it) {
          // PeakType peak;
          // peak.setMZ(map_it->first);
          // peak.setIntensity(map_it->second);
          // output.push_back(peak);
          // std::cout << map_it->first << " " << map_it->second << " snt: " << std::endl;
          // }
          // std::cout << "--------------------" << std::endl;

          const Size num_raw_points = peak_raw_data.size();

          std::vector<double> raw_mz_values;
          std::vector<double> raw_int_values;

          for (std::map<double, double>::const_iterator map_it = peak_raw_data.begin(); map_it != peak_raw_data.end(); ++map_it)
          {
            raw_mz_values.push_back(map_it->first);
            raw_int_values.push_back(map_it->second);
          }

          // setup gsl splines
          gsl_interp_accel * spline_acc = gsl_interp_accel_alloc();
          gsl_interp_accel * first_deriv_acc = gsl_interp_accel_alloc();
          gsl_spline * peak_spline = gsl_spline_alloc(gsl_interp_cspline, num_raw_points);
          gsl_spline_init(peak_spline, &(*raw_mz_values.begin()), &(*raw_int_values.begin()), num_raw_points);


          // calculate maximum by evaluating the spline's 1st derivative
          // (bisection method)
          double max_peak_mz = central_peak_mz, max_peak_int = central_peak_int;
          double threshold = 0.000001;
          double lefthand = left_neighbor_mz;
          double righthand = right_neighbor_mz;

          bool lefthand_sign = 1;
          double eps = std::numeric_limits<double>::epsilon();


          // bisection
          do
          {
            double mid = (lefthand + righthand) / 2;

            double midpoint_deriv_val = gsl_spline_eval_deriv(peak_spline, mid, first_deriv_acc);

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

            // TODO: #ifdef DEBUG_ ...
            // PeakType peak;
            // peak.setMZ(mid);
            // peak.setIntensity(gsl_spline_eval(peak_spline, mid, spline_acc));
            // output.push_back(peak);

          }
          while (std::fabs(lefthand - righthand) > threshold);

          // sanity check?
          max_peak_mz = (lefthand + righthand) / 2;
          max_peak_int = gsl_spline_eval(peak_spline, max_peak_mz, spline_acc);

          // save picked pick into output spectrum
          PeakType peak;
          peak.setMZ(max_peak_mz);
          peak.setIntensity(max_peak_int);
          output.push_back(peak);

          // free allocated gsl memory
          gsl_spline_free(peak_spline);
          gsl_interp_accel_free(spline_acc);
          gsl_interp_accel_free(first_deriv_acc);

          // jump over raw data points that have been considered already
          i = i + k - 1;
        }
      }

      return;
    }

    /**
      @brief Applies the peak-picking algorithm to a single chromatogram
      (MSChromatogram). The resulting picked peaks are written to the output
      chromatogram.
    */
    template <typename PeakType>
    void pick(const MSChromatogram<PeakType> & input, MSChromatogram<PeakType> & output) const
    {
      // copy meta data of the input chromatogram
      output.clear(true);
      output.ChromatogramSettings::operator=(input);
      output.MetaInfoInterface::operator=(input);
      output.setName(input.getName());

      MSSpectrum<PeakType> input_spectrum;
      MSSpectrum<PeakType> output_spectrum;
      for (typename MSChromatogram<PeakType>::const_iterator it = input.begin(); it != input.end(); ++it)
      {
        input_spectrum.push_back(*it);
      }
      pick(input_spectrum, output_spectrum);
      for (typename MSSpectrum<PeakType>::const_iterator it = output_spectrum.begin(); it != output_spectrum.end(); ++it)
      {
        output.push_back(*it);
      }

    }

    /**
      @brief Applies the peak-picking algorithm to a map (MSExperiment). This
      method picks peaks for each scan in the map consecutively. The resulting
      picked peaks are written to the output map.
    */
    template <typename PeakType, typename ChromatogramPeakT>
    void pickExperiment(const MSExperiment<PeakType, ChromatogramPeakT> & input, MSExperiment<PeakType, ChromatogramPeakT> & output) const
    {
      // make sure that output is clear
      output.clear(true);

      // copy experimental settings
      static_cast<ExperimentalSettings &>(output) = input;

      // resize output with respect to input
      output.resize(input.size());

      bool ms1_only = param_.getValue("ms1_only").toBool();
      Size progress = 0;

      startProgress(0, input.size() + input.getChromatograms().size(), "smoothing data");
      for (Size scan_idx = 0; scan_idx != input.size(); ++scan_idx)
      {
        if (ms1_only && (input[scan_idx].getMSLevel() != 1))
        {
          output[scan_idx] = input[scan_idx];
        }
        else
        {
          pick(input[scan_idx], output[scan_idx]);
        }
        setProgress(++progress);
      }
      for (Size i = 0; i < input.getChromatograms().size(); ++i)
      {
        MSChromatogram<ChromatogramPeakT> chromatogram;
        pick(input.getChromatograms()[i], chromatogram);
        output.addChromatogram(chromatogram);
        setProgress(++progress);
      }

      endProgress();

      return;
    }

protected:
    // signal-to-noise parameter
    double signal_to_noise_;

    // docu in base class
    void updateMembers_();

  }; // end PeakPickerHiRes

} // namespace OpenMS

#endif
