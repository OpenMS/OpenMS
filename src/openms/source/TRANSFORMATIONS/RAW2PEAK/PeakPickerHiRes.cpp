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
// $Maintainer: Timo Sachsenberg$
// --------------------------------------------------------------------------
//

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>

#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/MATH/MISC/SplineBisection.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>


using namespace std;

namespace OpenMS
{
  PeakPickerHiRes::PeakPickerHiRes() :
    DefaultParamHandler("PeakPickerHiRes"),
    ProgressLogger()
  {
    // set default parameter values
    defaults_.setValue("signal_to_noise", 0.0, "Minimal signal-to-noise ratio for a peak to be picked (0.0 disables SNT estimation!)");
    defaults_.setMinFloat("signal_to_noise", 0.0);

    defaults_.setValue("spacing_difference_gap", 4.0, "The extension of a peak is stopped if the spacing between two subsequent data points exceeds 'spacing_difference_gap * min_spacing'. 'min_spacing' is the smaller of the two spacings from the peak apex to its two neighboring points. '0' to disable the constraint. Not applicable to chromatograms.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("spacing_difference_gap", 0.0);

    defaults_.setValue("spacing_difference", 1.5, "Maximum allowed difference between points during peak extension, in multiples of the minimal difference between the peak apex and its two neighboring points. If this difference is exceeded a missing point is assumed (see parameter 'missing'). A higher value implies a less stringent peak definition, since individual signals within the peak are allowed to be further apart. '0' to disable the constraint. Not applicable to chromatograms.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("spacing_difference", 0.0);

    defaults_.setValue("missing", 1, "Maximum number of missing points allowed when extending a peak to the left or to the right. A missing data point occurs if the spacing between two subsequent data points exceeds 'spacing_difference * min_spacing'. 'min_spacing' is the smaller of the two spacings from the peak apex to its two neighboring points. Not applicable to chromatograms.", ListUtils::create<String>("advanced"));
    defaults_.setMinInt("missing", 0);

    defaults_.setValue("ms_levels", ListUtils::create<Int>(""), "List of MS levels for which the peak picking is applied. If empty, auto mode is enabled, all peaks which aren't picked yet will get picked. Other scans are copied to the output without changes.");
    defaults_.setMinInt("ms_levels", 1);

    defaults_.setValue("report_FWHM", "false", "Add metadata for FWHM (as floatDataArray named 'FWHM' or 'FWHM_ppm', depending on param 'report_FWHM_unit') for each picked peak.");
    defaults_.setValidStrings("report_FWHM", ListUtils::create<String>("true,false"));
    defaults_.setValue("report_FWHM_unit", "relative", "Unit of FWHM. Either absolute in the unit of input, e.g. 'm/z' for spectra, or relative as ppm (only sensible for spectra, not chromatograms).");
    defaults_.setValidStrings("report_FWHM_unit", ListUtils::create<String>("relative,absolute"));

    // parameters for STN estimator
    defaults_.insert("SignalToNoise:", SignalToNoiseEstimatorMedian< MSSpectrum >().getDefaults());

    // write defaults into Param object param_
    defaultsToParam_();
    updateMembers_();
  }

  PeakPickerHiRes::~PeakPickerHiRes()
  {
  }

  void PeakPickerHiRes::pick(const MSSpectrum& input, MSSpectrum& output) const
  {
    std::vector<PeakBoundary> boundaries;
    pick(input, output, boundaries);
  }

  void PeakPickerHiRes::pick(const MSChromatogram& input, MSChromatogram& output) const
  {
    std::vector<PeakBoundary> boundaries;
    pick(input, output, boundaries);
  }

  void PeakPickerHiRes::pick(const MSSpectrum& input, MSSpectrum& output, std::vector<PeakBoundary>& boundaries, bool check_spacings) const
  {
    // copy meta data of the input spectrum
    output.clear(true);
    output.SpectrumSettings::operator=(input);
    output.MetaInfoInterface::operator=(input);
    output.setRT(input.getRT());
    output.setMSLevel(input.getMSLevel());
    output.setName(input.getName());
    output.setType(SpectrumSettings::CENTROID);
    pick_(input, output, boundaries, check_spacings);
  }

  void PeakPickerHiRes::pick(const MSChromatogram& input, MSChromatogram& output, std::vector<PeakBoundary>& boundaries, bool check_spacings) const
  {
    // copy meta data of the input chromatogram
    output.clear(true);
    output.ChromatogramSettings::operator=(input);
    output.MetaInfoInterface::operator=(input);
    output.setName(input.getName());

    pick_(input, output, boundaries, check_spacings);
  }

  template <typename ContainerType>
  void PeakPickerHiRes::pick_(const ContainerType& input, ContainerType& output, std::vector<PeakBoundary>& boundaries, bool check_spacings) const
  {
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
    SignalToNoiseEstimatorMedian< ContainerType > snt;
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
        double threshold = 1e-6;
        OpenMS::Math::spline_bisection(peak_spline, left_neighbor_mz, right_neighbor_mz, max_peak_mz, max_peak_int, threshold);

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
          }
          else
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
            } while (fabs(int_mid - fwhm_int) > threshold);
          }
          const double fwhm_left_mz = mz_mid;

          // right ...
          double mz_right = peak_raw_data.rbegin()->first;
          mz_center = max_peak_mz;
          if (peak_spline.eval(mz_right) > fwhm_int)
          { // the spline ends before half max is reached -- take the rightmost point (probably an underestimation)
            mz_mid = mz_right;
          }
          else
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

            } while (fabs(int_mid - fwhm_int) > threshold);
          }
          const double fwhm_right_mz = mz_mid;
          const double fwhm_absolute = fwhm_right_mz - fwhm_left_mz;
          output.getFloatDataArrays()[0].push_back( report_FWHM_as_ppm_ ? fwhm_absolute / max_peak_mz  * 1e6 : fwhm_absolute);
        } // FWHM

        // save picked peak into output spectrum
        typename ContainerType::PeakType peak;
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

  void PeakPickerHiRes::pickExperiment(const PeakMap& input, PeakMap& output, const bool check_spectrum_type) const
  {
    std::vector<std::vector<PeakBoundary> > boundaries_spec;
    std::vector<std::vector<PeakBoundary> > boundaries_chrom;
    pickExperiment(input, output, boundaries_spec, boundaries_chrom, check_spectrum_type);
  }

  void PeakPickerHiRes::pickExperiment(const PeakMap& input, PeakMap& output, 
                                       std::vector<std::vector<PeakBoundary> >& boundaries_spec, 
                                       std::vector<std::vector<PeakBoundary> >& boundaries_chrom,
                                       const bool check_spectrum_type) const
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
        if (ms_levels_.empty()) // auto mode
        {
          SpectrumSettings::SpectrumType spectrum_type = input[scan_idx].getType();
          if (spectrum_type == SpectrumSettings::CENTROID)
          {
            output[scan_idx] = input[scan_idx];
          }
          else
          {
            std::vector<PeakBoundary> boundaries_s; // peak boundaries of a single spectrum

            pick(input[scan_idx], output[scan_idx], boundaries_s);
            boundaries_spec.push_back(boundaries_s);
          }
        }
        else if (!ListUtils::contains(ms_levels_, input[scan_idx].getMSLevel())) // manual mode
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

  void PeakPickerHiRes::pickExperiment(/* const */ OnDiscMSExperiment& input, PeakMap& output, const bool check_spectrum_type) const
  {
    // make sure that output is clear
    output.clear(true);

    // copy experimental settings
    static_cast<ExperimentalSettings &>(output) = *input.getExperimentalSettings();

    Size progress = 0;
    startProgress(0, input.size() + input.getNrChromatograms(), "picking peaks");

    // resize output with respect to input
    output.resize(input.size());

    if (input.getNrSpectra() > 0)
    {
      for (Size scan_idx = 0; scan_idx != input.size(); ++scan_idx)
      {
        if (ms_levels_.empty()) //auto mode
        {
          MSSpectrum s = input[scan_idx];
          s.sortByPosition();

          // determine type of spectral data (profile or centroided)
          SpectrumSettings::SpectrumType spectrumType = s.getType();
          if (spectrumType == SpectrumSettings::CENTROID)
          {
            output[scan_idx] = input[scan_idx];
          }
          else
          {
            pick(s, output[scan_idx]);
          }
        }
        else if (!ListUtils::contains(ms_levels_, input[scan_idx].getMSLevel())) // manual mode
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

  void PeakPickerHiRes::updateMembers_()
  {
    signal_to_noise_ = param_.getValue("signal_to_noise");
    spacing_difference_gap_ = param_.getValue("spacing_difference_gap");
    if (spacing_difference_gap_ == 0.0) spacing_difference_gap_ = std::numeric_limits<double>::infinity();
    spacing_difference_ = param_.getValue("spacing_difference");
    if (spacing_difference_ == 0.0) spacing_difference_ = std::numeric_limits<double>::infinity();
    missing_ = param_.getValue("missing");

    ms_levels_ = getParameters().getValue("ms_levels");
    report_FWHM_ = getParameters().getValue("report_FWHM").toBool();
    report_FWHM_as_ppm_ = getParameters().getValue("report_FWHM_unit")!="absolute";
  }

}
