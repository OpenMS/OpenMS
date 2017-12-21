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

#include <cmath>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>

#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/TwoDOptimization.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePick.h>
#include <OpenMS/FILTERING/TRANSFORMERS/TICFilter.h>
#include <boost/math/special_functions/fpclassify.hpp>

#ifdef _OPENMP
#ifdef OPENMS_WINDOWSPLATFORM
#include <omp.h>
#endif
#endif

#include <cmath>
#include <vector>
#include <algorithm>

//#define DEBUG_PEAK_PICKING2 1

using namespace std;

namespace OpenMS
{
  PeakPickerCWT::PeakPickerCWT() :
    DefaultParamHandler("PeakPickerCWT"),
    ProgressLogger(),
    radius_(0),
    scale_(0.0),
    peak_corr_bound_(0.0),
    noise_level_(0.0),
    optimization_(false)
  {
    defaults_.setValue("signal_to_noise", 1.0, "Minimal signal to noise ratio for a peak to be picked.");
    defaults_.setMinFloat("signal_to_noise", 0.0);
    defaults_.setValue("thresholds:peak_bound", 10.0, "Minimal peak intensity.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("thresholds:peak_bound", 0.0);
    defaults_.setValue("thresholds:peak_bound_ms2_level", 10.0, "Minimal peak intensity for MS/MS peaks.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("thresholds:peak_bound_ms2_level", 0.0);

    // if a peak picking parameter is missed in the param object the value should be substituted by a default value
    defaults_.setValue("centroid_percentage", 0.8, "Percentage of the maximum height that the raw data points must exceed to be taken into account for the calculation of the centroid. " \
                                                   "If it is 1 the centroid position corresponds to the position of the highest intensity.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("centroid_percentage", 0.0);
    defaults_.setMaxFloat("centroid_percentage", 1.0);
    defaults_.setValue("thresholds:correlation", 0.5, "minimal correlation of a peak and the raw signal. " \
                                                      "If a peak has a lower correlation it is skipped.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("thresholds:correlation", 0.0);
    defaults_.setMaxFloat("thresholds:correlation", 1.0);
    defaults_.setValue("peak_width", 0.15, "Approximate fwhm of the peaks.");
    defaults_.setMinFloat("peak_width", 0.0);
    defaults_.setValue("estimate_peak_width", "false", "Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.");
    std::vector<String> valid_opts;
    valid_opts.push_back("true");
    valid_opts.push_back("false");
    defaults_.setValidStrings("estimate_peak_width", valid_opts);

    defaults_.setValue("fwhm_lower_bound_factor", 0.7, "Factor that calculates the minimal fwhm value from the peak_width. All peaks with width smaller than fwhm_bound_factor * peak_width are discarded.", ListUtils::create<String>("advanced"));
    defaults_.setValue("fwhm_upper_bound_factor", 20., "Factor that calculates the maximal fwhm value from the peak_width. All peaks with width greater than fwhm_upper_bound_factor * peak_width are discarded.", ListUtils::create<String>("advanced"));

    defaults_.setMinFloat("fwhm_lower_bound_factor", 0.0);
    defaults_.setMinFloat("fwhm_upper_bound_factor", 0.0);

    defaults_.setValue("wavelet_transform:spacing", 0.001, "Spacing of the CWT. Note that the accuracy of the picked peak's centroid position depends in the Raw data spacing, i.e., 50% of raw peak distance at most.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("wavelet_transform:spacing", 0.0);
    defaults_.setValue("thresholds:noise_level", 0.1, "noise level for the search of the peak endpoints.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("thresholds:noise_level", 0.0);
    defaults_.setValue("thresholds:search_radius", 3, "search radius for the search of the maximum in the signal after a maximum in the cwt was found", ListUtils::create<String>("advanced"));
    defaults_.setMinInt("thresholds:search_radius", 0);

    //Optimization parameters
    defaults_.setValue("optimization", "no", "If the peak parameters position, intensity and left/right width" \
                                             "shall be optimized set optimization to one_dimensional or two_dimensional.", ListUtils::create<String>("advanced"));
    valid_opts.clear();
    valid_opts.push_back("no");
    valid_opts.push_back("one_dimensional");
    valid_opts.push_back("two_dimensional");
    defaults_.setValidStrings("optimization", valid_opts);
    defaults_.setValue("optimization:penalties:position", 0.0, "penalty term for the fitting of the position:" \
                                                               "If it differs too much from the initial one it can be penalized ", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("optimization:penalties:position", 0.0);
    defaults_.setValue("optimization:penalties:left_width", 1.0, "penalty term for the fitting of the left width:" \
                                                                 "If the left width differs too much from the initial one during the fitting it can be penalized.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("optimization:penalties:left_width", 0.0);
    defaults_.setValue("optimization:penalties:right_width", 1.0, "penalty term for the fitting of the right width:" \
                                                                  "If the right width differs too much from the initial one during the fitting it can be penalized.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("optimization:penalties:right_width", 0.0);
    defaults_.setValue("optimization:penalties:height", 1.0, "penalty term for the fitting of the intensity (only used in 2D Optimization):" \
                                                             "If it gets negative during the fitting it can be penalized.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("optimization:penalties:height", 0.0);
    defaults_.setValue("optimization:iterations", 400, "maximal number of iterations for the fitting step", ListUtils::create<String>("advanced"));
    defaults_.setMinInt("optimization:iterations", 1);
    // additional 2d optimization parameters
    defaults_.setValue("optimization:2d:tolerance_mz", 2.2, "mz tolerance for cluster construction", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("optimization:2d:tolerance_mz", 0.0);
    defaults_.setValue("optimization:2d:max_peak_distance", 1.2, "maximal peak distance in mz in a cluster", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("optimization:2d:max_peak_distance", 0.0);
    // deconvolution parameters
    defaults_.setValue("deconvolution:deconvolution", "false", "If you want heavily overlapping peaks to be separated set this value to \"true\"", ListUtils::create<String>("advanced"));
    valid_opts.clear();
    valid_opts.push_back("true");
    valid_opts.push_back("false");
    defaults_.setValidStrings("deconvolution:deconvolution", valid_opts);
    defaults_.setValue("deconvolution:asym_threshold", 0.3, "If the symmetry of a peak is smaller than asym_thresholds it is assumed that it consists of more than one peak and the deconvolution procedure is started.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("deconvolution:asym_threshold", 0.0);
    defaults_.setValue("deconvolution:left_width", 2.0, "1/left_width is the initial value for the left width of the peaks found in the deconvolution step.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("deconvolution:left_width", 0.0);
    defaults_.setValue("deconvolution:right_width", 2.0, "1/right_width is the initial value for the right width of the peaks found in the deconvolution step.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("deconvolution:right_width", 0.0);
    defaults_.setValue("deconvolution:scaling", 0.12, "Initial scaling of the cwt used in the separation of heavily overlapping peaks. The initial value is used for charge 1, for higher charges it is adapted to scaling/charge.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("deconvolution:scaling", 0.0);
    defaults_.setValue("deconvolution:fitting:penalties:position", 0.0, "penalty term for the fitting of the peak position:" \
                                                                        "If the position changes more than 0.5Da during the fitting it can be penalized as well as " \
                                                                        "discrepancies of the peptide mass rule.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("deconvolution:fitting:penalties:position", 0.0);
    defaults_.setValue("deconvolution:fitting:penalties:height", 1.0, "penalty term for the fitting of the intensity:" \
                                                                      "If it gets negative during the fitting it can be penalized.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("deconvolution:fitting:penalties:height", 0.0);
    defaults_.setValue("deconvolution:fitting:penalties:left_width", 0.0, "penalty term for the fitting of the left width:" \
                                                                          "If the left width gets too broad or negative during the fitting it can be penalized.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("deconvolution:fitting:penalties:left_width", 0.0);
    defaults_.setValue("deconvolution:fitting:penalties:right_width", 0.0, "penalty term for the fitting of the right width:" \
                                                                           "If the right width gets too broad or negative during the fitting it can be penalized.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("deconvolution:fitting:penalties:right_width", 0.0);
    defaults_.setValue("deconvolution:fitting:fwhm_threshold", 0.7, "If the FWHM of a peak is higher than 'fwhm_thresholds' it is assumed that it consists of more than one peak and the deconvolution procedure is started.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("deconvolution:fitting:fwhm_threshold", 0.0);
    defaults_.setValue("deconvolution:fitting:eps_abs", 1e-05f, "if the absolute error gets smaller than this value the fitting is stopped.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("deconvolution:fitting:eps_abs", 0.0);
    defaults_.setValue("deconvolution:fitting:eps_rel", 1e-05f, "if the relative error gets smaller than this value the fitting is stopped.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("deconvolution:fitting:eps_rel", 0.0);
    defaults_.setValue("deconvolution:fitting:max_iteration", 10, "maximal number of iterations for the fitting step", ListUtils::create<String>("advanced"));
    defaults_.setMinInt("deconvolution:fitting:max_iteration", 1);

    //this->subsections_.push_back("SignalToNoiseEstimationParameter");
    SignalToNoiseEstimatorMeanIterative<MSSpectrum> sne;      // make sure this is the same as in pick()!
    Param param_sne_defaults =  sne.getDefaults();
    Param::ParamIterator param_it = param_sne_defaults.begin();
    for (; param_it != param_sne_defaults.end(); ++param_it)
    {
      if (!param_sne_defaults.hasTag(param_it.getName(), "advanced"))
        param_sne_defaults.addTag(param_it.getName(), "advanced");
    }
    this->defaults_.insert("SignalToNoiseEstimationParameter:", param_sne_defaults);

    defaultsToParam_();
  }

  PeakPickerCWT::~PeakPickerCWT()
  {
  }

  void PeakPickerCWT::updateMembers_()
  {
    signal_to_noise_ = (float)param_.getValue("signal_to_noise");
    peak_bound_ = (float)param_.getValue("thresholds:peak_bound");
    peak_bound_ms2_level_ = (float)param_.getValue("thresholds:peak_bound_ms2_level");
    scale_ = (float)param_.getValue("peak_width");
    fwhm_bound_ = (float)param_.getValue("fwhm_lower_bound_factor") * scale_;
    peak_corr_bound_ = (float)param_.getValue("thresholds:correlation");
    String opt = param_.getValue("optimization").toString();
    if (opt == "one_dimensional")
    {
      optimization_ = true;
      two_d_optimization_ = false;
    }
    else if (opt == "two_dimensional")
    {
      two_d_optimization_ = true;
      optimization_ = false;
    }
    else
    {
      optimization_ = false;
      two_d_optimization_ = false;
    }

    noise_level_ = (float)param_.getValue("thresholds:noise_level");
    radius_ = (Int)param_.getValue("thresholds:search_radius");
    signal_to_noise_ = (float)param_.getValue("signal_to_noise");

    deconvolution_ = param_.getValue("deconvolution:deconvolution").toBool();
  }

  bool PeakPickerCWT::getMaxPosition_(
    const PeakIterator first,
    const PeakIterator last,
    const ContinuousWaveletTransform & wt,
    PeakArea_ & area,
    const Int distance_from_scan_border,
    const double fxn_param_peak_bound,
    const double peak_bound_cwt,
    const Int direction) const
  {
    // ATTENTION! It is assumed that the resolution==1 (no resolution higher than 1).
    // Comment: Who cares ??
    const double noise_level = fxn_param_peak_bound;
    const double noise_level_cwt = peak_bound_cwt;

    const Int zeros_left_index  = wt.getLeftPaddingIndex();
    const Int zeros_right_index = wt.getRightPaddingIndex();

    // Points to most intensive data point in the signal
    double max_value;

    // Given direction, start the search from left or right
    Int start = (direction > 0) ? ((zeros_left_index + 2) + distance_from_scan_border) : ((zeros_right_index - 2) - distance_from_scan_border);
    Int end   = (direction > 0) ? (zeros_right_index - 1) : zeros_left_index + 1;

    Int i = 0, j = 0, k, max_pos;
    for (i = start, k = 0; i != end; i += direction, ++k)
    {
#ifdef DEBUG_PEAK_PICKING
      std::cout << "Search for max pos in cwt " <<  noise_level_cwt << std::endl;
      std::cout << wt[i - 1] << "\t" <<  wt[i] << "\t" <<  wt[i + 1] << "\n";
#endif

      // Check for maximum in cwt at position i
      if (((wt[i - 1] - wt[i]) < 0)
         && ((wt[i] - wt[i + 1]) > 0)
         && (wt[i]  >  noise_level_cwt))
      {
        max_pos = (direction > 0) ? (i - distance_from_scan_border) : i;
        if (first + max_pos < first || first + max_pos >= last)
          break;
#ifdef DEBUG_PEAK_PICKING
        std::cout << "MAX in CWT at " << (first + max_pos)->getMZ() << " with " << wt[i]
                  << std::endl;
#endif

        max_value = (first + max_pos)->getIntensity();


        // search for the corresponding maximum in the signal (consider the radius left and right adjacent points)
        Int start_intervall = ((max_pos - (Int)radius_) < 0) ? 0 : (max_pos - (Int)radius_);
        Int end_intervall = ((max_pos + (Int)radius_) >= distance(first, last)) ? 0 : (max_pos + (Int)radius_);

        for (j = start_intervall; j <= end_intervall; ++j)
        {
          if ((first + j)->getIntensity() > max_value)
          {
            max_pos = j;
            max_value = (first + j)->getIntensity();
          }

        }

        // if the maximum position is high enough and isn't one of the border points, we return it
        if (((first + max_pos)->getIntensity() >= noise_level)
           && (((first + max_pos) != first)
              && (first + max_pos)   != (last - 1)))
        {
          area.max = first + max_pos;
#ifdef DEBUG_PEAK_PICKING

          std::cout << "_________Max in data at__________ " << area.max->getMZ()
                    << std::endl;
#endif

          return true;
        }

#ifdef DEBUG_PEAK_PICKING
        std::cout << "NO REAL MAX!!!!!!!" << std::endl;
#endif

      }
    }

    // No relevant peak was found
    return false;
  }

  bool PeakPickerCWT::getPeakEndPoints_(PeakIterator first,
                                        PeakIterator last,
                                        PeakArea_ & area,
                                        Int distance_from_scan_border,
                                        Int & peak_left_index,
                                        Int & peak_right_index, ContinuousWaveletTransformNumIntegration & wt) const
  {
    // the Maximum may neither be the first or last point in the signal
    if ((area.max <= first) || (area.max >= last - 1))
    {
      return false;
    }

    PeakIterator it_help = area.max - 1;
    DPosition<1> vec_pos;
    Int cwt_pos;
    Int ep_radius = 2;
    Int start;
    Int stop;
    bool monoton;

    Int zeros_left_index  = wt.getLeftPaddingIndex();
    if (it_help != first)
    {
      // search for the left endpoint
      while (((it_help - 1) > first) && (it_help->getIntensity() > noise_level_))
      {

#ifdef DEBUG_PEAK_PICKING
        std::cout << "while left endpoint " << std::endl;
#endif
        // if the values are still falling to the left, everything is ok.
        if ((it_help - 1)->getIntensity() < it_help->getIntensity())
        {
          --it_help;

#ifdef DEBUG_PEAK_PICKING
          std::cout << "it_help " << it_help->getMZ() << std::endl;
#endif

        }
        // if the values are _rising_, we have to check the cwt
        else
        {
          if ((it_help - 2) <= first)
          {
#ifdef DEBUG_PEAK_PICKING
            std::cout << "it_help-2) <= first"  << std::endl;
#endif
            break;
          }
          // now check the value to the left of the problematic value
          if ((it_help - 2)->getIntensity() > (it_help - 1)->getIntensity())                   // we probably ran into another peak
          {
#ifdef DEBUG_PEAK_PICKING
            std::cout << "((it_help-2)->getIntensity() > (it_help-1)->getIntensity()"  << std::endl;
#endif
            // if peaks are broad and asymmetric, the max in the cwt can be a
            // few data points away from the real max
            // --> check also for distance in m/z in comparison with the scale
            if ((area.max - 1)->getMZ() - (it_help - 2)->getMZ() > scale_ / 2)
            {
              break;
            }

          }


          // to the left, the values are falling again => let the cwt decide if we
          // are seeing a new peak or just noise

          // compute the position of the corresponding point in the cwt
          cwt_pos = distance(first, it_help);
          vec_pos = it_help->getMZ();

          // since the cwt is pretty smooth usually, we consider the point as noise
          // if the cwt is monotonous in this region
          // TODO: better monotonicity test... say two or three points more
          monoton = true;

          start   =   (cwt_pos < ep_radius)
                    ? (distance_from_scan_border + zeros_left_index) + 2
                    : cwt_pos - ep_radius + (distance_from_scan_border + zeros_left_index + 2);
          stop    =   ((cwt_pos + ep_radius) > distance(it_help, last))
                    ? (wt.getSize() - 2)
                    : cwt_pos + ep_radius + (distance_from_scan_border + zeros_left_index + 2);

#ifdef DEBUG_PEAK_PICKING
          std::cout << "start " << start << " stop " << stop << std::endl;
#endif
          for (; start < stop; ++start)
          {
            if ((wt[start - 1] - wt[start])
                * (wt[start]   - wt[start + 1]) < 0)
            {
              // different slopes at the sides => stop here
#ifdef DEBUG_PEAK_PICKING
              std::cout << "monoton test " << wt.getSignal()[start - 1].getMZ()
                        << " " <<  wt.getSignal()[start].getMZ()
                        << " " <<  wt.getSignal()[start + 1].getMZ() << std::endl;
#endif
              monoton = false;
              break;
            }
          }

          if (!monoton)
          {
            break;
          }
          --it_help;
        }
      }
    }
    area.left = it_help;

    it_help = area.max + 1;
    if (it_help != last)
    {
      // search for the right endpoint ???
      while (((it_help + 1) < last) && (it_help->getIntensity() > noise_level_))
      {
#ifdef DEBUG_PEAK_PICKING
        std::cout << "while right endpoint " << std::endl;
#endif
        // if the values are still falling to the right, everything is ok.
        if (it_help->getIntensity() > (it_help + 1)->getIntensity())
        {
          ++it_help; // i.e. rightmost point is inclusive (not a past-end iterator)
#ifdef DEBUG_PEAK_PICKING
          std::cout << "it_help " << it_help->getMZ() << std::endl;
#endif
        }
        // if the values are _rising_, we have to check the cwt
        else
        {
          if ((it_help + 2) >= last)
          {
#ifdef DEBUG_PEAK_PICKING
            std::cout << "it_help+2) <= first"  << std::endl;
#endif
            break;
          }
          // now check the value to the right of the problematic value
          if ((it_help + 2)->getIntensity() > (it_help + 1)->getIntensity())                   // we probably ran into another peak
          {
#ifdef DEBUG_PEAK_PICKING
            std::cout << "(it_help+2)->getIntensity() > (it_help+1)->getIntensity())"  << std::endl;
#endif
            // if peaks are broad and asymmetric, the max in the cwt can be a few datapoints away from the real max
            // --> check also for distance in m/z in comparison with the scale
            if ((it_help + 2)->getMZ() - (area.max - 1)->getMZ() > scale_ / 2)
            {
              break;
            }
          }

          // to the left, the values are falling again => let the cwt decide if we
          // are seeing a new peak or just noise
          // compute the position of the corresponding point in the cwt
          cwt_pos = distance(first, it_help);
          //cwt_pos = distance(first, it_help);
          vec_pos = it_help->getMZ();

          // since the cwt is pretty smooth usually, we consider the point as noise
          // if the cwt is monotonous in this region
          // TODO: better monotonicity test... say two or three points more
          monoton = true;

          start   =   (cwt_pos < ep_radius)
                    ? (distance_from_scan_border + zeros_left_index) + 2
                    : cwt_pos - ep_radius + (distance_from_scan_border + zeros_left_index + 2);
          stop    =   ((cwt_pos + ep_radius) > distance(it_help, last))
                    ? (wt.getSize() - 2)
                    : cwt_pos + ep_radius + (distance_from_scan_border + zeros_left_index + 2);

#ifdef DEBUG_PEAK_PICKING
          std::cout << "start " << start << " stop " << stop << std::endl;
#endif
          for (; start < stop; ++start)
          {
            if ((wt[start - 1] - wt[start])
                * (wt[start]  - wt[start + 1]) < 0)
            {
              // different slopes at the sides => stop here
#ifdef DEBUG_PEAK_PICKING
              std::cout << "monoton test " << wt.getSignal()[start - 1].getMZ()
                        << " " <<  wt.getSignal()[start].getMZ()
                        << " " <<  wt.getSignal()[start + 1].getMZ() << std::endl;
#endif
              monoton = false;
              break;
            }
          }

          if (!monoton)
          {
            break;
          }
          ++it_help;
        }
      }
    }
    area.right = it_help;

    peak_left_index = distance(first, area.left);
    peak_right_index = distance(first, area.right);

    // The minimal profile data points per peak should be 2
    if ((distance(area.left, area.max) > 0) && (distance(area.max, area.right) > 0))
    {
#ifdef DEBUG_PEAK_PICKING
      std::cout << "OK !!!! The endpoints are "
                << area.left->getMZ()
                << " and "
                << area.right->getMZ()
                << std::endl;
#endif

      return true;
    }
#ifdef DEBUG_PEAK_PICKING
    std::cout << "NO!!! The endpoints are "
              << area.left->getMZ()
              << " and "
              << area.right->getMZ()
              << std::endl;
#endif

    return false;
  }

  void PeakPickerCWT::getPeakCentroid_(PeakArea_& area) const
  {
    PeakIterator left_it = area.max - 1, right_it = area.max;
    double rel_peak_height = area.max->getIntensity() * (double)param_.getValue("centroid_percentage");
    double sum = 0., w = 0.;

    // compute the centroid position (use weighted mean)
    while ((left_it >= area.left) && (left_it->getIntensity() >= rel_peak_height))
    {
      w += left_it->getIntensity() * left_it->getMZ();
      sum += left_it->getIntensity();
      --left_it;
    }

    // right point is inclusive (but should normally not be reached due to intensity cutoff)
    while ((right_it <= area.right) && (right_it->getIntensity() >= rel_peak_height))
    {
      w += right_it->getIntensity() * right_it->getMZ();
      sum += right_it->getIntensity();
      ++right_it;
    }

    area.centroid_position = w / sum;

#ifdef DEBUG_PEAK_PICKING
    std::cout << "________Centroid is___________ " << area.centroid_position << std::endl;
#endif

  }

  void PeakPickerCWT::initializeWT_(ContinuousWaveletTransformNumIntegration& wt, const double peak_bound_in, double& peak_bound_ms_cwt) const
  {
#ifdef DEBUG_PEAK_PICKING
    std::cout << "PeakPickerCWT<D>::initialize_ peak_bound_" << peak_bound_ <<  std::endl;
#endif
    //initialize wavelet transformer
    double spacing = (double)param_.getValue("wavelet_transform:spacing");
    wt.init(scale_, spacing);

    //calculate peak bound in CWT

    // build a Lorentz peak of height peak_bound_
    // compute its cwt, and compute the resulting height
    // of the transformed peak

    //compute the peak in the interval [-2*scale,2*scale]
    Int n = (Int)(scale_ / spacing * 4) + 1;

    double lambda = 2. / scale_;
    // compute the width parameter using height=peak_bound_ and the peak endpoints should be -scale and +scale, so at
    // positions -scale and +scale the peak value should correspond to the noise_level_
    //double lambda = sqrt((-noise_level_*(-peak_bound_+noise_level_)))/(noise_level_*scale_);

    MSSpectrum lorentz_peak;
    lorentz_peak.reserve(n);

    // TODO: switch the type of the transform

    ContinuousWaveletTransformNumIntegration lorentz_cwt;

    lorentz_cwt.init(scale_, spacing);
    double start = -2 * scale_;
    for (Int i = 0; i < n; ++i)
    {
      double p = i * spacing + start;
      MSSpectrum::value_type peak(p, lorentz_(peak_bound_in, lambda, 0, p));
      lorentz_peak.push_back(peak);
    }

    const float resolution = 1.;
    lorentz_cwt.transform(lorentz_peak.begin(), lorentz_peak.end(), resolution);

    peak_bound_ms_cwt = 0;
    for (Int i = 0; i < lorentz_cwt.getSignalLength(); ++i)
    {
      if (lorentz_cwt[i] > peak_bound_ms_cwt)
      {
        peak_bound_ms_cwt = lorentz_cwt[i];
      }
    }

#ifdef DEBUG_PEAK_PICKING
    std::cout << "PEAK BOUND IN CWT " << peak_bound_ms_cwt << std::endl;
#endif

  }

  void PeakPickerCWT::getPeakArea_(const PeakPickerCWT::PeakArea_& area, double& area_left, double& area_right) const
  {
    area_left = 0.0;
    // this does not depend on equal peak spacing
    for (PeakIterator pi = area.left; pi < area.max; ++pi)
    {
      PeakIterator pi_next = pi + 1;
      area_left += (pi->getIntensity() + pi_next->getIntensity()) / 2 // average intensity 
                   * (pi_next->getMZ() - pi->getMZ()); // m/z diff
    }

    // same with right side (reverse order to sum up small numbers first (numerical stability))
    area_right = 0.0;
    for (PeakIterator pi = area.right; pi > area.max; --pi)
    {
      PeakIterator pi_prev = pi - 1;
      area_right += (pi->getIntensity() + pi_prev->getIntensity()) / 2 // average intensity 
                    * (pi->getMZ() - pi_prev->getMZ()); // m/z diff
    }
  }

  PeakShape PeakPickerCWT::fitPeakShape_(const PeakPickerCWT::PeakArea_& area) const
  {

#ifdef DEBUG_PEAK_PICKING
    std::cout << "Left end point: " << area.left->getMZ() << std::endl;
#endif

    double max_intensity   =   area.max->getIntensity();
    double left_intensity  =  area.left->getIntensity();
    double right_intensity = area.right->getIntensity();

#ifdef DEBUG_PEAK_PICKING
    std::cout << "fit at the peak maximum " << std::endl;
#endif
    // compute peak areas from the left and right minima
    double peak_area_left, peak_area_right;
    getPeakArea_(area, peak_area_left, peak_area_right);

    // first the Lorentz peak ...
    // (see equation 8.14 on p. 74 in Dissertation of Eva Lange -- the equation has a typo: it should say A_r instead of A, i.e. lambda_r = h/A_r*arctan(...) )
    double left_width = max_intensity / peak_area_left * atan(sqrt(max_intensity / left_intensity - 1.));
    double right_width = max_intensity / peak_area_right * atan(sqrt(max_intensity / right_intensity - 1.));

    PeakShape lorentz(max_intensity, area.max->getMZ(),
                      left_width, right_width, peak_area_left + peak_area_right,
                      PeakShape::LORENTZ_PEAK);

    lorentz.r_value = correlate_(lorentz, area);

    // now the Sech2 peak ...
    // (see equation 8.17 on p. 74 in Dissertation of Eva Lange -- the equation has a typo: it should say A_r instead of A, i.e. lambda_r = h/A_r*sqrt(...) )
    left_width  = max_intensity / peak_area_left * sqrt(1. - left_intensity / max_intensity);
    right_width  = max_intensity / peak_area_right * sqrt(1. - right_intensity / max_intensity);

    PeakShape sech(max_intensity, area.max->getMZ(),
                    left_width, right_width,
                    peak_area_left + peak_area_right,
                    PeakShape::SECH_PEAK);

    sech.r_value = correlate_(sech, area);

#ifdef DEBUG_PEAK_PICKING

    std::cout << "r: " << lorentz.r_value << " " << sech.r_value << std::endl;
    std::cout << "pos: " << lorentz.mz_position <<  " " << sech.mz_position << std::endl;
    std::cout << "w1, w2: " << lorentz.left_width << " " << lorentz.right_width << " "
              << sech.left_width << " " << sech.right_width << std::endl;
    std::cout << "h: " << lorentz.height << std::endl;
#endif

    // take shape with higher correlation (Sech2 can be NaN, so Lorentzian might be the only option)
    if ((lorentz.r_value > sech.r_value) || boost::math::isnan(sech.r_value))
    {
      return lorentz;
    }
    else
    {
      return sech;
    }
  }

  bool PeakPickerCWT::deconvolutePeak_(PeakShape & shape, std::vector<PeakShape> & peak_shapes, double peak_bound_cwt) const
  {
    // scaling for charge one
    float scaling_DC = (float) param_.getValue("deconvolution:scaling");
    double resolution = 10;
    // init and calculate the transform of the signal in the convoluted region
    // first take the scaling for charge 2
    ContinuousWaveletTransformNumIntegration wtDC;
    wtDC.init(scaling_DC / 2, (double)param_.getValue("wavelet_transform:spacing"));
    wtDC.transform(shape.getLeftEndpoint(), shape.getRightEndpoint(), resolution);
    
#ifdef DEBUG_DECONV
    std::cout << "------------------\n---------------------\nconvoluted area begin " << shape.getLeftEndpoint()->getMZ() << "\tend " << shape.getRightEndpoint()->getMZ() << std::endl;
#endif


    Int charge = 2;
    std::vector<double> peak_values;
    std::vector<PeakShape> peaks_DC;
    Int peaks = getNumberOfPeaks_(shape.getLeftEndpoint(), shape.getRightEndpoint(), peak_values, 1, resolution, wtDC, peak_bound_cwt);

#ifdef DEBUG_PEAK_PICKING
    std::cout << "Number of peaks: " << peaks << std::endl;
#endif
    OptimizePeakDeconvolution::Data data;
    // one peak needn't be deconvoluted
    if (peaks > 1 && charge > 0)
    {

      data.positions.clear();
      data.signal.clear();
      data.peaks.clear();

      // enter zero-intensity at the left margin
      data.positions.push_back((shape.getLeftEndpoint())->getMZ() - 0.2);
      data.signal.push_back(0);

      for (Size i = 0; shape.getLeftEndpoint() + i != shape.getRightEndpoint(); ++i)
      {
        data.positions.push_back((shape.getLeftEndpoint() + i)->getMZ());
        data.signal.push_back((shape.getLeftEndpoint() + i)->getIntensity());
      }
      data.positions.push_back((shape.getRightEndpoint())->getMZ());
      data.signal.push_back((shape.getRightEndpoint())->getIntensity());

      data.positions.push_back((shape.getRightEndpoint())->getMZ() + 0.2);
      data.signal.push_back(0);



      // initial parameters for the optimization
      double leftwidth = (float)param_.getValue("deconvolution:left_width");
      double rightwidth = (float)param_.getValue("deconvolution:right_width");
      std::vector<double> cwt_distances(peaks - 1);
      peaks_DC.resize(peaks);
      for (Int i = 0; i < peaks; ++i)
      {
        PeakShape peak(peak_values[2 * i], peak_values[2 * i + 1], leftwidth, rightwidth, 0,
                       PeakShape::SECH_PEAK);
        peaks_DC[i] = peak;
        if (i < (peaks - 1))
        {
          cwt_distances[i] = fabs(peak_values[2 * i + 1] - peak_values[2 * (i + 1) + 1]);
        }

#ifdef DEBUG_DECONV
        std::cout << "peak.mz_position " << peak.mz_position << std::endl;
#endif
      }

      OptimizePeakDeconvolution opt;
      opt.setParameters(param_.copy("deconvolution:fitting:", true));
      opt.setCharge(charge);

#ifdef DEBUG_DECONV
      std::cout << "OptimizationType: Levenberg-Marquardt mit " << peaks_DC.size()
                << "peaks\n";
#endif

      //std::vector<PeakShape> peaks_DC2 = peaks_DC;
      //OptimizePeakDeconvolution::Data data2 = data;
      opt.optimize(peaks_DC, data);
      for (Int i = 0; i < peaks; ++i)
      {
        if (i < (peaks - 1))
        {
          if ((cwt_distances[i] - fabs(peaks_DC[i].mz_position - peaks_DC[i + 1].mz_position)) > 0.1)
          {
#ifdef DEBUG_DECONV
            std::cout << cwt_distances[i] << " - " << fabs(peaks_DC[i].mz_position - peaks_DC[i + 1].mz_position)
                      << " ---> cwt_distances not good" << std::endl;
#endif
            return false;
          }
        }
      }




#ifdef DEBUG_DECONV

      std::cout << " \nLM results:\n";

      //print results:
      for (Int i = 0; i < (Int)peaks_DC.size(); ++i)
      {
        std::cout << "\nposLM(" << i + 1 << ")=" << peaks_DC[i].mz_position;
        std::cout << "\nheightLM(" << i + 1 << ")=" << peaks_DC[i].height << std::endl;
        //  if(peaks_DC[i].type == PeakShape::LORENTZ_PEAK)
        //                          std::cout<<"\npeak("<<i+1<<")=0\n";
        //                      else
        //                          std::cout<<"\npeak("<<i+1<<")=1\n";
      }
      //                std::cout<<"\nlw="<<peaks_DC[0].left_width;
      //                std::cout<<"\nrw="<<peaks_DC[0].right_width<<"\n\n";

#endif
      // add the new peaks
      for (Size curr_peak = 0; curr_peak < peaks_DC.size(); ++curr_peak)
      {
        peak_shapes.push_back(peaks_DC[curr_peak]);
      }

      data.peaks.clear();
      data.signal.clear();
      data.positions.clear();
      peaks_DC.clear();
    }
    else
      return false;

    return true;

  }

  void PeakPickerCWT::addPeak_(std::vector<PeakShape> & peaks_DC, PeakArea_ & area, double left_width, double right_width, OptimizePeakDeconvolution::Data & data) const
  {
    // just enter a peak using equally spaced peak positions

    double peak_width = area.right->getMZ() - area.left->getMZ();
    Size num_peaks = peaks_DC.size() + 1;

    double dist = peak_width / (num_peaks + 1);

    // put peak into peak vector using default values for the widths and peak type
    peaks_DC.push_back(PeakShape(0, 0, left_width, right_width, 0, PeakShape::SECH_PEAK));

    // adjust the positions and get their initial intensities from the profile data
    for (Size i = 0; i < num_peaks; ++i)
    {
      peaks_DC[i].mz_position = area.left->getMZ() + dist / 2 + i * dist;

      std::vector<double>::iterator it_help = lower_bound(data.positions.begin(),
                                                              data.positions.end(),
                                                              peaks_DC[i].mz_position);
      if (it_help != data.positions.end())
      {
        peaks_DC[i].height =
          data.signal[distance(data.positions.begin(), it_help)] / 10;
#ifdef DEBUG_DECONV
        std::cout << "height " << i << "   " <<  peaks_DC[i].height << "\t"
                  << distance(data.positions.begin(), it_help) << "\n";
#endif
      }
      else
      {
        peaks_DC[i].height =
          data.signal[data.positions.size() - 1];
#ifdef DEBUG_DECONV
        std::cout << "else height " << i << "   " << peaks_DC[i].height << "\n";
#endif
      }
    }


  }

  Int PeakPickerCWT::getNumberOfPeaks_(ConstPeakIterator first,
                                       ConstPeakIterator last,
                                       std::vector<double> & peak_values,
                                       Int direction,
                                       double resolution,
                                       ContinuousWaveletTransformNumIntegration & wt,
                                       double peak_bound_cwt) const
  {
    double noise_level = peak_bound_;
    double noise_level_cwt = peak_bound_cwt;

#ifdef DEBUG_DECONV
    std::cout << "noise_level = " << noise_level << ";\tnoise_level_cwt = " << noise_level_cwt << ";\n";
    std::cout << "resoltuion " << resolution << "\n";
    std::cout << "first und last " << first->getMZ() << "\t" << last->getMZ() << std::endl;

#endif

    Int found = 0;

    Int zeros_left_index  = wt.getLeftPaddingIndex();
    Int zeros_right_index = wt.getRightPaddingIndex();

    //double max_value;T
    Int start = (direction > 0) ? zeros_left_index + 2 : zeros_right_index - 2;
    Int end   = (direction > 0) ? zeros_right_index - 1 : zeros_left_index + 1;

    Int i = 0, max_pos;
    Int k = 0;

    //std::vector<double>::iterator checker;
    while (wt.getSignal()[start + 1].getMZ() <= first->getMZ())
      ++start;
    //k=i;
    Int offset = start;
    while (wt.getSignal()[end].getMZ() > last->getMZ())
      --end;
    for (i = start; i != end; i += direction, k += direction)
    {

      // Does a maximum in the wavelet transform at position i exists?
      if (((wt[i - 1] - wt[i]) < 0)
         && ((wt[i]   - wt[i + 1]) > 0)
         && (wt[i]   >  noise_level_cwt))
      {
#ifdef DEBUG_DECONV
        std::cout << "MAX in CWT at " << (first + (i - offset) / resolution)->getPosition()[0] << " with " << wt[i]
                  << std::endl;

#endif
        max_pos = Int((i - offset) / resolution);

        // if the maximum position is high enough and isn't one of the border points, we return it
        if (((first + max_pos)->getIntensity() >= noise_level)
           && (((first + max_pos) != first)
              && (first + max_pos)   != (last - 1))
            )
        {
#ifdef DEBUG_DECONV
          std::cout << "_________REAL Max in CWT at__________ " << (first + max_pos)->getPosition()[0] << "\n";
#endif
          peak_values.push_back((first + max_pos)->getIntensity());
          peak_values.push_back((first + max_pos)->getMZ());

          ++found;
        }

      }
    }
#ifdef DEBUG_DECONV
    std::cout << "number of peaks found: " << found << "\n";
#endif
    return found;
  }

  Int PeakPickerCWT::determineChargeState_(std::vector<double> & peak_values) const
  {
    Int charge;
    Int peaks = (Int)peak_values.size() / 2;
    if (peaks > 1)
    {
      double dif = 0;
      Int i = peaks - 1;
      while (i > 0)
      {
        dif += fabs(peak_values[2 * i + 1] - peak_values[2 * (i - 1) + 1]);
#ifdef DEBUG_DECONV
        std::cout << "peak_values[" << i << "] " << peak_values[2 * i + 1]
                  << " peak_values[" << i - 1 << "] " << peak_values[2 * (i - 1) + 1]
                  << " dif[" << i << "]=" << dif << std::endl;
#endif
        --i;
      }
      dif /= peaks - 1;
      charge = (Int) Math::round(1 / dif);
      if (boost::math::isnan((double)charge) || boost::math::isinf((double)charge))
        charge = 0;
#ifdef DEBUG_DECONV
      std::cout << "1/dif = " << 1 / dif << ";\tcharge = " << charge << std::endl;
#endif
    }
    else
      charge = 1;

    return charge;
  }

  double PeakPickerCWT::correlate_(const PeakShape & peak,
                                   const PeakPickerCWT::PeakArea_ & area,
                                   Int direction) const
  {
    // compute the averages
    double data_average = 0., fit_average = 0.;
    double data_sqr = 0., fit_sqr = 0.;
    double cross = 0.;

    Int number_of_points = 0;
    PeakIterator corr_begin = area.left;
    PeakIterator corr_end = area.right;

    // for separate overlapping peak correlate until the max position...
    if (direction > 0)
      corr_end = area.max;
    else if (direction < 0)
      corr_begin = area.max;

    for (PeakIterator pi = corr_begin; pi <= corr_end; ++pi)
    {
      double data_val = pi->getIntensity();
      double peak_val = peak(pi->getMZ());

      data_average += data_val;
      fit_average  += peak_val;

      data_sqr += data_val * data_val;
      fit_sqr  += peak_val * peak_val;

      cross += data_val * peak_val;

      ++number_of_points;
    }

    if (number_of_points == 0)
      return 0.;

    data_average /= number_of_points;
    fit_average  /= number_of_points;

    double SSxx = data_sqr - number_of_points * (data_average * data_average);
    double SSyy = fit_sqr - number_of_points * (fit_average * fit_average);
    double SSxy = cross - number_of_points * (data_average * fit_average);
    // using r^2 = cov(x,y) / (var(x) * var(y))
    return (SSxy * SSxy) / (SSxx * SSyy);
  }

  void PeakPickerCWT::pickExperiment(const PeakMap & input, PeakMap & output)
  {
    // if estimatePeakWidth-flag is set estimate it
    if (param_.getValue("estimate_peak_width") == "true")
    {
      double p_w = estimatePeakWidth(input);
      if (p_w == 0.)
      {
        throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "estimatePeakWidth()", "Peak width could not be determined from data!");
      }
      else
      {
        param_.setValue("peak_width", p_w);
        updateMembers_();
      }
    }
    // clear output container
    output.clear(true);

    // copy the experimental settings
    static_cast<ExperimentalSettings &>(output) = input;
    output.resize(input.size());
    // pick peaks on each scan
    startProgress(0, input.size(), "picking peaks");
    Size progress = 0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize i = 0; i < (SignedSize)input.size(); ++i)
    {
      // pick the peaks in scan i
      // this is needed to eliminate empty spectra in the end
      pick(input[i], output[i]);
#ifdef _OPENMP
#pragma omp critical (PeakPickerCWT_PickExperiment)
#endif
      {
        setProgress(++progress);             //do not use 'i' here, as each thread will be assigned different blocks
      }
    }
    //optimize peak positions
    if (two_d_optimization_ || optimization_)
    {
      TwoDOptimization my_2d;
      my_2d.setParameters(param_.copy("optimization:", true));
      my_2d.optimize(input.begin(), input.end(), output, two_d_optimization_);
    }
    endProgress();

  }

  void PeakPickerCWT::pick(const MSSpectrum & input, MSSpectrum & output) const
  {
    // copy the spectrum meta data
    output.clear(true);
    output.SpectrumSettings::operator=(input);
    output.MetaInfoInterface::operator=(input);
    output.setRT(input.getRT());
    output.setMSLevel(input.getMSLevel());
    output.setName(input.getName());
    //make sure the data type is set correctly
    output.setType(SpectrumSettings::CENTROID);

    // nearly empty spectra shouldn't be picked
    if (input.size() < 2)
      return;

    //set up meta data arrays
    output.getFloatDataArrays().clear();
    output.getFloatDataArrays().resize(7);
    output.getFloatDataArrays()[0].setName("rValue");
    output.getFloatDataArrays()[1].setName("maximumIntensity");
    output.getFloatDataArrays()[2].setName("fwhm");
    output.getFloatDataArrays()[3].setName("leftWidth");
    output.getFloatDataArrays()[4].setName("rightWidth");
    output.getFloatDataArrays()[5].setName("peakShape");
    output.getFloatDataArrays()[6].setName("SignalToNoise");

    /// The continuous wavelet "transformer"
    ContinuousWaveletTransformNumIntegration wt;
    /// The minimal height which defines a peak in the CWT
    double peak_bound_ms_cwt = 0.0;
    double bound = (input.getMSLevel() <= 1 ? peak_bound_ : peak_bound_ms2_level_);
    // now initialize every time as every spectrum is picked with its own cwt
    initializeWT_(wt, bound, peak_bound_ms_cwt);

    //create the peak shapes vector
    std::vector<PeakShape> peak_shapes;

#ifdef DEBUG_PEAK_PICKING
    std::cout << "****************** PICK ******************" << input.getRT()
              << " " << input.getMSLevel() << std::endl;
#endif

    // vector of peak endpoint positions
    // std::vector<double> peak_endpoints;

    // copy the profile data into a std::vector<Peak1D>
    MSSpectrum raw_peak_array;
    // signal to noise estimator
    SignalToNoiseEstimatorMeanIterative<MSSpectrum> sne;
    Param sne_param(param_.copy("SignalToNoiseEstimationParameter:", true));
    sne.setParameters(sne_param);

    raw_peak_array.insert(raw_peak_array.end(), input.begin(), input.end());

    PeakIterator it_pick_begin = raw_peak_array.begin();
    PeakIterator it_pick_end   = raw_peak_array.end();

    sne.init(it_pick_begin, it_pick_end);

    // Upper peak width bound
    double fwhm_upper_bound = (double)param_.getValue("fwhm_upper_bound_factor") * scale_;

    // thresholds for deconvolution
    double fwhm_threshold = (float)param_.getValue("deconvolution:fitting:fwhm_threshold");
    double symm_threshold = (float)param_.getValue("deconvolution:asym_threshold");

    // start the peak picking until no more maxima can be found in the wavelet transform
    UInt number_of_peaks = 0;

    do
    {
      number_of_peaks = 0;
      Int peak_left_index, peak_right_index;

      // compute the continuous wavelet transform with resolution 1
      const double resolution = 1;
      wt.transform(it_pick_begin, it_pick_end, resolution);
      PeakArea_ area;

      // search for maximum positions in the cwt and extract potential peaks
      Int direction = 1;
      Int distance_from_scan_border = 0;
      while ((distance(it_pick_begin, it_pick_end) > 3)
            && getMaxPosition_(it_pick_begin,
                               it_pick_end,
                               wt,
                               area,
                               distance_from_scan_border,
                               bound,
                               peak_bound_ms_cwt,
                               direction))
      {
        // if the signal to noise ratio at the max position is too small
        // the peak isn't considered
        if ((area.max != it_pick_end) && (sne.getSignalToNoise(area.max) < signal_to_noise_))
        {
          it_pick_begin = area.max;
          distance_from_scan_border = distance(raw_peak_array.begin(), it_pick_begin);
          continue;
        }
        else if (area.max >= it_pick_end)
        {
          break;
        }
        //search for the endpoints of the peak
        bool regular_endpoints = getPeakEndPoints_(it_pick_begin,
                                                  it_pick_end,
                                                  area,
                                                  distance_from_scan_border,
                                                  peak_left_index,
                                                  peak_right_index, wt);

        // compute the centroid position (area.centroid_position)
        getPeakCentroid_(area);

        // if the peak achieves a minimal width, start the peak fitting
        if (regular_endpoints)
        {
          // #ifdef DEBUG_PEAK_PICKING
          // std::cout << "The endpoints are "
          // << area.left->getPosition()
          // << " and "
          // << area.right->getPosition()
          // << std::endl;
          // #endif
          // determine the best fitting Lorentzian or sech2 function
          PeakShape shape = fitPeakShape_(area);
          shape.setLeftEndpoint((input.begin() + distance(raw_peak_array.begin(), area.left)));
          shape.setRightEndpoint((input.begin() + distance(raw_peak_array.begin(), area.right)));
          if (shape.getRightEndpoint() == input.end())
            shape.setRightEndpoint(input.end() - 1);
          // Use the centroid for Optimization
          shape.mz_position = area.centroid_position[0];
          if ((shape.r_value > peak_corr_bound_)
             && (shape.getFWHM() >= fwhm_bound_)
             && (shape.getFWHM() <= fwhm_upper_bound))
          {
            shape.signal_to_noise = sne.getSignalToNoise(area.max);
            peak_shapes.push_back(shape);
            ++number_of_peaks;
          }
          else
          {
            // #ifdef DEBUG_PEAK_PICKING
            // std::cout << "Corr: " << shape.r_value << " SN: " << sne.getSignalToNoise(area.max) << " FWHM: " << shape.getFWHM() << std::endl;
            // std::cout << "Bad fitting peak "<< std::endl;
            // #endif
          }
        }

        // remove the peak from the signal
        // TODO: does this work as expected???
        for (PeakIterator pi = area.left; pi != area.right + 1; ++pi)
        {
          pi->setIntensity(0.);
        }

        // search for the next peak
        it_pick_begin = area.right;
        distance_from_scan_border = distance(raw_peak_array.begin(), it_pick_begin);

      }               //end while (getMaxPosition_(it_pick_begin, it_pick_end, wt, area, distance_from_scan_border, ms_level, direction))
      it_pick_begin = raw_peak_array.begin();
    }
    while (number_of_peaks != 0);

    // start the nonlinear optimization for all peaks in split
#ifdef DEBUG_PEAK_PICKING
    std::cout << "Try the optimization run... with " << peak_shapes.size() << std::endl;
#endif

    if (!peak_shapes.empty())
    {
      // overlapping peaks are mostly broad or asymmetric
      // we distinguish them from broad or asymmetric isotopic peaks
      // (e.g. charge one peaks, or peaks in the high mass range)
      // by a simple heuristic: if the distances to adjacent peaks
      // are dissimilar, the FWHM is much broader than the FWHM of
      // adjacent peaks or if the peak has no near neighbors
      // we assume a convolved peak pattern and start the deconvolution.
      // sort the peaks according to their positions
      sort(peak_shapes.begin(), peak_shapes.end(), PeakShape::PositionLess());
      std::set<UInt> peaks_to_skip;
      // search for broad or asymmetric peaks
      if (deconvolution_)
      {
        UInt n = (UInt) peak_shapes.size();
        for (UInt i = 0; i < n; ++i)
        {
          if ((peak_shapes[i].getFWHM() > fwhm_threshold)
             || (peak_shapes[i].getSymmetricMeasure() < symm_threshold))
          {
#ifdef DEBUG_DECONV
            std::cout << "check " << peak_shapes[i].mz_position
                      << " with fwhm: " << peak_shapes[i].getFWHM()
                      << " and " << peak_shapes[i].left_width
                      << ' ' << peak_shapes[i].right_width
                      << ' ' << peak_shapes[i].getLeftEndpoint()->getMZ()
                      << ' ' << peak_shapes[i].getRightEndpoint()->getMZ()
                      << std::endl;
#endif
            // this might be a convolved peak pattern
            // and we check the distance to the neighboring peaks as well as their fwhm values:
            //float max_distance = 1.1;
            float dist_left = ((i > 0) && (fabs(peak_shapes[i].mz_position - peak_shapes[i - 1].mz_position) < 1.2)) ? fabs(peak_shapes[i].mz_position - peak_shapes[i - 1].mz_position) : -1;
            float dist_right = ((i < (n - 1)) && (fabs(peak_shapes[i].mz_position - peak_shapes[i + 1].mz_position) < 1.2)) ? fabs(peak_shapes[i].mz_position - peak_shapes[i + 1].mz_position) : -1;

            // left and right neighbor
            if ((dist_left > 0) && (dist_right > 0))
            {
              // if distances to left and right adjacent peaks is dissimilar deconvolute
              double ratio = (dist_left > dist_right) ? dist_right / dist_left : dist_left / dist_right;
#ifdef DEBUG_DECONV
              std::cout << "Ratio " << ratio << std::endl;
#endif
              if (ratio < 0.6)
              {
#ifdef DEBUG_DECONV
                std::cout << "deconvolute: dissimilar left and right neighbor "  << peak_shapes[i - 1].mz_position << ' ' << peak_shapes[i + 1].mz_position << std::endl;
#endif
                if (deconvolutePeak_(peak_shapes[i], peak_shapes, peak_bound_ms_cwt))
                  peaks_to_skip.insert(i);
              }
            }
            // has only one or no neighbor peak
            else
            {
              // only left neighbor
              if (dist_left > 0)
              {
                // check distance and compare fwhm
                double dist = 1.00235;
                //check charge 1 or 2
                bool dist_ok = ((fabs(dist - dist_left) < 0.21) || (fabs(dist / 2. - dist_left) < 0.11)) ? true : false;
                // distance complies peptide mass rule
                if (dist_ok)
                {
#ifdef DEBUG_DECONV
                  std::cout << "left neighbor " << peak_shapes[i - 1].mz_position << ' ' << peak_shapes[i - 1].getFWHM() << std::endl;
#endif
                  // if the left peak has a fwhm which is smaller than 60% of the fwhm of the broad peak deconvolute
                  if ((peak_shapes[i - 1].getFWHM() / peak_shapes[i].getFWHM()) < 0.6)
                  {
#ifdef DEBUG_DECONV
                    std::cout << " too small fwhm" << std::endl;
#endif
                    if (deconvolutePeak_(peak_shapes[i], peak_shapes, peak_bound_ms_cwt))
                      peaks_to_skip.insert(i);
                  }
                }
                else
                {
#ifdef DEBUG_DECONV
                  std::cout << "distance not ok" << dist_left << ' ' << peak_shapes[i - 1].mz_position << std::endl;
#endif
                  if (deconvolutePeak_(peak_shapes[i], peak_shapes, peak_bound_ms_cwt))
                    peaks_to_skip.insert(i);
                }
              }
              else
              {
                // only right neighbor
                if (dist_right > 0)
                {
                  // check distance and compare fwhm
                  double dist = 1.00235;
                  //check charge 1 or 2
                  bool dist_ok = ((fabs(dist - dist_right) < 0.21) || (fabs(dist / 2. - dist_right) < 0.11)) ? true : false;
                  // distance complies with peptide mass rule
                  if (dist_ok)
                  {
#ifdef DEBUG_DECONV
                    std::cout << "right neighbor " << peak_shapes[i + 1].mz_position << ' ' << peak_shapes[i + 1].getFWHM() << std::endl;
#endif
                    // if the left peak has a fwhm which is smaller than 60% of the fwhm of the broad peak deconvolute
                    if ((peak_shapes[i + 1].getFWHM() / peak_shapes[i].getFWHM()) < 0.6)
                    {
#ifdef DEBUG_DECONV
                      std::cout << "too small fwhm"  << std::endl;
#endif
                      if (deconvolutePeak_(peak_shapes[i], peak_shapes, peak_bound_ms_cwt))
                        peaks_to_skip.insert(i);
                    }
                  }
                  else
                  {
#ifdef DEBUG_DECONV
                    std::cout << "distance not ok" << dist_right << ' ' << peak_shapes[i + 1].mz_position << std::endl;
#endif
                    if (deconvolutePeak_(peak_shapes[i], peak_shapes, peak_bound_ms_cwt))
                      peaks_to_skip.insert(i);
                  }
                }
                // no neighbor
                else
                {
#ifdef DEBUG_DECONV
                  std::cout << "no neighbor" << std::endl;
#endif
                  if (deconvolutePeak_(peak_shapes[i], peak_shapes, peak_bound_ms_cwt))
                    peaks_to_skip.insert(i);
                }
              }
            }
          }
        }
      }

      //reserve space in the output container
      Size npeaks = peak_shapes.size() - peaks_to_skip.size();
      output.reserve(npeaks);
      output.getFloatDataArrays()[0].reserve(npeaks);
      output.getFloatDataArrays()[1].reserve(npeaks);
      output.getFloatDataArrays()[2].reserve(npeaks);
      output.getFloatDataArrays()[3].reserve(npeaks);
      output.getFloatDataArrays()[4].reserve(npeaks);
      output.getFloatDataArrays()[5].reserve(npeaks);
      output.getFloatDataArrays()[6].reserve(npeaks);

      // write the picked peaks to the output container
      for (Size i = 0; i < peak_shapes.size(); ++i)
      {
        // put it out only if the peak was not deconvoluted
        if (peaks_to_skip.find((UInt)i) == peaks_to_skip.end())
        {
          //store output peak
          Peak1D picked_peak;
          // store area as intensity, the maximal intensity is stored in the metadata arrays
          picked_peak.setIntensity(peak_shapes[i].area);
          picked_peak.setMZ(peak_shapes[i].mz_position);
          output.push_back(picked_peak);
          //store meta data
          output.getFloatDataArrays()[0].push_back(peak_shapes[i].r_value);
          output.getFloatDataArrays()[1].push_back(peak_shapes[i].height);
          output.getFloatDataArrays()[2].push_back(peak_shapes[i].getFWHM());
          output.getFloatDataArrays()[3].push_back(peak_shapes[i].left_width);
          output.getFloatDataArrays()[4].push_back(peak_shapes[i].right_width);
          output.getFloatDataArrays()[5].push_back(peak_shapes[i].type);
          output.getFloatDataArrays()[6].push_back(peak_shapes[i].signal_to_noise);

#ifdef DEBUG_PEAK_PICKING
          if (peak_shapes[i].getFWHM() > 1.)
          {
            std::cout << output.getRT() << "\t" << peak_shapes[i].mz_position << "\t"
                      << peak_shapes[i].getFWHM() << "\t" << peak_shapes[i].height << std::endl;
          }
#endif

        }
      }
    }         // if (peak_shapes.size() > 0)
  }

  double PeakPickerCWT::estimatePeakWidth(const PeakMap & input)
  {
    // The peak widths which are tested
    // Note that the limit is 1 Da in Peak Width. The Wavelet has a certain tolerance, but very broad peaks (e.g. top-down MS) will not be found
    // unless the search range is extended
    DoubleList test_widths = ListUtils::create<double>("1.,0.5,0.25,0.125,0.1,0.05,0.025,0.0125,0.01,0.005,0.0025,0.00125,0.0005,0.0001");
#ifdef DEBUG_PEAK_PICKING2
    std::cout << "calculating max tic" << std::endl;
#endif
    // make a backup of param .. restore later
    Param p_backup = param_;

    // use a high S/N to ensure we only look at high peaks.. low noise peaks might distort the FWHM
    param_.setValue("signal_to_noise", 10); 

    // determine spectrum used for estimation
    // use the one with the highest tic
    TICFilter tic_filter;
    typedef std::vector<std::pair<Size, double> > TicVec;
    TicVec index_tic_vec;
    for (Size s = 0; s < input.size(); ++s)
    {
      double tmp_tic = tic_filter.apply(input[s]);
      index_tic_vec.push_back(make_pair(s, tmp_tic));
    }

    // now sort the index_tic_vec according to tic and get the three highest spectra
    sort(index_tic_vec.begin(), index_tic_vec.end(), PairComparatorSecondElement<std::pair<Size, double> >());
    
    std::vector<double> best_FWHMs;
    for (Size vec_idx = 0; vec_idx < 3 && vec_idx < index_tic_vec.size(); ++vec_idx)
    {
      std::pair<Size, double> tic = index_tic_vec[vec_idx];
      // skip empty spectra
      if (tic.second == 0) continue;

      std::vector<Size> peak_counts;
#ifdef DEBUG_PEAK_PICKING2
      std::cout << "RT: " << input[tic.first].getRT() << "\n";
#endif
      // holds average FWHM of this spectrum
      std::vector<double> observed_FWHM;
      // ... for different peak widths
      for (Size w = 0; w < test_widths.size(); ++w)
      {
        MSSpectrum spec;
        param_.setValue("peak_width", test_widths[w]);
        updateMembers_();
#ifdef DEBUG_PEAK_PICKING2
        std::cout << "peak_width " << param_.getValue("peak_width") << "\tfwhm_bound_ " << fwhm_bound_ << "\t";
#endif
        pick(input[tic.first], spec);

        double fwhm = 0.;
        for (Size p = 0; p < spec.size(); ++p)
        {
          fwhm += spec.getFloatDataArrays()[2][p];
        }
        fwhm /= (double) spec.size();
        observed_FWHM.push_back(fwhm);

        peak_counts.push_back(spec.size());
      }

      // determine slope
      std::vector<double> slopes;
      double m_min = std::numeric_limits<double>::max();
      for (Size i = 0; i < peak_counts.size() - 1; ++i)
      {
        double m = (double)((SignedSize)peak_counts[i] - (SignedSize)peak_counts[i + 1]) / (test_widths[i] - test_widths[i + 1]);
        slopes.push_back(m);
#ifdef DEBUG_PEAK_PICKING2
        std::cout << m << "\t";
#endif
        if (m < m_min)
        {
          m_min = m;
        }
#ifdef DEBUG_PEAK_PICKING2
        std::cout << (test_widths[i] + test_widths[i + 1]) / 2 << "\t" << m << std::endl;
#endif
      }
      // determine point where slope starts decreasing , there the plateau should end
      bool min_found = false;
#ifdef DEBUG_PEAK_PICKING2
      std::cout << "m_min: " << m_min << std::endl;
#endif
      for (Size s = 0; s < slopes.size(); ++s)
      {
#ifdef DEBUG_PEAK_PICKING2
        std::cout << slopes[s] << std::endl;
#endif
        if (fabs(slopes[s] - m_min) < 0.01)
        {
          min_found = true;
        }
        if (min_found && slopes[s] / m_min < 0.5)
        {
          double pw = observed_FWHM[s+1];
#ifdef DEBUG_PEAK_PICKING2
          std::cout << "peak_width = " << pw << std::endl;
#endif
          best_FWHMs.push_back(pw);
          break;
        }
      }
    }
    // average width over the three tested spectra
    double avg_width = 0.;
    if (best_FWHMs.empty())
    {
      std::cout << "Could not estimate peak width!" << std::endl;
    } 
    else
    {
      for (Size s = 0; s < best_FWHMs.size(); ++s)
      {
        std::cout << "  Estimated FWHM of spectrum #" << s << ": " << best_FWHMs[s] << std::endl;
        avg_width += best_FWHMs[s];
      }
      avg_width /= (double)best_FWHMs.size();
      std::cout << "--> Average estimated FWHM: " << avg_width << std::endl;
    }

    // restore old param set
    param_ = p_backup; 

    return avg_width;
  }

}
