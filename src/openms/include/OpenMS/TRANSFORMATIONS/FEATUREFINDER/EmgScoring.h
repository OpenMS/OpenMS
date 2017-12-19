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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_EMGSCORING_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_EMGSCORING_H

#include <vector>
#include <boost/math/special_functions/fpclassify.hpp> // for isnan
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgModel.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>

#include <OpenMS/KERNEL/MRMFeature.h>
#include <OpenMS/KERNEL/MRMTransitionGroup.h>

#include <OpenMS/KERNEL/StandardTypes.h>


namespace OpenMS
{

  /**
    @brief Scoring of an elution peak using an exponentially modified gaussian distribution model.

    This class uses the original ideas from FeatureFinderAlgorithmMRM to
    construct an interface that allows scoring of chromatographic peaks.

  */
  class EmgScoring
  {

  public :

    EmgScoring() { }

    ~EmgScoring() { }

    void setFitterParam(Param param)
    {
      fitter_emg1D_.setParameters(param);
    }

    Param getDefaults()
    {
      return fitter_emg1D_.getDefaults();
    }

    /// calculate the elution profile fit score
    template<typename SpectrumType, class TransitionT>
    double calcElutionFitScore(MRMFeature & mrmfeature, MRMTransitionGroup<SpectrumType, TransitionT> & transition_group)
    {
      std::vector<double> fit_scores;
      double avg_score = 0;
      bool smooth_data = false;
      for (Size k = 0; k < transition_group.size(); k++)
      {
        // get the id, then find the corresponding transition and features within this peakgroup
        String native_id = transition_group.getChromatograms()[k].getNativeID();
        Feature f = mrmfeature.getFeature(native_id);
        OPENMS_PRECONDITION(f.getConvexHulls().size() == 1, "Convex hulls need to have exactly one hull point structure");

        // TODO what if score is -1 ?? e.g. if it is undefined
        double fscore = elutionModelFit(f.getConvexHulls()[0].getHullPoints(), smooth_data);
        fit_scores.push_back(fscore);
        avg_score += fscore;
      }

      avg_score /= transition_group.size();
      return avg_score;
    }

    // Fxn from FeatureFinderAlgorithmMRM
    // TODO: check whether we can leave out some of the steps here, e.g. gaussian smoothing
    double elutionModelFit(ConvexHull2D::PointArrayType current_section, bool smooth_data)
    {
      // We need at least 2 datapoints in order to create a fit
      if (current_section.size() < 2)
      {
        return -1;
      }

      // local PeakType is a small hack since here we *need* data of type
      // Peak1D, otherwise our fitter will not accept it.
      typedef Peak1D LocalPeakType;

      // -- cut line 301 of FeatureFinderAlgorithmMRM
      std::vector<LocalPeakType> data_to_fit;
      prepareFit_(current_section, data_to_fit, smooth_data);
      InterpolationModel * model_rt = nullptr;
      double quality = fitRT_(data_to_fit, model_rt);
      // cut line 354 of FeatureFinderAlgorithmMRM
      delete model_rt;

      return quality;

    }

  protected:
    template<class LocalPeakType>
    double fitRT_(std::vector<LocalPeakType> & rt_input_data, InterpolationModel * & model)
    {
      double quality;
      //Param param;

      /*EmgFitter
       param.setValue( "tolerance_stdev_bounding_box", tolerance_stdev_box_);
       param.setValue( "statistics:mean", rt_stat_.mean() );
       param.setValue( "statistics:variance", rt_stat_.variance() );
       param.setValue( "interpolation_step", interpolation_step_rt_ );
       param.setValue( "max_iteration", max_iteration_);
       param.setValue( "deltaAbsError", deltaAbsError_);
       param.setValue( "deltaRelError", deltaRelError_);
       */

      // Set parameter for fitter
      //fitter_emg1D.setParameters(param);
      // Construct model for rt
      quality = fitter_emg1D_.fit1d(rt_input_data, model);

      // Check quality
      if (boost::math::isnan(quality)) quality = -1.0;
      return quality;
    }

    // Fxn from FeatureFinderAlgorithmMRM
    // TODO: check whether we can leave out some of the steps here, e.g. gaussian smoothing
    template<class LocalPeakType>
    void prepareFit_(const ConvexHull2D::PointArrayType & current_section, std::vector<LocalPeakType> & data_to_fit, bool smooth_data)
    {
      // typedef Peak1D LocalPeakType;
      PeakSpectrum filter_spec;
      // first smooth the data to prevent outliers from destroying the fit
      for (ConvexHull2D::PointArrayType::const_iterator it = current_section.begin(); it != current_section.end(); it++)
      {
        LocalPeakType p;
        p.setMZ(it->getX());
        p.setIntensity(it->getY());
        filter_spec.push_back(p);
      }

      // add two peaks at the beginning and at the end for better fit
      // therefore calculate average distance first
      std::vector<double> distances;
      for (Size j = 1; j < filter_spec.size(); ++j)
      {
        distances.push_back(filter_spec[j].getMZ() - filter_spec[j - 1].getMZ());
      }
      double dist_average = std::accumulate(distances.begin(), distances.end(), 0.0) / (double) distances.size();

      // append peaks
      Peak1D new_peak;
      new_peak.setIntensity(0);
      new_peak.setMZ(filter_spec.back().getMZ() + dist_average);
      filter_spec.push_back(new_peak);
      new_peak.setMZ(filter_spec.back().getMZ() + dist_average);
      filter_spec.push_back(new_peak);
      new_peak.setMZ(filter_spec.back().getMZ() + dist_average);
      filter_spec.push_back(new_peak);

      // prepend peaks
      new_peak.setMZ(filter_spec.front().getMZ() - dist_average);
      filter_spec.insert(filter_spec.begin(), new_peak);
      new_peak.setMZ(filter_spec.front().getMZ() - dist_average);
      filter_spec.insert(filter_spec.begin(), new_peak);
      new_peak.setMZ(filter_spec.front().getMZ() - dist_average);
      filter_spec.insert(filter_spec.begin(), new_peak);

      // To get an estimate of the peak quality, we probably should not smooth
      // and/or transform the data.
      if (smooth_data)
      {
        GaussFilter filter;
        Param filter_param(filter.getParameters());
        filter.setParameters(filter_param);
        filter_param.setValue("gaussian_width", 4 * dist_average);
        filter.setParameters(filter_param);
        filter.filter(filter_spec);
      }

      // transform the data for fitting and fit RT profile
      for (Size j = 0; j != filter_spec.size(); ++j)
      {
        LocalPeakType p;
        p.setPosition(filter_spec[j].getMZ());
        p.setIntensity(filter_spec[j].getIntensity());
        data_to_fit.push_back(p);
      }
    }

    EmgFitter1D fitter_emg1D_;
  };
}

#endif /* EMGSCORING_H_ */
