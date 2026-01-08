// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FEATUREFINDER/EmgFitter1D.h>
#include <OpenMS/FEATUREFINDER/EmgModel.h>
#include <OpenMS/PROCESSING/SMOOTHING/GaussFilter.h>

#include <OpenMS/KERNEL/MRMFeature.h>
#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <OpenMS/KERNEL/StandardTypes.h>

#include <vector>
#include <cmath> // for isnan

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

    EmgScoring() = default;

    ~EmgScoring() = default;

    /// overwrites params for the Emg1DFitter. Unspecified params will stay default.
    /// use getDefaults to see what you can set.
    void setFitterParam(const Param& param)
    {
      fitter_emg1D_params_ = param;
    }

    /// Get default params for the Emg1D fitting
    Param getDefaults()
    {
      return EmgFitter1D().getDefaults();
    }

    /// calculate the elution profile fit score
    template<typename SpectrumType, class TransitionT>
    double calcElutionFitScore(MRMFeature & mrmfeature, MRMTransitionGroup<SpectrumType, TransitionT> & transition_group) const
    {
      double avg_score = 0;
      bool smooth_data = false;

      for (Size k = 0; k < transition_group.size(); k++)
      {
        // get the id, then find the corresponding transition and features within this peakgroup
        String native_id = transition_group.getChromatograms()[k].getNativeID();
        Feature f = mrmfeature.getFeature(native_id);
        OPENMS_PRECONDITION(f.getConvexHulls().size() == 1, "Convex hulls need to have exactly one hull point structure");

        //TODO think about penalizing aborted fits even more. Currently -1 is just the "lowest" pearson correlation to
        // a fit that you can have.
        double fscore = elutionModelFit(f.getConvexHulls()[0].getHullPoints(), smooth_data);
        avg_score += fscore;
      }

      avg_score /= transition_group.size();
      return avg_score;
    }

    // Fxn from FeatureFinderAlgorithmMRM
    // TODO: check whether we can leave out some of the steps here, e.g. gaussian smoothing
    double elutionModelFit(const ConvexHull2D::PointArrayType& current_section, bool smooth_data) const
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
      std::unique_ptr<InterpolationModel> model_rt;
      double quality = fitRT_(data_to_fit, model_rt);
      // cut line 354 of FeatureFinderAlgorithmMRM

      return quality;
    }

  protected:
    template<class LocalPeakType>
    double fitRT_(std::vector<LocalPeakType>& rt_input_data, std::unique_ptr<InterpolationModel>& model) const
    {
      EmgFitter1D fitter_emg1D;
      fitter_emg1D.setParameters(fitter_emg1D_params_);
      // Construct model for rt
      // NaN is checked in fit1d: if (std::isnan(quality)) quality = -1.0;
      return fitter_emg1D.fit1d(rt_input_data, model);
    }

    // Fxn from FeatureFinderAlgorithmMRM
    // TODO: check whether we can leave out some of the steps here, e.g. gaussian smoothing
    template<class LocalPeakType>
    void prepareFit_(const ConvexHull2D::PointArrayType & current_section, std::vector<LocalPeakType> & data_to_fit, bool smooth_data) const
    {
      // typedef Peak1D LocalPeakType;
      PeakSpectrum filter_spec;
      // first smooth the data to prevent outliers from destroying the fit
      for (const auto& pa : current_section)
      {
        LocalPeakType p;
        using IntensityType = typename LocalPeakType::IntensityType;
        p.setMZ(pa.getX());
        p.setIntensity(IntensityType(pa.getY()));
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

    Param fitter_emg1D_params_;
  };

}

