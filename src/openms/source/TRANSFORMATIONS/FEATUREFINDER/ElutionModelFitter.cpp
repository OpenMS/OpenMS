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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ElutionModelFitter.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLinear.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EGHTraceFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussTraceFitter.h>

using namespace OpenMS;
using namespace std;


ElutionModelFitter::ElutionModelFitter():
  DefaultParamHandler("ElutionModelFitter")
{
  vector<String> truefalse = ListUtils::create<String>("true,false");
  vector<String> advanced(1, "advanced");

  defaults_.setValue("asymmetric", "false", "Fit an asymmetric (exponential-Gaussian hybrid) model? By default a symmetric (Gaussian) model is used.");
  defaults_.setValidStrings("asymmetric", truefalse);

  defaults_.setValue("add_zeros", 0.2, "Add zero-intensity points outside the feature range to constrain the model fit. This parameter sets the weight given to these points during model fitting; '0' to disable.", advanced);
  defaults_.setMinFloat("add_zeros", 0.0);

  defaults_.setValue("unweighted_fit", "false", "Suppress weighting of mass traces according to theoretical intensities when fitting elution models", advanced);
  defaults_.setValidStrings("unweighted_fit", truefalse);

  defaults_.setValue("no_imputation", "false", "If fitting the elution model fails for a feature, set its intensity to zero instead of imputing a value from the initial intensity estimate", advanced);
  defaults_.setValidStrings("no_imputation", truefalse);

  defaults_.setValue("check:min_area", 1.0, "Lower bound for the area under the curve of a valid elution model", advanced);
  defaults_.setMinFloat("check:min_area", 0.0);

  defaults_.setValue("check:boundaries", 0.5, "Time points corresponding to this fraction of the elution model height have to be within the data region used for model fitting", advanced);
  defaults_.setMinFloat("check:boundaries", 0.0);
  defaults_.setMaxFloat("check:boundaries", 1.0);

  defaults_.setValue("check:width", 10.0, "Upper limit for acceptable widths of elution models (Gaussian or EGH), expressed in terms of modified (median-based) z-scores; '0' to disable", advanced);
  defaults_.setMinFloat("check:width", 0.0);

  defaults_.setValue("check:asymmetry", 10.0, "Upper limit for acceptable asymmetry of elution models (EGH only), expressed in terms of modified (median-based) z-scores; '0' to disable", advanced);
  defaults_.setMinFloat("check:asymmetry", 0.0);

  defaults_.setSectionDescription("check", "Parameters for checking the validity of elution models (and rejecting them if necessary)");

  defaultsToParam_();
}


ElutionModelFitter::~ElutionModelFitter() {}


double ElutionModelFitter::calculateFitQuality_(const TraceFitter* fitter,
                                                const MassTraces& traces)
{
  double mre = 0.0;
  double total_weights = 0.0;
  double rt_start = max(fitter->getLowerRTBound(), traces[0].peaks[0].first);
  double rt_end = min(fitter->getUpperRTBound(), 
                      traces[0].peaks.back().first);

  for (MassTraces::const_iterator tr_it = traces.begin();
       tr_it != traces.end(); ++tr_it)
  {
    for (vector<pair<double, const Peak1D*> >::const_iterator p_it = 
           tr_it->peaks.begin(); p_it != tr_it->peaks.end(); ++p_it)
    {
      double rt = p_it->first;
      if ((rt >= rt_start) && (rt <= rt_end))
      {
        double model_value = fitter->getValue(rt);
        double diff = fabs(model_value * tr_it->theoretical_int -
                           p_it->second->getIntensity());
        mre += diff / model_value;
        total_weights += tr_it->theoretical_int;
      }
    }
  }
  return mre / total_weights;
}


void ElutionModelFitter::fitElutionModels(FeatureMap& features)
{
  bool asymmetric = param_.getValue("asymmetric").toBool();
  double add_zeros = param_.getValue("add_zeros");
  bool weighted = !param_.getValue("unweighted_fit").toBool();
  bool impute = !param_.getValue("no_imputation").toBool();
  double check_boundaries = param_.getValue("check:boundaries");
  double area_limit = param_.getValue("check:min_area");
  double width_limit = param_.getValue("check:width");
  double asym_limit = (asymmetric ? 
                       double(param_.getValue("check:asymmetry")) : 0.0);

  TraceFitter* fitter;
  if (asymmetric)
  {
    fitter = new EGHTraceFitter();
  }
  else fitter = new GaussTraceFitter();
  if (weighted)
  {
    Param params = fitter->getDefaults();
    params.setValue("weighted", "true");
    fitter->setParameters(params);
  }

  // store model parameters to find outliers later; store values redundantly -
  // once aligned with the features in the map, once only for successful models:
  vector<double> widths_all, widths_good, asym_all, asym_good;
  if (width_limit > 0)
  {
    widths_all.resize(features.size(), numeric_limits<double>::quiet_NaN());
    widths_good.reserve(features.size());
  }
  if (asym_limit > 0)
  {
    asym_all.resize(features.size(), numeric_limits<double>::quiet_NaN());
    asym_good.reserve(features.size());
  }

  // collect peaks that constitute mass traces:
  LOG_DEBUG << "Fitting elution models to features:" << endl;
  Size index = 0;
  for (FeatureMap::Iterator feat_it = features.begin(); 
       feat_it != features.end(); ++feat_it, ++index)
  {
    // LOG_DEBUG << String(feat_it->getMetaValue("PeptideRef")) << endl;
    double region_start = double(feat_it->getMetaValue("leftWidth"));
    double region_end = double(feat_it->getMetaValue("rightWidth"));

    if (feat_it->getSubordinates().empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No subordinate features for mass traces available.");
    }
    const Feature& sub = feat_it->getSubordinates()[0];
    if (sub.getConvexHulls().empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No hull points for mass trace in subordinate feature available.");
    }

    vector<Peak1D> peaks;
    // reserve space once, to avoid copying and invalidating pointers:
    Size points_per_hull = sub.getConvexHulls()[0].getHullPoints().size();
    peaks.reserve(feat_it->getSubordinates().size() * points_per_hull +
                  (add_zeros > 0.0)); // don't forget additional zero point
    MassTraces traces;
    traces.max_trace = 0;
    // need a mass trace for every transition, plus maybe one for add. zeros:
    traces.reserve(feat_it->getSubordinates().size() + (add_zeros > 0.0));
    for (vector<Feature>::iterator sub_it = feat_it->getSubordinates().begin();
         sub_it != feat_it->getSubordinates().end(); ++sub_it)
    {
      MassTrace trace;
      trace.peaks.reserve(points_per_hull);
      trace.theoretical_int = sub_it->getMetaValue("isotope_probability");
      const ConvexHull2D& hull = sub_it->getConvexHulls()[0];
      for (ConvexHull2D::PointArrayTypeConstIterator point_it = 
             hull.getHullPoints().begin(); point_it !=
             hull.getHullPoints().end(); ++point_it)
      {
        double intensity = point_it->getY();
        if (intensity > 0.0) // only use non-zero intensities for fitting
        {
          Peak1D peak;
          peak.setMZ(sub_it->getMZ());
          peak.setIntensity(intensity);
          peaks.push_back(peak);
          trace.peaks.push_back(make_pair(point_it->getX(), &peaks.back()));
        }
      }
      trace.updateMaximum();
      if (!trace.peaks.empty()) traces.push_back(trace);
    }

    // find the trace with maximal intensity:
    Size max_trace = 0;
    double max_intensity = 0;
    for (Size i = 0; i < traces.size(); ++i)
    {
      if (traces[i].max_peak->getIntensity() > max_intensity)
      {
        max_trace = i;
        max_intensity = traces[i].max_peak->getIntensity();
      }
    }
    traces.max_trace = max_trace;
    traces.baseline = 0.0;

    if (add_zeros > 0.0)
    {
      MassTrace trace;
      trace.peaks.reserve(2);
      trace.theoretical_int = add_zeros;
      Peak1D peak;
      peak.setMZ(feat_it->getSubordinates()[0].getMZ());
      peak.setIntensity(0.0);
      peaks.push_back(peak);
      double offset = 0.2 * (region_start - region_end);
      trace.peaks.push_back(make_pair(region_start - offset, &peaks.back()));
      trace.peaks.push_back(make_pair(region_end + offset, &peaks.back()));
      traces.push_back(trace);
    }

    // fit the model:
    bool fit_success = true;
    try
    {
      fitter->fit(traces);
    }
    catch (Exception::UnableToFit& except)
    {
      LOG_ERROR << "Error fitting model to feature '" << feat_it->getUniqueId()
                << "': " << except.getName() << " - " << except.getMessage()
                << endl;
      fit_success = false;
    }

    // record model parameters:
    double center = fitter->getCenter(), height = fitter->getHeight();
    feat_it->setMetaValue("model_height", height);
    feat_it->setMetaValue("model_FWHM", fitter->getFWHM());
    feat_it->setMetaValue("model_center", center);
    feat_it->setMetaValue("model_lower", fitter->getLowerRTBound());
    feat_it->setMetaValue("model_upper", fitter->getUpperRTBound());
    if (asymmetric)
    {
      EGHTraceFitter* egh = static_cast<EGHTraceFitter*>(fitter);
      feat_it->setMetaValue("model_EGH_tau", egh->getTau());
      feat_it->setMetaValue("model_EGH_sigma", egh->getSigma());
    }
    else
    {
      GaussTraceFitter* gauss = static_cast<GaussTraceFitter*>(fitter);
      feat_it->setMetaValue("model_Gauss_sigma", gauss->getSigma());
    }

    // goodness of fit:
    double mre = -1.0; // mean relative error
    if (fit_success)
    {
      mre = calculateFitQuality_(fitter, traces);
    }
    feat_it->setMetaValue("model_error", mre);

    // check model validity:
    double area = fitter->getArea();
    feat_it->setMetaValue("model_area", area);
    if ((area != area) || (area <= area_limit)) // x != x: test for NaN
    {
      feat_it->setMetaValue("model_status", "1 (invalid area)");
    }
    else if ((center <= region_start) || (center >= region_end))
    {
      feat_it->setMetaValue("model_status", "2 (center out of bounds)");
    }
    else if (fitter->getValue(region_start) > check_boundaries * height)
    {
      feat_it->setMetaValue("model_status", "3 (left side out of bounds)");
    }
    else if (fitter->getValue(region_end) > check_boundaries * height)
    {
      feat_it->setMetaValue("model_status", "4 (right side out of bounds)");
    }
    else
    {
      feat_it->setMetaValue("model_status", "0 (valid)");
      // store model parameters to find outliers later:
      if (asymmetric)
      {
        double sigma = feat_it->getMetaValue("model_EGH_sigma");
        double abs_tau = fabs(double(feat_it->getMetaValue("model_EGH_tau")));
        if (width_limit > 0)
        {
          // see implementation of "EGHTraceFitter::getArea":
          double width = sigma * 0.6266571 + abs_tau;
          widths_all[index] = width;
          widths_good.push_back(width);
        }
        if (asym_limit > 0)
        {
          double asymmetry = abs_tau / sigma;
          asym_all[index] = asymmetry;
          asym_good.push_back(asymmetry);
        }
      }
      else if (width_limit > 0)
      {
        double width = feat_it->getMetaValue("model_Gauss_sigma");
        widths_all[index] = width;
        widths_good.push_back(width);
      }
    }
  }
  delete fitter;

  // find outliers in model parameters:
  if (width_limit > 0)
  {
    double median_width = Math::median(widths_good.begin(), widths_good.end());
    vector<double> abs_diffs(widths_good.size());
    for (Size i = 0; i < widths_good.size(); ++i)
    {
      abs_diffs[i] = fabs(widths_good[i] - median_width);
    }
    // median absolute deviation (constant factor to approximate std. dev.):
    double mad_width = 1.4826 * Math::median(abs_diffs.begin(), 
                                             abs_diffs.end());

    for (Size i = 0; i < features.size(); ++i)
    {
      double width = widths_all[i];
      if (width != width) continue; // NaN (failed model)
      double z_width = (width - median_width) / mad_width; // mod. z-score
      if (z_width > width_limit)
      {
        features[i].setMetaValue("model_status", "5 (width too large)");
        if (asym_limit > 0) // skip asymmetry check below
        {
          asym_all[i] = numeric_limits<double>::quiet_NaN();
        }
      }
      else if (z_width < -width_limit)
      {
        features[i].setMetaValue("model_status", "6 (width too small)");
        if (asym_limit > 0) // skip asymmetry check below
        {
          asym_all[i] = numeric_limits<double>::quiet_NaN();
        }
      }
    }
  }
  if (asym_limit > 0)
  {
    double median_asym = Math::median(asym_good.begin(), asym_good.end());
    vector<double> abs_diffs(asym_good.size());
    for (Size i = 0; i < asym_good.size(); ++i)
    {
      abs_diffs[i] = fabs(asym_good[i] - median_asym);
    }
    // median absolute deviation (constant factor to approximate std. dev.):
    double mad_asym = 1.4826 * Math::median(abs_diffs.begin(),
                                            abs_diffs.end());

    for (Size i = 0; i < features.size(); ++i)
    {
      double asym = asym_all[i];
      if (asym != asym) continue; // NaN (failed model)
      double z_asym = (asym - median_asym) / mad_asym; // mod. z-score
      if (z_asym > asym_limit)
      {
        features[i].setMetaValue("model_status", "7 (asymmetry too high)");
      }
      else if (z_asym < -asym_limit) // probably shouldn't happen in practice
      {
        features[i].setMetaValue("model_status", "8 (asymmetry too low)");
      }
    }
  }

  // impute approximate results for failed model fits (basically bring the
  // OpenSWATH intensity estimates to the same scale as the model-based ones):
  TransformationModel::DataPoints quant_values;
  vector<FeatureMap::Iterator> failed_models;
  Size model_successes = 0, model_failures = 0;

  for (FeatureMap::Iterator feat_it = features.begin(); 
       feat_it != features.end(); ++feat_it, ++index)
  {
    feat_it->setMetaValue("raw_intensity", feat_it->getIntensity());
    if (String(feat_it->getMetaValue("model_status"))[0] != '0')
    {
      if (impute) failed_models.push_back(feat_it);
      else feat_it->setIntensity(0.0);
      model_failures++;
    }
    else
    {
      double area = feat_it->getMetaValue("model_area");
      if (impute)
      { // apply log-transform to weight down high outliers:
        double raw_intensity = feat_it->getIntensity();
        LOG_DEBUG << "Successful model: x = " << raw_intensity << ", y = "
                  << area << "; log(x) = " << log(raw_intensity) 
                  << ", log(y) = " << log(area) << endl;
        quant_values.push_back(make_pair(log(raw_intensity), log(area)));
      }
      feat_it->setIntensity(area);
      model_successes++;
    }
  }
  LOG_INFO << "Model fitting: " << model_successes << " successes, "
           << model_failures << " failures" << endl;

  if (impute) // impute results for cases where the model fit failed
  {
    TransformationModelLinear lm(quant_values, Param());
    double slope, intercept;
    String x_weight, y_weight;
    double x_datum_min, x_datum_max, y_datum_min, y_datum_max;
    lm.getParameters(slope, intercept, x_weight, y_weight, x_datum_min, x_datum_max, y_datum_min, y_datum_max);
    LOG_DEBUG << "LM slope: " << slope << ", intercept: " << intercept << endl;
    for (vector<FeatureMap::Iterator>::iterator it = failed_models.begin();
         it != failed_models.end(); ++it)
    {
      double area = exp(lm.evaluate(log((*it)->getIntensity())));
      (*it)->setIntensity(area);
    }
  }
}
