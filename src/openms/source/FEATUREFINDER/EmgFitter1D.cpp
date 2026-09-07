// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/EmgFitter1D.h>
#include <OpenMS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/FEATUREFINDER/EmgModel.h>

#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <Eigen/Core>

#include <algorithm>
#include <cmath>

namespace OpenMS
{
  const EmgFitter1D::CoordinateType EmgFitter1D::EgmFitterFunctor::sqrt2pi = sqrt(2.0 * Constants::PI);
  const EmgFitter1D::CoordinateType EmgFitter1D::EgmFitterFunctor::emg_const = 2.4055;
  const EmgFitter1D::CoordinateType EmgFitter1D::EgmFitterFunctor::sqrt_2 = sqrt(2.0);
  const EmgFitter1D::CoordinateType EmgFitter1D::EgmFitterFunctor::c = -emg_const / sqrt_2;

  int EmgFitter1D::EgmFitterFunctor::operator()(const double* x, double* fvec) const
  {
    // Create Eigen::Map views for convenient indexing
    Eigen::Map<const Eigen::VectorXd> x_map(x, m_inputs);
    Eigen::Map<Eigen::VectorXd> fvec_map(fvec, m_values);

    Size n = m_data->n;
    const EmgFitter1D::RawDataArrayType& set = m_data->set;

    EmgFitter1D::CoordinateType h = x_map(0);
    EmgFitter1D::CoordinateType w = x_map(1);
    EmgFitter1D::CoordinateType s = x_map(2);
    EmgFitter1D::CoordinateType z = x_map(3);

    EmgFitter1D::CoordinateType Yi = 0.0;
    double prefix = (h * w / s) * sqrt2pi;
    double part1 = w*w / (2 * s*s);
    double part2 = w / s;

    // iterate over all points of the signal
    for (Size i = 0; i < n; i++)
    {
      double t = set[i].getPos();
      double diff = t - z;

      // Simplified EMG (doi=10.1.1.915.3568) Equation 9
      Yi = prefix * exp(part1 - (diff / s)) / (1 + exp(c * ((diff / w) - part2)));

      fvec_map(i) = Yi - set[i].getIntensity();
    }
    return 0;
  }

  // compute Jacobian matrix for the different parameters
  int EmgFitter1D::EgmFitterFunctor::df(const double* x, double* J) const
  {
    // Create Eigen::Map views for convenient indexing
    Eigen::Map<const Eigen::VectorXd> x_map(x, m_inputs);
    Eigen::Map<Eigen::MatrixXd> J_map(J, m_values, m_inputs);

    Size n =  m_data->n;
    const EmgFitter1D::RawDataArrayType& set = m_data->set;

    EmgFitter1D::CoordinateType h = x_map(0);
    EmgFitter1D::CoordinateType w = x_map(1);
    EmgFitter1D::CoordinateType w2 = w*w;
    EmgFitter1D::CoordinateType s = x_map(2);
    EmgFitter1D::CoordinateType s2 = s*s;
    EmgFitter1D::CoordinateType s3 = s2 * s;
    EmgFitter1D::CoordinateType z = x_map(3);

    EmgFitter1D::CoordinateType diff, exp1, exp2, exp3 = 0.0;
    EmgFitter1D::CoordinateType derivative_height, derivative_width, derivative_symmetry, derivative_retention = 0.0;


    // iterate over all points of the signal
    for (Size i = 0; i < n; i++)
    {
      EmgFitter1D::CoordinateType t = set[i].getPos();
      diff = t - z;
      exp1 = exp((w2 / (2 * s2)) - (diff / s));
      exp3 = exp((-emg_const / sqrt_2) * ((diff / w) - w / s));
      exp2 = 1 + exp3;

      // f'(h)
      derivative_height = w / s * sqrt2pi * exp1 / exp2;

      // f'(w)
      derivative_width = h / s * sqrt2pi * exp1 / exp2 + (h * w2) / s3 * sqrt2pi * exp1 / exp2 + (emg_const * h * w) / s * sqrt2pi * exp1 * (-diff / w2 - 1 / s) * exp3 / ((exp2 * exp2) * sqrt_2);

      // f'(s)
      derivative_symmetry = -h * w / s2 * sqrt2pi * exp1 / exp2 + h * w / s * sqrt2pi * (-(w * w) / s3 + diff / s2) * exp1 / exp2 + (emg_const * h * w2) / s3 * sqrt2pi * exp1 * exp3 / ((exp2 * exp2) * sqrt_2);

      // f'(z)
      derivative_retention = h * w / s2 * sqrt2pi * exp1 / exp2 - (emg_const * h) / s * sqrt2pi * exp1 * exp3 / ((exp2 * exp2) * sqrt_2);

      // set the jacobian matrix
      J_map(i, 0) = derivative_height;
      J_map(i, 1) = derivative_width;
      J_map(i, 2) = derivative_symmetry;
      J_map(i, 3) = derivative_retention;
    }
    return 0;
  }

  EmgFitter1D::EmgFitter1D() :
    LevMarqFitter1D()
  {
    setName("EmgFitter1D");
    defaults_.setValue("init_mom", "false", "Initialize parameters using method of moments estimators.", {"advanced"});
    defaults_.setValidStrings("init_mom", {"true","false"});
    defaults_.setValue("statistics:variance", 1.0, "Variance of the model.", {"advanced"});
    defaultsToParam_();
  }

  EmgFitter1D::EmgFitter1D(const EmgFitter1D& source) :
    LevMarqFitter1D(source)
  {
    setParameters(source.getParameters());
    updateMembers_();
  }

  EmgFitter1D::~EmgFitter1D() = default;

  EmgFitter1D& EmgFitter1D::operator=(const EmgFitter1D& source)
  {
    if (&source == this)
    {
      return *this;
    }
    LevMarqFitter1D::operator=(source);
    setParameters(source.getParameters());
    updateMembers_();

    return *this;
  }

  EmgFitter1D::QualityType EmgFitter1D::fit1d(const RawDataArrayType& set, std::unique_ptr<InterpolationModel>& model)
  {
    // Robustness (issue #6239): this method must never throw, never read out of bounds,
    // and never emit non-finite model parameters. A degenerate or failed fit yields a
    // stable, in-range fallback model and the quality sentinel -1.0. Well-formed fits are
    // handled exactly as before (no change to scores feeding OpenSWATH/FeatureFinderId.).

    // Builds an EmgModel from a parameter set; shared by the regular and fallback paths so
    // they can never diverge on parameter keys / interpolation step.
    auto build_model = [this](double height, double width, double symmetry, double retention,
                              double bb_min, double bb_max) -> std::unique_ptr<InterpolationModel>
    {
      std::unique_ptr<InterpolationModel> m(new EmgModel());
      m->setInterpolationStep((std::isfinite(interpolation_step_) && interpolation_step_ > 0.0) ? interpolation_step_ : 0.2);
      Param p;
      p.setValue("bounding_box:min", bb_min);
      p.setValue("bounding_box:max", bb_max);
      p.setValue("statistics:variance", statistics_.variance());
      p.setValue("statistics:mean", statistics_.mean());
      p.setValue("emg:height", height);
      p.setValue("emg:width", width);
      p.setValue("emg:symmetry", symmetry);
      p.setValue("emg:retention", retention);
      m->setParameters(p);
      return m;
    };

    const Size n = set.size();

    // Degenerate: fewer points than the 4 EMG parameters (LM requires N >= params).
    if (n < 4)
    {
      double lo = 0.0, hi = 1.0, apex = 0.0, height = 1.0;
      if (n > 0)
      {
        lo = hi = apex = set[0].getPos();
        height = set[0].getIntensity();
        for (Size i = 1; i < n; ++i)
        {
          const double pos = set[i].getPos();
          if (std::isfinite(pos)) { lo = std::min(lo, pos); hi = std::max(hi, pos); }
          if (set[i].getIntensity() > height) { height = set[i].getIntensity(); apex = pos; }
        }
      }
      if (!std::isfinite(lo) || !std::isfinite(hi) || hi <= lo) { lo = 0.0; hi = 1.0; }
      if (!std::isfinite(apex) || apex < lo || apex > hi) { apex = lo; }
      const double w = std::max((hi - lo) / 6.0, 1e-3);
      model = build_model(std::isfinite(height) ? height : 1.0, w, w, apex, lo, hi);
      return -1.0;
    }

    // Single pass: data range, apex (max-intensity position), intensity sum, finiteness.
    double min_pos = set[0].getPos(), max_pos = set[0].getPos();
    double apex_pos = set[0].getPos(), max_int = set[0].getIntensity();
    double weight_sum = 0.0;
    bool finite_input = std::isfinite(set[0].getPos()) && std::isfinite(set[0].getIntensity());
    for (Size i = 0; i < n; ++i)
    {
      const double pos = set[i].getPos();
      const double in = set[i].getIntensity();
      if (!std::isfinite(pos) || !std::isfinite(in)) { finite_input = false; }
      if (pos < min_pos) { min_pos = pos; }
      if (pos > max_pos) { max_pos = pos; }
      if (in > max_int) { max_int = in; apex_pos = pos; }
      weight_sum += in;
    }
    // Sanitize the fallback range so a fallback model is ALWAYS finite and in-range, even
    // when the input itself is non-finite (e.g. a NaN/Inf leading position leaves min_pos/
    // max_pos non-finite, since NaN comparisons never update them). For well-formed data
    // (finite positions, positive span) every condition below is false, so the regular path
    // is unchanged. (issue #6239 contract: never emit non-finite model parameters.)
    if (!std::isfinite(min_pos) || !std::isfinite(max_pos) || !(max_pos > min_pos))
    {
      const double center = std::isfinite(apex_pos) ? apex_pos : 0.0;
      min_pos = center - 0.5;
      max_pos = center + 0.5;
    }
    const double span = max_pos - min_pos;
    const double fb_w = std::max(span / 6.0, 1e-3);
    if (!std::isfinite(apex_pos) || apex_pos < min_pos || apex_pos > max_pos) { apex_pos = min_pos; }
    if (!std::isfinite(max_int) || !(max_int > 0.0)) { max_int = 1.0; }

    // Unusable input or fitter settings -> stable fallback + sentinel.
    if (!finite_input || !(span > 0.0) || !std::isfinite(weight_sum) || !(weight_sum > 0.0)
        || !std::isfinite(interpolation_step_) || !(interpolation_step_ > 0.0))
    {
      model = build_model(max_int, fb_w, fb_w, apex_pos, min_pos, max_pos);
      return -1.0;
    }

    // Bounding box enlarged by a few multiples of the standard deviation (unchanged).
    const CoordinateType stdev = sqrt(statistics_.variance()) * tolerance_stdev_box_;
    const CoordinateType min_bb = min_pos - stdev;
    const CoordinateType max_bb = max_pos + stdev;

    // Set advanced parameters for residual_  und jacobian_ method
    EmgFitter1D::Data d;
    d.n = set.size();
    d.set = set;

    // Compute start parameters (method-of-moments path is OOB-safe, see setInitialParametersMOM_)
    setInitialParameters_(set);

    // Optimize parameter with Levenberg-Marquardt algorithm
    std::vector<double> x_init(4);
    x_init[0] = height_;
    x_init[1] = width_;
    x_init[2] = symmetry_;
    x_init[3] = retention_;

    bool fit_ok = true;
    if (!symmetric_)
    {
      try
      {
        EgmFitterFunctor functor(4, &d);
        optimize_(x_init, functor);
      }
      catch (const Exception::UnableToFit&)
      {
        fit_ok = false;
      }
    }

    // Set optimized parameters
    height_ = x_init[0];
    width_ = x_init[1];
    symmetry_ = x_init[2];
    retention_ = x_init[3];

    // Reject non-finite optimizer output: emit a stable fallback model instead of garbage.
    if (!fit_ok || !std::isfinite(height_) || !std::isfinite(width_)
        || !std::isfinite(symmetry_) || !std::isfinite(retention_))
    {
      model = build_model(max_int, fb_w, fb_w, apex_pos, min_pos, max_pos);
      return -1.0;
    }

    // build model
    model = build_model(height_, width_, symmetry_, retention_, min_bb, max_bb);

    // calculate pearson correlation
    std::vector<float> real_data;
    real_data.reserve(set.size());
    std::vector<float> model_data;
    model_data.reserve(set.size());

    for (Size i = 0; i < set.size(); ++i)
    {
      real_data.push_back(set[i].getIntensity());
      model_data.push_back(model->getIntensity(DPosition<1>(set[i].getPosition())));
    }

    QualityType correlation = Math::pearsonCorrelationCoefficient(real_data.begin(), real_data.end(), model_data.begin(), model_data.end());

    // Non-finite correlation -> unusable fit: stable fallback model + sentinel.
    if (!std::isfinite(correlation))
    {
      model = build_model(max_int, fb_w, fb_w, apex_pos, min_pos, max_pos);
      return -1.0;
    }
    return correlation;
  }

  void EmgFitter1D::setInitialParametersMOM_(const RawDataArrayType& set)
  {
    std::vector<CoordinateType> weighted;
    weighted.reserve(set.size());
    CoordinateType weighted_sum = 0.;
    CoordinateType weight_sum = 0.;
    for (Size s = 0 ; s < set.size() ; ++s)
    {
      weighted_sum += set[s].getPos() * set[s].getIntensity();
      weight_sum += set[s].getIntensity();
    }
    CoordinateType weighted_mean = weighted_sum / weight_sum;

    int weighted_median_idx = 0;
    double sum = weight_sum - set[0].getIntensity(); // sum is the total weight of all `x[i] > x[k]`

    // Bound the index (issue #6239): with non-positive/degenerate weights `sum` may not
    // decrease, so guard against advancing past the last element (out-of-bounds read).
    while (sum > weight_sum/2. && weighted_median_idx + 1 < static_cast<int>(set.size()))
    {
      ++weighted_median_idx;
      sum -= set[weighted_median_idx].getIntensity();
    }
    CoordinateType weighted_median = set[weighted_median_idx].getPos();

    CoordinateType weighted_sd = 0.;
    for (Size s = 0 ; s < set.size() ; ++s)
    {
      weighted_sd += std::pow(weighted_mean - set[s].getPos(), 2) * set[s].getIntensity();
    }
    weighted_sd /= weight_sum;
    weighted_sd = std::sqrt(weighted_sd);
    CoordinateType weighted_skew = std::fabs(weighted_mean - weighted_median) / weighted_sd;

    CoordinateType max_peak_width = fabs(set[set.size() - 1].getPos() - set[weighted_median_idx].getPos()); // cannot be wider than this

    // calculate the height of the peak
    height_ = set[weighted_median_idx].getIntensity();

    // calculate retention time
    retention_ = weighted_mean - weighted_sd * std::pow(weighted_skew / 2., 1./3.);

    // default is an asymmetric peak
    symmetric_ = false;

    // calculate the symmetry (fronted peak: s<1 , tailed peak: s>1)
    symmetry_ = weighted_sd * std::pow(weighted_skew / 2., 1./3.);

    // check the symmetry
    if (std::isinf(symmetry_) || std::isnan(symmetry_))
    {
      symmetric_ = true;
      symmetry_ = 10.0;
    }

    // optimize the symmetry
    // The computations can lead to an overflow error at very low values of symmetry (s~0).
    // For s~5 the parameter can be approximated by the Levenberg-Marquardt algorithms.
    // (the other parameters are much greater than one)
    if (symmetry_ < 1)
    {
      symmetry_ += 5;
    }

    // Need to ensure that we do not go beyond the maximal width of the peak
    symmetry_ = std::min(symmetry_, max_peak_width);

    // calculate the width of the peak
    // rt-values with intensity zero are not allowed for calculation of the width
    // normally: width_ = fabs( set[set.size() - 1].getPos() - set[0].getPos() );
    // but its better for the emg function to proceed from narrow peaks
    width_ = symmetry_;
    //MOM estimator would be the following, but it is too large for the test
    //width_ = weighted_sd * std::sqrt(std::pow(1. - (weighted_skew / 2.), 2./3.));
  }

  void EmgFitter1D::setInitialParameters_(const RawDataArrayType& set)
  {
    if (param_.getValue("init_mom").toBool())
    {
      setInitialParametersMOM_(set);
      return;
    }

    // sum over all intensities
    CoordinateType sum = 0.0;
    for (Size i = 0; i < set.size(); ++i)
    {
      sum += set[i].getIntensity();
    }
    // calculate the median
    Size median = 0;
    float count = 0.0;
    for (Size i = 0; i < set.size(); ++i)
    {
      count += set[i].getIntensity();
      if (count <= sum / 2)
      {
        median = i;
      }
    }

    double max_peak_width = fabs(set[set.size() - 1].getPos() - set[median].getPos()); // cannot be wider than this

    // calculate the height of the peak
    height_ = set[median].getIntensity();

    // calculate retention time
    retention_ = set[median].getPos();

    // default is an asymmetric peak
    symmetric_ = false;

    // calculate the symmetry (fronted peak: s<1 , tailed peak: s>1)
    symmetry_ = fabs(set[set.size() - 1].getPos() - set[median].getPos()) / fabs(set[median].getPos() - set[0].getPos());

    // check the symmetry
    if (std::isinf(symmetry_) || std::isnan(symmetry_))
    {
      symmetric_ = true;
      symmetry_ = 10.0;
    }

    // optimize the symmetry
    // The computations can lead to an overflow error at very low values of symmetry (s~0).
    // For s~5 the parameter can be approximated by the Levenberg-Marquardt algorithms.
    // (the other parameters are much greater than one)
    if (symmetry_ < 1)
    {
      symmetry_ += 5;
    }

    // Need to ensure that we do not go beyond the maximal width of the peak
    symmetry_ = std::min(symmetry_, max_peak_width);

    // calculate the width of the peak
    // rt-values with intensity zero are not allowed for calculation of the width
    // normally: width_ = fabs( set[set.size() - 1].getPos() - set[0].getPos() );
    // but its better for the emg function to proceed from narrow peaks
    width_ = symmetry_;
  }

  void EmgFitter1D::updateMembers_()
  {
    LevMarqFitter1D::updateMembers_();
    statistics_.setVariance(param_.getValue("statistics:variance"));
  }

}
