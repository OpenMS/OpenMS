// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/MATH/MISC/EmgGradientDescent.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <numeric>

namespace OpenMS
{
  EmgGradientDescent::EmgGradientDescent() :
    DefaultParamHandler("EmgGradientDescent")
  {
    getDefaultParameters(defaults_);
    defaultsToParam_(); // write defaults into Param object param_
  }

  void EmgGradientDescent::getDefaultParameters(Param& defaults)
  {
    defaults.clear();

    defaults.setValue("print_debug", (UInt)0, "The level of debug information to print in the terminal. Valid values are: 0, 1, 2. Higher values mean more information.");
    defaults.setMinInt("print_debug", 0);
    defaults.setMaxInt("print_debug", 2);

    defaults.setValue("max_gd_iter", (UInt)100000, "The maximum number of iterations permitted to the gradient descent algorithm.");
    defaults.setMinInt("max_gd_iter", 0);

    defaults.setValue("compute_additional_points", "true", "Whether additional points should be added when fitting EMG peak model.");
    defaults.setValidStrings("compute_additional_points", {"true","false"});
  }

  void EmgGradientDescent::updateMembers_()
  {
    print_debug_ = (UInt)param_.getValue("print_debug");
    max_gd_iter_ = (UInt)param_.getValue("max_gd_iter");
    compute_additional_points_ = param_.getValue("compute_additional_points").toBool();
  }

  double EmgGradientDescent::compute_z(
    const double x,
    const double mu,
    const double sigma,
    const double tau
  ) const
  {
    return (1.0 / std::sqrt(2.0)) * (sigma / tau - (x - mu) / sigma);
  }

  double EmgGradientDescent::E_wrt_h(
    const std::vector<double>& xs,
    const std::vector<double>& ys,
    const double h,
    const double mu,
    const double sigma,
    const double tau
  ) const
  {
    const double u = mu;
    const double s = sigma;
    const double t = tau;
    std::vector<double> diffs(xs.size());
    for (Size i = 0; i < xs.size(); ++i)
    {
      const double x = xs[i];
      const double y = ys[i];
      const double z = compute_z(x, mu, sigma, tau);
      if (z < 0)
      {
        diffs[i] = ((s * std::exp((std::pow(s,2.0) + 2.0 * t * u - 4.0 * t * x)/(2.0 * std::pow(t,2.0))) * std::erfc((std::pow(s,2.0) + t * (u - x))/(std::sqrt(2.0) * s * t)) * (PI * h * s * std::exp((std::pow(s,2.0) + 2 * t * u)/(2 * std::pow(t,2.0))) * std::erfc((std::pow(s,2.0) + t * (u - x))/(std::sqrt(2.0) * s * t)) - std::sqrt(2.0 * PI) * t * y * std::exp(x/t)))/std::pow(t,2.0)) / static_cast<double>(xs.size());
      }
      else if (z <= 6.71e7)
      {
        diffs[i] = ((std::sqrt(2.0 * PI) * s * std::exp(1.0/2.0 * std::pow((s/t - (x - u)/s),2.0) - std::pow((x - u),2.0)/(2.0 * std::pow(s,2.0))) * std::erfc((s/t - (x - u)/s)/std::sqrt(2.0)) * ((std::sqrt(PI/2.0) * h * s * std::exp(1.0/2.0 * std::pow((s/t - (x - u)/s), 2.0) - std::pow((x - u),2.0)/(2 * std::pow(s,2.0))) * std::erfc((s/t - (x - u)/s)/std::sqrt(2.0)))/t - y))/t) / static_cast<double>(xs.size());
      }
      else
      {
        diffs[i] = ((2 * std::exp(-std::pow((x - u),2.0)/(2.0 * std::pow(s,2.0))) * ((h * std::exp(-std::pow((x - u),2.0)/(2 * std::pow(s,2.0))))/(1 - (t * (x - u))/std::pow(s,2.0)) - y))/(1 - (t * (x - u))/std::pow(s,2.0))) / static_cast<double>(xs.size());
      }
    }
    const double result = std::accumulate(diffs.begin(), diffs.end(), 0.0);
    if (print_debug_ == 2)
    {
      std::cout << std::endl << "E_wrt_h() diffs:" << std::endl;
      for (const double d : diffs)
      {
        std::cout << d << " ";
      }
      std::cout << std::endl;
      std::cout << "result=" << result << std::endl;
    }
    return result;
  }

  double EmgGradientDescent::E_wrt_mu(
    const std::vector<double>& xs,
    const std::vector<double>& ys,
    const double h,
    const double mu,
    const double sigma,
    const double tau
  ) const
  {
    const double u = mu;
    const double s = sigma;
    const double t = tau;
    std::vector<double> diffs(xs.size());
    for (Size i = 0; i < xs.size(); ++i)
    {
      const double x = xs[i];
      const double y = ys[i];
      const double z = compute_z(x, mu, sigma, tau);
      if (z < 0)
      {
        diffs[i] = (2 * ((std::sqrt(PI/2.0) * h * s * std::exp(std::pow(s,2.0)/(2.0 * std::pow(t,2.0)) - (x - u)/t) * std::erfc((s/t - (x - u)/s)/std::sqrt(2.0)))/std::pow(t,2.0) - (h * std::exp(std::pow(s,2.0)/(2.0 * std::pow(t,2.0)) - 1.0/2.0 * std::pow((s/t - (x - u)/s),2.0) - (x - u)/t))/t) * ((std::sqrt(PI/2.0) * h * s * std::exp(std::pow(s,2.0)/(2.0 * std::pow(t,2.0)) - (x - u)/t) * std::erfc((s/t - (x - u)/s)/std::sqrt(2.0)))/t - y)) / static_cast<double>(xs.size());
      }
      else if (z <= 6.71e7)
      {
        diffs[i] = (2 * ((std::sqrt(PI/2.0) * h * s * std::exp(1.0/2.0 * std::pow((s/t - (x - u)/s),2.0) - std::pow((x - u),2.0)/(2.0 * std::pow(s,2.0))) * ((x - u)/std::pow(s,2.0) + (s/t - (x - u)/s)/s) * std::erfc((s/t - (x - u)/s)/std::sqrt(2.0)))/t - (h * std::exp(-std::pow((x - u),2.0)/(2.0 * std::pow(s,2.0))))/t) * ((std::sqrt(PI/2.0) * h * s * std::exp(1.0/2.0 * std::pow((s/t - (x - u)/s),2.0) - std::pow((x - u),2.0)/(2.0 * std::pow(s,2.0))) * std::erfc((s/t - (x - u)/s)/std::sqrt(2.0)))/t - y)) / static_cast<double>(xs.size());
      }
      else
      {
        diffs[i] = (2.0 * ((h * (x - u) * std::exp(-std::pow((x - u),2.0)/(2.0 * std::pow(s,2.0))))/(std::pow(s,2.0) * (1.0 - (t * (x - u))/std::pow(s,2.0))) - (h * t * std::exp(-std::pow((x - u),2.0)/(2.0 * std::pow(s,2.0))))/(std::pow(s,2.0) * std::pow((1.0 - (t * (x - u))/std::pow(s,2.0)),2.0))) * ((h * std::exp(-std::pow((x - u),2.0)/(2.0 * std::pow(s,2.0))))/(1.0 - (t * (x - u))/std::pow(s,2.0)) - y)) / static_cast<double>(xs.size());
      }
    }
    const double result = std::accumulate(diffs.begin(), diffs.end(), 0.0);
    if (print_debug_ == 2)
    {
      std::cout << std::endl << "E_wrt_mu() diffs:" << std::endl;
      for (const double d : diffs)
      {
        std::cout << d << " ";
      }
      std::cout << std::endl;
      std::cout << "result=" << result << std::endl;
    }
    return result;
  }

  double EmgGradientDescent::E_wrt_sigma(
    const std::vector<double>& xs,
    const std::vector<double>& ys,
    const double h,
    const double mu,
    const double sigma,
    const double tau
  ) const
  {
    const double u = mu;
    const double s = sigma;
    const double t = tau;
    std::vector<double> diffs(xs.size());
    for (Size i = 0; i < xs.size(); ++i)
    {
      const double x = xs[i];
      const double y = ys[i];
      const double z = compute_z(x, mu, sigma, tau);
      if (z < 0)
      {
        diffs[i] = (2.0 * ((std::sqrt(PI/2.0) * h * std::exp(std::pow(s,2.0)/(2.0 * std::pow(t,2.0)) - (x - u)/t) * std::erfc((s/t - (x - u)/s)/std::sqrt(2.0)))/t + (std::sqrt(PI/2.0) * h * std::pow(s,2.0) * std::exp(std::pow(s,2.0)/(2.0 * std::pow(t,2.0)) - (x - u)/t) * std::erfc((s/t - (x - u)/s)/std::sqrt(2.0)))/std::pow(t,3.0) - (h * s * std::exp(std::pow(s,2.0)/(2.0 * std::pow(t,2.0)) - 1.0/2.0 * std::pow((s/t - (x - u)/s),2.0) - (x - u)/t) * ((x - u)/std::pow(s,2.0) + 1.0/t))/t) * ((std::sqrt(PI/2.0) * h * s * std::exp(std::pow(s,2.0)/(2.0 * std::pow(t,2.0)) - (x - u)/t) * std::erfc((s/t - (x - u)/s)/std::sqrt(2.0)))/t - y)) / static_cast<double>(xs.size());
      }
      else if (z <= 6.71e7)
      {
        diffs[i] = (2.0 * ((std::sqrt(PI/2.0) * h * std::exp(1.0/2.0 * std::pow((s/t - (x - u)/s),2.0) - std::pow((x - u),2.0)/(2.0 * std::pow(s,2.0))) * std::erfc((s/t - (x - u)/s)/std::sqrt(2.0)))/t + (std::sqrt(PI/2.0) * h * s * std::exp(1.0/2.0 * std::pow((s/t - (x - u)/s),2.0) - std::pow((x - u),2.0)/(2.0 * std::pow(s,2.0))) * (std::pow((x - u),2.0)/std::pow(s,3.0) + ((x - u)/std::pow(s,2.0) + 1.0/t) * (s/t - (x - u)/s)) * std::erfc((s/t - (x - u)/s)/std::sqrt(2.0)))/t - (h * s * std::exp(-std::pow((x - u),2.0)/(2.0 * std::pow(s,2.0))) * ((x - u)/std::pow(s,2.0) + 1.0/t))/t) * ((std::sqrt(PI/2.0) * h * s * std::exp(1.0/2.0 * std::pow((s/t - (x - u)/s),2.0) - std::pow((x - u),2.0)/(2.0 * std::pow(s,2.0))) * std::erfc((s/t - (x - u)/s)/std::sqrt(2.0)))/t - y)) / static_cast<double>(xs.size());
      }
      else
      {
        diffs[i] = (2.0 * ((h * std::pow((x - u),2.0) * std::exp(-std::pow((x - u),2.0)/(2.0 * std::pow(s,2.0))))/(std::pow(s,3.0) * (1.0 - (t * (x - u))/std::pow(s,2.0))) - (2.0 * h * t * (x - u) * std::exp(-std::pow((x - u),2.0)/(2 * std::pow(s,2.0))))/(std::pow(s,3.0) * std::pow((1.0 - (t * (x - u))/std::pow(s,2.0)),2.0))) * ((h * std::exp(-std::pow(x-u,2.0)/(2 * std::pow(s,2.0))))/(1 - (t * (x - u))/std::pow(s,2.0)) - y)) / static_cast<double>(xs.size());
      }
    }
    const double result = std::accumulate(diffs.begin(), diffs.end(), 0.0);
    if (print_debug_ == 2)
    {
      std::cout << std::endl << "E_wrt_sigma() diffs:" << std::endl;
      for (const double d : diffs)
      {
        std::cout << d << " ";
      }
      std::cout << std::endl;
      std::cout << "result=" << result << std::endl;
    }
    return result;
  }

  double EmgGradientDescent::E_wrt_tau(
    const std::vector<double>& xs,
    const std::vector<double>& ys,
    const double h,
    const double mu,
    const double sigma,
    const double tau
  ) const
  {
    const double u = mu;
    const double s = sigma;
    const double t = tau;
    std::vector<double> diffs(xs.size());
    const double PI = OpenMS::Constants::PI;
    for (Size i = 0; i < xs.size(); ++i)
    {
      const double x = xs[i];
      const double y = ys[i];
      const double z = compute_z(x, mu, sigma, tau);
      if (z < 0)
      {
        diffs[i] = (2 * (-(std::sqrt(PI/2.0) * h * s * std::exp(std::pow(s,2.0)/(2 * std::pow(t,2.0)) - (x - u)/t) * std::erfc((s/t - (x - u)/s)/std::sqrt(2.0)))/std::pow(t,2.0) + (std::sqrt(PI/2.0) * h * s * std::exp(std::pow(s,2.0)/(2.0 * std::pow(t,2.0)) - (x - u)/t) * ((x - u)/std::pow(t,2.0) - std::pow(s,2.0)/std::pow(t,3.0)) * std::erfc((s/t - (x - u)/s)/std::sqrt(2.0)))/t + (h * std::pow(s,2.0) * std::exp(std::pow(s,2.0)/(2.0 * std::pow(t,2.0)) - 1.0/2.0 * std::pow((s/t - (x - u)/s),2.0) - (x - u)/t))/std::pow(t,3.0)) * ((std::sqrt(PI/2.0) * h * s * std::exp(std::pow(s,2.0)/(2.0 * std::pow(t,2.0)) - (x - u)/t) * std::erfc((s/t - (x - u)/s)/std::sqrt(2.0)))/t - y)) / static_cast<double>(xs.size());
      }
      else if (z <= 6.71e7)
      {
        diffs[i] = (2 * (-(std::sqrt(PI/2.0) * h * std::pow(s,2.0) * std::exp(1.0/2.0 * std::pow((s/t - (x - u)/s),2.0) - std::pow(x-u,2.0)/(2.0 * std::pow(s,2.0))) * (s/t - (x - u)/s) * std::erfc((s/t - (x - u)/s)/std::sqrt(2.0)))/std::pow(t,3.0) - (std::sqrt(PI/2.0) * h * s * std::exp(1.0/2.0 * std::pow((s/t - (x - u)/s),2.0) - std::pow(x-u,2.0)/(2.0 * std::pow(s,2.0))) * std::erfc((s/t - (x - u)/s)/std::sqrt(2.0)))/std::pow(t,2.0) + (h * std::pow(s,2.0) * std::exp(-std::pow(x-u,2.0)/(2 * std::pow(s,2.0))))/std::pow(t,3.0)) * ((std::sqrt(PI/2.0) * h * s * std::exp(1.0/2.0 * std::pow((s/t - (x - u)/s),2.0) - std::pow(x-u,2.0)/(2 * std::pow(s,2.0))) * std::erfc((s/t - (x - u)/s)/std::sqrt(2.0)))/t - y)) / static_cast<double>(xs.size());
      }
      else
      {
        diffs[i] = ((2.0 * h * (x - u) * std::exp(-std::pow(x-u,2.0)/(2.0 * std::pow(s,2.0))) * ((h * std::exp(-std::pow(x-u,2.0)/(2.0 * std::pow(s,2.0))))/(1.0 - (t * (x - u))/std::pow(s,2.0)) - y))/(std::pow(s,2.0) * std::pow((1.0 - (t * (x - u))/std::pow(s,2.0)),2.0))) / static_cast<double>(xs.size());
      }
    }
    const double result = std::accumulate(diffs.begin(), diffs.end(), 0.0);
    if (print_debug_ == 2)
    {
      std::cout << std::endl << "E_wrt_tau() diffs:" << std::endl;
      for (const double d : diffs)
      {
        std::cout << d << " ";
      }
      std::cout << std::endl;
      std::cout << "result=" << result << std::endl;
    }
    return result;
  }

  double EmgGradientDescent::Loss_function(
    const std::vector<double>& xs,
    const std::vector<double>& ys,
    const double h,
    const double mu,
    const double sigma,
    const double tau
  ) const
  {
    std::vector<double> diffs(xs.size());
    for (Size i = 0; i < xs.size(); ++i)
    {
      diffs[i] = std::pow(emg_point(xs[i], h, mu, sigma, tau) - ys[i], 2.0) / xs.size();
    }
    const double result = std::accumulate(diffs.begin(), diffs.end(), 0.0);
    if (print_debug_ == 2)
    {
      std::cout << std::endl << "Loss_function() diffs:" << std::endl;
      for (const double d : diffs)
      {
        std::cout << d << " ";
      }
      std::cout << std::endl;
      std::cout << "result=" << result << std::endl;
    }
    return result;
  }

  void EmgGradientDescent::applyEstimatedParameters(
    const std::vector<double>& xs,
    const double h,
    const double mu,
    const double sigma,
    const double tau,
    std::vector<double>& out_xs,
    std::vector<double>& out_ys
  ) const
  {
    out_xs = xs; // Copy all positions to output
    out_ys.clear();
    for (const double x : out_xs) // For each x, estimate y
    {
      out_ys.push_back(emg_point(x, h, mu, sigma, tau));
    }

    if (!compute_additional_points_)
    {
      return;
    }
    // Compute the sampling step for the additional points
    double avg_sampling { 0.0 };
    for (Size i = 1; i < xs.size(); ++i)
    {
      avg_sampling += xs[i] - xs[i - 1];
    }
    avg_sampling /= xs.size() - 1;

    // Stop adding points if the estimated y <= `est_y_threshold`
    const double est_y_threshold { 1e-3 };

    // Stop adding points if the peak becomes too large
    std::vector<double>::const_iterator apex_pos_it = std::max_element(out_ys.cbegin(), out_ys.cend());
    const double apex_pos = out_xs[std::distance(out_ys.cbegin(), apex_pos_it)];

    // Decide on which side the eventual additional points should be added
    // The loop stops once the last added point's intensity is:
    // - lower than the intensity on the opposite boundary
    // - lower than `est_y_threshold`
    // The loop stops if the cutoff side becomes 3 times larger than the other side
    if (out_ys.front() > out_ys.back())
    {
      const double pos_boundary = apex_pos - (out_xs.back() - apex_pos) * 3;
      const double target_intensity = out_ys.back();
      while (out_ys.front() > target_intensity && out_ys.front() > est_y_threshold)
      {
        const double position = out_xs.front() - avg_sampling;
        if (position < pos_boundary)
        {
          break;
        }
        out_xs.insert(out_xs.begin(), position);
        out_ys.insert(out_ys.begin(), emg_point(position, h, mu, sigma, tau));
      }
    }
    else
    {
      const double pos_boundary = apex_pos + (apex_pos - out_xs.front()) * 3;
      const double target_intensity = out_ys.front();
      while (out_ys.back() > target_intensity && out_ys.back() > est_y_threshold)
      {
        const double position = out_xs.back() + avg_sampling;
        if (position > pos_boundary)
        {
          break;
        }
        out_xs.push_back(position);
        out_ys.push_back(emg_point(position, h, mu, sigma, tau));
      }
    }
  }

  double EmgGradientDescent::computeInitialMean(
    const std::vector<double>& xs,
    const std::vector<double>& ys
  ) const
  {
    if (xs.empty())
    {
      throw Exception::SizeUnderflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 0);
    }
    const double max_intensity = *std::max_element(ys.begin(), ys.end());
    // The intensity levels at which the mean candidates are computed
    const std::vector<double> percentages = { 0.6, 0.65, 0.7, 0.75, 0.8, 0.85 };
    Size i = 0;
    Size j = xs.size() - 1;
    // Make sure left and right positions have an initial value
    // This is to avoid situations (eg. cutoff peaks) where `max_intensity_threshold`
    // is higher than the first point on a boundary of the peak. In such a case,
    // the following nested loops would not get a chance to execute and the
    // algorithm would fail.
    // The avoidance of using the highest points of the peak apex also provides
    // robustness to spurious points or random fluctuations in detector sampling
    // from inflating the maximum peak height.
    double left_pos = xs.front();
    double right_pos = xs.back();
    std::vector<double> mean_candidates;
    for (const double height_percentage : percentages)
    {
      const double max_intensity_threshold = max_intensity * height_percentage;
      for (; i < xs.size() - 1 && ys[i] <= max_intensity_threshold; ++i)
      {
        left_pos = xs[i];
      }
      for (; j >= 1 && ys[j] <= max_intensity_threshold; --j)
      {
        right_pos = xs[j];
      }
      mean_candidates.push_back( (left_pos + right_pos) / 2.0 );
    }
    // Return the average of all middle RTs
    return std::accumulate(mean_candidates.begin(), mean_candidates.end(), 0.0) / mean_candidates.size();
  }

  void EmgGradientDescent::iRpropPlus(
    const double prev_diff_E_param,
    double& diff_E_param,
    double& param_lr,
    double& param_update,
    double& param,
    const double current_E,
    const double previous_E
  ) const
  {
    if (prev_diff_E_param * diff_E_param > 0.0)
    {
      // Using value 2000 as upper bound (iRprop+ paper recommends a value of 50)
      param_lr = std::min(param_lr * 1.2, 2000.0);
      param_update = - ( diff_E_param / std::fabs(diff_E_param) ) * param_lr;
      param += param_update;
    }
    else if (prev_diff_E_param * diff_E_param < 0.0)
    {
      param_lr = std::max(param_lr * 0.5, 0.0);
      if (current_E > previous_E)
      {
        param -= param_update;
      }
      diff_E_param = 0.0;
    }
    else
    {
      if (diff_E_param)
      {
        param_update = - ( diff_E_param / std::fabs(diff_E_param) ) * param_lr;
      }
      else
      {
        param_update = - param_lr;
      }
      param += param_update;
    }
  }

  double EmgGradientDescent::emg_point(
    const double x,
    const double h,
    const double mu,
    const double sigma,
    const double tau
  ) const
  {
    const double z = compute_z(x, mu, sigma, tau);
    const double u = mu;
    const double s = sigma;
    const double t = tau;
    if (z < 0)
    {
      return ((h*s)/t) * std::sqrt(PI/2.0) * std::exp((1.0/2.0)*(std::pow(s/t,2.0))-(x-u)/t) * std::erfc((1.0/std::sqrt(2.0)) * (s/t - (x-u)/s));
    }
    else if (z <= 6.71e7)
    {
      return h * std::exp(-(1.0/2.0) * std::pow(((x - u)/s),2.0)) * (s/t) * std::sqrt(PI/2.0) * std::exp(std::pow((1.0/std::sqrt(2.0) * (s/t - (x - u)/s)),2.0)) * std::erfc(1.0/std::sqrt(2.0) * (s/t - (x - u)/s));
    }
    else
    {
      return (h * std::exp(-(1.0/2.0) * (std::pow(((x-u) / s),2.0)))) / (1.0 - (((x-u) * t) / (std::pow(s,2.0))));
    }
  }

  void EmgGradientDescent::extractTrainingSet(
    const std::vector<double>& xs,
    const std::vector<double>& ys,
    std::vector<double>& TrX,
    std::vector<double>& TrY
  ) const
  {
    if (xs.size() < 2) // A valid training set cannot be computed
    {
      throw Exception::SizeUnderflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, xs.size());
    }

    const double intensity_threshold = *std::max_element(ys.begin(), ys.end()) * 0.8;
    std::vector<std::pair<double,double>> points;

    // Add points from the LEFT side, until `intensity_threshold` is reached
    points.emplace_back(xs.front(), ys.front()); // Add FIRST point, no matter the threshold
    Size i = 1;
    for (; i < xs.size() - 1 && ys[i] < intensity_threshold; ++i)
    {
      points.emplace_back(xs[i], ys[i]);
    }

    // Add points from the RIGHT side, until `intensity_threshold` is reached
    points.emplace_back(xs.back(), ys.back()); // Add LAST point, no matter the threshold
    Size j = xs.size() - 2;
    for (; i <= j && ys[j] < intensity_threshold; --j)
    {
      points.emplace_back(xs[j], ys[j]);
    }

    // Compute the derivative for points of intensity greater than `intensity_threshold`
    // According to the value of the highest derivative, it will be decided if
    // a given point is to be added or to be skipped
    // `derivatives` contains the information for both directions
    std::vector<double> derivatives(xs.size() + 1); // One more element to account for derivatives from right to left
    derivatives.front() = 1.0;
    derivatives.back() = -1.0;
    for (Size k = i - 1; k < xs.size() && k <= j + 1 && i > 1; ++k)
    {
      derivatives[k] = (ys[k] - ys[k - 1]) / (xs[k] - xs[k - 1]);
    }

    const double max_derivative = *std::max_element(
      derivatives.begin() + i,
      derivatives.begin() + j + 2,
      [](const double a, const double b)
      {
        return std::fabs(a) < std::fabs(b);
      }
    );
    const double derivative_percent { 0.3 };
    const double derivative_threshold = std::fabs(max_derivative) * derivative_percent;

    // Starting from `i` and proceeding toward the RIGHT side,
    // add points until the derivative conditions are satisfied
    for (
      ; i < xs.size() - 1 &&
      i <= j &&
      derivatives[i] > 0.0 &&
      (std::fabs(derivatives[i]) >= derivative_threshold ||
       derivatives[i] / derivatives[i - 1] >= 0.6);
      ++i
    )
    {
      points.emplace_back(xs[i], ys[i]);
    }

    // Starting from `j` and proceeding toward the LEFT side,
    // add points until the derivative conditions are satisfied
    for (
      ; j > 0 &&
      i <= j &&
      derivatives[j + 1] < 0.0 &&
      (std::fabs(derivatives[j + 1]) >= derivative_threshold ||
       derivatives[j + 1] / derivatives[j + 2] >= 0.6);
      --j
    )
    {
      points.emplace_back(xs[j], ys[j]);
    }

    // Create the output vectors containing the training set
    TrX.clear();
    TrY.clear();
    for (const std::pair<double,double>& point : points)
    {
      TrX.push_back(point.first);
      TrY.push_back(point.second);
    }
  }

  double EmgGradientDescent::computeMuMaxDistance(const std::vector<double>& xs) const
  {
    const std::pair<
      std::vector<double>::const_iterator,
      std::vector<double>::const_iterator
    > p = std::minmax_element(xs.begin(), xs.end());
    if (p.first == xs.end() || p.second == xs.end()) return 0.0;
    const double min_pos = *p.first;
    const double max_pos = *p.second;
    // Return the maximum distance permitted for the Mean parameter, to avoid
    // diverging from the optimal solution in gradient descent
    return (max_pos - min_pos) * 0.35;
  }

  UInt EmgGradientDescent::estimateEmgParameters(
    const std::vector<double>& xs,
    const std::vector<double>& ys,
    double& best_h,
    double& best_mu,
    double& best_sigma,
    double& best_tau
  ) const
  {
    // Initial parameters
    double h { *std::max_element(ys.begin(), ys.end()) };
    double mu { computeInitialMean(xs, ys) };
    double sigma { mu * 1e-2 };
    double tau { sigma * 2.0 };

    const double h_lower_boundary { h }; // Parameter `h` won't decrease below this value

    std::vector<double> TrX, TrY; // Training set (positions and intensities)
    extractTrainingSet(xs, ys, TrX, TrY);

    // Variables containing the "previous" differentials
    double prev_diff_E_h, prev_diff_E_mu, prev_diff_E_sigma, prev_diff_E_tau, previous_E;
    prev_diff_E_h = prev_diff_E_mu = prev_diff_E_sigma = prev_diff_E_tau = 0.0;

    // Part of computation in iRpropPlus()
    // The parameter will change as much as these terms between iterations
    double term_h, term_mu, term_sigma, term_tau;
    term_h = term_mu = term_sigma = term_tau = 0.0;

    // These variables will contain the values obtained at the best iteration
    // The best iteration is decided by the smallest E found
    // Therefore, here `best_E` is initialized with the maximum value for type `double`
    double best_E;
    previous_E = best_h = best_mu = best_sigma = best_tau = best_E = std::numeric_limits<double>::max();

    // Keep track of the current iteration index, and the best iteration index
    UInt iter_idx, best_iter;
    iter_idx = best_iter = 0;

    // This value will increase according to the number of iterations occurred,
    // to avoid spamming the terminal with too much debug information
    UInt info_iter_threshold { 1 };

    // Learning rates (used in gradient descent and iRprop+)
    double lr_h, lr_mu, lr_sigma, lr_tau;
    lr_h = lr_mu = lr_sigma = lr_tau = 0.0125; // iRprop+ paper recommends 0.0125

    // Variables to limit the change in position `mu`
    const double mu_max_dist = computeMuMaxDistance(TrX);
    const double mu_left_boundary { mu - mu_max_dist };
    const double mu_right_boundary { mu + mu_max_dist };

    // The standard deviation between a selection of the precedent Es is computed.
    // If said standard deviation is lower than a certain value,
    // the computation of gradient descent is terminated
    const Size last_few_Es_dim { 10 };
    std::vector<double> last_few_Es(last_few_Es_dim, 0.0);
    Size last_few_Es_idx = 0;
    const double Es_std_dev_min = 1.0; // NOTE: magic value

    if (print_debug_ == 1)
    {
      std::cout << "GRADIENT DESCENT\nInput vectors size: " << xs.size() << "; Training set size: " << TrX.size() << std::endl;
      std::cout << "The possible mu range is [" << mu_left_boundary << " " << mu_right_boundary << "]" << std::endl;
    }

    while (++iter_idx <= max_gd_iter_)
    {
      // Break if parameters are `nan` or `inf`
      if (
        std::isnan(h) || std::isnan(mu) || std::isnan(sigma) || std::isnan(tau) ||
        std::isinf(h) || std::isinf(mu) || std::isinf(sigma) || std::isinf(tau)
      )
      {
        std::cout << std::endl << "[" << iter_idx << "]" << std::endl;
        std::cout << "One or more parameters are invalid." << std::endl;
        std::cout << "Bad: h=" << h << " mu=" << mu << " sigma=" << sigma << " tau=" << tau << std::endl;
        break;
      }

      // Compute the cost given the current parameters
      const double current_E = Loss_function(TrX, TrY, h, mu, sigma, tau);

      // Break if the computed cost is an invalid value
      if (std::isnan(current_E) || std::isinf(current_E))
      {
        std::cout << std::endl << "[" << iter_idx << "]" << std::endl;
        std::cout << "Bad: E value is invalid. current_E=" << current_E << std::endl;
        break;
      }

      // If the current iteration is the best one, save the relevant values
      if (current_E < best_E)
      {
        best_h = h;
        best_mu = mu;
        best_sigma = sigma;
        best_tau = tau;
        best_E = current_E;
        best_iter = iter_idx;
      }

      // Compute the partial derivatives given the current parameters
      double diff_E_h = E_wrt_h(TrX, TrY, h, mu, sigma, tau);
      double diff_E_mu = E_wrt_mu(TrX, TrY, h, mu, sigma, tau);
      double diff_E_sigma = E_wrt_sigma(TrX, TrY, h, mu, sigma, tau);
      double diff_E_tau = E_wrt_tau(TrX, TrY, h, mu, sigma, tau);

      // Logging info to the terminal
      if (print_debug_ == 1 && iter_idx % info_iter_threshold == 0)
      {
        std::cout << std::endl << "[" << iter_idx << "] [prev. E=" << current_E << "]" << std::endl;
        std::cout << "[diff_E_h=" << diff_E_h << "] [diff_E_mu=" << diff_E_mu << "] [diff_E_sigma=" << diff_E_sigma << "] [diff_E_tau=" << diff_E_tau << "]" << std::endl;
        std::cout << "[h=" << h << "] \t[mu=" << mu << "] \t[sigma=" << sigma << "] \t[tau=" << tau << "]" << std::endl;
        std::cout << "[lr_h=" << lr_h << "] \t[lr_mu=" << lr_mu << "] \t[lr_sigma=" << lr_sigma << "] \t[lr_tau=" << lr_tau << "]" << std::endl;
        // Avoid spamming the terminal: increase `info_iter_threshold` dynamically (until 10k)
        if (iter_idx < 10000 && iter_idx / info_iter_threshold >= 10) info_iter_threshold *= 10;
      }

      // If the cost function doesn't change enough, the gradient descent algorithm is terminated
      // This is decided by computing the standard deviation between a selection of the last few Es
      if (iter_idx % 50 == 0) // NOTE: magic value
      {
        last_few_Es[last_few_Es_idx++ % last_few_Es_dim] = current_E;
        const double mean = std::accumulate(last_few_Es.begin(), last_few_Es.end(), 0.0) / last_few_Es_dim;
        double squared_diffs {0.0};
        for (const double current_E : last_few_Es)
        {
          squared_diffs += std::pow(current_E - mean, 2.0);
        }
        const double Es_std_dev = std::sqrt(squared_diffs / static_cast<double>(last_few_Es_dim));
        if (Es_std_dev < Es_std_dev_min)
        {
          if (print_debug_ == 1)
          {
            std::cout << std::endl << "[" << iter_idx << "] The cost function is not changing enough, anymore. Breaking.";
            std::cout << std::endl << "[" << iter_idx << "] [mean=" << mean << "] [Es_std_dev=" << Es_std_dev << "]" << std::endl;
          }
          break;
        }
      }

      // Simultaneous update of all parameters for gradient descent
      iRpropPlus(prev_diff_E_h, diff_E_h, lr_h, term_h, h, current_E, previous_E);
      iRpropPlus(prev_diff_E_mu, diff_E_mu, lr_mu, term_mu, mu, current_E, previous_E);
      iRpropPlus(prev_diff_E_sigma, diff_E_sigma, lr_sigma, term_sigma, sigma, current_E, previous_E);
      iRpropPlus(prev_diff_E_tau, diff_E_tau, lr_tau, term_tau, tau, current_E, previous_E);

      // Apply the parameters' constraints
      h = std::max(h_lower_boundary, h);
      if (mu < mu_left_boundary || mu_right_boundary < mu)
      {
        mu = mu < mu_left_boundary ? mu_left_boundary : mu_right_boundary;
      }
      sigma = std::min(std::max(1e-4, sigma), 20.0); // NOTE: magic value
      tau = std::min(std::max(sigma, tau), sigma * 15.0); // NOTE: magic value

      // Saving values to be used at the next iteration
      prev_diff_E_h = diff_E_h;
      prev_diff_E_mu = diff_E_mu;
      prev_diff_E_sigma = diff_E_sigma;
      prev_diff_E_tau = diff_E_tau;
      previous_E = current_E;
    }
    if (print_debug_ == 1)
    {
      std::cout << std::endl << "[" << best_iter << "] RESULT: best_E=" << best_E << std::endl;
      // TODO: Remove the following "GEOGEBRA" line
      std::cout << "[" << best_iter << "] GEOGEBRA: Execute[{\"h = " << best_h << "\", \"mu = " << best_mu << "\",\"sigma = " << best_sigma << "\", \"tau = " << best_tau << "\"}]" << std::endl;
    }
    // The method has a maximum number of iterations permitted
    // (see class parameter `max_gd_iter`).
    // Said limit is rarely reached, and instead the method will finish after
    // a lower number of iterations. The method returns such number.
    return iter_idx;
  }

  template <typename PeakContainerT>
  void EmgGradientDescent::fitEMGPeakModel(
    const PeakContainerT& input_peak,
    PeakContainerT& output_peak,
    const double left_pos,
    const double right_pos
  ) const
  {
    // Extract points
    typename PeakContainerT::const_iterator start_it = left_pos ? input_peak.PosBegin(left_pos) : input_peak.begin();
    typename PeakContainerT::const_iterator end_it = right_pos ? input_peak.PosEnd(right_pos) : input_peak.end();
    std::vector<double> xs, ys;
    for (typename PeakContainerT::const_iterator it = start_it; it != end_it; ++it)
    {
      xs.push_back(it->getPos());
      ys.push_back(it->getIntensity());
    }

    // EMG parameter estimation with gradient descent
    double h, mu, sigma, tau;
    estimateEmgParameters(xs, ys, h, mu, sigma, tau);

    // Estimate the intensities for each point
    std::vector<double> out_xs;
    std::vector<double> out_ys;
    applyEstimatedParameters(xs, h, mu, sigma, tau, out_xs, out_ys);

    // Prepare the output peak
    output_peak = input_peak;
    output_peak.clear(false); // Remove the points, but keep the metadata
    for (Size i = 0; i < out_xs.size(); ++i)
    {
      // NOTE: casting to avoid -Wnarrowing compiler warning/error
      // TODO: remove cast once issue #3379 is solved
      // https://github.com/OpenMS/OpenMS/issues/3379
      typename PeakContainerT::PeakType point { out_xs[i], static_cast<float>(out_ys[i]) };
      output_peak.push_back(point);
    }

    // Add the EMG parameters as metadata
    typename PeakContainerT::FloatDataArray fda;
    fda.setName("emg_parameters");
    fda.push_back(h);
    fda.push_back(mu);
    fda.push_back(sigma);
    fda.push_back(tau);
    output_peak.getFloatDataArrays().push_back(fda);

    if (print_debug_ == 1)
    {
      std::cout << std::endl << "Input size: " << input_peak.size() << ". ";
      std::cout << "Number of additional points: " << (output_peak.size() - input_peak.size()) << "\n\n" << std::endl;
    }
  }

  template void OPENMS_DLLAPI EmgGradientDescent::fitEMGPeakModel<MSChromatogram>(
    const MSChromatogram& input_peak,
    MSChromatogram& output_peak,
    const double left_pos,
    const double right_pos
  ) const;

  template void OPENMS_DLLAPI EmgGradientDescent::fitEMGPeakModel<MSSpectrum>(
    const MSSpectrum& input_peak,
    MSSpectrum& output_peak,
    const double left_pos,
    const double right_pos
  ) const;
}
