// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/EGHTraceFitter.h>

#include <unsupported/Eigen/NonLinearOptimization>
#include <Eigen/Core>

#include <OpenMS/CONCEPT/LogStream.h>

#include <numeric> // for "accumulate"

namespace OpenMS
{
  // from table 1 in the Lan & Jorgenson paper:
  const double EGHTraceFitter::EPSILON_COEFS_[] =
  {4.0, -6.293724, 9.232834, -11.342910, 9.123978, -4.173753, 0.827797};


  const Size EGHTraceFitter::NUM_PARAMS_ = 4;

  EGHTraceFitter::EGHTraceFunctor::EGHTraceFunctor(int dimensions,
                                                   const TraceFitter::ModelData* data) :
    TraceFitter::GenericFunctor(dimensions, data->traces_ptr->getPeakCount()), m_data(data)
  {
  }

  EGHTraceFitter::EGHTraceFunctor::~EGHTraceFunctor() = default;

  int EGHTraceFitter::EGHTraceFunctor::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec)
  {
    double H  = x(0);
    double tR = x(1);
    double sigma = x(2);
    double tau = x(3);

    double t_diff, t_diff2, denominator = 0.0;

    double fegh = 0.0;

    UInt count = 0;
    for (Size t = 0; t < m_data->traces_ptr->size(); ++t)
    {
      const FeatureFinderAlgorithmPickedHelperStructs::MassTrace& trace = m_data->traces_ptr->at(t);
      double weight = m_data->weighted ? trace.theoretical_int : 1.0;
      for (Size i = 0; i < trace.peaks.size(); ++i)
      {
        double rt = trace.peaks[i].first;

        t_diff = rt - tR;
        t_diff2 = t_diff * t_diff; // -> (t - t_R)^2

        denominator = 2 * sigma * sigma + tau * t_diff; // -> 2\sigma_{g}^{2} + \tau \left(t - t_R\right)

        if (denominator > 0.0)
        {
          fegh =  m_data->traces_ptr->baseline + trace.theoretical_int * H * exp(-t_diff2 / denominator);
        }
        else
        {
          fegh = 0.0;
        }

        fvec(count) = (fegh - trace.peaks[i].second->getIntensity()) * weight;
        ++count;
      }
    }
    return 0;
  }

  int EGHTraceFitter::EGHTraceFunctor::df(const Eigen::VectorXd& x, Eigen::MatrixXd& J)
  {
    double H  = x(0);
    double tR = x(1);
    double sigma = fabs(x(2)); // must be non-negative!
    double tau = x(3);

    double derivative_H, derivative_tR, derivative_sigma, derivative_tau = 0.0;
    double t_diff, t_diff2, exp1, denominator = 0.0;

    UInt count = 0;
    for (Size t = 0; t < m_data->traces_ptr->size(); ++t)
    {
      const FeatureFinderAlgorithmPickedHelperStructs::MassTrace& trace = m_data->traces_ptr->at(t);
      double weight = m_data->weighted ? trace.theoretical_int : 1.0;
      for (Size i = 0; i < trace.peaks.size(); ++i)
      {
        double rt = trace.peaks[i].first;

        t_diff = rt - tR;
        t_diff2 = t_diff * t_diff; // -> (t - t_R)^2

        denominator = 2 * sigma * sigma + tau * t_diff; // -> 2\sigma_{g}^{2} + \tau \left(t - t_R\right)

        if (denominator > 0)
        {
          exp1 = exp(-t_diff2 / denominator);

          // \partial H f_{egh}(t) = \exp\left( \frac{-\left(t-t_R \right)}{2\sigma_{g}^{2} + \tau \left(t - t_R\right)} \right)
          derivative_H = trace.theoretical_int * exp1;

          // \partial t_R f_{egh}(t) &=& H \exp \left( \frac{-\left(t-t_R \right)}{2\sigma_{g}^{2} + \tau \left(t - t_R\right)} \right) \left( \frac{\left( 4 \sigma_{g}^{2} + \tau \left(t-t_R \right) \right) \left(t-t_R \right)}{\left( 2\sigma_{g}^{2} + \tau \left(t - t_R\right) \right)^2} \right)
          derivative_tR = trace.theoretical_int * H * exp1 * ((4 * sigma * sigma + tau * t_diff) * t_diff) / (denominator * denominator);

          // \partial \sigma_{g}^{2} f_{egh}(t) &=& H \exp \left( \frac{-\left(t-t_R \right)^2}{2\sigma_{g}^{2} + \tau \left(t - t_R\right)} \right) \left( \frac{ 2 \left(t - t_R\right)^2}{\left( 2\sigma_{g}^{2} + \tau \left(t - t_R\right) \right)^2} \right)
          // // \partial \sigma_{g}^{2} f_{egh}(t) &=& H \exp \left( \frac{-\left(t-t_R \right)^2}{2\sigma_{g}^{2} + \tau \left(t - t_R\right)} \right) \left( \frac{ 2 \left(t - t_R\right)^2}{\left( 2\sigma_{g}^{2} + \tau \left(t - t_R\right) \right)^2} \right)
          // derivative_sigma_square = trace.theoretical_int * H * exp1 * 2 * t_diff2 / (denominator * denominator));

          // \partial \sigma_{g} f_{egh}(t) &=& H \exp \left( \frac{-\left(t-t_R \right)^2}{2\sigma_{g}^{2} + \tau \left(t - t_R\right)} \right) \left( \frac{ 4 \sigma_{g} \left(t - t_R\right)^2}{\left( 2\sigma_{g}^{2} + \tau \left(t - t_R\right) \right)^2} \right)
          derivative_sigma = trace.theoretical_int * H * exp1 * 4 * sigma * t_diff2 / (denominator * denominator);

          // \partial \tau f_{egh}(t) &=& H \exp \left( \frac{-\left(t-t_R \right)^2}{2\sigma_{g}^{2} + \tau \left(t - t_R\right)} \right) \left( \frac{ \left(t - t_R\right)^3}{\left( 2\sigma_{g}^{2} + \tau \left(t - t_R\right) \right)^2} \right)
          derivative_tau = trace.theoretical_int * H * exp1 * t_diff * t_diff2 / (denominator * denominator);
        }
        else
        {
          derivative_H = 0.0;
          derivative_tR = 0.0;
          derivative_sigma = 0.0;
          derivative_tau = 0.0;
        }

        // set the jacobian matrix
        J(count, 0) = derivative_H * weight;
        J(count, 1) = derivative_tR * weight;
        J(count, 2) = derivative_sigma * weight;
        J(count, 3) = derivative_tau * weight;
        ++count;
      }
    }
    return 0;
  }

  EGHTraceFitter::EGHTraceFitter() = default;

  EGHTraceFitter::EGHTraceFitter(const EGHTraceFitter& other) :
    TraceFitter(other)
  {
    this->height_ = other.height_;
    this->apex_rt_ = other.apex_rt_;
    this->sigma_ = other.sigma_;
    this->tau_ = other.tau_;
    this->region_rt_span_ = other.region_rt_span_;

    this->sigma_5_bound_ = other.sigma_5_bound_;

    updateMembers_();
  }

  EGHTraceFitter& EGHTraceFitter::operator=(const EGHTraceFitter& source)
  {
    TraceFitter::operator=(source);

    this->height_ = source.height_;
    this->apex_rt_ = source.apex_rt_;
    this->sigma_ = source.sigma_;
    this->tau_ = source.tau_;
    this->region_rt_span_ = source.region_rt_span_;

    this->sigma_5_bound_ = source.sigma_5_bound_;

    updateMembers_();

    return *this;
  }

  EGHTraceFitter::~EGHTraceFitter() = default;

  void EGHTraceFitter::fit(FeatureFinderAlgorithmPickedHelperStructs::MassTraces& traces)
  {
    setInitialParameters_(traces);

    Eigen::VectorXd x_init(NUM_PARAMS_);
    x_init(0) = height_;
    x_init(1) = apex_rt_;
    x_init(2) = sigma_;
    x_init(3) = tau_;

    TraceFitter::ModelData data{};
    data.traces_ptr = &traces;
    data.weighted = this->weighted_;
    EGHTraceFunctor functor(NUM_PARAMS_, &data);

    TraceFitter::optimize_(x_init, functor);
  }

  double EGHTraceFitter::getLowerRTBound() const
  {
    return sigma_5_bound_.first;
  }

  double EGHTraceFitter::getTau() const
  {
    return tau_;
  }

  double EGHTraceFitter::getUpperRTBound() const
  {
    return sigma_5_bound_.second;
  }

  double EGHTraceFitter::getHeight() const
  {
    return height_;
  }

  double EGHTraceFitter::getSigma() const
  {
    return sigma_;
  }

  double EGHTraceFitter::getCenter() const
  {
    return apex_rt_;
  }

  bool EGHTraceFitter::checkMaximalRTSpan(const double max_rt_span)
  {
    return (sigma_5_bound_.second - sigma_5_bound_.first) > max_rt_span * region_rt_span_;
  }

  bool EGHTraceFitter::checkMinimalRTSpan(const std::pair<double, double>& rt_bounds, const double min_rt_span)
  {
    return (rt_bounds.second - rt_bounds.first) < min_rt_span * (sigma_5_bound_.second - sigma_5_bound_.first);
  }

  double EGHTraceFitter::getValue(double rt) const
  {
    // equation 12 from Lan & Jorgenson paper:
    double fegh = 0.0;
    double t_diff = rt - apex_rt_;
    double denominator = 2 * sigma_ * sigma_ + tau_ * t_diff;
    if (denominator > 0.0)
    {
      fegh = height_ * exp(-t_diff * t_diff / denominator);
    }

    return fegh;
  }

  double EGHTraceFitter::getArea()
  {
    // equation 21 from Lan & Jorgenson paper:
    double abs_tau = fabs(tau_);
    double abs_sigma = fabs(sigma_);
    double phi = atan(abs_tau / abs_sigma);
    double epsilon = EPSILON_COEFS_[0];
    double phi_pow = phi;
    for (Size i = 1; i < 7; ++i)
    {
      epsilon += phi_pow * EPSILON_COEFS_[i];
      phi_pow *= phi;
    }
    // 0.62... is approx. sqrt(pi / 8):
    return height_ * (abs_sigma * 0.6266571 + abs_tau) * epsilon;
  }

  double EGHTraceFitter::getFWHM() const
  {
    std::pair<double, double> bounds = getAlphaBoundaries_(0.5);
    return bounds.second - bounds.first;
  }

  String EGHTraceFitter::getGnuplotFormula(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace& trace, const char function_name, const double baseline, const double rt_shift)
  {
    std::stringstream s;
    s << String(function_name)  << "(x)= " << baseline << " + ";
    s << "("; // the overall bracket
    s << "((" << 2 * sigma_ * sigma_ << " + " << tau_ << " * (x - " << (rt_shift + apex_rt_) << " )) > 0) ? "; // condition
    s <<  (trace.theoretical_int *  height_) << " * exp(-1 * (x - " << (rt_shift + apex_rt_) << ")**2 " <<
      "/" <<
      " ( " << 2 * sigma_ * sigma_ << " + " << tau_ << " * (x - " << (rt_shift + apex_rt_) << " )))";
    s << " : 0)";
    return String(s.str());
  }

  std::pair<double, double> EGHTraceFitter::getAlphaBoundaries_(const double alpha) const
  {
    std::pair<double, double> bounds;
    // solved equations A.2 and A.3 from the Lan & Jorgenson paper (Appendix
    // A) for the boundaries A_alpha and B_alpha:
    double L = log(alpha);
    double s = sqrt(((L * tau_) * (L * tau_) / 4) - 2 * L * sigma_ * sigma_);

    double s1, s2;
    s1 = (-1 * (L * tau_) / 2) + s;
    s2 = (-1 * (L * tau_) / 2) - s;

    // the smaller one (should be < 0) = lower bound
    bounds.first = apex_rt_ + std::min(s1, s2);
    // bigger one (should be > 0) = upper bound
    bounds.second = apex_rt_ + std::max(s1, s2);

    return bounds;
  }

  void EGHTraceFitter::getOptimizedParameters_(const Eigen::VectorXd& x_init)
  {
    height_ =  x_init(0);
    apex_rt_ =  x_init(1);
    sigma_ =  x_init(2);
    tau_ =  x_init(3);

    // we set alpha to 0.04 which is conceptually equal to
    // 2.5 sigma for lower and upper bound
    sigma_5_bound_ = getAlphaBoundaries_(0.043937);
  }

  void EGHTraceFitter::setInitialParameters_(FeatureFinderAlgorithmPickedHelperStructs::MassTraces& traces)
  {
    OPENMS_LOG_DEBUG << "EGHTraceFitter->setInitialParameters(...)" << std::endl;
    OPENMS_LOG_DEBUG << "Number of traces: " << traces.size() << std::endl;

    // aggregate data; some peaks (where intensity is zero) can be missing!
    // mapping: RT -> total intensity over all mass traces
    std::list<std::pair<double, double> > total_intensities;
    traces.computeIntensityProfile(total_intensities);

    // compute moving average for smoothing:
    const Size N = total_intensities.size();
    const Size LEN = 2;   // window size: 2 * LEN + 1
    std::vector<double> totals(N + 2 * LEN);   // pad with zeros at ends
    Int index = LEN;
    // OPENMS_LOG_DEBUG << "Summed intensities:\n";
    for (std::list<std::pair<double, double> >::iterator it =
           total_intensities.begin(); it != total_intensities.end(); ++it)
    {
      totals[index++] = it->second;
      // OPENMS_LOG_DEBUG << it->second << std::endl;
    }

    std::vector<double> smoothed(N);
    Size max_index = 0;   // index of max. smoothed intensity
    // OPENMS_LOG_DEBUG << "Smoothed intensities:\n";
    double sum = std::accumulate(&totals[LEN], &totals[2 * LEN], 0.0);
    for (Size i = 0; i < N; ++i)
    {
      sum += totals[i + 2 * LEN];
      smoothed[i] = sum / (2 * LEN + 1);
      sum -= totals[i];
      if (smoothed[i] > smoothed[max_index])
      {
        max_index = i;
      }
      // OPENMS_LOG_DEBUG << smoothed[i] << std::endl;
    }
    OPENMS_LOG_DEBUG << "Maximum at index " << max_index << std::endl;
    height_ = smoothed[max_index] - traces.baseline;
    OPENMS_LOG_DEBUG << "height: " << height_ << std::endl;
    std::list<std::pair<double, double> >::iterator it = total_intensities.begin();
    std::advance(it, max_index);
    apex_rt_ = it->first;
    OPENMS_LOG_DEBUG << "apex_rt: " << apex_rt_ << std::endl;
    region_rt_span_ = (total_intensities.rbegin()->first -
                       total_intensities.begin()->first);
    OPENMS_LOG_DEBUG << "region_rt_span: " << region_rt_span_ << std::endl;

    // find RT values where intensity is at half-maximum:
    index = static_cast<Int>(max_index);
    while ((index > 0) && (smoothed[index] > height_ * 0.5))
    {
      --index;
    }
    double left_height = smoothed[index];
    it = total_intensities.begin();
    std::advance(it, index);
    double left_rt = it->first;
    OPENMS_LOG_DEBUG << "Left half-maximum at index " << index << ", RT " << left_rt
              << std::endl;
    index = static_cast<Int>(max_index);
    while ((index < Int(N - 1)) && (smoothed[index] > height_ * 0.5))
    {
      ++index;
    }
    double right_height = smoothed[index];
    it = total_intensities.end();
    std::advance(it, index - Int(N));
    double right_rt = it->first;
    OPENMS_LOG_DEBUG << "Right half-maximum at index " << index << ", RT "
              << right_rt << std::endl;

    double A = apex_rt_ - left_rt;
    double B = right_rt - apex_rt_;
    //OPENMS_LOG_DEBUG << "A: " << A << std::endl;
    //OPENMS_LOG_DEBUG << "B: " << B << std::endl;

    // compute estimates for tau / sigma based on A and B:
    double alpha = (left_height + right_height) * 0.5 / height_;   // ~0.5
    double log_alpha = log(alpha);

    tau_ = -1 / log_alpha * (B - A);
    //EGH function fails when tau==0
    if (tau_ == 0)
    {
      tau_ = std::numeric_limits<double>::epsilon();
    }
    OPENMS_LOG_DEBUG << "tau: " << tau_ << std::endl;
    sigma_ = sqrt(-0.5 / log_alpha * B * A);
    OPENMS_LOG_DEBUG << "sigma: " << sigma_ << std::endl;
  }

  void EGHTraceFitter::updateMembers_()
  {
    TraceFitter::updateMembers_();
  }

} // namespace OpenMS
