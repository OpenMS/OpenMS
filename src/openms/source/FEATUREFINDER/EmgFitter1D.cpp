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

#include <unsupported/Eigen/NonLinearOptimization>

namespace OpenMS
{
  const EmgFitter1D::CoordinateType EmgFitter1D::EgmFitterFunctor::sqrt2pi = sqrt(2.0 * Constants::PI);
  const EmgFitter1D::CoordinateType EmgFitter1D::EgmFitterFunctor::emg_const = 2.4055;
  const EmgFitter1D::CoordinateType EmgFitter1D::EgmFitterFunctor::sqrt_2 = sqrt(2.0);
  const EmgFitter1D::CoordinateType EmgFitter1D::EgmFitterFunctor::c = -emg_const / sqrt_2;

  int EmgFitter1D::EgmFitterFunctor::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const
  {
    Size n = m_data->n;
    EmgFitter1D::RawDataArrayType set = m_data->set;

    EmgFitter1D::CoordinateType h = x(0);
    EmgFitter1D::CoordinateType w = x(1);
    EmgFitter1D::CoordinateType s = x(2);
    EmgFitter1D::CoordinateType z = x(3);

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

      fvec(i) = Yi - set[i].getIntensity();
    }
    return 0;
  }

  // compute Jacobian matrix for the different parameters
  int EmgFitter1D::EgmFitterFunctor::df(const Eigen::VectorXd& x, Eigen::MatrixXd& J) const
  {
    Size n =  m_data->n;
    EmgFitter1D::RawDataArrayType set = m_data->set;

    EmgFitter1D::CoordinateType h = x(0);
    EmgFitter1D::CoordinateType w = x(1);
    EmgFitter1D::CoordinateType w2 = w*w;
    EmgFitter1D::CoordinateType s = x(2);
    EmgFitter1D::CoordinateType s2 = s*s;
    EmgFitter1D::CoordinateType s3 = s2 * s;
    EmgFitter1D::CoordinateType z = x(3);

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
      J(i, 0) = derivative_height;
      J(i, 1) = derivative_width;
      J(i, 2) = derivative_symmetry;
      J(i, 3) = derivative_retention;
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
    // Calculate bounding box
    CoordinateType min_bb = set[0].getPos(), max_bb = set[0].getPos();
    for (Size pos = 1; pos < set.size(); ++pos)
    {
      CoordinateType tmp = set[pos].getPos();
      if (min_bb > tmp)
      {
        min_bb = tmp;
      }
      if (max_bb < tmp)
      {
        max_bb = tmp;
      }
    }

    // Enlarge the bounding box by a few multiples of the standard deviation
    const CoordinateType stdev = sqrt(statistics_.variance()) * tolerance_stdev_box_;
    min_bb -= stdev;
    max_bb += stdev;


    // Set advanced parameters for residual_  und jacobian_ method
    EmgFitter1D::Data d;
    d.n = set.size();
    d.set = set;

    // Compute start parameters
    setInitialParameters_(set);

    // Optimize parameter with Levenberg-Marquardt algorithm
    Eigen::VectorXd x_init(4);
    x_init(0) = height_;
    x_init(1) = width_;
    x_init(2) = symmetry_;
    x_init(3) = retention_;

    if (!symmetric_)
    {
        EgmFitterFunctor functor(4, &d);
        optimize_(x_init, functor);
    }

    // Set optimized parameters
    height_ = x_init[0];
    width_ = x_init[1];
    symmetry_ = x_init[2];
    retention_ = x_init[3];

    // build model
    model = std::unique_ptr<InterpolationModel>(new EmgModel());
    model->setInterpolationStep(interpolation_step_);

    Param tmp;
    tmp.setValue("bounding_box:min", min_bb);
    tmp.setValue("bounding_box:max", max_bb);
    tmp.setValue("statistics:variance", statistics_.variance());
    tmp.setValue("statistics:mean", statistics_.mean());
    tmp.setValue("emg:height", height_);
    tmp.setValue("emg:width", width_);
    tmp.setValue("emg:symmetry", symmetry_);
    tmp.setValue("emg:retention", retention_);
    model->setParameters(tmp);

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
    if (std::isnan(correlation))
    {
      correlation = -1.0;
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

    while (sum > weight_sum/2.)
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
