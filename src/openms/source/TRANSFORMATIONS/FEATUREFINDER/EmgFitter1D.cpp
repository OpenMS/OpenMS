// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Factory.h>

namespace OpenMS
{
  int EmgFitter1D::EgmFitterFunctor::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec)
  {
    Size n = m_data->n;
    EmgFitter1D::RawDataArrayType set = m_data->set;

    EmgFitter1D::CoordinateType h = x(0);
    EmgFitter1D::CoordinateType w = x(1);
    EmgFitter1D::CoordinateType s = x(2);
    EmgFitter1D::CoordinateType z = x(3);

    EmgFitter1D::CoordinateType Yi = 0.0;

    // iterate over all points of the signal
    for (Size i = 0; i < n; i++)
    {
      double t = set[i].getPos();

      // Simplified EMG
      Yi = (h * w / s) * sqrt(2.0 * Constants::PI) * exp((pow(w, 2) / (2 * pow(s, 2))) - ((t - z) / s)) / (1 + exp((-2.4055 / sqrt(2.0)) * (((t - z) / w) - w / s)));

      fvec(i) = Yi - set[i].getIntensity();
    }
    return 0;
  }

  // compute Jacobian matrix for the different parameters
  int EmgFitter1D::EgmFitterFunctor::df(const Eigen::VectorXd& x, Eigen::MatrixXd& J)
  {
    Size n =  m_data->n;
    EmgFitter1D::RawDataArrayType set = m_data->set;

    EmgFitter1D::CoordinateType h = x(0);
    EmgFitter1D::CoordinateType w = x(1);
    EmgFitter1D::CoordinateType s = x(2);
    EmgFitter1D::CoordinateType z = x(3);

    const EmgFitter1D::CoordinateType emg_const = 2.4055;
    const EmgFitter1D::CoordinateType sqrt_2pi = sqrt(2 * Constants::PI);
    const EmgFitter1D::CoordinateType sqrt_2 = sqrt(2.0);

    EmgFitter1D::CoordinateType exp1, exp2, exp3 = 0.0;
    EmgFitter1D::CoordinateType derivative_height, derivative_width, derivative_symmetry, derivative_retention = 0.0;

    // iterate over all points of the signal
    for (Size i = 0; i < n; i++)
    {
      EmgFitter1D::CoordinateType t = set[i].getPos();

      exp1 = exp(((w * w) / (2 * s * s)) - ((t - z) / s));
      exp2 = (1 + exp((-emg_const / sqrt_2) * (((t - z) / w) - w / s)));
      exp3 = exp((-emg_const / sqrt_2) * (((t - z) / w) - w / s));

      // f'(h)
      derivative_height = w / s * sqrt_2pi * exp1 / exp2;

      // f'(w)
      derivative_width = h / s * sqrt_2pi * exp1 / exp2 + (h * w * w) / (s * s * s) * sqrt_2pi * exp1 / exp2 + (emg_const * h * w) / s * sqrt_2pi * exp1 * (-(t - z) / (w * w) - 1 / s) * exp3 / ((exp2 * exp2) * sqrt_2);

      // f'(s)
      derivative_symmetry = -h * w / (s * s) * sqrt_2pi * exp1 / exp2 + h * w / s * sqrt_2pi * (-(w * w) / (s * s * s) + (t - z) / (s * s)) * exp1 / exp2 + (emg_const * h * w * w) / (s * s * s) * sqrt_2pi * exp1 * exp3 / ((exp2 * exp2) * sqrt_2);

      // f'(z)
      derivative_retention = h * w / (s * s) * sqrt_2pi * exp1 / exp2 - (emg_const * h) / s * sqrt_2pi * exp1 * exp3 / ((exp2 * exp2) * sqrt_2);

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
    setName(getProductName());
    defaults_.setValue("statistics:variance", 1.0, "Variance of the model.", ListUtils::create<String>("advanced"));
    defaultsToParam_();
  }

  EmgFitter1D::EmgFitter1D(const EmgFitter1D& source) :
    LevMarqFitter1D(source)
  {
    setParameters(source.getParameters());
    updateMembers_();
  }

  EmgFitter1D::~EmgFitter1D()
  {
  }

  EmgFitter1D& EmgFitter1D::operator=(const EmgFitter1D& source)
  {
    if (&source == this)
      return *this;

    LevMarqFitter1D::operator=(source);
    setParameters(source.getParameters());
    updateMembers_();

    return *this;
  }

  EmgFitter1D::QualityType EmgFitter1D::fit1d(const RawDataArrayType& set, InterpolationModel*& model)
  {
    // Calculate bounding box
    CoordinateType min_bb = set[0].getPos(), max_bb = set[0].getPos();
    for (Size pos = 1; pos < set.size(); ++pos)
    {
      CoordinateType tmp = set[pos].getPos();
      if (min_bb > tmp)
        min_bb = tmp;
      if (max_bb < tmp)
        max_bb = tmp;
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
//    CoordinateType x_init[4] = { height_, width_, symmetry_, retention_ };
    Eigen::VectorXd x_init(4);
    x_init(0) = height_;
    x_init(1) = width_;
    x_init(2) = symmetry_;
    x_init(3) = retention_;
    if (symmetric_ == false)
    {
      EgmFitterFunctor functor(4, &d);
      optimize_(x_init, functor);
    }

    // Set optimized parameters
    height_ = x_init[0];
    width_ = x_init[1];
    symmetry_ = x_init[2];
    retention_ = x_init[3];

#ifdef DEBUG_FEATUREFINDER
    if (getGslStatus_() != "success")
    {
      std::cout << "status: " << getGslStatus_() << std::endl;
    }
#endif

    // build model
    model = static_cast<InterpolationModel*>(Factory<BaseModel<1> >::create("EmgModel"));
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

  void EmgFitter1D::setInitialParameters_(const RawDataArrayType& set)
  {
    // sum over all intensities
    CoordinateType sum = 0.0;
    for (Size i = 0; i < set.size(); ++i)
      sum += set[i].getIntensity();

    // calculate the median
    Size median = 0;
    float count = 0.0;
    for (Size i = 0; i < set.size(); ++i)
    {
      count += set[i].getIntensity();
      if (count <= sum / 2)
        median = i;
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
