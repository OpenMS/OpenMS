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
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <numeric>
#include <cmath>

#include <OpenMS/SIMULATION/EGHModel.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

using std::endl;

namespace OpenMS
{
  EGHModel::EGHModel() :
    InterpolationModel()
  {
    setName(getProductName());

    // Since the interpolation table is (re-)initialized after setting
    // parameters, we set an empty bounding_box to avoid silly computations.
    defaults_.setValue("interpolation_step", 0.1, "Sampling rate for the interpolation of the model function.", ListUtils::create<String>("advanced"));

    defaults_.setValue("statistics:mean", 0.0f, "Centroid position of the model.", ListUtils::create<String>("advanced"));
    defaults_.setValue("statistics:variance", 1.0f, "The variance of the model.", ListUtils::create<String>("advanced"));

    defaults_.setValue("egh:height", 1000.0f, "Height of the exponential-Gaussian hybrid.");
    defaults_.setValue("egh:retention", 1200.0f, "Retention time of the exponential-Gaussian hybrid.", ListUtils::create<String>("advanced"));

    defaults_.setValue("egh:guess_parameter", "true", "If set to true, the EGHModel will try to estimate the model parameters (tau and sigma-square) based on A,B, and alpha. If set to false, it will use the parameters tau and sigma-square directly.");
    defaults_.setValidStrings("egh:guess_parameter", ListUtils::create<String>("true,false"));

    defaults_.setValue("egh:A", 100.0f, "Horizontal distance between the vertical line at the peak maximum (H) and the leading line where the EGH has H*alpha (e.g. the left half-max for alpha=0.5).");
    defaults_.setValue("egh:B", 100.0f, "Horizontal distance between the vertical line at the peak maximum (H) and the trailing line where the EGH has H*alpha (e.g. the right half-max for alpha=0.5).");
    defaults_.setValue("egh:alpha", 0.5, "See egh:A and egh:B.");
    defaults_.setMinFloat("egh:alpha", 0.0);
    defaults_.setMaxFloat("egh:alpha", 1.0);


    defaults_.setValue("egh:tau", 0.0, "Time constant of the exponential decay (tau is zero for gaussian peaks).", ListUtils::create<String>("advanced"));
    defaults_.setValue("egh:sigma_square", 1803.4, "Standard deviation of the peak.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("egh:sigma_square", 0.0);

    defaults_.setValue("bounding_box:compute", "true", "If true, the EGHModel will compute its own bounding box.");
    defaults_.setValidStrings("bounding_box:compute", ListUtils::create<String>("true,false"));

    defaults_.setValue("bounding_box:min", 0.0, "Lower end of bounding box enclosing the data used to fit the model.", ListUtils::create<String>("advanced"));
    defaults_.setValue("bounding_box:max", 0.0, "Upper end of bounding box enclosing the data used to fit the model.", ListUtils::create<String>("advanced"));

    defaultsToParam_();
  }

  EGHModel::EGHModel(const EGHModel & source) :
    InterpolationModel(source)
  {
    setParameters(source.getParameters());
    updateMembers_();
  }

  EGHModel::~EGHModel()
  {
  }

  EGHModel & EGHModel::operator=(const EGHModel & source)
  {
    if (&source == this)
      return *this;

    InterpolationModel::operator=(source);
    setParameters(source.getParameters());
    updateMembers_();

    return *this;
  }

  void EGHModel::setSamples()
  {
    ContainerType & data = interpolation_.getData();
    data.clear();
    if (max_ == min_)
      return;

    data.reserve(UInt((max_ - min_) / interpolation_step_ + 1));
    CoordinateType pos = min_;

    // EGH computation
    for (UInt i = 0; pos < max_; ++i)
    {
      pos = min_ + i * interpolation_step_;
      CoordinateType t_diff = pos - apex_rt_;
      CoordinateType egh_value;

      evaluateEGH_(t_diff, egh_value);

      data.push_back(egh_value);
    }

//      B_a = 150;
//       A_a = 250;
//       a = 0.5;
//
//       sig_2 = ( -1 / (2 * log(a)) ) * (B_a * A_a);
//       tau =  ( -1 / log(a) ) * (B_a - A_a);
//
//       t_R = 0.0;
//       H = 1000;
//
//       t_diff = (t - t_R);
//
//       if((2*sig_2 + tau*t_diff) > 0)
//           numerator = - t_diff^2;
//           denominator = 2 * sig_2 + tau * t_diff;
//
//           f = H * exp(numerator / denominator);
//       else
//           f = 0;
//       end

    interpolation_.setScale(interpolation_step_);
    interpolation_.setOffset(min_);
    return;
  }

  void EGHModel::setOffset(CoordinateType offset)
  {
    double diff = offset - getInterpolation().getOffset();
    min_ += diff;
    max_ += diff;

    // sync with params
    param_.setValue("bounding_box:min", min_);
    param_.setValue("bounding_box:max", max_);

    statistics_.setMean(statistics_.mean() + diff);

    InterpolationModel::setOffset(offset);
    param_.setValue("statistics:mean", statistics_.mean());
  }

  EGHModel::CoordinateType EGHModel::getCenter() const
  {
    return statistics_.mean();
  }

  void EGHModel::updateMembers_()
  {
    InterpolationModel::updateMembers_();

    statistics_.setMean(param_.getValue("statistics:mean"));
    statistics_.setVariance(param_.getValue("statistics:variance"));

    height_ = param_.getValue("egh:height");
    apex_rt_ = param_.getValue("egh:retention");

    if (param_.getValue("egh:guess_parameter") == "true")
    {
      A_ = param_.getValue("egh:A");
      B_ = param_.getValue("egh:B");

      CoordinateType alpha = param_.getValue("egh:alpha");
      CoordinateType log_alpha = log(alpha);

      tau_ = (-1 / log_alpha) * (B_ - A_);
      sigma_square_ = (-1 / (2 * log_alpha)) * (B_ * A_);

      // sync with params
      param_.setValue("egh:sigma_square", sigma_square_);
      param_.setValue("egh:tau", tau_);
    }
    else
    {
      tau_ = param_.getValue("egh:tau");
      sigma_square_ = param_.getValue("egh:sigma_square");

      // these values are needed to estimate the bounding box
      A_ = sqrt(sigma_square_);
      B_ = A_;
    }
    sigma_square_2_ = 2 * sigma_square_;

    if (param_.getValue("bounding_box:compute") == "true")
    {
      computeBoundaries_();
      // sync with params
      param_.setValue("bounding_box:min", min_);
      param_.setValue("bounding_box:max", max_);

    }
    else
    {
      min_ = param_.getValue("bounding_box:min");
      max_ = param_.getValue("bounding_box:max");
    }

//      LOG_DEBUG << "Computed EGH-Parameters:\n";
//      LOG_DEBUG << "A:            " << A_ << "\n";
//      LOG_DEBUG << "B:            " << B_ << "\n";
//      LOG_DEBUG << "retention:    " << apex_rt_ << "\n";
//      LOG_DEBUG << "sigma_square: " << sigma_square_ << "\n";
//      LOG_DEBUG << "tau:          " << tau_ << "\n";
//      LOG_DEBUG << "Feature bounding box is <" << min_ << "," << max_ << ">" << std::endl;

    setSamples();
  }

  void EGHModel::computeBoundaries_()
  {
    // reset boundaries
    min_ = max_ = 0.0;

    CoordinateType egh_value;
    CoordinateType threshold = height_ / 1000.0;

    // go left .. A_ defines the step width
    egh_value = height_;
    min_ = -A_;

    while (egh_value > threshold)
    {
      min_ -= A_;
//        LOG_DEBUG << "Decreased feature (" << apex_rt_ << ") min_ to " << (min_ + apex_rt_) << "\n";

      evaluateEGH_(min_, egh_value);

//        LOG_DEBUG << "egh(" << min_ << ")=" << egh_value << " (H = " << height_ << ")" << endl;
    }

    // go right .. B_ defines the step width
    egh_value = height_;
    max_ = B_;
    while (egh_value > threshold)
    {
      max_ += B_;
//        LOG_DEBUG << "Increased feature (" << apex_rt_ << ") max_ to " << (max_ + apex_rt_) << "\n";

      evaluateEGH_(max_, egh_value);

//        LOG_DEBUG << "egh(" << max_ << ")=" << egh_value << " (H = " << height_ << ")" << endl;
    }

    // set boundaries at the correct position on the RT scale
    max_ += apex_rt_;
    min_ += apex_rt_;

    if (min_ < 0.0)  // check if we are below the absolute lower scan limit -> 0.0
    {
      min_ = 0.0;
    }
  }

}


//////////////////////////////////////////////////////////////////////////////
