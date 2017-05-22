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
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLinear.h>

#include <Wm5Vector2.h>
#include <Wm5ApprLineFit2.h>

namespace OpenMS
{

  TransformationModelLinear::TransformationModelLinear(const TransformationModel::DataPoints& data, const Param& params)
  {
    params_ = params;
    data_given_ = !data.empty();
    if (!data_given_ && params.exists("slope") && (params.exists("intercept")))
    {
      // don't estimate parameters, use given values
      slope_ = params.getValue("slope");
      intercept_ = params.getValue("intercept");
    }
    else // estimate parameters from data
    {
      Param defaults;
      getDefaultParameters(defaults);
      params_.setDefaults(defaults);
      symmetric_ = params_.getValue("symmetric_regression") == "true";

      size_t size = data.size();
      std::vector<Wm5::Vector2d> points;
      if (size == 0) // no data
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "no data points for 'linear' model");
      }
      else if (size == 1) // degenerate case, but we can still do something
      {
        slope_ = 1.0;
        intercept_ = data[0].second - data[0].first;
      }
      else // compute least-squares fit
      {
        for (size_t i = 0; i < size; ++i)
        {
          points.push_back(Wm5::Vector2d(data[i].first, data[i].second));
        }
        if (!Wm5::HeightLineFit2<double>(static_cast<int>(size), &points.front(), slope_, intercept_))
        {
          throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "TransformationModelLinear", "Unable to fit linear transformation to data points.");
        }
      }
      // update params
      params_.setValue("slope", slope_);
      params_.setValue("intercept", intercept_);
    }
  }

  TransformationModelLinear::~TransformationModelLinear()
  {
  }

  double TransformationModelLinear::evaluate(double value) const
  {
    return slope_ * value + intercept_;
  }

  void TransformationModelLinear::invert()
  {
    if (slope_ == 0)
    {
      throw Exception::DivisionByZero(__FILE__, __LINE__,
                                      OPENMS_PRETTY_FUNCTION);
    }
    intercept_ = -intercept_ / slope_;
    slope_ = 1.0 / slope_;
    // update parameters:
    params_.setValue("slope", slope_);
    params_.setValue("intercept", intercept_);
  }

  void TransformationModelLinear::getParameters(double& slope, double& intercept) const
  {
    slope = slope_;
    intercept = intercept_;
  }

  void TransformationModelLinear::getDefaultParameters(Param& params)
  {
    params.clear();
    params.setValue("symmetric_regression", "false", "Perform linear regression"
                                                     " on 'y - x' vs. 'y + x', instead of on 'y' vs. 'x'.");
    params.setValidStrings("symmetric_regression",
                           ListUtils::create<String>("true,false"));
  }

} // namespace
