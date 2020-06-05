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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BiGaussModel.h>

namespace OpenMS
{
  BiGaussModel::BiGaussModel() :
    InterpolationModel(), statistics1_(), statistics2_()
  {
    setName(getProductName());

    defaults_.setValue("bounding_box:min", 0.0, "Lower end of bounding box enclosing the data used to fit the model.", ListUtils::create<String>("advanced"));
    defaults_.setValue("bounding_box:max", 1.0, "Upper end of bounding box enclosing the data used to fit the model.", ListUtils::create<String>("advanced"));
    defaults_.setValue("statistics:mean", 0.0, "Centroid position of the model, this also separates both halves of the model.", ListUtils::create<String>("advanced"));
    defaults_.setValue("statistics:variance1", 1.0, "Variance of the first gaussian, used for the lower half of the model.", ListUtils::create<String>("advanced"));
    defaults_.setValue("statistics:variance2", 1.0, "Variance of the second gaussian, used for the upper half of the model.", ListUtils::create<String>("advanced"));

    defaultsToParam_();
  }

  BiGaussModel::BiGaussModel(const BiGaussModel & source) :
    InterpolationModel(source)
  {
    setParameters(source.getParameters());
    updateMembers_();
  }

  BiGaussModel::~BiGaussModel()
  {
  }

  BiGaussModel & BiGaussModel::operator=(const BiGaussModel & source)
  {
    if (&source == this)
      return *this;

    InterpolationModel::operator=(source);
    setParameters(source.getParameters());
    updateMembers_();

    return *this;
  }

  void BiGaussModel::setSamples()
  {
    LinearInterpolation::container_type & data = interpolation_.getData();
    data.clear();
    if (max_ == min_)
      return;

    data.reserve(UInt((max_ - min_) / interpolation_step_ + 1));
    CoordinateType pos = min_;

    for (UInt i = 0; pos < max_; ++i)
    {
      pos = min_ + i * interpolation_step_;
      if (pos < statistics1_.mean())
      {
        data.push_back(statistics1_.normalDensity_sqrt2pi(pos));
      }
      else
      {
        data.push_back(statistics2_.normalDensity_sqrt2pi(pos));
      }
    }
    // scale data so that integral over distribution equals one
    // multiply sum by interpolation_step_ -> rectangular approximation of integral
    IntensityType factor = scaling_ / interpolation_step_ /
                           std::accumulate(data.begin(), data.end(), IntensityType(0));

    for (LinearInterpolation::container_type::iterator it = data.begin(); it != data.end(); ++it)
    {
      *it *= factor;
    }

    interpolation_.setScale(interpolation_step_);
    interpolation_.setOffset(min_);
  }

  void BiGaussModel::updateMembers_()
  {
    InterpolationModel::updateMembers_();

    min_ = param_.getValue("bounding_box:min");
    max_ = param_.getValue("bounding_box:max");
    statistics1_.setMean(param_.getValue("statistics:mean"));
    statistics2_.setMean(param_.getValue("statistics:mean"));
    statistics1_.setVariance(param_.getValue("statistics:variance1"));
    statistics2_.setVariance(param_.getValue("statistics:variance2"));

    setSamples();
  }

  void BiGaussModel::setOffset(CoordinateType offset)
  {
    double diff = offset - getInterpolation().getOffset();
    min_ += diff;
    max_ += diff;
    statistics1_.setMean(statistics1_.mean() + diff);
    statistics2_.setMean(statistics2_.mean() + diff);

    InterpolationModel::setOffset(offset);

    param_.setValue("bounding_box:min", min_);
    param_.setValue("bounding_box:max", max_);
    param_.setValue("statistics:mean", statistics1_.mean());
  }

  BiGaussModel::CoordinateType BiGaussModel::getCenter() const
  {
    return statistics2_.mean();
  }

}
