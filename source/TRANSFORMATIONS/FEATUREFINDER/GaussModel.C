// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussModel.h>
#include <numeric>

namespace OpenMS
{
  GaussModel::GaussModel() :
    InterpolationModel(),
    statistics_()
  {
    setName(getProductName());

    defaults_.setValue("bounding_box:min", 0.0, "Lower end of bounding box enclosing the data used to fit the model.", StringList::create("advanced"));
    defaults_.setValue("bounding_box:max", 1.0, "Upper end of bounding box enclosing the data used to fit the model.", StringList::create("advanced"));
    defaults_.setValue("statistics:mean", 0.0, "Centroid position of the model (Gaussian).", StringList::create("advanced"));
    defaults_.setValue("statistics:variance", 1.0, "The variance of the Gaussian.", StringList::create("advanced"));

    defaultsToParam_();
  }

  GaussModel::GaussModel(const GaussModel & source) :
    InterpolationModel(source)
  {
    setParameters(source.getParameters());
    updateMembers_();
  }

  GaussModel::~GaussModel()
  {
  }

  GaussModel & GaussModel::operator=(const GaussModel & source)
  {
    if (&source == this)
      return *this;

    setParameters(source.getParameters());
    InterpolationModel::operator=(source);
    updateMembers_();

    return *this;
  }

  void GaussModel::setSamples()
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
      data.push_back(statistics_.normalDensity_sqrt2pi(pos));
    }
    // scale data so that integral over distribution equals one
    // multiply sum by interpolation_step_ -> rectangular approximation of integral
    IntensityType factor = scaling_ / interpolation_step_ /
                           std::accumulate(data.begin(), data.end(), IntensityType(0));

    for (LinearInterpolation::container_type::iterator it = data.begin(); it != data.end(); ++it)
      *it *= factor;
    interpolation_.setScale(interpolation_step_);
    interpolation_.setOffset(min_);
  }

  void GaussModel::updateMembers_()
  {
    InterpolationModel::updateMembers_();

    min_ = param_.getValue("bounding_box:min");
    max_ = param_.getValue("bounding_box:max");
    statistics_.setMean(param_.getValue("statistics:mean"));
    statistics_.setVariance(param_.getValue("statistics:variance"));

    setSamples();
  }

  void GaussModel::setOffset(CoordinateType offset)
  {
    DoubleReal diff = offset - getInterpolation().getOffset();
    min_ += diff;
    max_ += diff;
    statistics_.setMean(statistics_.mean() + diff);

    InterpolationModel::setOffset(offset);

    param_.setValue("bounding_box:min", min_);
    param_.setValue("bounding_box:max", max_);
    param_.setValue("statistics:mean", statistics_.mean());
  }

  GaussModel::CoordinateType GaussModel::getCenter() const
  {
    return statistics_.mean();
  }

}
