// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BiGaussFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CONCEPT/Factory.h>
#include <boost/math/special_functions/fpclassify.hpp>

namespace OpenMS
{
  BiGaussFitter1D::BiGaussFitter1D() :
    MaxLikeliFitter1D()
  {
    setName(getProductName());

    defaults_.setValue("statistics:variance1", 1.0, "Variance of the first gaussian, used for the lower half of the model.", ListUtils::create<String>("advanced"));
    defaults_.setValue("statistics:variance2", 1.0, "Variance of the second gaussian, used for the upper half of the model.", ListUtils::create<String>("advanced"));

    defaultsToParam_();
  }

  BiGaussFitter1D::BiGaussFitter1D(const BiGaussFitter1D & source) :
    MaxLikeliFitter1D(source)
  {
    updateMembers_();
  }

  BiGaussFitter1D::~BiGaussFitter1D()
  {
  }

  BiGaussFitter1D & BiGaussFitter1D::operator=(const BiGaussFitter1D & source)
  {
    if (&source == this)
      return *this;

    MaxLikeliFitter1D::operator=(source);
    updateMembers_();

    return *this;
  }

  BiGaussFitter1D::QualityType BiGaussFitter1D::fit1d(const RawDataArrayType & set, InterpolationModel * & model)
  {
    // Calculate bounding box
    min_ = max_ = set[0].getPos();
    for (UInt pos = 1; pos < set.size(); ++pos)
    {
      CoordinateType tmp = set[pos].getPos();
      if (min_ > tmp)
        min_ = tmp;
      if (max_ < tmp)
        max_ = tmp;
    }

    // Enlarge the bounding box by a few multiples of the standard deviation
    {
      stdev1_ = sqrt(statistics1_.variance()) * tolerance_stdev_box_;
      stdev2_ = sqrt(statistics2_.variance()) * tolerance_stdev_box_;
      min_ -= stdev1_;
      max_ += stdev2_;
    }

    // build model
    model = static_cast<InterpolationModel *>(Factory<BaseModel<1> >::create("BiGaussModel"));
    model->setInterpolationStep(interpolation_step_);
    Param tmp;
    tmp.setValue("bounding_box:min", min_);
    tmp.setValue("bounding_box:max", max_);
    tmp.setValue("statistics:mean", statistics1_.mean());
    tmp.setValue("statistics:variance1", statistics1_.variance());
    tmp.setValue("statistics:variance2", statistics2_.variance());
    model->setParameters(tmp);

    // fit offset
    QualityType quality;
    quality = fitOffset_(model, set, stdev1_, stdev2_, interpolation_step_);
    if (boost::math::isnan(quality))
      quality = -1.0;

    return quality;
  }

  void BiGaussFitter1D::updateMembers_()
  {
    MaxLikeliFitter1D::updateMembers_();
    statistics1_.setMean(param_.getValue("statistics:mean"));
    statistics1_.setVariance(param_.getValue("statistics:variance1"));
    statistics2_.setMean(param_.getValue("statistics:mean"));
    statistics2_.setVariance(param_.getValue("statistics:variance2"));
  }

}
