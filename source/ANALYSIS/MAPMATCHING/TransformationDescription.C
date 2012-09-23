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
// $Authors: Clemens Groepl, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/CONCEPT/LogStream.h>

using namespace std;

namespace OpenMS
{
  TransformationDescription::TransformationDescription() :
    data_(TransformationDescription::DataPoints()), model_type_("none"),
    model_(new TransformationModel())
  {
  }

  TransformationDescription::TransformationDescription(
    const TransformationDescription::DataPoints & data) :
    data_(data), model_type_("none"), model_(new TransformationModel())
  {
  }

  TransformationDescription::~TransformationDescription()
  {
    delete model_;
  }

  TransformationDescription::TransformationDescription(
    const TransformationDescription & rhs)
  {
    data_ = rhs.data_;
    model_type_ = "none";
    model_ = 0;     // initialize this before the "delete" call in "fitModel"!
    Param params;
    rhs.getModelParameters(params);
    fitModel(rhs.model_type_, params);
  }

  TransformationDescription & TransformationDescription::operator=(
    const TransformationDescription & rhs)
  {
    if (this == &rhs)
      return *this;

    data_ = rhs.data_;
    model_type_ = "none";
    Param params;
    rhs.getModelParameters(params);
    fitModel(rhs.model_type_, params);

    return *this;
  }

  void TransformationDescription::fitModel(const String & model_type,
                                           const Param & params)
  {
    // if the transformation is the identity, don't fit another model:
    if (model_type_ == "identity")
      return;

    delete model_;
    model_ = 0;     // avoid segmentation fault in case of exception
    if ((model_type == "none") || (model_type == "identity"))
    {
      model_ = new TransformationModel();
    }
    else if (model_type == "linear")
    {
      model_ = new TransformationModelLinear(data_, params);
      // // debug output:
      // DoubleReal slope, intercept;
      // TransformationModelLinear* lm = dynamic_cast<TransformationModelLinear*>(model_);
      // lm->getParameters(slope, intercept);
      // cout << "slope: " << slope << ", intercept: " << intercept << endl;
    }
    else if (model_type == "interpolated")
    {
      model_ = new TransformationModelInterpolated(data_, params);
    }
    else if (model_type == "b_spline")
    {
      model_ = new TransformationModelBSpline(data_, params);
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "unknown model type '" + model_type + "'");
    }
    model_type_ = model_type;
  }

  DoubleReal TransformationDescription::apply(DoubleReal value) const
  {
    return model_->evaluate(value);
  }

  const String & TransformationDescription::getModelType() const
  {
    return model_type_;
  }

  void TransformationDescription::getModelTypes(StringList & result)
  {
    result = StringList::create("linear,b_spline,interpolated");
    // "none" and "identity" don't count
  }

  void TransformationDescription::setDataPoints(const DataPoints & data)
  {
    data_ = data;
    model_type_ = "none";     // reset the model even if it was "identity"
    delete model_;
    model_ = new TransformationModel();
  }

  const TransformationDescription::DataPoints &
  TransformationDescription::getDataPoints() const
  {
    return data_;
  }

  void TransformationDescription::getModelParameters(Param & params) const
  {
    model_->getParameters(params);
  }

  void TransformationDescription::invert()
  {
    for (TransformationDescription::DataPoints::iterator it = data_.begin();
         it != data_.end(); ++it)
    {
      *it = make_pair(it->second, it->first);
    }
    // ugly hack for linear model with explicit slope/intercept parameters:
    if ((model_type_ == "linear") && data_.empty())
    {
      TransformationModelLinear * lm =
        dynamic_cast<TransformationModelLinear *>(model_);
      lm->invert();
    }
    else
    {
      Param params;
      getModelParameters(params);
      fitModel(model_type_, params);
    }
  }

} // end of namespace OpenMS
