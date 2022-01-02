// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MaxLikeliFitter1D.h>

#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>

namespace OpenMS
{

  void MaxLikeliFitter1D::updateMembers_()
  {
    Fitter1D::updateMembers_();
  }

  MaxLikeliFitter1D::QualityType MaxLikeliFitter1D::fitOffset_(std::unique_ptr<InterpolationModel>& model,
                                                               const RawDataArrayType & set,
                                                               const CoordinateType stdev1,
                                                               const CoordinateType stdev2,
                                                               const CoordinateType offset_step) const
  {
    const CoordinateType offset_min = model->getInterpolation().supportMin() - stdev1;
    const CoordinateType offset_max = model->getInterpolation().supportMin() + stdev2;

    CoordinateType offset;
    QualityType correlation;

    //test model with default offset
    std::vector<float> real_data;
    real_data.reserve(set.size());
    std::vector<float> model_data;
    model_data.reserve(set.size());

    for (Size i = 0; i < set.size(); ++i)
    {
      real_data.push_back(set[i].getIntensity());
      model_data.push_back(model->getIntensity(DPosition<1>(set[i].getPosition())));
    }

    CoordinateType max_offset = model->getInterpolation().getOffset();
    QualityType max_correlation = Math::pearsonCorrelationCoefficient(real_data.begin(), real_data.end(), model_data.begin(), model_data.end());

    //test different offsets
    for (offset = offset_min; offset <= offset_max; offset += offset_step)
    {
      // set offset
      model->setOffset(offset);

      // get samples
      model_data.clear();
      for (Size i = 0; i < set.size(); ++i)
      {
        model_data.push_back(model->getIntensity(DPosition<1>(set[i].getPosition())));
      }

      correlation = Math::pearsonCorrelationCoefficient(real_data.begin(), real_data.end(), model_data.begin(), model_data.end());

      if (correlation > max_correlation)
      {
        max_correlation = correlation;
        max_offset = offset;
      }
    }

    model->setOffset(max_offset);

    return max_correlation;
  }

} // namespace OpenMS
