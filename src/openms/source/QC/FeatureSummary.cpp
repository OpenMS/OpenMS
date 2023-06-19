// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Axel Walter $
// $Authors: Axel Walter $
// --------------------------------------------------------------------------


#include <OpenMS/QC/FeatureSummary.h>

using namespace std;

namespace OpenMS
{
  FeatureSummary::Result FeatureSummary::compute(const FeatureMap& feature_map)
  {
    FeatureSummary::Result result;
    float sum_rt_deviations = 0;
    UInt rt_count = 0;
    result.feature_count = feature_map.size();
    for (const auto& f : feature_map)
    {
      if (f.metaValueExists("rt_deviation"))
      {
        sum_rt_deviations += (float)f.getMetaValue("rt_deviation");
        rt_count += 1;
      }
    }

    // calculate mean rt shift (sec)
    if (rt_count != 0)
    {
      result.rt_shift_mean = sum_rt_deviations / rt_count;
    }
    else
    {
      result.rt_shift_mean = 0;
    }

    return result;
  }

  bool FeatureSummary::Result::operator==(const Result& rhs) const
  {
    return feature_count == rhs.feature_count && rt_shift_mean == rhs.rt_shift_mean;
  }

  /// Returns the name of the metric
  const String& FeatureSummary::getName() const
  {
    return name_;
  }

  /// Returns required file input i.e. MzML.
  /// This is encoded as a bit in a Status object.
  QCBase::Status FeatureSummary::requirements() const
  {
    return QCBase::Status(QCBase::Requires::PREFDRFEAT);
  }
} // namespace OpenMS
