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

#include <OpenMS/ANALYSIS/MAPMATCHING/SimplePairFinder.h>

namespace OpenMS
{
  SimplePairFinder::SimplePairFinder() :
    Base()
  {
    //set the name for DefaultParamHandler error messages
    setName(getProductName());

    defaults_.setValue("similarity:diff_intercept:RT", 1.0, "This parameter controls the asymptotic decay rate for large differences (for more details see the similarity measurement).", ListUtils::create<String>("advanced"));
    defaults_.setValue("similarity:diff_intercept:MZ", 0.1, "This parameter controls the asymptotic decay rate for large differences (for more details see the similarity measurement).", ListUtils::create<String>("advanced"));
    defaults_.setValue("similarity:diff_exponent:RT", 2.0, "This parameter is important for small differences (for more details see the similarity measurement).", ListUtils::create<String>("advanced"));
    defaults_.setValue("similarity:diff_exponent:MZ", 1.0, "This parameter is important for small differences (for more details see the similarity measurement).", ListUtils::create<String>("advanced"));
    defaults_.setValue("similarity:pair_min_quality", 0.01, "Minimum required pair quality.", ListUtils::create<String>("advanced"));

    Base::defaultsToParam_();
  }

  void SimplePairFinder::run(const std::vector<ConsensusMap> & input_maps, ConsensusMap & result_map)
  {
    if (input_maps.size() != 2)
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "exactly two input maps required");
    checkIds_(input_maps);

    // progress dots
    Int progress_dots = 0;
    if (this->param_.exists("debug::progress_dots"))
    {
      progress_dots = (Int) this->param_.getValue("debug:progress_dots");
    }
    Int number_of_considered_element_pairs = 0;

    // For each element in map 0, find its best friend in map 1
    std::vector<UInt> best_companion_index_0(input_maps[0].size(), UInt(-1));
    std::vector<double> best_companion_quality_0(input_maps[0].size(), 0);
    for (UInt fi0 = 0; fi0 < input_maps[0].size(); ++fi0)
    {
      double best_quality = -std::numeric_limits<double>::max();
      for (UInt fi1 = 0; fi1 < input_maps[1].size(); ++fi1)
      {
        double quality = similarity_(input_maps[0][fi0], input_maps[1][fi1]);
        if (quality > best_quality)
        {
          best_quality = quality;
          best_companion_index_0[fi0] = fi1;
        }

        ++number_of_considered_element_pairs;
        if (progress_dots && !(number_of_considered_element_pairs % progress_dots))
        {
          std::cout << '-' << std::flush;
        }

      }
      best_companion_quality_0[fi0] = best_quality;
    }

    // For each element in map 1, find its best friend in map 0
    std::vector<UInt> best_companion_index_1(input_maps[1].size(), UInt(-1));
    std::vector<double> best_companion_quality_1(input_maps[1].size(), 0);
    for (UInt fi1 = 0; fi1 < input_maps[1].size(); ++fi1)
    {
      double best_quality = -std::numeric_limits<double>::max();
      for (UInt fi0 = 0; fi0 < input_maps[0].size(); ++fi0)
      {
        double quality = similarity_(input_maps[0][fi0], input_maps[1][fi1]);
        if (quality > best_quality)
        {
          best_quality = quality;
          best_companion_index_1[fi1] = fi0;
        }

        ++number_of_considered_element_pairs;
        if (progress_dots && !(number_of_considered_element_pairs % progress_dots))
        {
          std::cout << '+' << std::flush;
        }

      }
      best_companion_quality_1[fi1] = best_quality;
    }

    // And if both like each other, they become a pair.
    // element_pairs_->clear();
    for (UInt fi0 = 0; fi0 < input_maps[0].size(); ++fi0)
    {
      // fi0 likes someone ...
      if (best_companion_quality_0[fi0] > pair_min_quality_)
      {
        // ... who likes him too ...
        UInt best_companion_of_fi0 = best_companion_index_0[fi0];
        if (best_companion_index_1[best_companion_of_fi0] == fi0 &&
            best_companion_quality_1[best_companion_of_fi0] > pair_min_quality_
            )
        {
          ConsensusFeature f;
          f.insert(input_maps[0][fi0]);
          f.insert(input_maps[1][best_companion_of_fi0]);
          f.computeConsensus();
          f.setQuality(best_companion_quality_0[fi0] + best_companion_quality_1[best_companion_of_fi0]);
          result_map.push_back(f);
        }
      }
    }
    return;
  }

  void SimplePairFinder::updateMembers_()
  {
    diff_intercept_[Peak2D::RT] = (double)param_.getValue("similarity:diff_intercept:RT");
    if (diff_intercept_[Peak2D::RT] <= 0)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "intercept for RT must be > 0");
    }

    diff_intercept_[Peak2D::MZ] = (double)param_.getValue("similarity:diff_intercept:MZ");
    if (diff_intercept_[Peak2D::MZ] <= 0)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "intercept for MZ must be > 0");
    }

    diff_exponent_[Peak2D::RT] = (double)param_.getValue("similarity:diff_exponent:RT");
    diff_exponent_[Peak2D::MZ] = (double)param_.getValue("similarity:diff_exponent:MZ");
    pair_min_quality_ = (double)param_.getValue("similarity:pair_min_quality");
  }

  double SimplePairFinder::similarity_(ConsensusFeature const & left, ConsensusFeature const & right) const
  {
    double right_intensity(right.getIntensity());
    if (right_intensity == 0)
      return 0;

    double intensity_ratio = left.getIntensity() / right_intensity;
    if (intensity_ratio > 1.)
      intensity_ratio = 1. / intensity_ratio;

    // if the right map is the transformed map, take the transformed right position
    DPosition<2> position_difference = left.getPosition() - right.getPosition();

    for (UInt dimension = 0; dimension < 2; ++dimension)
    {
      // the formula is explained in class doc
      if (position_difference[dimension] < 0)
      {
        position_difference[dimension] = -position_difference[dimension];
      }
      position_difference[dimension] *= diff_intercept_[dimension];
      position_difference[dimension] += 1.0;
      position_difference[dimension] = pow(position_difference[dimension], diff_exponent_[dimension]);
    }

    return intensity_ratio / position_difference[Peak2D::RT] / position_difference[Peak2D::MZ];
  }

}
