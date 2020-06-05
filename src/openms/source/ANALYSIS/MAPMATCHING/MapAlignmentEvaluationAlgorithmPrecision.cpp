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
// $Authors: Katharina Albers $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmPrecision.h>

namespace OpenMS
{


  MapAlignmentEvaluationAlgorithmPrecision::MapAlignmentEvaluationAlgorithmPrecision() :
    MapAlignmentEvaluationAlgorithm()
  {
  }

  MapAlignmentEvaluationAlgorithmPrecision::~MapAlignmentEvaluationAlgorithmPrecision()
  {
  }

  void MapAlignmentEvaluationAlgorithmPrecision::evaluate(const ConsensusMap & consensus_map_in, const ConsensusMap & consensus_map_gt, const double & rt_dev, const double & mz_dev, const Peak2D::IntensityType & int_dev, const bool use_charge, double & out)
  {
    //Precision = 1/N * sum ( gt_subtend_tilde_tool_i / tilde_tool_i )

    ConsensusMap cons_map_gt;     /* = consensus_map_gt; */

    for (Size i = 0; i < consensus_map_gt.size(); ++i)
    {
      if (consensus_map_gt[i].size() >= 2)
      {
        cons_map_gt.push_back(consensus_map_gt[i]);
      }
    }

    ConsensusMap cons_map_tool = consensus_map_in;

    std::vector<Size> gt_subtend_tilde_tool;        //holds the numerators of the sum
    std::vector<Size> tilde_tool;               //holds the denominators of the sum

    Size gt_subtend_tilde_tool_i = 0;       //filling material for the vectors
    Size tilde_tool_i = 0;

    Size cons_tool_size = 0;            //size  of the actual consensus feature of the tool
    Size gt_i_subtend_tool_j = 0;       //size of the intersection of the actual cons. feat. of the tool with the c.f. of GT

    double precision = 0;       //holds the output
    double sum = 0;         //intermediate step: the sum

    //loop over all consensus features of the ground truth
    for (Size i = 0; i < cons_map_gt.size(); ++i)      //N = cons_map_gt.size()
    {

      ConsensusFeature & gt_elem = cons_map_gt[i];

      //for every i = 1, ..., N:
      gt_subtend_tilde_tool_i = 0;
      tilde_tool_i = 0;

      //loop over all consensus features of the tool's consensus map
      for (Size j = 0; j < cons_map_tool.size(); ++j)
      {
        ConsensusFeature & tool_elem = cons_map_tool[j];
        cons_tool_size = cons_map_tool[j].size();

        gt_i_subtend_tool_j = 0;

        //loop over all features in the ith consensus feature of the gt
        for (HandleIterator gt_it = gt_elem.begin(); gt_it != gt_elem.end(); ++gt_it)
        {
          //loop over all features in the jth consensus feature of the tool's map
          for (HandleIterator tool_it = tool_elem.begin(); tool_it != tool_elem.end(); ++tool_it)
          {
            //++cons_tool_size;

            if (isSameHandle(*tool_it, *gt_it, rt_dev, mz_dev, int_dev, use_charge))
            {
              ++gt_i_subtend_tool_j;
              break;
            }
          }

        }
        if ((cons_tool_size >= 2) && (gt_i_subtend_tool_j > 0))
        {
          gt_subtend_tilde_tool_i += gt_i_subtend_tool_j;
          tilde_tool_i += cons_tool_size;
        }
      }

      gt_subtend_tilde_tool.push_back(gt_subtend_tilde_tool_i);
      tilde_tool.push_back(tilde_tool_i);
    }
    for (Size k = 0; k < gt_subtend_tilde_tool.size(); ++k)
    {
      double fraction = 0;        //intermediate step: the fraction

      if (gt_subtend_tilde_tool[k] != 0)
      {
        fraction = double(gt_subtend_tilde_tool[k]) / double(tilde_tool[k]);
      }
      sum += fraction;
    }
    precision = (1.0 / double(cons_map_gt.size())) * sum;
    out = precision;
  }

} // namespace OpenMS
