// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Katharina Albers $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmRecall.h>

namespace OpenMS
{


  MapAlignmentEvaluationAlgorithmRecall::MapAlignmentEvaluationAlgorithmRecall() :
    MapAlignmentEvaluationAlgorithm()
  {
  }

  MapAlignmentEvaluationAlgorithmRecall::~MapAlignmentEvaluationAlgorithmRecall() = default;

  void MapAlignmentEvaluationAlgorithmRecall::evaluate(const ConsensusMap & consensus_map_in, const ConsensusMap & consensus_map_gt, const double & rt_dev, const double & mz_dev, const Peak2D::IntensityType & int_dev, const bool use_charge, double & out)
  {
    //Recall = 1/N * sum( gt_subtend_tilde_tool_i / ( m_i * gt_i ) )

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
    std::vector<Size> m;                //holds the denominators of the sum
    std::vector<Size> gt;               //holds the denominators of the sum

    //loop over all consensus features of the ground truth
    for (Size i = 0; i < cons_map_gt.size(); ++i)     //N = cons_map_gt.size()
    {

      ConsensusFeature & gt_elem = cons_map_gt[i];

      //for every i = 1, ..., N:
      Size gt_subtend_tilde_tool_i = 0; // holds the numerators of the sum
      Size m_i = 0;
      Size gt_i = 0; // size  of the actual consensus feature of the GT

      //loop over all consensus features of the tool's consensus map
      for (Size j = 0; j < cons_map_tool.size(); ++j)
      {
        ConsensusFeature & tool_elem = cons_map_tool[j];
        Size gt_i_subtend_tool_j = 0;                    //size of the intersection of the actual cons. feat. of the tool with the c.f. of GT
        Size cons_tool_size = cons_map_tool[j].size();   //size  of the actual consensus feature of the tool_subtend_tilde_tool

        //loop over all features in the ith consensus feature of the gt
        for (HandleIterator gt_it = gt_elem.begin(); gt_it != gt_elem.end(); ++gt_it)
        {
          ++gt_i;

          //loop over all features in the jth consensus feature of the tool's map
          for (HandleIterator tool_it = tool_elem.begin(); tool_it != tool_elem.end(); ++tool_it)
          {
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
          ++m_i;
        }
      }

      gt_subtend_tilde_tool.push_back(gt_subtend_tilde_tool_i);
      m.push_back(m_i);
      gt.push_back(gt_i / cons_map_tool.size());

    }

    double sum = 0;    // intermediate step: the sum
    for (Size k = 0; k < gt_subtend_tilde_tool.size(); ++k)
    {
      double fraction = 0;

      if (gt_subtend_tilde_tool[k] != 0)
      {
        fraction = double(gt_subtend_tilde_tool[k]) / (m[k] * gt[k]);
      }
      sum += fraction;
    }
    double recall = (1.0 / double(cons_map_gt.size())) * sum;
    out = recall;
  }

} // namespace OpenMS
