// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Clemens Groepl $
// $Authors: Katharina Albers $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmPrecision.h>

#include <vector>

namespace OpenMS
{


	MapAlignmentEvaluationAlgorithmPrecision::MapAlignmentEvaluationAlgorithmPrecision()
		: MapAlignmentEvaluationAlgorithm()
	{
	}

	MapAlignmentEvaluationAlgorithmPrecision::~MapAlignmentEvaluationAlgorithmPrecision()
	{
	}

	void MapAlignmentEvaluationAlgorithmPrecision::evaluate(const ConsensusMap& consensus_map_in, const ConsensusMap& consensus_map_gt, const DoubleReal& rt_dev, const DoubleReal& mz_dev, const Peak2D::IntensityType& int_dev, const bool use_charge, DoubleReal& out)
	{
		//Precision = 1/N * sum ( gt_subtend_tilde_tool_i / tilde_tool_i )

		ConsensusMap cons_map_gt; /* = consensus_map_gt; */

		for ( Size i = 0; i < consensus_map_gt.size(); ++i )
		{
			if (consensus_map_gt[i].size() >= 2 )
			{
				cons_map_gt.push_back(consensus_map_gt[i]);
			}
		}

		ConsensusMap cons_map_tool = consensus_map_in;

		std::vector<Size> gt_subtend_tilde_tool;	//holds the numerators of the sum
		std::vector<Size> tilde_tool;			//holds the denominators of the sum

		Size gt_subtend_tilde_tool_i = 0;	//filling material for the vectors
		Size tilde_tool_i = 0;

		Size cons_tool_size = 0;		//size  of the actual consensus feature of the tool
		Size gt_i_subtend_tool_j = 0;	//size of the intersection of the actual cons. feat. of the tool with the c.f. of GT

		DoubleReal precision = 0;	//holds the output
		DoubleReal fraction = 0;	//intermediate step: the fraction
		DoubleReal sum = 0;		//intermediate step: the sum

		//loop over all consensus features of the ground truth
		for ( Size i = 0; i < cons_map_gt.size(); ++i) //N = cons_map_gt.size()
		{

			ConsensusFeature& gt_elem = cons_map_gt[i];

			//for every i = 1, ..., N:
			gt_subtend_tilde_tool_i = 0;
			tilde_tool_i = 0;

			//loop over all consensus features of the tool's consensus map
			for (Size j = 0; j < cons_map_tool.size(); ++j)
			{
				ConsensusFeature& tool_elem = cons_map_tool[j];
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
				if((cons_tool_size >= 2) && (gt_i_subtend_tool_j > 0))
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
			fraction = 0;

			if (gt_subtend_tilde_tool[k] != 0)
			{
				fraction = DoubleReal(gt_subtend_tilde_tool[k]) / DoubleReal(tilde_tool[k]);
			}
			sum += fraction;
		}
		precision = (1.0 / DoubleReal(cons_map_gt.size()) ) * sum;
		out = precision;
	}

} // namespace OpenMS

