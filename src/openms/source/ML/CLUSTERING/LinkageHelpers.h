// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/BinaryTreeNode.h>
#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>

#include <set>
#include <utility>
#include <vector>

namespace OpenMS
{
  namespace Internal
  {
    /**
      @brief Private helpers shared by the hierarchical cluster methods AverageLinkage and CompleteLinkage. Not installed.
    */
    namespace LinkageHelpers
    {
      /**
        @brief Agglomerative hierarchical clustering with a Lance-Williams update of the distance matrix

        Repeatedly merges the two clusters i,j of minimal distance until that distance reaches @p threshold or
        only two clusters remain; the tree is then filled with dummy nodes. The clustering is done on the indices
        pointing to elements in the input matrix.

        @param[in,out] original_distance The distances; MUST have >= 2 elements (reduced during clustering)
        @param[out] cluster_tree The resulting tree
        @param[in] threshold Distance at which merging stops
        @param[in] logger Used to report the progress
        @param[in] make_update make_update(m_i, m_j) is called before the clusters i,j (with m_i, m_j elements) are merged
                   and returns a function (d(i,k), d(j,k)) -> d((i,j),k), the distance of the merged cluster to a cluster k
      */
      template <typename MakeUpdate>
      void lanceWilliamsClustering(DistanceMatrix<float>& original_distance, std::vector<BinaryTreeNode>& cluster_tree, const float threshold,
                                   const ProgressLogger& logger, MakeUpdate&& make_update)
      {
        std::vector<std::set<Size> > clusters(original_distance.dimensionsize());
        for (Size i = 0; i < original_distance.dimensionsize(); ++i)
        {
          clusters[i].insert(i);
        }

        cluster_tree.clear();
        cluster_tree.reserve(original_distance.dimensionsize() - 1);

        // Initial minimum-distance pair
        original_distance.updateMinElement();
        std::pair<Size, Size> min = original_distance.getMinElementCoordinates();

        Size overall_cluster_steps(original_distance.dimensionsize());
        logger.startProgress(0, original_distance.dimensionsize(), "clustering data");

        while (original_distance(min.first, min.second) < threshold)
        {
          //grow the tree
          cluster_tree.emplace_back(*(clusters[min.second].begin()), *(clusters[min.first].begin()), original_distance(min.first, min.second));
          if (cluster_tree.back().left_child > cluster_tree.back().right_child)
          {
            std::swap(cluster_tree.back().left_child, cluster_tree.back().right_child);
          }

          if (original_distance.dimensionsize() > 2)
          {
            //pick minimum-distance pair i,j and merge them

            //lance-williams update for d((i,j),k)
            const auto update = make_update(clusters[min.first].size(), clusters[min.second].size());

            //pushback elements of second to first (and then erase second)
            clusters[min.second].insert(clusters[min.first].begin(), clusters[min.first].end());
            // erase first one
            clusters.erase(clusters.begin() + min.first);

            //update original_distance matrix
            for (Size k = 0; k < min.second; ++k)
            {
              float dik = original_distance.getValue(min.first, k);
              float djk = original_distance.getValue(min.second, k);
              original_distance.setValueQuick(min.second, k, update(dik, djk));
            }
            for (Size k = min.second + 1; k < original_distance.dimensionsize(); ++k)
            {
              float dik = original_distance.getValue(min.first, k);
              float djk = original_distance.getValue(min.second, k);
              original_distance.setValueQuick(k, min.second, update(dik, djk));
            }

            //reduce
            original_distance.reduce(min.first);

            //update minimum-distance pair
            original_distance.updateMinElement();

            //get new min-pair
            min = original_distance.getMinElementCoordinates();
          }
          else
          {
            break;
          }
          logger.setProgress(overall_cluster_steps - original_distance.dimensionsize());

          //repeat until only two cluster remains or threshold exceeded, last step skips matrix operations
        }
        //fill tree with dummy nodes
        Size sad(*clusters.front().begin());
        for (Size i = 1; i < clusters.size() && (cluster_tree.size() < cluster_tree.capacity()); ++i)
        {
          cluster_tree.emplace_back(sad, *clusters[i].begin(), -1.0);
        }

        logger.endProgress();
      }
    } // namespace LinkageHelpers
  } // namespace Internal
} // namespace OpenMS
