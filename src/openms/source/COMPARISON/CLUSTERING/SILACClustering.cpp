// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Lars Nilse $
// $Authors: Bastian Blank $
// --------------------------------------------------------------------------

#include <OpenMS/COMPARISON/CLUSTERING/SILACClustering.h>

using namespace OpenMS;

void SILACClustering::cluster()
{
  HierarchicalClustering<SILACPattern *>::cluster();
  if (rt_min)
    removeSmall_();
  if (rt_max_spacing)
    joinLarge_();
}

void SILACClustering::removeSmall_()
{
  for (Grid::grid_iterator cell_it = grid.grid_begin();
       cell_it != grid.grid_end();
       ++cell_it)
  {
    Grid::cell_iterator cluster_it = cell_it->second.begin();
    while (cluster_it != cell_it->second.end())
    {
      Grid::cell_iterator cluster_mod_it = cluster_it++;

      // Remove clusters with a size smaller then the given threshold
      if (cluster_mod_it->second.bbox.size()[0] < rt_min)
      {
        cell_it->second.erase(cluster_mod_it);
      }
    }
  }
}

void SILACClustering::joinLarge_()
{
  typedef std::multimap<double, std::pair<PointCoordinate, Cluster> > Sort;
  Sort sort_mz;

  // Sort all input clusters according to MZ
  for (Grid::grid_iterator cell_it = grid.grid_begin();
       cell_it != grid.grid_end();
       ++cell_it)
  {
    for (Grid::cell_iterator cluster_it = cell_it->second.begin();
         cluster_it != cell_it->second.end();
         ++cluster_it)
    {
      sort_mz.insert(std::make_pair(cluster_it->first[1], *cluster_it));
    }
  }

  grid.clear();

  typedef std::vector<std::pair<double, Sort> > SortRT;
  SortRT sort_rt;
  SortRT::iterator sort_rt_check = sort_rt.end();

  // Bucketize clusters according to MZ
  for (Sort::iterator sort_it = sort_mz.begin();
       sort_it != sort_mz.end();
       ++sort_it)
  {
    // Generate new bucket if we got none yet
    if (sort_rt_check == sort_rt.end() ||
        // or if the new point is too far from the first one
        (sort_it->first - sort_rt_check->first) > grid.cell_dimension[1])
    {
      sort_rt_check = sort_rt.insert(sort_rt_check, std::make_pair(sort_it->first, Sort()));
    }

    // Insert new point
    sort_rt_check->second.insert(std::make_pair(sort_it->second.first[0], sort_it->second));
  }

  typedef std::vector<Cluster> Clusters;

  for (SortRT::iterator mz_it = sort_rt.begin();
       mz_it != sort_rt.end();
       ++mz_it)
  {
    Clusters cluster;
    Clusters::iterator cluster_check = cluster.end();

    for (Sort::iterator sort_it = mz_it->second.begin();
         sort_it != mz_it->second.end();
         ++sort_it)
    {
      Cluster & cluster_sort = sort_it->second.second;

      // Generate new output cluster if we got none yet
      if (cluster_check == cluster.end() ||
          // or if the spacing is too high
          (cluster_sort.bbox.first[0] - cluster_check->bbox.second[0]) > rt_max_spacing)
      {
        cluster_check = cluster.insert(cluster_check, cluster_sort);
      }
      // Otherwise merge the clusters
      else
      {
        // Merge the points
        for (Cluster::iterator it = cluster_sort.begin(); it != cluster_sort.end(); ++it)
        {
          cluster_check->insert(*it);
        }
        // Update bounding box
        cluster_check->bbox |= cluster_sort.bbox;
      }
    }

    // Re-add clusters into the grid
    for (Clusters::const_iterator it = cluster.begin(); it != cluster.end(); ++it)
    {
      grid.insert(std::make_pair(it->bbox, *it));
    }
  }
}
