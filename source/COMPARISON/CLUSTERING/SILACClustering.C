// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:expandtab
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2011 -- Bastian Blank
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
// $Maintainer: Lars Nilse $
// $Authors: Bastian Blank $
// --------------------------------------------------------------------------

#include <OpenMS/COMPARISON/CLUSTERING/SILACClustering.h>

using namespace OpenMS;

void SILACClustering::cluster()
{
  HierarchicalClustering<SILACPattern *>::cluster();
  if (rt_min) removeSmall_();
  if (rt_max_spacing) joinLarge_();
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
  typedef std::multimap<DoubleReal, std::pair<PointCoordinate, Cluster> > Sort;
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

  typedef std::vector<std::pair<DoubleReal, Sort> > SortRT;
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
      Cluster &cluster_sort = sort_it->second.second;

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
