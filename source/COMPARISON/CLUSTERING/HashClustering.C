// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Steffen Sass $
// $Authors: $
// --------------------------------------------------------------------------


#include <OpenMS/COMPARISON/CLUSTERING/HashClustering.h>
#include <map>
#include <list>
#include <limits>
#include <cmath>
#include <algorithm>

namespace OpenMS
{

  typedef std::map<std::pair<int,int>, std::list<GridElement*> > GridElements;

  HashClustering::HashClustering(std::vector<DataPoint>& data, DoubleReal rt_threshold, DoubleReal mz_threshold, ClusteringMethod& method_) : grid(HashGrid(rt_threshold,mz_threshold))
  {
    if(data.size() < 2)
    {
      throw InsufficientInput(__FILE__, __LINE__, __PRETTY_FUNCTION__, "The data set contains not enough elements");
    }

    method=&method_;

    for (std::vector<DataPoint>::iterator it = data.begin(); it != data.end(); ++it)
    {
      grid.insert(new DataSubset(*it));
    }
    init();
  }

  DoubleReal HashClustering::getDistance(DataSubset& subset1, DataSubset& subset2)
  {
    return method->getDistance(subset1, subset2);
  }

  DoubleReal HashClustering::getDistance(DataPoint& point1, DataPoint& point2)
  {
    return method->getDistance(point1, point2);
  }

  // Calculate initial distances
  void HashClustering::init()
  {
    min_distance=std::numeric_limits<DoubleReal>::max();
    min_distance_subsets=std::pair<DataSubset*, DataSubset*>(0, 0);
    // Iterate over all cells in the grid
    for (GridElements::iterator it = grid.begin();it != grid.end(); ++it)
    {
      std::pair<int,int> act_coords = it->first;
      std::list<GridElement*>& elements = it->second;
      int x=act_coords.first;
      int y=act_coords.second;

      // Iterate over the upper right cells and calculate the distances between the elements of the current cell and this five neighbor cells
      for (int i = x - 1; i <= x + 1; ++i)
      {
        if (i < 0 || i > grid.getGridSizeX())
        {
          continue;
        }

        for (int j=y;j<=y+1;++j)
        {
          if (j < 0 || j > grid.getGridSizeY() || (i == x - 1 && j == y))
          {
            continue;
          }

          GridElements::iterator act_pos = grid.find(std::make_pair(i, j));

          if (act_pos == grid.end())
          {
            continue;
          }

          std::list<GridElement*>& neighbor_elements = act_pos->second;

          for (std::list<GridElement*>::iterator current_element = elements.begin(); current_element != elements.end(); ++current_element)
          {
            for (std::list<GridElement*>::iterator neighbor_element = neighbor_elements.begin(); neighbor_element != neighbor_elements.end(); ++neighbor_element)
            {
              if ((i == x && j == y && (*current_element)->getID()>(*neighbor_element)->getID()) || *current_element == *neighbor_element)
              {
                continue;
              }

              DataSubset* element_ptr = dynamic_cast<DataSubset*> (*current_element);
              DataSubset* neighbor_ptr = dynamic_cast<DataSubset*> (*neighbor_element);
              DoubleReal act_distance=getDistance(*element_ptr, *neighbor_ptr);

              // Store distances in distance map
              std::pair<DistanceSet::iterator,bool> position = distances.insert(DistanceEntry(element_ptr, neighbor_ptr,act_distance));

              if (position.second)
              {
                // If distance was inserted, insert a pointer to this distance in the current DataSubset
                element_ptr->distance_iterators.insert(std::make_pair(neighbor_ptr, position.first));
                // Remember minimal distance and corresponding DataSubsets
                if (act_distance < min_distance)
                {
                  min_distance = act_distance;
                  min_distance_subsets=std::make_pair(element_ptr, neighbor_ptr);
                }
              }
            }
          }
        }
      }
    }
  }

  typedef std::map<GridElement*,DistanceSet::iterator> IteratorMap;
  typedef std::map<std::pair<int, int>, std::list<GridElement*> > ElementMap;

  void HashClustering::merge()
  {
    // Merge the two subtrees of the minimal distance DataSubset, append a new node and insert it to the first DataSubset

    // check if one of the two pointer is a NULL pointer
    if (min_distance_subsets.first == 0 || min_distance_subsets.second == 0)
    {
      return;
    }

    // Determine that the id of the first element is always smaller than the id of the second element
    if (min_distance_subsets.first->data_points.front() > min_distance_subsets.second->data_points.front())
    {
      std::swap(min_distance_subsets.first, min_distance_subsets.second);
    }

    std::vector<BinaryTreeNode>& tree1 = min_distance_subsets.first->tree;
    std::vector<BinaryTreeNode>& tree2 = min_distance_subsets.second->tree;
    tree1.insert(tree1.end(), tree2.begin(), tree2.end());
    BinaryTreeNode act_node = BinaryTreeNode(min_distance_subsets.first->data_points.front(), min_distance_subsets.second->data_points.front(), min_distance);

    // Append the new node
    tree1.insert(tree1.end(), 1, act_node);
    DataSubset& subset1=*min_distance_subsets.first;
    DataSubset& subset2=*min_distance_subsets.second;

    // Remember old position
    int x = subset1.mz / grid.getMZThreshold();
    int y = subset1.rt / grid.getRTThreshold();

    // Calculate new centroid
    DoubleReal mz_centroid = (subset1.mz * subset1.size()) + (subset2.mz * subset2.size());
    mz_centroid /= subset1.size() + subset2.size();
    DoubleReal rt_centroid = (subset1.rt * subset1.size()) + (subset2.rt * subset2.size());
    rt_centroid /= subset1.size() + subset2.size();
    subset1.rt = rt_centroid;
    subset1.mz = mz_centroid;

    // Insert all data points from subset2 to subset1
    subset1.data_points.insert(subset1.data_points.end(), subset2.data_points.begin(), subset2.data_points.end());
    grid.removeElement(&subset2);

    // Calculate new position
    int x_new = subset1.mz / grid.getMZThreshold();
    int y_new = subset1.rt / grid.getRTThreshold();

    // Iterate over all distances from subset1 and recalculate them
    for (IteratorMap::iterator dist_it = subset1.distance_iterators.begin(); dist_it != subset1.distance_iterators.end();)
    {
      DistanceEntry entry =* (dist_it->second);
      int n_x = entry.data_point->mz / grid.getMZThreshold();
      int n_y = entry.data_point->rt / grid.getRTThreshold();

      DoubleReal new_distance = getDistance(subset1, *(entry.data_point));

      // Check, if the distance is valid, i.e. the DataSubset, to which the distance points, lies in one of the four upper right cells and the distance does not point from subset1 to itself or subset2. Delete it instead.
      if (abs(x_new - n_x) >= 2 || abs(y_new - n_y) >= 2 || (y_new > n_y) || (y_new == n_y && x_new == n_x + 1) || (x_new == n_x && y_new == n_y) || entry.data_point->getID() == entry.owner->getID() || entry.data_point == &subset2)
      {
        distances.erase(dist_it->second);
        subset1.distance_iterators.erase(dist_it++);
      }
      else
      {
        entry.distance = new_distance;
        distances.replace(dist_it->second,entry);
        ++dist_it;
      }
    }
    // Iterate over all distances from subset2. Recalculate and insert them, if they are not already contained in subset1
    for (IteratorMap::iterator dist_it = subset2.distance_iterators.begin(); dist_it != subset2.distance_iterators.end(); ++dist_it)
    {
      DistanceEntry entry = *(dist_it->second);
      int n_x = entry.data_point->mz /  grid.getMZThreshold();
      int n_y = entry.data_point->rt /  grid.getRTThreshold();

      DoubleReal new_distance=getDistance(subset1, *(entry.data_point));

      // Check, if the distance is valid, i.e. the DataSubset, to which the distance points, lies in one of the four upper right cells and the distance does not point from subset2 to itself or subset1. Delete it instead.
      if (abs(x_new - n_x) >= 2 || abs(y_new - n_y) >= 2 || (y_new > n_y) || (y_new == n_y && x_new == n_x + 1) || (x_new == n_x && y_new == n_y) || entry.data_point->getID() == entry.owner->getID() || entry.data_point == &subset1)
      {
        distances.erase(dist_it->second);
      }      
      else
      {
        entry.owner = &subset1;
        entry.distance = new_distance;
        bool replaced = distances.replace(dist_it->second,entry);

        if (replaced)
        {
          subset1.distance_iterators.insert(make_pair(entry.data_point,dist_it->second));
        }
        else
        {
          distances.erase(dist_it->second);
        }
      }
    }
    // Make a set of all coordinates of the five lower left cells (including the current cell itself) of subset1 centroid, subset2 centroid and new calculated centroid
    std::set<std::pair<int, int> > traverse_set;
    if (x != x_new || y != y_new)
    {
      grid.removeElement(&subset1, x, y);
      grid.insert(&subset1);
      for (int i = x_new - 1; i <= x_new + 1; ++i)
      {
        if (i < 0 || i > grid.getGridSizeX())
          continue;
        for (int j = y_new - 1; j <= y_new; ++j)
        {
          if (j < 0 || j > grid.getGridSizeY() || (i == x_new + 1 && j == y_new))
            continue;
          traverse_set.insert(std::make_pair(i, j));
        }
      }
    }
    for (int i = x - 1; i <= x + 1; ++i)
    {
      if (i < 0 || i > grid.getGridSizeX())
        continue;

      for (int j = y - 1; j <= y; ++j)
      {
        if (j < 0 || j > grid.getGridSizeY() || (i == x + 1 && j == y))
          continue;
        traverse_set.insert(std::make_pair(i, j));
      }
    }
    x = subset2.mz / grid.getMZThreshold();
    y = subset2.rt / grid.getRTThreshold();
    for (int i = x - 1; i <= x + 1; ++i)
    {
      if (i < 0 || i > grid.getGridSizeX())
        continue;
      for (int j = y - 1; j <= y; ++j)
      {
        if (j < 0 || j > grid.getGridSizeY() || (i == x + 1 && j == y))
          continue;
        traverse_set.insert(std::make_pair(i, j));
      }
    }

    // Iterate over all possible neighbor cells and insert/delete/recalculate distances
    for (std::set<std::pair<int, int> >::iterator set_it=traverse_set.begin();set_it!=traverse_set.end();++set_it)
    {
      x = set_it->first;
      y = set_it->second;
      ElementMap::iterator act_pos = grid.find(std::make_pair(x, y));
      if (act_pos == grid.end())
        continue;
      std::list<GridElement*>& neighbor_elements=act_pos->second;


      for (std::list<GridElement*>::iterator lit = neighbor_elements.begin(); lit != neighbor_elements.end(); ++lit)
      {
        if (*lit == &subset1)
          continue;
        DataSubset* neighbor_ptr = dynamic_cast<DataSubset*> (*lit);
        DoubleReal act_distance = getDistance(subset1, *neighbor_ptr);
        IteratorMap::iterator pos = neighbor_ptr->distance_iterators.find(&subset2);
        if (pos!=neighbor_ptr->distance_iterators.end())
        {
          distances.erase(pos->second);
          neighbor_ptr->distance_iterators.erase(pos);
        }
        pos=neighbor_ptr->distance_iterators.find(&subset1);

        // Check, if the distance is valid, i.e. the DataSubset, to which the distance points(in this case subset1), lies in one of the four upper right cells of the neighbor element and the distance does not point from subset1 to itself. Delete it instead.
        if (std::abs(x_new - x) >= 2 || std::abs(y_new - y) >= 2 || (y > y_new) || (y == y_new && x == x_new + 1) || act_distance <= 0)
        {
          if (pos != neighbor_ptr->distance_iterators.end())
          {
            distances.erase(pos->second);
            neighbor_ptr->distance_iterators.erase(pos);
          }
        }

        // If it is valid, insert it to the distance set
        else
        {
          if (pos != neighbor_ptr->distance_iterators.end())
          {
            distances.replace(pos->second, DistanceEntry(neighbor_ptr, &subset1, act_distance));
          }
          else
          {
            std::pair<DistanceSet::iterator, bool> new_pos = distances.insert(DistanceEntry(neighbor_ptr, &subset1, act_distance));
            neighbor_ptr->distance_iterators.insert(std::make_pair(&subset1, new_pos.first));
          }
        }
      }
    }
    delete(min_distance_subsets.second);
  }


  void HashClustering::updateMinElements()
  {
    DistanceSet::index<Dist>::type& dist_elements = distances.get<Dist>();

    // Since the distance set is sorted by the distances, it is sufficient to take the first element as the minimal distance element
    DistanceSet::index<Dist>::type::iterator dist_pos = dist_elements.begin();
    if (dist_pos!=dist_elements.end())
    {
      DistanceEntry entry =* dist_pos;
      min_distance = entry.distance;
      min_distance_subsets=std::make_pair(entry.owner, entry.data_point);
    }
  }

  void HashClustering::performClustering()
  {
    do
    {
      // Merge the two subsets
      merge();

      // Find the two new minimal distance DataSubsets
      updateMinElements();
    }
    while(distances.size() > 0);
  }


  std::vector< Real > HashClustering::averageSilhouetteWidth(DataSubset& subset)
  {
    std::vector<BinaryTreeNode>& tree=subset.tree;

    std::vector< Real > average_silhouette_widths; //for each step from the average silhouette widths of the clusters
    std::map<DataPoint*, Real > interdist_i; //for each element i holds the min. average intercluster distance in cluster containing i
    std::map<DataPoint*, DataPoint* > cluster_with_interdist; //for each element i holds which cluster originated the min. intercluster distance
    std::map<DataPoint*, Real > intradist_i; //for each element i holds the average intracluster distance in [i]

    // Initial leafs
    std::set<DataPoint*> leafs;
    for (std::vector<BinaryTreeNode>::iterator it=tree.begin(); it!=tree.end(); ++it)
    {
      leafs.insert(it->data1);
      leafs.insert(it->data2);
      interdist_i.insert(std::make_pair( it->data1, std::numeric_limits<Real>::max() ));
      interdist_i.insert(std::make_pair( it->data2, std::numeric_limits<Real>::max() ));
      cluster_with_interdist.insert(std::make_pair(it->data1, new DataPoint()));
      cluster_with_interdist.insert(std::make_pair(it->data2, new DataPoint()));
      intradist_i.insert(std::make_pair(it->data1, (Real)0.));
      intradist_i.insert(std::make_pair(it->data2, (Real)0.));
      if(it->distance == -1)
      {
        break;
      }
    }

    // Inital values for interdis_i and cluster_with_interdist
    std::set<DataPoint*>::iterator it = leafs.begin();
    ++it;
    for (; it != leafs.end(); ++it)
    {
      std::set<DataPoint*>::iterator jt = leafs.begin();
      for ( ; *jt !=*it; ++jt)
      {
        if(getDistance(**it, **jt)<interdist_i[*it])
        {
          interdist_i[*it] = getDistance(**it, **jt);
          cluster_with_interdist[*it] = *jt;
        }
        if(getDistance(**it, **jt)<interdist_i[*jt])
        {
          interdist_i[*jt] = getDistance(**it, **jt);
          cluster_with_interdist[*jt] = *it;
        }
      }
    }

    // Initial cluster state
    std::map<DataPoint*, std::vector < DataPoint* > > clusters;


    for (std::set<DataPoint*>::iterator it = leafs.begin(); it != leafs.end(); ++it)
    {
      std::vector < DataPoint* > v;
      v.push_back(*it);
      clusters.insert(make_pair(*it, v));

    }

    // Subsequent cluster states after silhouette calc
    for (std::vector<BinaryTreeNode>::iterator tree_it = tree.begin(); tree_it != tree.end(); ++tree_it)
    {

      if (*tree_it == tree.back()) //last steps silhouettes would be all 0 respectively not defined
        break;
      for (std::set<DataPoint*>::iterator it = leafs.begin(); it != leafs.end(); ++it)
      {
        std::vector<DataPoint*>::iterator in_left = std::find(clusters[tree_it->data1].begin(), clusters[tree_it->data1].end(), *it);
        std::vector<DataPoint*>::iterator in_right = std::find(clusters[tree_it->data2].begin(), clusters[tree_it->data2].end(), *it);


        if(in_left == clusters[tree_it->data1].end() && in_right == clusters[tree_it->data2].end()) //*it (!element_of) left or right
        {
          // Intradist_i is always kept
          // Handle interdist:
          if(tree_it->data1 != cluster_with_interdist[*it] && tree_it->data2 != cluster_with_interdist[*it]) //s(i)_nr (!element_of) left or right
          {
            Real interdist_merged = 0;
            for (unsigned int j = 0; j < clusters[tree_it->data1].size(); ++j)
            {
              interdist_merged += getDistance(**it, *clusters[tree_it->data1][j]);
            }
            for (unsigned int j = 0; j < clusters[tree_it->data2].size(); ++j)
            {
              interdist_merged += getDistance(**it, *clusters[tree_it->data2][j]);
            }
            interdist_merged /= (Real)(clusters[tree_it->data1].size() + clusters[tree_it->data2].size());

            if (interdist_merged < interdist_i[*it])
            {
              interdist_i[*it] = interdist_merged;
              cluster_with_interdist[*it] = tree_it->data1;
            }
          }
          else      // s(i)_nr (element_of) left or right
          {
            // Calculate interdist_i to merged
            DataPoint* k; //the one cluster of the two merged which does NOT contain s(i)_nr
            if(tree_it->data2 != cluster_with_interdist[*it] )
            {
              k = tree_it->data2;
            }
            else
            {
              k = tree_it->data1;
            }
            Real interdist_merged = 0;
            for (unsigned int j = 0; j < clusters[k].size(); ++j)
            {
              interdist_merged += getDistance(**it, *clusters[k][j]);
            }
            interdist_merged += (clusters[cluster_with_interdist[*it]].size()*interdist_i[*it]);
            interdist_merged /= (Real)(clusters[k].size()+clusters[cluster_with_interdist[*it]].size());

            // if new inderdist is smaller that old min. nothing else has to be done
            if (interdist_merged <= interdist_i[*it])
            {
              interdist_i[*it] = interdist_merged;
              cluster_with_interdist[*it] = tree_it->data1;
            }

            // else find min av. dist from other clusters to i
            else
            {
              interdist_i[*it] = interdist_merged;
              cluster_with_interdist[*it] = tree_it->data1;

              for (std::map<DataPoint*, std::vector < DataPoint* > >::iterator uit = clusters.begin(); uit != clusters.end(); ++uit)
              {
                if(uit->first != tree_it->data1 && uit->first != tree_it->data2 && !uit->second.empty() && uit->second.end() == std::find(uit->second.begin(), uit->second.end(), *it))
                {
                  Real min_interdist_i = 0;
                  for (unsigned int v = 0; v < uit->second.size(); ++v)
                  {
                    min_interdist_i += getDistance(*uit->second[v], **it);
                  }
                  min_interdist_i /= (Real)uit->second.size();
                  if (min_interdist_i < interdist_i[*it])
                  {
                    interdist_i[*it] = min_interdist_i;
                    cluster_with_interdist[*it] = uit->first;
                  }
                }
              }
            }
          }

        }
        else      // i (element_of) left or right
        {
          DataPoint* k;
          DataPoint* l;     // k is the cluster that is one of the merged but not the one containing i, l the cluster containing i
          if(in_left == clusters[tree_it->data1].end())
          {
            l = tree_it->data2;
            k = tree_it->data1;
          }
          else
          {
            l = tree_it->data1;
            k = tree_it->data2;
          }
          if(k != cluster_with_interdist[*it]) //s(i)_nr (!element_of) left or right cluster
          {
            // Interdist_i is kept
            // But intradist_i has to be updated
            intradist_i[*it] *= clusters[l].size() - 1;
            for (unsigned int j = 0; j < clusters[k].size(); ++j)
            {
              intradist_i[*it] += getDistance(**it, *clusters[k][j]);
            }
            intradist_i[*it] /= (Real)(clusters[k].size()+(clusters[l].size()-1));
          }
          else      //s(i)_nr (element_of) left or right
          {
            // Intradist_i has to be updated
            intradist_i[*it] *= clusters[l].size() - 1;
            intradist_i[*it] += (clusters[k].size() * interdist_i[*it]);
            intradist_i[*it] /= (Real)(clusters[k].size()+(clusters[l].size() - 1));

            // Find new min av. interdist_i
            interdist_i[*it]=std::numeric_limits<Real>::max();
            for (std::map<DataPoint*, std::vector < DataPoint* > >::iterator uit = clusters.begin(); uit != clusters.end(); ++uit)
            {
              if(uit->first!=l && uit->first!=k && !uit->second.empty())
              {
                Real av_interdist_i=0;
                for (unsigned int v = 0; v < uit->second.size(); ++v)
                {
                  av_interdist_i += getDistance(*uit->second[v],**it);
                }
                av_interdist_i /= (Real)uit->second.size();
                if (av_interdist_i < interdist_i[*it])
                {
                  interdist_i[*it]=av_interdist_i;
                  cluster_with_interdist[*it] = uit->first;
                }
              }
            }
          }
        }
      }

      // Redo clustering following tree
      // Pushback elements of data2 to data1 (and then erase second)


      std::map<DataPoint*, std::vector < DataPoint* > >::iterator clusterPos1 = clusters.find(tree_it->data1);
      std::map<DataPoint*, std::vector < DataPoint* > >::iterator clusterPos2 = clusters.find(tree_it->data2);

      std::vector < DataPoint* >& data = clusterPos1->second;

      data.insert(clusterPos1->second.begin(), clusterPos2->second.begin(), clusterPos2->second.end());
      clusterPos2->second.clear();
      // Erase second one

      // Calculate average silhouette width for clusters and then overall average silhouette width for cluster step
      Real average_overall_silhouette = 0; // from cluster step

      for (std::map<DataPoint*, std::vector < DataPoint* > >::iterator git=clusters.begin(); git != clusters.end(); ++git)
      {
        if(git->second.size() > 1)
        {
          // Collect silhouettes clusterwise so that average cluster silhouettes will be easily accessible
          for (unsigned int h = 0; h < git->second.size(); ++h)
          {
            if(interdist_i[git->second[h]] != 0)
            {
              average_overall_silhouette += (interdist_i[git->second[h]] - intradist_i[git->second[h]]) / std::max(interdist_i[git->second[h]], intradist_i[git->second[h]]);

            }
          }
        }
      }
      average_silhouette_widths.push_back(average_overall_silhouette / (Real)(tree.size() + 1));
    }
    average_silhouette_widths.push_back(0.0);
    return average_silhouette_widths;
  }

  bool dataPointPtrCompare(DataPoint* d1, DataPoint* d2)
  {
    return d1->getID() < d2->getID();
  }


  void HashClustering::cut(int cluster_quantity, std::vector< std::vector<DataPoint*> >& clusters, std::vector<BinaryTreeNode>& tree)
  {
    std::set<DataPoint*> leafs;
    for (std::vector<BinaryTreeNode>::iterator it = tree.begin(); it != tree.end(); ++it)
    {
      leafs.insert(it->data1).second;
      leafs.insert(it->data2).second;
      if(it->distance == -1)
      {
        break;
      }
    }

    std::map<DataPoint*, std::vector<DataPoint*> > cluster_map;
    std::set<DataPoint*>::iterator sit = leafs.begin();
    // ++it;
    for ( ; sit != leafs.end(); ++sit)

    {
      cluster_map[*sit] = std::vector<DataPoint*>(1, *sit);
    }
    // Redo clustering till step (original.dimensionsize()-cluster_quantity)
    std::vector<BinaryTreeNode>::iterator it = tree.begin();
    for (unsigned int cluster_step = 0; cluster_step < tree.size() + 1 - cluster_quantity; ++cluster_step)
    {
      // Pushback elements of data2 to data1 (and then erase second)
      cluster_map[it->data1].insert(cluster_map[it->data1].end(), cluster_map[it->data2].begin(), cluster_map[it->data2].end());

      // Erase second one
      cluster_map[it->data2].clear();
      ++it;
    }
    // clusters.reserve(tree.size()+1);
    std::map<DataPoint*,std::vector<DataPoint*> >::iterator iter;
    for( iter = cluster_map.begin(); iter != cluster_map.end(); ++iter )
    {
      std::vector<DataPoint*> actCluster = iter->second;
      clusters.push_back(actCluster);
    }
    // ~ sorts by first element contained!!
    for (std::vector< std::vector<DataPoint*> >::iterator it = clusters.begin();it != clusters.end(); ++it)
    {
      std::sort(it->begin(),it->end(), dataPointPtrCompare);
    }
    std::sort(clusters.begin(), clusters.end());
    std::reverse(clusters.begin(), clusters.end());
    clusters.resize(cluster_quantity);
    std::sort(clusters.begin(), clusters.end());
  }

  void HashClustering::getSubtrees(std::vector<std::vector<BinaryTreeNode> >& subtrees)
  {
    // Extract the subtrees and append them to the subtree vector
    for (ElementMap::iterator it = grid.begin(); it != grid.end(); ++it)
    {
      std::list<GridElement*>& elements = it->second;
      for (std::list<GridElement*>::iterator lit = elements.begin(); lit != elements.end(); ++lit)
      {
        DataSubset* subset_ptr = dynamic_cast<DataSubset*> (*lit);
        if (subset_ptr->size() < 2)
          continue;
        std::vector<BinaryTreeNode> tree;
        tree.insert(tree.begin(), subset_ptr->tree.begin(), subset_ptr->tree.end());
        sort(tree.begin(), tree.end());
        subtrees.push_back(tree);
      }
    }
  }



  typedef std::vector<DataPoint*> Cluster;

  void HashClustering::createClusters(std::vector<Cluster>& clusters)
  {
    Size cluster_id = 0;

    // Run silhoutte optimization for all subtrees and find the appropriate best_n
    for (ElementMap::iterator it = grid.begin(); it != grid.end(); ++it)
    {
      std::list<GridElement*>& elements=it->second;
      for (std::list<GridElement*>::iterator lit = elements.begin(); lit != elements.end(); ++lit)
      {
        DataSubset* subset_ptr = dynamic_cast<DataSubset*> (*lit);
        if (subset_ptr->size() < 2)
        {
          if (subset_ptr->data_points.begin() != subset_ptr->data_points.end())
          {
            Cluster act_cluster;
            DataPoint* act_element = subset_ptr->data_points.front();
            act_element->cluster_id = cluster_id++;
            act_cluster.push_back(act_element);
            clusters.push_back(act_cluster);
          }
          continue;
        }
        sort(subset_ptr->tree.begin(), subset_ptr->tree.end());
        std::vector< Real > asw = averageSilhouetteWidth(*subset_ptr);
        silhouettes.push_back(asw);

        // Look only in the front area of the silhoutte values to avoid getting the wrong number
        std::vector< Real >::iterator max_el(max_element(asw.begin() + 0.9 * asw.size(), asw.end()));
        std::vector<Cluster> act_clusters;
        Size best_n = (Size)subset_ptr->tree.size();
        if (*max_el < 0.75)
        {
          std::vector<DataPoint*> actCluster;
          actCluster.insert(actCluster.begin(), subset_ptr->data_points.begin(), subset_ptr->data_points.end());
          act_clusters.push_back(actCluster);
        }
        else
        {
          for (Size i = 0; i < asw.size(); ++i)
          {
            if(std::fabs(asw[i]-(*max_el)) <= 0)
            {
              best_n = Size(subset_ptr->tree.size() - i);
              break;
            }
          }

          // Get clusters from the subtree
          cut(best_n, act_clusters, subset_ptr->tree);
        }



        // Run through all elements in each cluster and update cluster number
        for (std::vector<Cluster>::iterator cluster_it = act_clusters.begin(); cluster_it != act_clusters.end(); ++cluster_it)
        {
          for (Cluster::iterator element_it = cluster_it->begin(); element_it != cluster_it->end(); ++element_it)
          {
            (*element_it)->cluster_id = cluster_id;
            (*element_it)->cluster_size = cluster_it->size();
          }
          ++cluster_id;
          clusters.push_back(*cluster_it);
        }

      }
    }

  }

  std::vector<std::vector<Real> > HashClustering::getSilhouetteValues()
  {
    return silhouettes;
  }

  HashClustering::InsufficientInput::InsufficientInput(const char* file, int line, const char* function, const char* message) throw()
    : BaseException(file, line, function, "ClusterFunctor::InsufficentInput", message)
  {

  }

  HashClustering::InsufficientInput::~InsufficientInput() throw()
  {

  }
}
