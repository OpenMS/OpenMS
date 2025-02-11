// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/ML/CLUSTERING/ClusterAnalyzer.h>

#include <OpenMS/DATASTRUCTURES/String.h>

#include <map>
#include <set>
#include <list>

//using namespace std;

namespace OpenMS
{
  ClusterAnalyzer::ClusterAnalyzer() = default;

  ClusterAnalyzer::~ClusterAnalyzer() = default;

  ClusterAnalyzer & ClusterAnalyzer::operator=(const ClusterAnalyzer & source)
  {
    //ALWAYS CHECK FOR SELF ASSIGNMENT!
    if (this == &source)
    {
      return *this;
    }
    //...
    return *this;
  }

  std::vector<float> ClusterAnalyzer::averageSilhouetteWidth(const std::vector<BinaryTreeNode> & tree, const DistanceMatrix<float> & original)
  {
    //throw exception if cannot be legal clustering
    if (tree.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "tree is empty but minimal clustering hirachy has at least one level");
    }

    std::vector<float> average_silhouette_widths;       //for each step from the average silhouette widths of the clusters
    std::map<Size, float> interdist_i;      //for each element i holds the min. average intercluster distance in cluster containing i
    std::map<Size, Size> cluster_with_interdist;      //for each element i holds which cluster originated the min. intercluster distance
    std::map<Size, float> intradist_i;      //for each element i holds the average intracluster distance in [i]

    //initial leafs
    std::set<Size> leafs;
    for (Size i = 0; i < tree.size(); ++i)
    {
      leafs.insert(tree[i].left_child);
      leafs.insert(tree[i].right_child);
      interdist_i.insert(std::make_pair(tree[i].left_child, std::numeric_limits<float>::max()));
      interdist_i.insert(std::make_pair(tree[i].right_child, std::numeric_limits<float>::max()));
      cluster_with_interdist.insert(std::make_pair(tree[i].left_child, 0));
      cluster_with_interdist.insert(std::make_pair(tree[i].right_child, 0));
      intradist_i.insert(std::make_pair(tree[i].left_child, (float)0.));
      intradist_i.insert(std::make_pair(tree[i].right_child, (float)0.));
      if (tree[i].distance == -1)
      {
        break;
      }
    }


    //initial values for interdis_i and cluster_with_interdist
    std::set<Size>::iterator leafs_it = leafs.begin();
    ++leafs_it;
    for (; leafs_it != leafs.end(); ++leafs_it)
    {
      std::set<Size>::iterator jt = leafs.begin();
      for (; *jt < *leafs_it; ++jt)
      {
        if (original.getValue(*leafs_it, *jt) < interdist_i[*leafs_it])
        {
          interdist_i[*leafs_it] = original.getValue(*leafs_it, *jt);
          cluster_with_interdist[*leafs_it] = *jt;
        }
        if (original.getValue(*leafs_it, *jt) < interdist_i[*jt])
        {
          interdist_i[*jt] = original.getValue(*leafs_it, *jt);
          cluster_with_interdist[*jt] = *leafs_it;
        }
      }
    }

    /* to manually retrace
    for (Size i = 0; i < original.dimensionsize(); ++i)
    {
        std::cout << interdist_i[i] << " | " << cluster_with_interdist[i] << " | " << intradist_i[i] << std::endl;
    }
    */

    //initial cluster state
    std::map<Size, std::vector<Size> > clusters;
    for (std::set<Size>::iterator it = leafs.begin(); it != leafs.end(); ++it)
    {
      clusters[*it].push_back(*it);
    }

    //subsequent cluster states after silhouette calc
    for (Size t = 0; t < tree.size() - 1; ++t)   //last steps silhouettes would be all 0 respectively not defined
    {

      for (std::set<Size>::iterator it = leafs.begin(); it != leafs.end(); ++it)
      {

        std::vector<Size>::iterator in_left = std::find(clusters[tree[t].left_child].begin(), clusters[tree[t].left_child].end(), *it);
        std::vector<Size>::iterator in_right = std::find(clusters[tree[t].right_child].begin(), clusters[tree[t].right_child].end(), *it);

        if (in_left == clusters[tree[t].left_child].end() && in_right == clusters[tree[t].right_child].end())    //*it (!element_of) left or right
        {
          //intradist_i is always kept
          //handle interdist:
          if (tree[t].left_child != cluster_with_interdist[*it] && tree[t].right_child != cluster_with_interdist[*it])          //s(i)_nr (!element_of) left or right
          {
            float interdist_merged(0);
            for (Size j = 0; j < clusters[tree[t].left_child].size(); ++j)
            {
              interdist_merged += original.getValue(*it, clusters[tree[t].left_child][j]);
            }
            for (Size j = 0; j < clusters[tree[t].right_child].size(); ++j)
            {
              interdist_merged += original.getValue(*it, clusters[tree[t].right_child][j]);
            }
            interdist_merged /= (float)(clusters[tree[t].left_child].size() + clusters[tree[t].right_child].size());
            if (interdist_merged < interdist_i[*it])
            {
              interdist_i[*it] = interdist_merged;
              cluster_with_interdist[*it] = tree[t].left_child;
            }
          }
          else           //s(i)_nr (element_of) left or right
          {
            //calculate interdist_i to merged
            Size k;             //the one cluster of the two merged which does NOT contain s(i)_nr
            if (tree[t].right_child != cluster_with_interdist[*it])
            {
              k = tree[t].right_child;
            }
            else
            {
              k = tree[t].left_child;
            }
            float interdist_merged(0);
            for (Size j = 0; j < clusters[k].size(); ++j)
            {
              interdist_merged += original.getValue(*it, clusters[k][j]);
            }
            interdist_merged += (clusters[cluster_with_interdist[*it]].size() * interdist_i[*it]);
            interdist_merged /= (float)(clusters[k].size() + clusters[cluster_with_interdist[*it]].size());
            //if new inderdist is smaller that old min. nothing else has to be done
            if (interdist_merged <= interdist_i[*it])
            {
              interdist_i[*it] = interdist_merged;
              cluster_with_interdist[*it] = tree[t].left_child;
            }
            // else find min av. dist from other clusters to i
            else
            {
              interdist_i[*it] = interdist_merged;
              cluster_with_interdist[*it] = tree[t].left_child;

              for (Size u = 0; u < clusters.size(); ++u)
              {
                if (u != tree[t].left_child && u != tree[t].right_child && !clusters[u].empty() && clusters[u].end() == std::find(clusters[u].begin(), clusters[u].end(), *it))
                {
                  float min_interdist_i(0);
                  for (Size v = 0; v < clusters[u].size(); ++v)
                  {
                    min_interdist_i += original.getValue(clusters[u][v], *it);
                  }
                  min_interdist_i /= (float)clusters[u].size();
                  if (min_interdist_i < interdist_i[*it])
                  {
                    interdist_i[*it] = min_interdist_i;
                    cluster_with_interdist[*it] = u;
                  }
                }
              }
            }
          }

        }
        else         //i (element_of) left or right
        {
          Size k, l;          //k is the cluster that is one of the merged but not the one containing i, l the cluster containing i
          if (in_left == clusters[tree[t].left_child].end())
          {
            l = tree[t].right_child;
            k = tree[t].left_child;
          }
          else
          {
            l = tree[t].left_child;
            k = tree[t].right_child;
          }

          if (k != cluster_with_interdist[*it])          //s(i)_nr (!element_of) left or right cluster
          {
            //interdist_i is kept
            //but intradist_i has to be updated
            intradist_i[*it] *= clusters[l].size() - 1;
            for (Size j = 0; j < clusters[k].size(); ++j)
            {
              intradist_i[*it] += original.getValue(*it, clusters[k][j]);
            }
            intradist_i[*it] /= (float)(clusters[k].size() + (clusters[l].size() - 1));
          }
          else           //s(i)_nr (element_of) left or right
          {
            //intradist_i has to be updated
            intradist_i[*it] *= clusters[l].size() - 1;
            intradist_i[*it] += (clusters[k].size() * interdist_i[*it]);
            intradist_i[*it] /= (float)(clusters[k].size() + (clusters[l].size() - 1));
            //find new min av. interdist_i
            interdist_i[*it] = std::numeric_limits<float>::max();
            for (Size u = 0; u < clusters.size(); ++u)
            {
              if (u != l && u != k && !clusters[u].empty())
              {
                float av_interdist_i(0);
                for (Size v = 0; v < clusters[u].size(); ++v)
                {
                  av_interdist_i += original.getValue(clusters[u][v], *it);
                }
                av_interdist_i /= (float)clusters[u].size();
                if (av_interdist_i < interdist_i[*it])
                {
                  interdist_i[*it] = av_interdist_i;
                  cluster_with_interdist[*it] = u;
                }
              }
            }
          }
        }
      }
      //redo clustering following tree
      //pushback elements of right_child to left_child (and then erase second)
      clusters[tree[t].left_child].insert(clusters[tree[t].left_child].end(), clusters[tree[t].right_child].begin(), clusters[tree[t].right_child].end());

      //erase second one
      clusters[tree[t].right_child].clear();

      //~ //adept the cluster indices in clusters with interdist
      //~ for (Size x = 0; x < cluster_with_interdist.size();++x)
      //~ {
      //~ if(cluster_with_interdist[x]>tree[t].right_child)
      //~ {
      //~ --cluster_with_interdist[x];
      //~ }
      //~ }

      /* to manually retrace
      for (Size x = 0; x < clusters.size();++x)
      {
          for (Size y = 0; y < clusters[x].size();++y)
          {
              std::cout << clusters[x][y] << " ";
          }
          std::cout << " | ";
      }
      std::cout << std::endl;
      std::cout << "---------" << std::endl;
      for (Size z = 0; z < original.dimensionsize(); ++z)
      {
          std::cout << interdist_i[z] << " , " << intradist_i[z] << " , " << cluster_with_interdist[z] << " , ";
          std::cout << ((interdist_i[z] - intradist_i[z]) / std::max(interdist_i[z],intradist_i[z])) << std::endl;

      }
      std::cout << "---------" << std::endl;
      */

      //calculate average silhouette width for clusters and then overall average silhouette width for cluster step
      float average_overall_silhouette(0);       // from cluster step
      /* to manually retrace
      std::vector<float> silhouettes(original.dimensionsize(),0.0);
      */
      for (Size g = 0; g < clusters.size(); ++g)
      {
        if (clusters[g].size() > 1)
        {
          //collect silhouettes clusterwise so that average cluster silhouettes will be easily accessible
          for (Size h = 0; h < clusters[g].size(); ++h)
          {
            if (interdist_i[clusters[g][h]] != 0)
            {
              average_overall_silhouette += (interdist_i[clusters[g][h]] - intradist_i[clusters[g][h]]) / std::max(interdist_i[clusters[g][h]], intradist_i[clusters[g][h]]);
              /* to manually retrace
                  silhouettes[clusters[g][h]] = (interdist_i[clusters[g][h]] - intradist_i[clusters[g][h]]) / std::max(interdist_i[clusters[g][h]],intradist_i[clusters[g][h]]);
              */
            }
          }
        }
      }
      /* to manually retrace
          for (Size i = 0; i < silhouettes.size(); ++i)
          {
              std::cout << "s(" <<  (i) << ") = " << silhouettes[i] << std::endl;
          }
          std::cout << "---------" << std::endl;
      */
      average_silhouette_widths.push_back(average_overall_silhouette / (float)(tree.size() + 1));
    }
    average_silhouette_widths.push_back(0.0);
    return average_silhouette_widths;
  }

  std::vector<float> ClusterAnalyzer::dunnIndices(const std::vector<BinaryTreeNode> & tree, const DistanceMatrix<float> & original, const bool tree_from_singlelinkage)
  {
    //throw exception if cannot be legal clustering
    if (tree.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "tree is empty but minimal clustering hirachy has at least one level");
    }

    std::vector<float> all_dunn_indices;
    all_dunn_indices.reserve(tree.size() + 1);

    std::set<Size> leafs;
    for (Size i = 0; i < tree.size(); ++i)
    {
      leafs.insert(tree[i].left_child);
      leafs.insert(tree[i].right_child);
    }

    //initial cluster state
    //~ Size sz = *(leafs.rbegin())+1;
    std::vector<std::vector<Size> > clusters(original.dimensionsize());
    std::vector<std::pair<float, Size> > min_intercluster_distances(original.dimensionsize(), std::make_pair<float, Size>(-1, 0));
    for (std::set<Size>::iterator it = leafs.begin(); it != leafs.end(); ++it)
    {
      clusters[*it].push_back(*it);
      std::set<Size>::iterator it_2 = leafs.begin();
      for (; it_2 != it; ++it_2)
      {
        float d = original.getValue(*it, *it_2);
        if (d < min_intercluster_distances[*it].first || min_intercluster_distances[*it].first == -1)
        {
          min_intercluster_distances[*it].first = d;
          min_intercluster_distances[*it].second = *it_2;
        }
      }
      it_2 = it;
      ++it_2;
      for (; it_2 != leafs.end(); ++it_2)
      {
        float d = original.getValue(*it, *it_2);
        if (d < min_intercluster_distances[*it].first || min_intercluster_distances[*it].first == -1)
        {
          min_intercluster_distances[*it].first = d;
          min_intercluster_distances[*it].second = *it_2;
        }
      }
    }
    Size min_intercluster_distance_index(0);
    for (Size i = min_intercluster_distance_index + 1; i < min_intercluster_distances.size(); ++i)
    {
      if (min_intercluster_distances[min_intercluster_distance_index].first == -1)
      {
        min_intercluster_distance_index = i;
      }
      else if (min_intercluster_distances[i].first != -1 && min_intercluster_distances[i].first < min_intercluster_distances[min_intercluster_distance_index].first)
      {
        min_intercluster_distance_index = i;
      }
    }

    //initial state for min inter and max intra distances
    float max_intracluster_distance(0);
    for (Size cluster_step = 0; cluster_step < tree.size() - 1; ++cluster_step)
    {

      //max intracluster distance changed?
      for (Size x = 0; x < clusters[tree[cluster_step].left_child].size(); ++x)
      {
        for (Size y = 0; y < clusters[tree[cluster_step].right_child].size(); ++y)
        {
          if (original.getValue(clusters[tree[cluster_step].left_child][x], clusters[tree[cluster_step].right_child][y]) > max_intracluster_distance)
          {
            max_intracluster_distance = original.getValue(clusters[tree[cluster_step].left_child][x], clusters[tree[cluster_step].right_child][y]);
          }
        }
      }

      //redo clustering following tree
      //pushback elements of right_child to left_child (and then erase second)
      clusters[tree[cluster_step].left_child].insert(clusters[tree[cluster_step].left_child].end(), clusters[tree[cluster_step].right_child].begin(), clusters[tree[cluster_step].right_child].end());

      //erase second one
      clusters[tree[cluster_step].right_child].clear();

      //min intercluster distance changed?
      if (!tree_from_singlelinkage)
      {
        min_intercluster_distances[tree[cluster_step].right_child].first = -1;
        min_intercluster_distances[tree[cluster_step].right_child].second = 0;

        if ((min_intercluster_distance_index == tree[cluster_step].right_child && min_intercluster_distances[min_intercluster_distance_index].second == tree[cluster_step].left_child)
           ||
            (min_intercluster_distance_index == tree[cluster_step].left_child  && min_intercluster_distances[min_intercluster_distance_index].second == tree[cluster_step].right_child))
        {
          //find new min intercluster distance
          min_intercluster_distances[tree[cluster_step].left_child].first = std::numeric_limits<float>::max();

          for (Size j = 0; j < clusters[tree[cluster_step].left_child].size(); ++j)
          {
            Size k(0);
            for (; k < tree[cluster_step].left_child; ++k)
            {
              for (Size l = 0; l < clusters[k].size(); ++l)
              {
                if (original.getValue(clusters[tree[cluster_step].left_child][j], clusters[k][l]) < min_intercluster_distances[tree[cluster_step].left_child].first)
                {
                  min_intercluster_distances[tree[cluster_step].left_child].first = original.getValue(clusters[tree[cluster_step].left_child][j], clusters[k][l]);
                  min_intercluster_distances[tree[cluster_step].left_child].second = k;
                }
              }
            }
            ++k;
            for (; k < clusters.size(); ++k)
            {
              for (Size l = 0; l < clusters[k].size(); ++l)
              {
                if (original.getValue(clusters[tree[cluster_step].left_child][j], clusters[k][l]) < min_intercluster_distances[tree[cluster_step].left_child].first)
                {
                  min_intercluster_distances[tree[cluster_step].left_child].first = original.getValue(clusters[tree[cluster_step].left_child][j], clusters[k][l]);
                  min_intercluster_distances[tree[cluster_step].left_child].second = k;
                }
              }
            }
          }

          min_intercluster_distance_index = 0;
          for (Size i = min_intercluster_distance_index + 1; i < min_intercluster_distances.size(); ++i)
          {
            if (min_intercluster_distances[min_intercluster_distance_index].first == -1)
            {
              min_intercluster_distance_index = i;
            }
            else if (min_intercluster_distances[i].first != -1 && min_intercluster_distances[i].first < min_intercluster_distances[min_intercluster_distance_index].first)
            {
              min_intercluster_distance_index = i;
            }
          }

        }
        else
        {
          if (min_intercluster_distances[tree[cluster_step].right_child].first < min_intercluster_distances[tree[cluster_step].left_child].first)
          {
            min_intercluster_distances[tree[cluster_step].left_child].first = min_intercluster_distances[tree[cluster_step].right_child].first;
            min_intercluster_distances[tree[cluster_step].left_child].second = min_intercluster_distances[tree[cluster_step].right_child].second;
          }
        }

        for (Size k = 0; k < min_intercluster_distances.size(); ++k)
        {
          if (min_intercluster_distances[k].second == tree[cluster_step].right_child)
          {
            min_intercluster_distances[k].second = tree[cluster_step].left_child;
          }
        }
      }

      //shortcut for single linkage generated hirachy as merging criterion is min intercluster distance
      if (tree_from_singlelinkage)
      {
        float dunn_index(0);
        if (max_intracluster_distance > 0)
        {
          dunn_index = tree[cluster_step + 1].distance / max_intracluster_distance;
        }
        all_dunn_indices.push_back(dunn_index);
      }
      else
      {
        //find max dunn index and deduct the corresponding cluster step
        float dunn_index(0);
        if (max_intracluster_distance > 0)
        {
          dunn_index = min_intercluster_distances[min_intercluster_distance_index].first / max_intracluster_distance;
        }
        all_dunn_indices.push_back(dunn_index);
      }

      /* to manually retrace
          std::cout << min_intercluster_distance << std::endl;
          std::cout << clusters_with_min_intercluster_dist.first << " , " << clusters_with_min_intercluster_dist.second << std::endl;
          std::cout << max_intracluster_distance << std::endl;
      */

    }
    all_dunn_indices.push_back(0.0);     //last one is clearly 0
    return all_dunn_indices;
  }

  void ClusterAnalyzer::cut(const Size cluster_quantity, const std::vector<BinaryTreeNode> & tree, std::vector<std::vector<Size> > & clusters)
  {
    if (cluster_quantity == 0)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "You requested 0 clusters. Minimal partitioning contains one cluster, not zero.");
    }
    if (cluster_quantity > tree.size() + 1)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Cluster count exceeds node count. No partitioning exists.");
    }

    std::set<Size> leafs;
    for (Size i = 0; i < tree.size(); ++i)
    {
      leafs.insert(tree[i].left_child);
      leafs.insert(tree[i].right_child);
    }

    std::map<Size, std::vector<Size> > cluster_map;
    std::set<Size>::iterator it = leafs.begin();
    for (; it != leafs.end(); ++it)
    {
      cluster_map[*it] = std::vector<Size>(1, *it);
    }

    //redo clustering till step (original.dimensionsize()-cluster_quantity)
    for (Size cluster_step = 0; cluster_step < tree.size() + 1 - cluster_quantity; ++cluster_step)
    {
      if (tree[cluster_step].distance == -1)
      {
        break;
      }
      //pushback elements of right_child to left_child (and then erase second)
      cluster_map[tree[cluster_step].left_child].insert(cluster_map[tree[cluster_step].left_child].end(), cluster_map[tree[cluster_step].right_child].begin(), cluster_map[tree[cluster_step].right_child].end());

      // erase second one
      cluster_map[tree[cluster_step].right_child].clear();
    }

    // convert Map to Vector
    std::map<Size, std::vector<Size> >::iterator iter;
    for (iter = cluster_map.begin(); iter != cluster_map.end(); ++iter)
    {
      if (iter->second.empty())
        continue;
      std::vector<Size> actCluster = iter->second;
      clusters.push_back(actCluster);
    }
    //~ sorts by first element contained!!
    for (Size cluster_num = 0; cluster_num < clusters.size(); ++cluster_num)
    {
      std::sort(clusters[cluster_num].begin(), clusters[cluster_num].end());
    }
    std::sort(clusters.begin(), clusters.end());
    std::reverse(clusters.begin(), clusters.end());
    clusters.resize(cluster_quantity);
    std::sort(clusters.begin(), clusters.end());
  }

  void ClusterAnalyzer::cut(const Size cluster_quantity, const std::vector<BinaryTreeNode> & tree, std::vector<std::vector<BinaryTreeNode> > & subtrees)
  {
    if (cluster_quantity == 0)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "minimal partition contains one cluster, not zero");
    }
    if (cluster_quantity >= tree.size() + 1)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "maximal partition contains singleton clusters, further separation is not possible");
    }
    subtrees.clear();
    subtrees.resize(cluster_quantity);

    std::vector<std::vector<Size> > clusters;
    cut(cluster_quantity, tree, clusters);

    //~ unused nodes are discarded, (tree.begin()+tree.size()+1-cluster_quantity) is maximal tree.end() since cluster_quantity is always > 1! (tree.end()==tree.begin()+tree.size())
    std::list<BinaryTreeNode> tc(tree.begin(), (tree.begin() + (tree.size() + 1 - cluster_quantity)));
    for (Size cluster = 0; cluster < clusters.size(); ++cluster)
    {
      std::sort(clusters[cluster].begin(), clusters[cluster].end());
      std::list<BinaryTreeNode>::iterator it = tc.begin();
      while (it != tc.end())
      {
        std::vector<Size>::iterator left = std::find(clusters[cluster].begin(), clusters[cluster].end(), it->left_child);
        std::vector<Size>::iterator right = std::find(clusters[cluster].begin(), clusters[cluster].end(), it->right_child);
        if ((left != clusters[cluster].end() || right != clusters[cluster].end()))
        {
          subtrees[cluster].push_back(*it);
          it = tc.erase(it);
        }
        else
        {
          ++it;
        }
      }
    }
  }

  float ClusterAnalyzer::averagePopulationAberration(Size cluster_quantity, std::vector<BinaryTreeNode> & tree)
  {
    if (cluster_quantity == 0)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "minimal partition contains one cluster, not zero");
    }
    if (cluster_quantity >= tree.size() + 1)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "maximal partition contains singleton clusters, further separation is not possible");
    }
    std::vector<float> average_sizes;
    average_sizes.reserve(tree.size() + 1);
    std::vector<std::vector<Size> > clusters;
    clusters.reserve(tree.size() + 1);

    clusters.clear();
    clusters.reserve(tree.size() + 1);
    for (Size i = 0; i < tree.size() + 1; ++i)
    {
      clusters.emplace_back(1, i);
    }
    //redo clustering till step (original.dimensionsize()-cluster_quantity)
    for (Size cluster_step = 0; cluster_step < tree.size() + 1 - cluster_quantity; ++cluster_step)
    {
      //pushback elements of right_child to left_child (and then erase second)
      clusters[tree[cluster_step].left_child].insert(clusters[tree[cluster_step].left_child].end(), clusters[tree[cluster_step].right_child].begin(), clusters[tree[cluster_step].right_child].end());

      // clear second one
      clusters[tree[cluster_step].right_child].clear();
    }

    float average = (float)(tree.size() + 1) / (float)cluster_quantity;
    float aberration(0);
    float cluster_number(0);
    for (Size i = 0; i < clusters.size(); ++i)
    {
      if (!clusters[i].empty())
      {
        aberration += std::fabs((float)clusters[i].size() - average);
        ++cluster_number;
      }
    }
    aberration /= cluster_number;

    return aberration;
  }

  std::vector<float> ClusterAnalyzer::cohesion(const std::vector<std::vector<Size> > & clusters, const DistanceMatrix<float> & original)
  {
    if (clusters.empty() || clusters.size() > original.dimensionsize())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "invalid clustering");
    }

    float av_dist(0);     // average of all pairwise distances
    for (Size i = 0; i < original.dimensionsize(); ++i)
    {
      for (Size j = i + 1; j < original.dimensionsize(); ++j)
      {
        av_dist += original.getValue(i, j);
      }
    }
    av_dist /= (((float)original.dimensionsize() * (float)(original.dimensionsize() - 1.0)) / 2.0f);

    std::vector<float> cohesions;
    cohesions.reserve(clusters.size());
    for (Size i = 0; i < clusters.size(); ++i)
    {
      float av_c_dist(0);       // all pairwise distances in cluster i
      for (Size j = 0; j < clusters[i].size(); ++j)
      {
        for (Size k = 0; k < j; ++k)
        {
          av_c_dist += original.getValue(clusters[i][j], clusters[i][k]);
        }
      }

      av_c_dist /= (((float)clusters[i].size() * (float)(clusters[i].size() - 1.0)) / 2.0f); //now av. intra cluster distance
      if (clusters[i].size() == 1)
      {
        av_c_dist = av_dist;
      }
      //~ std::cout << " av clu i " << av_c_dist << std::endl;
      cohesions.push_back(av_c_dist);
    }
    return cohesions;
  }

  String ClusterAnalyzer::newickTree(const std::vector<BinaryTreeNode> & tree, const bool include_distance)
  {
    std::set<Size> leafs;
    for (Size i = 0; i < tree.size(); ++i)
    {
      leafs.insert(tree[i].left_child);
      leafs.insert(tree[i].right_child);
    }

    std::vector<String> clusters(*(leafs.rbegin()) + 1, "");
    for (std::set<Size>::iterator it = leafs.begin(); it != leafs.end(); ++it)
    {
      clusters[*it] = String(*it);
    }

    //redo clustering till step (original.dimensionsize()-1)
    for (Size cluster_step = 0; cluster_step < tree.size(); ++cluster_step)
    {
      //append string right_child to left_child
      clusters[tree[cluster_step].left_child].insert(0, "( ");
      if (include_distance)
      {
        clusters[tree[cluster_step].left_child] += ":";
        clusters[tree[cluster_step].left_child] += String(tree[cluster_step].distance);
      }
      clusters[tree[cluster_step].left_child] += " , ";
      clusters[tree[cluster_step].left_child] += clusters[tree[cluster_step].right_child];
      if (include_distance)
      {
        clusters[tree[cluster_step].left_child] += ":";
        clusters[tree[cluster_step].left_child] += String(tree[cluster_step].distance);
      }
      clusters[tree[cluster_step].left_child] += " )";

      clusters[tree[cluster_step].right_child] = String("");
    }

    Size first_filled(0);
    for (Size i = 0; i < clusters.size(); ++i)
    {
      if (!clusters[i].empty())
      {
        first_filled = i;
        break;
      }
    }
    for (Size i = first_filled + 1; i < clusters.size(); ++i)
    {
      if (!clusters[i].empty())
      {
        clusters[first_filled].insert(0, "( ");
        if (include_distance)
        {
          clusters[first_filled] += ":";
          clusters[first_filled] += String("1");
        }
        clusters[first_filled] += " , ";
        clusters[first_filled] += clusters[i];
        if (include_distance)
        {
          clusters[first_filled] += ":";
          clusters[first_filled] += String("1");
        }
        clusters[first_filled] += " )";
      }
    }
    return clusters[first_filled];
    //~example inspectable with: http://cgi-www.daimi.au.dk/cgi-chili/phyfi/go [BMC Bioinformatics 2006, 7:315]
  }

  bool compareBinaryTreeNode(const BinaryTreeNode & x, const BinaryTreeNode & y)
  {
    return x.distance < y.distance;
  }

}
