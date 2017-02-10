// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Johannes Veit $
// $Authors: Johannes Veit $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmKD.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmKD.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

using namespace std;

namespace OpenMS
{

  FeatureGroupingAlgorithmKD::FeatureGroupingAlgorithmKD() :
    ProgressLogger()
  {
    setName("FeatureGroupingAlgorithmKD");
    defaults_.setValue("rt_tol", 60.0, "width of RT tolerance window (sec)");
    defaults_.setValue("mz_tol", 15.0, "m/z tolerance (in ppm or Da)");
    defaults_.setValue("mz_unit", "ppm", "unit of m/z tolerance");
    defaults_.setValidStrings("mz_unit", ListUtils::create<String>("ppm,Da"));
    defaults_.setValue("warp", "true", "whether or not to internally warp feature RTs using LOWESS transformation before linking (reported RTs in results will always be the original RTs)");
    defaults_.setValidStrings("warp", ListUtils::create<String>("true,false"));
    defaults_.setValue("min_rel_cc_size", 0.5, "only relevant during RT transformation: only connected components containing at least (warp_min_occur * number_of_input_maps) features are considered for computing the warping function");
    defaults_.setMinFloat("min_rel_cc_size", 0.0);
    defaults_.setMaxFloat("min_rel_cc_size", 1.0);
    defaults_.setValue("max_nr_conflicts", 0, "only relevant during RT transformation: allow up to this many conflicts (features from the same map) per connected component to be used for alignment (-1 means allow any number of conflicts)");
    defaults_.setMinInt("max_nr_conflicts", -1);
    defaults_.setValue("nr_partitions", 100, "number of partitions in m/z space");
    defaults_.setMinInt("nr_partitions", 1);
    defaultsToParam_();
    setLogType(CMD);
  }

  FeatureGroupingAlgorithmKD::~FeatureGroupingAlgorithmKD()
  {
  }

  template <typename MapType>
  void FeatureGroupingAlgorithmKD::group_(const vector<MapType>& input_maps,
                                          ConsensusMap& out)
  {
    rt_tol_secs_ = (double)(param_.getValue("rt_tol"));
    mz_tol_ = (double)(param_.getValue("mz_tol"));
    mz_ppm_ = param_.getValue("mz_unit").toString() == "ppm";

    // check that the number of maps is ok:
    if (input_maps.size() < 2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "At least two maps must be given!");
    }

    out.clear(false);

    vector<double> massrange;
    for (typename vector<MapType>::const_iterator map_it = input_maps.begin();
         map_it != input_maps.end(); ++map_it)
    {
      for (typename MapType::const_iterator feat_it = map_it->begin();
          feat_it != map_it->end(); feat_it++)
      {
        massrange.push_back(feat_it->getMZ());
      }
    }
    sort(massrange.begin(), massrange.end());

    // partition at boundaries -> this should be safe because there cannot be
    // any cluster reaching across boundaries

    int pts_per_partition = massrange.size() / (int)(param_.getValue("nr_partitions"));

    // compute partition boundaries
    vector<double> partition_boundaries;
    partition_boundaries.push_back(massrange.front());
    for (size_t j = 0; j < massrange.size()-1; j++)
    {
      // minimal differences between two m/z values
      double massrange_diff = mz_ppm_ ? mz_tol_ * 1e-6 * massrange[j+1] : mz_tol_;

      if (fabs(massrange[j] - massrange[j+1]) > massrange_diff)
      {
        if (j >= (partition_boundaries.size() ) * pts_per_partition  )
        {
          partition_boundaries.push_back((massrange[j] + massrange[j+1])/2.0);
        }
      }
    }
    // add last partition (a bit more since we use "smaller than" below)
    partition_boundaries.push_back(massrange.back() + 1.0);

    // ------------ compute RT transformation models ------------

    MapAlignmentAlgorithmKD aligner(input_maps.size(), param_);
    bool align = param_.getValue("warp").toString() == "true";
    if (align)
    {
      Size progress = 0;
      startProgress(0, partition_boundaries.size(), "computing RT transformations");
      for (size_t j = 0; j < partition_boundaries.size()-1; j++)
      {
        double partition_start = partition_boundaries[j];
        double partition_end = partition_boundaries[j+1];

        std::vector<MapType> tmp_input_maps(input_maps.size());
        for (size_t k = 0; k < input_maps.size(); k++)
        {
          // iterate over all features in the current input map and append
          // matching features (within the current partition) to the temporary
          // map
          for (size_t m = 0; m < input_maps[k].size(); m++)
          {
            if (input_maps[k][m].getMZ() >= partition_start &&
                input_maps[k][m].getMZ() < partition_end)
            {
              tmp_input_maps[k].push_back(input_maps[k][m]);
            }
          }
          tmp_input_maps[k].updateRanges();
        }

        // set up kd-tree
        KDTreeFeatureMaps kd_data(tmp_input_maps, param_);
        aligner.addRTFitData(kd_data);
        setProgress(progress++);
      }

      // fit LOWESS on RT fit data collected across all partitions
      try
      {
        aligner.fitLOWESS();
      }
      catch (Exception::BaseException& e)
      {
        LOG_ERROR << "Error: " << e.what() << endl;
        return;
      }

      endProgress();
    }

    // ------------ run alignment + feature linking on individual partitions ------------
    Size progress = 0;
    startProgress(0, partition_boundaries.size(), "linking features");
    for (size_t j = 0; j < partition_boundaries.size()-1; j++)
    {
      double partition_start = partition_boundaries[j];
      double partition_end = partition_boundaries[j+1];

      std::vector<MapType> tmp_input_maps(input_maps.size());
      for (size_t k = 0; k < input_maps.size(); k++)
      {
        // iterate over all features in the current input map and append
        // matching features (within the current partition) to the temporary
        // map
        for (size_t m = 0; m < input_maps[k].size(); m++)
        {
          if (input_maps[k][m].getMZ() >= partition_start &&
              input_maps[k][m].getMZ() < partition_end)
          {
            tmp_input_maps[k].push_back(input_maps[k][m]);
          }
        }
        tmp_input_maps[k].updateRanges();
      }

      // set up kd-tree
      KDTreeFeatureMaps kd_data(tmp_input_maps, param_);

      // alignment
      if (align)
      {
        aligner.transform(kd_data);
      }

      // link features
      runClustering_(kd_data, out);
      setProgress(progress++);
    }
    endProgress();

    // add protein IDs and unassigned peptide IDs to the result map here,
    // to keep the same order as the input maps (useful for output later):
    for (typename vector<MapType>::const_iterator map_it = input_maps.begin();
         map_it != input_maps.end(); ++map_it)
    {
      // add protein identifications to result map:
      out.getProteinIdentifications().insert(
        out.getProteinIdentifications().end(),
        map_it->getProteinIdentifications().begin(),
        map_it->getProteinIdentifications().end());

      // add unassigned peptide identifications to result map:
      out.getUnassignedPeptideIdentifications().insert(
        out.getUnassignedPeptideIdentifications().end(),
        map_it->getUnassignedPeptideIdentifications().begin(),
        map_it->getUnassignedPeptideIdentifications().end());
    }

    // canonical ordering for checking the results:
    startProgress(0, 3, String("sorting results"));
    out.sortByQuality();
    setProgress(1);
    out.sortByMaps();
    setProgress(2);
    out.sortBySize();
    endProgress();
    return;
  }

  void FeatureGroupingAlgorithmKD::group(const std::vector<FeatureMap>& maps,
                                         ConsensusMap& out)
  {
    group_(maps, out);
  }

  void FeatureGroupingAlgorithmKD::group(const std::vector<ConsensusMap>& maps,
                                         ConsensusMap& out)
  {
    group_(maps, out);
  }

  void FeatureGroupingAlgorithmKD::runClustering_(const KDTreeFeatureMaps& kd_data, ConsensusMap& out)
  {
    Size n = kd_data.size();

    // pass 1: initialize best potential clusters for all possible cluster centers
    set<Size> update_these;
    for (Size i = 0; i < kd_data.size(); ++i)
    {
      update_these.insert(i);
    }
    set<ClusterProxyKD> potential_clusters;
    vector<ClusterProxyKD> cluster_for_idx(n);
    vector<Int> assigned(n, false);
    updateClusterProxies_(potential_clusters, cluster_for_idx, update_these, assigned, kd_data);

    // pass 2: construct consensus features until all points assigned.
    while (!potential_clusters.empty())
    {
      // get index of current best cluster center (as defined by ClusterProxyKD::operator<)
      Size i = potential_clusters.begin()->getCenterIndex();

      // compile the actual list of sub feature indices for cluster with center i
      vector<Size> cf_indices;
      computeBestClusterForCenter_(i, cf_indices, assigned, kd_data);

      // add consensus feature
      addConsensusFeature_(cf_indices, kd_data, out);

      // mark selected sub features assigned and delete them from potential_clusters
      for (vector<Size>::const_iterator f_it = cf_indices.begin(); f_it != cf_indices.end(); ++f_it)
      {
        assigned[*f_it] = true;
        potential_clusters.erase(cluster_for_idx[*f_it]);
      }

      // compile set of all points whose neighborhoods will need updating
      update_these = set<Size>();
      for (vector<Size>::const_iterator f_it = cf_indices.begin(); f_it != cf_indices.end(); ++f_it)
      {
        vector<Size> f_neighbors;
        kd_data.getNeighborhood(*f_it, f_neighbors, true);
        for (vector<Size>::const_iterator it = f_neighbors.begin(); it != f_neighbors.end(); ++it)
        {
          if (!assigned[*it])
          {
            update_these.insert(*it);
          }
        }
      }

      // now that the points are marked assigned, update the neighborhoods of their neighbors
      updateClusterProxies_(potential_clusters, cluster_for_idx, update_these, assigned, kd_data);
    }
  }



  void FeatureGroupingAlgorithmKD::updateClusterProxies_(set<ClusterProxyKD>& potential_clusters,
                                                         vector<ClusterProxyKD>& cluster_for_idx,
                                                         const set<Size>& update_these,
                                                         const vector<Int>& assigned,
                                                         const KDTreeFeatureMaps& kd_data)
  {
    for (set<Size>::const_iterator it = update_these.begin(); it != update_these.end(); ++it)
    {
      Size i = *it;
      const ClusterProxyKD& old_proxy = cluster_for_idx[i];
      vector<Size> unused;
      ClusterProxyKD new_proxy = computeBestClusterForCenter_(i, unused, assigned, kd_data);

      // only need to update if size and/or average distance have changed
      if (new_proxy != old_proxy)
      {
        potential_clusters.erase(old_proxy);
        cluster_for_idx[i] = new_proxy;
        potential_clusters.insert(new_proxy);
      }
    }
  }

  ClusterProxyKD FeatureGroupingAlgorithmKD::computeBestClusterForCenter_(Size i, vector<Size>& cf_indices, const vector<Int>& assigned, const KDTreeFeatureMaps& kd_data) const
  {
    // for scaling distances relative to tolerance windows below
    pair<double, double> dummy_mz_win = Math::getTolWindow(1000, mz_tol_, mz_ppm_);
    double max_rt_dist = rt_tol_secs_;
    double max_mz_dist = (dummy_mz_win.second - dummy_mz_win.first) / 2;

    // compute i's neighborhood, together with a look-up table
    // map index -> corresponding points
    map<Size, vector<Size> > points_for_map_index;
    vector<Size> neighbors;
    kd_data.getNeighborhood(i, neighbors, true);
    Int charge_i = kd_data.charge(i);
    for (vector<Size>::const_iterator it = neighbors.begin(); it != neighbors.end(); ++it)
    {
      if (!assigned[*it] && kd_data.charge(*it) == charge_i)
      {
        points_for_map_index[kd_data.mapIndex(*it)].push_back(*it);
      }
    }
    // center i is always part of CF, no other points from i's map can be contained
    points_for_map_index[kd_data.mapIndex(i)] = vector<Size>(1, i);

    // compile list of sub feature indices and corresponding average distance
    double avg_distance = 0.0;
    for (map<Size, vector<Size> >::const_iterator it = points_for_map_index.begin();
         it != points_for_map_index.end();
         ++it)
    {
      const vector<Size>& candidates = it->second;

      // choose a point j with minimal distance to center i
      double min_dist = numeric_limits<double>::max();
      Size best_index = numeric_limits<double>::max();
      for (vector<Size>::const_iterator c_it = candidates.begin(); c_it != candidates.end(); ++c_it)
      {
        double dist = distance_(kd_data.mz(*c_it) / max_mz_dist,
                                kd_data.rt(*c_it) / max_rt_dist,
                                kd_data.mz(i) / max_mz_dist,
                                kd_data.rt(i) / max_rt_dist);
        if (dist < min_dist)
        {
          min_dist = dist;
          best_index = *c_it;
        }
      }
      cf_indices.push_back(best_index);
      avg_distance += min_dist;
    }
    avg_distance /= cf_indices.size();

    return ClusterProxyKD(cf_indices.size(), avg_distance, i);
  }

  void FeatureGroupingAlgorithmKD::addConsensusFeature_(const vector<Size>& indices, const KDTreeFeatureMaps& kd_data, ConsensusMap& out) const
  {
    ConsensusFeature cf;
    float avg_quality = 0;
    for (vector<Size>::const_iterator it = indices.begin(); it != indices.end(); ++it)
    {
      Size i = *it;
      cf.insert(kd_data.mapIndex(i), *(kd_data.feature(i)));
      avg_quality += kd_data.feature(i)->getQuality();
    }
    avg_quality /= indices.size();
    cf.setQuality(avg_quality);
    cf.computeConsensus();
    out.push_back(cf);
  }

  double FeatureGroupingAlgorithmKD::distance_(double mz_1, double rt_1, double mz_2, double rt_2) const
  {
    return sqrt(pow((mz_1 - mz_2), 2) + pow((rt_1 - rt_2), 2));
  }

} // namespace OpenMS
