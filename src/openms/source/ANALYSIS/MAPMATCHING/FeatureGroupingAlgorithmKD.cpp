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
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>

using namespace std;

namespace OpenMS
{

  FeatureGroupingAlgorithmKD::FeatureGroupingAlgorithmKD() :
    ProgressLogger(),
    kd_data_()
  {
    setName("FeatureGroupingAlgorithmKD");
    defaults_.setValue("rt_tol", 60.0, "width of RT tolerance window (sec)");
    defaults_.setValue("mz_tol", 15.0, "m/z tolerance (in ppm or Da)");
    defaults_.setValue("mz_unit", "ppm", "unit of m/z tolerance");
    defaults_.setValidStrings("mz_unit", ListUtils::create<String>("ppm,Da"));
    defaultsToParam_();
    setLogType(CMD);
  }

  FeatureGroupingAlgorithmKD::~FeatureGroupingAlgorithmKD()
  {
  }

  template <typename MapType>
  void FeatureGroupingAlgorithmKD::group_(const vector<MapType>& maps,
                                          ConsensusMap& out)
  {
    rt_tol_secs_ = (double)(param_.getValue("rt_tol"));
    mz_tol_ = (double)(param_.getValue("mz_tol"));
    mz_ppm_ = param_.getValue("mz_unit").toString() == "ppm";

    // check that the number of maps is ok:
    if (maps.size() < 2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                       "At least two maps must be given!");
    }

    // set up kd-tree
    setUpTree_(maps);

    // link features
    runClustering_(out);

    // add protein IDs and unassigned peptide IDs to the result map here,
    // to keep the same order as the input maps (useful for output later):
    for (typename vector<MapType>::const_iterator map_it = maps.begin();
         map_it != maps.end(); ++map_it)
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
    startProgress(0, 1, String("sorting results"));
    out.sortByQuality();
    out.sortByMaps();
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

  template <typename MapType>
  void FeatureGroupingAlgorithmKD::setUpTree_(const std::vector<MapType>& maps)
  {
    kd_data_.clear();
    kd_data_.setParameters(param_);

    progress_ = 0;
    startProgress(0, maps.size(), String("building kd-tree"));
    for (Size i = 0; i < maps.size(); ++i)
    {
      const MapType& m = maps[i];
      for (typename MapType::const_iterator it = m.begin(); it != m.end(); ++it)
      {
        kd_data_.addFeature(i, &(*it));
      }
      setProgress(++progress_);
    }
    endProgress();

    LOG_INFO << "Added " << kd_data_.size() << " features from " << maps.size() << " maps." << endl;

    startProgress(0, 1, String("optimizing k-d tree"));
    kd_data_.optimizeTree();
    endProgress();
  }

  void FeatureGroupingAlgorithmKD::runClustering_(ConsensusMap& out)
  {
    Size n = kd_data_.size();

    // pass 1: initialize best potential clusters for all possible cluster centers
    progress_ = 0;
    startProgress(0, n, String("linking features (pass 1/2)"));
    set<Size> update_these;
    for (Size i = 0; i < kd_data_.size(); ++i)
    {
      update_these.insert(i);
    }
    set<ClusterProxyKD> potential_clusters;
    vector<ClusterProxyKD> cluster_for_idx(n);
    vector<Int> assigned(n, false);
    updateClusterProxies_(potential_clusters, cluster_for_idx, update_these, assigned, true);
    endProgress();

    // pass 2: construct consensus features until all points assigned.
    progress_ = 0;
    startProgress(0, n, String("linking features (pass 2/2)"));
    while(!potential_clusters.empty())
    {
      // get index of current best cluster center (as defined by ClusterProxyKD::operator<)
      Size i = potential_clusters.begin()->getCenterIndex();

      // compile the actual list of sub feature indices for cluster with center i
      vector<Size> cf_indices;
      computeBestClusterForCenter(i, cf_indices, assigned);

      // add consensus feature
      addConsensusFeature_(cf_indices, out);

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
        kd_data_.getNeighborhood(*f_it, f_neighbors, true);
        for (vector<Size>::const_iterator it = f_neighbors.begin(); it != f_neighbors.end(); ++it)
        {
          if (!assigned[*it])
          {
            update_these.insert(*it);
          }
        }
      }

      // now that the points are marked assigned, update the neighborhoods of their neighbors
      updateClusterProxies_(potential_clusters, cluster_for_idx, update_these, assigned);
      progress_ += cf_indices.size();
      setProgress(progress_);
    }

    endProgress();
  }



  void FeatureGroupingAlgorithmKD::updateClusterProxies_(set<ClusterProxyKD>& potential_clusters,
                                                         vector<ClusterProxyKD>& cluster_for_idx,
                                                         const set<Size>& update_these,
                                                         const vector<Int>& assigned,
                                                         bool update_progress)
  {
    for (set<Size>::const_iterator it = update_these.begin(); it != update_these.end(); ++it)
    {
      Size i = *it;
      const ClusterProxyKD& old_proxy = cluster_for_idx[i];
      vector<Size> unused;
      ClusterProxyKD new_proxy = computeBestClusterForCenter(i, unused, assigned);

      // only need to update if size and/or average distance have changed
      if (new_proxy != old_proxy)
      {
        potential_clusters.erase(old_proxy);
        cluster_for_idx[i] = new_proxy;
        potential_clusters.insert(new_proxy);
      }
      if (update_progress)
      {
        setProgress(++progress_);
      }
    }
  }

  ClusterProxyKD FeatureGroupingAlgorithmKD::computeBestClusterForCenter(Size i, vector<Size>& cf_indices, const vector<Int>& assigned) const
  {
    // for scaling distances relative to tolerance windows below
    pair<double, double> dummy_mz_win = kd_data_.getTolWindow(1000, mz_tol_, mz_ppm_);
    double max_rt_dist = rt_tol_secs_;
    double max_mz_dist = (dummy_mz_win.second - dummy_mz_win.first) / 2;

    // compute i's neighborhood, together with a look-up table
    // map index -> corresponding points
    map<Size, vector<Size> > points_for_map_index;
    vector<Size> neighbors;
    kd_data_.getNeighborhood(i, neighbors, true);
    Int charge_i = kd_data_.charge(i);
    for (vector<Size>::const_iterator it = neighbors.begin(); it != neighbors.end(); ++it)
    {
      if (!assigned[*it] && kd_data_.charge(*it) == charge_i)
      {
        points_for_map_index[kd_data_.mapIndex(*it)].push_back(*it);
      }
    }
    // center i is always part of CF, no other points from i's map can be contained
    points_for_map_index[kd_data_.mapIndex(i)] = vector<Size>(1, i);

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
        double dist = distance_(kd_data_.mz(*c_it) / max_mz_dist,
                                kd_data_.rt(*c_it) / max_rt_dist,
                                kd_data_.mz(i) / max_mz_dist,
                                kd_data_.rt(i) / max_rt_dist);
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

  void FeatureGroupingAlgorithmKD::addConsensusFeature_(const vector<Size>& indices, ConsensusMap& out) const
  {
    ConsensusFeature cf;
    float avg_quality = 0;
    for (vector<Size>::const_iterator it = indices.begin(); it != indices.end(); ++it)
    {
      Size i = *it;
      cf.insert(kd_data_.mapIndex(i), *(kd_data_.feature(i)));
      avg_quality += kd_data_.feature(i)->getQuality();
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
