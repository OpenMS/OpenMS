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

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmKD.h>
#include <queue>

using namespace std;

namespace OpenMS
{

MapAlignmentAlgorithmKD::MapAlignmentAlgorithmKD(Size num_maps) :
  fit_data_(num_maps),
  transformations_(num_maps)
{
  setLogType(CMD);
}

MapAlignmentAlgorithmKD::~MapAlignmentAlgorithmKD()
{
}

void MapAlignmentAlgorithmKD::addRTFitData(const KDTreeData& kd_data)
{
  Size num_maps = kd_data.numMaps();
  Size min_cc_size = num_maps / 2; //TODO

  // compute connected components
  map<Size, vector<Size> > ccs;
  getCCs_(kd_data, ccs);
  // keep only conflict-free CCs of sufficient size
  map<Size, vector<Size> > filtered_ccs;
  filterCCs_(kd_data, filtered_ccs, ccs, min_cc_size);
  // save some memory
  ccs.clear();

  // compute average RTs for all CCs
  map<Size, double> avg_rts;
  for (map<Size, vector<Size> >::const_iterator it = filtered_ccs.begin(); it != filtered_ccs.end(); ++it)
  {
    double avg_rt = 0;
    Size cc_index = it->first;
    const vector<Size>& cc = it->second;
    for (vector<Size>::const_iterator cc_it = cc.begin(); cc_it != cc.end(); ++cc_it)
    {
      Size i = *cc_it;
      avg_rt += kd_data.rt(i);
    }
    avg_rt /= cc.size();
    avg_rts[cc_index] = avg_rt;
  }

  // generate fit data for each map, add to fit_data_
  for (map<Size, vector<Size> >::const_iterator it = filtered_ccs.begin(); it != filtered_ccs.end(); ++it)
  {
    Size cc_index = it->first;
    const vector<Size>& cc = it->second;
    for (vector<Size>::const_iterator cc_it = cc.begin(); cc_it != cc.end(); ++cc_it)
    {
      Size i = *cc_it;
      double rt = kd_data.rt(i);
      double avg_rt = avg_rts[cc_index];
      fit_data_[kd_data.mapIndex(i)].push_back(make_pair(rt, avg_rt));
    }
  }
}

void MapAlignmentAlgorithmKD::fitLOWESS()
{
  //SignedSize model_progress = 0;
  Size num_maps = fit_data_.size();
  //startProgress(0, num_maps, String("fitting LOWESS transformation models"));
  for (Size i = 0; i < num_maps; ++i)
  {
    transformations_[i] = new TransformationModelLowess(fit_data_[i], Param());
    //setProgress(++model_progress);
  }
  //endProgress();
}

void MapAlignmentAlgorithmKD::transform(KDTreeData& kd_data) const
{
  // apply transformations to kd_data
  //startProgress(0, 1, String("applying LOWESS transformations"));
  kd_data.applyTransformations(transformations_);
  //endProgress();

  // re-optimize kd-tree
  //startProgress(0, 1, String("re-optimizing kd-tree"));
  kd_data.optimizeTree();
  //endProgress();
}

Size MapAlignmentAlgorithmKD::computeCCs_(const KDTreeData& kd_data, vector<Size>& result) const
{
  //compute CCs by means of repeated BFS (without actually storing the graph (edges) in memory)

  Size num_nodes = kd_data.size();
  //startProgress(0, num_nodes, String("computing connected components"));

  //clear CC indices
  result.clear();
  result.resize(num_nodes, numeric_limits<Size>::max());

  //set up data structures
  queue<Size> bfs_queue;
  vector<Int> bfs_visited(num_nodes, false);
  Size search_pos = 0;
  Size cc_index = 0;
  //SignedSize cc_progress = 0;

  //BFS until every node has been visited
  while (true)
  {
    bool finished = true;
    for (Size i = search_pos; i < num_nodes; ++i)
    {
      if (!bfs_visited[i])
      {
        bfs_queue.push(i);
        bfs_visited[i] = true;
        finished = false;
        search_pos = i + 1;
        break;
      }
    }
    if (finished) break;

    while (!bfs_queue.empty())
    {
      Size i = bfs_queue.front();
      bfs_queue.pop();
      result[i] = cc_index;
      //setProgress(++cc_progress);

      vector<Size> compatible_features;
      kd_data.getNeighborhood(i, compatible_features);
      for (vector<Size>::const_iterator it = compatible_features.begin();
           it != compatible_features.end();
           ++it)
      {
        Size j = *it;
        if (!bfs_visited[j])
        {
          bfs_queue.push(j);
          bfs_visited[j] = true;
        }
      }
    }
    ++cc_index;
  }
  //endProgress();

  return cc_index;
}

void MapAlignmentAlgorithmKD::getCCs_(const KDTreeData& kd_data, map<Size, vector<Size> >& result) const
{
  vector<Size> cc_index;
  computeCCs_(kd_data, cc_index);

  result.clear();
  for (Size i = 0; i < kd_data.size(); ++i)
  {
    result[cc_index[i]].push_back(i);
  }
}

void MapAlignmentAlgorithmKD::filterCCs_(const KDTreeData& kd_data, map<Size, vector<Size> >& filtered_ccs, const map<Size, vector<Size> >& ccs, Size min_size) const
{
  //SignedSize filter_progress = 0;
  //startProgress(0, ccs.size(), String("filtering connected components"));

  Size max_size = fit_data_.size();
  filtered_ccs.clear();

  for (map<Size, vector<Size> >::const_iterator it = ccs.begin(); it != ccs.end(); ++it)
  {
    //setProgress(filter_progress++);
    const vector<Size>& cc = it->second;
    // size OK?
    if (cc.size() < min_size || cc.size() > max_size)
    {
      continue;
    }
    // map indices unique and diameter within user-specified tolerances?
    set<Size> map_indices;
    double rt_low = numeric_limits<double>::max();
    double rt_high = 0;
    double mz_low = numeric_limits<double>::max();
    double mz_high = 0;
    double rt_tol = kd_data.rtTolerance();
    double mz_tol;
    if (!kd_data.mzPPM())
    {
      mz_tol = kd_data.mzTolerance();
    }
    else
    {
      pair<double, double> win = kd_data.getTolWindow(kd_data.mz(cc[0]),
                                                      kd_data.mzTolerance(),
                                                      true);
      mz_tol = (win.second - win.first) / 2.0;
    }
    bool passes = true;
    for (vector<Size>::const_iterator idx_it = cc.begin(); idx_it != cc.end(); ++idx_it)
    {
      // filter out if reoccuring index found
      Size i = *idx_it;
      if (map_indices.find(i) != map_indices.end())
      {
        passes = false;
        break;
      }
      map_indices.insert(kd_data.mapIndex(i));
//      // filter out if diameter too large
//      double rt = kd_data.rt(i);
//      double mz = kd_data.mz(i);
//      rt_low = min(rt_low, rt);
//      rt_high = max(rt_high, rt);
//      mz_low = min(mz_low, mz);
//      mz_high = max(mz_high, mz);
//      if (rt_high - rt_low > rt_tol || mz_high - mz_low > mz_tol)
//      {
//        passes = false;
//        break;
//      }
    }
    if (passes)
    {
      filtered_ccs[it->first] = cc;
    }
  }
  //endProgress();
}

}
