// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

#include <algorithm>

using namespace std;

namespace OpenMS
{

  MapAlignmentAlgorithmPoseClustering::MapAlignmentAlgorithmPoseClustering() :
    DefaultParamHandler("MapAlignmentAlgorithmPoseClustering"), 
    ProgressLogger(), max_num_peaks_considered_(0)
  {
    defaults_.insert("superimposer:", PoseClusteringAffineSuperimposer().getParameters());
    defaults_.insert("pairfinder:", StablePairFinder().getParameters());
    defaults_.setValue("max_num_peaks_considered", 1000, "The maximal number of peaks/features to be considered per map. To use all, set to '-1'.");
    defaults_.setMinInt("max_num_peaks_considered", -1);

    defaultsToParam_();
  }

  void MapAlignmentAlgorithmPoseClustering::updateMembers_()
  {
    superimposer_.setParameters(param_.copy("superimposer:", true));
    superimposer_.setLogType(getLogType());

    pairfinder_.setParameters(param_.copy("pairfinder:", true));
    pairfinder_.setLogType(getLogType());

    max_num_peaks_considered_ = param_.getValue("max_num_peaks_considered");
  }

  MapAlignmentAlgorithmPoseClustering::~MapAlignmentAlgorithmPoseClustering() = default;

  void MapAlignmentAlgorithmPoseClustering::selectPeaks_(UInt64 map_index, const PeakMap& map, ConsensusMap& output, Size n)
  {
    output.clear(true);

    std::vector<Peak2D> peaks;
    peaks.reserve(map.getSize()); // an upper bound only, see below

    // TODO Avoid tripling the memory consumption by this call
    map.get2DData(peaks);

    // Clamp n only now, against the number of peaks actually collected.
    // map.getSize() is the wrong bound: it counts the peaks of every
    // spectrum at every MS level plus all chromatogram points, whereas
    // get2DData() collects MS1 peaks only. With n > peaks.size() the middle
    // iterator of the partial_sort and the copy loop below would run past
    // the end of peaks (out-of-bounds reads and writes, consensus features
    // built from garbage).
    if (n > peaks.size())
    {
      n = peaks.size();
    }
    output.reserve(n);

    // most intense first; equal intensities are ordered by RT, then m/z, so that the selected peaks and their
    // order do not depend on the standard library's partial_sort implementation
    std::partial_sort(peaks.begin(),
                      peaks.begin() + n,
                      peaks.end(),
                      [](const Peak2D& left, const Peak2D& right)
                      {
                        if (left.getIntensity() != right.getIntensity()) return left.getIntensity() > right.getIntensity();
                        if (left.getRT() != right.getRT()) return left.getRT() < right.getRT();
                        return left.getMZ() < right.getMZ();
                      });

    for (Size element_index = 0; element_index < n; ++element_index)
    {
      output.push_back(ConsensusFeature(map_index, peaks[element_index], element_index));
    }

    output.getColumnHeaders()[map_index].size = n;
    output.updateRanges();
  }

  void MapAlignmentAlgorithmPoseClustering::setReference(const FeatureMap& map)
  {
    MapConversion::convert(0, map, reference_, max_num_peaks_considered_);
  }

  void MapAlignmentAlgorithmPoseClustering::setReference(const PeakMap& map)
  {
    selectPeaks_(0, map, reference_, max_num_peaks_considered_);
  }

  void MapAlignmentAlgorithmPoseClustering::align(const FeatureMap& map, TransformationDescription& trafo)
  {
    ConsensusMap map_scene;
    MapConversion::convert(1, map, map_scene, max_num_peaks_considered_);
    align(map_scene, trafo);
  }

  void MapAlignmentAlgorithmPoseClustering::align(const PeakMap& map, TransformationDescription& trafo)
  {
    ConsensusMap map_scene;
    selectPeaks_(1, map, map_scene, max_num_peaks_considered_);
    align(map_scene, trafo);
  }

  void MapAlignmentAlgorithmPoseClustering::align(const ConsensusMap& map, TransformationDescription& trafo)
  {
    // TODO: move this to updateMembers_? (if ConsensusMap prevails)
    // TODO: why does superimposer work on consensus map???
    const ConsensusMap & map_model = reference_;
    ConsensusMap map_scene = map;

    // run superimposer to find the global transformation
    TransformationDescription si_trafo;
    superimposer_.run(map_model, map_scene, si_trafo);

    // apply transformation to consensus features and contained feature
    // handles
    for (Size j = 0; j < map_scene.size(); ++j)
    {
      // Calculate new RT
      double rt = map_scene[j].getRT();
      rt = si_trafo.apply(rt);
      // Set RT of consensus feature centroid
      map_scene[j].setRT(rt);
      // Set RT of consensus feature handles
      map_scene[j].begin()->asMutable().setRT(rt);
    }

    // run pairfinder to find pairs
    ConsensusMap result;
    // TODO: add another 2-map interface to pairfinder?
    std::vector<ConsensusMap> input(2);
    input[0] = map_model;
    input[1] = map_scene;
    pairfinder_.run(input, result);

    // calculate the local transformation
    si_trafo.invert(); // to undo the transformation applied above
    TransformationDescription::DataPoints data;
    for (ConsensusFeature& cfeature : result)
    {
      if (cfeature.size() == 2) // two matching features
      {
        ConsensusFeature::iterator feat_it = cfeature.begin();
        double y = feat_it->getRT();
        double x = si_trafo.apply((++feat_it)->getRT());
        // one feature should be from the reference map:
        if (feat_it->getMapIndex() != 0)
        {
          data.push_back(make_pair(x, y));
        }
        else
        {
          data.push_back(make_pair(y, x));
        }
      }
    }
    trafo = TransformationDescription(data);
    trafo.fitModel("linear");
  }

} // namespace
