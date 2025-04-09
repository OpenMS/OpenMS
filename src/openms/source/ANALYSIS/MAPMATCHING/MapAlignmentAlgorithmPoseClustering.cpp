// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>

#include <OpenMS/FORMAT/FileHandler.h>

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

  void MapAlignmentAlgorithmPoseClustering::align(const FeatureMap& map, TransformationDescription& trafo)
  {
    ConsensusMap map_scene;
    MapConversion::convert(1, map, map_scene, max_num_peaks_considered_);
    align(map_scene, trafo);
  }

  void MapAlignmentAlgorithmPoseClustering::align(const PeakMap& map, TransformationDescription& trafo)
  {
    ConsensusMap map_scene;
    PeakMap map2(map);
    MapConversion::convert(1, map2, map_scene, max_num_peaks_considered_); // copy MSExperiment here, since it is sorted internally by intensity
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
