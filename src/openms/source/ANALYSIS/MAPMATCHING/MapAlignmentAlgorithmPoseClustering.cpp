// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>

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

  namespace
  {
    // most intense first; equal intensities are ordered by RT, then m/z, so that the selection and its order do not
    // depend on the standard library. A strict weak ordering only without NaN intensities (rejected by the callers).
    template <typename T>
    bool moreIntense(const T& left, const T& right)
    {
      if (left.getIntensity() != right.getIntensity()) return left.getIntensity() > right.getIntensity();
      if (left.getRT() != right.getRT()) return left.getRT() < right.getRT();
      return left.getMZ() < right.getMZ();
    }

    void throwIfNaN(float intensity, const char* what)
    {
      if (std::isnan(intensity))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      std::string(what) + " with a NaN intensity cannot be aligned.", "NaN");
      }
    }
  } // namespace

  void MapAlignmentAlgorithmPoseClustering::selectPeaks_(UInt64 map_index, const PeakMap& map, ConsensusMap& output, Size n)
  {
    output.clear(true);

    // n is capped by the number of MS1 peaks (not by map.getSize(), which also counts the peaks of other MS levels
    // and all chromatogram points)
    Size ms1_peaks = 0;
    for (const MSSpectrum& spectrum : map)
    {
      if (spectrum.getMSLevel() == 1)
      {
        ms1_peaks += spectrum.size();
      }
    }
    n = std::min(n, ms1_peaks);

    // keep the n most intense peaks in a heap whose front is the least intense of them, instead of copying all peaks
    const auto more_intense = [](const Peak2D& left, const Peak2D& right) { return moreIntense(left, right); };
    std::vector<Peak2D> selected;
    selected.reserve(n);
    for (const MSSpectrum& spectrum : map)
    {
      if (spectrum.getMSLevel() != 1)
      {
        continue;
      }
      for (const Peak1D& peak : spectrum)
      {
        throwIfNaN(peak.getIntensity(), "An MS1 peak");
        if (n == 0)
        {
          continue;
        }
        Peak2D candidate;
        candidate.setRT(spectrum.getRT());
        candidate.setMZ(peak.getMZ());
        candidate.setIntensity(peak.getIntensity());
        if (selected.size() < n)
        {
          selected.push_back(candidate);
          std::push_heap(selected.begin(), selected.end(), more_intense);
        }
        else if (more_intense(candidate, selected.front()))
        {
          std::pop_heap(selected.begin(), selected.end(), more_intense);
          selected.back() = candidate;
          std::push_heap(selected.begin(), selected.end(), more_intense);
        }
      }
    }
    std::sort_heap(selected.begin(), selected.end(), more_intense); // most intense first

    output.reserve(selected.size());
    for (Size element_index = 0; element_index < selected.size(); ++element_index)
    {
      output.push_back(ConsensusFeature(map_index, selected[element_index], element_index));
    }

    output.getColumnHeaders()[map_index].size = selected.size();
    output.updateRanges();
  }

  void MapAlignmentAlgorithmPoseClustering::selectFeatures_(UInt64 map_index, const FeatureMap& map, ConsensusMap& output, Size n)
  {
    for (const Feature& feature : map)
    {
      throwIfNaN(feature.getIntensity(), "A feature");
    }
    n = std::min(n, map.size());

    // indices of the n most intense features (ties by RT, m/z and index), in map order
    std::vector<Size> selected(map.size());
    std::iota(selected.begin(), selected.end(), 0);
    if (n < map.size())
    {
      std::nth_element(selected.begin(), selected.begin() + n, selected.end(), [&map](Size left, Size right)
      {
        if (moreIntense(map[left], map[right])) return true;
        if (moreIntense(map[right], map[left])) return false;
        return left < right;
      });
      selected.resize(n);
      std::sort(selected.begin(), selected.end());
    }

    // as MapConversion::convert(UInt64, const FeatureMap&, ConsensusMap&, Size)
    output.clear(true);
    output.reserve(n);
    output.setUniqueId(map.getUniqueId());
    for (Size index : selected)
    {
      output.push_back(ConsensusFeature(map_index, map[index]));
    }
    output.getColumnHeaders()[map_index].size = map.size();
    output.setProteinIdentifications(map.getProteinIdentifications());
    output.setUnassignedPeptideIdentifications(map.getUnassignedPeptideIdentifications());
    output.updateRanges();
  }

  void MapAlignmentAlgorithmPoseClustering::setReference(const FeatureMap& map)
  {
    selectFeatures_(0, map, reference_, max_num_peaks_considered_);
  }

  void MapAlignmentAlgorithmPoseClustering::setReference(const PeakMap& map)
  {
    selectPeaks_(0, map, reference_, max_num_peaks_considered_);
  }

  void MapAlignmentAlgorithmPoseClustering::align(const FeatureMap& map, TransformationDescription& trafo)
  {
    ConsensusMap map_scene;
    selectFeatures_(1, map, map_scene, max_num_peaks_considered_);
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
