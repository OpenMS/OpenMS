// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Clemens Groepl $
// $Authors: Eva Lange, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#include <iostream>

using namespace std;

namespace OpenMS
{

  MapAlignmentAlgorithmPoseClustering::MapAlignmentAlgorithmPoseClustering() :
    MapAlignmentAlgorithm(),
    max_num_peaks_considered_(0)
  {
    setName("MapAlignmentAlgorithmPoseClustering");

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

  MapAlignmentAlgorithmPoseClustering::~MapAlignmentAlgorithmPoseClustering()
  {
  }

  void MapAlignmentAlgorithmPoseClustering::align(const FeatureMap& map, TransformationDescription& trafo)
  {
    ConsensusMap map_scene;
    MapConversion::convert(1, map, map_scene, max_num_peaks_considered_);
    align(map_scene, trafo);
  }

  void MapAlignmentAlgorithmPoseClustering::align(const MSExperiment<>& map, TransformationDescription& trafo)
  {
    ConsensusMap map_scene;
    MSExperiment<> map2(map);
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
    for (ConsensusMap::Iterator it = result.begin(); it != result.end(); ++it)
    {
      if (it->size() == 2) // two matching features
      {
        ConsensusFeature::iterator feat_it = it->begin();
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
