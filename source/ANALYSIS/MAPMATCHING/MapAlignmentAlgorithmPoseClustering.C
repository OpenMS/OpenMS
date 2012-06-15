// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
  {}


  void MapAlignmentAlgorithmPoseClustering::align( const FeatureMap<>& map, TransformationDescription& trafo )
  {
    ConsensusMap map_scene;
    ConsensusMap::convert(1, map, map_scene, max_num_peaks_considered_);
    align(map_scene, trafo);
  }

  void MapAlignmentAlgorithmPoseClustering::align( const MSExperiment<>& map, TransformationDescription& trafo )
  {
    ConsensusMap map_scene;
	MSExperiment<> map2(map);
    ConsensusMap::convert(1, map2, map_scene, max_num_peaks_considered_); // copy MSExperiment here, since it is sorted internally by intensity
    align(map_scene, trafo);
  }

  void MapAlignmentAlgorithmPoseClustering::align( const ConsensusMap& map, TransformationDescription& trafo )
  {
    // TODO: move this to updateMembers_? (if consensusMap prevails)
    // TODO: why does superimposer work on consensus map???
    const ConsensusMap& map_model = reference_;
    ConsensusMap map_scene = map;

    // run superimposer to find the global transformation
    TransformationDescription si_trafo;
    superimposer_.run(map_model, map_scene, si_trafo);

    // apply transformation to consensus features and contained feature
    // handles
    for (Size j = 0; j < map_scene.size(); ++j)
    {
      //Calculate new RT
      DoubleReal rt = map_scene[j].getRT();
      rt = si_trafo.apply(rt);
      //Set RT of consensus feature centroid
      map_scene[j].setRT(rt);
      //Set RT of consensus feature handles
      map_scene[j].begin()->asMutable().setRT(rt);
    }

    //run pairfinder to find pairs
    ConsensusMap result;
    //TODO: add another 2map interface to pairfinder?
    std::vector< ConsensusMap > input(2);
    input[0] = map_model;
    input[1] = map_scene;
    pairfinder_.run(input, result);

    // calculate the local transformation
    si_trafo.invert();         // to undo the transformation applied above
    TransformationDescription::DataPoints data;
    for (ConsensusMap::Iterator it = result.begin(); it != result.end();
      ++it)
    {
      if (it->size() == 2)           // two matching features
      {
        ConsensusFeature::iterator feat_it = it->begin();
        DoubleReal y = feat_it->getRT();
        DoubleReal x = si_trafo.apply((++feat_it)->getRT());
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

} //namespace
