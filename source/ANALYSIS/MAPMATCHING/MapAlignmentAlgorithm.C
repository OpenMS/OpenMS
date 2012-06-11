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
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithm.h>

// Helper class to apply transformations to all sorts of KERNEL types
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>

// Derived classes are included here
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmSpectrumAlignment.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>

using namespace std;

namespace OpenMS
{
  //register products here
  void MapAlignmentAlgorithm::registerChildren()
  {
    Factory<MapAlignmentAlgorithm>::registerProduct(MapAlignmentAlgorithmIdentification::getProductName(),    &MapAlignmentAlgorithmIdentification::create);
    Factory<MapAlignmentAlgorithm>::registerProduct(MapAlignmentAlgorithmPoseClustering::getProductName(),    &MapAlignmentAlgorithmPoseClustering::create);
    Factory<MapAlignmentAlgorithm>::registerProduct(MapAlignmentAlgorithmSpectrumAlignment::getProductName(), &MapAlignmentAlgorithmSpectrumAlignment::create);
  }

  MapAlignmentAlgorithm::MapAlignmentAlgorithm() :
    DefaultParamHandler("MapAlignmentAlgorithm"),
    ProgressLogger()
  {
  }

  MapAlignmentAlgorithm::~MapAlignmentAlgorithm()
  {
  }

  void MapAlignmentAlgorithm::alignPeakMaps(vector<MSExperiment<> > &, vector<TransformationDescription> &)
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
  }

  void MapAlignmentAlgorithm::alignCompactFeatureMaps(vector<std::vector<Peak2D> > &, vector<TransformationDescription> &)
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
  }

  void MapAlignmentAlgorithm::alignFeatureMaps(vector<FeatureMap<> > &, vector<TransformationDescription> &)
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
  }

  void MapAlignmentAlgorithm::alignConsensusMaps(vector<ConsensusMap> & cms, vector<TransformationDescription> & tf)
  {
    LOG_WARN << "MapAlignmentAlgorithm::alignConsensusMaps() does not support ConsensusMaps directly. Converting to FeatureMaps." << endl;

    vector<FeatureMap<> > maps_f;
    for (Size i = 0; i < cms.size(); ++i)
    {
      FeatureMap<> fm;
      ConsensusMap::convert(cms[i], true, fm);
      maps_f.push_back(fm);
    }
    // call FeatureMap version of group()
    alignFeatureMaps(maps_f, tf);
    // apply transform
    MapAlignmentTransformer::transformConsensusMaps(cms, tf);
  }

  void MapAlignmentAlgorithm::alignPeptideIdentifications(vector<vector<PeptideIdentification> > &, vector<TransformationDescription> &)
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
  }

  void MapAlignmentAlgorithm::setReference(Size reference_index,
                                           const String & reference_file)
  {
    if (reference_index || !reference_file.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "This algorithm does not support a reference for the alignment.");
    }
  }

  void MapAlignmentAlgorithm::fitModel(const String & model_type, const Param & params, vector<TransformationDescription> & trafos)
  {
    for (vector<TransformationDescription>::iterator it = trafos.begin();
         it != trafos.end(); ++it)
    {
      it->fitModel(model_type, params);
    }
  }

}
