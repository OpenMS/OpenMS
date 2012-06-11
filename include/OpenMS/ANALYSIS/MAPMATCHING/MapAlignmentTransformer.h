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
//  version 2.1 of the License, or (at your option) any later version
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Stephan Aiche $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTTRANSFORMER_H
#define OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTTRANSFORMER_H

#include <vector>
#include <OpenMS/config.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenMS
{

  class TransformationDescription;
  class ConsensusMap;
  class PeptideIdentification;
  class ConsensusFeature;

  /**
   * @brief The MapAlignmentTransformer class
   */
  class OPENMS_DLLAPI MapAlignmentTransformer
  {

public:
    /// Applies the <i>given</i> transformations to peak maps
    static void transformPeakMaps(std::vector<MSExperiment<> > & maps, const std::vector<TransformationDescription> & given_trafos);

    /// Applies the <i>given</i> transformations to feature maps
    static void transformFeatureMaps(std::vector<FeatureMap<> > & maps, const std::vector<TransformationDescription> & given_trafos);

    /// Applies the <i>given</i> transformations to consensus maps
    static void transformConsensusMaps(std::vector<ConsensusMap> & maps, const std::vector<TransformationDescription> & given_trafos);

    /// Applies the <i>given</i> transformations to peptide identifications
    static void transformPeptideIdentifications(std::vector<std::vector<PeptideIdentification> > & maps, const std::vector<TransformationDescription> & given_trafos);


    /// Applies the <i>given</i> transformations to a single peak map
    static void transformSinglePeakMap(MSExperiment<> & msexp, const TransformationDescription & trafo);

    /// Applies the <i>given</i> transformations to a single feature map
    static void transformSingleFeatureMap(FeatureMap<> & fmap, const TransformationDescription & trafo);

    /// Applies the <i>given</i> transformations to a single consensus map
    static void transformSingleConsensusMap(ConsensusMap & cmap, const TransformationDescription & trafo);

    /// Applies the <i>given</i> transformations to a single peptide identification
    static void transformSinglePeptideIdentification(std::vector<PeptideIdentification> & pepids, const TransformationDescription & trafo);

private:

    /// apply a transformation to a feature
    static void applyToFeature_(Feature & feature, const TransformationDescription & trafo);

    /// apply a transformation to a basic feature
    static void applyToBaseFeature_(BaseFeature & feature, const TransformationDescription & trafo);

    /// apply a transformation to a consensus feature
    static void applyToConsensusFeature_(ConsensusFeature & feature, const TransformationDescription & trafo);

  };
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTTRANSFORMER_H
