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
// $Maintainer: Clemens Groepl $
// $Authors: Eva Lange, Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMPOSECLUSTERING_H
#define OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMPOSECLUSTERING_H

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithm.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/StablePairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringAffineSuperimposer.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

namespace OpenMS
{
  /**
    @brief A map alignment algorithm based on pose clustering.

    Pose clustering analyzes pair distances to find the most probable transformation of retention times.
		The algorithm choses the x most intensity peaks/features per map.
    This is modeled via the parameter 'max_num_peaks_considered', which in turn influences the runtime and stability of the results.
		Bigger values prolong computation, smaller values might lead to no or unstable trafos. Set to -1 to use all features (might take very
		long for large maps).
	
    For further details see:
    @n Eva Lange et.al
    @n A Geometric Approach for the Alignment of Liquid Chromatography-Mass Spectrometry Data
    @n ISMB/ECCB 2007

    @htmlinclude OpenMS_MapAlignmentAlgorithmPoseClustering.parameters

    @ingroup MapAlignment

  */
  class OPENMS_DLLAPI MapAlignmentAlgorithmPoseClustering :
    public MapAlignmentAlgorithm
  {
public:
    /// Default constructor
    MapAlignmentAlgorithmPoseClustering();

    /// Destructor
    virtual ~MapAlignmentAlgorithmPoseClustering();

    void align(const FeatureMap<>& map, TransformationDescription& trafo);
    void align(const MSExperiment<>& map, TransformationDescription& trafo);
    void align(const ConsensusMap& map, TransformationDescription& trafo);

    template <typename MapType> void setReference( const MapType& map )
    {
      MapType map2 = map; // todo: avoid copy (MSExperiment version of convert() demands non-const version)
      ConsensusMap::convert(0, map2, reference_, max_num_peaks_considered_);
    }

    /// Creates a new instance of this class (for Factory)
    static MapAlignmentAlgorithm * create()
    {
      return new MapAlignmentAlgorithmPoseClustering();
    }

    /// Returns the product name (for the Factory)
    static String getProductName()
    {
      return "pose_clustering";
    }

protected:

    virtual void updateMembers_();

    PoseClusteringAffineSuperimposer superimposer_;

    StablePairFinder pairfinder_;

    ConsensusMap reference_;

    Int max_num_peaks_considered_;

private:

    /// Copy constructor intentionally not implemented -> private
    MapAlignmentAlgorithmPoseClustering(const MapAlignmentAlgorithmPoseClustering &);
    ///Assignment operator intentionally not implemented -> private
    MapAlignmentAlgorithmPoseClustering & operator=(const MapAlignmentAlgorithmPoseClustering &);
  };
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMPOSECLUSTERING_H
