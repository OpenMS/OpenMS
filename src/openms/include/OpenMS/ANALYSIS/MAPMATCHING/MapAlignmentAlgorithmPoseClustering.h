// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange, Clemens Groepl $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/StablePairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringAffineSuperimposer.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

#include <OpenMS/KERNEL/ConversionHelper.h>

namespace OpenMS
{
  /**
    @brief A map alignment algorithm based on pose clustering.

    Pose clustering analyzes pair distances to find the most probable
    transformation of retention times.

    The algorithm chooses the x most intensity peaks/features per map.  This is
    modeled via the parameter 'max_num_peaks_considered', which in turn
    influences the runtime and stability of the results.  Bigger values prolong
    computation, smaller values might lead to no or unstable trafos. Set to -1
    to use all features (might take very long for large maps).

    For further details see:
    @n Eva Lange et al.
    @n A Geometric Approach for the Alignment of Liquid Chromatography-Mass Spectrometry Data
    @n ISMB/ECCB 2007

    @htmlinclude OpenMS_MapAlignmentAlgorithmPoseClustering.parameters

    @ingroup MapAlignment

  */
  class OPENMS_DLLAPI MapAlignmentAlgorithmPoseClustering :
    public DefaultParamHandler,
    public ProgressLogger
  {
public:
    /// Default constructor
    MapAlignmentAlgorithmPoseClustering();

    /// Destructor
    ~MapAlignmentAlgorithmPoseClustering() override;

    void align(const FeatureMap& map, TransformationDescription& trafo);
    void align(const PeakMap& map, TransformationDescription& trafo);
    void align(const ConsensusMap& map, TransformationDescription& trafo);

    /// Sets the reference for the alignment
    template <typename MapType>
    void setReference(const MapType& map)
    {
      MapType map2 = map; // todo: avoid copy (MSExperiment version of convert() demands non-const version)
      MapConversion::convert(0, map2, reference_, max_num_peaks_considered_);
    }

protected:

    void updateMembers_() override;

    PoseClusteringAffineSuperimposer superimposer_;

    StablePairFinder pairfinder_;

    ConsensusMap reference_;

    Int max_num_peaks_considered_;

private:

    /// Copy constructor intentionally not implemented -> private
    MapAlignmentAlgorithmPoseClustering(const MapAlignmentAlgorithmPoseClustering&);
    /// Assignment operator intentionally not implemented -> private
    MapAlignmentAlgorithmPoseClustering& operator=(const MapAlignmentAlgorithmPoseClustering&);
  };
} // namespace OpenMS

