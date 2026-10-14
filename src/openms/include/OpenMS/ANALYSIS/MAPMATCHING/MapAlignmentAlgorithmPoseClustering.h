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

    /// Sets the reference for the alignment: the first 'max_num_peaks_considered' features of @p map
    void setReference(const FeatureMap& map);

    /// Sets the reference for the alignment: the 'max_num_peaks_considered' most intense MS1 peaks of @p map
    void setReference(const PeakMap& map);

protected:

    void updateMembers_() override;

    /**
      @brief Copies the @p n most intense MS1 peaks of @p map into @p output, as consensus features of map @p map_index

      @p output is cleared first. Only the peaks of MS level 1 spectra are used; spectra of other MS levels and
      chromatograms do not contribute. The peaks are ordered by descending intensity; equal intensities are ordered
      by RT, then m/z. @p n is capped by the number of MS1 peaks, so @c Size(-1) keeps all of them. The column header
      @c size for @p map_index is the number of peaks copied.

      The consensus features only carry the positions and intensities the alignment needs: they have no unique ids,
      and the single handle of each stores the index of the peak in @p output.
    */
    static void selectPeaks_(UInt64 map_index, const PeakMap& map, ConsensusMap& output, Size n);

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

