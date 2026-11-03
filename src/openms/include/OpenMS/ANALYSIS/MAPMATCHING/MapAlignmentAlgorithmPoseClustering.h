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
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{
  /**
    @brief A map alignment algorithm based on pose clustering.

    Pose clustering analyzes pair distances to find the most probable
    transformation of retention times.

    The algorithm uses up to 'max_num_peaks_considered' peaks/features per
    map: the most intense features of a feature map, or MS1 peaks of a peak
    map selected as set by 'peak_selection:method' (see below). This parameter
    influences the runtime and stability of the results.  Bigger values prolong
    computation, smaller values might lead to no or unstable trafos. Set to -1
    to use all features (might take very long for large maps).

    For peak maps, 'peak_selection:method' 'spread' (the default) selects at
    most one peak per analyte (instead of many peaks of the elution profile and
    isotope pattern of a few intense analytes), spread over the RT and m/z
    range; see selectSpreadPeaks_(). 'intensity' selects the most intense
    peaks; see selectPeaks_().

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

    /// Sets the reference for the alignment: the 'max_num_peaks_considered' most intense features of @p map
    void setReference(const FeatureMap& map);

    /// Sets the reference for the alignment: 'max_num_peaks_considered' MS1 peaks of @p map, selected as set by 'peak_selection:method'
    void setReference(const PeakMap& map);

protected:

    void updateMembers_() override;

    /**
      @brief Copies the @p n most intense MS1 peaks of @p map into @p output, as consensus features of map @p map_index

      @p output is cleared first. Only the peaks of MS level 1 spectra are used; spectra of other MS levels and
      chromatograms do not contribute. The peaks are ordered by descending intensity; equal intensities are ordered
      by RT, then m/z. @p n is capped by the number of MS1 peaks, so @c Size(-1) keeps all of them. Only the
      selected peaks are held in memory besides @p map. The column header @c size for @p map_index is the number of
      peaks copied.

      The consensus features only carry the positions and intensities the alignment needs: they have no unique ids,
      and the single handle of each stores the index of the peak in @p output.

      @exception Exception::InvalidValue if an MS1 peak has a NaN intensity or a non-finite m/z, or a non-empty MS1
                 spectrum has a non-finite RT
    */
    static void selectPeaks_(UInt64 map_index, const PeakMap& map, ConsensusMap& output, Size n);

    /**
      @brief Copies up to @p n MS1 peaks of @p map into @p output that represent different analytes, spread over the
             RT and m/z range of @p map

      Many peaks belong to the same analyte: the peaks of its elution profile (the same m/z in consecutive spectra)
      and of its isotope pattern. Taking the most intense peaks (selectPeaks_()) therefore picks a few intense
      analytes many times, which the superimposer and the pair finder cannot use. This selection keeps one peak per
      analyte and spreads the peaks over the map:

      1. The MS1 peaks are binned into a grid of @p rt_bins x @p mz_bins cells over the RT and m/z range of the
         MS1 peaks. Each cell keeps only its most intense candidates, at first a multiple of its share of @p n, so
         that memory stays proportional to @p n. Steps 1 to 3 are repeated with more candidates for each cell whose
         candidates were cut off where that was too few: step 3 would take further peaks from the cell (e.g. because
         other cells are empty, or its candidates belong to only a few analytes), or a cut-off peak of the cell may
         belong to the analyte of a peak selected from another cell (e.g. an elution profile that crosses the cell
         border, with the apex cut off). The selection is thus practically the one with all MS1 peaks as candidates,
         at a fraction of the memory.
      2. Going from the most intense candidate down, a candidate is accepted unless an accepted peak lies within
         @p rt_window (in seconds) of it at the same m/z, or at an isotope position m/z +/- k * 1.0033548378 / z
         (k = 1..3, z = 1..@p max_charge), within @p mz_tolerance_ppm. The accepted peaks are the apices of
         different analytes.
      3. The accepted peaks are taken round-robin over the grid cells (the most intense of each cell first, then the
         second most intense of each cell, ...), until @p n peaks are selected.

      Fewer than @p n peaks are selected if fewer analytes are found. The selected peaks are ordered by descending
      intensity (equal intensities by RT, then m/z), and @p output is otherwise set up as by selectPeaks_().

      @exception Exception::InvalidValue if an MS1 peak has a NaN intensity or a non-finite m/z, or a non-empty MS1
                 spectrum has a non-finite RT
    */
    static void selectSpreadPeaks_(UInt64 map_index, const PeakMap& map, ConsensusMap& output, Size n,
                                   double rt_window, double mz_tolerance_ppm, Size max_charge, Size rt_bins, Size mz_bins);

    /// Selects the peaks of @p map to align, by selectPeaks_() or selectSpreadPeaks_() as set by 'peak_selection:method'
    void selectPeakLandmarks_(UInt64 map_index, const PeakMap& map, ConsensusMap& output) const;

    /**
      @brief Copies the @p n most intense features of @p map into @p output, as consensus features of map @p map_index

      Like MapConversion::convert(UInt64, const FeatureMap&, ConsensusMap&, Size), but selects the most intense
      features instead of the first ones (ties are broken by RT, then m/z, then position in @p map). The selected
      features keep their order in @p map. With @p n at or above the number of features (e.g. @c Size(-1)), the
      result equals that conversion.

      @exception Exception::InvalidValue if a feature has a NaN intensity or a non-finite RT or m/z
    */
    static void selectFeatures_(UInt64 map_index, const FeatureMap& map, ConsensusMap& output, Size n);

    PoseClusteringAffineSuperimposer superimposer_;

    StablePairFinder pairfinder_;

    ConsensusMap reference_;

    Int max_num_peaks_considered_;

    /// 'peak_selection:' parameters
    bool spread_peaks_ = true;
    double spread_rt_window_ = 20.0;
    double spread_mz_tolerance_ppm_ = 10.0;
    Size spread_max_charge_ = 4;
    Size spread_rt_bins_ = 10;
    Size spread_mz_bins_ = 10;

private:

    /// Copy constructor intentionally not implemented -> private
    MapAlignmentAlgorithmPoseClustering(const MapAlignmentAlgorithmPoseClustering&);
    /// Assignment operator intentionally not implemented -> private
    MapAlignmentAlgorithmPoseClustering& operator=(const MapAlignmentAlgorithmPoseClustering&);
  };
} // namespace OpenMS

