// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <OpenMS/FEATUREFINDER/MultiplexDeltaMasses.h>
#include <OpenMS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteredMSExperiment.h>
#include <OpenMS/PROCESSING/MISC/SplinePackage.h>
#include <OpenMS/ML/CLUSTERING/GridBasedCluster.h>

#include <vector>
#include <fstream>
#include <map>

namespace OpenMS
{

class OPENMS_DLLAPI FeatureFinderMultiplexAlgorithm :
  public DefaultParamHandler, public ProgressLogger
{
public:
  /// default constructor
  FeatureFinderMultiplexAlgorithm();

  /// main method for feature detection
  void run(MSExperiment& exp, bool progress);

  /// get methods
  FeatureMap& getFeatureMap();
  ConsensusMap& getConsensusMap();
  MSExperiment& getBlacklist();

protected:

  // experimental data
  MSExperiment exp_profile_;
  MSExperiment exp_centroid_;

  bool centroided_;

  ProgressLogger prog_log_;

  bool progress_;

  unsigned charge_min_;
  unsigned charge_max_;

  unsigned isotopes_per_peptide_min_;
  unsigned isotopes_per_peptide_max_;


  // mass shift names and their values
  std::map<String, double> label_mass_shift_;

  // final results, maps of detected features
  FeatureMap feature_map_;
  ConsensusMap consensus_map_;

  // blacklist
  MSExperiment exp_blacklist_;

  /**
   * @brief generate list of m/z shifts
   *
   * @param charge_min    minimum charge
   * @param charge_max    maximum charge
   * @param peaks_per_peptide_max    maximum number of isotopes in peptide
   * @param mass_pattern_list    mass shifts due to labelling
   *
   * @return list of m/z shifts
   */
  std::vector<MultiplexIsotopicPeakPattern> generatePeakPatterns_(int charge_min, int charge_max, int peaks_per_peptide_max, const std::vector<MultiplexDeltaMasses>& mass_pattern_list);

  /**
   * @brief determine ratios through linear regression and correct peptide intensities
   *
   * In most labelled mass spectrometry experiments, the fold change i.e. ratio and not the individual peptide intensities
   * are of primary interest. For that reason, we determine the ratios from interpolated chromatogram data points directly,
   * and then correct the current ones.
   *
   */
  void correctPeptideIntensities_(const MultiplexIsotopicPeakPattern& pattern, std::map<size_t, SplinePackage>& spline_chromatograms, const std::vector<double>& rt_peptide, std::vector<double>& intensity_peptide) const;

  /**
   * @brief calculate peptide intensities
   *
   * @param pattern
   * @param satellites
   *
   * @return vector with intensities for each of the peptides
   */
  std::vector<double> determinePeptideIntensitiesCentroided_(const MultiplexIsotopicPeakPattern& pattern, const std::multimap<size_t, MultiplexSatelliteCentroided >& satellites);

  /**
   * @brief calculate peptide intensities
   *
   * @param pattern
   * @param satellites
   *
   * @return vector with intensities for each of the peptides
   */
  std::vector<double> determinePeptideIntensitiesProfile_(const MultiplexIsotopicPeakPattern& pattern, const std::multimap<size_t, MultiplexSatelliteProfile >& satellites);

  /**
   * @brief generates consensus and feature maps containing all peptide multiplets
   *
   * @param patterns    patterns of isotopic peaks we have been searching for
   * @param filter_results    filter results for each of the patterns
   * @param cluster_results    clusters of filter results
   */
  void generateMapsCentroided_(const std::vector<MultiplexIsotopicPeakPattern>& patterns, const std::vector<MultiplexFilteredMSExperiment>& filter_results, std::vector<std::map<int, GridBasedCluster> >& cluster_results);

  /**
   * @brief generates consensus and feature maps containing all peptide multiplets
   *
   * @param patterns    patterns of isotopic peaks we have been searching for
   * @param filter_results    filter results for each of the patterns
   * @param cluster_results    clusters of filter results
   */
  void generateMapsProfile_(const std::vector<MultiplexIsotopicPeakPattern>& patterns, const std::vector<MultiplexFilteredMSExperiment>& filter_results, const std::vector<std::map<int, GridBasedCluster> >& cluster_results);

};

}
