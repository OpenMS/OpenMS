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

  /**
  FeatureFinderMultiplexAlgorithm is a tool for the fully automated analysis of quantitative proteomics data. It detects pairs of isotopic envelopes with fixed m/z separation.
It requires no prior sequence identification of the peptides and works on both profile or centroided spectra. In what follows we outline the algorithm.

<b>Algorithm</b>
The algorithm is divided into three parts: filtering, clustering and linear fitting, see Fig. (d), (e) and (f).
In the following discussion let us consider a particular mass spectrum at retention time 1350 s, see Fig. (a).
It contains a peptide of mass 1492 Da and its 6 Da heavier labelled counterpart. Both are doubly charged in this instance.
Their isotopic envelopes therefore appear at 746 and 749 in the spectrum. The isotopic peaks within each envelope are separated by 0.5.
The spectrum was recorded at finite intervals. In order to read accurate intensities at arbitrary m/z we spline-fit over the data, see Fig. (b).
We would like to search for such peptide pairs in our LC-MS data set. As a warm-up let us consider a standard intensity cut-off filter, see Fig. (c).

Scanning through the entire m/z range (red dot) only data points with intensities above a certain threshold pass the filter. 
Unlike such a local filter, the filter used in our algorithm takes intensities at a range of m/z positions into account, see Fig. (d). A data point (red dot) passes if
- all six intensities at m/z, m/z+0.5, m/z+1, m/z+3, m/z+3.5 and m/z+4 lie above a certain threshold,
- the intensity profiles in neighbourhoods around all six m/z positions show a good correlation and
- the relative intensity ratios within a peptide agree up to a factor with the ratios of a theoretic averagine model.

Let us now filter not only a single spectrum but all spectra in our data set. Data points that pass the filter form clusters in the t-m/z plane, see Fig. (e).
Each cluster corresponds to the mono-isotopic mass trace of the lightest peptide of a SILAC pattern. We now use hierarchical clustering methods to assign each data point to a specific cluster.
The optimum number of clusters is determined by maximizing the silhouette width of the partitioning.
Each data point in a cluster corresponds to three pairs of intensities (at [m/z, m/z+3], [m/z+0.5, m/z+3.5] and [m/z+1, m/z+4]).
A plot of all intensity pairs in a cluster shows a clear linear correlation, see Fig. (f).
Using linear regression we can determine the relative amounts of labelled and unlabelled peptides in the sample.

@image html SILACAnalyzer_algorithm.png

  */

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
