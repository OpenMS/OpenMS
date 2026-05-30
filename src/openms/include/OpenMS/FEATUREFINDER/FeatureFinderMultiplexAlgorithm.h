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
    @brief Identification-free quantitative feature finder for SILAC / Dimethyl / ICPL-style label-pair experiments.

    Detects pairs (or larger multiplets) of isotopic envelopes separated by fixed label-induced
    m/z shifts. Works on both profile and centroided MS1 spectra; no prior sequence
    identification is needed. Configurable via @ref DefaultParamHandler (parameter sections
    @c "algorithm:*" and @c "labels:*"; defaults are set by the constructor).

    <b>Algorithm</b>
    Three stages — filtering, clustering, and linear fitting:
      - Filtering retains data points where all isotope-pattern peaks (e.g. at m/z, m/z+0.5,
        m/z+1, m/z+3, m/z+3.5, m/z+4 for a doubly charged pair) clear an intensity threshold,
        show good profile correlation across the m/z range, and whose within-peptide
        intensity ratios match an averagine model within a per-call tolerance.
      - Clustering uses hierarchical methods on the surviving points in the t-m/z plane.
        Each cluster corresponds to the monoisotopic mass trace of the lightest peptide of
        a label pattern. The number of clusters is chosen by maximising the silhouette
        width.
      - Linear fitting on the in-cluster intensity pairs (e.g. [m/z, m/z+3], [m/z+0.5,
        m/z+3.5], [m/z+1, m/z+4] for a doubly charged pair) yields the relative amounts of
        labelled and unlabelled peptide.

    @image html SILACAnalyzer_algorithm.png

    @ingroup Quantitation
  */
class OPENMS_DLLAPI FeatureFinderMultiplexAlgorithm :
  public DefaultParamHandler, public ProgressLogger
{
public:
  /// Construct with built-in defaults; parameter sections @c "algorithm" and @c "labels" are registered. Call @c setParameters to override before @ref run.
  FeatureFinderMultiplexAlgorithm();

  /**
    @brief Run the full detection pipeline on @p exp and populate the internal feature / consensus / blacklist maps.

    The pipeline:
      -# Parses the configured charge and isotopes-per-peptide ranges from the
         @c "algorithm:charge" and @c "algorithm:isotopes_per_peptide" Param values
         (formatted as @c "min:max"; swapped silently if the input is reversed).
      -# Loads the label mass shifts from the @c "labels:*" Param section into the
         in-class @c label_mass_shift_ map.
      -# Clears @p exp's chromatograms, calls @c updateRanges and @c sortSpectra on it,
         and then **moves** @p exp into either the internal centroided or profile
         working buffer via @c MSExperiment::swap — after this call returns, the
         caller's @p exp is left empty.
      -# Determines the spectrum type from @c "algorithm:spectrum_type":
         @c "automatic" reads @c exp[0].getType(true) and treats UNKNOWN as profile;
         @c "centroid" / @c "profile" force the choice.
      -# When the input is in profile mode, runs @ref PeakPickerHiRes::pickExperiment
         on @c "MS1" only with @c signal_to_noise @c = @c 0 (S/N estimation disabled).
      -# Builds the label delta-mass list with @ref MultiplexDeltaMassesGenerator,
         optionally extending it with knock-out delta masses when
         @c "algorithm:knock_out" is @c "true".
      -# Generates isotopic peak patterns, runs the filter / cluster / linear-fit
         stages described in the class brief, and populates @c feature_map_,
         @c consensus_map_, and @c exp_blacklist_ accessible through the @c get*
         methods.

    @param[in,out] exp      MS1 experiment to analyse. Chromatograms are cleared, ranges are recomputed, spectra sorted, and the container is moved into the internal working buffer — the caller's @p exp is empty on return.
    @param[in]     progress When @c true, progress is reported through the inherited @ref ProgressLogger.

    @throws OpenMS::Exception::FileEmpty when @p exp has no spectra.
  */
  void run(MSExperiment& exp, bool progress);

  /// Return the FeatureMap populated by the most recent @ref run call (empty before @ref run).
  FeatureMap& getFeatureMap();
  /// Return the ConsensusMap of detected multiplets populated by the most recent @ref run call (empty before @ref run).
  ConsensusMap& getConsensusMap();
  /// Return the per-pattern blacklist MSExperiment produced during @ref run (peak regions consumed by detected multiplets).
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
   * @param[in] charge_min    minimum charge
   * @param[in] charge_max    maximum charge
   * @param[in] peaks_per_peptide_max    maximum number of isotopes in peptide
   * @param[in] mass_pattern_list    mass shifts due to labelling
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
   * @param[in] pattern Isotopic peak pattern
   * @param[in,out] spline_chromatograms Spline chromatograms to be used/modified
   * @param[in] rt_peptide Retention times of peptides
   * @param[out] intensity_peptide Corrected peptide intensities
   */
  void correctPeptideIntensities_(const MultiplexIsotopicPeakPattern& pattern, std::map<size_t, SplinePackage>& spline_chromatograms, const std::vector<double>& rt_peptide, std::vector<double>& intensity_peptide) const;

  /**
   * @brief calculate peptide intensities
   *
   * @param[in] pattern Isotopic peak pattern
   * @param[in] satellites Satellite peaks
   *
   * @return vector with intensities for each of the peptides
   */
  std::vector<double> determinePeptideIntensitiesCentroided_(const MultiplexIsotopicPeakPattern& pattern, const std::multimap<size_t, MultiplexSatelliteCentroided >& satellites);

  /**
   * @brief calculate peptide intensities
   *
   * @param[in] pattern Isotopic peak pattern
   * @param[in] satellites Satellite peaks
   *
   * @return vector with intensities for each of the peptides
   */
  std::vector<double> determinePeptideIntensitiesProfile_(const MultiplexIsotopicPeakPattern& pattern, const std::multimap<size_t, MultiplexSatelliteProfile >& satellites);

  /**
   * @brief generates consensus and feature maps containing all peptide multiplets
   *
   * @param[in] patterns    patterns of isotopic peaks we have been searching for
   * @param[in] filter_results    filter results for each of the patterns
   * @param[in,out] cluster_results    clusters of filter results
   */
  void generateMapsCentroided_(const std::vector<MultiplexIsotopicPeakPattern>& patterns, const std::vector<MultiplexFilteredMSExperiment>& filter_results, std::vector<std::map<int, GridBasedCluster> >& cluster_results);

  /**
   * @brief generates consensus and feature maps containing all peptide multiplets
   *
   * @param[in] patterns    patterns of isotopic peaks we have been searching for
   * @param[in] filter_results    filter results for each of the patterns
   * @param[in] cluster_results    clusters of filter results
   */
  void generateMapsProfile_(const std::vector<MultiplexIsotopicPeakPattern>& patterns, const std::vector<MultiplexFilteredMSExperiment>& filter_results, const std::vector<std::map<int, GridBasedCluster> >& cluster_results);

};

}
