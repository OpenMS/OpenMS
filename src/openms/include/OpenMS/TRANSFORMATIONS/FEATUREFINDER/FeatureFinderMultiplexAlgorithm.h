// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMasses.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteredMSExperiment.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SplinePackage.h>
#include <OpenMS/COMPARISON/CLUSTERING/GridBasedCluster.h>

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
   * @param intensity_peptide    peptide intensities to be corrected
   */
  void correctPeptideIntensities_(const MultiplexIsotopicPeakPattern& pattern, std::map<size_t, SplinePackage>& spline_chromatograms, const std::vector<double>& rt_peptide, std::vector<double>& intensity_peptide);
  
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
