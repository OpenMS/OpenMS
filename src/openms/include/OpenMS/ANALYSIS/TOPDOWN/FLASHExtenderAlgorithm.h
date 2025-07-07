// Copyright (c) 2002-2024, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <boost/dynamic_bitset.hpp>
#include <iomanip>
#include <iostream>

namespace OpenMS
{
/**
@brief Extend between tags found by FLASHTaggerAlgorithm. In practice, the proteoform characterization is
 performed here. Only blind modification is implemented, and variable modification search is under develolpment.
@ingroup Topdown
*/

class OPENMS_DLLAPI FLASHExtenderAlgorithm : public DefaultParamHandler, public ProgressLogger
{
public:
  /// constructor
  FLASHExtenderAlgorithm();

  /// destructor
  ~FLASHExtenderAlgorithm() override = default;

  /// copy constructor
  FLASHExtenderAlgorithm(const FLASHExtenderAlgorithm&) = default;

  /// move constructor
  FLASHExtenderAlgorithm(FLASHExtenderAlgorithm&& other) = default;

  /// assignment operator
  FLASHExtenderAlgorithm& operator=(const FLASHExtenderAlgorithm& other);

  /**
   * The main run function to perform extension algorithm. Take the candidate protein hits and perform extension for each protein with the tags.
   * @param hits the candidate protein hits
   * @param tags the sequence tags from FLASHTaggerAlgorithm
   * @param dspec the deconvolved spectrum
   * @param ppm mass ppm tolerance
   * @param multiple_hits_per_spec should multiple proteins be considered per spectrum or only the best protein should be considered?
   */
  void run(std::vector<ProteinHit>& hits, const std::vector<FLASHHelperClasses::Tag>& tags,
           const DeconvolvedSpectrum& dspec, double ppm, bool multiple_hits_per_spec);

  /// fill the characterized proteoforms in @p hits
  void fillProteoforms(std::vector<ProteinHit>& hits) const
  {
    for (const auto& hit : proteoform_hits_)
    {
      hits.push_back(hit);
    }
  }

  /// set modification map to be considered. Note that they are not for variable modification search. The mass shifts found by blind modifications matching to the masses of these modifications are annotated.
  void setModificationMap(const std::map<double, std::vector<ResidueModification>>& mod_map)
  {
    mod_map_ = mod_map;
  }

protected:
  void updateMembers_() override;
  /// implemented for DefaultParamHandler
  void setDefaultParams_();

private:
  /// Protein hit information that should be retained throughout the run
  struct OPENMS_DLLAPI HitInformation
  {
  public:
    /// if a vertex is already considered (visited) when building a DAG.
    boost::dynamic_bitset<> visited_;
    /// The DAG consisting of paths representing proteoforms
    FLASHHelperClasses::DAG dag_;
    /// We perform three times of extensions - for suffix, for prefix, and for both termini
    int mode_;
    /// For each mode, information on each mass should be kept. These retain the information of each mass
    std::map<int, MSSpectrum> node_spec_map_, tol_spec_map_;
    /// For each mode, protein suffix or prefix masses are stored here.
    std::map<int, std::vector<double>> pro_mass_map_;
    /// protein start and end position for protein truncation, calculated after extension (after proteoform characterization).
    int protein_start_position_ = -1, protein_end_position_ = -1;
    /// calculated precursor mass from fragment mass pairing (representing complementary ion pairs)
    double calculated_precursor_mass_ = -1;
  };

  /// modification mass to modification index. To use find nearest function
  std::map<double, std::vector<ResidueModification>> mod_map_;

  /// get protein prefix or suffix masses (depending on the mode)
  static void getProMasses_(const ProteinHit& hit, std::vector<double>& pro_masses, int mode);

  /// calculate precursor mass with fragment mass pairs
  void calculatePrecursorMass_(const ProteinHit& hit,
                               const std::map<int, std::vector<Size>>& best_path_map,
                               HitInformation& hi);

  /// define nodes for DAG
  void defineNodes_(const DeconvolvedSpectrum& dspec, HitInformation& hi,
                   double max_mass);

  /// main run function per protein hit. The results are stored in @p all_paths_per_mode.
  void run_(const ProteinHit& hit, HitInformation& hi,
            const std::vector<FLASHHelperClasses::Tag>& matched_tags,
            std::map<int, std::vector<Size>>& all_paths_per_mode, int max_mod_cntr_for_last_mode); // per hit

  Size getVertex_(int node_index, int pro_index, int score, int num_mod, Size pro_mass_size) const;
  int getNodeIndex_(Size vertex, Size pro_mass_size) const;
  int getProIndex_(Size vertex, Size pro_mass_size) const;
  int getModNumber_(Size vertex) const;
  int getScore_(Size vertex) const;
  void constructDAG_(std::set<Size>& sinks,
                     HitInformation& hi,
                     const std::vector<std::vector<int>>& tag_edges,
                     int max_mod_cntr_for_last_mode,
                     bool use_tags);

  void connectBetweenTags_(std::set<Size>& visited_tag_edges,
                           HitInformation& hi,
                           std::map<Size, std::set<std::pair<double, double>>>& sinks,
                           Size vertex,
                           double truncation_mass,
                           double cumulative_shift,
                           std::map<Size, std::map<Size, int>>& node_max_score_map,
                           const std::vector<std::vector<int>>& tag_edges,
                           int max_mod_cntr_for_last_mode,
                           bool use_tags);

  void extendBetweenTags_(std::map<Size, std::set<std::pair<double, double>>>& sinks,
                          HitInformation& hi,
                          Size start_vertex,
                          int end_node_index,
                          int end_pro_index,
                          int diagonal_counter,
                          double truncation_mass,
                          double cumulative_mod_mass,
                          std::map<Size, std::map<Size, int>>& node_max_score_map,
                          int max_mod_cntr_for_last_mode);

  int getProteinLength_(const std::vector<Size>& path, const std::vector<double>& pro_masses) const;
  double getSpecMassSpan_(const std::vector<Size>& path, const MSSpectrum& node_spec, int pro_mass_size) const;
  double getProteinMassSpan_(const std::vector<Size>& path, const std::vector<double>& pro_masses) const;
  int getModifiedAACount_(const std::vector<Size>& path) const;

  std::vector<std::string> ion_types_str_;
  std::vector<double> prefix_shifts_;
  std::vector<double> suffix_shifts_;
  std::vector<ProteinHit> proteoform_hits_;
  double tol_;

  int max_mod_cntr_ = 0;
  int allowed_isotope_error_ = 1;
  const int max_path_score_ = 1200;
  const int min_path_score_ = -20;
  const int max_extension_stretch_ = 50;
  double max_mod_mass_ = 500.0;
  double given_precursor_mass_ = -1;
};
} // namespace OpenMS