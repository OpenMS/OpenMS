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
@brief
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

  void run(std::vector<ProteinHit>& hits, const std::vector<FLASHHelperClasses::Tag>& tags,
           const DeconvolvedSpectrum& dspec,
           const MSSpectrum& spec, double ppm, bool multiple_hits_per_spec);

  void getProteoforms(std::vector<ProteinHit>& hits) const
  {
    for (const auto& hit : proteoform_hits_)
    {
      hits.push_back(hit);
    }
  }

  bool hasProteoforms() const
  {
    return ! proteoform_hits_.empty();
  }

  void setModificationMap(const std::map<double, std::vector<ResidueModification>>& mod_map)
  {
    mod_map_ = mod_map;
  }

  const static int multi_ion_score = 1;

protected:
  void updateMembers_() override;
  /// implemented for DefaultParamHandler
  void setDefaultParams_();

private:
  struct OPENMS_DLLAPI HitInformation
  {
  public:
    boost::dynamic_bitset<> visited_;
    FLASHHelperClasses::DAG dag_;
    std::map<int, MSSpectrum> node_spec_map_, tol_spec_map_;
    std::map<int, std::vector<double>> pro_mass_map_;
    int mode_;
    int protein_start_position_ = -1, protein_end_position_ = -1;
    double calculated_precursor_mass_ = -1;
    std::vector<int> start_pro_indices_;
  };

  std::map<double, std::vector<ResidueModification>> mod_map_; // modification mass to modification index. To use find nearest function

  static void getProMasses_(const ProteinHit& hit, std::vector<double>& pro_masses, int mode);

  void calculatePrecursorMass_(const ProteinHit& hit,
                               const std::map<int, std::vector<Size>>& best_path_map,
                               HitInformation& hi);

  void defineNodes_(const DeconvolvedSpectrum& dspec, HitInformation& hi,
                   double max_mass);
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
  std::vector<FLASHHelperClasses::Tag> tags_;
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