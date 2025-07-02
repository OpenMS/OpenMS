// Copyright (c) 2002-2024, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <boost/dynamic_bitset.hpp>
#include <iomanip>
#include <iostream>

namespace OpenMS
{
/**
@brief De novo sequence tag finding algorithm for Top Down proteomics. The sequence tags are
 generated on deconvolved spectrum (DeconvolvedSpectrum instance) quickly in the descending order of
 scores. The typical length of a tag ranges from 3 to 5 (user specified) and the tags are used to
 filter out protein candidates from the input fasta entries.
@ingroup Topdown
*/

class OPENMS_DLLAPI FLASHTaggerAlgorithm : public DefaultParamHandler, public ProgressLogger
{
public:
  /// constructor
  FLASHTaggerAlgorithm();

  /// destructor
  ~FLASHTaggerAlgorithm() override = default;

  /// copy constructor
  FLASHTaggerAlgorithm(const FLASHTaggerAlgorithm&);

  /// move constructor
  FLASHTaggerAlgorithm(FLASHTaggerAlgorithm&& other) = default;

  /// assignment operator
  FLASHTaggerAlgorithm& operator=(const FLASHTaggerAlgorithm& other);

  /**
    @brief Generate the tags from the input deconvolved spectrum with given ppm tolerance
    @param deconvolved_spectrum deconvolved spectrum from FLASHDeconv
    @param ppm The acceptable ppm tolerance for mass

  */
  void run(const DeconvolvedSpectrum& deconvolved_spectrum, double ppm);

  /**
   *@brief  Match the tags against protein sequences in the input fasta entry.
   * The maximum modification mass is used to skip protein sequences that do not match with
   * tag flanking masses.
   * @param fasta_entry fasta entries from the input proteome database in fasta format
   * @param deconvolved_spectrum deconvolved spectrum from FLASHDeconv
   * @param max_mod_mass maximum modification mass (a positive number)
   */
  void runMatching(const std::vector<FASTAFile::FASTAEntry>& fasta_entry,
                   const DeconvolvedSpectrum& deconvolved_spectrum,
                   double max_mod_mass = 0);

  /**
   * @brief fill tags with the length of @p tag_length in @p tags
   * @param tags
   * @param tag_length
   */
  void fillTags(std::vector<FLASHHelperClasses::Tag>& tags, int tag_length = 0) const;

  /// get the node score from peakgroup
  static int getNodeScore(const PeakGroup& peak_group);


  /**
  * Fill matched protein sequence positions and corresponding flankiing masses
  * for the input protein sequence (matched by tags)
  * @param positions positions of the protein sequence matches
  * @param masses flanking masses of the protein sequence matches
  * @param flanking_mass_tol up to how large mass flanking mass is allowed?
  * @param seq the input protein sequence
  * @param tag the tags to be matched against sequence
  */
  static void fillMatchedPositionsAndFlankingMassDiffs(std::vector<int>& positions,
                                                      std::vector<double>& masses,
                                                      double flanking_mass_tol,
                                                      const String& seq,
                                                      const FLASHHelperClasses::Tag& tag);

  /**
   * fill the protein hits, up to @p max_target_count hits. The hits are selected
   * from the high to low scoring ones.
   * After collecting max_target_count hits, if more hits are found with the same score,
   * they are also filled.
   * @param hits
   * @param max_target_count
   */
  void fillProteinHits(std::vector<ProteinHit>& hits, int max_target_count) const;

  /// maximum node score for tag generation and extension
  const static int max_node_score = 8;

protected:
  void updateMembers_() override;
  /// implemented for DefaultParamHandler
  void setDefaultParams_();

private:
  /**
    @brief generate sequence tags on the deconvolved spectrum.
    the generated tags are stored internally.
    @param dspec deconvolved spectrum
    @param ppm The acceptable ppm tolerance for mass.
    @param mode process mode : 0 - process common n c ion shift. 1 - n 2 - c term exclusive ion shift
  */

  void generateTags_(const DeconvolvedSpectrum& dspec, double ppm, int mode);
  /**
    @brief makes three vectors containing monoisotopic mass and score and scan numbers of each peakgroups.
    @param mzs mass
    @param scores score of each mass
    @param scan scan number
    @param ppm The acceptable ppm tolerance for mass.
    @param mode process mode : 0 - process common n c ion shift. 1 - n 2 - c term exclusive ion shift
  */
  void getTags_(const std::vector<double>& mzs, const std::vector<int>& scores, int scan, double ppm, int mode);
  void constructDAG_(FLASHHelperClasses::DAG& dag, const std::vector<double>& mzs, const std::vector<int>& scores, int length, double tol, int mode);
  std::vector<Residue> getAA_(double l, double r, double tol, int consider_ion_diff, int mode) const;
  std::vector<std::vector<Residue>> getGap_(double l, double r, double tol, int iso_offset) const;
  void updateEdgeMasses_();
  Size getVertex_(int index, int path_score, int level, int iso_level, int gap_level) const;
  int getIndex_(Size vertex) const;

  void getScoreAndMatchCount_(const std::vector<int>& spec_vec,
                              const boost::dynamic_bitset<>& pro_vec,
                              //const boost::dynamic_bitset<>& mask_pro_vec,
                              const std::set<int>& spec_pro_diffs,
                              const std::vector<int>& spec_scores,
                              int& max_score, int& match_cntr) const;


  void updateTagSet_(std::set<FLASHHelperClasses::Tag>& tag_set,
                     std::map<String, std::vector<FLASHHelperClasses::Tag>>& seq_tag,
                     const std::vector<Size>& path,
                     const std::vector<double>& mzs,
                     const std::vector<int>& scores,
                     int scan,
                     double ppm, int mode);

  //static std::vector<boost::dynamic_bitset<>> vectorized_fasta_entry_, rev_vectorized_fasta_entry_;
  //static std::vector<std::map<int, double>> mass_map_, rev_mass_map_;

  static void vectorizeFasta_(const std::vector<FASTAFile::FASTAEntry>& fasta_entry, bool reverse);

  static Size find_with_X_(const std::string_view& A, const String& B, Size pos = 0);

  std::set<const Residue*> aas_ = ResidueDB::getInstance()->getResidues("Natural20");
  std::map<double, std::vector<Residue>> aa_mass_map_;
  std::map<double, std::vector<std::vector<Residue>>> gap_mass_map_;
  std::map<int, std::map<int, std::vector<String>>> edge_aa_map_;

  std::vector<FLASHHelperClasses::Tag> tags_; // from scan to tags
  std::vector<ProteinHit> protein_hits_;

  std::set<double> common_shifts_;
  std::set<double> n_term_shifts_;
  std::set<double> c_term_shifts_;

  std::map<String, std::set<Size>> indexed_fasta_;
  bool consider_diff_ion_jumps_ = false;
  std::vector<Size> max_tag_counts_ {0, 0, 0, 100, 500, 1000}; // tag count for length 0, 1, 2, 3, 4, 5

  int min_tag_length_ = 0;
  int max_tag_length_ = 0;
  //int max_iso_in_tag_ = 0;
  int max_path_score_ = 0;
  int min_path_score_ = 0;
  int max_gap_count_ = 0;
  int max_aa_in_gap_ = 2;
  int min_cov_aa_ = 3;
  double max_edge_mass_ = 0;
};
} // namespace OpenMS