// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeon $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#pragma once

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
  @brief
  @ingroup Topdown
  */

  class OPENMS_DLLAPI TopDownTagger : public DefaultParamHandler
  {
  public:
    /// constructor
    TopDownTagger();

    /// destructor
    ~TopDownTagger() override = default;

    /// copy constructor
    TopDownTagger(const TopDownTagger&);

    /// move constructor
    TopDownTagger(TopDownTagger&& other) = default;

    /// assignment operator
    TopDownTagger& operator=(const TopDownTagger& other);

    /// Find sequence tags from @p mzs and @p intensities then store them in @p tags.
    void run(const std::vector<double>& mzs, const std::vector<int>& scores, double ppm, const std::function<int(int, int)>& edge_score);
    void run(const std::vector<double>& mzs, const std::vector<int>& scores, double ppm);
    void run(const DeconvolvedSpectrum& dspec, double ppm);
    void runMatching(const String& fasta_file);

    const std::vector<ProteinHit>& getProteinHits() const;
    const std::vector<ProteinHit> getProteinHits(const FLASHDeconvHelperStructs::Tag& tag) const;
    std::vector<FLASHDeconvHelperStructs::Tag> getTags(const ProteinHit& hit) const;
    const std::vector<FLASHDeconvHelperStructs::Tag>& getTags() const;
    int getProteinIndex(const ProteinHit& hit) const;
    int getTagIndex(const FLASHDeconvHelperStructs::Tag& tag) const;

  protected:
    void updateMembers_() override;
    /// implemented for DefaultParamHandler
    void setDefaultParams_();

  private:
    class DAC_;

    void constructDAC_(TopDownTagger::DAC_& dac, const std::vector<double>& mzs, const std::vector<int>& scores, int length, double tol);
    std::vector<Residue> getAA_(double l, double r, double tol, int iso_offset = 0) const;
    // std::vector<std::vector<Residue>> getGap_(double l, double r, double tol, int z) const;
    void updateEdgeMasses_();
    int getVertex_(int index, int path_score, int level, int iso_level) const;
    int getIndex_(int vertex) const;
    void updateTagSet_(std::set<FLASHDeconvHelperStructs::Tag>& tag_set, std::map<String, std::vector<FLASHDeconvHelperStructs::Tag>>& seq_tag, const std::vector<int>& path, const std::vector<double>& mzs, const std::vector<int>& scores, double ppm);

    static int edgeScore_(int vertex_score1, int vertex_score2);
    bool connectEdge_(TopDownTagger::DAC_& dac, int vertex1, int vertex2, boost::dynamic_bitset<>& visited);

    static Size find_with_X_(const String& A, const String& B);

    std::function<int(int, int)> edge_score_;

    std::set<const Residue*> aas_ = ResidueDB::getInstance()->getResidues("Natural20");
    std::map<double, std::vector<Residue>> aa_mass_map_;
    std::map<int, std::map<int, std::vector<String>>> edge_aa_map_;

    std::vector<FLASHDeconvHelperStructs::Tag> tags_;
    std::vector<ProteinHit> protein_hits_;
    std::vector<std::vector<int>> matching_tags_indices_; // from protein hit to tag index
    std::vector<std::vector<int>> matching_hits_indices_; // from tag to protein hit index

    int max_tag_count_ = 0;
    int min_tag_length_ = 0;
    int max_tag_length_ = 0;
    int max_iso_in_tag_ = 0;
    int max_path_score_ = 0;
    int min_path_score_ = 0;
    int min_cov_aa_ = 5;
    double fdr_ = 1.0;
    bool keep_decoy_ = false;
    double max_edge_mass_ = 0;
    double flanking_mass_tol_ = 200.0;
  };
} // namespace OpenMS