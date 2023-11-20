// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeon $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <iomanip>
#include <iostream>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <boost/dynamic_bitset.hpp>

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

    ///Find sequence tags from @p mzs and @p intensities then store them in @p tags.
    void run(const std::vector<double>& mzs, const std::vector<int>& scores, double ppm, std::vector<FLASHDeconvHelperStructs::Tag>& tags, const std::function<int(int, int)>& edge_score);
    void run(const std::vector<double>& mzs, const std::vector<int>& scores, double ppm, std::vector<FLASHDeconvHelperStructs::Tag>& tags);
    void run(const DeconvolvedSpectrum& dspec, double ppm, std::vector<FLASHDeconvHelperStructs::Tag>& tags);


  protected:
    void updateMembers_() override;
    /// implemented for DefaultParamHandler
    void setDefaultParams_();

  private:
    class DAC_;

    void constructDAC_(TopDownTagger::DAC_& dac, const std::vector<double>& mzs, const std::vector<int>& scores, int z, double tol);
    std::vector<Residue> getAA_(double l, double r, double tol, int z) const;
    std::vector<std::vector<Residue>> getGap_(double l, double r, double tol, int z) const;
    void updateEdgeMasses_();
    int getVertex_(int index, int path_score, int level, int gap_level) const;
    int getIndex_(int vertex) const;
    void updateTagSet_(std::set<FLASHDeconvHelperStructs::Tag>& tag_set, const std::vector<int>& path, const std::vector<double>& mzs, int z, int score);

    static int edgeScore_(int vertex_score1, int vertex_score2);
    bool connectEdge_(TopDownTagger::DAC_& dac, int vertex1, int vertex2, boost::dynamic_bitset<>& visited);

    std::function<int(int, int)> edge_score_;

    std::set<const Residue*> aas_ = ResidueDB::getInstance()->getResidues("Natural20");
    std::map<double, std::vector<Residue>> aa_mass_map_;
    std::map<double, std::vector<std::vector<Residue>>> gap_mass_map_;

    std::map<int, std::map<int, std::vector<String>>> edge_aa_map_;

    int allowed_isotope_error_ = 0;
    int max_tag_count_ = 0;
    int tag_length_ = 0;
    int max_aa_in_gap_ = 0;
    int max_gap_count_ = 0;
    int max_path_score_ = 0;
    int min_path_score_ = 0;

    double max_edge_mass_ = 0;
    int min_charge_ = 0;
    int max_charge_ = 0;
  };
} // namespace OpenMS