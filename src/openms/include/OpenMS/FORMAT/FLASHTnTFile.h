// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------


#pragma once

#include <OpenMS/config.h>
#include <iomanip>
#include <iostream>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTAlgorithm.h>

namespace OpenMS
{
  class OPENMS_DLLAPI FLASHTnTFile
{
  public:
    /// write header line for tag output file
    static void writeTagHeader(std::fstream& fs);

    /// write header line for PrSM output file
    static void writePrSMHeader(std::fstream& fs);

    /// write header line for Proteoform output file
    static void writeProHeader(std::fstream& fs);

    /// write the features in regular file output
    static void writeTags(const FLASHTnTAlgorithm& tnt, double flanking_mass_tol, std::fstream& fs);

    static void writePrSMs(const std::vector<ProteinHit>& hits, std::fstream& fs);

    static String generateProFormaString_(const String& sequence, int seq_start, int seq_end,
                                  const std::vector<double>& mod_masses,
                                  const std::vector<int>& mod_starts,
                                  const std::vector<int>& mod_ends,
                                  const std::vector<String>& mod_ids);

    static void writeProteoforms(const std::vector<ProteinHit>& hits, std::fstream& fs, double pro_fdr);
  };
} // namespace OpenMS