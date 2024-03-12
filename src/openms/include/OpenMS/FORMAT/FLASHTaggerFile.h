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

namespace OpenMS
{
  class OPENMS_DLLAPI FLASHTaggerFile
  {
  public:
    /// write header line for regular file output
    static void writeTagHeader(std::fstream& fs);

    /// write header line for topFD feature file
    static void writeProteinHeader(std::fstream& fs);

    /// write the features in regular file output
    static void writeTags(const FLASHTaggerAlgorithm& tagger, std::fstream& fs);

    static void writeProteins(const FLASHTaggerAlgorithm& tagger, std::fstream& fs);
  };
} // namespace OpenMS