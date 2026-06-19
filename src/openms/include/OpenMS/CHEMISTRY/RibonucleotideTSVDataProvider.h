// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/RibonucleotideDataProvider.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

namespace OpenMS
{
  /**
    @brief Loads Ribonucleotide data from a tab-separated values (TSV) file.

    Parses the legacy TSV format used for custom RNA modifications and
    returns ribonucleotide entries as RibonucleotideEntry objects.

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI RibonucleotideTSVDataProvider : public RibonucleotideDataProvider
  {
  public:
    /**
      @param filename Path to a TSV file (resolved via File::find)
    */
    explicit RibonucleotideTSVDataProvider(const std::string& filename);

    std::vector<RibonucleotideEntry> loadRibonucleotides() override;

  private:
    std::string filename_;
  };

} // namespace OpenMS
