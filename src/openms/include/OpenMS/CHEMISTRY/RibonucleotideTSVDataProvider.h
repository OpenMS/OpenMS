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

    An optional column named @c baseloss_formula specifies the formula retained
    after nucleobase loss, including the sugar modification (e.g. @c C6H10O5
    for LNA or @c C5H9FO4 for 2'-fluorodeoxyribose). This is the hydrated
    base-loss formula used by Ribonucleotide::getBaselossFormula(), before the
    water losses needed to form an a-B ion. Append this column after the
    optional @c alternatives column, leaving alternatives empty for unambiguous
    modifications. An explicit formula overrides inference from the code;
    an absent, empty, or @c - value retains the legacy inference.

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
