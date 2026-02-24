// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/RibonucleotideDataProvider.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{
  /**
    @brief Loads Ribonucleotide data from a Modomics JSON file.

    Parses the JSON format used by the Modomics database
    (https://www.genesilico.pl/modomics/api/modifications) and returns
    ribonucleotide entries as RibonucleotideEntry objects.

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI ModomicsJSONDataProvider : public RibonucleotideDataProvider
  {
  public:
    /**
      @param filename Path to a Modomics JSON file (resolved via File::find)
    */
    explicit ModomicsJSONDataProvider(const String& filename);

    std::vector<RibonucleotideEntry> loadRibonucleotides() override;

  private:
    String filename_;
  };

} // namespace OpenMS
