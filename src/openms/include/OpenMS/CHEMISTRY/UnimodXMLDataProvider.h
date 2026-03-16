// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/ModificationDataProvider.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{
  /**
    @brief Loads ResidueModification data from a Unimod XML file.

    Wraps UnimodXMLFile::load() behind the ModificationDataProvider interface.
    Each modification gets its FullId set after parsing.

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI UnimodXMLDataProvider : public ModificationDataProvider
  {
  public:
    /// @param filename Path to a Unimod XML file (resolved via File::find)
    explicit UnimodXMLDataProvider(const String& filename);

    std::vector<std::unique_ptr<ResidueModification>> loadModifications() override;

  private:
    String filename_;
  };

} // namespace OpenMS
