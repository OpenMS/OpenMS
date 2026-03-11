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
    @brief Loads ResidueModification data from an OBO ontology file (PSI-MOD or XLMOD).

    Parses OBO format files and returns modifications as unique_ptr.
    The @p cross_links_only flag controls filtering:
    - false (default): returns all modifications EXCEPT cross-linkers (reactionSites==2).
      Used by ModificationsDB for PSI-MOD.obo and XLMOD.obo (mono-links only).
    - true: returns only cross-linkers (reactionSites==2), skipping mono-links.
      Used by CrossLinksDB for XLMOD.obo.

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI OBODataProvider : public ModificationDataProvider
  {
  public:
    /**
      @param filename Path to an OBO file (resolved via File::find)
      @param cross_links_only If true, return only cross-linkers; if false, exclude them
    */
    explicit OBODataProvider(const String& filename, bool cross_links_only = false);

    std::vector<std::unique_ptr<ResidueModification>> loadModifications() override;

  private:
    String filename_;
    bool cross_links_only_;
  };

} // namespace OpenMS
