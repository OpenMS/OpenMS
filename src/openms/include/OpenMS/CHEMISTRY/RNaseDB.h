// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/DigestionEnzymeDB.h>
#include <OpenMS/CHEMISTRY/DigestionEnzymeRNA.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <memory>
#include <vector>

namespace OpenMS
{
  /**
    @ingroup Chemistry

    @brief Database for enzymes that digest RNA (RNases)

    The enzymes stored in this DB are defined in an XML file under "share/CHEMISTRY/Enzymes_RNA.xml".

    Supports dependency injection via a provider-based constructor for testing
    and custom enzyme sources.
  */
  class OPENMS_DLLAPI RNaseDB: public DigestionEnzymeDB<DigestionEnzymeRNA, RNaseDB>
  {
    // allow access to constructor in DigestionEnzymeDB::getInstance():
    friend class DigestionEnzymeDB<DigestionEnzymeRNA, RNaseDB>;

  protected:
    /// default constructor: loads enzymes from XML file
    RNaseDB();

  public:
    /// @brief Construct from custom data providers (for testing / dependency injection)
    /// @param[in] providers Data providers supplying enzyme definitions
    explicit RNaseDB(std::vector<std::unique_ptr<DigestionEnzymeDataProvider<DigestionEnzymeRNA>>> providers);
  };
}
