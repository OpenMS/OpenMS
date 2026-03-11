// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Xiao Liang $
// $Authors: Xiao Liang, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/DigestionEnzymeDB.h>
#include <OpenMS/CHEMISTRY/DigestionEnzymeProtein.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <memory>
#include <vector>

namespace OpenMS
{
  /**
    @ingroup Chemistry

    @brief Database for enzymes that digest proteins (proteases)

    The enzymes stored in this DB are defined as built-in defaults. Additional
    user-defined enzymes can be loaded from share/CHEMISTRY/Enzymes.xml if present.

    Supports dependency injection via a provider-based constructor for testing
    and custom enzyme sources.
  */
  class OPENMS_DLLAPI ProteaseDB: public DigestionEnzymeDB<DigestionEnzymeProtein, ProteaseDB>
  {
    // allow access to constructor in DigestionEnzymeDB::getInstance():
    friend class DigestionEnzymeDB<DigestionEnzymeProtein, ProteaseDB>;

  protected:
    /// default constructor: loads built-in enzymes and optional XML file
    ProteaseDB();

  public:
    /// @brief Construct from custom data providers (for testing / dependency injection)
    /// @param[in] providers Data providers supplying enzyme definitions
    explicit ProteaseDB(std::vector<std::unique_ptr<DigestionEnzymeDataProvider<DigestionEnzymeProtein>>> providers);

    /// returns all the enzyme names available for XTandem
    void getAllXTandemNames(std::vector<String>& all_names) const;

    /// returns all the enzyme names available for Comet
    void getAllCometNames(std::vector<String>& all_names) const;

     /// returns all the enzyme names available for OMSSA
    void getAllOMSSANames(std::vector<String>& all_names) const;

    /// returns all the enzyme names available for MSGFPlus
    void getAllMSGFNames(std::vector<String>& all_names) const;

    /// writes the full names to a TSV file
    void writeTSV(const String& filename);
  };
}
