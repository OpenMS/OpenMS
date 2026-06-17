// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/OpenMSConfig.h>

#include <memory>
#include <vector>

namespace OpenMS
{
  /**
    @brief Interface for providing ResidueModification data to ModificationsDB.

    Implementations of this interface abstract the source of modification data,
    enabling dependency injection. File-based providers (UnimodXMLDataProvider,
    OBODataProvider) handle I/O; InMemoryDataProvider supports testing.

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI ModificationDataProvider
  {
  public:
    virtual ~ModificationDataProvider() = default;

    /**
      @brief Load modifications from whatever source this provider wraps.
      @return Vector of modifications with ownership transferred to caller.
      @note Providers may only be called once. Subsequent calls may return empty results.
    */
    virtual std::vector<std::unique_ptr<ResidueModification>> loadModifications() = 0;
  };

  /**
    @brief Data provider that serves pre-built modifications from memory.

    Useful for unit testing ModificationsDB without requiring files on disk.

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI InMemoryDataProvider : public ModificationDataProvider
  {
  public:
    explicit InMemoryDataProvider(std::vector<std::unique_ptr<ResidueModification>> mods)
      : mods_(std::move(mods))
    {
    }

    InMemoryDataProvider(const InMemoryDataProvider&) = delete;
    InMemoryDataProvider& operator=(const InMemoryDataProvider&) = delete;
    InMemoryDataProvider(InMemoryDataProvider&&) = default;
    InMemoryDataProvider& operator=(InMemoryDataProvider&&) = default;

    std::vector<std::unique_ptr<ResidueModification>> loadModifications() override
    {
      return std::move(mods_);
    }

  private:
    std::vector<std::unique_ptr<ResidueModification>> mods_;
  };

} // namespace OpenMS
