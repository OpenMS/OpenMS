// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/Ribonucleotide.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OpenMSConfig.h>

#include <memory>
#include <vector>

namespace OpenMS
{
  /**
    @brief Holds a Ribonucleotide together with optional ambiguity codes.

    When a ribonucleotide represents an ambiguous modification (e.g., a methyl
    group whose localization cannot be determined), the two alternative codes
    identify the concrete modifications it could represent.

    @ingroup Chemistry
  */
  struct OPENMS_DLLAPI RibonucleotideEntry
  {
    std::unique_ptr<Ribonucleotide> ribo;
    String alternative_code_1; ///< code of first alternative (empty if unambiguous)
    String alternative_code_2; ///< code of second alternative (empty if unambiguous)

    RibonucleotideEntry() = default;
    ~RibonucleotideEntry() = default;
    RibonucleotideEntry(const RibonucleotideEntry&) = delete;
    RibonucleotideEntry& operator=(const RibonucleotideEntry&) = delete;
    RibonucleotideEntry(RibonucleotideEntry&&) = default;
    RibonucleotideEntry& operator=(RibonucleotideEntry&&) = default;

    /// Returns true if this entry represents an ambiguous modification
    bool isAmbiguous() const { return !alternative_code_1.empty(); }
  };

  /**
    @brief Interface for providing Ribonucleotide data to RibonucleotideDB.

    Implementations of this interface abstract the source of ribonucleotide data,
    enabling dependency injection. File-based providers (ModomicsJSONDataProvider,
    RibonucleotideTSVDataProvider) handle I/O; InMemoryRibonucleotideDataProvider
    supports testing.

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI RibonucleotideDataProvider
  {
  public:
    virtual ~RibonucleotideDataProvider() = default;

    /**
      @brief Load ribonucleotides from whatever source this provider wraps.
      @return Vector of ribonucleotide entries with ownership transferred to caller.
      @note Providers may only be called once. Subsequent calls may return empty results.
    */
    virtual std::vector<RibonucleotideEntry> loadRibonucleotides() = 0;
  };

  /**
    @brief Data provider that serves pre-built ribonucleotides from memory.

    Useful for unit testing RibonucleotideDB without requiring files on disk.

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI InMemoryRibonucleotideDataProvider : public RibonucleotideDataProvider
  {
  public:
    explicit InMemoryRibonucleotideDataProvider(std::vector<RibonucleotideEntry> entries)
      : entries_(std::move(entries))
    {
    }

    InMemoryRibonucleotideDataProvider(const InMemoryRibonucleotideDataProvider&) = delete;
    InMemoryRibonucleotideDataProvider& operator=(const InMemoryRibonucleotideDataProvider&) = delete;
    InMemoryRibonucleotideDataProvider(InMemoryRibonucleotideDataProvider&&) = default;
    InMemoryRibonucleotideDataProvider& operator=(InMemoryRibonucleotideDataProvider&&) = default;

    std::vector<RibonucleotideEntry> loadRibonucleotides() override
    {
      return std::move(entries_);
    }

  private:
    std::vector<RibonucleotideEntry> entries_;
  };

} // namespace OpenMS
