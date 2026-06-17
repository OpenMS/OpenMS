// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Sachsenberg $
// $Authors: Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <memory>
#include <vector>

namespace OpenMS
{
  /**
    @ingroup Chemistry

    @brief Abstract interface for providing digestion enzyme data

    Implementations of this interface supply enzyme definitions to DigestionEnzymeDB
    subclasses. This decouples data sourcing (file I/O, built-in definitions, etc.)
    from the database indexing logic.

    Template parameter @p EnzymeType should be a subclass of DigestionEnzyme
    (e.g. DigestionEnzymeProtein, DigestionEnzymeRNA).
  */
  template <typename EnzymeType>
  class DigestionEnzymeDataProvider
  {
  public:
    virtual ~DigestionEnzymeDataProvider() = default;

    /// @brief Load enzyme definitions. Ownership is transferred to the caller.
    /// @note Providers may only support a single call. Callers should invoke loadEnzymes() exactly once per provider instance.
    virtual std::vector<std::unique_ptr<EnzymeType>> loadEnzymes() = 0;
  };

  /**
    @ingroup Chemistry

    @brief In-memory data provider for digestion enzymes (primarily for testing)

    Accepts pre-built enzyme objects and returns them via loadEnzymes().
  */
  template <typename EnzymeType>
  class InMemoryDigestionEnzymeDataProvider : public DigestionEnzymeDataProvider<EnzymeType>
  {
  public:
    /// @brief Construct from a vector of pre-built enzymes (takes ownership)
    /// @param[in] enzymes Pre-built enzyme objects; ownership is transferred to this provider
    explicit InMemoryDigestionEnzymeDataProvider(std::vector<std::unique_ptr<EnzymeType>> enzymes)
      : enzymes_(std::move(enzymes))
    {
    }

    /// Load enzymes. This moves from internal storage; subsequent calls return an empty vector.
    std::vector<std::unique_ptr<EnzymeType>> loadEnzymes() override
    {
      return std::move(enzymes_);
    }

  private:
    std::vector<std::unique_ptr<EnzymeType>> enzymes_;
  };

} // namespace OpenMS
