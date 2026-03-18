// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CONCEPT/HashUtils.h>
#include <functional>
#include <vector>
#include <map>
#include <set>

namespace OpenMS
{  
struct OPENMS_DLLAPI NuXLFragmentAdductDefinition
{
  EmpiricalFormula formula; // formula
  String name;  // name used in annotation
  double mass = 0;

  NuXLFragmentAdductDefinition() = default;

  NuXLFragmentAdductDefinition(const NuXLFragmentAdductDefinition&) = default;

  NuXLFragmentAdductDefinition(NuXLFragmentAdductDefinition&&) = default;

  NuXLFragmentAdductDefinition(const EmpiricalFormula& f, const String& n, double m) : formula(f), name(n), mass(m) {}

  NuXLFragmentAdductDefinition& operator=(const NuXLFragmentAdductDefinition&) = default;

  NuXLFragmentAdductDefinition& operator=(NuXLFragmentAdductDefinition&&) = default;

  bool operator<(const NuXLFragmentAdductDefinition& other) const;

  bool operator==(const NuXLFragmentAdductDefinition& other) const;

};

} // namespace OpenMS

// Hash function specialization for NuXLFragmentAdductDefinition
// Placed in std namespace to allow use with std::unordered_map/set
namespace std
{
  /**
   * @brief Hash function for OpenMS::NuXLFragmentAdductDefinition.
   *
   * Computes a hash based on all fields used in operator==:
   * formula (EmpiricalFormula) and name (String).
   *
   * @note Hash is consistent with operator==.
   */
  template<>
  struct hash<OpenMS::NuXLFragmentAdductDefinition>
  {
    std::size_t operator()(const OpenMS::NuXLFragmentAdductDefinition& fad) const noexcept
    {
      std::size_t seed = 0;
      // Hash formula using EmpiricalFormula's std::hash specialization
      OpenMS::hash_combine(seed, std::hash<OpenMS::EmpiricalFormula>{}(fad.formula));
      // Hash name using fnv1a_hash_string (String inherits from std::string)
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(fad.name));
      return seed;
    }
  };
} // namespace std

