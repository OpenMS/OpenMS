// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/ID/DBSearchParam.h>

#include <tuple>

namespace OpenMS::IdentificationDataInternal
{
  bool DBSearchParam::operator<(const DBSearchParam& other) const
  {
    return (std::tie(molecule_type, mass_type,
                     database, database_version, taxonomy,
                     charges, fixed_mods, variable_mods,
                     precursor_mass_tolerance, fragment_mass_tolerance,
                     precursor_tolerance_ppm, fragment_tolerance_ppm,
                     digestion_enzyme, enzyme_term_specificity, missed_cleavages,
                     min_length, max_length) <
            std::tie(other.molecule_type, other.mass_type,
                     other.database, other.database_version, other.taxonomy,
                     other.charges, other.fixed_mods, other.variable_mods,
                     other.precursor_mass_tolerance, other.fragment_mass_tolerance,
                     other.precursor_tolerance_ppm, other.fragment_tolerance_ppm,
                     other.digestion_enzyme, other.enzyme_term_specificity, other.missed_cleavages,
                     other.min_length, other.max_length));
  }

  bool DBSearchParam::operator==(const DBSearchParam& other) const
  {
    return (std::tie(molecule_type, mass_type, database,
                     database_version, taxonomy, charges, fixed_mods,
                     variable_mods, fragment_mass_tolerance,
                     precursor_mass_tolerance, fragment_tolerance_ppm,
                     precursor_tolerance_ppm, digestion_enzyme, enzyme_term_specificity,
                     missed_cleavages, min_length, max_length) ==
            std::tie(other.molecule_type, other.mass_type,
                     other.database, other.database_version, other.taxonomy,
                     other.charges, other.fixed_mods, other.variable_mods,
                     other.fragment_mass_tolerance,
                     other.precursor_mass_tolerance,
                     other.fragment_tolerance_ppm,
                     other.precursor_tolerance_ppm,
                     other.digestion_enzyme, other.enzyme_term_specificity,
                     other.missed_cleavages,
                     other.min_length, other.max_length));
  }
} // namespace OpenMS::IdentificationDataInternal
