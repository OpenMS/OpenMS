// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/METADATA/ID/MetaData.h>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /** @brief Parameters specific to a database search step.
    */
    struct DBSearchParam: public MetaInfoInterface
    {
      enum MoleculeType molecule_type;
      enum MassType mass_type;

      String database;
      String database_version;
      String taxonomy;

      std::set<Int> charges;

      std::set<String> fixed_mods;
      std::set<String> variable_mods;

      double precursor_mass_tolerance;
      double fragment_mass_tolerance;
      bool precursor_tolerance_ppm;
      bool fragment_tolerance_ppm;

      // allow for either "DigestionEnzymeProtein" or "DigestionEnzymeRNA":
      const DigestionEnzyme* digestion_enzyme;
      EnzymaticDigestion::Specificity enzyme_term_specificity;
      Size missed_cleavages;
      Size min_length;
      Size max_length;

      DBSearchParam():
          molecule_type(MoleculeType::PROTEIN),
          mass_type(MassType::MONOISOTOPIC),
          precursor_mass_tolerance(0.0), fragment_mass_tolerance(0.0),
          precursor_tolerance_ppm(false), fragment_tolerance_ppm(false),
          digestion_enzyme(nullptr), enzyme_term_specificity(EnzymaticDigestion::SPEC_UNKNOWN),
          missed_cleavages(0), min_length(0), max_length(0)
      {
      }

      DBSearchParam(const DBSearchParam& other) = default;

      bool operator<(const DBSearchParam& other) const
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

      bool operator==(const DBSearchParam& other) const
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
    };

    typedef std::set<DBSearchParam> DBSearchParams;
    typedef IteratorWrapper<DBSearchParams::iterator> SearchParamRef;
    typedef std::map<ProcessingStepRef, SearchParamRef> DBSearchSteps;

  }
}
