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
#include <OpenMS/METADATA/ID/ProcessingStep.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

#include <map>
#include <set>
#include <string>
#include <tuple>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /** @brief Parameters specific to a database search step.
    */
    struct OPENMS_DLLAPI DBSearchParam: public MetaInfoInterface
    {
      enum MoleculeType molecule_type;
      enum MassType mass_type;

      std::string database;
      std::string database_version;
      std::string taxonomy;

      std::set<Int> charges;

      std::set<std::string> fixed_mods;
      std::set<std::string> variable_mods;

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

      /// Lexicographic order over all parameters (defined out of line: comparing the
      /// 17-member tuple is expensive to instantiate in every including translation unit)
      bool operator<(const DBSearchParam& other) const;

      bool operator==(const DBSearchParam& other) const;
    };

    typedef std::set<DBSearchParam> DBSearchParams;
    typedef IteratorWrapper<DBSearchParams::iterator> SearchParamRef;
    typedef std::map<ProcessingStepRef, SearchParamRef> DBSearchSteps;

  }
}
