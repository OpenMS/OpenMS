// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/NUXL/NuXLFragmentAdductDefinition.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLModificationsGenerator.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <vector>
#include <map>
#include <set>
#include <iostream>

namespace OpenMS
{  

// fast (flat) data structure to store feasible x-,y-,a-ion fragment adducts and observable marker ions
using NucleotideToFeasibleFragmentAdducts = std::pair<char, std::vector<NuXLFragmentAdductDefinition> >;

// stores the fragment adducts and marker ions for a given precursor adduct
struct OPENMS_DLLAPI MS2AdductsOfSinglePrecursorAdduct
{
  std::vector<NucleotideToFeasibleFragmentAdducts> feasible_adducts;
  std::vector<NuXLFragmentAdductDefinition> marker_ions;
};

// helper struct to facilitate parsing of parameters (modifications, nucleotide adducts, ...)
struct OPENMS_DLLAPI NuXLParameterParsing
{
  /// Query ResidueModifications (given as strings) from ModificationsDB
  static std::vector<ResidueModification> getModifications(StringList modNames);

  // Map a nucleotide (e.g. U to all possible fragment adducts)
  using NucleotideToFragmentAdductMap = std::map<char, std::set<NuXLFragmentAdductDefinition> >;
  // @brief Parse tool parameter to create map from target nucleotide to all its fragment adducts
  // It maps a single letter nucleotide (e.g., 'T', 'C', ...)
  // to the maximum set of fragment adducts that may arise if the nucleotide is cross-linked.
  // Losses, that might reduce this set, are not considered in this data structure and handled later
  // when specific precursor adducts are considered.
  static NucleotideToFragmentAdductMap getTargetNucleotideToFragmentAdducts(StringList fragment_adducts);

  // @brief Determines the fragment adducts and marker ions for a given precursor.
  // The precursor adduct (the oligo including losses, e.g.: "TC-H3PO4") is mapped to all contained nucleotides
  // and their marker ions. In addition, each cross-linkable nucleotide is mapped to its chemically feasible fragment adducts.
  // Chemical feasible means in this context, that the fragment or marker ion adduct is a subformula of the precursor adduct.
  static MS2AdductsOfSinglePrecursorAdduct getFeasibleFragmentAdducts(
    const String& exp_pc_adduct,
    const String& exp_pc_formula,
    const NucleotideToFragmentAdductMap& nucleotide_to_fragment_adducts,
    const std::set<char>& can_xl,
    const bool always_add_default_marker_ions,
    const bool default_marker_ions_RNA
  );

  // Maps a precursor adduct (e.g.: "UU-H2O") to all chemically feasible fragment adducts.
  using PrecursorsToMS2Adducts = std::map<std::string, MS2AdductsOfSinglePrecursorAdduct>;

  // @brief extract all marker ions into a vector and make it unique according to mass. (e.g., used for matching agains all possible marker ions for MIC calculation)
  static std::vector<NuXLFragmentAdductDefinition> getMarkerIonsMassSet(const PrecursorsToMS2Adducts& pc2adducts);

  // @brief Calculate all chemically feasible fragment adducts for all possible precursor adducts
  // Same as getFeasibleFragmentAdducts but calculated from all precursor adducts
  static PrecursorsToMS2Adducts getAllFeasibleFragmentAdducts(
    const NuXLModificationMassesResult& precursor_adducts,
    const NucleotideToFragmentAdductMap& nucleotide_to_fragment_adducts,
    const std::set<char>& can_xl,
    const bool always_add_default_marker_ions,
    const bool default_marker_ions_RNA
  );
};

}

