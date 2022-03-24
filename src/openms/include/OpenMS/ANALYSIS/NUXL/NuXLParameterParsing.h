// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

