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
#include <vector>
#include <map>
#include <set>
#include <iostream>

namespace OpenMS
{
  class AASequence;

  /**
    @brief Result of @ref NuXLModificationsGenerator::initModificationMassesNA — empirical formulas with their monoisotopic mass and the disambiguating nucleotide-composition lookup.

    Two parallel maps keyed by the same empirical-formula string:
      - @c formula2mass — empirical formula → monoisotopic mass.
      - @c mod_combinations — empirical formula → set of nucleotide-composition strings that
        produce that formula. A set rather than a single value because different oligos can
        share a chemical formula, e.g. @c C18H22N4O16P2 → @c {CU-H3N1, UU-H2O1}.

    The nested @ref MyStringLengthCompare comparator orders the composition strings primarily
    by length and lexicographically within the same length so the shortest (most likely)
    composition appears first when iterating a value.

    @ingroup Analysis_ID
  */
  struct OPENMS_DLLAPI NuXLModificationMassesResult
  {
    /// Comparator that orders strings primarily by length (shortest first) and lexicographically within the same length.
    struct MyStringLengthCompare
    {
      bool operator () (const std::string & p_lhs, const std::string & p_rhs) const
      {
        const size_t lhsLength = p_lhs.length() ;
        const size_t rhsLength = p_rhs.length() ;
        if(lhsLength == rhsLength)
        {
            return (p_lhs < p_rhs) ; // when two strings have the same
                                     // length, defaults to the normal
                                     // string comparison
        }
        return (lhsLength < rhsLength) ; // compares with the length
      }
    };
    std::map<String, double> formula2mass; ///< Empirical formula → monoisotopic mass.

    /// Length-ordered set of nucleotide-composition strings producing the same empirical formula.
    using NucleotideFormulas = std::set<String, MyStringLengthCompare>;
    /// Empirical formula → @ref NucleotideFormulas.
    using MapSumFormulaToNucleotideFormulas = std::map<String, NucleotideFormulas>;
    MapSumFormulaToNucleotideFormulas mod_combinations; ///< Empirical formula → all nucleotide compositions yielding that formula (kept as a set because mutations / losses can make the mapping ambiguous).
  };

  /**
    @brief Enumerator of precursor-adduct masses for NuXL cross-link searches.

    Builds the table of empirical formulas + monoisotopic masses + disambiguating
    nucleotide compositions that NuXL searches against. Driven by:
      - the list of target nucleotides with their monophosphate empirical formulas
        (e.g. @c "U=C9H13N2O9P" entries in @p target_nucleotides),
      - the list of source→target mappings used to express modifications
        (e.g. @c "C->T" entries in @p mappings),
      - the list of @c "nucleotide:+formula-formula" modifications applied to each
        target nucleotide (e.g. @c "U:+H2O-H2O" entries in @p modifications),
      - and an optional @p sequence_restriction limiting which oligomer compositions
        survive (only adducts whose nucleotide string is a substring of the restriction
        sequence, sort-insensitive at each window, are kept).

    Optional @p cysteine_adduct toggles a hardcoded DTT-derived @c C4H8S2O2 entry
    (commonly referred to as the @em 152 modification) into the result.

    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI NuXLModificationsGenerator
  {
    public:
      /**
        @brief Build the full set of precursor adducts (formula + mass + nucleotide composition) for the given configuration.

        Iterates source-nucleotide combinations up to @p max_length, applies the
        per-nucleotide modifications, and (optionally) keeps only those compositions that
        are substrings of @p sequence_restriction. When @p sequence_restriction is empty,
        every length-@c k combination of source nucleotides for @c k @c in @c [1, max_length]
        is considered.

        @param[in] target_nucleotides    Entries in the form @c "N=Empirical_formula"
                                         giving the monophosphate formula for each target
                                         nucleotide (e.g. @c "U=C9H13N2O9P"). The entries
                                         are split on @c '='.
        @param[in] nt_groups             Group identifiers consulted during the
                                         combinatorial expansion (forwarded verbatim into
                                         the per-modification application loop).
        @param[in] can_xl                Set of nucleotides eligible to carry a cross-link;
                                         compositions without at least one cross-linkable
                                         nucleotide are dropped.
        @param[in] mappings              Source→target nucleotide mappings in the form
                                         @c "X->Y" (split on @c "->"). Used to expand
                                         source-nucleotide sequences into target sequences
                                         and to derive the @em source alphabet for the
                                         combinatorial enumeration.
        @param[in] modifications         Per-nucleotide modification descriptors in the
                                         form @c "N:+formula-formula" (e.g.
                                         @c "U:+H2O-H2O"). The second character must be
                                         @c ':' — otherwise the call throws.
        @param[in] sequence_restriction  Optional NA reference sequence; when set, only
                                         oligo compositions that appear as a length-matched
                                         windowed permutation of the reference survive.
                                         Default empty (no restriction → all combinations
                                         up to @p max_length).
        @param[in] cysteine_adduct       If @c true, insert the hardcoded DTT-derived
                                         @c C4H8S2O2 adduct (the @em 152 modification) into
                                         the result.
        @param[in] max_length            Maximum oligomer length to enumerate when
                                         @p sequence_restriction is empty. Default @c 4.

        @return Populated @ref NuXLModificationMassesResult.

        @throws OpenMS::Exception::MissingInformation when a @p modifications entry does
                not follow the @c "N:+formula-formula" format (specifically when the
                second character of the entry is not @c ':').
      */
      static NuXLModificationMassesResult initModificationMassesNA(const StringList& target_nucleotides,
                                                                     const StringList& nt_groups,
                                                                     const std::set<char>& can_xl,
                                                                     const StringList& mappings,
                                                                     const StringList& modifications,
                                                                     String sequence_restriction = "",
                                                                     bool cysteine_adduct = false,
                                                                     Int max_length = 4);
    private:
      /**
        @brief Test whether @p query is @em not present as a sorted-window permutation of @p res_seq.

        Returns @c true iff no contiguous window of length @c query.size() in @p res_seq
        has the same multi-set of characters as @p query. The window comparison sorts both
        sides first, so the test is permutation-insensitive. Empty @p query short-circuits
        to @c false (treated as always present).
      */
      static bool notInSeq(const String& res_seq, const String& query);

      /**
        @brief Recursively expand @p res_seq into every target-substituted variant according to @p map_source2target.

        Starting at @p param_pos, walks @p res_seq character by character. Whenever the
        current character is a key in @p map_source2target with more than one candidate
        target, recurses once per candidate that differs from the current character; the
        unchanged path falls through. After the recursive walk, the resulting sequence is
        appended to @p target_sequences only when every position is either a non-mapped
        character or a character that is simultaneously a source and a valid target
        nucleotide (i.e. no @em pure source nucleotide is left). @p target_sequences is
        appended to; existing entries are preserved.
      */
      static void generateTargetSequences(const String& res_seq, Size param_pos, const std::map<char, std::vector<char> >& map_source2target, StringList& target_sequences);
    };
}

