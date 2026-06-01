// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLParameterParsing.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLAnnotatedHit.h>

#include <set>
#include <map>
#include <vector>
#include <algorithm>

namespace OpenMS
{
/**
  @brief Post-scoring step of the OpenNuXL pipeline: annotates fragment-ion peaks of each candidate hit and localises the NA cross-link to a residue.

  Consumes the per-scan vectors of @ref NuXLAnnotatedHit produced upstream by the
  OpenNuXL search, reconstructs each hit's fully-modified peptide sequence, generates
  the theoretical fragment spectrum (b/a/y ions plus precursor + immonium variants),
  aligns it against the experimental spectrum, and writes the alignment results plus a
  per-residue localisation score back into the hit. Single-purpose helper class — only
  the static @ref annotateAndLocate_ entry point is exposed.

  @ingroup Analysis_ID
*/
class OPENMS_DLLAPI NuXLAnnotateAndLocate
{
  public:

  /**
    @brief Annotate every candidate hit in @p annotated_hits and localise the NA cross-link site.

    Iterates over @p annotated_hits in parallel (OpenMP @c parallel @c for over the scan
    dimension). For each non-empty @c annotated_hits[scan_index]:
      -# Reads the experimental @c PeakSpectrum and its first precursor charge from
         @p exp[scan_index].
      -# For each candidate hit @c a:
         - Rebuilds @c a.sequence's fully-modified form by reapplying
           @ref ModifiedPeptideGenerator::applyFixedModifications and
           @ref ModifiedPeptideGenerator::applyVariableModifications (the upstream search
           only stored modification indices to save memory) and selects the entry at
           @c a.peptide_mod_index.
         - Picks the precursor NA adduct from @p mm via @c a.NA_mod_index (using the
           first entry of the composition set at that key).
         - Generates a total-loss theoretical spectrum (b/a/y + precursor + manually
           added abundant immoniums for charge 1 and precursor peaks for every charge)
           via @ref TheoreticalSpectrumGenerator with hardcoded settings
           (@c add_first_prefix_ion, @c add_metainfo, @c add_a_ions, @c add_b_ions,
           @c add_y_ions enabled; @c add_c_ions, @c add_x_ions, @c add_z_ions,
           @c add_internal_fragments, @c add_abundant_immonium_ions,
           @c add_all_precursor_charges disabled — the latter two are added manually).
         - Aligns the theoretical and experimental spectra under @p fragment_mass_tolerance
           (in @em ppm when @p fragment_mass_tolerance_unit_ppm is @c true, else Da),
           with charge matching enforced.
         - Annotates the matched peaks as b/a/y ions, immoniums, marker ions and
           precursor variants (the immoniums / precursors are added by the NuXL fragment
           generator using the per-adduct definitions in @p all_feasible_adducts).
         - Computes per-residue site scores from the aligned shifted-ion evidence,
           lowercases every residue that ties for the maximum score in the resulting
           @c best_localization string, and stores the last-tied position in
           @c best_localization_position (the cpp's @c "Note: check if there are
           situations where multiple have the same score" caveat applies — ties resolve
           to the last index, not the first).
         - Writes @c localization_scores (comma-separated @c "0..100" percentages),
           @c best_localization, @c best_localization_score, and
           @c best_localization_position back into the hit.

    All output goes into @p annotated_hits — no exceptions are thrown by this function;
    malformed alignments emit @c OPENMS_LOG_WARN. Hits whose first precursor adduct
    cannot be resolved fall back to an unmodified-sequence localisation
    (@c best_localization @c == @c unmodified_sequence, @c best_localization_score @c = @c 0).

    The trailing underscore on the method name is a historical artefact from the
    OpenNuXL feature branch; the entry point is still public.

    @param[in]     exp                              MS2 experiment whose spectra are aligned against the theoretical spectra of the candidate hits.
    @param[in,out] annotated_hits                   Per-scan vectors of candidate hits to annotate and localise in place.
    @param[in]     mm                               Precursor-adduct mass table (see @ref NuXLModificationsGenerator::initModificationMassesNA); indexed by @c NuXLAnnotatedHit::NA_mod_index.
    @param[in]     fixed_modifications              Cached fixed-modification table consumed by @ref ModifiedPeptideGenerator.
    @param[in]     variable_modifications           Cached variable-modification table consumed by @ref ModifiedPeptideGenerator.
    @param[in]     max_variable_mods_per_peptide    Cap on the number of variable modifications placed per peptide.
    @param[in]     fragment_mass_tolerance          Half-width of the fragment-alignment tolerance window.
    @param[in]     fragment_mass_tolerance_unit_ppm Interpret @p fragment_mass_tolerance as ppm (true) or Da (false).
    @param[in]     all_feasible_adducts             Per-precursor-adduct fragment-adduct definitions used to generate marker and shifted immonium ions.
  */
  static void annotateAndLocate_(
    const PeakMap& exp,
    std::vector<std::vector<NuXLAnnotatedHit>>& annotated_hits,
    const NuXLModificationMassesResult& mm,
    const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications,
    const ModifiedPeptideGenerator::MapToResidueType& variable_modifications,
    Size max_variable_mods_per_peptide,
    double fragment_mass_tolerance,
    bool fragment_mass_tolerance_unit_ppm,
    const NuXLParameterParsing::PrecursorsToMS2Adducts & all_feasible_adducts);

};
} // namespace OpenMS


