// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <optional>
#include <string>
#include <unordered_map>

namespace OpenMS
{
  /**
    @brief Establish and validate canonical precursor/transition identifiers for OpenSWATH libraries.

    OpenSWATH algorithms use string identifiers in @ref OpenSwath::LightTargetedExperiment,
    while persistent PQP/OSWPQ schemas use integer foreign keys. This helper defines the
    application-level invariant that the LightTargetedExperiment carries those canonical
    integer identifiers encoded as decimal strings:

      - @c LightCompound.id is the canonical precursor ID.
      - @c LightTransition.peptide_ref references that canonical precursor ID.
      - @c LightTransition.transition_name is the canonical transition ID.

    Arbitrary source identifiers (for example from TSV or TraML) are normalized once at the
    OpenSWATH library-loading boundary. Database-backed inputs (PQP/OSWPQ) already carry
    persistent numeric IDs and are validated instead of renumbered.

    @ingroup TargetedQuantitation
  */
  class OPENMS_DLLAPI OpenSwathLibraryIDNormalizer
  {
  public:
    /**
      @brief Source-to-canonical precursor provenance captured during library loading.

      The operational LightTargetedExperiment contains only canonical numeric IDs.
      This mapping retains the original precursor identifiers long enough for callers
      that still need source-ID semantics (for example configurable target/decoy
      pairing) without reusing source identifiers as operational foreign keys.
    */
    struct SourceIDMapping
    {
      /// Original/source precursor ID -> canonical operational precursor ID.
      std::unordered_map<std::string, std::string> precursor_source_to_canonical;

      /// Canonical operational precursor ID -> original/source precursor ID.
      std::unordered_map<std::string, std::string> precursor_canonical_to_source;

      /// Canonical operational transition ID -> original/source transition ID.
      /// Source transition IDs are not required to be unique, so only the canonical-to-source
      /// direction is retained.
      std::unordered_map<std::string, std::string> transition_canonical_to_source;
    };

    /**
      @brief Replace arbitrary source precursor/transition identifiers with deterministic canonical IDs.

      Precursor IDs are assigned by lexicographically sorting all unique source compound IDs and
      numbering them from zero, matching the persistent PQP convention. Compounds that are not
      referenced by any transition are then omitted from the operational LightTargetedExperiment
      without compressing the remaining IDs, so sparse precursor IDs are expected. This keeps direct
      source loading consistent with source-to-PQP round-trips. Transition IDs are assigned from their
      order in @p exp.transitions, also starting at zero. Existing transition @c peptide_ref values
      are rewritten to the new precursor IDs.

      Normalization rebuilds the LightTargetedExperiment rather than changing compound IDs in place,
      ensuring that its internal compound-reference lookup cache cannot retain source-ID keys after
      canonicalization.

      This function is intended for source-oriented formats such as TSV and TraML. It must not
      be used to renumber libraries that already contain persistent canonical IDs.

      @param[in,out] exp Source-oriented experiment to canonicalize in place.
      @return Source-ID provenance for the canonicalized experiment. The precursor maps cover
              every source compound, including compounds omitted from the operational @p exp.
      @throws Exception::InvalidValue If a compound ID is empty/duplicated or a transition
              references an unknown/empty compound ID.
    */
    static SourceIDMapping normalizeSourceIDs(OpenSwath::LightTargetedExperiment& exp);

    /**
      @brief Materialize decoy flags from a configurable source precursor prefix.

      Canonicalization replaces source precursor references and transition names with numeric
      IDs. This helper uses the provenance returned by normalizeSourceIDs() (or reconstructed
      by a persistent-format reader) to preserve prefix-based decoy semantics before filtering.
      A configured prefix on either the source precursor ID or source transition ID materializes
      the explicit transition decoy flag. Existing explicit decoy flags are retained.

      @param[in,out] exp Canonical experiment whose transition decoy flags may be updated.
      @param[in] source_ids Source precursor/transition provenance for the canonical experiment.
      @param[in] decoy_prefix Configured source precursor decoy prefix. Empty disables
                 additional prefix materialization.
    */
    static void materializeDecoyPrefix(OpenSwath::LightTargetedExperiment& exp,
                                       const SourceIDMapping& source_ids,
                                       const std::string& decoy_prefix);

    /**
      @brief Resolve the canonical target precursor paired with a source-prefix decoy.

      @p decoy_ref is normally a canonical precursor ID. The source precursor ID is
      recovered from @p source_ids, @p decoy_prefix is removed there, and the matching
      target source ID is mapped back to its canonical operational ID. If no provenance
      mapping is available, source-ID input is supported as a compatibility path.

      @return The matching target precursor ID in the same operational domain as the
              loaded experiment, or std::nullopt if the decoy cannot be paired.
    */
    static std::optional<std::string> canonicalTargetForDecoyPrecursor(
      const std::string& decoy_ref,
      const SourceIDMapping& source_ids,
      const std::string& decoy_prefix);

    /**
      @brief Report whether a LightTargetedExperiment satisfies the canonical-ID invariant.

      @param[in] exp Experiment to inspect.
      @return True if @c validateCanonicalIDs() would accept @p exp.
    */
    static bool hasCanonicalIDs(const OpenSwath::LightTargetedExperiment& exp);

    /**
      @brief Report whether all operational ID fields use canonical decimal syntax.

      This checks only the textual ID representation. It deliberately does not check uniqueness
      or whether transition precursor references resolve. Writer compatibility overloads use it
      to distinguish source-style identifiers from malformed canonical-looking input.

      @param[in] exp Experiment to inspect.
      @return True if every compound ID, transition ID and transition precursor reference is a
              non-negative Int64 in canonical decimal form.
    */
    static bool hasCanonicalIDFormat(const OpenSwath::LightTargetedExperiment& exp) noexcept;

    /**
      @brief Validate that a LightTargetedExperiment already uses canonical numeric IDs.

      IDs must be unique, non-negative Int64 values written in canonical decimal form
      (e.g. @c "7", not @c "007", @c "+7", or @c "-0"). IDs may be sparse; filtering must not force
      renumbering. Every transition @c peptide_ref must exactly match an existing compound ID.

      @param[in] exp Experiment to validate.
      @throws Exception::InvalidValue If any canonical-ID invariant is violated.
    */
    static void validateCanonicalIDs(const OpenSwath::LightTargetedExperiment& exp);
  };
} // namespace OpenMS
