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

      @throws Exception::InvalidValue If a compound ID is empty/duplicated or a transition
              references an unknown/empty compound ID.
    */
    static void normalizeSourceIDs(OpenSwath::LightTargetedExperiment& exp);

    /**
      @brief Validate that a LightTargetedExperiment already uses canonical numeric IDs.

      IDs must be unique, non-negative Int64 values written in canonical decimal form
      (e.g. @c "7", not @c "007" or @c "+7"). IDs may be sparse; filtering must not force
      renumbering. Every transition @c peptide_ref must exactly match an existing compound ID.

      @throws Exception::InvalidValue If any canonical-ID invariant is violated.
    */
    static void validateCanonicalIDs(const OpenSwath::LightTargetedExperiment& exp);
  };
} // namespace OpenMS
