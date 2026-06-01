// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Axel Walter $
// $Authors: Axel Walter $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/ConsensusMap.h>

namespace OpenMS
{
  /**
    @brief Helpers for the Ion Identity Molecular Networking (IIMN) workflow on GNPS.

    IIMN extends GNPS Feature-Based Molecular Networking by grouping consensus features
    that represent the same compound seen as different adducts (e.g. [M+H]+ and [M+Na]+).
    The grouping information is usually produced upstream by an adduct decharger (e.g.
    @ref TOPP_MetaboliteAdductDecharger), which records a set of adduct-group IDs in the
    @c Constants::UserParam::IIMN_LINKED_GROUPS meta value on each consensus feature.

    Two stages live in this class:

      -# annotateConsensusMap() consumes @c IIMN_LINKED_GROUPS, builds a bipartite graph
         (consensus feature ↔ adduct group), runs Boost connected-components on it, and
         writes back three IIMN meta values: a 1-based row id, the network number, and
         the list of partner-feature row ids. The @c IIMN_LINKED_GROUPS meta value is
         consumed and removed.
      -# writeSupplementaryPairTable() turns those annotations into the CSV file required
         by GNPS IIMN as the "supplementary pairs" upload.

    Output of stage 1 is also a precondition for the meta-value / quantification tables
    written by @ref GNPSMetaValueFile and @ref GNPSQuantificationFile, and is invoked
    by the @ref TOPP_GNPSExport tool.

    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI IonIdentityMolecularNetworking
  {
    public:
      /**
        @brief Populate the IIMN meta values on @p consensus_map in-place.

        Reads @c Constants::UserParam::IIMN_LINKED_GROUPS (a list of adduct-group ids,
        typically from @ref TOPP_MetaboliteAdductDecharger), constructs a bipartite graph
        whose nodes are either consensus features or adduct groups, computes its
        connected components, and writes three meta values per feature:

          - @c IIMN_ROW_ID — 1-based index of the feature in @p consensus_map (assigned
            to every feature, including those with no group annotation).
          - @c IIMN_ANNOTATION_NETWORK_NUMBER — 1-based connected-component id;
            features that are direct or transitive adduct partners share this id.
            Set only on features that carried @c IIMN_LINKED_GROUPS.
          - @c IIMN_ADDUCT_PARTNERS — semicolon-separated list of partner @c IIMN_ROW_ID
            values; set only on features that share at least one adduct group.

        The @c IIMN_LINKED_GROUPS meta value is consumed and removed from every feature
        when the routine returns.

        @param[in,out] consensus_map Map to annotate. Features that lack @c IIMN_LINKED_GROUPS
                                     receive only @c IIMN_ROW_ID.
      */
      static void annotateConsensusMap(ConsensusMap& consensus_map);

      /**
        @brief Write the IIMN supplementary pairs CSV consumed by GNPS.

        Emits one row per pair of features that share at least one adduct group, with
        columns:

          - @c ID1, @c ID2 — @c IIMN_ROW_ID of the two features (from @ref annotateConsensusMap).
          - @c EdgeType   — fixed string "MS1 annotation".
          - @c Score      — number of direct partners across both sides of the pair
                            (sum of each side's partner counts, minus 2 to remove the
                            self-counts).
          - @c Annotation — best-ion adduct string for each side (defaulting to "default"
                            if no @c IIMN_BEST_ION was set) plus the absolute delta m/z
                            between the two precursors.

        Reverse edges are deduplicated as the file is written (only one row per unordered
        pair). The schema matches https://ccms-ucsd.github.io/GNPSDocumentation/fbmn-iin/#supplementary-pairs.

        @note Returns silently if the first feature has no @c IIMN_ROW_ID (i.e.
              @ref annotateConsensusMap has not been run on @p consensus_map).

        @param[in] consensus_map Annotated consensus map (call @ref annotateConsensusMap first).
        @param[in] output_file   Destination CSV path (overwritten if it exists).
      */
      static void writeSupplementaryPairTable(const ConsensusMap& consensus_map, const String& output_file);
  };
} // closing namespace OpenMS