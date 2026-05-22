// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenMS
{

  /**
    @brief Static helpers that convert between the three OpenMS container
           types (@c PeakMap, @c FeatureMap, @c ConsensusMap).

    Each overload replaces the contents of the destination, gives it a
    container-level unique id, and leaves it in a state where the standard
    range queries reflect the new contents without further bookkeeping by
    the caller.

    The conversions are not symmetric: which fields are preserved
    (container / element unique ids, protein and unassigned peptide
    identifications) depends on the source and target container, so see
    each overload for the exact contract.

    @ingroup Kernel
  */
  class OPENMS_DLLAPI MapConversion
  {
public:

    /**
      @brief Copy the most intense peaks of a @c PeakMap into a
             @c ConsensusMap.

      The output's previous contents are dropped and it is given a fresh
      container unique id (@c PeakMap has no container-level unique id, so
      one is generated). The @p n peaks with the highest intensity are
      written to the output as @c ConsensusFeature entries tagged with
      @p input_map_index; their order in the output is by descending
      intensity. The column header @c size for @p input_map_index reflects
      the number of peaks written.

      @param[in]     input_map_index Index assigned to the input map in the
                                     resulting @c ConsensusMap column headers.
      @param[in,out] input_map       Source peaks; its range queries are
                                     left consistent with its contents as a
                                     side effect.
      @param[out]    output_map      Resulting @c ConsensusMap; previous
                                     contents are replaced.
      @param[in]     n               Maximum number of peaks to copy. The
                                     default (@c Size(-1)) keeps all peaks.
    */
    static void convert(UInt64 const input_map_index,
                        PeakMap& input_map,
                        ConsensusMap& output_map,
                        Size n = -1);

    /**
      @brief Convert a @c ConsensusMap to a @c FeatureMap.

      Every element of @p input_map is converted to a @c Feature in
      @p output_map (the @c BaseFeature portion is copied; positional and
      meta data are preserved). The document identifier and the protein /
      unassigned peptide identifications are preserved on the output.

      @param[in]  input_map  Source @c ConsensusMap.
      @param[in]  keep_uids  If @c true, the container unique id and every
                             element's unique id are preserved from
                             @p input_map; otherwise they are replaced with
                             fresh ones.
      @param[out] output_map Resulting @c FeatureMap; previous contents are
                             replaced.
    */
    static void convert(ConsensusMap const& input_map,
                        const bool keep_uids,
                        FeatureMap& output_map);

    /**
      @brief Convert a @c FeatureMap to a @c ConsensusMap.

      The first @p n features of @p input_map (in input order; no sorting
      is performed) are written to @p output_map as @c ConsensusFeature
      entries tagged with @p input_map_index. The output's container unique
      id is taken from @p input_map (an intentional design choice -- callers
      that need a fresh id must overwrite it afterwards). The protein and
      unassigned peptide identifications are preserved on the output.

      @note Because features are taken in input order, @p n is mainly useful
            after pre-sorting @p input_map (e.g. by intensity); it exists for
            symmetry with the @c PeakMap overload above.

      @note The column header @c size for @p input_map_index is set to
            @c input_map.size() -- the full input size -- even when @p n
            truncates the actual copy. Inspect the output container's size
            for the number of features that were actually written.

      @param[in]  input_map_index Index assigned to the input map in the
                                  resulting @c ConsensusMap column headers.
      @param[in]  input_map       Source @c FeatureMap.
      @param[out] output_map      Resulting @c ConsensusMap; previous
                                  contents are replaced.
      @param[in]  n               Maximum number of features to copy. The
                                  default (@c Size(-1)) keeps all features.
    */
    static void convert(UInt64 const input_map_index,
                        FeatureMap const& input_map,
                        ConsensusMap& output_map,
                        Size n = -1);
  };
} // namespace OpenMS

