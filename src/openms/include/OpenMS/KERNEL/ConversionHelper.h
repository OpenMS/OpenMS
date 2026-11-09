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
#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenMS
{

  /**
    @brief Static helpers that convert between @c FeatureMap and @c ConsensusMap.

    There is no conversion from a @c PeakMap: the handle of a consensus
    feature built from a peak could only refer to the peak's index, not to
    an element with a unique id. MapAlignmentAlgorithmPoseClustering selects
    the peaks it aligns itself.

    Each overload replaces the contents of the destination and leaves it in
    a state where the standard range queries reflect the new contents
    without further bookkeeping by the caller. The container unique id of
    the destination is copied from the source or freshly drawn, depending
    on the overload.

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
      @brief Convert a @c ConsensusMap to a @c FeatureMap.

      Every element of @p input_map is converted to a @c Feature in
      @p output_map (the @c BaseFeature portion is copied; positional and
      meta data are preserved). The document identifier, the protein /
      unassigned peptide identifications and the identification data with
      the elements' links to it are preserved on the output.

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
      unassigned peptide identifications are preserved on the output. The
      peptide identifications of the copied features are tagged with
      @p input_map_index (meta value @c map_index); those of features that
      are not copied and of subordinates are dropped. The identification
      data follows the same rule: identifications linked by copied features
      get @c map_index, identifications linked only by other features or
      subordinates are removed, unlinked ones are kept.

      @note Because features are taken in input order, @p n is mainly useful
            after pre-sorting @p input_map (e.g. by intensity).

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

