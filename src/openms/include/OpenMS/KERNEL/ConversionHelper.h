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
    @brief Static helpers that convert between the three OpenMS container types
           (@c PeakMap, @c FeatureMap, @c ConsensusMap).

    Every overload clears the destination container before writing, sets a
    unique id on the output, and calls @c updateRanges() on the output so
    range queries return the expected values immediately after the conversion.

    The conversions are not symmetric: which fields are preserved (unique ids,
    protein / unassigned peptide identifications) depends on the source and
    target container -- see each overload for the exact behaviour.

    @ingroup Kernel
  */
  class OPENMS_DLLAPI MapConversion
  {
public:

    /**
      @brief Copy the @p n most intense peaks of a @c PeakMap into a
             @c ConsensusMap.

      The output is cleared and assigned a fresh unique id (@c PeakMap has no
      container-level unique id, so one is generated). The input's ranges are
      updated (the underlying 2D peak data is then materialized via
      @c get2DData) and the result is partial-sorted by intensity in descending
      order so the top @p n peaks are kept. Each surviving peak is appended as
      a @c ConsensusFeature carrying @p input_map_index, the peak, and its
      rank in the sorted list. The column header @c size for
      @p input_map_index is set to @p n, and the output's ranges are updated.

      @param[in]     input_map_index Index assigned to the input map in the
                                     resulting @c ConsensusMap column headers.
      @param[in,out] input_map       Source peaks; @c updateRanges() is called
                                     on it as a side effect.
      @param[out]    output_map      Resulting @c ConsensusMap; previous
                                     contents are cleared.
      @param[in]     n               Maximum number of peaks to copy. The
                                     default (@c Size(-1)) is clamped to
                                     @c input_map.getSize(), i.e. all peaks
                                     are copied.
    */
    static void convert(UInt64 const input_map_index,
                        PeakMap& input_map,
                        ConsensusMap& output_map,
                        Size n = -1);

    /**
      @brief Convert a @c ConsensusMap to a @c FeatureMap.

      The output is cleared and resized to @c input_map.size(). The input's
      @c DocumentIdentifier is copied to the output, then either the
      container-level unique id is copied as well (when @p keep_uids is
      @c true) or a new one is generated (when @c false). Protein
      identifications and unassigned peptide identifications are copied
      verbatim from input to output. For each element, the @c BaseFeature
      portion of the @c ConsensusFeature is copied into the corresponding
      @c Feature; when @p keep_uids is @c false the per-feature unique id is
      regenerated. Finally @c updateRanges() is called on the output.

      @param[in]  input_map  Source @c ConsensusMap.
      @param[in]  keep_uids  If @c true, preserve the container-level and
                             per-element unique ids from @p input_map;
                             otherwise mint new ones.
      @param[out] output_map Resulting @c FeatureMap; previous contents are
                             cleared.
    */
    static void convert(ConsensusMap const& input_map,
                        const bool keep_uids,
                        FeatureMap& output_map);

    /**
      @brief Convert a @c FeatureMap to a @c ConsensusMap.

      The output is cleared, its container-level unique id is overwritten with
      the @c FeatureMap's unique id (an intentional design choice -- callers
      that need a fresh id must overwrite it afterwards), and the first @p n
      features (in input order; no sorting is performed) are appended as
      @c ConsensusFeature entries tagged with @p input_map_index. Protein and
      unassigned peptide identifications are copied verbatim from input to
      output. The column header @c size for @p input_map_index is set to
      @c input_map.size() (i.e. the full input size, even when @p n truncates
      the copy), and @c updateRanges() is called on the output.

      @note Because elements are taken in input order, @p n is mainly useful
            after pre-sorting @p input_map (e.g. by intensity); it exists for
            symmetry with the @c PeakMap overload above.

      @param[in]  input_map_index Index assigned to the input map in the
                                  resulting @c ConsensusMap column headers.
      @param[in]  input_map       Source @c FeatureMap.
      @param[out] output_map      Resulting @c ConsensusMap; previous contents
                                  are cleared.
      @param[in]  n               Maximum number of features to copy. The
                                  default (@c Size(-1)) is clamped to
                                  @c input_map.size(), i.e. all features are
                                  copied.
    */
    static void convert(UInt64 const input_map_index,
                        FeatureMap const& input_map,
                        ConsensusMap& output_map,
                        Size n = -1);
  };
} // namespace OpenMS

