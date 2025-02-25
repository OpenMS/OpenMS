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

  class OPENMS_DLLAPI MapConversion
  {
public:

    /**
      @brief Similar to @p convert for FeatureMaps.

      Only the @p n most intense elements are copied.

      Currently PeakMap does not have a unique id but ConsensusMap has
      one, so we assign a new one here.

      @param input_map_index The index of the input map.
      @param input_map The input map to be converted.
      @param output_map The resulting ConsensusMap.
      @param n The maximum number of elements to be copied.
    */
    static void convert(UInt64 const input_map_index,
                        PeakMap& input_map,
                        ConsensusMap& output_map,
                        Size n = -1);

    /**
      @brief Convert a ConsensusMap to a FeatureMap (of any feature type).

      The previous content of output_map is cleared. UID's of the elements and
      the container is copied if the @p keep_uids flag is set.

      @param input_map The container to be converted.
      @param keep_uids Shall the UID's of the elements and the container be kept or created anew
      @param output_map The resulting ConsensusMap.
    */
    static void convert(ConsensusMap const& input_map,
                        const bool keep_uids,
                        FeatureMap& output_map);

    /**
      @brief Convert a FeatureMap (of any feature type) to a ConsensusMap.

      Each ConsensusFeature contains a map index, so this has to be given as
      well. The previous content of @p output_map is cleared. An arguable
      design decision is that the unique id of the FeatureMap is copied (!) to
      the ConsensusMap, because that is the way it is meant to be used in the
      algorithms.

      Only the first (!) @p n elements are copied. (This parameter exists
      mainly for compatibility with @p convert for MSExperiments. To use it in
      a meaningful way, apply one of the sorting methods to @p input_map
      beforehand.)

      @param input_map_index The index of the input map.
      @param input_map The container to be converted.
      @param output_map The resulting ConsensusMap.
      @param n The maximum number of elements to be copied.
    */
    static void convert(UInt64 const input_map_index,
                        FeatureMap const& input_map,
                        ConsensusMap& output_map,
                        Size n = -1);
  };
} // namespace OpenMS

