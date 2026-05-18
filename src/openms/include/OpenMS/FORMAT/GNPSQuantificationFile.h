// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Axel Walter $
// $Authors: Axel Walter $

#pragma once

#include <OpenMS/KERNEL/ConsensusMap.h>

namespace OpenMS
{
  /**
    @brief Writer for the feature-quantification (.TXT) table required by GNPS Feature-Based Molecular Networking.

    For each consensus feature in a @ref ConsensusMap, this writer emits one row giving
    the consensus-level @c rt_cf, @c mz_cf, @c intensity_cf, @c charge_cf, @c width_cf,
    @c quality_cf (in that order), optionally followed by IIMN columns if
    @c Constants::UserParam::IIMN_ROW_ID is set on the first feature, followed by per-map
    sub-feature columns. The file starts with @c \#MAP and @c \#CONSENSUS header lines and
    one @c MAP row per @c ConsensusMap::getColumnHeaders entry (id, filename, label, size).
    The format matches the schema at
    https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/#feature-quantification-table.

    The file is **required** for FBMN (in contrast to the optional meta-value table).
    Produced by the @ref TOPP_GNPSExport tool.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI GNPSQuantificationFile
  {
    public:
      /**
        @brief Write the FBMN feature-quantification table.

        Output structure (tab-separated, written via @ref SVOutStream):
          - @c "#MAP\\tid\\tfilename\\tlabel\\tsize" header.
          - @c "#CONSENSUS\\trt_cf\\tmz_cf\\tintensity_cf\\tcharge_cf\\twidth_cf\\tquality_cf"
            header, followed by IIMN columns (@c IIMN_ROW_ID, @c IIMN_BEST_ION,
            @c IIMN_ADDUCT_PARTNERS, @c IIMN_ANNOTATION_NETWORK_NUMBER) iff the first
            consensus feature carries @c IIMN_ROW_ID, then per-map sub-feature columns
            @c "rt_&lt;i&gt;\\tmz_&lt;i&gt;\\tintensity_&lt;i&gt;\\tcharge_&lt;i&gt;\\twidth_&lt;i&gt;".
          - One @c MAP row per @p consensus_map.getColumnHeaders entry.
          - One @c CONSENSUS row per consensus feature; missing per-map sub-features
            are filled with empty strings.

        @param[in] consensus_map Consensus map whose column headers and features drive the rows.
        @param[in] output_file   Destination path for the TXT file (overwritten if it exists).
      */
      static void store(const ConsensusMap& consensus_map, const String& output_file);
  };
}
