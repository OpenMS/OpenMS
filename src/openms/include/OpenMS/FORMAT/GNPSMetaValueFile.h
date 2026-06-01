// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Axel Walter $
// $Authors: Axel Walter $

#include <OpenMS/KERNEL/ConsensusMap.h>

#pragma once

namespace OpenMS
{
  /**
    @brief Writer for the meta-value (.TSV) table consumed by GNPS Feature-Based Molecular Networking.

    GNPS FBMN can pick up per-sample metadata so the web platform can group / colour
    network nodes by sample attribute. The writer here produces the minimal seed of
    such a table: one row per mzML file used to build the input @ref ConsensusMap,
    populated with that file's basename and a sequential @c MAP&lt;i&gt; identifier. Users
    typically extend the file with additional attribute columns before uploading to
    GNPS (see https://ccms-ucsd.github.io/GNPSDocumentation/metadata/).

    The file is optional for FBMN — supplying it lets users group/colour by attribute
    on the GNPS side. Produced by the @ref TOPP_GNPSExport tool.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI GNPSMetaValueFile
  {
    public:
      /**
        @brief Write the FBMN meta-value TSV seed.

        Pulls the primary mzML run paths from @p consensus_map via
        @c ConsensusMap::getPrimaryMSRunPath and emits one row per path. The file
        starts with a header row @c "\\tfilename\\tATTRIBUTE_MAPID" (the first
        column header is empty); each data row is
        @c "&lt;index&gt;\\t&lt;basename(path)&gt;\\tMAP&lt;index&gt;" with @c index
        running from 0.

        @param[in] consensus_map Consensus map whose primary MS-run paths provide the file list.
        @param[in] output_file   Destination path for the TSV file (overwritten if it exists).
      */
      static void store(const ConsensusMap& consensus_map, const String& output_file);
  };
}
