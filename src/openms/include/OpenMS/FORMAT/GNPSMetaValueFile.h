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
  class OPENMS_DLLAPI GNPSMetaValueFile
  {
    public:
      /// Generate a meta value table (tsv file) for GNPS FBMN with information on the input mzML files extracted from ConsensusMap.
      static void store(const ConsensusMap& consensus_map, const String& output_file);
  };
}
