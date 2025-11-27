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
  class OPENMS_DLLAPI GNPSQuantificationFile
  {
    public:
      /// Write feature quantification table (txt file) from a consensusXML file. Required for GNPS FBMN.
      /// The table contains map information on the featureXML files from which the consensusXML file was generated as well as
      /// a row for every consensus feature with information on rt, mz, intensity, width and quality. The same information is
      /// added for each original feature in the consensus feature.
      static void store(const ConsensusMap& consensus_map, const String& output_file);
  };
}
