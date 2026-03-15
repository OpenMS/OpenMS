// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: OpenMS Team $
// $Authors: OpenMS Team $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>

namespace OpenMS
{
  /**
    @brief File adapter for bedRMod v2 files.

    Exports RNA modification sites from IdentificationData into the BED-style
    bedRMod v2 format.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI BedRModFile
  {
  public:
    /// Default constructor
    BedRModFile() = default;

    /// Destructor
    ~BedRModFile() = default;

    /**
      @brief Store RNA modification sites as bedRMod v2.

      @param[in] out_file Output path for the bedRMod file.
      @param[in] id_data Identification data containing RNA observation matches.
      @param[in] chebi_mapping_file Optional CSV mapping file with columns
                 "mod" (or "name") and "chebi_id" (or "chebi id").
    */
    void store(const String& out_file,
               const IdentificationData& id_data,
               const String& chebi_mapping_file = "") const;
  };
}
