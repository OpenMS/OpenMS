// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: OpenMS Team $
// $Authors: OpenMS Team $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <string>
#include <map>
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
    /// Copy constructor
    BedRModFile(const BedRModFile&) = default;

    /// Destructor
    ~BedRModFile() = default;

    /**
      @brief Store RNA modification sites as bedRMod v2.

      Canonical residues contribute to the frequency denominator but do not
      produce modification rows. Frequencies count distinct observations at each
      site, including alternative hypotheses and overlapping parent mappings.
      MODOMICS short names identify modifications without a positive ChEBI mapping.

      @param[in] out_file Output path for the bedRMod file.
      @param[in] id_data Identification data containing RNA observation matches.
      @param[in] chebi_mapping_file Optional CSV mapping file with columns
                 "mod" (or "name") and "chebi_id" (or "chebi id").
      @param[in] reference_metadata Optional organism, assembly, annotation_source
                 and annotation_version values. Unspecified values are written as NA.
      @throws Exception::InvalidParameter if reference metadata contains a line break.
      @throws Exception::FileNotWritable if the output file cannot be opened.
    */
    void store(const std::string& out_file,
               const IdentificationData& id_data,
               const std::string& chebi_mapping_file = "",
               const std::map<std::string, std::string>& reference_metadata = {}) const;
  };
}
