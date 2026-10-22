// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Sam Wein $
// $Authors: Sam Wein $
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
      canonical-base site, including alternative hypotheses and overlapping parent
      mappings. One row aggregates each target parent, position and modification;
      its score is the minimum supporting PSM-level q-value (lower is better).
      Decoy parents and matches without assigned q-values are omitted. Callers
      must apply their chosen FDR threshold and ranking before export.
      Fragment columns list distinct supporting intervals; the mapping-count
      column is the maximum number of target mappings among supporting oligos.
      MODOMICS short names identify modifications without a positive ChEBI mapping.

      @param[in] out_file Output path for the bedRMod file.
      @param[in] id_data Identification data containing RNA observation matches.
      @param[in] chebi_mapping_file Optional CSV mapping file with columns
                 "mod" (or "name") and "chebi_id" (or "chebi id") and/or
                 "modomics_short_name". The latter explicitly maps custom codes
                 to registered MODOMICS names with the same formula and origin.
      @param[in] reference_metadata Optional organism, assembly, annotation_source
                 annotation_version, sequencing_platform, basecalling,
                 bioinformatics_workflow, experiment and external_source values.
                 Unspecified values are written as NA.
      @throws Exception::InvalidParameter if reference metadata contains a line
              break, or a modification cannot be represented by a valid MODOMICS
              dictionary name. Output is not opened before these checks succeed.
      @throws Exception::FileNotWritable if the output file cannot be opened.
    */
    void store(const std::string& out_file,
               const IdentificationData& id_data,
               const std::string& chebi_mapping_file = "",
               const std::map<std::string, std::string>& reference_metadata = {}) const;

    /**
      @brief Write the observation-level bedRMod layout used by the RNOME stage-8 pipeline.

      Includes canonical and decoy rows and preserves the historical column names,
      coverage, frequency, mapping, and fragment conventions. This is a compatibility
      format for existing RNOME downstream tooling; @ref store writes site-level calls.
    */
    void storeLegacy(const std::string& out_file,
                     const IdentificationData& id_data,
                     const std::string& chebi_mapping_file = "",
                     const std::map<std::string, std::string>& reference_metadata = {}) const;
  };
}
