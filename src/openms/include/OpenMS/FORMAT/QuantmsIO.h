// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

#include <vector>
#include <set>

namespace OpenMS
{
  /**
    @brief File adapter for writing PSM (Peptide Spectrum Match) data to parquet or TSV files

    This class converts OpenMS ProteinIdentification and PeptideIdentification 
    objects to parquet or TSV format following the quantms.io PSM specification.
    
    The output contains columns following the quantms.io PSM specification:
    - sequence: unmodified peptide sequence
    - peptidoform: peptide sequence with modifications
    - modifications: peptide modifications (null for now)
    - precursor_charge: precursor charge
    - posterior_error_probability: PEP score from metavalues (nullable)
    - is_decoy: decoy flag (0=target, 1=decoy) based on target_decoy metavalue
    - calculated_mz: theoretical m/z from sequence
    - observed_mz: experimental precursor m/z
    - additional_scores: additional scores (null for now)
    - protein_accessions: protein accessions (null for now)
    - predicted_rt: predicted retention time (null for now)
    - reference_file_name: reference file name
    - cv_params: CV parameters (null for now)
    - scan: scan identifier
    - rt: retention time in seconds (nullable)
    - ion_mobility: ion mobility value (nullable, null for now)
    - number_peaks: number of peaks (nullable, null for now)
    - mz_array: m/z values array (null for now)
    - intensity_array: intensity values array (null for now)
    - file_metadata: file-level metadata with quantmsio_version (1.0), creator (OpenMS), file_type (psm), creation_date (actual timestamp), uuid (generated), scan_format (scan), software_provider (OpenMS) - only preserved in parquet format

    Only the first peptide hit per peptide identification is processed by default (no rank field).
    When export_all_psms is enabled, all peptide hits are processed with a rank field.
    PEP scores are automatically detected from metavalues using known PEP score names.
    Optional meta value columns can be added for specific keys.
    
    Output format is determined by file extension:
    - .parquet: Apache Parquet format with embedded metadata
    - .tsv, .csv, .txt: Tab-separated values format (metadata not preserved)

    @ingroup FileIO
  */
  class OPENMS_DLLAPI QuantmsIO :
    public ProgressLogger
  {
  public:
    /// Default constructor
    QuantmsIO() = default;

    /// Destructor
    ~QuantmsIO();

    /**
      @brief Store peptide and protein identifications in parquet or TSV format

      @param filename Output filename (format determined by extension: .parquet, .tsv, .csv, .txt)
      @param protein_identifications Vector of protein identifications
      @param peptide_identifications Vector of peptide identifications

      @throws Exception::UnableToCreateFile if file cannot be created
      @throws Exception::InvalidParameter if unsupported file format is specified
    */
    void store(const String& filename,
               const std::vector<ProteinIdentification>& protein_identifications,
               const PeptideIdentificationList& peptide_identifications);

    /**
      @brief Store peptide and protein identifications in parquet or TSV format with all PSMs

      @param filename Output filename (format determined by extension: .parquet, .tsv, .csv, .txt)
      @param protein_identifications Vector of protein identifications
      @param peptide_identifications Vector of peptide identifications
      @param export_all_psms If true, export all PSMs per spectrum with rank column

      @throws Exception::UnableToCreateFile if file cannot be created
      @throws Exception::InvalidParameter if unsupported file format is specified
    */
    void store(const String& filename,
               const std::vector<ProteinIdentification>& protein_identifications,
               const PeptideIdentificationList& peptide_identifications,
               bool export_all_psms);

    /**
      @brief Store peptide and protein identifications in parquet or TSV format with enhanced options

      @param filename Output filename (format determined by extension: .parquet, .tsv, .csv, .txt)
      @param protein_identifications Vector of protein identifications
      @param peptide_identifications Vector of peptide identifications
      @param export_all_psms If true, export all PSMs per spectrum with rank column. If false, export only first PSM
      @param meta_value_keys Set of meta value keys to export as additional columns

      @throws Exception::UnableToCreateFile if file cannot be created
      @throws Exception::InvalidParameter if unsupported file format is specified
    */
    void store(const String& filename,
               const std::vector<ProteinIdentification>& protein_identifications,
               const PeptideIdentificationList& peptide_identifications,
               bool export_all_psms,
               const std::set<String>& meta_value_keys);

  };

} // namespace OpenMS
