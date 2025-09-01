// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

#include <vector>

// Forward declarations for Arrow types
namespace arrow {
  class Schema;
  class Table;
}

namespace OpenMS
{
  /**
    @brief File adapter for writing PSM (Peptide Spectrum Match) data to parquet files

    This class converts OpenMS ProteinIdentification and PeptideIdentification 
    objects to parquet format following the quantms.io PSM specification.
    
    The parquet output contains columns following the quantms.io PSM specification:
    - sequence: unmodified peptide sequence
    - peptidoform: peptide sequence with modifications
    - modifications: peptide modifications (null for now)
    - precursor_charge: precursor charge
    - posterior_error_probability: PEP score from metavalues (nullable)
    - is_decoy: decoy flag (0=target, 1=decoy) based on target_decoy metavalue
    - calculated_mz: theoretical m/z from sequence
    - observed_mz: experimental precursor m/z
    - additional_scores: additional scores (null for now)
    - mp_accessions: protein accessions (null for now)
    - predicted_rt: predicted retention time (null for now)
    - reference_file_name: reference file name
    - cv_params: CV parameters (null for now)
    - scan: scan identifier
    - rt: retention time in seconds (nullable)
    - ion_mobility: ion mobility value (nullable, null for now)
    - num_peaks: number of peaks (nullable, null for now)
    - mz_array: m/z values array (null for now)
    - intensity_array: intensity values array (null for now)

    Only the first peptide hit per peptide identification is processed (no rank field).
    PEP scores are automatically detected from metavalues using known PEP score names.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI QuantmsIO :
    public ProgressLogger
  {
  public:
    /// Default constructor
    QuantmsIO();

    /// Destructor
    ~QuantmsIO();

    /**
      @brief Store peptide and protein identifications in parquet format

      @param filename Output filename (should end with .parquet)
      @param protein_identifications Vector of protein identifications
      @param peptide_identifications Vector of peptide identifications

      @throws Exception::UnableToCreateFile if file cannot be created
    */
    void store(const String& filename,
               const std::vector<ProteinIdentification>& protein_identifications,
               const PeptideIdentificationList& peptide_identifications);

  protected:
    
    /**
      @brief Create Arrow schema for PSM data following quantms.io specification
      
      @return Shared pointer to Arrow schema
    */
    std::shared_ptr<arrow::Schema> createPSMSchema_();

    /**
      @brief Convert peptide identifications to Arrow table
      
      @param protein_identifications Vector of protein identifications (for metadata)
      @param peptide_identifications Vector of peptide identifications
      @return Shared pointer to Arrow table
    */
    std::shared_ptr<arrow::Table> convertToArrowTable_(
      const std::vector<ProteinIdentification>& protein_identifications,
      const PeptideIdentificationList& peptide_identifications);

    /**
      @brief Write Arrow table to parquet file
      
      @param table Arrow table to write
      @param filename Output filename
    */
    void writeParquetFile_(const std::shared_ptr<arrow::Table>& table, 
                          const String& filename);
  };

} // namespace OpenMS
