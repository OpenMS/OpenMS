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

// Forward declarations for Arrow types (only when WITH_PARQUET is enabled)
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
    
    The parquet output contains columns for:
    - sequence: peptide sequence
    - spectrum_reference: spectrum identifier
    - charge: precursor charge
    - retention_time: RT in seconds  
    - mass_to_charge: precursor m/z
    - score: identification score
    - rank: hit rank
    - protein_accessions: protein accessions (comma-separated)
    - modifications: peptide modifications
    - is_decoy: decoy flag
    - search_engine: search engine name
    - search_engine_score_name: score type name

    @note This class requires OpenMS to be compiled with parquet support (WITH_PARQUET=ON).
          When WITH_PARQUET is not enabled, store() will throw Exception::NotImplemented.

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
      @throws Exception::NotImplemented if WITH_PARQUET is not defined
    */
    void store(const String& filename,
               const std::vector<ProteinIdentification>& protein_identifications,
               const PeptideIdentificationList& peptide_identifications);

  protected:
    
#ifdef WITH_PARQUET
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
#endif // WITH_PARQUET
  };

} // namespace OpenMS
