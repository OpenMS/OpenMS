// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/QuantmsIO.h>

#ifdef WITH_PARQUET

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideEvidence.h>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>

#include <algorithm>

using namespace std;

namespace OpenMS
{
  QuantmsIO::QuantmsIO() = default;

  QuantmsIO::~QuantmsIO() = default;

  void QuantmsIO::store(const String& filename,
                       const std::vector<ProteinIdentification>& protein_identifications,
                       const PeptideIdentificationList& peptide_identifications)
  {
    // Convert data to Arrow table
    auto table = convertToArrowTable_(protein_identifications, peptide_identifications);
    
    // Write to parquet file
    writeParquetFile_(table, filename);
  }

  std::shared_ptr<arrow::Schema> QuantmsIO::createPSMSchema_()
  {
    std::vector<std::shared_ptr<arrow::Field>> fields = {
      arrow::field("sequence", arrow::utf8()),
      arrow::field("spectrum_reference", arrow::utf8()),
      arrow::field("charge", arrow::int32()),
      arrow::field("retention_time", arrow::float64()),
      arrow::field("mass_to_charge", arrow::float64()),
      arrow::field("score", arrow::float64()),
      arrow::field("rank", arrow::int32()),
      arrow::field("protein_accessions", arrow::utf8()),
      arrow::field("modifications", arrow::utf8()),
      arrow::field("is_decoy", arrow::boolean()),
      arrow::field("search_engine", arrow::utf8()),
      arrow::field("search_engine_score_name", arrow::utf8())
    };

    return arrow::schema(fields);
  }

  std::shared_ptr<arrow::Table> QuantmsIO::convertToArrowTable_(
    const std::vector<ProteinIdentification>& protein_identifications,
    const PeptideIdentificationList& peptide_identifications)
  {
    // Create builders for each column
    arrow::StringBuilder sequence_builder;
    arrow::StringBuilder spectrum_reference_builder;
    arrow::Int32Builder charge_builder;
    arrow::DoubleBuilder retention_time_builder;
    arrow::DoubleBuilder mass_to_charge_builder;
    arrow::DoubleBuilder score_builder;
    arrow::Int32Builder rank_builder;
    arrow::StringBuilder protein_accessions_builder;
    arrow::StringBuilder modifications_builder;
    arrow::BooleanBuilder is_decoy_builder;
    arrow::StringBuilder search_engine_builder;
    arrow::StringBuilder search_engine_score_name_builder;

    // Find associated protein identification for metadata
    std::map<String, const ProteinIdentification*> protein_id_map;
    for (const auto& protein_id : protein_identifications)
    {
      protein_id_map[protein_id.getIdentifier()] = &protein_id;
    }

    // Process each peptide identification
    for (const auto& peptide_id : peptide_identifications)
    {
      const ProteinIdentification* protein_id = nullptr;
      auto it = protein_id_map.find(peptide_id.getIdentifier());
      if (it != protein_id_map.end())
      {
        protein_id = it->second;
      }

      // Get spectrum reference
      String spectrum_ref = peptide_id.getSpectrumReference();
      if (spectrum_ref.empty())
      {
        // Use a fallback format if spectrum reference is not available
        spectrum_ref = peptide_id.getBaseName() + "_RT_" + String(peptide_id.getRT());
      }

      // Process each peptide hit
      const auto& hits = peptide_id.getHits();
      for (size_t rank = 0; rank < hits.size(); ++rank)
      {
        const PeptideHit& hit = hits[rank];

        // Sequence
        String sequence = hit.getSequence().toString();
        auto status = sequence_builder.Append(sequence.c_str());
        if (!status.ok()) {
          throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                          "Failed to append sequence: " + String(status.ToString()));
        }

        // Spectrum reference
        auto status2 = spectrum_reference_builder.Append(spectrum_ref.c_str());
        if (!status2.ok()) {
          throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                          "Failed to append spectrum reference: " + String(status2.ToString()));
        }

        // Charge
        auto status3 = charge_builder.Append(hit.getCharge());
        if (!status3.ok()) {
          throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                          "Failed to append charge: " + String(status3.ToString()));
        }

        // Retention time
        double rt = peptide_id.getRT();
        auto status4 = retention_time_builder.Append(rt);
        if (!status4.ok()) {
          throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                          "Failed to append retention time: " + String(status4.ToString()));
        }

        // Mass to charge
        double mz = peptide_id.getMZ();
        auto status5 = mass_to_charge_builder.Append(mz);
        if (!status5.ok()) {
          throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                          "Failed to append mass to charge: " + String(status5.ToString()));
        }

        // Score
        double score = hit.getScore();
        auto status6 = score_builder.Append(score);
        if (!status6.ok()) {
          throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                          "Failed to append score: " + String(status6.ToString()));
        }

        // Rank (1-based)
        auto status7 = rank_builder.Append(static_cast<int32_t>(rank + 1));
        if (!status7.ok()) {
          throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                          "Failed to append rank: " + String(status7.ToString()));
        }

        // Protein accessions
        String protein_accessions;
        const auto& peptide_evidences = hit.getPeptideEvidences();
        for (size_t i = 0; i < peptide_evidences.size(); ++i)
        {
          if (i > 0) protein_accessions += ",";
          protein_accessions += peptide_evidences[i].getProteinAccession();
        }
        auto status8 = protein_accessions_builder.Append(protein_accessions.c_str());
        if (!status8.ok()) {
          throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                          "Failed to append protein accessions: " + String(status8.ToString()));
        }

        // Modifications
        String modifications = hit.getSequence().toString(); // TODO: Extract actual modifications
        auto status9 = modifications_builder.Append(modifications.c_str());
        if (!status9.ok()) {
          throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                          "Failed to append modifications: " + String(status9.ToString()));
        }

        // Is decoy
        bool is_decoy = false;
        if (protein_id)
        {
          // Check if any protein accession is marked as decoy
          for (const auto& evidence : peptide_evidences)
          {
            String accession = evidence.getProteinAccession();
            // Common decoy prefixes
            if (accession.hasPrefix("DECOY_") || accession.hasPrefix("REV_") || 
                accession.hasPrefix("decoy_") || accession.hasPrefix("rev_"))
            {
              is_decoy = true;
              break;
            }
          }
        }
        auto status10 = is_decoy_builder.Append(is_decoy);
        if (!status10.ok()) {
          throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                          "Failed to append is_decoy: " + String(status10.ToString()));
        }

        // Search engine
        String search_engine = "Unknown";
        if (protein_id)
        {
          search_engine = protein_id->getSearchEngine();
          if (search_engine.empty()) search_engine = "Unknown";
        }
        auto status11 = search_engine_builder.Append(search_engine.c_str());
        if (!status11.ok()) {
          throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                          "Failed to append search engine: " + String(status11.ToString()));
        }

        // Search engine score name
        String score_name = peptide_id.getScoreType();
        if (score_name.empty()) score_name = "Unknown";
        auto status12 = search_engine_score_name_builder.Append(score_name.c_str());
        if (!status12.ok()) {
          throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                          "Failed to append search engine score name: " + String(status12.ToString()));
        }
      }
    }

    // Finish builders and create arrays
    std::shared_ptr<arrow::Array> sequence_array;
    auto status = sequence_builder.Finish(&sequence_array);
    if (!status.ok()) {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish sequence array: " + String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> spectrum_reference_array;
    status = spectrum_reference_builder.Finish(&spectrum_reference_array);
    if (!status.ok()) {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish spectrum reference array: " + String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> charge_array;
    status = charge_builder.Finish(&charge_array);
    if (!status.ok()) {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish charge array: " + String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> retention_time_array;
    status = retention_time_builder.Finish(&retention_time_array);
    if (!status.ok()) {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish retention time array: " + String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> mass_to_charge_array;
    status = mass_to_charge_builder.Finish(&mass_to_charge_array);
    if (!status.ok()) {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish mass to charge array: " + String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> score_array;
    status = score_builder.Finish(&score_array);
    if (!status.ok()) {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish score array: " + String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> rank_array;
    status = rank_builder.Finish(&rank_array);
    if (!status.ok()) {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish rank array: " + String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> protein_accessions_array;
    status = protein_accessions_builder.Finish(&protein_accessions_array);
    if (!status.ok()) {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish protein accessions array: " + String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> modifications_array;
    status = modifications_builder.Finish(&modifications_array);
    if (!status.ok()) {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish modifications array: " + String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> is_decoy_array;
    status = is_decoy_builder.Finish(&is_decoy_array);
    if (!status.ok()) {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish is_decoy array: " + String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> search_engine_array;
    status = search_engine_builder.Finish(&search_engine_array);
    if (!status.ok()) {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish search engine array: " + String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> search_engine_score_name_array;
    status = search_engine_score_name_builder.Finish(&search_engine_score_name_array);
    if (!status.ok()) {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish search engine score name array: " + String(status.ToString()));
    }

    // Create table
    auto schema = createPSMSchema_();
    auto table = arrow::Table::Make(schema, {
      sequence_array,
      spectrum_reference_array,
      charge_array,
      retention_time_array,
      mass_to_charge_array,
      score_array,
      rank_array,
      protein_accessions_array,
      modifications_array,
      is_decoy_array,
      search_engine_array,
      search_engine_score_name_array
    });

    return table;
  }

  void QuantmsIO::writeParquetFile_(const std::shared_ptr<arrow::Table>& table, 
                                   const String& filename)
  {
    std::shared_ptr<arrow::io::FileOutputStream> outfile;
    auto result = arrow::io::FileOutputStream::Open(filename.c_str());
    if (!result.ok()) {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                         filename, "Failed to open parquet file: " + String(result.status().ToString()));
    }
    outfile = result.ValueOrDie();

    auto status = parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), 
                                           outfile, table->num_rows());
    if (!status.ok()) {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                         filename, "Failed to write parquet file: " + String(status.ToString()));
    }
  }

} // namespace OpenMS

#else // WITH_PARQUET

namespace OpenMS
{
  QuantmsIO::QuantmsIO() = default;
  QuantmsIO::~QuantmsIO() = default;

  void QuantmsIO::store(const String& /*filename*/,
                       const std::vector<ProteinIdentification>& /*protein_identifications*/,
                       const PeptideIdentificationList& /*peptide_identifications*/)
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
  }
}

#endif // WITH_PARQUET