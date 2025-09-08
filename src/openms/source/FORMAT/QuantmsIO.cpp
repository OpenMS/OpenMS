// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/QuantmsIO.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>

#include <algorithm>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <functional>

using namespace std;

namespace
{
  // Helper functions moved to anonymous namespace to hide Arrow types from header

  std::shared_ptr<arrow::Schema> createPSMSchema()
  {
    std::vector<std::shared_ptr<arrow::Field>> fields = {
      arrow::field("sequence", arrow::utf8()),
      arrow::field("peptidoform", arrow::utf8()),
      arrow::field("modifications", arrow::null(), true), // nullable - null for now
      arrow::field("precursor_charge", arrow::int32()),
      arrow::field("posterior_error_probability", arrow::float32(), true), // nullable
      arrow::field("is_decoy", arrow::int32()),
      arrow::field("calculated_mz", arrow::float32()),
      arrow::field("observed_mz", arrow::float32()),
      arrow::field("additional_scores", arrow::null(), true), // nullable - null for now
      arrow::field("protein_accessions", arrow::null(), true), // nullable - null for now
      arrow::field("predicted_rt", arrow::float32(), true), // nullable
      arrow::field("reference_file_name", arrow::utf8()),
      arrow::field("cv_params", arrow::null(), true), // nullable - null for now
      arrow::field("scan", arrow::utf8()),
      arrow::field("rt", arrow::float32(), true), // nullable
      arrow::field("ion_mobility", arrow::float32(), true), // nullable
      arrow::field("number_peaks", arrow::int32(), true), // nullable
      arrow::field("mz_array", arrow::null(), true), // nullable - null for now
      arrow::field("intensity_array", arrow::null(), true) // nullable - null for now
    };

    return arrow::schema(fields);
  }

  void writeParquetFile(const std::shared_ptr<arrow::Table>& table, 
                       const OpenMS::String& filename,
                       const std::map<std::string, std::string>& file_metadata = {})
  {
    std::shared_ptr<arrow::io::FileOutputStream> outfile;
    auto result = arrow::io::FileOutputStream::Open(filename.c_str());
    if (!result.ok()) {
      throw OpenMS::Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                         filename, "Failed to open parquet file: " + OpenMS::String(result.status().ToString()));
    }
    outfile = result.ValueOrDie();

    // For now, let's use the basic WriteTable approach and add metadata support later
    // when proper Parquet metadata headers are available
    auto status = parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), 
                                           outfile, table->num_rows());
    if (!status.ok()) {
      throw OpenMS::Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                         filename, "Failed to write parquet file: " + OpenMS::String(status.ToString()));
    }
    
    // TODO: Add file metadata to parquet when appropriate headers are available
    // The metadata is currently passed but not yet written to the file
    // This will need parquet::WriterProperties or similar API
  }

  std::shared_ptr<arrow::Table> convertToArrowTable(
    const std::vector<OpenMS::ProteinIdentification>& protein_identifications,
    const OpenMS::PeptideIdentificationList& peptide_identifications)
  {
    // Create builders for each column - using simpler types for now
    arrow::StringBuilder sequence_builder;
    arrow::StringBuilder peptidoform_builder;
    arrow::Int32Builder precursor_charge_builder;
    arrow::FloatBuilder posterior_error_probability_builder;
    arrow::Int32Builder is_decoy_builder;
    arrow::FloatBuilder calculated_mz_builder;
    arrow::FloatBuilder observed_mz_builder;
    arrow::StringBuilder mp_accessions_builder; // Use comma-separated string for now
    arrow::FloatBuilder predicted_rt_builder;
    arrow::StringBuilder reference_file_name_builder;
    arrow::StringBuilder scan_builder;
    arrow::FloatBuilder rt_builder;
    arrow::FloatBuilder ion_mobility_builder;
    arrow::Int32Builder num_peaks_builder;

    // Find associated protein identification for metadata
    std::map<OpenMS::String, const OpenMS::ProteinIdentification*> protein_id_map;
    for (const auto& protein_id : protein_identifications)
    {
      protein_id_map[protein_id.getIdentifier()] = &protein_id;
    }

    // Find PEP score if available using IDScoreSwitcherAlgorithm
    OpenMS::String pep_score_name;
    bool is_main_score_pep = false;
    if (!peptide_identifications.empty())
    {
      OpenMS::IDScoreSwitcherAlgorithm score_switcher;
      const auto& first_peptide_id = peptide_identifications[0];
      auto pep_result = score_switcher.findScoreType(first_peptide_id, OpenMS::IDScoreSwitcherAlgorithm::ScoreType::PEP);
      
      pep_score_name = pep_result.score_name;
      is_main_score_pep = pep_result.is_main_score_type;
    }

    // Process each peptide identification (only first hit per peptide identification)
    for (const auto& peptide_id : peptide_identifications)
    {
      const auto& hits = peptide_id.getHits();
      if (hits.empty()) continue; // Skip if no hits
      
      // Only process the first hit (rank 1)
      const OpenMS::PeptideHit& hit = hits[0];

      // Sequence
      OpenMS::String sequence = hit.getSequence().toUnmodifiedString();
      auto status = sequence_builder.Append(sequence.c_str());
      if (!status.ok()) {
        throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append sequence: " + OpenMS::String(status.ToString()));
      }

      // Peptidoform (sequence with modifications in ProForma format)
      OpenMS::String peptidoform = hit.getSequence().toBracketString(true, false);
      // Convert round brackets to square brackets for ProForma format
      peptidoform.substitute("(", "[");
      peptidoform.substitute(")", "]");
      status = peptidoform_builder.Append(peptidoform.c_str());
      if (!status.ok()) {
        throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append peptidoform: " + OpenMS::String(status.ToString()));
      }

      // Precursor charge
      status = precursor_charge_builder.Append(hit.getCharge());
      if (!status.ok()) {
        throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append precursor charge: " + OpenMS::String(status.ToString()));
      }

      // Posterior error probability
      if (is_main_score_pep)
      {
        // Use main score as PEP value
        double pep_value = hit.getScore();
        status = posterior_error_probability_builder.Append(static_cast<float>(pep_value));
      }
      else if (!pep_score_name.empty() && hit.metaValueExists(pep_score_name))
      {
        double pep_value = hit.getMetaValue(pep_score_name);
        status = posterior_error_probability_builder.Append(static_cast<float>(pep_value));
      }
      else
      {
        status = posterior_error_probability_builder.AppendNull();
      }
      if (!status.ok()) {
        throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append posterior error probability: " + OpenMS::String(status.ToString()));
      }

      // Is decoy - use the isDecoy() method from PeptideHit
      int is_decoy = hit.isDecoy() ? 1 : 0;
      status = is_decoy_builder.Append(is_decoy);
      if (!status.ok()) {
        throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append is_decoy: " + OpenMS::String(status.ToString()));
      }

      // Calculated m/z (theoretical)
      double theoretical_mz = hit.getSequence().getMonoWeight(OpenMS::Residue::Full, hit.getCharge());
      status = calculated_mz_builder.Append(static_cast<float>(theoretical_mz));
      if (!status.ok()) {
        throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append calculated mz: " + OpenMS::String(status.ToString()));
      }

      // Observed m/z (experimental)
      double observed_mz = peptide_id.getMZ();
      status = observed_mz_builder.Append(static_cast<float>(observed_mz));
      if (!status.ok()) {
        throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append observed mz: " + OpenMS::String(status.ToString()));
      }

      // Protein accessions (comma-separated string for now)
      OpenMS::String protein_accessions;
      const auto& peptide_evidences = hit.getPeptideEvidences();
      for (size_t i = 0; i < peptide_evidences.size(); ++i)
      {
        if (i > 0) protein_accessions += ",";
        protein_accessions += peptide_evidences[i].getProteinAccession();
      }
      status = mp_accessions_builder.Append(protein_accessions.c_str());
      if (!status.ok()) {
        throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append protein accessions: " + OpenMS::String(status.ToString()));
      }

      // Predicted RT (null for now)
      status = predicted_rt_builder.AppendNull();
      if (!status.ok()) {
        throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append predicted rt: " + OpenMS::String(status.ToString()));
      }

      // Reference file name
      OpenMS::String file_name = peptide_id.getSpectrumReference();
      if (file_name.empty())
      {
        file_name = peptide_id.getBaseName();
      }
      if (file_name.empty())
      {
        file_name = "unknown";
      }
      status = reference_file_name_builder.Append(file_name.c_str());
      if (!status.ok()) {
        throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append reference file name: " + OpenMS::String(status.ToString()));
      }

      // Scan
      OpenMS::String scan = peptide_id.getSpectrumReference();
      if (scan.empty())
      {
        // Generate scan from RT if available
        std::ostringstream scan_stream;
        scan_stream << "RT_" << peptide_id.getRT();
        scan = scan_stream.str();
      }
      status = scan_builder.Append(scan.c_str());
      if (!status.ok()) {
        throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append scan: " + OpenMS::String(status.ToString()));
      }

      // RT (retention time)
      double rt = peptide_id.getRT();
      if (rt >= 0)
      {
        status = rt_builder.Append(static_cast<float>(rt));
      }
      else
      {
        status = rt_builder.AppendNull();
      }
      if (!status.ok()) {
        throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append rt: " + OpenMS::String(status.ToString()));
      }

      // Ion mobility 
      if (hit.metaValueExists(OpenMS::Constants::UserParam::IM))
      {
        double ion_mobility = hit.getMetaValue(OpenMS::Constants::UserParam::IM);
        status = ion_mobility_builder.Append(ion_mobility);
      }
      else
      {
        status = ion_mobility_builder.AppendNull();
      }
      if (!status.ok()) {
        throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append ion mobility: " + OpenMS::String(status.ToString()));
      }

      // Num peaks
      if (hit.metaValueExists(OpenMS::Constants::UserParam::NUM_PEAKS))
      {
        int num_peaks = hit.getMetaValue(OpenMS::Constants::UserParam::NUM_PEAKS);
        status = num_peaks_builder.Append(num_peaks);
      }
      else
      {
        status = num_peaks_builder.AppendNull();
      }
      if (!status.ok()) {
        throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append num peaks: " + OpenMS::String(status.ToString()));
      }
    }

    // Finish builders and create arrays
    std::shared_ptr<arrow::Array> sequence_array;
    auto status = sequence_builder.Finish(&sequence_array);
    if (!status.ok()) {
      throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish sequence array: " + OpenMS::String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> peptidoform_array;
    status = peptidoform_builder.Finish(&peptidoform_array);
    if (!status.ok()) {
      throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish peptidoform array: " + OpenMS::String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> precursor_charge_array;
    status = precursor_charge_builder.Finish(&precursor_charge_array);
    if (!status.ok()) {
      throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish precursor charge array: " + OpenMS::String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> posterior_error_probability_array;
    status = posterior_error_probability_builder.Finish(&posterior_error_probability_array);
    if (!status.ok()) {
      throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish posterior error probability array: " + OpenMS::String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> is_decoy_array;
    status = is_decoy_builder.Finish(&is_decoy_array);
    if (!status.ok()) {
      throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish is_decoy array: " + OpenMS::String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> calculated_mz_array;
    status = calculated_mz_builder.Finish(&calculated_mz_array);
    if (!status.ok()) {
      throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish calculated mz array: " + OpenMS::String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> observed_mz_array;
    status = observed_mz_builder.Finish(&observed_mz_array);
    if (!status.ok()) {
      throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish observed mz array: " + OpenMS::String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> mp_accessions_array;
    status = mp_accessions_builder.Finish(&mp_accessions_array);
    if (!status.ok()) {
      throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish mp_accessions array: " + OpenMS::String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> predicted_rt_array;
    status = predicted_rt_builder.Finish(&predicted_rt_array);
    if (!status.ok()) {
      throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish predicted rt array: " + OpenMS::String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> reference_file_name_array;
    status = reference_file_name_builder.Finish(&reference_file_name_array);
    if (!status.ok()) {
      throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish reference file name array: " + OpenMS::String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> scan_array;
    status = scan_builder.Finish(&scan_array);
    if (!status.ok()) {
      throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish scan array: " + OpenMS::String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> rt_array;
    status = rt_builder.Finish(&rt_array);
    if (!status.ok()) {
      throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish rt array: " + OpenMS::String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> ion_mobility_array;
    status = ion_mobility_builder.Finish(&ion_mobility_array);
    if (!status.ok()) {
      throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish ion mobility array: " + OpenMS::String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> num_peaks_array;
    status = num_peaks_builder.Finish(&num_peaks_array);
    if (!status.ok()) {
      throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish num peaks array: " + OpenMS::String(status.ToString()));
    }

    // Create simplified schema for now - will add complex nested types later
    auto schema = createPSMSchema();
    auto table = arrow::Table::Make(schema, {
      sequence_array,
      peptidoform_array,
      arrow::MakeArrayOfNull(arrow::null(), sequence_array->length()).ValueOrDie(), // modifications - null for now
      precursor_charge_array,
      posterior_error_probability_array,
      is_decoy_array,
      calculated_mz_array,
      observed_mz_array,
      arrow::MakeArrayOfNull(arrow::null(), sequence_array->length()).ValueOrDie(), // additional_scores - null for now
      arrow::MakeArrayOfNull(arrow::null(), sequence_array->length()).ValueOrDie(), // protein_accessions - using mp_accessions_array would need list type
      predicted_rt_array,
      reference_file_name_array,
      arrow::MakeArrayOfNull(arrow::null(), sequence_array->length()).ValueOrDie(), // cv_params - null for now
      scan_array,
      rt_array,
      ion_mobility_array,
      num_peaks_array,
      arrow::MakeArrayOfNull(arrow::null(), sequence_array->length()).ValueOrDie(), // mz_array - null for now
      arrow::MakeArrayOfNull(arrow::null(), sequence_array->length()).ValueOrDie() // intensity_array - null for now
    });

    return table;
  }

} // anonymous namespace

namespace OpenMS
{
  QuantmsIO::~QuantmsIO() = default;

  void QuantmsIO::store(const String& filename,
                       const std::vector<ProteinIdentification>& protein_identifications,
                       const PeptideIdentificationList& peptide_identifications)
  {
    // Generate file metadata
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    std::ostringstream creation_date_stream;
    creation_date_stream << std::put_time(std::gmtime(&time_t), "%Y-%m-%dT%H:%M:%SZ");
    std::string creation_date_str = creation_date_stream.str();
    
    // Generate a simple UUID based on current time and process
    std::ostringstream uuid_stream;
    uuid_stream << std::hex << std::hash<std::string>{}(creation_date_str) << "-0000-4000-8000-" 
                << std::setfill('0') << std::setw(12) << (std::hash<void*>{}(&protein_identifications) & 0xFFFFFFFFFFFF);
    std::string uuid_str = uuid_stream.str();

    // Create file metadata map
    std::map<std::string, std::string> file_metadata = {
      {"quantmsio_version", "1.0"},
      {"creator", "OpenMS"},
      {"file_type", "psm"},
      {"creation_date", creation_date_str},
      {"uuid", uuid_str},
      {"scan_format", "scan"},
      {"software_provider", "OpenMS"}
    };

    // Convert data to Arrow table
    auto table = convertToArrowTable(protein_identifications, peptide_identifications);
    
    // Write to parquet file with metadata
    writeParquetFile(table, filename, file_metadata);
  }

} // namespace OpenMS