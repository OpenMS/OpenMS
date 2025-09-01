// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/QuantmsIO.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>

#include <algorithm>
#include <sstream>

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
      arrow::field("peptidoform", arrow::utf8()),
      arrow::field("modifications", arrow::null(), true), // nullable - null for now
      arrow::field("precursor_charge", arrow::int32()),
      arrow::field("posterior_error_probability", arrow::float32(), true), // nullable
      arrow::field("is_decoy", arrow::int32()),
      arrow::field("calculated_mz", arrow::float32()),
      arrow::field("observed_mz", arrow::float32()),
      arrow::field("additional_scores", arrow::null(), true), // nullable - null for now
      arrow::field("mp_accessions", arrow::null(), true), // nullable - null for now
      arrow::field("predicted_rt", arrow::float32(), true), // nullable
      arrow::field("reference_file_name", arrow::utf8()),
      arrow::field("cv_params", arrow::null(), true), // nullable - null for now
      arrow::field("scan", arrow::utf8()),
      arrow::field("rt", arrow::float32(), true), // nullable
      arrow::field("ion_mobility", arrow::float32(), true), // nullable
      arrow::field("num_peaks", arrow::int32(), true), // nullable
      arrow::field("mz_array", arrow::null(), true), // nullable - null for now
      arrow::field("intensity_array", arrow::null(), true) // nullable - null for now
    };

    return arrow::schema(fields);
  }

  std::shared_ptr<arrow::Table> QuantmsIO::convertToArrowTable_(
    const std::vector<ProteinIdentification>& protein_identifications,
    const PeptideIdentificationList& peptide_identifications)
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
    std::map<String, const ProteinIdentification*> protein_id_map;
    for (const auto& protein_id : protein_identifications)
    {
      protein_id_map[protein_id.getIdentifier()] = &protein_id;
    }

    // Find PEP score if available using ScoreSwitcher logic
    String pep_score_name;
    bool is_main_score_pep = false;
    if (!peptide_identifications.empty() && !peptide_identifications[0].getHits().empty())
    {
      const auto& first_peptide_id = peptide_identifications[0];
      const auto& first_hit = first_peptide_id.getHits()[0];
      
      // First check if main score is already a PEP score
      const String& score_type = first_peptide_id.getScoreType();
      const std::set<String> pep_score_types = {"Posterior Error Probability", "pep", "PEP", "MS:1001493"};
      for (const auto& pep_type : pep_score_types)
      {
        if (score_type.hasSubstring(pep_type))
        {
          is_main_score_pep = true;
          break;
        }
      }
      
      // If main score is not PEP, look for PEP score in metavalues
      if (!is_main_score_pep)
      {
        const std::set<String> pep_names = {"pep", "PEP", "Posterior Error Probability", "posterior_error_probability", "MS:1001493"};
        for (const auto& name : pep_names)
        {
          if (first_hit.metaValueExists(name))
          {
            pep_score_name = name;
            break;
          }
          if (first_hit.metaValueExists(name + "_score"))
          {
            pep_score_name = name + "_score";
            break;
          }
        }
      }
    }

    // Process each peptide identification (only first hit per peptide identification)
    for (const auto& peptide_id : peptide_identifications)
    {
      const auto& hits = peptide_id.getHits();
      if (hits.empty()) continue; // Skip if no hits
      
      // Only process the first hit (rank 1)
      const PeptideHit& hit = hits[0];

      // Sequence
      String sequence = hit.getSequence().toUnmodifiedString();
      auto status = sequence_builder.Append(sequence.c_str());
      if (!status.ok()) {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append sequence: " + String(status.ToString()));
      }

      // Peptidoform (sequence with modifications)
      String peptidoform = hit.getSequence().toString();
      status = peptidoform_builder.Append(peptidoform.c_str());
      if (!status.ok()) {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append peptidoform: " + String(status.ToString()));
      }

      // Precursor charge
      status = precursor_charge_builder.Append(hit.getCharge());
      if (!status.ok()) {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append precursor charge: " + String(status.ToString()));
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
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append posterior error probability: " + String(status.ToString()));
      }

      // Is decoy - check target_decoy metavalue
      int is_decoy = 0;
      if (hit.metaValueExists("target_decoy"))
      {
        String target_decoy = hit.getMetaValue("target_decoy").toString();
        is_decoy = (target_decoy == "decoy" || target_decoy == "DECOY") ? 1 : 0;
      }
      status = is_decoy_builder.Append(is_decoy);
      if (!status.ok()) {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append is_decoy: " + String(status.ToString()));
      }

      // Calculated m/z (theoretical)
      double theoretical_mz = hit.getSequence().getMonoWeight(Residue::Full, hit.getCharge());
      status = calculated_mz_builder.Append(static_cast<float>(theoretical_mz));
      if (!status.ok()) {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append calculated mz: " + String(status.ToString()));
      }

      // Observed m/z (experimental)
      double observed_mz = peptide_id.getMZ();
      status = observed_mz_builder.Append(static_cast<float>(observed_mz));
      if (!status.ok()) {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append observed mz: " + String(status.ToString()));
      }

      // Protein accessions (comma-separated string for now)
      String protein_accessions;
      const auto& peptide_evidences = hit.getPeptideEvidences();
      for (size_t i = 0; i < peptide_evidences.size(); ++i)
      {
        if (i > 0) protein_accessions += ",";
        protein_accessions += peptide_evidences[i].getProteinAccession();
      }
      status = mp_accessions_builder.Append(protein_accessions.c_str());
      if (!status.ok()) {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append protein accessions: " + String(status.ToString()));
      }

      // Predicted RT (null for now)
      status = predicted_rt_builder.AppendNull();
      if (!status.ok()) {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append predicted rt: " + String(status.ToString()));
      }

      // Reference file name
      String file_name = peptide_id.getSpectrumReference();
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
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append reference file name: " + String(status.ToString()));
      }

      // Scan
      String scan = peptide_id.getSpectrumReference();
      if (scan.empty())
      {
        // Generate scan from RT if available
        std::ostringstream scan_stream;
        scan_stream << "RT_" << peptide_id.getRT();
        scan = scan_stream.str();
      }
      status = scan_builder.Append(scan.c_str());
      if (!status.ok()) {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append scan: " + String(status.ToString()));
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
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append rt: " + String(status.ToString()));
      }

      // Ion mobility (null for now)
      status = ion_mobility_builder.AppendNull();
      if (!status.ok()) {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append ion mobility: " + String(status.ToString()));
      }

      // Num peaks (null for now)
      status = num_peaks_builder.AppendNull();
      if (!status.ok()) {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to append num peaks: " + String(status.ToString()));
      }
    }

    // Finish builders and create arrays
    std::shared_ptr<arrow::Array> sequence_array;
    auto status = sequence_builder.Finish(&sequence_array);
    if (!status.ok()) {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish sequence array: " + String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> peptidoform_array;
    status = peptidoform_builder.Finish(&peptidoform_array);
    if (!status.ok()) {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish peptidoform array: " + String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> precursor_charge_array;
    status = precursor_charge_builder.Finish(&precursor_charge_array);
    if (!status.ok()) {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish precursor charge array: " + String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> posterior_error_probability_array;
    status = posterior_error_probability_builder.Finish(&posterior_error_probability_array);
    if (!status.ok()) {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish posterior error probability array: " + String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> is_decoy_array;
    status = is_decoy_builder.Finish(&is_decoy_array);
    if (!status.ok()) {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish is_decoy array: " + String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> calculated_mz_array;
    status = calculated_mz_builder.Finish(&calculated_mz_array);
    if (!status.ok()) {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish calculated mz array: " + String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> observed_mz_array;
    status = observed_mz_builder.Finish(&observed_mz_array);
    if (!status.ok()) {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish observed mz array: " + String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> mp_accessions_array;
    status = mp_accessions_builder.Finish(&mp_accessions_array);
    if (!status.ok()) {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish mp_accessions array: " + String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> predicted_rt_array;
    status = predicted_rt_builder.Finish(&predicted_rt_array);
    if (!status.ok()) {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish predicted rt array: " + String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> reference_file_name_array;
    status = reference_file_name_builder.Finish(&reference_file_name_array);
    if (!status.ok()) {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish reference file name array: " + String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> scan_array;
    status = scan_builder.Finish(&scan_array);
    if (!status.ok()) {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish scan array: " + String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> rt_array;
    status = rt_builder.Finish(&rt_array);
    if (!status.ok()) {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish rt array: " + String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> ion_mobility_array;
    status = ion_mobility_builder.Finish(&ion_mobility_array);
    if (!status.ok()) {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish ion mobility array: " + String(status.ToString()));
    }

    std::shared_ptr<arrow::Array> num_peaks_array;
    status = num_peaks_builder.Finish(&num_peaks_array);
    if (!status.ok()) {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Failed to finish num peaks array: " + String(status.ToString()));
    }

    // Create simplified schema for now - will add complex nested types later
    auto schema = createPSMSchema_();
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
      arrow::MakeArrayOfNull(arrow::null(), sequence_array->length()).ValueOrDie(), // mp_accessions - using mp_accessions_array would need list type
      predicted_rt_array,
      reference_file_name_array,
      arrow::MakeArrayOfNull(arrow::null(), sequence_array->length()).ValueOrDie(), // cv_params - null for now
      scan_array,
      rt_array,
      ion_mobility_array,
      num_peaks_array,
      arrow::MakeArrayOfNull(arrow::null(), sequence_array->length()).ValueOrDie(), // mz_array - null for now
      arrow::MakeArrayOfNull(arrow::null(), sequence_array->length()).ValueOrDie()  // intensity_array - null for now
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