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
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>

#include <algorithm>
#include <sstream>
#include <chrono>
#include <format>
#include <functional>

using namespace std;

namespace
{
  // Helper functions moved to anonymous namespace to hide Arrow types from header

  // Helper function to map DataValue type to Arrow field type
  std::shared_ptr<arrow::DataType> dataValueTypeToArrowType(OpenMS::DataValue::DataType data_type)
  {
    switch (data_type)
    {
      case OpenMS::DataValue::STRING_VALUE:
        return arrow::utf8();
      case OpenMS::DataValue::INT_VALUE:
        return arrow::int64(); // SignedSize maps to int64
      case OpenMS::DataValue::DOUBLE_VALUE:
        return arrow::float64();
      case OpenMS::DataValue::STRING_LIST:
        return arrow::list(arrow::utf8());
      case OpenMS::DataValue::INT_LIST:
        return arrow::list(arrow::int64());
      case OpenMS::DataValue::DOUBLE_LIST:
        return arrow::list(arrow::float64());
      case OpenMS::DataValue::EMPTY_VALUE:
      default:
        return arrow::utf8(); // Default to string for unknown or empty types
    }
  }

  // Helper function to determine meta value types by scanning all peptide hits
  std::map<OpenMS::String, OpenMS::DataValue::DataType> determineMetaValueTypes(
    const OpenMS::PeptideIdentificationList& peptide_identifications,
    const std::set<OpenMS::String>& meta_value_keys,
    bool export_all_psms)
  {
    std::map<OpenMS::String, OpenMS::DataValue::DataType> meta_value_types;
    
    // Initialize all keys as EMPTY_VALUE
    for (const auto& key : meta_value_keys)
    {
      meta_value_types[key] = OpenMS::DataValue::EMPTY_VALUE;
    }
    
    // Scan all peptide identifications and hits to determine types
    for (const auto& peptide_id : peptide_identifications)
    {
      const auto& hits = peptide_id.getHits();
      if (hits.empty()) continue;
      
      size_t num_hits_to_process = export_all_psms ? hits.size() : 1;
      
      for (size_t hit_index = 0; hit_index < num_hits_to_process; ++hit_index)
      {
        const OpenMS::PeptideHit& hit = hits[hit_index];
        
        for (const auto& key : meta_value_keys)
        {
          if (hit.metaValueExists(key))
          {
            OpenMS::DataValue meta_value = hit.getMetaValue(key);
            OpenMS::DataValue::DataType current_type = meta_value.valueType();
            
            if (current_type != OpenMS::DataValue::EMPTY_VALUE)
            {
              if (meta_value_types[key] == OpenMS::DataValue::EMPTY_VALUE)
              {
                // First non-empty value found, set the type
                meta_value_types[key] = current_type;
              }
              else if (meta_value_types[key] != current_type)
              {
                // Type conflict - throw an exception as requested
                throw OpenMS::Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                  "Meta value type conflict for key '" + key + "': found both " + 
                  OpenMS::DataValue::NamesOfDataType[meta_value_types[key]] + " and " + 
                  OpenMS::DataValue::NamesOfDataType[current_type] + ". All values for a meta key must have the same type.");
              }
            }
          }
        }
      }
    }
    
    // Convert any remaining EMPTY_VALUE types to STRING_VALUE
    for (auto& pair : meta_value_types)
    {
      if (pair.second == OpenMS::DataValue::EMPTY_VALUE)
      {
        pair.second = OpenMS::DataValue::STRING_VALUE;
      }
    }
    
    return meta_value_types;
  }

  // Base class for meta value builders to handle different types
  class MetaValueBuilderBase 
  {
  public:
    virtual ~MetaValueBuilderBase() = default;
    virtual arrow::Status AppendValue(const OpenMS::DataValue& value) = 0;
    virtual arrow::Status AppendNull() = 0;
    virtual arrow::Status Finish(std::shared_ptr<arrow::Array>* out) = 0;
  };

  // String meta value builder
  class StringMetaValueBuilder : public MetaValueBuilderBase
  {
  private:
    arrow::StringBuilder builder_;
  public:
    arrow::Status AppendValue(const OpenMS::DataValue& value) override
    {
      return builder_.Append(value.toString().c_str());
    }
    arrow::Status AppendNull() override
    {
      return builder_.AppendNull();
    }
    arrow::Status Finish(std::shared_ptr<arrow::Array>* out) override
    {
      return builder_.Finish(out);
    }
  };

  // Integer meta value builder
  class IntMetaValueBuilder : public MetaValueBuilderBase
  {
  private:
    arrow::Int64Builder builder_;
  public:
    arrow::Status AppendValue(const OpenMS::DataValue& value) override
    {
      try
      {
        return builder_.Append(static_cast<int64_t>(value));
      }
      catch (...)
      {
        return builder_.AppendNull(); // Fall back to null if conversion fails
      }
    }
    arrow::Status AppendNull() override
    {
      return builder_.AppendNull();
    }
    arrow::Status Finish(std::shared_ptr<arrow::Array>* out) override
    {
      return builder_.Finish(out);
    }
  };

  // Double meta value builder
  class DoubleMetaValueBuilder : public MetaValueBuilderBase
  {
  private:
    arrow::DoubleBuilder builder_;
  public:
    arrow::Status AppendValue(const OpenMS::DataValue& value) override
    {
      try
      {
        return builder_.Append(static_cast<double>(value));
      }
      catch (...)
      {
        return builder_.AppendNull(); // Fall back to null if conversion fails
      }
    }
    arrow::Status AppendNull() override
    {
      return builder_.AppendNull();
    }
    arrow::Status Finish(std::shared_ptr<arrow::Array>* out) override
    {
      return builder_.Finish(out);
    }
  };

  // String list meta value builder
  class StringListMetaValueBuilder : public MetaValueBuilderBase
  {
  private:
    arrow::ListBuilder builder_;
    arrow::StringBuilder value_builder_;
  public:
    StringListMetaValueBuilder() : builder_(arrow::default_memory_pool(), std::make_shared<arrow::StringBuilder>())
    {
    }
    
    arrow::Status AppendValue(const OpenMS::DataValue& value) override
    {
      try
      {
        auto string_list = value.toStringList();
        ARROW_RETURN_NOT_OK(builder_.Append());
        for (const auto& str : string_list)
        {
          ARROW_RETURN_NOT_OK(static_cast<arrow::StringBuilder*>(builder_.value_builder())->Append(str.c_str()));
        }
        return arrow::Status::OK();
      }
      catch (...)
      {
        return builder_.AppendNull(); // Fall back to null if conversion fails
      }
    }
    arrow::Status AppendNull() override
    {
      return builder_.AppendNull();
    }
    arrow::Status Finish(std::shared_ptr<arrow::Array>* out) override
    {
      return builder_.Finish(out);
    }
  };

  // Integer list meta value builder
  class IntListMetaValueBuilder : public MetaValueBuilderBase
  {
  private:
    arrow::ListBuilder builder_;
  public:
    IntListMetaValueBuilder() : builder_(arrow::default_memory_pool(), std::make_shared<arrow::Int64Builder>())
    {
    }
    
    arrow::Status AppendValue(const OpenMS::DataValue& value) override
    {
      try
      {
        auto int_list = value.toIntList();
        ARROW_RETURN_NOT_OK(builder_.Append());
        for (const auto& val : int_list)
        {
          ARROW_RETURN_NOT_OK(static_cast<arrow::Int64Builder*>(builder_.value_builder())->Append(static_cast<int64_t>(val)));
        }
        return arrow::Status::OK();
      }
      catch (...)
      {
        return builder_.AppendNull(); // Fall back to null if conversion fails
      }
    }
    arrow::Status AppendNull() override
    {
      return builder_.AppendNull();
    }
    arrow::Status Finish(std::shared_ptr<arrow::Array>* out) override
    {
      return builder_.Finish(out);
    }
  };

  // Double list meta value builder
  class DoubleListMetaValueBuilder : public MetaValueBuilderBase
  {
  private:
    arrow::ListBuilder builder_;
  public:
    DoubleListMetaValueBuilder() : builder_(arrow::default_memory_pool(), std::make_shared<arrow::DoubleBuilder>())
    {
    }
    
    arrow::Status AppendValue(const OpenMS::DataValue& value) override
    {
      try
      {
        auto double_list = value.toDoubleList();
        ARROW_RETURN_NOT_OK(builder_.Append());
        for (const auto& val : double_list)
        {
          ARROW_RETURN_NOT_OK(static_cast<arrow::DoubleBuilder*>(builder_.value_builder())->Append(val));
        }
        return arrow::Status::OK();
      }
      catch (...)
      {
        return builder_.AppendNull(); // Fall back to null if conversion fails
      }
    }
    arrow::Status AppendNull() override
    {
      return builder_.AppendNull();
    }
    arrow::Status Finish(std::shared_ptr<arrow::Array>* out) override
    {
      return builder_.Finish(out);
    }
  };

  // Factory function to create appropriate meta value builder
  std::unique_ptr<MetaValueBuilderBase> createMetaValueBuilder(OpenMS::DataValue::DataType data_type)
  {
    switch (data_type)
    {
      case OpenMS::DataValue::STRING_VALUE:
        return std::make_unique<StringMetaValueBuilder>();
      case OpenMS::DataValue::INT_VALUE:
        return std::make_unique<IntMetaValueBuilder>();
      case OpenMS::DataValue::DOUBLE_VALUE:
        return std::make_unique<DoubleMetaValueBuilder>();
      case OpenMS::DataValue::STRING_LIST:
        return std::make_unique<StringListMetaValueBuilder>();
      case OpenMS::DataValue::INT_LIST:
        return std::make_unique<IntListMetaValueBuilder>();
      case OpenMS::DataValue::DOUBLE_LIST:
        return std::make_unique<DoubleListMetaValueBuilder>();
      case OpenMS::DataValue::EMPTY_VALUE:
      default:
        return std::make_unique<StringMetaValueBuilder>(); // Default to string
    }
  }

  std::shared_ptr<arrow::Schema> createPSMSchema(bool export_all_psms = false, 
                                                const std::set<OpenMS::String>& meta_value_keys = {},
                                                const std::map<OpenMS::String, OpenMS::DataValue::DataType>& meta_value_types = {})
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

    // Add rank column if exporting all PSMs
    if (export_all_psms)
    {
      fields.insert(fields.begin() + 4, arrow::field("rank", arrow::int32())); // Insert after precursor_charge
    }

    // Add meta value columns with appropriate types
    for (const auto& key : meta_value_keys)
    {
      auto type_it = meta_value_types.find(key);
      std::shared_ptr<arrow::DataType> arrow_type;
      
      if (type_it != meta_value_types.end())
      {
        arrow_type = dataValueTypeToArrowType(type_it->second);
      }
      else
      {
        arrow_type = arrow::utf8(); // Default to string
      }
      
      fields.push_back(arrow::field(key.c_str(), arrow_type, true)); // nullable
    }

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

    // Create a writer with parquet::arrow::FileWriter
    auto writer_result = parquet::arrow::FileWriter::Open(*table->schema(), 
                                                       arrow::default_memory_pool(),
                                                       outfile);
    if (!writer_result.ok()) {
      throw OpenMS::Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      filename, "Failed to create parquet writer: " + OpenMS::String(writer_result.status().ToString()));
    }
    std::unique_ptr<parquet::arrow::FileWriter> writer = std::move(writer_result.ValueOrDie());
    
    // Add metadata to the parquet file using AddKeyValueMetadata
    if (!file_metadata.empty()) {
      std::vector<std::string> md_keys, md_values;
      for (const auto& kv : file_metadata) {
        md_keys.push_back(kv.first);
        md_values.push_back(kv.second);
      }
      const auto metadata = arrow::key_value_metadata(md_keys,md_values);
      auto status = writer->AddKeyValueMetadata(metadata);
      if (!status.ok()) {
        throw OpenMS::Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        filename, "Failed to add metadata: " + OpenMS::String(status.ToString()));
      }
    }
    
    
    // Write table using the FileWriter
    auto status = writer->WriteTable(*table, table->num_rows());
    if (!status.ok()) {
      throw OpenMS::Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                         filename, "Failed to write parquet file: " + OpenMS::String(status.ToString()));
    }
    // Close writer
    status = writer->Close();
    if (!status.ok()) {
      throw OpenMS::Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      filename, "Failed to close writer: " + OpenMS::String(status.ToString()));
    }
  }

  std::shared_ptr<arrow::Table> convertToArrowTable(
    const std::vector<OpenMS::ProteinIdentification>& protein_identifications,
    const OpenMS::PeptideIdentificationList& peptide_identifications,
    bool export_all_psms = false,
    const std::set<OpenMS::String>& meta_value_keys = {})
  {
    // First, determine meta value types by scanning all data
    auto meta_value_types = determineMetaValueTypes(peptide_identifications, meta_value_keys, export_all_psms);

    // Create builders for each column - using simpler types for now
    arrow::StringBuilder sequence_builder;
    arrow::StringBuilder peptidoform_builder;
    arrow::Int32Builder precursor_charge_builder;
    arrow::Int32Builder rank_builder; // Only used if export_all_psms is true
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

    // Create typed builders for meta value columns
    std::map<OpenMS::String, std::unique_ptr<MetaValueBuilderBase>> meta_value_builders;
    for (const auto& key : meta_value_keys)
    {
      auto type_it = meta_value_types.find(key);
      OpenMS::DataValue::DataType data_type = (type_it != meta_value_types.end()) ? 
                                             type_it->second : OpenMS::DataValue::STRING_VALUE;
      meta_value_builders[key] = createMetaValueBuilder(data_type);
    }

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

    // Process each peptide identification
    for (const auto& peptide_id : peptide_identifications)
    {
      const auto& hits = peptide_id.getHits();
      if (hits.empty()) continue; // Skip if no hits
      
      // Determine how many hits to process
      size_t num_hits_to_process = export_all_psms ? hits.size() : 1;
      
      for (size_t hit_index = 0; hit_index < num_hits_to_process; ++hit_index)
      {
        const OpenMS::PeptideHit& hit = hits[hit_index];

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

        // Rank (if exporting all PSMs)
        if (export_all_psms)
        {
          status = rank_builder.Append(static_cast<int>(hit_index + 1)); // rank is 1-based
          if (!status.ok()) {
            throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                            "Failed to append rank: " + OpenMS::String(status.ToString()));
          }
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
        scan = std::format("RT_{}", peptide_id.getRT());
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

      // Process meta values
      for (const auto& key : meta_value_keys)
      {
        auto& builder = meta_value_builders[key];
        if (hit.metaValueExists(key))
        {
          OpenMS::DataValue meta_value = hit.getMetaValue(key);
          status = builder->AppendValue(meta_value);
        }
        else
        {
          status = builder->AppendNull();
        }
        if (!status.ok()) {
          throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                          "Failed to append meta value " + key + ": " + OpenMS::String(status.ToString()));
        }
      }
      } // End hit processing loop
    } // End peptide identification loop

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

    // Finish rank array if needed
    std::shared_ptr<arrow::Array> rank_array;
    if (export_all_psms)
    {
      status = rank_builder.Finish(&rank_array);
      if (!status.ok()) {
        throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to finish rank array: " + OpenMS::String(status.ToString()));
      }
    }

    // Finish meta value arrays
    std::map<OpenMS::String, std::shared_ptr<arrow::Array>> meta_value_arrays;
    for (const auto& key : meta_value_keys)
    {
      std::shared_ptr<arrow::Array> meta_array;
      status = meta_value_builders[key]->Finish(&meta_array);
      if (!status.ok()) {
        throw OpenMS::Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Failed to finish meta value array for " + key + ": " + OpenMS::String(status.ToString()));
      }
      meta_value_arrays[key] = meta_array;
    }

    // Create simplified schema for now - will add complex nested types later
    auto schema = createPSMSchema(export_all_psms, meta_value_keys, meta_value_types);
    
    // Build the arrays vector in the correct order to match the schema
    std::vector<std::shared_ptr<arrow::Array>> arrays = {
      sequence_array,
      peptidoform_array,
      arrow::MakeArrayOfNull(arrow::null(), sequence_array->length()).ValueOrDie(), // modifications - null for now
      precursor_charge_array
    };
    
    // Add rank array if exporting all PSMs (inserted after precursor_charge)
    if (export_all_psms)
    {
      arrays.push_back(rank_array);
    }
    
    // Continue with the rest of the standard columns
    arrays.insert(arrays.end(), {
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
    
    // Add meta value arrays in the same order as in the schema
    for (const auto& key : meta_value_keys)
    {
      arrays.push_back(meta_value_arrays[key]);
    }
    
    auto table = arrow::Table::Make(schema, arrays);

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
    store(filename, protein_identifications, peptide_identifications, false, {});
  }

  void QuantmsIO::store(const String& filename,
                       const std::vector<ProteinIdentification>& protein_identifications,
                       const PeptideIdentificationList& peptide_identifications,
                       bool export_all_psms)
  {
    store(filename, protein_identifications, peptide_identifications, export_all_psms, {});
  }

  void QuantmsIO::store(const String& filename,
                       const std::vector<ProteinIdentification>& protein_identifications,
                       const PeptideIdentificationList& peptide_identifications,
                       bool export_all_psms,
                       const std::set<String>& meta_value_keys)
  {
    // Generate file metadata
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    std::string creation_date_str = std::format("{:%Y-%m-%dT%H:%M:%SZ}", 
                                                std::chrono::system_clock::from_time_t(time_t));
    
    // Generate a simple UUID based on current time and process
    std::string uuid_str = std::format("{:x}-0000-4000-8000-{:012x}", 
                                       std::hash<std::string>{}(creation_date_str),
                                       std::hash<const void*>{}(&protein_identifications) & 0xFFFFFFFFFFFF);

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
    auto table = convertToArrowTable(protein_identifications, peptide_identifications, export_all_psms, meta_value_keys);
    
    // Write to parquet file with metadata
    writeParquetFile(table, filename, file_metadata);
  }

} // namespace OpenMS