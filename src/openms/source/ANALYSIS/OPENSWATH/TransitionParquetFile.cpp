// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionParquetFile.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/SYSTEM/File.h>

#include <fstream>
#ifdef WITH_PARQUET
#include <arrow/api.h>
#include <arrow/io/file.h>
#include <parquet/arrow/writer.h>
#include <parquet/arrow/reader.h>
#include <OpenMS/SYSTEM/ExternalProcess.h>
#include <QtCore/QStringList>
#endif

#include <unordered_map>
#include <memory>

namespace
{
#ifdef WITH_PARQUET
  void appendOrThrow_(const arrow::Status& status, const std::string& context)
  {
    if (!status.ok())
    {
      throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Failed to append value for " + context, status.ToString());
    }
  }

  template <typename BuilderT>
  std::shared_ptr<arrow::Array> finishArray_(BuilderT& builder, const std::string& context)
  {
    auto result = builder.Finish();
    if (!result.ok())
    {
      throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Failed to finish column for " + context, result.status().ToString());
    }
    return result.ValueOrDie();
  }

  ::arrow::Status writeParquetTable_(const std::shared_ptr<arrow::Table>& table, const OpenMS::String& filename)
  {
    auto outfile_result = arrow::io::FileOutputStream::Open(std::string(filename));
    if (!outfile_result.ok())
    {
      return outfile_result.status();
    }
    auto outfile = outfile_result.ValueOrDie();
    return parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), outfile, 1024);
  }

  std::string joinProteinAccessions_(const std::vector<std::string>& accessions)
  {
    OpenMS::String joined;
    for (OpenMS::Size i = 0; i < accessions.size(); ++i)
    {
      if (i > 0) joined += ";";
      joined += accessions[i];
    }
    return std::string(joined);
  }

  void writeLibraryMetadata_(const OpenMS::String& library_dir, const OpenMS::String& library_name)
  {
    const OpenMS::String metadata_path = library_dir + "/metadata.json";
    std::ofstream out(metadata_path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.is_open())
    {
      throw OpenMS::Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, metadata_path);
    }

    out << "{\n"
        << "  \"mzspec_lib\": {\n"
        << "    \"format_version\": \"1.0\",\n"
        << "    \"attributes\": [\n"
        << "      {\"accession\": \"MS:1003186\", \"name\": \"library format version\", \"value\": \"1.0\"},\n"
        << "      {\"accession\": \"MS:1003188\", \"name\": \"library name\", \"value\": \"" << library_name << "\"},\n"
        << "      {\"accession\": \"MS:1003207\", \"name\": \"library creation software\", \"value\": \"OpenMS\"}\n"
        << "    ]\n"
        << "  },\n"
        << "  \"openms\": {\n"
        << "    \"schema_version\": 1,\n"
        << "    \"generator\": \"OpenMS TransitionParquetFile\"\n"
        << "  }\n"
        << "}\n";
  }

  struct PrecursorInfo
  {
    double precursor_mz = 0.0;
    double library_rt = 0.0;
    double drift_time = -1.0;
    int charge = 0;
    bool decoy = false;
    std::string traml_id;
    std::string modified_sequence;
    std::string unmodified_sequence;
    std::vector<std::string> protein_accessions;
  };

  std::shared_ptr<arrow::Table> readParquetTable_(const OpenMS::String& filename)
  {
    auto infile_result = arrow::io::ReadableFile::Open(std::string(filename));
    if (!infile_result.ok())
    {
      throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Failed to open parquet file", filename);
    }
    std::shared_ptr<arrow::io::ReadableFile> infile = *infile_result;

    auto reader_result = parquet::arrow::OpenFile(infile, arrow::default_memory_pool());
    if (!reader_result.ok())
    {
      throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Failed to create parquet reader", filename);
    }
    std::unique_ptr<parquet::arrow::FileReader> reader = std::move(reader_result.ValueOrDie());

    std::shared_ptr<arrow::Table> table;
    auto read_status = reader->ReadTable(&table);
    if (!read_status.ok())
    {
      throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Failed to read parquet table", filename);
    }

    auto combined = table->CombineChunks(arrow::default_memory_pool());
    if (!combined.ok())
    {
      throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Failed to combine parquet chunks", filename);
    }

    return *combined;
  }

  std::shared_ptr<arrow::Array> getColumn_(const std::shared_ptr<arrow::Table>& table,
                                           const std::string& name)
  {
    auto column = table->GetColumnByName(name);
    if (!column)
    {
      throw OpenMS::Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                                  "Missing required column '" + name + "'");
    }
    if (column->num_chunks() == 0)
    {
      throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Column has no chunks", name);
    }
    return column->chunk(0);
  }

  std::shared_ptr<arrow::Array> getOptionalColumn_(const std::shared_ptr<arrow::Table>& table,
                                                   const std::string& name)
  {
    auto column = table->GetColumnByName(name);
    if (!column)
    {
      return nullptr;
    }
    if (column->num_chunks() == 0)
    {
      throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Optional column has no chunks", name);
    }
    return column->chunk(0);
  }

#ifdef WITH_PARQUET
  OpenMS::String prepareParquetDirectory_(const OpenMS::String& input_path,
                                          std::unique_ptr<OpenMS::File::TempDir>& temp_dir)
  {
    if (OpenMS::File::isDirectory(input_path))
    {
      return input_path;
    }

    if (!OpenMS::File::readable(input_path))
    {
      throw OpenMS::Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, input_path);
    }

    temp_dir = std::make_unique<OpenMS::File::TempDir>();
    const OpenMS::String unpack_dir = temp_dir->getPath() + "/pqp_parquet_unpacked";
    OpenMS::File::makeDir(unpack_dir);

    OpenMS::ExternalProcess unzip_process;
    QStringList args;
    args << "-q" << input_path.toQString() << "-d" << unpack_dir.toQString();
    auto status = unzip_process.run("unzip", args, "", false);
    if (status != OpenMS::ExternalProcess::RETURNSTATE::SUCCESS)
    {
      throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Failed to unzip parquet library", input_path);
    }

    return unpack_dir;
  }

  void zipParquetDirectory_(const OpenMS::String& directory_path, const OpenMS::String& output_zip)
  {
    const OpenMS::String output_zip_abs = OpenMS::File::absolutePath(output_zip);
    if (OpenMS::File::exists(output_zip_abs))
    {
      OpenMS::File::remove(output_zip_abs);
    }

    OpenMS::ExternalProcess zip_process;
    QStringList args;
    args << "-r" << "-q" << output_zip_abs.toQString() << ".";
    auto status = zip_process.run("zip", args, directory_path.toQString(), false);
    if (status != OpenMS::ExternalProcess::RETURNSTATE::SUCCESS)
    {
      throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Failed to zip parquet library", output_zip_abs);
    }
  }
#endif

  int64_t getInt64_(const std::shared_ptr<arrow::Array>& array, int64_t row, int64_t default_value, bool allow_null)
  {
    if (!array)
    {
      return default_value;
    }
    if (array->IsNull(row))
    {
      if (allow_null) return default_value;
      throw OpenMS::Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                                  "Required integer value is null");
    }
    switch (array->type_id())
    {
      case arrow::Type::INT64:
        return std::static_pointer_cast<arrow::Int64Array>(array)->Value(row);
      case arrow::Type::INT32:
        return static_cast<int64_t>(std::static_pointer_cast<arrow::Int32Array>(array)->Value(row));
      case arrow::Type::INT16:
        return static_cast<int64_t>(std::static_pointer_cast<arrow::Int16Array>(array)->Value(row));
      case arrow::Type::INT8:
        return static_cast<int64_t>(std::static_pointer_cast<arrow::Int8Array>(array)->Value(row));
      case arrow::Type::UINT64:
        return static_cast<int64_t>(std::static_pointer_cast<arrow::UInt64Array>(array)->Value(row));
      case arrow::Type::UINT32:
        return static_cast<int64_t>(std::static_pointer_cast<arrow::UInt32Array>(array)->Value(row));
      case arrow::Type::UINT16:
        return static_cast<int64_t>(std::static_pointer_cast<arrow::UInt16Array>(array)->Value(row));
      case arrow::Type::UINT8:
        return static_cast<int64_t>(std::static_pointer_cast<arrow::UInt8Array>(array)->Value(row));
      default:
        throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                              "Unsupported integer column type", "");
    }
    return default_value;
  }

  double getDouble_(const std::shared_ptr<arrow::Array>& array, int64_t row, double default_value, bool allow_null)
  {
    if (!array)
    {
      return default_value;
    }
    if (array->IsNull(row))
    {
      if (allow_null) return default_value;
      throw OpenMS::Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                                  "Required numeric value is null");
    }
    switch (array->type_id())
    {
      case arrow::Type::DOUBLE:
        return std::static_pointer_cast<arrow::DoubleArray>(array)->Value(row);
      case arrow::Type::FLOAT:
        return static_cast<double>(std::static_pointer_cast<arrow::FloatArray>(array)->Value(row));
      case arrow::Type::INT64:
      case arrow::Type::INT32:
      case arrow::Type::INT16:
      case arrow::Type::INT8:
      case arrow::Type::UINT64:
      case arrow::Type::UINT32:
      case arrow::Type::UINT16:
      case arrow::Type::UINT8:
        return static_cast<double>(getInt64_(array, row, 0, false));
      default:
        throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                              "Unsupported numeric column type", "");
    }
    return default_value;
  }

  bool getBool_(const std::shared_ptr<arrow::Array>& array, int64_t row, bool default_value, bool allow_null)
  {
    if (!array)
    {
      return default_value;
    }
    if (array->IsNull(row))
    {
      if (allow_null) return default_value;
      throw OpenMS::Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                                  "Required boolean value is null");
    }
    switch (array->type_id())
    {
      case arrow::Type::BOOL:
        return std::static_pointer_cast<arrow::BooleanArray>(array)->Value(row);
      case arrow::Type::INT8:
      case arrow::Type::INT16:
      case arrow::Type::INT32:
      case arrow::Type::INT64:
      case arrow::Type::UINT8:
      case arrow::Type::UINT16:
      case arrow::Type::UINT32:
      case arrow::Type::UINT64:
        return getInt64_(array, row, 0, false) != 0;
      default:
        throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                              "Unsupported boolean column type", "");
    }
    return default_value;
  }

  std::string getString_(const std::shared_ptr<arrow::Array>& array, int64_t row)
  {
    if (!array)
    {
      return "";
    }
    if (array->IsNull(row)) return "";
    switch (array->type_id())
    {
      case arrow::Type::STRING:
        return std::static_pointer_cast<arrow::StringArray>(array)->GetString(row);
      case arrow::Type::LARGE_STRING:
        return std::static_pointer_cast<arrow::LargeStringArray>(array)->GetString(row);
      default:
        throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                              "Unsupported string column type", "");
    }
    return "";
  }

  std::vector<std::string> getStringList_(const std::shared_ptr<arrow::Array>& array, int64_t row)
  {
    std::vector<std::string> values;
    if (!array) return values;
    if (array->IsNull(row)) return values;

    if (array->type_id() == arrow::Type::STRING || array->type_id() == arrow::Type::LARGE_STRING)
    {
      OpenMS::String raw = getString_(array, row);
      if (!raw.empty())
      {
        std::vector<OpenMS::String> parts;
        raw.split(';', parts);
        values.reserve(parts.size());
        for (const auto& part : parts)
        {
          values.push_back(part);
        }
      }
      return values;
    }

    if (array->type_id() == arrow::Type::LIST || array->type_id() == arrow::Type::LARGE_LIST)
    {
      if (array->type_id() == arrow::Type::LIST)
      {
        auto list_array = std::static_pointer_cast<arrow::ListArray>(array);
        auto values_array = list_array->values();
        auto start = list_array->value_offset(row);
        auto end = list_array->value_offset(row + 1);
        for (int64_t i = start; i < end; ++i)
        {
          values.push_back(getString_(values_array, i));
        }
        return values;
      }
      auto list_array = std::static_pointer_cast<arrow::LargeListArray>(array);
      auto values_array = list_array->values();
      auto start = list_array->value_offset(row);
      auto end = list_array->value_offset(row + 1);
      for (int64_t i = start; i < end; ++i)
      {
        values.push_back(getString_(values_array, i));
      }
      return values;
    }

    throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Unsupported list column type", "");
  }

  void addModifications_(const std::string& sequence, OpenSwath::LightCompound& compound)
  {
    if (sequence.empty()) return;
    try
    {
      OpenMS::AASequence aa_sequence = OpenMS::AASequence::fromString(sequence);
      if (aa_sequence.hasNTerminalModification())
      {
        OpenSwath::LightModification mod;
        mod.location = -1;
        mod.unimod_id = aa_sequence.getNTerminalModification()->getUniModRecordId();
        compound.modifications.push_back(mod);
      }
      if (aa_sequence.hasCTerminalModification())
      {
        OpenSwath::LightModification mod;
        mod.location = static_cast<int>(aa_sequence.size());
        mod.unimod_id = aa_sequence.getCTerminalModification()->getUniModRecordId();
        compound.modifications.push_back(mod);
      }
      for (OpenMS::Size i = 0; i != aa_sequence.size(); i++)
      {
        if (aa_sequence[i].isModified())
        {
          OpenSwath::LightModification mod;
          mod.location = static_cast<int>(i);
          mod.unimod_id = aa_sequence.getResidue(i).getModification()->getUniModRecordId();
          compound.modifications.push_back(mod);
        }
      }
    }
    catch (OpenMS::Exception::InvalidValue&)
    {
      OPENMS_LOG_DEBUG << "Could not parse modifications from sequence: " << sequence << '\n';
    }
  }
#endif
} // namespace

namespace OpenMS
{
  void TransitionParquetFile::convertParquetToTargetedExperiment(
    const String& pqp_parquet_dir, OpenSwath::LightTargetedExperiment& targeted_exp) const
  {
#ifdef WITH_PARQUET
    std::unique_ptr<File::TempDir> temp_dir;
    const String base_dir = prepareParquetDirectory_(pqp_parquet_dir, temp_dir);
    const String library_dir = base_dir + "/library";
    const String precursors_path = library_dir + "/precursors.parquet";
    const String transitions_path = library_dir + "/transitions.parquet";

    if (!File::exists(precursors_path) || !File::exists(transitions_path))
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Missing required parquet files in '" + pqp_parquet_dir + "'");
    }

    auto precursors_table = readParquetTable_(precursors_path);
    auto transitions_table = readParquetTable_(transitions_path);

    auto precursor_id_col = getColumn_(precursors_table, "precursor_id");
    auto precursor_mz_col = getColumn_(precursors_table, "precursor_mz");
    auto charge_col = getColumn_(precursors_table, "charge");
    auto library_rt_col = getColumn_(precursors_table, "library_rt");
    auto drift_time_col = getOptionalColumn_(precursors_table, "library_drift_time");
    auto traml_id_col = getOptionalColumn_(precursors_table, "traml_id");
    auto decoy_col = getOptionalColumn_(precursors_table, "decoy");
    auto modified_sequence_col = getOptionalColumn_(precursors_table, "modified_sequence");
    auto unmodified_sequence_col = getOptionalColumn_(precursors_table, "unmodified_sequence");
    auto protein_accessions_col = getOptionalColumn_(precursors_table, "protein_accessions");

    std::unordered_map<int64_t, PrecursorInfo> precursor_map;
    precursor_map.reserve(precursors_table->num_rows());

    for (int64_t row = 0; row < precursors_table->num_rows(); ++row)
    {
      const int64_t precursor_id = getInt64_(precursor_id_col, row, 0, false);
      PrecursorInfo info;
      info.precursor_mz = getDouble_(precursor_mz_col, row, 0.0, false);
      info.charge = static_cast<int>(getInt64_(charge_col, row, 0, false));
      info.library_rt = getDouble_(library_rt_col, row, 0.0, true);
      info.drift_time = getDouble_(drift_time_col, row, -1.0, true);
      info.decoy = getBool_(decoy_col, row, false, true);
      info.traml_id = getString_(traml_id_col, row);
      info.modified_sequence = getString_(modified_sequence_col, row);
      info.unmodified_sequence = getString_(unmodified_sequence_col, row);
      info.protein_accessions = getStringList_(protein_accessions_col, row);

      precursor_map.emplace(precursor_id, std::move(info));
    }

    std::unordered_map<std::string, int> compound_map;
    std::unordered_map<std::string, int> protein_map;
    compound_map.reserve(precursor_map.size());

    for (const auto& entry : precursor_map)
    {
      const int64_t precursor_id = entry.first;
      const PrecursorInfo& info = entry.second;
      const String precursor_id_str(precursor_id);

      OpenSwath::LightCompound compound;
      compound.id = precursor_id_str;
      compound.drift_time = info.drift_time;
      compound.rt = info.library_rt;
      compound.charge = info.charge;
      compound.sequence = info.modified_sequence.empty() ? info.unmodified_sequence : info.modified_sequence;
      compound.protein_refs = info.protein_accessions;
      if (!compound.sequence.empty())
      {
        addModifications_(compound.sequence, compound);
      }

      targeted_exp.compounds.push_back(std::move(compound));
      compound_map[precursor_id_str] = 0;

      for (const auto& accession : info.protein_accessions)
      {
        if (protein_map.find(accession) == protein_map.end())
        {
          OpenSwath::LightProtein protein;
          protein.id = accession;
          protein.sequence = "";
          targeted_exp.proteins.push_back(std::move(protein));
          protein_map[accession] = 0;
        }
      }
    }

    auto transition_id_col = getColumn_(transitions_table, "transition_id");
    auto transition_traml_id_col = getOptionalColumn_(transitions_table, "traml_id");
    auto transition_precursor_id_col = getColumn_(transitions_table, "precursor_id");
    auto product_mz_col = getColumn_(transitions_table, "product_mz");
    auto fragment_charge_col = getColumn_(transitions_table, "charge");
    auto fragment_type_col = getColumn_(transitions_table, "type");
    auto fragment_annotation_col = getOptionalColumn_(transitions_table, "annotation");
    auto fragment_ordinal_col = getColumn_(transitions_table, "ordinal");
    auto detecting_col = getColumn_(transitions_table, "detecting");
    auto identifying_col = getColumn_(transitions_table, "identifying");
    auto quantifying_col = getColumn_(transitions_table, "quantifying");
    auto transition_intensity_col = getColumn_(transitions_table, "library_intensity");
    auto transition_decoy_col = getColumn_(transitions_table, "decoy");

    for (int64_t row = 0; row < transitions_table->num_rows(); ++row)
    {
      const int64_t precursor_id = getInt64_(transition_precursor_id_col, row, 0, false);
      auto precursor_it = precursor_map.find(precursor_id);
      if (precursor_it == precursor_map.end())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Transition references unknown precursor_id");
      }

      const int64_t transition_id = getInt64_(transition_id_col, row, 0, false);
      const std::string traml_id = getString_(transition_traml_id_col, row);
      const std::string transition_name = traml_id.empty() ? String(transition_id) : String(traml_id);
      const std::string fragment_type = getString_(fragment_type_col, row);
      const std::string annotation = getString_(fragment_annotation_col, row);

      OpenSwath::LightTransition transition;
      transition.transition_name = transition_name;
      transition.peptide_ref = String(precursor_id);
      transition.library_intensity = getDouble_(transition_intensity_col, row, 0.0, true);
      transition.precursor_mz = precursor_it->second.precursor_mz;
      transition.product_mz = getDouble_(product_mz_col, row, 0.0, false);
      transition.precursor_im = precursor_it->second.drift_time;
      transition.fragment_charge = static_cast<int8_t>(getInt64_(fragment_charge_col, row, 0, true));
      transition.fragment_nr = static_cast<int16_t>(getInt64_(fragment_ordinal_col, row, -1, true));
      transition.setFragmentType(fragment_type.empty() ? annotation : fragment_type);
      transition.setDecoy(getBool_(transition_decoy_col, row, false, true));
      transition.setDetectingTransition(getBool_(detecting_col, row, true, true));
      transition.setIdentifyingTransition(getBool_(identifying_col, row, false, true));
      transition.setQuantifyingTransition(getBool_(quantifying_col, row, true, true));

      targeted_exp.transitions.push_back(std::move(transition));
    }
#else
    (void)pqp_parquet_dir;
    (void)targeted_exp;
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
#endif
  }

  void TransitionParquetFile::convertLightTargetedExperimentToParquet(
    const String& pqp_parquet_path, const OpenSwath::LightTargetedExperiment& targeted_exp) const
  {
#ifdef WITH_PARQUET
    const bool output_is_dir = File::isDirectory(pqp_parquet_path);
    std::unique_ptr<File::TempDir> temp_dir;
    String base_dir = pqp_parquet_path;
    if (!output_is_dir)
    {
      temp_dir = std::make_unique<File::TempDir>();
      base_dir = temp_dir->getPath() + "/pqp_parquet_output";
      File::makeDir(base_dir);
    }

    const String library_dir = base_dir + "/library";
    File::makeDir(library_dir);
    String library_name = File::basename(pqp_parquet_path);
    if (library_name.empty())
    {
      library_name = "openms_library";
    }
    writeLibraryMetadata_(library_dir, library_name);

    std::unordered_map<std::string, int64_t> compound_to_precursor;
    compound_to_precursor.reserve(targeted_exp.compounds.size());

    std::unordered_map<std::string, bool> compound_decoy;
    compound_decoy.reserve(targeted_exp.compounds.size());
    for (const auto& transition : targeted_exp.transitions)
    {
      if (transition.getDecoy())
      {
        compound_decoy[transition.peptide_ref] = true;
      }
    }

    int64_t next_precursor_id = 1;
    for (const auto& compound : targeted_exp.compounds)
    {
      if (compound_to_precursor.find(compound.id) != compound_to_precursor.end())
      {
        continue;
      }

      int64_t precursor_id = 0;
      try
      {
        precursor_id = OpenMS::String(compound.id).toInt64();
      }
      catch (OpenMS::Exception::ConversionError&)
      {
        precursor_id = next_precursor_id++;
      }

      if (precursor_id >= next_precursor_id)
      {
        next_precursor_id = precursor_id + 1;
      }

      compound_to_precursor[compound.id] = precursor_id;
    }

    std::unordered_map<std::string, double> precursor_mz;
    for (const auto& transition : targeted_exp.transitions)
    {
      if (precursor_mz.find(transition.peptide_ref) == precursor_mz.end())
      {
        precursor_mz[transition.peptide_ref] = transition.precursor_mz;
      }
    }

    arrow::Int64Builder precursor_id_builder;
    arrow::DoubleBuilder precursor_mz_builder;
    arrow::Int32Builder precursor_charge_builder;
    arrow::DoubleBuilder library_rt_builder;
    arrow::DoubleBuilder drift_time_builder;
    arrow::BooleanBuilder decoy_builder;
    arrow::StringBuilder traml_id_builder;
    arrow::StringBuilder modified_sequence_builder;
    arrow::StringBuilder unmodified_sequence_builder;
    arrow::StringBuilder protein_accessions_builder;

    for (const auto& compound : targeted_exp.compounds)
    {
      const int64_t precursor_id = compound_to_precursor[compound.id];
      const bool is_decoy = compound_decoy[compound.id] ||
        OpenMS::String(compound.id).hasPrefix("DECOY_");
      appendOrThrow_(precursor_id_builder.Append(precursor_id), "precursor_id");
      appendOrThrow_(precursor_mz_builder.Append(precursor_mz[compound.id]), "precursor_mz");
      appendOrThrow_(precursor_charge_builder.Append(compound.charge), "charge");
      appendOrThrow_(library_rt_builder.Append(compound.rt), "library_rt");
      appendOrThrow_(drift_time_builder.Append(compound.drift_time), "library_drift_time");
      appendOrThrow_(decoy_builder.Append(is_decoy), "decoy");
      appendOrThrow_(traml_id_builder.Append(compound.id), "traml_id");
      appendOrThrow_(modified_sequence_builder.Append(compound.sequence), "modified_sequence");

      std::string unmodified_sequence;
      if (!compound.sequence.empty())
      {
        try
        {
          unmodified_sequence = OpenMS::AASequence::fromString(compound.sequence).toUnmodifiedString();
        }
        catch (OpenMS::Exception::InvalidValue&)
        {
          unmodified_sequence = "";
        }
      }
      appendOrThrow_(unmodified_sequence_builder.Append(unmodified_sequence), "unmodified_sequence");
      appendOrThrow_(protein_accessions_builder.Append(joinProteinAccessions_(compound.protein_refs)), "protein_accessions");
    }

    auto precursor_schema = arrow::schema({
      arrow::field("precursor_id", arrow::int64()),
      arrow::field("precursor_mz", arrow::float64()),
      arrow::field("charge", arrow::int32()),
      arrow::field("library_rt", arrow::float64()),
      arrow::field("library_drift_time", arrow::float64()),
      arrow::field("decoy", arrow::boolean()),
      arrow::field("traml_id", arrow::utf8()),
      arrow::field("modified_sequence", arrow::utf8()),
      arrow::field("unmodified_sequence", arrow::utf8()),
      arrow::field("protein_accessions", arrow::utf8())
    });

    auto precursors_table = arrow::Table::Make(
      precursor_schema,
      {finishArray_(precursor_id_builder, "precursor_id"),
       finishArray_(precursor_mz_builder, "precursor_mz"),
       finishArray_(precursor_charge_builder, "charge"),
       finishArray_(library_rt_builder, "library_rt"),
       finishArray_(drift_time_builder, "library_drift_time"),
       finishArray_(decoy_builder, "decoy"),
       finishArray_(traml_id_builder, "traml_id"),
       finishArray_(modified_sequence_builder, "modified_sequence"),
       finishArray_(unmodified_sequence_builder, "unmodified_sequence"),
       finishArray_(protein_accessions_builder, "protein_accessions")});

    auto precursors_status = writeParquetTable_(precursors_table, library_dir + "/precursors.parquet");
    if (!precursors_status.ok())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Failed to write precursors parquet", precursors_status.ToString());
    }

    arrow::Int64Builder transition_id_builder;
    arrow::Int64Builder transition_precursor_id_builder;
    arrow::StringBuilder transition_traml_id_builder;
    arrow::DoubleBuilder product_mz_builder;
    arrow::Int32Builder fragment_charge_builder;
    arrow::StringBuilder fragment_type_builder;
    arrow::StringBuilder annotation_builder;
    arrow::Int32Builder ordinal_builder;
    arrow::BooleanBuilder detecting_builder;
    arrow::BooleanBuilder identifying_builder;
    arrow::BooleanBuilder quantifying_builder;
    arrow::DoubleBuilder transition_intensity_builder;
    arrow::BooleanBuilder transition_decoy_builder;

    int64_t transition_id = 1;
    for (const auto& transition : targeted_exp.transitions)
    {
      const int64_t precursor_ref = compound_to_precursor[transition.peptide_ref];
      appendOrThrow_(transition_id_builder.Append(transition_id++), "transition_id");
      appendOrThrow_(transition_precursor_id_builder.Append(precursor_ref), "precursor_id");
      appendOrThrow_(transition_traml_id_builder.Append(transition.transition_name), "traml_id");
      appendOrThrow_(product_mz_builder.Append(transition.product_mz), "product_mz");
      appendOrThrow_(fragment_charge_builder.Append(static_cast<int32_t>(transition.fragment_charge)), "charge");
      appendOrThrow_(fragment_type_builder.Append(transition.getFragmentType()), "type");
      appendOrThrow_(annotation_builder.Append(transition.getAnnotation()), "annotation");
      appendOrThrow_(ordinal_builder.Append(transition.fragment_nr), "ordinal");
      appendOrThrow_(detecting_builder.Append(transition.isDetectingTransition()), "detecting");
      appendOrThrow_(identifying_builder.Append(transition.isIdentifyingTransition()), "identifying");
      appendOrThrow_(quantifying_builder.Append(transition.isQuantifyingTransition()), "quantifying");
      appendOrThrow_(transition_intensity_builder.Append(transition.library_intensity), "library_intensity");
      appendOrThrow_(transition_decoy_builder.Append(transition.getDecoy()), "decoy");
    }

    auto transition_schema = arrow::schema({
      arrow::field("transition_id", arrow::int64()),
      arrow::field("precursor_id", arrow::int64()),
      arrow::field("traml_id", arrow::utf8()),
      arrow::field("product_mz", arrow::float64()),
      arrow::field("charge", arrow::int32()),
      arrow::field("type", arrow::utf8()),
      arrow::field("annotation", arrow::utf8()),
      arrow::field("ordinal", arrow::int32()),
      arrow::field("detecting", arrow::boolean()),
      arrow::field("identifying", arrow::boolean()),
      arrow::field("quantifying", arrow::boolean()),
      arrow::field("library_intensity", arrow::float64()),
      arrow::field("decoy", arrow::boolean())
    });

    auto transitions_table = arrow::Table::Make(
      transition_schema,
      {finishArray_(transition_id_builder, "transition_id"),
       finishArray_(transition_precursor_id_builder, "precursor_id"),
       finishArray_(transition_traml_id_builder, "traml_id"),
       finishArray_(product_mz_builder, "product_mz"),
       finishArray_(fragment_charge_builder, "charge"),
       finishArray_(fragment_type_builder, "type"),
       finishArray_(annotation_builder, "annotation"),
       finishArray_(ordinal_builder, "ordinal"),
       finishArray_(detecting_builder, "detecting"),
       finishArray_(identifying_builder, "identifying"),
       finishArray_(quantifying_builder, "quantifying"),
       finishArray_(transition_intensity_builder, "library_intensity"),
       finishArray_(transition_decoy_builder, "decoy")});

    auto transitions_status = writeParquetTable_(transitions_table, library_dir + "/transitions.parquet");
    if (!transitions_status.ok())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Failed to write transitions parquet", transitions_status.ToString());
    }

    if (!output_is_dir)
    {
      zipParquetDirectory_(base_dir, pqp_parquet_path);
    }
#else
    (void)pqp_parquet_path;
    (void)targeted_exp;
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
#endif
  }
} // namespace OpenMS
