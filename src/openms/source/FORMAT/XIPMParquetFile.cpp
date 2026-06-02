// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/XIPMParquetFile.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/MSNumpressCoder.h>
#include <OpenMS/FORMAT/ZlibCompression.h>
#include <OpenMS/SYSTEM/File.h>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/reader.h>

#include <cstring>
#include <limits>
#include <memory>
#include <string>
#include <unordered_set>

namespace OpenMS
{
  namespace
  {
    std::shared_ptr<arrow::Table> readParquetTable_(const String& filename)
    {
      auto infile_result = arrow::io::ReadableFile::Open(std::string(filename));
      if (!infile_result.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to open parquet file", filename);
      }
      const std::shared_ptr<arrow::io::ReadableFile>& infile = *infile_result;

      auto reader_result = parquet::arrow::OpenFile(infile, arrow::default_memory_pool());
      if (!reader_result.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to create parquet reader", filename);
      }
      std::unique_ptr<parquet::arrow::FileReader> reader = std::move(*reader_result);

      auto table_result = reader->ReadTable();
      if (!table_result.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to read parquet table: " + table_result.status().ToString(), filename);
      }
      std::shared_ptr<arrow::Table> table = *table_result;

      auto combined = table->CombineChunks(arrow::default_memory_pool());
      if (!combined.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to combine parquet chunks", filename);
      }

      return *combined;
    }

    std::shared_ptr<arrow::Table> readParquetTable_(const std::vector<String>& filenames)
    {
      if (filenames.empty())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "No parquet files provided", "");
      }
      if (filenames.size() == 1)
      {
        return readParquetTable_(filenames.front());
      }

      std::vector<std::shared_ptr<arrow::Table>> tables;
      tables.reserve(filenames.size());
      for (const auto& filename : filenames)
      {
        tables.push_back(readParquetTable_(filename));
      }

      auto concat_result = arrow::ConcatenateTables(tables);
      if (!concat_result.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to concatenate parquet tables", concat_result.status().ToString());
      }

      auto combined = (*concat_result)->CombineChunks(arrow::default_memory_pool());
      if (!combined.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to combine parquet chunks", combined.status().ToString());
      }
      return *combined;
    }

    std::shared_ptr<arrow::Schema> readParquetSchema_(const String& filename)
    {
      auto infile_result = arrow::io::ReadableFile::Open(std::string(filename));
      if (!infile_result.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to open parquet file", filename);
      }
      const std::shared_ptr<arrow::io::ReadableFile>& infile = *infile_result;

      auto reader_result = parquet::arrow::OpenFile(infile, arrow::default_memory_pool());
      if (!reader_result.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to create parquet reader", filename);
      }
      std::unique_ptr<parquet::arrow::FileReader> reader = std::move(*reader_result);

      std::shared_ptr<arrow::Schema> schema;
      auto status = reader->GetSchema(&schema);
      if (!status.ok() || !schema)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to read parquet schema", filename);
      }
      return schema;
    }

    std::shared_ptr<arrow::Schema> readParquetSchemaAllFiles_(const std::vector<String>& filenames)
    {
      if (filenames.empty())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "No parquet files provided for schema validation", "");
      }
      auto base = readParquetSchema_(filenames.front());
      for (Size i = 1; i < filenames.size(); ++i)
      {
        auto other = readParquetSchema_(filenames[i]);
        if (!other->Equals(*base))
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Parquet input files have incompatible schemas: '" + filenames.front() +
                                        "' vs '" + filenames[i] + "'", "");
        }
      }
      return base;
    }

    std::shared_ptr<arrow::Array> getColumn_(const std::shared_ptr<arrow::Table>& table,
                                             const std::string& name,
                                             bool required = true)
    {
      auto column = table->GetColumnByName(name);
      if (!column)
      {
        if (required)
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                              "Missing required column '" + name + "'");
        }
        return nullptr;
      }
      if (column->num_chunks() == 0)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Column has no chunks", name);
      }
      if (column->num_chunks() == 1)
      {
        return column->chunk(0);
      }
      auto combined = arrow::Concatenate(column->chunks(), arrow::default_memory_pool());
      if (!combined.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to combine column chunks", combined.status().ToString());
      }
      return *combined;
    }

    bool getOptionalInt_(const std::shared_ptr<arrow::Array>& array, int64_t row, Int64& value)
    {
      if (!array || array->IsNull(row))
      {
        return false;
      }
      auto typed = std::static_pointer_cast<arrow::Int64Array>(array);
      value = typed->Value(row);
      return true;
    }

    bool getOptionalDouble_(const std::shared_ptr<arrow::Array>& array, int64_t row, double& value)
    {
      if (!array || array->IsNull(row))
      {
        return false;
      }
      auto typed = std::static_pointer_cast<arrow::DoubleArray>(array);
      value = typed->Value(row);
      return true;
    }

    bool getOptionalString_(const std::shared_ptr<arrow::Array>& array, int64_t row, String& value)
    {
      if (!array || array->IsNull(row))
      {
        return false;
      }
      auto typed = std::static_pointer_cast<arrow::StringArray>(array);
      value = typed->GetString(row);
      return true;
    }

    String getBinaryView_(const std::shared_ptr<arrow::Array>& array, int64_t row)
    {
      auto typed = std::static_pointer_cast<arrow::BinaryArray>(array);
      const auto view = typed->GetView(row);
      return String(view.data(), static_cast<Int>(view.size()));
    }

    void decodeBinary_(const String& data, Int64 compression, std::vector<double>& output)
    {
      output.clear();
      if (data.empty())
      {
        return;
      }

      auto decodeNumpress = [&](MSNumpressCoder::NumpressCompression type, const String& input)
      {
        MSNumpressCoder::NumpressConfig config;
        config.np_compression = type;
        MSNumpressCoder().decodeNPRaw(input, output, config);
      };

      switch (compression)
      {
        case 0:
        {
          if (data.size() % sizeof(double) != 0)
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Invalid binary data size (not divisible by 8)", String(data.size()));
          }
          const size_t count = data.size() / sizeof(double);
          output.resize(count);
          std::memcpy(output.data(), data.c_str(), count * sizeof(double));
          break;
        }
        case 1:
        {
          String decoded;
          ZlibCompression::uncompressString(data, decoded);
          if (decoded.size() % sizeof(double) != 0)
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Invalid decompressed binary data size", String(decoded.size()));
          }
          const size_t count = decoded.size() / sizeof(double);
          output.resize(count);
          std::memcpy(output.data(), decoded.c_str(), count * sizeof(double));
          break;
        }
        case 2:
          decodeNumpress(MSNumpressCoder::LINEAR, data);
          break;
        case 3:
          decodeNumpress(MSNumpressCoder::SLOF, data);
          break;
        case 4:
          decodeNumpress(MSNumpressCoder::PIC, data);
          break;
        case 5:
        {
          String decoded;
          ZlibCompression::uncompressString(data, decoded);
          decodeNumpress(MSNumpressCoder::LINEAR, decoded);
          break;
        }
        case 6:
        {
          String decoded;
          ZlibCompression::uncompressString(data, decoded);
          decodeNumpress(MSNumpressCoder::SLOF, decoded);
          break;
        }
        default:
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Unknown compression code", String(compression));
      }
    }
  } // namespace

  XIPMParquetFile::XIPMParquetFile(const String& filename) :
    filename_(filename),
    filenames_({filename})
  {
    if (!File::exists(filename))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
  }

  XIPMParquetFile::XIPMParquetFile(const std::vector<String>& filenames) :
    filename_(filenames.empty() ? "" : filenames.front()),
    filenames_(filenames)
  {
    if (filenames.empty())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No .xipm files given");
    }
    for (const auto& filename : filenames)
    {
      if (!File::exists(filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
    }
  }

  const String& XIPMParquetFile::getFilename() const
  {
    return filename_;
  }

  const std::vector<String>& XIPMParquetFile::getFilenames() const
  {
    return filenames_;
  }

  void XIPMParquetFile::load(std::vector<XIPMPeakMap>& output) const
  {
    getPeakMaps(output);
  }

  void XIPMParquetFile::getPeakMaps(std::vector<XIPMPeakMap>& output,
                                    Int64 precursor_id,
                                    Int64 transition_id,
                                    const String& modified_sequence,
                                    Int64 precursor_charge,
                                    Int64 product_charge,
                                    Int64 ms_level,
                                    Int64 run_id,
                                    const String& peakmap_type) const
  {
    output.clear();
    auto table = readParquetTable_(filenames_);

    auto run_id_col = getColumn_(table, "RUN_ID");
    auto source_file_col = getColumn_(table, "SOURCE_FILE", false);
    auto ms_level_col = getColumn_(table, "MS_LEVEL");
    auto peakmap_type_col = getColumn_(table, "PEAKMAP_TYPE", false);
    auto precursor_id_col = getColumn_(table, "PRECURSOR_ID", false);
    auto transition_id_col = getColumn_(table, "TRANSITION_ID", false);
    auto modified_sequence_col = getColumn_(table, "MODIFIED_SEQUENCE", false);
    auto precursor_charge_col = getColumn_(table, "PRECURSOR_CHARGE", false);
    auto product_charge_col = getColumn_(table, "PRODUCT_CHARGE", false);
    auto detecting_transition_col = getColumn_(table, "DETECTING_TRANSITION", false);
    auto precursor_decoy_col = getColumn_(table, "PRECURSOR_DECOY", false);
    auto product_decoy_col = getColumn_(table, "PRODUCT_DECOY", false);
    auto transition_ordinal_col = getColumn_(table, "TRANSITION_ORDINAL", false);
    auto transition_type_col = getColumn_(table, "TRANSITION_TYPE", false);
    auto annotation_col = getColumn_(table, "ANNOTATION", false);
    auto target_mz_col = getColumn_(table, "TARGET_MZ");
    auto target_rt_col = getColumn_(table, "TARGET_RT", false);
    auto target_im_col = getColumn_(table, "TARGET_ION_MOBILITY", false);
    auto rt_start_col = getColumn_(table, "RT_START", false);
    auto rt_end_col = getColumn_(table, "RT_END", false);
    auto mz_data_col = getColumn_(table, "MZ_DATA");
    auto rt_data_col = getColumn_(table, "RT_DATA");
    auto mobility_data_col = getColumn_(table, "MOBILITY_DATA");
    auto intensity_data_col = getColumn_(table, "INTENSITY_DATA");
    auto mz_compression_col = getColumn_(table, "MZ_COMPRESSION");
    auto rt_compression_col = getColumn_(table, "RT_COMPRESSION");
    auto mobility_compression_col = getColumn_(table, "MOBILITY_COMPRESSION");
    auto intensity_compression_col = getColumn_(table, "INTENSITY_COMPRESSION");

    const int64_t rows = table->num_rows();
    output.reserve(rows);

    for (int64_t row = 0; row < rows; ++row)
    {
      XIPMPeakMap peak_map;
      getOptionalInt_(run_id_col, row, peak_map.run_id);
      getOptionalString_(source_file_col, row, peak_map.source_file);
      getOptionalString_(peakmap_type_col, row, peak_map.peakmap_type);
      peak_map.has_precursor_id = getOptionalInt_(precursor_id_col, row, peak_map.precursor_id);
      peak_map.has_transition_id = getOptionalInt_(transition_id_col, row, peak_map.transition_id);
      getOptionalString_(modified_sequence_col, row, peak_map.modified_sequence);
      peak_map.has_precursor_charge = getOptionalInt_(precursor_charge_col, row, peak_map.precursor_charge);
      peak_map.has_product_charge = getOptionalInt_(product_charge_col, row, peak_map.product_charge);
      peak_map.has_detecting_transition = getOptionalInt_(detecting_transition_col, row, peak_map.detecting_transition);
      peak_map.has_precursor_decoy = getOptionalInt_(precursor_decoy_col, row, peak_map.precursor_decoy);
      peak_map.has_product_decoy = getOptionalInt_(product_decoy_col, row, peak_map.product_decoy);
      peak_map.has_transition_ordinal = getOptionalInt_(transition_ordinal_col, row, peak_map.transition_ordinal);
      getOptionalString_(transition_type_col, row, peak_map.transition_type);
      getOptionalString_(annotation_col, row, peak_map.annotation);
      getOptionalInt_(ms_level_col, row, peak_map.ms_level);
      getOptionalDouble_(target_mz_col, row, peak_map.target_mz);
      peak_map.has_target_rt = getOptionalDouble_(target_rt_col, row, peak_map.target_rt);
      peak_map.has_target_ion_mobility = getOptionalDouble_(target_im_col, row, peak_map.target_ion_mobility);
      peak_map.has_rt_start = getOptionalDouble_(rt_start_col, row, peak_map.rt_start);
      peak_map.has_rt_end = getOptionalDouble_(rt_end_col, row, peak_map.rt_end);
      if (!peak_map.has_target_rt && peak_map.has_rt_start && peak_map.has_rt_end && peak_map.rt_end > peak_map.rt_start)
      {
        peak_map.target_rt = (peak_map.rt_start + peak_map.rt_end) / 2.0;
        peak_map.has_target_rt = true;
      }

      if (precursor_id >= 0 && (!peak_map.has_precursor_id || peak_map.precursor_id != precursor_id)) continue;
      if (transition_id >= 0 && (!peak_map.has_transition_id || peak_map.transition_id != transition_id)) continue;
      if (!modified_sequence.empty() && peak_map.modified_sequence != modified_sequence) continue;
      if (precursor_charge >= 0 && (!peak_map.has_precursor_charge || peak_map.precursor_charge != precursor_charge)) continue;
      if (product_charge >= 0 && (!peak_map.has_product_charge || peak_map.product_charge != product_charge)) continue;
      if (ms_level >= 0 && peak_map.ms_level != ms_level) continue;
      if (run_id >= 0 && peak_map.run_id != run_id) continue;
      if (!peakmap_type.empty() && peak_map.peakmap_type != peakmap_type) continue;

      Int64 mz_compression = 0;
      Int64 rt_compression = 0;
      Int64 mobility_compression = 0;
      Int64 intensity_compression = 0;
      getOptionalInt_(mz_compression_col, row, mz_compression);
      getOptionalInt_(rt_compression_col, row, rt_compression);
      getOptionalInt_(mobility_compression_col, row, mobility_compression);
      getOptionalInt_(intensity_compression_col, row, intensity_compression);

      decodeBinary_(getBinaryView_(mz_data_col, row), mz_compression, peak_map.mz);
      decodeBinary_(getBinaryView_(rt_data_col, row), rt_compression, peak_map.rt);
      decodeBinary_(getBinaryView_(mobility_data_col, row), mobility_compression, peak_map.ion_mobility);
      decodeBinary_(getBinaryView_(intensity_data_col, row), intensity_compression, peak_map.intensity);

      if (peak_map.mz.size() != peak_map.rt.size() ||
          peak_map.mz.size() != peak_map.ion_mobility.size() ||
          peak_map.mz.size() != peak_map.intensity.size())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Peak-map array size mismatch", String(row));
      }

      output.push_back(std::move(peak_map));
    }
  }

  void XIPMParquetFile::getRuns(std::vector<XIPMRunInfo>& output) const
  {
    output.clear();
    auto table = readParquetTable_(filenames_);

    auto run_id_col = getColumn_(table, "RUN_ID");
    auto source_file_col = getColumn_(table, "SOURCE_FILE", false);

    std::unordered_set<std::string> seen;
    const int64_t rows = table->num_rows();
    output.reserve(rows);

    for (int64_t row = 0; row < rows; ++row)
    {
      XIPMRunInfo info;
      getOptionalInt_(run_id_col, row, info.run_id);
      getOptionalString_(source_file_col, row, info.source_file);

      const std::string key = std::to_string(info.run_id) + '\t' + std::string(info.source_file);
      if (seen.insert(key).second)
      {
        output.push_back(std::move(info));
      }
    }
  }

  void XIPMParquetFile::getColumns(std::vector<String>& output) const
  {
    output.clear();
    std::shared_ptr<arrow::Schema> schema = readParquetSchemaAllFiles_(filenames_);
    output.reserve(schema->num_fields());
    for (const auto& field : schema->fields())
    {
      output.emplace_back(field->name());
    }
  }
} // namespace OpenMS
