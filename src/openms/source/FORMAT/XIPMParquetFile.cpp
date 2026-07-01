// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/XIPMParquetFile.h>

#include <OpenMS/CONCEPT/Exception.h>
#include "DATAACCESS/ParquetReaderHelpers.h"
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/FORMAT/MSNumpressCoder.h>
#include <OpenMS/FORMAT/ZlibCompression.h>
#include <OpenMS/SYSTEM/File.h>

#include <arrow/api.h>

#include <cstring>
#include <memory>
#include <string>
#include <unordered_set>

namespace OpenMS
{
  namespace
  {
    namespace PRH = Internal::ParquetReaderHelpers;

    void validateXIPMTable_(const std::shared_ptr<arrow::Table>& table)
    {
      auto validation = ArrowSchemaValidation::validate(
        table, XIPMSchema::schema(), ArrowSchemaValidation::Mode::Subset);
      if (!validation.valid)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Invalid XIPM parquet schema", validation.toString());
      }
    }

    std::shared_ptr<arrow::Table> readParquetTable_(const std::string& filename)
    {
      auto table = PRH::readParquetTable(filename);
      validateXIPMTable_(table);
      return table;
    }

    std::shared_ptr<arrow::Table> readParquetTable_(const std::vector<std::string>& filenames)
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

      auto table = PRH::readParquetTable(filenames);
      validateXIPMTable_(table);
      return table;
    }

    std::shared_ptr<arrow::Table> readParquetTableColumns_(const std::vector<std::string>& filenames,
                                                           const std::vector<std::string>& columns)
    {
      if (filenames.empty())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "No parquet files provided", "");
      }
      auto table = PRH::readParquetTableColumns(filenames, columns);
      validateXIPMTable_(table);
      return table;
    }

    void decodeBinary_(const std::string& data, Int64 compression, std::vector<double>& output)
    {
      output.clear();
      if (data.empty())
      {
        return;
      }

      auto decodeNumpress = [&](MSNumpressCoder::NumpressCompression type, const std::string& input)
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
                                          "Invalid binary data size (not divisible by 8)", StringUtils::toStr(data.size()));
          }
          const size_t count = data.size() / sizeof(double);
          output.resize(count);
          std::memcpy(output.data(), data.data(), count * sizeof(double));
          break;
        }
        case 1:
        {
          std::string decoded;
          ZlibCompression::uncompressString(data, decoded);
          if (decoded.size() % sizeof(double) != 0)
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Invalid decompressed binary data size", StringUtils::toStr(decoded.size()));
          }
          const size_t count = decoded.size() / sizeof(double);
          output.resize(count);
          std::memcpy(output.data(), decoded.data(), count * sizeof(double));
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
          std::string decoded;
          ZlibCompression::uncompressString(data, decoded);
          decodeNumpress(MSNumpressCoder::LINEAR, decoded);
          break;
        }
        case 6:
        {
          std::string decoded;
          ZlibCompression::uncompressString(data, decoded);
          decodeNumpress(MSNumpressCoder::SLOF, decoded);
          break;
        }
        default:
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Unknown compression code", StringUtils::toStr(compression));
      }
    }
  } // namespace

  XIPMParquetFile::XIPMParquetFile(const std::string& filename) :
    filename_(filename),
    filenames_({filename})
  {
    if (!File::exists(filename))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
  }

  XIPMParquetFile::XIPMParquetFile(const std::vector<std::string>& filenames) :
    filename_(filenames.empty() ? std::string() : filenames.front()),
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

  const std::string& XIPMParquetFile::getFilename() const
  {
    return filename_;
  }

  const std::vector<std::string>& XIPMParquetFile::getFilenames() const
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
                                    const std::string& modified_sequence,
                                    Int64 precursor_charge,
                                    Int64 product_charge,
                                    Int64 ms_level,
                                    Int64 run_id,
                                    const std::string& peakmap_type) const
  {
    output.clear();
    auto table = readParquetTable_(filenames_);
    const int64_t rows = table->num_rows();
    if (rows == 0)
    {
      return;
    }

    auto run_id_col = PRH::getColumn(table, "RUN_ID");
    auto source_file_col = PRH::getColumn(table, "SOURCE_FILE", false);
    auto ms_level_col = PRH::getColumn(table, "MS_LEVEL");
    auto peakmap_type_col = PRH::getColumn(table, "PEAKMAP_TYPE", false);
    auto precursor_id_col = PRH::getColumn(table, "PRECURSOR_ID", false);
    auto transition_id_col = PRH::getColumn(table, "TRANSITION_ID", false);
    auto modified_sequence_col = PRH::getColumn(table, "MODIFIED_SEQUENCE", false);
    auto precursor_charge_col = PRH::getColumn(table, "PRECURSOR_CHARGE", false);
    auto product_charge_col = PRH::getColumn(table, "PRODUCT_CHARGE", false);
    auto detecting_transition_col = PRH::getColumn(table, "DETECTING_TRANSITION", false);
    auto precursor_decoy_col = PRH::getColumn(table, "PRECURSOR_DECOY", false);
    auto product_decoy_col = PRH::getColumn(table, "PRODUCT_DECOY", false);
    auto transition_ordinal_col = PRH::getColumn(table, "TRANSITION_ORDINAL", false);
    auto transition_type_col = PRH::getColumn(table, "TRANSITION_TYPE", false);
    auto annotation_col = PRH::getColumn(table, "ANNOTATION", false);
    auto target_mz_col = PRH::getColumn(table, "TARGET_MZ");
    auto target_rt_col = PRH::getColumn(table, "TARGET_RT", false);
    auto target_im_col = PRH::getColumn(table, "TARGET_ION_MOBILITY", false);
    auto rt_start_col = PRH::getColumn(table, "RT_START", false);
    auto rt_end_col = PRH::getColumn(table, "RT_END", false);
    auto mz_data_col = PRH::getColumn(table, "MZ_DATA");
    auto rt_data_col = PRH::getColumn(table, "RT_DATA");
    auto mobility_data_col = PRH::getColumn(table, "MOBILITY_DATA");
    auto intensity_data_col = PRH::getColumn(table, "INTENSITY_DATA");
    auto mz_compression_col = PRH::getColumn(table, "MZ_COMPRESSION");
    auto rt_compression_col = PRH::getColumn(table, "RT_COMPRESSION");
    auto mobility_compression_col = PRH::getColumn(table, "MOBILITY_COMPRESSION");
    auto intensity_compression_col = PRH::getColumn(table, "INTENSITY_COMPRESSION");

    output.reserve(rows);

    for (int64_t row = 0; row < rows; ++row)
    {
      XIPMPeakMap peak_map;
      PRH::getOptionalInt(run_id_col, row, peak_map.run_id);
      PRH::getOptionalString(source_file_col, row, peak_map.source_file);
      PRH::getOptionalString(peakmap_type_col, row, peak_map.peakmap_type);
      peak_map.has_precursor_id = PRH::getOptionalInt(precursor_id_col, row, peak_map.precursor_id);
      peak_map.has_transition_id = PRH::getOptionalInt(transition_id_col, row, peak_map.transition_id);
      PRH::getOptionalString(modified_sequence_col, row, peak_map.modified_sequence);
      peak_map.has_precursor_charge = PRH::getOptionalInt(precursor_charge_col, row, peak_map.precursor_charge);
      peak_map.has_product_charge = PRH::getOptionalInt(product_charge_col, row, peak_map.product_charge);
      peak_map.has_detecting_transition = PRH::getOptionalInt(detecting_transition_col, row, peak_map.detecting_transition);
      peak_map.has_precursor_decoy = PRH::getOptionalInt(precursor_decoy_col, row, peak_map.precursor_decoy);
      peak_map.has_product_decoy = PRH::getOptionalInt(product_decoy_col, row, peak_map.product_decoy);
      peak_map.has_transition_ordinal = PRH::getOptionalInt(transition_ordinal_col, row, peak_map.transition_ordinal);
      PRH::getOptionalString(transition_type_col, row, peak_map.transition_type);
      PRH::getOptionalString(annotation_col, row, peak_map.annotation);
      PRH::getOptionalInt(ms_level_col, row, peak_map.ms_level);
      PRH::getOptionalDouble(target_mz_col, row, peak_map.target_mz);
      peak_map.has_target_rt = PRH::getOptionalDouble(target_rt_col, row, peak_map.target_rt);
      peak_map.has_target_ion_mobility = PRH::getOptionalDouble(target_im_col, row, peak_map.target_ion_mobility);
      peak_map.has_rt_start = PRH::getOptionalDouble(rt_start_col, row, peak_map.rt_start);
      peak_map.has_rt_end = PRH::getOptionalDouble(rt_end_col, row, peak_map.rt_end);
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
      PRH::getOptionalInt(mz_compression_col, row, mz_compression);
      PRH::getOptionalInt(rt_compression_col, row, rt_compression);
      PRH::getOptionalInt(mobility_compression_col, row, mobility_compression);
      PRH::getOptionalInt(intensity_compression_col, row, intensity_compression);

      decodeBinary_(PRH::getBinaryView(mz_data_col, row), mz_compression, peak_map.mz);
      decodeBinary_(PRH::getBinaryView(rt_data_col, row), rt_compression, peak_map.rt);
      decodeBinary_(PRH::getBinaryView(mobility_data_col, row), mobility_compression, peak_map.ion_mobility);
      decodeBinary_(PRH::getBinaryView(intensity_data_col, row), intensity_compression, peak_map.intensity);

      if (peak_map.mz.size() != peak_map.rt.size() ||
          peak_map.mz.size() != peak_map.ion_mobility.size() ||
          peak_map.mz.size() != peak_map.intensity.size())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Peak-map array size mismatch", StringUtils::toStr(row));
      }

      output.push_back(std::move(peak_map));
    }
  }

  void XIPMParquetFile::getRuns(std::vector<XIPMRunInfo>& output) const
  {
    output.clear();
    auto table = readParquetTableColumns_(filenames_, {"RUN_ID", "SOURCE_FILE"});
    const int64_t rows = table->num_rows();
    if (rows == 0)
    {
      return;
    }

    auto run_id_col = PRH::getColumn(table, "RUN_ID");
    auto source_file_col = PRH::getColumn(table, "SOURCE_FILE", false);

    std::unordered_set<std::string> seen;
    output.reserve(rows);

    for (int64_t row = 0; row < rows; ++row)
    {
      XIPMRunInfo info;
      PRH::getOptionalInt(run_id_col, row, info.run_id);
      PRH::getOptionalString(source_file_col, row, info.source_file);

      const std::string key = std::to_string(info.run_id) + '\t' + std::string(info.source_file);
      if (seen.insert(key).second)
      {
        output.push_back(std::move(info));
      }
    }
  }

  void XIPMParquetFile::getColumns(std::vector<std::string>& output) const
  {
    output.clear();
    std::shared_ptr<arrow::Schema> schema = PRH::readParquetSchemaAllFiles(filenames_);
    output.reserve(schema->num_fields());
    for (const auto& field : schema->fields())
    {
      output.emplace_back(field->name());
    }
  }
} // namespace OpenMS
