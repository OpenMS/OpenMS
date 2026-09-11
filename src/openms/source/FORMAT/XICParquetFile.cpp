// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include "DATAACCESS/ParquetReaderHelpers.h"
#include <OpenMS/FORMAT/MSNumpressCoder.h>
#include <OpenMS/FORMAT/XICParquetFile.h>
#include <OpenMS/FORMAT/ZlibCompression.h>
#include <OpenMS/SYSTEM/File.h>


#include <cstring>
#include <cctype>
#include <exception>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

namespace OpenMS
{
  namespace
  {
    namespace PRH = Internal::ParquetReaderHelpers;

    const PRH::FilterPruningOptions& filterPruningOptions_()
    {
      static const PRH::FilterPruningOptions options{
        {"RT", "INTENSITY"},
        false,
        false
      };
      return options;
    }

    /// Decode RT/intensity arrays stored as compressed binary blobs.
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
        case 0: // no compression
        {
          if (data.size() % sizeof(double) != 0)
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Invalid RT/intensity data size (not divisible by 8)", StringUtils::toStr(data.size()));
          }
          const size_t count = data.size() / sizeof(double);
          output.resize(count);
          std::memcpy(output.data(), data.c_str(), count * sizeof(double));
          break;
        }
        case 1: // zlib only
        {
          std::string decoded;
          ZlibCompression::uncompressString(data, decoded);
          if (decoded.size() % sizeof(double) != 0)
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Invalid decompressed RT/intensity data size", StringUtils::toStr(decoded.size()));
          }
          const size_t count = decoded.size() / sizeof(double);
          output.resize(count);
          std::memcpy(output.data(), decoded.c_str(), count * sizeof(double));
          break;
        }
        case 2: // np-linear
        {
          decodeNumpress(MSNumpressCoder::LINEAR, data);
          break;
        }
        case 3: // np-slof
        {
          decodeNumpress(MSNumpressCoder::SLOF, data);
          break;
        }
        case 4: // np-pic
        {
          decodeNumpress(MSNumpressCoder::PIC, data);
          break;
        }
        case 5: // np-linear + zlib
        {
          std::string decoded;
          ZlibCompression::uncompressString(data, decoded);
          decodeNumpress(MSNumpressCoder::LINEAR, decoded);
          break;
        }
        case 6: // np-slof + zlib
        {
          std::string decoded;
          ZlibCompression::uncompressString(data, decoded);
          decodeNumpress(MSNumpressCoder::SLOF, decoded);
          break;
        }
        case 7: // np-pic + zlib
        {
          std::string decoded;
          ZlibCompression::uncompressString(data, decoded);
          decodeNumpress(MSNumpressCoder::PIC, decoded);
          break;
        }
        default:
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Unsupported compression type", StringUtils::toStr(compression));
      }
    }

    std::string upper_(const std::string& input)
    {
      std::string out = input;
      for (Size i = 0; i < out.size(); ++i)
      {
        out[i] = static_cast<char>(std::toupper(static_cast<unsigned char>(out[i])));
      }
      return out;
    }

    FilterExpression buildFilterFromArgs_(Int64 precursor_id,
                                          Int64 transition_id,
                                          const std::string& modified_sequence,
                                          Int64 precursor_charge,
                                          Int64 product_charge,
                                          Int64 ms_level,
                                          Int64 run_id)
    {
      ParquetFilter builder;
      if (precursor_id >= 0) builder.eq("PRECURSOR_ID", precursor_id);
      if (transition_id >= 0) builder.eq("TRANSITION_ID", transition_id);
      if (!modified_sequence.empty()) builder.eq("MODIFIED_SEQUENCE", modified_sequence);
      if (precursor_charge >= 0) builder.eq("PRECURSOR_CHARGE", precursor_charge);
      if (product_charge >= 0) builder.eq("PRODUCT_CHARGE", product_charge);
      if (ms_level >= 0) builder.eq("MS_LEVEL", ms_level);
      if (run_id >= 0) builder.eq("RUN_ID", run_id);
      return builder.expression();
    }

    void appendKeyPart_(std::ostringstream& oss, bool has_value, Int64 value)
    {
      if (has_value)
      {
        oss << value;
      }
      else
      {
        oss << "NULL";
      }
      oss << '|';
    }

    void appendKeyPart_(std::ostringstream& oss, const std::string& value)
    {
      oss << value << '|';
    }
  } // namespace

  XICParquetFile::XICParquetFile(const std::string& filename)
    : filename_(filename)
  {
    if (!File::exists(filename_))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename_);
    }
    filenames_.push_back(filename_);
  }

  XICParquetFile::XICParquetFile(const std::vector<std::string>& filenames)
    : filenames_(filenames)
  {
    if (filenames_.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "No parquet files provided", "");
    }
    for (const auto& filename : filenames_)
    {
      if (!File::exists(filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
    }
    filename_ = filenames_.front();
  }

  const std::string& XICParquetFile::getFilename() const
  {
    return filename_;
  }

  const std::vector<std::string>& XICParquetFile::getFilenames() const
  {
    return filenames_;
  }

  void XICParquetFile::load(std::vector<XICChromatogram>& output) const
  {
    getChromatograms(output);
  }

  void XICParquetFile::getChromatograms_(std::vector<XICChromatogram>& output,
                                         const FilterExpression& extra_filter,
                                         Int64 precursor_id,
                                         Int64 transition_id,
                                         const std::string& modified_sequence,
                                         Int64 precursor_charge,
                                         Int64 product_charge,
                                         Int64 ms_level,
                                         Int64 run_id,
                                         const std::string& filter) const
  {
    output.clear();

    std::shared_ptr<arrow::Table> table;

    std::unordered_map<std::string, ColumnType> column_types = {
      {"RUN_ID", ColumnType::INT},
      {"MS_LEVEL", ColumnType::INT},
      {"PRECURSOR_ID", ColumnType::INT},
      {"TRANSITION_ID", ColumnType::INT},
      {"PRECURSOR_CHARGE", ColumnType::INT},
      {"PRODUCT_CHARGE", ColumnType::INT},
      {"DETECTING_TRANSITION", ColumnType::INT},
      {"PRECURSOR_DECOY", ColumnType::INT},
      {"PRODUCT_DECOY", ColumnType::INT},
      {"TRANSITION_ORDINAL", ColumnType::INT},
      {"RT_COMPRESSION", ColumnType::INT},
      {"INTENSITY_COMPRESSION", ColumnType::INT},
      {"SOURCE_FILE", ColumnType::STRING},
      {"MODIFIED_SEQUENCE", ColumnType::STRING},
      {"TRANSITION_TYPE", ColumnType::STRING},
      {"ANNOTATION", ColumnType::STRING}
    };

    FilterExpression args_filter = buildFilterFromArgs_(precursor_id,
                                                       transition_id,
                                                       modified_sequence,
                                                       precursor_charge,
                                                       product_charge,
                                                       ms_level,
                                                       run_id);
    FilterExpression string_filter = PRH::parseFilter(filter, column_types);
    FilterExpression combined_filter = ParquetFilter::merge(args_filter, string_filter, "AND");
    combined_filter = ParquetFilter::merge(combined_filter, extra_filter, "AND");

    const std::string filter_context = filter.empty() ? std::string("<typed/args filter>") : filter;

    bool used_dataset = false;
#ifdef WITH_ARROW_DATASET
    if (!combined_filter.conditions.empty())
    {
      try
      {
        // Prefer dataset pushdown when available; fall back to compute filtering on failure.
        table = PRH::readParquetTableDataset(filenames_, combined_filter, filter_context, filterPruningOptions_());
        used_dataset = true;
      }
      catch (const std::exception& e)
      {
        OPENMS_LOG_WARN << "Arrow Dataset filter failed, falling back to compute filter: "
                        << e.what() << '\n';
        table = nullptr;
        used_dataset = false;
      }
      catch (...)
      {
        OPENMS_LOG_WARN << "Arrow Dataset filter failed, falling back to compute filter.\n";
        table = nullptr;
        used_dataset = false;
      }
    }
#endif
    if (!table)
    {
      table = PRH::readParquetTable(filenames_);
    }
    if (!used_dataset)
    {
      table = PRH::applyArrowFilter(table, combined_filter, filter_context, filterPruningOptions_());
    }

    auto run_id_col = PRH::getColumn(table, "RUN_ID");
    auto source_file_col = PRH::getColumn(table, "SOURCE_FILE", false);
    auto ms_level_col = PRH::getColumn(table, "MS_LEVEL", false);
    auto precursor_id_col = PRH::getColumn(table, "PRECURSOR_ID");
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
    auto rt_data_col = PRH::getColumn(table, "RT_DATA");
    auto intensity_data_col = PRH::getColumn(table, "INTENSITY_DATA");
    auto rt_compression_col = PRH::getColumn(table, "RT_COMPRESSION");
    auto intensity_compression_col = PRH::getColumn(table, "INTENSITY_COMPRESSION");

    const int64_t rows = table->num_rows();
    output.reserve(rows);

    for (int64_t row = 0; row < rows; ++row)
    {
      XICChromatogram chrom;

      // Defensive: RUN_ID is a required column but individual cells may be null in malformed files.
      // Use the shared typed extractor and throw a descriptive error if a required cell is null.
      if (!PRH::getOptionalInt(run_id_col, row, chrom.run_id))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "NULL value in required RUN_ID column",StringUtils::toStr(row));
      }
      chrom.has_precursor_id = PRH::getOptionalInt(precursor_id_col, row, chrom.precursor_id);
      chrom.has_transition_id = PRH::getOptionalInt(transition_id_col, row, chrom.transition_id);
      chrom.has_precursor_charge = PRH::getOptionalInt(precursor_charge_col, row, chrom.precursor_charge);
      chrom.has_product_charge = PRH::getOptionalInt(product_charge_col, row, chrom.product_charge);
      chrom.has_detecting_transition = PRH::getOptionalInt(detecting_transition_col, row, chrom.detecting_transition);
      chrom.has_precursor_decoy = PRH::getOptionalInt(precursor_decoy_col, row, chrom.precursor_decoy);
      chrom.has_product_decoy = PRH::getOptionalInt(product_decoy_col, row, chrom.product_decoy);
      chrom.has_transition_ordinal = PRH::getOptionalInt(transition_ordinal_col, row, chrom.transition_ordinal);

      PRH::getOptionalString(source_file_col, row, chrom.source_file);
      PRH::getOptionalString(modified_sequence_col, row, chrom.modified_sequence);
      PRH::getOptionalString(transition_type_col, row, chrom.transition_type);
      PRH::getOptionalString(annotation_col, row, chrom.annotation);

      Int64 ms_level_value = 0;
      if (PRH::getOptionalInt(ms_level_col, row, ms_level_value))
      {
        chrom.ms_level = ms_level_value;
      }

      if (precursor_id >= 0 && (!chrom.has_precursor_id || chrom.precursor_id != precursor_id)) continue;
      if (transition_id >= 0 && (!chrom.has_transition_id || chrom.transition_id != transition_id)) continue;
      if (!modified_sequence.empty() && chrom.modified_sequence != modified_sequence) continue;
      if (precursor_charge >= 0 && (!chrom.has_precursor_charge || chrom.precursor_charge != precursor_charge)) continue;
      if (product_charge >= 0 && (!chrom.has_product_charge || chrom.product_charge != product_charge)) continue;
      if (ms_level >= 0 && chrom.ms_level != ms_level) continue;
      if (run_id >= 0 && chrom.run_id != run_id) continue;

      Int64 rt_compression = 0;
      Int64 intensity_compression = 0;
      PRH::getOptionalInt(rt_compression_col, row, rt_compression);
      PRH::getOptionalInt(intensity_compression_col, row, intensity_compression);

      const std::string rt_data = PRH::getBinaryView(rt_data_col, row);
      const std::string intensity_data = PRH::getBinaryView(intensity_data_col, row);

      decodeBinary_(rt_data, rt_compression, chrom.rt);
      decodeBinary_(intensity_data, intensity_compression, chrom.intensity);
      if (chrom.rt.size() != chrom.intensity.size())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "RT and intensity array size mismatch",StringUtils::toStr(row));
      }

      output.push_back(std::move(chrom));
    }
  }

  void XICParquetFile::getRuns(std::vector<XICRunInfo>& output) const
  {
    output.clear();

    const std::vector<std::string> columns = {"RUN_ID", "SOURCE_FILE"};
    auto table = PRH::readParquetTableColumns(filenames_, columns);

    auto run_id_col = PRH::getColumn(table, "RUN_ID");
    auto source_file_col = PRH::getColumn(table, "SOURCE_FILE", false);

    std::unordered_set<std::string> seen;
    const int64_t rows = table->num_rows();
    output.reserve(rows);

    for (int64_t row = 0; row < rows; ++row)
    {
      XICRunInfo info;
      if (!PRH::getOptionalInt(run_id_col, row, info.run_id))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "NULL value in required RUN_ID column",StringUtils::toStr(row));
      }
      PRH::getOptionalString(source_file_col, row, info.source_file);

      std::ostringstream oss;
      appendKeyPart_(oss, true, info.run_id);
      appendKeyPart_(oss, info.source_file);
      if (seen.insert(oss.str()).second)
      {
        output.push_back(std::move(info));
      }
    }
  }

  void XICParquetFile::getChromatograms(std::vector<XICChromatogram>& output,
                                        Int64 precursor_id,
                                        Int64 transition_id,
                                        const std::string& modified_sequence,
                                        Int64 precursor_charge,
                                        Int64 product_charge,
                                        Int64 ms_level,
                                        Int64 run_id,
                                        const std::string& filter) const
  {
    FilterExpression empty_filter;
    getChromatograms_(output, empty_filter, precursor_id, transition_id, modified_sequence,
                      precursor_charge, product_charge, ms_level, run_id, filter);
  }

  void XICParquetFile::getChromatograms(std::vector<XICChromatogram>& output,
                                        const XICQuery& query) const
  {
    // forward to the positional overload; named fields make the column mapping unambiguous
    getChromatograms(output, query.precursor_id, query.transition_id, query.modified_sequence,
                      query.precursor_charge, query.product_charge, query.ms_level, query.run_id,
                      query.filter);
  }

  void XICParquetFile::getChromatograms(std::vector<XICChromatogram>& output,
                                        const ParquetFilter& filter) const
  {
    getChromatograms_(output, filter.expression(), -1, -1, "", -1, -1, -1, -1, "");
  }

  void XICParquetFile::getChromatograms(std::vector<XICChromatogram>& output,
                                        const ParquetFilterBuilder& filter) const
  {
    getChromatograms(output, filter.filter());
  }

  void XICParquetFile::getAnalytes(std::vector<XICAnalyte>& output,
                                   const std::vector<std::string>& columns,
                                   bool nest_transitions) const
  {
    output.clear();

    const std::vector<std::string> default_columns = {
      "PRECURSOR_ID",
      "MODIFIED_SEQUENCE",
      "PRECURSOR_CHARGE",
      "PRECURSOR_DECOY",
      "TRANSITION_ID",
      "PRODUCT_CHARGE",
      "TRANSITION_ORDINAL",
      "DETECTING_TRANSITION",
      "PRODUCT_DECOY",
      "TRANSITION_TYPE",
      "ANNOTATION"
    };
    const std::vector<std::string>& requested_columns = columns.empty() ? default_columns : columns;

    std::shared_ptr<arrow::Schema> schema = PRH::readParquetSchemaAllFiles(filenames_);
    std::unordered_set<std::string> schema_columns;
    schema_columns.reserve(schema->num_fields());
    for (const auto& field : schema->fields())
    {
      schema_columns.insert(upper_(std::string(field->name())));
    }

    std::unordered_set<std::string> allowed_columns = {
      "PRECURSOR_ID",
      "MODIFIED_SEQUENCE",
      "PRECURSOR_CHARGE",
      "PRECURSOR_DECOY",
      "TRANSITION_ID",
      "PRODUCT_CHARGE",
      "TRANSITION_ORDINAL",
      "DETECTING_TRANSITION",
      "PRODUCT_DECOY",
      "TRANSITION_TYPE",
      "ANNOTATION"
    };

    std::vector<std::string> normalized_columns;
    normalized_columns.reserve(requested_columns.size());
    for (const auto& name : requested_columns)
    {
      const std::string upper_name = upper_(name);
      if (!allowed_columns.contains(upper_name))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Unsupported analyte column", name);
      }
      if (!schema_columns.contains(upper_name))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Column not found in parquet schema", name);
      }
      normalized_columns.push_back(upper_name);
    }

    std::unordered_set<std::string> requested_set;
    requested_set.reserve(normalized_columns.size());
    for (const auto& name : normalized_columns)
    {
      requested_set.insert(name);
    }

    auto want = [&](const std::string& name) -> bool
    {
      return requested_set.contains(name);
    };

    auto table = PRH::readParquetTableColumns(filenames_, normalized_columns);

    auto precursor_id_col = PRH::getColumn(table, "PRECURSOR_ID", want("PRECURSOR_ID"));
    auto modified_sequence_col = PRH::getColumn(table, "MODIFIED_SEQUENCE", want("MODIFIED_SEQUENCE"));
    auto precursor_charge_col = PRH::getColumn(table, "PRECURSOR_CHARGE", want("PRECURSOR_CHARGE"));
    auto precursor_decoy_col = PRH::getColumn(table, "PRECURSOR_DECOY", want("PRECURSOR_DECOY"));
    auto transition_id_col = PRH::getColumn(table, "TRANSITION_ID", want("TRANSITION_ID"));
    auto product_charge_col = PRH::getColumn(table, "PRODUCT_CHARGE", want("PRODUCT_CHARGE"));
    auto transition_ordinal_col = PRH::getColumn(table, "TRANSITION_ORDINAL", want("TRANSITION_ORDINAL"));
    auto detecting_transition_col = PRH::getColumn(table, "DETECTING_TRANSITION", want("DETECTING_TRANSITION"));
    auto product_decoy_col = PRH::getColumn(table, "PRODUCT_DECOY", want("PRODUCT_DECOY"));
    auto transition_type_col = PRH::getColumn(table, "TRANSITION_TYPE", want("TRANSITION_TYPE"));
    auto annotation_col = PRH::getColumn(table, "ANNOTATION", want("ANNOTATION"));

    const int64_t rows = table->num_rows();
    output.reserve(rows);

    if (!nest_transitions)
    {
      std::unordered_set<std::string> seen;
      for (int64_t row = 0; row < rows; ++row)
      {
        XICAnalyte analyte;

        if (want("PRECURSOR_ID"))
        {
          analyte.has_precursor_id = PRH::getOptionalInt(precursor_id_col, row, analyte.precursor_id);
        }
        if (want("MODIFIED_SEQUENCE"))
        {
          PRH::getOptionalString(modified_sequence_col, row, analyte.modified_sequence);
        }
        if (want("PRECURSOR_CHARGE"))
        {
          analyte.has_precursor_charge = PRH::getOptionalInt(precursor_charge_col, row, analyte.precursor_charge);
        }
        if (want("PRECURSOR_DECOY"))
        {
          analyte.has_precursor_decoy = PRH::getOptionalInt(precursor_decoy_col, row, analyte.precursor_decoy);
        }

        if (want("TRANSITION_ID"))
        {
          analyte.has_transition_id = PRH::getOptionalInt(transition_id_col, row, analyte.transition_id);
        }
        if (want("PRODUCT_CHARGE"))
        {
          analyte.has_product_charge = PRH::getOptionalInt(product_charge_col, row, analyte.product_charge);
        }
        if (want("TRANSITION_ORDINAL"))
        {
          analyte.has_transition_ordinal = PRH::getOptionalInt(transition_ordinal_col, row, analyte.transition_ordinal);
        }
        if (want("DETECTING_TRANSITION"))
        {
          analyte.has_detecting_transition = PRH::getOptionalInt(detecting_transition_col, row, analyte.detecting_transition);
        }
        if (want("PRODUCT_DECOY"))
        {
          analyte.has_product_decoy = PRH::getOptionalInt(product_decoy_col, row, analyte.product_decoy);
        }
        if (want("TRANSITION_TYPE"))
        {
          PRH::getOptionalString(transition_type_col, row, analyte.transition_type);
        }
        if (want("ANNOTATION"))
        {
          PRH::getOptionalString(annotation_col, row, analyte.annotation);
        }

        std::ostringstream oss;
        if (want("PRECURSOR_ID"))
        {
          appendKeyPart_(oss, analyte.has_precursor_id, analyte.precursor_id);
        }
        if (want("MODIFIED_SEQUENCE"))
        {
          appendKeyPart_(oss, analyte.modified_sequence);
        }
        if (want("PRECURSOR_CHARGE"))
        {
          appendKeyPart_(oss, analyte.has_precursor_charge, analyte.precursor_charge);
        }
        if (want("PRECURSOR_DECOY"))
        {
          appendKeyPart_(oss, analyte.has_precursor_decoy, analyte.precursor_decoy);
        }
        if (want("TRANSITION_ID"))
        {
          appendKeyPart_(oss, analyte.has_transition_id, analyte.transition_id);
        }
        if (want("PRODUCT_CHARGE"))
        {
          appendKeyPart_(oss, analyte.has_product_charge, analyte.product_charge);
        }
        if (want("TRANSITION_ORDINAL"))
        {
          appendKeyPart_(oss, analyte.has_transition_ordinal, analyte.transition_ordinal);
        }
        if (want("DETECTING_TRANSITION"))
        {
          appendKeyPart_(oss, analyte.has_detecting_transition, analyte.detecting_transition);
        }
        if (want("PRODUCT_DECOY"))
        {
          appendKeyPart_(oss, analyte.has_product_decoy, analyte.product_decoy);
        }
        if (want("TRANSITION_TYPE"))
        {
          appendKeyPart_(oss, analyte.transition_type);
        }
        if (want("ANNOTATION"))
        {
          appendKeyPart_(oss, analyte.annotation);
        }

        if (seen.insert(oss.str()).second)
        {
          output.push_back(std::move(analyte));
        }
      }
      return;
    }

    std::unordered_map<std::string, Size> precursor_index;
    std::unordered_map<std::string, std::unordered_set<std::string>> transition_seen;

    for (int64_t row = 0; row < rows; ++row)
    {
      XICAnalyte analyte;
      analyte.has_precursor_id = PRH::getOptionalInt(precursor_id_col, row, analyte.precursor_id);
      PRH::getOptionalString(modified_sequence_col, row, analyte.modified_sequence);
      analyte.has_precursor_charge = PRH::getOptionalInt(precursor_charge_col, row, analyte.precursor_charge);
      analyte.has_precursor_decoy = PRH::getOptionalInt(precursor_decoy_col, row, analyte.precursor_decoy);

      std::ostringstream precursor_key;
      if (want("PRECURSOR_ID"))
      {
        appendKeyPart_(precursor_key, analyte.has_precursor_id, analyte.precursor_id);
      }
      if (want("MODIFIED_SEQUENCE"))
      {
        appendKeyPart_(precursor_key, analyte.modified_sequence);
      }
      if (want("PRECURSOR_CHARGE"))
      {
        appendKeyPart_(precursor_key, analyte.has_precursor_charge, analyte.precursor_charge);
      }
      if (want("PRECURSOR_DECOY"))
      {
        appendKeyPart_(precursor_key, analyte.has_precursor_decoy, analyte.precursor_decoy);
      }
      const std::string precursor_key_str = precursor_key.str();

      Size index = 0;
      auto it = precursor_index.find(precursor_key_str);
      if (it == precursor_index.end())
      {
        output.push_back(std::move(analyte));
        index = output.size() - 1;
        precursor_index[precursor_key_str] = index;
      }
      else
      {
        index = it->second;
      }

      Int64 transition_id = 0;
      bool has_transition_id = false;
      Int64 product_charge = 0;
      bool has_product_charge = false;
      Int64 transition_ordinal = 0;
      bool has_transition_ordinal = false;
      Int64 detecting_transition = 0;
      bool has_detecting_transition = false;
      Int64 product_decoy = 0;
      bool has_product_decoy = false;
      std::string transition_type;
      std::string annotation;

      if (want("TRANSITION_ID"))
      {
        has_transition_id = PRH::getOptionalInt(transition_id_col, row, transition_id);
      }
      if (want("PRODUCT_CHARGE"))
      {
        has_product_charge = PRH::getOptionalInt(product_charge_col, row, product_charge);
      }
      if (want("TRANSITION_ORDINAL"))
      {
        has_transition_ordinal = PRH::getOptionalInt(transition_ordinal_col, row, transition_ordinal);
      }
      if (want("DETECTING_TRANSITION"))
      {
        has_detecting_transition = PRH::getOptionalInt(detecting_transition_col, row, detecting_transition);
      }
      if (want("PRODUCT_DECOY"))
      {
        has_product_decoy = PRH::getOptionalInt(product_decoy_col, row, product_decoy);
      }
      if (want("TRANSITION_TYPE"))
      {
        PRH::getOptionalString(transition_type_col, row, transition_type);
      }
      if (want("ANNOTATION"))
      {
        PRH::getOptionalString(annotation_col, row, annotation);
      }

      std::ostringstream transition_key;
      if (want("TRANSITION_ID"))
      {
        appendKeyPart_(transition_key, has_transition_id, transition_id);
      }
      if (want("PRODUCT_CHARGE"))
      {
        appendKeyPart_(transition_key, has_product_charge, product_charge);
      }
      if (want("TRANSITION_ORDINAL"))
      {
        appendKeyPart_(transition_key, has_transition_ordinal, transition_ordinal);
      }
      if (want("DETECTING_TRANSITION"))
      {
        appendKeyPart_(transition_key, has_detecting_transition, detecting_transition);
      }
      if (want("PRODUCT_DECOY"))
      {
        appendKeyPart_(transition_key, has_product_decoy, product_decoy);
      }
      if (want("TRANSITION_TYPE"))
      {
        appendKeyPart_(transition_key, transition_type);
      }
      if (want("ANNOTATION"))
      {
        appendKeyPart_(transition_key, annotation);
      }
      const std::string transition_key_str = transition_key.str();

      if (transition_seen[precursor_key_str].insert(transition_key_str).second)
      {
        if (want("TRANSITION_ID"))
        {
          output[index].transition_ids.push_back(has_transition_id ? transition_id : -1);
        }
        if (want("PRODUCT_CHARGE"))
        {
          output[index].product_charges.push_back(has_product_charge ? product_charge : -1);
        }
        if (want("TRANSITION_ORDINAL"))
        {
          output[index].transition_ordinals.push_back(has_transition_ordinal ? transition_ordinal : -1);
        }
        if (want("DETECTING_TRANSITION"))
        {
          output[index].detecting_transitions.push_back(has_detecting_transition ? detecting_transition : -1);
        }
        if (want("PRODUCT_DECOY"))
        {
          output[index].product_decoys.push_back(has_product_decoy ? product_decoy : -1);
        }
        if (want("TRANSITION_TYPE"))
        {
          output[index].transition_types.push_back(transition_type);
        }
        if (want("ANNOTATION"))
        {
          output[index].annotations.push_back(annotation);
        }
      }
    }
  }

  void XICParquetFile::getColumns(std::vector<std::string>& output) const
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
