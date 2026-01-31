// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/XICParquetFile.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/MSNumpressCoder.h>
#include <OpenMS/FORMAT/ZlibCompression.h>
#include <OpenMS/SYSTEM/File.h>

#ifdef WITH_PARQUET
#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/reader.h>
#endif

#include <cstring>
#include <cctype>
#include <unordered_map>

namespace OpenMS
{
  namespace
  {
#ifdef WITH_PARQUET
    std::shared_ptr<arrow::Table> readParquetTable_(const String& filename)
    {
      auto infile_result = arrow::io::ReadableFile::Open(std::string(filename));
      if (!infile_result.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to open parquet file", filename);
      }
      std::shared_ptr<arrow::io::ReadableFile> infile = *infile_result;

      auto reader_result = parquet::arrow::OpenFile(infile, arrow::default_memory_pool());
      if (!reader_result.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to create parquet reader", filename);
      }
      std::unique_ptr<parquet::arrow::FileReader> reader = std::move(reader_result.ValueOrDie());

      std::shared_ptr<arrow::Table> table;
      auto read_status = reader->ReadTable(&table);
      if (!read_status.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to read parquet table", filename);
      }

      auto combined = table->CombineChunks(arrow::default_memory_pool());
      if (!combined.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to combine parquet chunks", filename);
      }

      return *combined;
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
      return column->chunk(0);
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
        case 0: // no compression
        {
          const size_t count = data.size() / sizeof(double);
          output.resize(count);
          std::memcpy(output.data(), data.c_str(), count * sizeof(double));
          break;
        }
        case 1: // zlib only
        {
          String decoded;
          ZlibCompression::uncompressString(data, decoded);
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
          String decoded;
          ZlibCompression::uncompressString(data, decoded);
          decodeNumpress(MSNumpressCoder::LINEAR, decoded);
          break;
        }
        case 6: // np-slof + zlib
        {
          String decoded;
          ZlibCompression::uncompressString(data, decoded);
          decodeNumpress(MSNumpressCoder::SLOF, decoded);
          break;
        }
        case 7: // np-pic + zlib
        {
          String decoded;
          ZlibCompression::uncompressString(data, decoded);
          decodeNumpress(MSNumpressCoder::PIC, decoded);
          break;
        }
        default:
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Unsupported compression type", String(compression));
      }
    }

    enum class ColumnType
    {
      INT,
      STRING
    };

    struct Condition
    {
      String column;
      String op;
      std::vector<String> values;
      ColumnType type{ColumnType::INT};
    };

    struct FilterExpression
    {
      std::vector<Condition> conditions;
      std::vector<String> connectors; // "AND" / "OR"
    };

    String upper_(const String& input)
    {
      String out = input;
      for (Size i = 0; i < out.size(); ++i)
      {
        out[i] = static_cast<char>(std::toupper(out[i]));
      }
      return out;
    }

    std::vector<String> tokenize_(const String& expr)
    {
      std::vector<String> tokens;
      String current;
      bool in_string = false;
      char quote_char = '\0';

      for (Size i = 0; i < expr.size(); ++i)
      {
        const char c = expr[i];
        if (in_string)
        {
          if (c == quote_char)
          {
            tokens.push_back(current);
            current.clear();
            in_string = false;
            quote_char = '\0';
          }
          else
          {
            current += c;
          }
          continue;
        }

        if (c == '"' || c == '\'')
        {
          if (!current.empty())
          {
            tokens.push_back(current);
            current.clear();
          }
          in_string = true;
          quote_char = c;
          continue;
        }

        if (std::isspace(static_cast<unsigned char>(c)))
        {
          if (!current.empty())
          {
            tokens.push_back(current);
            current.clear();
          }
          continue;
        }

        if (c == '[' || c == ']' || c == ',' || c == '(' || c == ')')
        {
          if (!current.empty())
          {
            tokens.push_back(current);
            current.clear();
          }
          tokens.emplace_back(String(c));
          continue;
        }

        if (c == '=' || c == '!' || c == '<' || c == '>')
        {
          if (!current.empty())
          {
            tokens.push_back(current);
            current.clear();
          }
          String op;
          op += c;
          if (i + 1 < expr.size())
          {
            const char next = expr[i + 1];
            if (next == '=')
            {
              op += next;
              ++i;
            }
          }
          tokens.push_back(op);
          continue;
        }

        current += c;
      }

      if (!current.empty())
      {
        tokens.push_back(current);
      }

      return tokens;
    }

    FilterExpression parseFilter_(const String& filter,
                                  const std::unordered_map<String, ColumnType>& column_types)
    {
      FilterExpression expr;
      if (filter.empty())
      {
        return expr;
      }

      auto tokens = tokenize_(filter);
      Size i = 0;
      while (i < tokens.size())
      {
        String column = upper_(tokens[i++]);
        if (i >= tokens.size()) break;
        String op = upper_(tokens[i++]);

        Condition cond;
        cond.column = column;
        cond.op = op;
        auto type_it = column_types.find(column);
        cond.type = (type_it == column_types.end()) ? ColumnType::INT : type_it->second;

        if (op == "IN")
        {
          if (i < tokens.size() && tokens[i] == "[") ++i;
          while (i < tokens.size() && tokens[i] != "]")
          {
            if (tokens[i] != ",")
            {
              cond.values.push_back(tokens[i]);
            }
            ++i;
          }
          if (i < tokens.size() && tokens[i] == "]") ++i;
        }
        else
        {
          if (i < tokens.size())
          {
            cond.values.push_back(tokens[i++]);
          }
        }

        expr.conditions.push_back(cond);
        if (i < tokens.size())
        {
          String connector = upper_(tokens[i]);
          if (connector == "AND" || connector == "OR")
          {
            expr.connectors.push_back(connector);
            ++i;
          }
        }
      }

      return expr;
    }

    bool matchCondition_(const Condition& cond,
                         int64_t int_value,
                         bool has_value,
                         const String& str_value,
                         bool has_str)
    {
      if (cond.type == ColumnType::INT)
      {
        if (!has_value)
        {
          return false;
        }
        if (cond.op == "IN")
        {
          for (const auto& value : cond.values)
          {
            if (value.toInt64() == int_value) return true;
          }
          return false;
        }
        if (cond.values.empty()) return false;
        Int64 rhs = cond.values.front().toInt64();
        if (cond.op == "=") return int_value == rhs;
        if (cond.op == "!=") return int_value != rhs;
        if (cond.op == "<") return int_value < rhs;
        if (cond.op == "<=") return int_value <= rhs;
        if (cond.op == ">") return int_value > rhs;
        if (cond.op == ">=") return int_value >= rhs;
        return false;
      }
      else
      {
        if (!has_str)
        {
          return false;
        }
        if (cond.op == "IN")
        {
          for (const auto& value : cond.values)
          {
            if (str_value == value) return true;
          }
          return false;
        }
        if (cond.values.empty()) return false;
        const String& rhs = cond.values.front();
        if (cond.op == "=") return str_value == rhs;
        if (cond.op == "!=") return str_value != rhs;
        return false;
      }
    }

    bool evaluateFilter_(const FilterExpression& expr,
                         const std::unordered_map<String, std::pair<Int64, bool>>& int_values,
                         const std::unordered_map<String, std::pair<String, bool>>& str_values)
    {
      if (expr.conditions.empty())
      {
        return true;
      }

      bool result = true;
      for (Size i = 0; i < expr.conditions.size(); ++i)
      {
        const Condition& cond = expr.conditions[i];
        bool cond_result = false;
        if (cond.type == ColumnType::INT)
        {
          auto it = int_values.find(cond.column);
          if (it != int_values.end())
          {
            cond_result = matchCondition_(cond, it->second.first, it->second.second, "", false);
          }
        }
        else
        {
          auto it = str_values.find(cond.column);
          if (it != str_values.end())
          {
            cond_result = matchCondition_(cond, 0, false, it->second.first, it->second.second);
          }
        }

        if (i == 0)
        {
          result = cond_result;
        }
        else if (i - 1 < expr.connectors.size())
        {
          if (expr.connectors[i - 1] == "AND")
          {
            result = result && cond_result;
          }
          else
          {
            result = result || cond_result;
          }
        }
      }
      return result;
    }
#endif
  } // namespace

  XICParquetFile::XICParquetFile(const String& filename)
    : filename_(filename)
  {
    if (!File::exists(filename_))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename_);
    }
  }

  const String& XICParquetFile::getFilename() const
  {
    return filename_;
  }

  void XICParquetFile::load(std::vector<XICChromatogram>& output) const
  {
    getChromatograms(output);
  }

  void XICParquetFile::getChromatograms(std::vector<XICChromatogram>& output,
                                        Int64 precursor_id,
                                        Int64 transition_id,
                                        const String& modified_sequence,
                                        Int64 precursor_charge,
                                        Int64 product_charge,
                                        Int64 ms_level,
                                        Int64 run_id,
                                        const String& filter) const
  {
#ifndef WITH_PARQUET
    (void)output;
    (void)precursor_id;
    (void)transition_id;
    (void)modified_sequence;
    (void)precursor_charge;
    (void)product_charge;
    (void)ms_level;
    (void)run_id;
    (void)filter;
    throw Exception::MissingFeature(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "OpenMS was built without Parquet support");
#else
    output.clear();

    auto table = readParquetTable_(filename_);

    auto run_id_col = getColumn_(table, "RUN_ID");
    auto source_file_col = getColumn_(table, "SOURCE_FILE", false);
    auto ms_level_col = getColumn_(table, "MS_LEVEL", false);
    auto precursor_id_col = getColumn_(table, "PRECURSOR_ID");
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
    auto rt_data_col = getColumn_(table, "RT_DATA");
    auto intensity_data_col = getColumn_(table, "INTENSITY_DATA");
    auto rt_compression_col = getColumn_(table, "RT_COMPRESSION");
    auto intensity_compression_col = getColumn_(table, "INTENSITY_COMPRESSION");

    std::unordered_map<String, ColumnType> column_types = {
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

    FilterExpression parsed_filter = parseFilter_(filter, column_types);

    const int64_t rows = table->num_rows();
    output.reserve(rows);

    for (int64_t row = 0; row < rows; ++row)
    {
      XICChromatogram chrom;

      chrom.run_id = std::static_pointer_cast<arrow::Int64Array>(run_id_col)->Value(row);
      chrom.has_precursor_id = getOptionalInt_(precursor_id_col, row, chrom.precursor_id);
      chrom.has_transition_id = getOptionalInt_(transition_id_col, row, chrom.transition_id);
      chrom.has_precursor_charge = getOptionalInt_(precursor_charge_col, row, chrom.precursor_charge);
      chrom.has_product_charge = getOptionalInt_(product_charge_col, row, chrom.product_charge);
      chrom.has_detecting_transition = getOptionalInt_(detecting_transition_col, row, chrom.detecting_transition);
      chrom.has_precursor_decoy = getOptionalInt_(precursor_decoy_col, row, chrom.precursor_decoy);
      chrom.has_product_decoy = getOptionalInt_(product_decoy_col, row, chrom.product_decoy);
      chrom.has_transition_ordinal = getOptionalInt_(transition_ordinal_col, row, chrom.transition_ordinal);

      getOptionalString_(source_file_col, row, chrom.source_file);
      getOptionalString_(modified_sequence_col, row, chrom.modified_sequence);
      getOptionalString_(transition_type_col, row, chrom.transition_type);
      getOptionalString_(annotation_col, row, chrom.annotation);

      Int64 ms_level_value = 0;
      if (getOptionalInt_(ms_level_col, row, ms_level_value))
      {
        chrom.ms_level = ms_level_value;
      }

      std::unordered_map<String, std::pair<Int64, bool>> int_values = {
        {"RUN_ID", {chrom.run_id, true}},
        {"MS_LEVEL", {chrom.ms_level, ms_level_col != nullptr}},
        {"PRECURSOR_ID", {chrom.precursor_id, chrom.has_precursor_id}},
        {"TRANSITION_ID", {chrom.transition_id, chrom.has_transition_id}},
        {"PRECURSOR_CHARGE", {chrom.precursor_charge, chrom.has_precursor_charge}},
        {"PRODUCT_CHARGE", {chrom.product_charge, chrom.has_product_charge}},
        {"DETECTING_TRANSITION", {chrom.detecting_transition, chrom.has_detecting_transition}},
        {"PRECURSOR_DECOY", {chrom.precursor_decoy, chrom.has_precursor_decoy}},
        {"PRODUCT_DECOY", {chrom.product_decoy, chrom.has_product_decoy}},
        {"TRANSITION_ORDINAL", {chrom.transition_ordinal, chrom.has_transition_ordinal}}
      };

      std::unordered_map<String, std::pair<String, bool>> str_values = {
        {"SOURCE_FILE", {chrom.source_file, !chrom.source_file.empty()}},
        {"MODIFIED_SEQUENCE", {chrom.modified_sequence, !chrom.modified_sequence.empty()}},
        {"TRANSITION_TYPE", {chrom.transition_type, !chrom.transition_type.empty()}},
        {"ANNOTATION", {chrom.annotation, !chrom.annotation.empty()}}
      };

      if (precursor_id >= 0 && (!chrom.has_precursor_id || chrom.precursor_id != precursor_id)) continue;
      if (transition_id >= 0 && (!chrom.has_transition_id || chrom.transition_id != transition_id)) continue;
      if (!modified_sequence.empty() && chrom.modified_sequence != modified_sequence) continue;
      if (precursor_charge >= 0 && (!chrom.has_precursor_charge || chrom.precursor_charge != precursor_charge)) continue;
      if (product_charge >= 0 && (!chrom.has_product_charge || chrom.product_charge != product_charge)) continue;
      if (ms_level >= 0 && chrom.ms_level != ms_level) continue;
      if (run_id >= 0 && chrom.run_id != run_id) continue;
      if (!evaluateFilter_(parsed_filter, int_values, str_values)) continue;

      Int64 rt_compression = 0;
      Int64 intensity_compression = 0;
      getOptionalInt_(rt_compression_col, row, rt_compression);
      getOptionalInt_(intensity_compression_col, row, intensity_compression);

      const String rt_data = getBinaryView_(rt_data_col, row);
      const String intensity_data = getBinaryView_(intensity_data_col, row);

      decodeBinary_(rt_data, rt_compression, chrom.rt);
      decodeBinary_(intensity_data, intensity_compression, chrom.intensity);

      output.push_back(std::move(chrom));
    }
#endif
  }
} // namespace OpenMS
