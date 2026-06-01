// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/XIMParquetFile.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/ParquetFilter.h>
#include <OpenMS/FORMAT/MSNumpressCoder.h>
#include <OpenMS/FORMAT/ZlibCompression.h>
#include <OpenMS/SYSTEM/File.h>

#include <arrow/api.h>
#include <arrow/compute/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/reader.h>
#ifdef WITH_ARROW_DATASET
#include <arrow/dataset/api.h>
#include <arrow/filesystem/api.h>
#endif

#include <cstring>
#include <cctype>
#include <cmath>
#include <limits>
#include <exception>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

namespace OpenMS
{
  namespace
  {
    /// Read a single parquet file into an Arrow table.
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

      std::shared_ptr<arrow::Table> table;
      auto read_status = reader->ReadTable(&table);
      if (!read_status.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to read parquet table: " + read_status.ToString(), filename);
      }

      auto combined = table->CombineChunks(arrow::default_memory_pool());
      if (!combined.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to combine parquet chunks", filename);
      }

      return *combined;
    }

    /// Read the parquet schema from a single file.
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

    /// Read parquet schema from all files and ensure they are compatible.
    std::shared_ptr<arrow::Schema> readParquetSchemaAllFiles_(const std::vector<String>& filenames)
    {
      if (filenames.empty())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "No parquet files provided for schema validation", "");
      }
      auto base = readParquetSchema_(filenames.front());
      for (size_t i = 1; i < filenames.size(); ++i)
      {
        auto other = readParquetSchema_(filenames[i]);
        // Use Arrow Schema equality check to ensure fields and types match.
        if (!other->Equals(*base))
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Parquet input files have incompatible schemas: '" + filenames.front() + "' vs '" + filenames[i] + "'", "");
        }
      }
      return base;
    }

    /// Read selected columns from a parquet file into an Arrow table.
    std::shared_ptr<arrow::Table> readParquetTableColumns_(const String& filename,
                                                           const std::vector<String>& columns)
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
      auto schema_status = reader->GetSchema(&schema);
      if (!schema_status.ok() || !schema)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to read parquet schema", filename);
      }

      std::vector<int> column_indices;
      column_indices.reserve(columns.size());
      for (const auto& name : columns)
      {
        const int idx = schema->GetFieldIndex(std::string(name));
        if (idx >= 0)
        {
          column_indices.push_back(idx);
        }
      }
      if (column_indices.empty())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "No requested columns found in parquet file '" +
                                            filename + "'");
      }

      std::shared_ptr<arrow::Table> table;
      auto read_status = reader->ReadTable(column_indices, &table);
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

    /// Concatenate multiple Arrow tables (multiple input files).
    std::shared_ptr<arrow::Table> concatenateTables_(const std::vector<std::shared_ptr<arrow::Table>>& tables)
    {
      if (tables.empty())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "No parquet tables provided", "");
      }
      if (tables.size() == 1)
      {
        return tables.front();
      }

      // Concatenate multiple file tables into a single logical view.
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

    /// Read and concatenate multiple parquet files.
    std::shared_ptr<arrow::Table> readParquetTable_(const std::vector<String>& filenames)
    {
      std::vector<std::shared_ptr<arrow::Table>> tables;
      tables.reserve(filenames.size());
      for (const auto& filename : filenames)
      {
        tables.push_back(readParquetTable_(filename));
      }
      return concatenateTables_(tables);
    }

    /// Read and concatenate selected columns from multiple parquet files.
    std::shared_ptr<arrow::Table> readParquetTableColumns_(const std::vector<String>& filenames,
                                                           const std::vector<String>& columns)
    {
      std::vector<std::shared_ptr<arrow::Table>> tables;
      tables.reserve(filenames.size());
      for (const auto& filename : filenames)
      {
        tables.push_back(readParquetTableColumns_(filename, columns));
      }
      return concatenateTables_(tables);
    }

    /// Fetch a named column or throw if required and missing.
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

    /// Read optional Int64 value at row; returns false if null or absent.
    // Forward declarations for helpers defined later in this translation unit.
    bool getIntValueFromArray_(const std::shared_ptr<arrow::Array>& array, int64_t row, Int64& value);
    bool getStringValueFromArray_(const std::shared_ptr<arrow::Array>& array, int64_t row, String& value);
    bool getOptionalInt_(const std::shared_ptr<arrow::Array>& array, int64_t row, Int64& value)
    {
      // Delegate to the type-checked helper which also validates integer column types.
      try
      {
        return getIntValueFromArray_(array, row, value);
      }
      catch (const Exception::BaseException&)
      {
        // Re-throw to preserve original error behaviour for unsupported types.
        throw;
      }
    }

    /// Read optional Double value at row; returns false if null or absent.
    bool getOptionalDouble_(const std::shared_ptr<arrow::Array>& array, int64_t row, double& value)
    {
      if (!array || array->IsNull(row))
      {
        return false;
      }
      switch (array->type_id())
      {
        case arrow::Type::DOUBLE:
          value = static_cast<arrow::DoubleArray*>(array.get())->Value(row);
          return true;
        case arrow::Type::FLOAT:
          value = static_cast<double>(static_cast<arrow::FloatArray*>(array.get())->Value(row));
          return true;
        default:
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Unsupported floating-point column type", String(array->type()->ToString()));
      }
    }

    /// Read optional String value at row; returns false if null or absent.
    bool getOptionalString_(const std::shared_ptr<arrow::Array>& array, int64_t row, String& value)
    {
      // Reuse the validated helper that supports STRING and LARGE_STRING
      try
      {
        return getStringValueFromArray_(array, row, value);
      }
      catch (const Exception::BaseException&)
      {
        throw;
      }
    }

    /// Read a binary cell as a String view.
    String getBinaryView_(const std::shared_ptr<arrow::Array>& array, int64_t row)
    {
      if (!array || array->IsNull(row))
      {
        return String();
      }
      switch (array->type_id())
      {
        case arrow::Type::BINARY:
        {
          auto typed = static_cast<arrow::BinaryArray*>(array.get());
          const auto view = typed->GetView(row);
          {
            const size_t vsize = view.size();
            if (vsize > static_cast<size_t>(std::numeric_limits<Size>::max()))
            {
              throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Binary data too large", String(vsize));
            }
            return String(view.data(), static_cast<Size>(vsize));
          }
        }
        case arrow::Type::LARGE_BINARY:
        {
          auto typed = static_cast<arrow::LargeBinaryArray*>(array.get());
          const auto view = typed->GetView(row);
          {
            const size_t vsize = view.size();
            if (vsize > static_cast<size_t>(std::numeric_limits<Size>::max()))
            {
              throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Binary data too large", String(vsize));
            }
            return String(view.data(), static_cast<Size>(vsize));
          }
        }
        default:
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Unsupported binary column type", String(array->type()->ToString()));
      }
    }

    /// Decode mobility/intensity arrays stored as compressed binary blobs.
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
          if (data.size() % sizeof(double) != 0)
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Invalid mobility/intensity data size (not divisible by 8)", String(data.size()));
          }
          const size_t count = data.size() / sizeof(double);
          output.resize(count);
          std::memcpy(output.data(), data.c_str(), count * sizeof(double));
          break;
        }
        case 1: // zlib only
        {
          String decoded;
          ZlibCompression::uncompressString(data, decoded);
          if (decoded.size() % sizeof(double) != 0)
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Invalid decompressed mobility/intensity data size", String(decoded.size()));
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

    /// Uppercase helper for column normalization.
    String upper_(const String& input)
    {
      String out = input;
      for (Size i = 0; i < out.size(); ++i)
      {
        out[i] = static_cast<char>(std::toupper(out[i]));
      }
      return out;
    }

    /// Tokenize a filter string into operators/identifiers/values.
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
          tokens.emplace_back(c);
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

        if (c == '&' || c == '|')
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
            if (next == c)
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

      if (in_string)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Unclosed quote in filter expression", expr);
      }

      return tokens;
    }

    /// Parse filter text into a structured expression.
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
        if (i >= tokens.size())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Missing filter operator for column", column);
        }
        String op = upper_(tokens[i++]);
        if (op == "==")
        {
          op = "=";
        }
        if (op != "=" && op != "!=" && op != "<" && op != "<=" && op != ">" && op != ">=" && op != "IN")
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Unsupported filter operator", op);
        }

        Condition cond;
        cond.column = column;
        cond.op = op;
        auto type_it = column_types.find(column);
        cond.type = (type_it == column_types.end()) ? ColumnType::INT : type_it->second;

        if (op == "IN")
        {
          String close_token;
          if (i < tokens.size() && (tokens[i] == "[" || tokens[i] == "("))
          {
            close_token = (tokens[i] == "[") ? "]" : ")";
            ++i;
          }
          else
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "IN operator expects a bracketed value list", filter);
          }
          while (i < tokens.size() && tokens[i] != close_token)
          {
            if (tokens[i] != ",")
            {
              cond.values.push_back(tokens[i]);
            }
            ++i;
          }
          if (cond.values.empty())
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "IN operator expects at least one value", filter);
          }
          if (i >= tokens.size() || tokens[i] != close_token)
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Unclosed IN list", filter);
          }
          ++i;
        }
        else
        {
          if (i >= tokens.size())
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Missing filter value for operator", op);
          }
          cond.values.push_back(tokens[i++]);
        }

        expr.conditions.push_back(cond);
        if (i < tokens.size())
        {
          String connector = upper_(tokens[i]);
          bool consumed_connector = false;
          if (connector == "AND" || connector == "&&" || connector == "&")
          {
            expr.connectors.push_back("AND");
            ++i;
            consumed_connector = true;
          }
          else if (connector == "OR" || connector == "||" || connector == "|")
          {
            expr.connectors.push_back("OR");
            ++i;
            consumed_connector = true;
          }
          if (consumed_connector && i >= tokens.size())
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Dangling connector in filter expression", filter);
          }
        }
      }

      return expr;
    }

    FilterExpression buildFilterFromArgs_(Int64 precursor_id,
                                          Int64 transition_id,
                                          const String& modified_sequence,
                                          Int64 precursor_charge,
                                          Int64 product_charge,
                                          Int64 ms_level,
                                          Int64 run_id,
                                          const String& mobilogram_type,
                                          Int64 feature_id,
                                          double feature_rt)
    {
      ParquetFilter builder;
      if (precursor_id >= 0) builder.eq("PRECURSOR_ID", precursor_id);
      if (transition_id >= 0) builder.eq("TRANSITION_ID", transition_id);
      if (!modified_sequence.empty()) builder.eq("MODIFIED_SEQUENCE", modified_sequence);
      if (precursor_charge >= 0) builder.eq("PRECURSOR_CHARGE", precursor_charge);
      if (product_charge >= 0) builder.eq("PRODUCT_CHARGE", product_charge);
      if (ms_level >= 0) builder.eq("MS_LEVEL", ms_level);
      if (run_id >= 0) builder.eq("RUN_ID", run_id);
      if (!mobilogram_type.empty()) builder.eq("MOBILOGRAM_TYPE", mobilogram_type);
      if (feature_id >= 0) builder.eq("FEATURE_ID", feature_id);
      (void)feature_rt;
      return builder.expression();
    }

    /// Map upper-case column names to their schema names.
    std::unordered_map<String, String> buildSchemaNameMap_(const std::shared_ptr<arrow::Schema>& schema)
    {
      std::unordered_map<String, String> name_map;
      if (!schema)
      {
        return name_map;
      }
      name_map.reserve(schema->num_fields());
      for (const auto& field : schema->fields())
      {
        name_map[upper_(String(field->name()))] = String(field->name());
      }
      return name_map;
    }

    /// Join columns for diagnostics/log output.
    String joinColumns_(const std::vector<String>& columns)
    {
      std::ostringstream oss;
      for (Size i = 0; i < columns.size(); ++i)
      {
        if (i > 0) oss << ", ";
        oss << columns[i];
      }
      return String(oss.str());
    }

    /// Normalize column names to schema and drop unsupported columns.
    void normalizeAndPruneColumns_(FilterExpression& expr,
                                   const std::shared_ptr<arrow::Schema>& schema,
                                   std::vector<String>& dropped_columns)
    {
      dropped_columns.clear();
      if (expr.conditions.empty())
      {
        return;
      }

      auto name_map = buildSchemaNameMap_(schema);
      // Unsupported columns that cannot be filtered on.
      // MOBILITY and INTENSITY are stored as binary blobs in parquet, so we cannot filter on them here.
      // FEATURE_RT is a floating-point column and currently not supported by the integer/string filter DSL.
      const std::unordered_set<String> unsupported = {"MOBILITY", "INTENSITY", "FEATURE_RT"};

      std::vector<Condition> kept;
      kept.reserve(expr.conditions.size());
      std::vector<String> new_connectors;
      new_connectors.reserve(expr.connectors.size());

      for (Size i = 0; i < expr.conditions.size(); ++i)
      {
        Condition cond = expr.conditions[i];
        const String original = cond.column;
        auto it = name_map.find(upper_(cond.column));
        if (it != name_map.end())
        {
          cond.column = it->second;
        }

        const bool unsupported_col = unsupported.count(upper_(cond.column)) > 0;
        const bool exists_in_schema = name_map.find(upper_(cond.column)) != name_map.end();
        if (unsupported_col || !exists_in_schema)
        {
          dropped_columns.push_back(original);
          continue;
        }

        if (!kept.empty())
        {
          const Size prev_index = i - 1;
          if (prev_index < expr.connectors.size())
          {
            new_connectors.push_back(expr.connectors[prev_index]);
          }
          else
          {
            new_connectors.push_back("AND");
          }
        }
        kept.push_back(cond);
      }

      expr.conditions.swap(kept);
      expr.connectors.swap(new_connectors);
    }

    /// Append an integer part to a composite key.
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

    /// Append a string part to a composite key.
    void appendKeyPart_(std::ostringstream& oss, const String& value)
    {
      oss << value << '|';
    }

    /// Build a single Arrow condition expression from one clause.
    arrow::compute::Expression buildConditionExpression_(const Condition& cond)
    {
      auto field = arrow::compute::field_ref(std::string(cond.column));

      if (cond.op == "IN")
      {
        if (cond.type == ColumnType::INT)
        {
          arrow::Int64Builder builder;
          for (const auto& value : cond.values)
          {
            auto status = builder.Append(value.toInt64());
            if (!status.ok())
            {
              throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Failed to build INT64 value set", status.ToString());
            }
          }
          std::shared_ptr<arrow::Array> value_set;
          auto status = builder.Finish(&value_set);
          if (!status.ok())
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Failed to finalize INT64 value set", status.ToString());
          }
          auto options = std::make_shared<arrow::compute::SetLookupOptions>(arrow::Datum(value_set));
          return arrow::compute::call("is_in", {field}, options);
        }
        else
        {
          arrow::StringBuilder builder;
          for (const auto& value : cond.values)
          {
            auto status = builder.Append(std::string(value));
            if (!status.ok())
            {
              throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Failed to build STRING value set", status.ToString());
            }
          }
          std::shared_ptr<arrow::Array> value_set;
          auto status = builder.Finish(&value_set);
          if (!status.ok())
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Failed to finalize STRING value set", status.ToString());
          }
          auto options = std::make_shared<arrow::compute::SetLookupOptions>(arrow::Datum(value_set));
          return arrow::compute::call("is_in", {field}, options);
        }
      }

      if (cond.values.empty())
      {
        return arrow::compute::literal(false);
      }

      arrow::compute::Expression literal;
      if (cond.type == ColumnType::INT)
      {
        literal = arrow::compute::literal(static_cast<int64_t>(cond.values.front().toInt64()));
      }
      else
      {
        literal = arrow::compute::literal(std::string(cond.values.front()));
      }

      if (cond.op == "=") return arrow::compute::equal(field, literal);
      if (cond.op == "!=") return arrow::compute::not_equal(field, literal);
      if (cond.op == "<") return arrow::compute::less(field, literal);
      if (cond.op == "<=") return arrow::compute::less_equal(field, literal);
      if (cond.op == ">") return arrow::compute::greater(field, literal);
      if (cond.op == ">=") return arrow::compute::greater_equal(field, literal);
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Unsupported filter operator", cond.op);
    }

    /// Combine conditions into a single Arrow expression.
    arrow::compute::Expression buildFilterExpression_(const FilterExpression& expr)
    {
      if (expr.conditions.empty())
      {
        return arrow::compute::literal(true);
      }

      std::vector<arrow::compute::Expression> cond_exprs;
      cond_exprs.reserve(expr.conditions.size());
      for (const auto& cond : expr.conditions)
      {
        cond_exprs.push_back(buildConditionExpression_(cond));
      }

      arrow::compute::Expression result = cond_exprs.front();
      for (Size i = 1; i < cond_exprs.size(); ++i)
      {
        if (i - 1 < expr.connectors.size() && expr.connectors[i - 1] == "OR")
        {
          result = arrow::compute::or_(result, cond_exprs[i]);
        }
        else
        {
          result = arrow::compute::and_(result, cond_exprs[i]);
        }
      }
      return result;
    }

    struct ManualCondition
    {
      Condition cond;
      int index{-1};
      std::vector<Int64> int_values;
    };

    bool getIntValueFromArray_(const std::shared_ptr<arrow::Array>& array, int64_t row, Int64& value)
    {
      if (!array || array->IsNull(row))
      {
        return false;
      }

      switch (array->type_id())
      {
        case arrow::Type::INT64:
          value = static_cast<arrow::Int64Array*>(array.get())->Value(row);
          return true;
        case arrow::Type::INT32:
          value = static_cast<arrow::Int32Array*>(array.get())->Value(row);
          return true;
        case arrow::Type::INT16:
          value = static_cast<arrow::Int16Array*>(array.get())->Value(row);
          return true;
        case arrow::Type::INT8:
          value = static_cast<arrow::Int8Array*>(array.get())->Value(row);
          return true;
        case arrow::Type::UINT64:
          value = static_cast<arrow::UInt64Array*>(array.get())->Value(row);
          return true;
        case arrow::Type::UINT32:
          value = static_cast<arrow::UInt32Array*>(array.get())->Value(row);
          return true;
        case arrow::Type::UINT16:
          value = static_cast<arrow::UInt16Array*>(array.get())->Value(row);
          return true;
        case arrow::Type::UINT8:
          value = static_cast<arrow::UInt8Array*>(array.get())->Value(row);
          return true;
        case arrow::Type::BOOL:
          value = static_cast<arrow::BooleanArray*>(array.get())->Value(row) ? 1 : 0;
          return true;
        default:
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Unsupported integer column type", String(array->type()->ToString()));
      }
    }

    bool getStringValueFromArray_(const std::shared_ptr<arrow::Array>& array, int64_t row, String& value)
    {
      if (!array || array->IsNull(row))
      {
        return false;
      }

      switch (array->type_id())
      {
        case arrow::Type::STRING:
          value = static_cast<arrow::StringArray*>(array.get())->GetString(row);
          return true;
        case arrow::Type::LARGE_STRING:
          value = static_cast<arrow::LargeStringArray*>(array.get())->GetString(row);
          return true;
        default:
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Unsupported string column type", String(array->type()->ToString()));
      }
    }

    bool evaluateCondition_(const ManualCondition& cond,
                            const std::shared_ptr<arrow::RecordBatch>& batch,
                            int64_t row)
    {
      if (!batch || cond.index < 0 || cond.index >= batch->num_columns())
      {
        return false;
      }

      const auto& array = batch->column(cond.index);
      if (cond.cond.type == ColumnType::INT)
      {
        Int64 lhs = 0;
        if (!getIntValueFromArray_(array, row, lhs))
        {
          return false;
        }

        if (cond.cond.op == "IN")
        {
          for (const auto& rhs : cond.int_values)
          {
            if (lhs == rhs) return true;
          }
          return false;
        }

        if (cond.int_values.empty())
        {
          return false;
        }
        const Int64 rhs = cond.int_values.front();
        if (cond.cond.op == "=") return lhs == rhs;
        if (cond.cond.op == "!=") return lhs != rhs;
        if (cond.cond.op == "<") return lhs < rhs;
        if (cond.cond.op == "<=") return lhs <= rhs;
        if (cond.cond.op == ">") return lhs > rhs;
        if (cond.cond.op == ">=") return lhs >= rhs;
      }
      else
      {
        String lhs;
        if (!getStringValueFromArray_(array, row, lhs))
        {
          return false;
        }

        if (cond.cond.op == "IN")
        {
          for (const auto& rhs : cond.cond.values)
          {
            if (lhs == rhs) return true;
          }
          return false;
        }

        if (cond.cond.values.empty())
        {
          return false;
        }
        const String& rhs = cond.cond.values.front();
        if (cond.cond.op == "=") return lhs == rhs;
        if (cond.cond.op == "!=") return lhs != rhs;
        if (cond.cond.op == "<") return lhs < rhs;
        if (cond.cond.op == "<=") return lhs <= rhs;
        if (cond.cond.op == ">") return lhs > rhs;
        if (cond.cond.op == ">=") return lhs >= rhs;
      }

      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Unsupported filter operator", cond.cond.op);
    }

    bool evaluateFilterExpression_(const std::vector<ManualCondition>& conditions,
                                   const std::vector<String>& connectors,
                                   const std::shared_ptr<arrow::RecordBatch>& batch,
                                   int64_t row)
    {
      if (conditions.empty())
      {
        return true;
      }

      bool result = evaluateCondition_(conditions.front(), batch, row);
      for (Size i = 1; i < conditions.size(); ++i)
      {
        const String connector = (i - 1 < connectors.size()) ? connectors[i - 1] : "AND";
        const bool next = evaluateCondition_(conditions[i], batch, row);
        if (connector == "OR")
        {
          result = result || next;
        }
        else
        {
          result = result && next;
        }
      }
      return result;
    }

    std::shared_ptr<arrow::Table> applyManualFilter_(const std::shared_ptr<arrow::Table>& table,
                                                     const FilterExpression& parsed,
                                                     const String& filter)
    {
      std::vector<ManualCondition> conditions;
      conditions.reserve(parsed.conditions.size());

      auto schema = table->schema();
      for (const auto& cond : parsed.conditions)
      {
        ManualCondition manual;
        manual.cond = cond;
        manual.index = schema->GetFieldIndex(std::string(cond.column));
        if (manual.index < 0)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Column not found in parquet schema", cond.column);
        }
        if (cond.type == ColumnType::INT)
        {
          manual.int_values.reserve(cond.values.size());
          for (const auto& value : cond.values)
          {
            try
            {
              manual.int_values.push_back(value.toInt64());
            }
            catch (const std::exception&)
            {
              throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Invalid integer filter value", value);
            }
          }
        }
        conditions.push_back(std::move(manual));
      }

      arrow::TableBatchReader reader(*table);
      std::vector<std::shared_ptr<arrow::RecordBatch>> batches;
      std::shared_ptr<arrow::RecordBatch> batch;

      while (true)
      {
        auto read_status = reader.ReadNext(&batch);
        if (!read_status.ok())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Failed to read record batch", filter);
        }
        if (!batch)
        {
          break;
        }

        const int64_t num_rows = batch->num_rows();
        for (int64_t row = 0; row < num_rows; ++row)
        {
          if (evaluateFilterExpression_(conditions, parsed.connectors, batch, row))
          {
            batches.push_back(batch->Slice(row, 1));
          }
        }
      }

      if (batches.empty())
      {
        std::vector<std::shared_ptr<arrow::Array>> empty_columns;
        empty_columns.reserve(schema->num_fields());
        for (const auto& field : schema->fields())
        {
          auto array_result = arrow::MakeArrayOfNull(field->type(), 0);
          if (!array_result.ok())
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Failed to create empty column", array_result.status().ToString());
          }
          empty_columns.push_back(*array_result);
        }
        return arrow::Table::Make(schema, empty_columns);
      }

      auto table_result = arrow::Table::FromRecordBatches(schema, batches);
      if (!table_result.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to build filtered table", filter);
      }
      auto combined = (*table_result)->CombineChunks(arrow::default_memory_pool());
      if (!combined.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to combine filtered chunks", combined.status().ToString());
      }
      return *combined;
    }

#ifdef WITH_ARROW_DATASET
    /// Read parquet files via Arrow Dataset with predicate pushdown.
    std::shared_ptr<arrow::Table> readParquetTableDataset_(const std::vector<String>& filenames,
                                                           const FilterExpression& parsed,
                                                           const String& filter_context)
    {
      (void)filter_context;
      // Arrow Dataset enables predicate pushdown across multiple parquet files.
      auto filesystem = std::make_shared<arrow::fs::LocalFileSystem>();
      auto format = std::make_shared<arrow::dataset::ParquetFileFormat>();
      arrow::dataset::FileSystemFactoryOptions options;

      std::vector<std::string> paths;
      paths.reserve(filenames.size());
      for (const auto& filename : filenames)
      {
        paths.emplace_back(filename);
      }
      auto factory_result = arrow::dataset::FileSystemDatasetFactory::Make(
        filesystem, paths, format, options);
      if (!factory_result.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to create dataset factory", "dataset");
      }

      auto dataset_result = (*factory_result)->Finish();
      if (!dataset_result.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to create dataset", "dataset");
      }
      std::shared_ptr<arrow::dataset::Dataset> dataset = *dataset_result;

      auto builder = std::make_shared<arrow::dataset::ScannerBuilder>(dataset);
      if (!parsed.conditions.empty())
      {
        FilterExpression pruned = parsed;
        std::vector<String> dropped;
        normalizeAndPruneColumns_(pruned, dataset->schema(), dropped);
        if (!dropped.empty())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Unsupported filter columns after pruning", joinColumns_(dropped));
        }
        if (pruned.conditions.empty())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Filter expression became empty after pruning unsupported columns", filter_context);
        }
        // Convert the parsed filter into an Arrow expression for pushdown.
        arrow::compute::Expression expr = buildFilterExpression_(pruned);
        auto status = builder->Filter(expr);
        if (!status.ok())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Failed to bind dataset filter expression", status.ToString());
        }
      }

      auto scanner_result = builder->Finish();
      if (!scanner_result.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to create dataset scanner", "dataset");
      }
      auto table_result = (*scanner_result)->ToTable();
      if (!table_result.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to read parquet table via dataset scanner", "dataset");
      }

      auto combined = (*table_result)->CombineChunks(arrow::default_memory_pool());
      if (!combined.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to combine parquet chunks", combined.status().ToString());
      }
      return *combined;
    }
#endif

    /// Read parquet files via Arrow Dataset with predicate pushdown.
    std::shared_ptr<arrow::Table> applyArrowFilter_(const std::shared_ptr<arrow::Table>& table,
                                                    const FilterExpression& parsed,
                                                    const String& filter_context)
    {
      if (parsed.conditions.empty())
      {
        return table;
      }

      FilterExpression pruned = parsed;
      std::vector<String> dropped;
      normalizeAndPruneColumns_(pruned, table->schema(), dropped);
      if (!dropped.empty())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Unsupported filter columns after pruning", joinColumns_(dropped));
      }
      if (pruned.conditions.empty())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Filter expression became empty after pruning unsupported columns", filter_context);
      }

      // Execute the expression in memory when dataset pushdown is unavailable.
      arrow::compute::Expression expr = buildFilterExpression_(pruned);
      auto bound_result = expr.Bind(*table->schema());
      if (!bound_result.ok())
      {
        OPENMS_LOG_WARN << "Arrow compute filter unavailable, falling back to manual filter: "
                        << bound_result.status().ToString() << '\n';
        return applyManualFilter_(table, pruned, filter_context);
      }
      const arrow::compute::Expression& bound_expr = *bound_result;

      arrow::TableBatchReader reader(*table);
      std::vector<std::shared_ptr<arrow::RecordBatch>> batches;
      std::shared_ptr<arrow::RecordBatch> batch;

      while (true)
      {
        auto read_status = reader.ReadNext(&batch);
        if (!read_status.ok())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Failed to read record batch", filter_context);
        }
        if (!batch)
        {
          break;
        }

        auto mask_result = arrow::compute::ExecuteScalarExpression(
          bound_expr, *batch->schema(), arrow::Datum(batch));
        if (!mask_result.ok())
        {
          OPENMS_LOG_WARN << "Arrow compute filter failed, falling back to manual filter: "
                          << mask_result.status().ToString() << '\n';
          return applyManualFilter_(table, pruned, filter_context);
        }

        auto filtered_result = arrow::compute::Filter(arrow::Datum(batch), *mask_result);
        if (!filtered_result.ok())
        {
          OPENMS_LOG_WARN << "Arrow compute filter failed, falling back to manual filter: "
                          << filtered_result.status().ToString() << '\n';
          return applyManualFilter_(table, pruned, filter_context);
        }

        const auto& filtered_batch = filtered_result->record_batch();
        if (filtered_batch && filtered_batch->num_rows() > 0)
        {
          batches.push_back(filtered_batch);
        }
      }

      auto table_result = arrow::Table::FromRecordBatches(table->schema(), batches);
      if (!table_result.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to build filtered table", filter_context);
      }
      auto combined = (*table_result)->CombineChunks(arrow::default_memory_pool());
      if (!combined.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to combine filtered chunks", combined.status().ToString());
      }
      return *combined;
    }
  } // namespace

  XIMParquetFile::XIMParquetFile(const String& filename)
    : filename_(filename)
  {
    if (!File::exists(filename_))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename_);
    }
    filenames_.push_back(filename_);
  }

  XIMParquetFile::XIMParquetFile(const std::vector<String>& filenames)
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

  const String& XIMParquetFile::getFilename() const
  {
    return filename_;
  }

  const std::vector<String>& XIMParquetFile::getFilenames() const
  {
    return filenames_;
  }

  void XIMParquetFile::load(std::vector<XIMMobilogram>& output) const
  {
    getMobilograms(output);
  }

  void XIMParquetFile::getMobilograms_(std::vector<XIMMobilogram>& output,
                                         const FilterExpression& extra_filter,
                                         Int64 precursor_id,
                                         Int64 transition_id,
                                         const String& modified_sequence,
                                         Int64 precursor_charge,
                                         Int64 product_charge,
                                         Int64 ms_level,
                                         Int64 run_id,
                                         const String& mobilogram_type,
                                         Int64 feature_id,
                                         double feature_rt,
                                         const String& filter) const
  {
    output.clear();

    std::shared_ptr<arrow::Table> table;

    std::unordered_map<String, ColumnType> column_types = {
      {"RUN_ID", ColumnType::INT},
      {"MS_LEVEL", ColumnType::INT},
      {"FEATURE_ID", ColumnType::INT},
      {"PRECURSOR_ID", ColumnType::INT},
      {"TRANSITION_ID", ColumnType::INT},
      {"PRECURSOR_CHARGE", ColumnType::INT},
      {"PRODUCT_CHARGE", ColumnType::INT},
      {"DETECTING_TRANSITION", ColumnType::INT},
      {"PRECURSOR_DECOY", ColumnType::INT},
      {"PRODUCT_DECOY", ColumnType::INT},
      {"TRANSITION_ORDINAL", ColumnType::INT},
      {"MOBILITY_COMPRESSION", ColumnType::INT},
      {"INTENSITY_COMPRESSION", ColumnType::INT},
      {"SOURCE_FILE", ColumnType::STRING},
      {"MOBILOGRAM_TYPE", ColumnType::STRING},
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
                                                       run_id,
                                                       mobilogram_type,
                                                       feature_id,
                                                       feature_rt);
    FilterExpression string_filter = parseFilter_(filter, column_types);
    FilterExpression combined_filter = ParquetFilter::merge(args_filter, string_filter, "AND");
    combined_filter = ParquetFilter::merge(combined_filter, extra_filter, "AND");

    const String filter_context = filter.empty() ? String("<typed/args filter>") : filter;

    bool used_dataset = false;
#ifdef WITH_ARROW_DATASET
    if (!combined_filter.conditions.empty())
    {
      try
      {
        // Prefer dataset pushdown when available; fall back to compute filtering on failure.
        table = readParquetTableDataset_(filenames_, combined_filter, filter_context);
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
      table = readParquetTable_(filenames_);
    }
    if (!used_dataset)
    {
      table = applyArrowFilter_(table, combined_filter, filter_context);
    }

    auto run_id_col = getColumn_(table, "RUN_ID");
    auto source_file_col = getColumn_(table, "SOURCE_FILE", false);
    auto ms_level_col = getColumn_(table, "MS_LEVEL", false);
    auto mobilogram_type_col = getColumn_(table, "MOBILOGRAM_TYPE", false);
    auto precursor_id_col = getColumn_(table, "PRECURSOR_ID");
    auto transition_id_col = getColumn_(table, "TRANSITION_ID", false);
    auto feature_id_col = getColumn_(table, "FEATURE_ID", false);
    auto feature_rt_col = getColumn_(table, "FEATURE_RT", false);
    auto modified_sequence_col = getColumn_(table, "MODIFIED_SEQUENCE", false);
    auto precursor_charge_col = getColumn_(table, "PRECURSOR_CHARGE", false);
    auto product_charge_col = getColumn_(table, "PRODUCT_CHARGE", false);
    auto detecting_transition_col = getColumn_(table, "DETECTING_TRANSITION", false);
    auto precursor_decoy_col = getColumn_(table, "PRECURSOR_DECOY", false);
    auto product_decoy_col = getColumn_(table, "PRODUCT_DECOY", false);
    auto transition_ordinal_col = getColumn_(table, "TRANSITION_ORDINAL", false);
    auto transition_type_col = getColumn_(table, "TRANSITION_TYPE", false);
    auto annotation_col = getColumn_(table, "ANNOTATION", false);
    auto mobility_data_col = getColumn_(table, "MOBILITY_DATA");
    auto intensity_data_col = getColumn_(table, "INTENSITY_DATA");
    auto mobility_compression_col = getColumn_(table, "MOBILITY_COMPRESSION");
    auto intensity_compression_col = getColumn_(table, "INTENSITY_COMPRESSION");

    const int64_t rows = table->num_rows();
    output.reserve(rows);

    for (int64_t row = 0; row < rows; ++row)
    {
      XIMMobilogram mobilogram;

      // Defensive: RUN_ID is a required column but individual cells may be null in malformed files.
      // Use getOptionalInt_ and throw a descriptive error if a required cell is null.
      if (!getOptionalInt_(run_id_col, row, mobilogram.run_id))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "NULL value in required RUN_ID column", String(row));
      }
      mobilogram.has_precursor_id = getOptionalInt_(precursor_id_col, row, mobilogram.precursor_id);
      mobilogram.has_transition_id = getOptionalInt_(transition_id_col, row, mobilogram.transition_id);
      mobilogram.has_feature_id = getOptionalInt_(feature_id_col, row, mobilogram.feature_id);
      mobilogram.has_feature_rt = getOptionalDouble_(feature_rt_col, row, mobilogram.feature_rt);
      mobilogram.has_precursor_charge = getOptionalInt_(precursor_charge_col, row, mobilogram.precursor_charge);
      mobilogram.has_product_charge = getOptionalInt_(product_charge_col, row, mobilogram.product_charge);
      mobilogram.has_detecting_transition = getOptionalInt_(detecting_transition_col, row, mobilogram.detecting_transition);
      mobilogram.has_precursor_decoy = getOptionalInt_(precursor_decoy_col, row, mobilogram.precursor_decoy);
      mobilogram.has_product_decoy = getOptionalInt_(product_decoy_col, row, mobilogram.product_decoy);
      mobilogram.has_transition_ordinal = getOptionalInt_(transition_ordinal_col, row, mobilogram.transition_ordinal);

      getOptionalString_(source_file_col, row, mobilogram.source_file);
      getOptionalString_(mobilogram_type_col, row, mobilogram.mobilogram_type);
      getOptionalString_(modified_sequence_col, row, mobilogram.modified_sequence);
      getOptionalString_(transition_type_col, row, mobilogram.transition_type);
      getOptionalString_(annotation_col, row, mobilogram.annotation);

      Int64 ms_level_value = 0;
      if (getOptionalInt_(ms_level_col, row, ms_level_value))
      {
        mobilogram.ms_level = ms_level_value;
      }

      if (precursor_id >= 0 && (!mobilogram.has_precursor_id || mobilogram.precursor_id != precursor_id)) continue;
      if (transition_id >= 0 && (!mobilogram.has_transition_id || mobilogram.transition_id != transition_id)) continue;
      if (!modified_sequence.empty() && mobilogram.modified_sequence != modified_sequence) continue;
      if (precursor_charge >= 0 && (!mobilogram.has_precursor_charge || mobilogram.precursor_charge != precursor_charge)) continue;
      if (product_charge >= 0 && (!mobilogram.has_product_charge || mobilogram.product_charge != product_charge)) continue;
      if (ms_level >= 0 && mobilogram.ms_level != ms_level) continue;
      if (run_id >= 0 && mobilogram.run_id != run_id) continue;
      if (!mobilogram_type.empty() && mobilogram.mobilogram_type != mobilogram_type) continue;
      if (feature_id >= 0 && (!mobilogram.has_feature_id || mobilogram.feature_id != feature_id)) continue;
      if (feature_rt >= 0.0 && (!mobilogram.has_feature_rt || std::abs(mobilogram.feature_rt - feature_rt) > 1e-9)) continue;

      Int64 mobility_compression = 0;
      Int64 intensity_compression = 0;
      getOptionalInt_(mobility_compression_col, row, mobility_compression);
      getOptionalInt_(intensity_compression_col, row, intensity_compression);

      const String mobility_data = getBinaryView_(mobility_data_col, row);
      const String intensity_data = getBinaryView_(intensity_data_col, row);

      decodeBinary_(mobility_data, mobility_compression, mobilogram.mobility);
      decodeBinary_(intensity_data, intensity_compression, mobilogram.intensity);
      if (mobilogram.mobility.size() != mobilogram.intensity.size())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Mobility and intensity array size mismatch", String(row));
      }

      output.push_back(std::move(mobilogram));
    }
  }

  void XIMParquetFile::getRuns(std::vector<XIMRunInfo>& output) const
  {
    output.clear();

    const std::vector<String> columns = {"RUN_ID", "SOURCE_FILE"};
    auto table = readParquetTableColumns_(filenames_, columns);

    auto run_id_col = getColumn_(table, "RUN_ID");
    auto source_file_col = getColumn_(table, "SOURCE_FILE", false);

    std::unordered_set<std::string> seen;
    const int64_t rows = table->num_rows();
    output.reserve(rows);

    for (int64_t row = 0; row < rows; ++row)
    {
      XIMRunInfo info;
      if (!getOptionalInt_(run_id_col, row, info.run_id))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "NULL value in required RUN_ID column", String(row));
      }
      getOptionalString_(source_file_col, row, info.source_file);

      std::ostringstream oss;
      appendKeyPart_(oss, true, info.run_id);
      appendKeyPart_(oss, info.source_file);
      if (seen.insert(oss.str()).second)
      {
        output.push_back(std::move(info));
      }
    }
  }

  void XIMParquetFile::getMobilograms(std::vector<XIMMobilogram>& output,
                                        Int64 precursor_id,
                                        Int64 transition_id,
                                        const String& modified_sequence,
                                        Int64 precursor_charge,
                                        Int64 product_charge,
                                        Int64 ms_level,
                                        Int64 run_id,
                                        const String& mobilogram_type,
                                        Int64 feature_id,
                                        double feature_rt,
                                        const String& filter) const
  {
    FilterExpression empty_filter;
    getMobilograms_(output, empty_filter, precursor_id, transition_id, modified_sequence,
                      precursor_charge, product_charge, ms_level, run_id,
                      mobilogram_type, feature_id, feature_rt, filter);
  }

  void XIMParquetFile::getMobilograms(std::vector<XIMMobilogram>& output,
                                        const ParquetFilter& filter) const
  {
    getMobilograms_(output, filter.expression(), -1, -1, "", -1, -1, -1, -1, "", -1, -1.0, "");
  }

  void XIMParquetFile::getMobilograms(std::vector<XIMMobilogram>& output,
                                        const ParquetFilterBuilder& filter) const
  {
    getMobilograms(output, filter.filter());
  }

  void XIMParquetFile::getAnalytes(std::vector<XIMAnalyte>& output,
                                   const std::vector<String>& columns,
                                   bool nest_transitions) const
  {
    output.clear();

    const std::vector<String> default_columns = {
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
    const std::vector<String>& requested_columns = columns.empty() ? default_columns : columns;

  std::shared_ptr<arrow::Schema> schema = readParquetSchemaAllFiles_(filenames_);
    std::unordered_set<String> schema_columns;
    schema_columns.reserve(schema->num_fields());
    for (const auto& field : schema->fields())
    {
      schema_columns.insert(upper_(String(field->name())));
    }

    std::unordered_set<String> allowed_columns = {
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

    std::vector<String> normalized_columns;
    normalized_columns.reserve(requested_columns.size());
    for (const auto& name : requested_columns)
    {
      const String upper_name = upper_(name);
      if (allowed_columns.find(upper_name) == allowed_columns.end())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Unsupported analyte column", name);
      }
      if (schema_columns.find(upper_name) == schema_columns.end())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Column not found in parquet schema", name);
      }
      normalized_columns.push_back(upper_name);
    }

    std::unordered_set<String> requested_set;
    requested_set.reserve(normalized_columns.size());
    for (const auto& name : normalized_columns)
    {
      requested_set.insert(name);
    }

    if (nest_transitions)
    {
      // Ensure at least one precursor-level discriminator is requested when nesting transitions,
      // otherwise all rows will collapse into a single precursor key.
      static const std::unordered_set<String> precursor_discriminators = {
        "PRECURSOR_ID", "MODIFIED_SEQUENCE", "PRECURSOR_CHARGE", "PRECURSOR_DECOY"
      };
      bool has_precursor_discriminator = false;
      for (const auto& key : precursor_discriminators)
      {
        if (requested_set.find(key) != requested_set.end())
        {
          has_precursor_discriminator = true;
          break;
        }
      }
      if (!has_precursor_discriminator)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "nest_transitions=true requires at least one precursor discriminator column",
                                      "PRECURSOR_ID, MODIFIED_SEQUENCE, PRECURSOR_CHARGE, or PRECURSOR_DECOY");
      }
    }

    auto want = [&](const String& name) -> bool
    {
      return requested_set.find(name) != requested_set.end();
    };

    auto table = readParquetTableColumns_(filenames_, normalized_columns);

    auto precursor_id_col = getColumn_(table, "PRECURSOR_ID", want("PRECURSOR_ID"));
    auto modified_sequence_col = getColumn_(table, "MODIFIED_SEQUENCE", want("MODIFIED_SEQUENCE"));
    auto precursor_charge_col = getColumn_(table, "PRECURSOR_CHARGE", want("PRECURSOR_CHARGE"));
    auto precursor_decoy_col = getColumn_(table, "PRECURSOR_DECOY", want("PRECURSOR_DECOY"));
    auto transition_id_col = getColumn_(table, "TRANSITION_ID", want("TRANSITION_ID"));
    auto product_charge_col = getColumn_(table, "PRODUCT_CHARGE", want("PRODUCT_CHARGE"));
    auto transition_ordinal_col = getColumn_(table, "TRANSITION_ORDINAL", want("TRANSITION_ORDINAL"));
    auto detecting_transition_col = getColumn_(table, "DETECTING_TRANSITION", want("DETECTING_TRANSITION"));
    auto product_decoy_col = getColumn_(table, "PRODUCT_DECOY", want("PRODUCT_DECOY"));
    auto transition_type_col = getColumn_(table, "TRANSITION_TYPE", want("TRANSITION_TYPE"));
    auto annotation_col = getColumn_(table, "ANNOTATION", want("ANNOTATION"));

    const int64_t rows = table->num_rows();
    output.reserve(rows);

    if (!nest_transitions)
    {
      std::unordered_set<std::string> seen;
      for (int64_t row = 0; row < rows; ++row)
      {
        XIMAnalyte analyte;

        if (want("PRECURSOR_ID"))
        {
          analyte.has_precursor_id = getOptionalInt_(precursor_id_col, row, analyte.precursor_id);
        }
        if (want("MODIFIED_SEQUENCE"))
        {
          getOptionalString_(modified_sequence_col, row, analyte.modified_sequence);
        }
        if (want("PRECURSOR_CHARGE"))
        {
          analyte.has_precursor_charge = getOptionalInt_(precursor_charge_col, row, analyte.precursor_charge);
        }
        if (want("PRECURSOR_DECOY"))
        {
          analyte.has_precursor_decoy = getOptionalInt_(precursor_decoy_col, row, analyte.precursor_decoy);
        }

        if (want("TRANSITION_ID"))
        {
          analyte.has_transition_id = getOptionalInt_(transition_id_col, row, analyte.transition_id);
        }
        if (want("PRODUCT_CHARGE"))
        {
          analyte.has_product_charge = getOptionalInt_(product_charge_col, row, analyte.product_charge);
        }
        if (want("TRANSITION_ORDINAL"))
        {
          analyte.has_transition_ordinal = getOptionalInt_(transition_ordinal_col, row, analyte.transition_ordinal);
        }
        if (want("DETECTING_TRANSITION"))
        {
          analyte.has_detecting_transition = getOptionalInt_(detecting_transition_col, row, analyte.detecting_transition);
        }
        if (want("PRODUCT_DECOY"))
        {
          analyte.has_product_decoy = getOptionalInt_(product_decoy_col, row, analyte.product_decoy);
        }
        if (want("TRANSITION_TYPE"))
        {
          getOptionalString_(transition_type_col, row, analyte.transition_type);
        }
        if (want("ANNOTATION"))
        {
          getOptionalString_(annotation_col, row, analyte.annotation);
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
      XIMAnalyte analyte;
      analyte.has_precursor_id = getOptionalInt_(precursor_id_col, row, analyte.precursor_id);
      getOptionalString_(modified_sequence_col, row, analyte.modified_sequence);
      analyte.has_precursor_charge = getOptionalInt_(precursor_charge_col, row, analyte.precursor_charge);
      analyte.has_precursor_decoy = getOptionalInt_(precursor_decoy_col, row, analyte.precursor_decoy);

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
      String transition_type;
      String annotation;

      if (want("TRANSITION_ID"))
      {
        has_transition_id = getOptionalInt_(transition_id_col, row, transition_id);
      }
      if (want("PRODUCT_CHARGE"))
      {
        has_product_charge = getOptionalInt_(product_charge_col, row, product_charge);
      }
      if (want("TRANSITION_ORDINAL"))
      {
        has_transition_ordinal = getOptionalInt_(transition_ordinal_col, row, transition_ordinal);
      }
      if (want("DETECTING_TRANSITION"))
      {
        has_detecting_transition = getOptionalInt_(detecting_transition_col, row, detecting_transition);
      }
      if (want("PRODUCT_DECOY"))
      {
        has_product_decoy = getOptionalInt_(product_decoy_col, row, product_decoy);
      }
      if (want("TRANSITION_TYPE"))
      {
        getOptionalString_(transition_type_col, row, transition_type);
      }
      if (want("ANNOTATION"))
      {
        getOptionalString_(annotation_col, row, annotation);
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

  void XIMParquetFile::getColumns(std::vector<String>& output) const
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
