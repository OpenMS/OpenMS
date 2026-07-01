// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include "ParquetReaderHelpers.h"

#include <arrow/io/api.h>
#ifdef WITH_ARROW_DATASET
#include <arrow/dataset/api.h>
#include <arrow/filesystem/api.h>
#endif
#include <parquet/arrow/reader.h>

#include <cctype>
#include <exception>
#include <limits>
#include <sstream>
#include <utility>

namespace OpenMS
{
  namespace Internal
  {
    namespace ParquetReaderHelpers
    {
      namespace
      {
        std::unique_ptr<parquet::arrow::FileReader> openReader_(const std::string& filename)
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

          return std::move(*reader_result);
        }

        std::shared_ptr<arrow::Table> combineChunks_(const std::shared_ptr<arrow::Table>& table,
                                                     const std::string& filename,
                                                     const std::string& message)
        {
          auto combined = table->CombineChunks(arrow::default_memory_pool());
          if (!combined.ok())
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          message, filename);
          }

          return *combined;
        }

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

          auto concat_result = arrow::ConcatenateTables(tables);
          if (!concat_result.ok())
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Failed to concatenate parquet tables",
                                          concat_result.status().ToString());
          }

          auto combined = (*concat_result)->CombineChunks(arrow::default_memory_pool());
          if (!combined.ok())
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Failed to combine parquet chunks",
                                          combined.status().ToString());
          }

          return *combined;
        }

        std::shared_ptr<arrow::Table> readTableFromReader_(parquet::arrow::FileReader& reader,
                                                           const std::string& filename)
        {
          std::shared_ptr<arrow::Table> table;
          const auto status = reader.ReadTable(&table);
          if (!status.ok() || !table)
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Failed to read parquet table",
                                          filename + ": " + status.ToString());
          }
          return table;
        }

        std::shared_ptr<arrow::Table> readTableFromReader_(parquet::arrow::FileReader& reader,
                                                           const std::vector<int>& column_indices,
                                                           const std::string& filename)
        {
          std::shared_ptr<arrow::Table> table;
          const auto status = reader.ReadTable(column_indices, &table);
          if (!status.ok() || !table)
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Failed to read parquet table",
                                          filename + ": " + status.ToString());
          }
          return table;
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

        std::vector<std::string> tokenize_(const std::string& expr)
        {
          std::vector<std::string> tokens;
          std::string current;
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
              tokens.emplace_back(StringUtils::toStr(c));
              continue;
            }

            if (c == '=' || c == '!' || c == '<' || c == '>')
            {
              if (!current.empty())
              {
                tokens.push_back(current);
                current.clear();
              }
              std::string op;
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
              std::string op;
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

        std::unordered_map<std::string, std::string> buildSchemaNameMap_(const std::shared_ptr<arrow::Schema>& schema)
        {
          std::unordered_map<std::string, std::string> name_map;
          if (!schema)
          {
            return name_map;
          }

          name_map.reserve(schema->num_fields());
          for (const auto& field : schema->fields())
          {
            name_map[upper_(std::string(field->name()))] = std::string(field->name());
          }

          return name_map;
        }

        struct PrunedFilterColumns
        {
          std::vector<std::string> unsupported_columns;
          std::vector<std::string> unknown_columns;
        };

        void normalizeAndPruneColumns_(FilterExpression& expr,
                                       const std::shared_ptr<arrow::Schema>& schema,
                                       const std::unordered_set<std::string>& unsupported_columns,
                                       PrunedFilterColumns& pruned_columns)
        {
          pruned_columns.unsupported_columns.clear();
          pruned_columns.unknown_columns.clear();
          if (expr.conditions.empty())
          {
            return;
          }

          const auto name_map = buildSchemaNameMap_(schema);

          std::vector<Condition> kept;
          kept.reserve(expr.conditions.size());
          std::vector<std::string> new_connectors;
          new_connectors.reserve(expr.connectors.size());

          for (Size i = 0; i < expr.conditions.size(); ++i)
          {
            Condition cond = expr.conditions[i];
            const std::string original = cond.column;
            const std::string upper_column = upper_(cond.column);

            auto it = name_map.find(upper_column);
            if (it == name_map.end())
            {
              pruned_columns.unknown_columns.push_back(original);
              continue;
            }
            cond.column = it->second;

            if (unsupported_columns.contains(upper_column))
            {
              pruned_columns.unsupported_columns.push_back(original);
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
            kept.push_back(std::move(cond));
          }

          expr.conditions.swap(kept);
          expr.connectors.swap(new_connectors);
        }

        std::string joinColumns_(const std::vector<std::string>& columns)
        {
          std::ostringstream oss;
          for (Size i = 0; i < columns.size(); ++i)
          {
            if (i > 0)
            {
              oss << ", ";
            }
            oss << columns[i];
          }
          return oss.str();
        }

        bool normalizeAndValidateFilter_(FilterExpression& expr,
                                         const std::shared_ptr<arrow::Schema>& schema,
                                         const std::string& filter_context,
                                         const FilterPruningOptions& pruning_options)
        {
          PrunedFilterColumns pruned_columns;
          normalizeAndPruneColumns_(expr, schema, pruning_options.unsupported_columns, pruned_columns);

          if (!pruned_columns.unknown_columns.empty())
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Unknown filter columns",
                                          joinColumns_(pruned_columns.unknown_columns));
          }

          if (!pruned_columns.unsupported_columns.empty())
          {
            const std::string dropped_text = joinColumns_(pruned_columns.unsupported_columns);
            if (pruning_options.reject_unsupported_after_pruning)
            {
              throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Unsupported filter columns after pruning", dropped_text);
            }
            OPENMS_LOG_WARN << "Dropping unsupported filter columns: "
                            << dropped_text << '\n';
          }

          if (expr.conditions.empty())
          {
            if (pruning_options.reject_empty_after_pruning)
            {
              throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Filter expression became empty after pruning unsupported columns",
                                            filter_context);
            }
            return false;
          }

          return true;
        }

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
                auto status = builder.Append(StringUtils::toInt64(value));
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

          if (cond.values.empty())
          {
            return arrow::compute::literal(false);
          }

          arrow::compute::Expression literal;
          if (cond.type == ColumnType::INT)
          {
            literal = arrow::compute::literal(static_cast<int64_t>(StringUtils::toInt64(cond.values.front())));
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

        arrow::compute::Expression buildFilterExpression_(const FilterExpression& expr)
        {
          if (expr.conditions.empty())
          {
            return arrow::compute::literal(true);
          }

          std::vector<arrow::compute::Expression> condition_expressions;
          condition_expressions.reserve(expr.conditions.size());
          for (const auto& cond : expr.conditions)
          {
            condition_expressions.push_back(buildConditionExpression_(cond));
          }

          arrow::compute::Expression result = condition_expressions.front();
          for (Size i = 1; i < condition_expressions.size(); ++i)
          {
            if (i - 1 < expr.connectors.size() && expr.connectors[i - 1] == "OR")
            {
              result = arrow::compute::or_(result, condition_expressions[i]);
            }
            else
            {
              result = arrow::compute::and_(result, condition_expressions[i]);
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
            {
              const uint64_t raw_value = static_cast<arrow::UInt64Array*>(array.get())->Value(row);
              if (raw_value > static_cast<uint64_t>(std::numeric_limits<Int64>::max()))
              {
                throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                              "UINT64 value exceeds Int64 range",
                                              StringUtils::toStr(raw_value));
              }
              value = static_cast<Int64>(raw_value);
              return true;
            }
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
                                            "Unsupported integer column type",
                                            std::string(array->type()->ToString()));
          }
        }

        bool getStringValueFromArray_(const std::shared_ptr<arrow::Array>& array, int64_t row, std::string& value)
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
                                            "Unsupported string column type",
                                            std::string(array->type()->ToString()));
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
                if (lhs == rhs)
                {
                  return true;
                }
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
            std::string lhs;
            if (!getStringValueFromArray_(array, row, lhs))
            {
              return false;
            }

            if (cond.cond.op == "IN")
            {
              for (const auto& rhs : cond.cond.values)
              {
                if (lhs == rhs)
                {
                  return true;
                }
              }
              return false;
            }

            if (cond.cond.values.empty())
            {
              return false;
            }

            const std::string& rhs = cond.cond.values.front();
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
                                       const std::vector<std::string>& connectors,
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
            const std::string connector = (i - 1 < connectors.size()) ? connectors[i - 1] : "AND";
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
      } // namespace

      std::shared_ptr<arrow::Table> readParquetTable(const std::string& filename)
      {
        auto reader = openReader_(filename);
        return combineChunks_(readTableFromReader_(*reader, filename), filename,
                              "Failed to combine parquet chunks");
      }

      std::shared_ptr<arrow::Table> readParquetTable(const std::vector<std::string>& filenames)
      {
        std::vector<std::shared_ptr<arrow::Table>> tables;
        tables.reserve(filenames.size());
        for (const auto& filename : filenames)
        {
          tables.push_back(readParquetTable(filename));
        }
        return concatenateTables_(tables);
      }

      std::shared_ptr<arrow::Schema> readParquetSchema(const std::string& filename)
      {
        auto reader = openReader_(filename);

        std::shared_ptr<arrow::Schema> schema;
        auto status = reader->GetSchema(&schema);
        if (!status.ok() || !schema)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Failed to read parquet schema", filename);
        }
        return schema;
      }

      std::shared_ptr<arrow::Schema> readParquetSchemaAllFiles(const std::vector<std::string>& filenames)
      {
        if (filenames.empty())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "No parquet files provided for schema validation", "");
        }

        auto base = readParquetSchema(filenames.front());
        for (size_t i = 1; i < filenames.size(); ++i)
        {
          auto other = readParquetSchema(filenames[i]);
          if (!other->Equals(*base))
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Parquet input files have incompatible schemas: '" +
                                            filenames.front() + "' vs '" + filenames[i] + "'",
                                          "");
          }
        }
        return base;
      }

      std::shared_ptr<arrow::Table> readParquetTableColumns(const std::string& filename,
                                                            const std::vector<std::string>& columns)
      {
        auto reader = openReader_(filename);

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
                                              "No requested columns found in parquet file '" + filename + "'");
        }

        return combineChunks_(readTableFromReader_(*reader, column_indices, filename), filename,
                              "Failed to combine parquet chunks");
      }

      std::shared_ptr<arrow::Table> readParquetTableColumns(const std::vector<std::string>& filenames,
                                                            const std::vector<std::string>& columns)
      {
        std::vector<std::shared_ptr<arrow::Table>> tables;
        tables.reserve(filenames.size());
        for (const auto& filename : filenames)
        {
          tables.push_back(readParquetTableColumns(filename, columns));
        }
        return concatenateTables_(tables);
      }

      std::shared_ptr<arrow::Array> getColumn(const std::shared_ptr<arrow::Table>& table,
                                              const std::string& name,
                                              bool required)
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

      bool getOptionalInt(const std::shared_ptr<arrow::Array>& array, int64_t row, Int64& value)
      {
        return getIntValueFromArray_(array, row, value);
      }

      bool getOptionalDouble(const std::shared_ptr<arrow::Array>& array, int64_t row, double& value)
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
                                          "Unsupported floating-point column type",
                                          std::string(array->type()->ToString()));
        }
      }

      bool getOptionalString(const std::shared_ptr<arrow::Array>& array, int64_t row, std::string& value)
      {
        return getStringValueFromArray_(array, row, value);
      }

      std::string getBinaryView(const std::shared_ptr<arrow::Array>& array, int64_t row)
      {
        if (!array || array->IsNull(row))
        {
          return std::string();
        }

        switch (array->type_id())
        {
          case arrow::Type::BINARY:
          {
            auto typed = static_cast<arrow::BinaryArray*>(array.get());
            const auto view = typed->GetView(row);
            const size_t view_size = view.size();
            if (view_size > static_cast<size_t>(std::numeric_limits<Size>::max()))
            {
              throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Binary data too large", StringUtils::toStr(view_size));
            }
            return std::string(view.data(), view_size);
          }
          case arrow::Type::LARGE_BINARY:
          {
            auto typed = static_cast<arrow::LargeBinaryArray*>(array.get());
            const auto view = typed->GetView(row);
            const size_t view_size = view.size();
            if (view_size > static_cast<size_t>(std::numeric_limits<Size>::max()))
            {
              throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Binary data too large", StringUtils::toStr(view_size));
            }
            return std::string(view.data(), view_size);
          }
          default:
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Unsupported binary column type",
                                          std::string(array->type()->ToString()));
        }
      }

      FilterExpression parseFilter(const std::string& filter,
                                   const std::unordered_map<std::string, ColumnType>& column_types)
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
          std::string column = upper_(tokens[i++]);
          if (i >= tokens.size())
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Missing filter operator for column", column);
          }

          std::string op = upper_(tokens[i++]);
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
            std::string close_token;
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

          expr.conditions.push_back(std::move(cond));
          if (i < tokens.size())
          {
            std::string connector = upper_(tokens[i]);
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

#ifdef WITH_ARROW_DATASET
      std::shared_ptr<arrow::Table> readParquetTableDataset(const std::vector<std::string>& filenames,
                                                            const FilterExpression& parsed,
                                                            const std::string& filter_context,
                                                            const FilterPruningOptions& pruning_options)
      {
        auto filesystem = std::make_shared<arrow::fs::LocalFileSystem>();
        auto format = std::make_shared<arrow::dataset::ParquetFileFormat>();
        arrow::dataset::FileSystemFactoryOptions options;

        std::vector<std::string> paths;
        paths.reserve(filenames.size());
        for (const auto& filename : filenames)
        {
          paths.emplace_back(filename);
        }

        auto factory_result = arrow::dataset::FileSystemDatasetFactory::Make(filesystem, paths, format, options);
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
          if (!normalizeAndValidateFilter_(pruned, dataset->schema(), filter_context, pruning_options))
          {
            return readParquetTable(filenames);
          }

          auto status = builder->Filter(buildFilterExpression_(pruned));
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

      std::shared_ptr<arrow::Table> applyManualFilter(const std::shared_ptr<arrow::Table>& table,
                                                      const FilterExpression& parsed,
                                                      const std::string& filter_context)
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
                manual.int_values.push_back(StringUtils::toInt64(value));
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
                                          "Failed to read record batch", filter_context);
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

      std::shared_ptr<arrow::Table> applyArrowFilter(const std::shared_ptr<arrow::Table>& table,
                                                     const FilterExpression& parsed,
                                                     const std::string& filter_context,
                                                     const FilterPruningOptions& pruning_options)
      {
        if (parsed.conditions.empty())
        {
          return table;
        }

        FilterExpression pruned = parsed;
        if (!normalizeAndValidateFilter_(pruned, table->schema(), filter_context, pruning_options))
        {
          return table;
        }

        auto bound_result = buildFilterExpression_(pruned).Bind(*table->schema());
        if (!bound_result.ok())
        {
          OPENMS_LOG_WARN << "Arrow compute filter unavailable, falling back to manual filter: "
                          << bound_result.status().ToString() << '\n';
          return applyManualFilter(table, pruned, filter_context);
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
            return applyManualFilter(table, pruned, filter_context);
          }

          auto filtered_result = arrow::compute::Filter(arrow::Datum(batch), *mask_result);
          if (!filtered_result.ok())
          {
            OPENMS_LOG_WARN << "Arrow compute filter failed, falling back to manual filter: "
                            << filtered_result.status().ToString() << '\n';
            return applyManualFilter(table, pruned, filter_context);
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
    } // namespace ParquetReaderHelpers
  } // namespace Internal
} // namespace OpenMS
