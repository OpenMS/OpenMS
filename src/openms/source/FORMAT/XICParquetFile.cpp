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
#include <arrow/compute/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/reader.h>
#endif
#ifdef WITH_ARROW_DATASET
#include <arrow/dataset/api.h>
#include <arrow/filesystem/api.h>
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
      return arrow::compute::literal(false);
    }

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

#ifdef WITH_ARROW_DATASET
    std::shared_ptr<arrow::Table> readParquetTableDataset_(const String& filename,
                                                           const String& filter,
                                                           const std::unordered_map<String, ColumnType>& column_types)
    {
      auto filesystem = std::make_shared<arrow::fs::LocalFileSystem>();
      auto format = std::make_shared<arrow::dataset::ParquetFileFormat>();
      arrow::dataset::FileSystemFactoryOptions options;

      auto factory_result = arrow::dataset::FileSystemDatasetFactory::Make(
        filesystem, std::vector<std::string>{std::string(filename)}, format, options);
      if (!factory_result.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to create dataset factory", filename);
      }

      auto dataset_result = (*factory_result)->Finish();
      if (!dataset_result.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to create dataset", filename);
      }
      std::shared_ptr<arrow::dataset::Dataset> dataset = *dataset_result;

      auto builder = std::make_shared<arrow::dataset::ScannerBuilder>(dataset);
      if (!filter.empty())
      {
        FilterExpression parsed = parseFilter_(filter, column_types);
        arrow::compute::Expression expr = buildFilterExpression_(parsed);
        auto status = builder->Filter(expr);
        if (!status.ok())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Failed to bind dataset filter expression", filter);
        }
      }

      auto scanner_result = builder->Finish();
      if (!scanner_result.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to create dataset scanner", filename);
      }
      auto table_result = (*scanner_result)->ToTable();
      if (!table_result.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to read parquet table via dataset scanner", filename);
      }

      auto combined = (*table_result)->CombineChunks(arrow::default_memory_pool());
      if (!combined.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to combine parquet chunks", filename);
      }
      return *combined;
    }
#endif

    std::shared_ptr<arrow::Table> applyArrowFilter_(const std::shared_ptr<arrow::Table>& table,
                                                    const String& filter,
                                                    const std::unordered_map<String, ColumnType>& column_types)
    {
      if (filter.empty())
      {
        return table;
      }

      FilterExpression parsed = parseFilter_(filter, column_types);
      if (parsed.conditions.empty())
      {
        return table;
      }

      arrow::compute::Expression expr = buildFilterExpression_(parsed);
      auto bound_result = expr.Bind(*table->schema());
      if (!bound_result.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to bind filter expression", filter);
      }
      arrow::compute::Expression bound_expr = *bound_result;

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

        auto mask_result = arrow::compute::ExecuteScalarExpression(
          bound_expr, *batch->schema(), arrow::Datum(batch));
        if (!mask_result.ok())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Failed to evaluate filter expression", filter);
        }

        auto filtered_result = arrow::compute::Filter(arrow::Datum(batch), *mask_result);
        if (!filtered_result.ok())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Failed to apply filter expression", filter);
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
                                      "Failed to build filtered table", filter);
      }
      return *table_result;
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

    bool used_dataset = false;
#ifdef WITH_ARROW_DATASET
    if (!filter.empty())
    {
      table = readParquetTableDataset_(filename_, filter, column_types);
      used_dataset = true;
    }
#endif
    if (!used_dataset)
    {
      table = applyArrowFilter_(table, filter, column_types);
    }

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

      if (precursor_id >= 0 && (!chrom.has_precursor_id || chrom.precursor_id != precursor_id)) continue;
      if (transition_id >= 0 && (!chrom.has_transition_id || chrom.transition_id != transition_id)) continue;
      if (!modified_sequence.empty() && chrom.modified_sequence != modified_sequence) continue;
      if (precursor_charge >= 0 && (!chrom.has_precursor_charge || chrom.precursor_charge != precursor_charge)) continue;
      if (product_charge >= 0 && (!chrom.has_product_charge || chrom.product_charge != product_charge)) continue;
      if (ms_level >= 0 && chrom.ms_level != ms_level) continue;
      if (run_id >= 0 && chrom.run_id != run_id) continue;

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
