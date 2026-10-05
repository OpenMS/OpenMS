// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $
#pragma once

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/IdentificationDataFile.h>
#include <arrow/api.h>
#include <arrow/io/file.h>
#include <filesystem>
#include <limits>
#include <map>
#include <nlohmann/json.hpp>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>
#include <tuple>
#include <type_traits>

namespace OpenMS::Internal::IdentificationDataIO
{
using ID = IdentificationData;
using Json = nlohmann::ordered_json;
using Options = IdentificationDataFile::Options;
[[noreturn]] void invalid(const std::string& message);
void check(const arrow::Status& status);
template<class T>
T value(arrow::Result<T> result)
{
  if (! result.ok()) invalid(result.status().ToString());
  return std::move(result).ValueOrDie();
}
/// JSON numeric conversions otherwise silently narrow or accept fractional values.
template<class T>
T integer(const Json& input)
{
  static_assert(std::is_integral_v<T>);
  if (input.is_number_unsigned())
  {
    const auto v = input.get<uint64_t>();
    if (v > static_cast<uint64_t>(std::numeric_limits<T>::max())) invalid("Integer outside the declared range");
    return static_cast<T>(v);
  }
  if (input.is_number_integer())
  {
    const auto v = input.get<int64_t>();
    if constexpr (std::is_unsigned_v<T>)
    {
      if (v < 0 || static_cast<uint64_t>(v) > std::numeric_limits<T>::max()) invalid("Integer outside the declared range");
    }
    else if (v < std::numeric_limits<T>::min() || v > std::numeric_limits<T>::max())
      invalid("Integer outside the declared range");
    return static_cast<T>(v);
  }
  invalid("Expected an integer, not a floating-point or other JSON value");
}
void validateText(const std::string& text);
void validateOptions(const Options& options);
std::filesystem::path tablePath(const std::filesystem::path& root, const std::string& relative);

template<class B, class T>
void append(arrow::ArrayBuilder& builder, const T& item)
{ check(static_cast<B&>(builder).Append(item)); }
void appendText(arrow::ArrayBuilder& builder, const std::string& text);
template<class A>
typename A::value_type number(const arrow::Array& array, int64_t row)
{
  if (array.IsNull(row)) invalid("Unexpected null in required column");
  return static_cast<const A&>(array).Value(row);
}
std::string text(const arrow::Array& array, int64_t row);
template<class A>
std::optional<typename A::value_type> optionalNumber(const arrow::Array& array, int64_t row)
{
  if (array.IsNull(row)) return std::nullopt;
  return number<A>(array, row);
}

// Descriptor ID identifies name, DataValue type and unit; row cells hold typed values only.
class Dictionary
{
public:
  UInt32 add(const std::string& name, const DataValue& item);
  UInt32 find(const std::string& name, const DataValue& item) const;
  void collect(const MetaInfoInterface& metadata);
  Json toJson() const;
  void load(const Json& descriptors);
  struct Descriptor
  {
    std::string name;
    DataValue::DataType type;
    DataValue::UnitType unit_type;
    Int32 unit;
  };
  const Descriptor& at(UInt32 index) const;

private:
  using Key = std::tuple<std::string, unsigned, unsigned, Int32>;
  std::map<Key, UInt32> ids_;
  std::vector<Descriptor> values_;
};
std::shared_ptr<arrow::DataType> metadataType();
void appendMetadata(arrow::ArrayBuilder& builder, const MetaInfoInterface& metadata, const Dictionary& dictionary);
void readMetadata(const arrow::Array& array, int64_t row, MetaInfoInterface& metadata, const Dictionary& dictionary);
Size metadataBytes(const MetaInfoInterface& metadata);
Json metadataJson(const MetaInfoInterface& metadata);
void readMetadataJson(const Json& json, MetaInfoInterface& metadata);
Json scoreJson(const ID::ScoreDefinition& score);
ID::ScoreDefinition readScoreJson(const Json& json);

// Record builders are reused until either the byte target or row limit is reached.
class TableWriter
{
public:
  TableWriter(const std::filesystem::path& path, std::shared_ptr<arrow::Schema> schema, const Options& options);
  ~TableWriter();
  arrow::ArrayBuilder& column(Size index)
  { return *builders_.at(index); }
  void finishRow(Size bytes);
  void close();
  UInt64 rows() const
  { return total_; }

private:
  void flush_();
  std::shared_ptr<arrow::Schema> schema_;
  Options options_;
  std::shared_ptr<arrow::io::FileOutputStream> sink_;
  std::unique_ptr<parquet::arrow::FileWriter> writer_;
  std::vector<std::unique_ptr<arrow::ArrayBuilder>> builders_;
  Size buffered_ = 0;
  Size bytes_ = 0;
  UInt64 total_ = 0;
};

// Keeps one Arrow batch and a current-row cursor. Projection uses top-level field names.
class TableReader
{
public:
  TableReader(const std::filesystem::path& path,
              const std::shared_ptr<arrow::Schema>& expected,
              const Options& options,
              const std::vector<std::string>& columns = {});
  bool next();
  const arrow::Array& column(const std::string& name) const;
  const arrow::Array& column(Size index) const
  { return *batch_->column(static_cast<int>(index)); }
  int64_t row() const
  { return row_; }
  UInt64 rows() const
  { return total_rows_; }
  bool hasColumn(const std::string& name) const;

private:
  std::shared_ptr<arrow::io::ReadableFile> source_;
  std::unique_ptr<parquet::arrow::FileReader> reader_;
  std::shared_ptr<arrow::RecordBatchReader> batches_;
  std::shared_ptr<arrow::RecordBatch> batch_;
  int64_t row_ = -1;
  UInt64 total_rows_ = 0;
};

// Configuration JSON excludes hits/groups. Per-record protein values are typed tables.
Json processingJson(const ProteinIdentification& processing);
ProteinIdentification readProcessingJson(const Json& json);
Json writeInference(const std::filesystem::path& directory, const ID::InferenceResult& result, const Options& options);
ID::InferenceResult readInference(const std::filesystem::path& directory, const Json& descriptor, const Options& options);
void validateInferenceTables(const std::filesystem::path& directory,
                             const Json& descriptor,
                             const Options& options,
                             const std::map<std::string, UInt64>& counters = {});
void validateParents(const std::filesystem::path& path, const Dictionary& dictionary, const Options& options);
void writeParents(const std::filesystem::path& path, const std::vector<ID::ParentRecord>& parents, Dictionary& dictionary, const Options& options);
std::vector<ID::ParentRecord> readParents(const std::filesystem::path& path, const Dictionary& dictionary, const Options& options);
} // namespace OpenMS::Internal::IdentificationDataIO
