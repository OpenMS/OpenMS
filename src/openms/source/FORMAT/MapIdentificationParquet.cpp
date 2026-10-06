// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
#include "MapIdentificationParquet.h"

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/CONCEPT/UniqueIdInterface.h>
#include <OpenMS/FORMAT/IdentificationDataFile.h>

#include <arrow/api.h>
#include <arrow/io/file.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>

#include <chrono>
#include <map>
#include <optional>
#include <set>
#include <system_error>
#include <thread>

namespace OpenMS::Internal::MapIdentificationParquet
{
namespace
{
  namespace fs = std::filesystem;
  using ID = IdentificationData;

  [[noreturn]] void invalid(const std::string& message)
  { throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message, "identification links"); }
  [[noreturn]] void fileError(const fs::path& path, const std::string& message)
  { throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path.string(), message); }
  void check(const arrow::Status& status)
  {
    if (! status.ok()) invalid(status.ToString());
  }
  template<class T>
  T value(arrow::Result<T> result)
  {
    check(result.status());
    return std::move(result).ValueUnsafe();
  }

  // One row per link; "primary" rows carry the molecular identity, the others a record reference.
  const std::string PRIMARY = "primary", QUERY = "query", MATCH = "match";
  // Run UUIDs are dictionary-encoded strings, as in every table that carries them: one entry per run.
  std::shared_ptr<arrow::DataType> runUuidType()
  { return arrow::dictionary(arrow::int32(), arrow::utf8()); }
  std::shared_ptr<arrow::Schema> linkSchema(const std::shared_ptr<arrow::DataType>& run_uuid_type = runUuidType())
  {
    return arrow::schema({arrow::field("feature_unique_id", arrow::uint64(), false), arrow::field("link", arrow::utf8(), false),
                          arrow::field("run_uuid", run_uuid_type, true), arrow::field("record_id", arrow::uint64(), true),
                          arrow::field("encoding", arrow::uint8(), true), arrow::field("representation", arrow::utf8(), true)});
  }

  /// Run UUID column of a link batch; tables written by other tools may store it as a plain string.
  class RunUuids
  {
  public:
    explicit RunUuids(const std::shared_ptr<arrow::Array>& array)
    {
      if (array->type_id() == arrow::Type::DICTIONARY)
      {
        dictionary_ = std::static_pointer_cast<arrow::DictionaryArray>(array);
        values_ = std::static_pointer_cast<arrow::StringArray>(dictionary_->dictionary());
      }
      else
        values_ = std::static_pointer_cast<arrow::StringArray>(array);
    }
    bool isNull(int64_t row) const
    { return dictionary_ ? dictionary_->IsNull(row) : values_->IsNull(row); }
    std::string at(int64_t row) const
    { return dictionary_ ? values_->GetString(dictionary_->GetValueIndex(row)) : values_->GetString(row); }

  private:
    std::shared_ptr<arrow::DictionaryArray> dictionary_;
    std::shared_ptr<arrow::StringArray> values_;
  };

  std::error_code renameWithRetry(const fs::path& from, const fs::path& to)
  {
    // On Windows, a scanner or indexer may still hold a handle for a moment.
    std::error_code ec;
    for (int attempt = 0; attempt < 5; ++attempt)
    {
      fs::rename(from, to, ec);
      if (! ec) return ec;
      std::this_thread::sleep_for(std::chrono::milliseconds(50 * (attempt + 1)));
    }
    return ec;
  }
} // namespace

bool hasNativeIdentifications(const IdentificationData& data, const std::vector<const BaseFeature*>& features)
{
  if (! data.empty()) return true;
  for (const auto* feature : features)
    if (feature->hasPrimaryID() || ! feature->getIDMatches().empty() || ! feature->getIDQueries().empty()) return true;
  return false;
}

void validateLinks(const IdentificationData& data, const std::vector<const BaseFeature*>& features)
{
  for (const auto* feature : features)
  {
    for (const auto& reference : feature->getIDQueries())
    {
      const auto* run = data.findRunByUuid(reference.run_uuid);
      if (! run || ! run->findIdentification(reference.query))
        invalid("Feature " + std::to_string(feature->getUniqueId()) + " links to a query that does not exist in the map's identification data");
    }
    for (const auto& reference : feature->getIDMatches())
    {
      const auto* run = data.findRunByUuid(reference.run_uuid);
      if (! run || ! run->findMatch(reference.match))
        invalid("Feature " + std::to_string(feature->getUniqueId()) + " links to a match that does not exist in the map's identification data");
    }
  }
}

void store(const fs::path& directory, const IdentificationData& data, const std::vector<const BaseFeature*>& features)
{
  if (! hasNativeIdentifications(data, features)) return;
  validateLinks(data, features);
  arrow::UInt64Builder feature_ids, record_ids;
  arrow::StringBuilder links, representations;
  arrow::StringDictionary32Builder uuids;
  arrow::UInt8Builder encodings;
  std::set<UInt64> linked;
  for (const auto* feature : features)
  {
    const bool has_links = feature->hasPrimaryID() || ! feature->getIDQueries().empty() || ! feature->getIDMatches().empty();
    if (! has_links) continue;
    const UInt64 unique_id = feature->getUniqueId();
    // Links are keyed by unique ID, like psms.parquet, so linked features need distinct valid IDs.
    if (! UniqueIdInterface::isValid(unique_id)) invalid("A feature with identification links has no valid unique ID");
    if (! linked.insert(unique_id).second) invalid("Features with identification links share unique ID " + std::to_string(unique_id));
    const auto add = [&](const std::string& link, const std::string* uuid, std::optional<UInt64> record, std::optional<ID::MoleculeIdentity> identity) {
      check(feature_ids.Append(unique_id));
      check(links.Append(link));
      check(uuid ? uuids.Append(*uuid) : uuids.AppendNull());
      check(record ? record_ids.Append(*record) : record_ids.AppendNull());
      check(identity ? encodings.Append(static_cast<uint8_t>(identity->encoding)) : encodings.AppendNull());
      check(identity ? representations.Append(identity->representation) : representations.AppendNull());
    };
    if (feature->hasPrimaryID()) add(PRIMARY, nullptr, std::nullopt, feature->getPrimaryID());
    for (const auto& reference : feature->getIDQueries())
      add(QUERY, &reference.run_uuid, reference.query.value, std::nullopt);
    for (const auto& reference : feature->getIDMatches())
      add(MATCH, &reference.run_uuid, reference.match.value, std::nullopt);
  }
  std::vector<std::shared_ptr<arrow::Array>> columns(6);
  check(feature_ids.Finish(&columns[0]));
  check(links.Finish(&columns[1]));
  check(uuids.Finish(&columns[2]));
  check(record_ids.Finish(&columns[3]));
  check(encodings.Finish(&columns[4]));
  check(representations.Finish(&columns[5]));
  const auto table = arrow::Table::Make(linkSchema(), columns);
  check(table->ValidateFull());
  const auto properties = parquet::WriterProperties::Builder().compression(parquet::Compression::ZSTD)->build();
  // Store the Arrow schema so that readers get the dictionary type back.
  const auto arrow_properties = parquet::ArrowWriterProperties::Builder().store_schema()->build();
  const auto file = directory / LINKS_FILE;
  auto sink = value(arrow::io::FileOutputStream::Open(file.string()));
  check(parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), sink, std::max<int64_t>(1, table->num_rows()), properties,
                                   arrow_properties));
  check(sink->Close());
  // The nested bundle is written even when empty so that its presence states the map's model.
  IdentificationDataFile::store((directory / IDENTIFICATIONS_DIRECTORY).string(), data);
}

void load(const fs::path& directory, IdentificationData& data, const std::vector<BaseFeature*>& features)
{
  std::error_code ec;
  const auto nested = directory / IDENTIFICATIONS_DIRECTORY;
  const auto links_file = directory / LINKS_FILE;
  const bool has_data = fs::exists(nested, ec);
  const bool has_links = fs::exists(links_file, ec);
  if (! has_data && ! has_links) return;
  if (! has_data) invalid("Identification links exist without the map's identification data");
  IdentificationData loaded;
  IdentificationDataFile::load(nested.string(), loaded);

  // Decode the links first; apply them only after everything validated.
  struct Links
  {
    std::optional<ID::MoleculeIdentity> primary;
    std::set<ID::QueryReference> queries;
    std::set<ID::MatchReference> matches;
  };
  std::map<UInt64, Links> by_feature;
  if (has_links)
  {
    auto input = value(arrow::io::ReadableFile::Open(links_file.string()));
    auto reader = value(parquet::arrow::OpenFile(input, arrow::default_memory_pool()));
    std::shared_ptr<arrow::Table> table;
#if ARROW_VERSION_MAJOR >= 24
    table = value(reader->ReadTable());
#else
    check(reader->ReadTable(&table));
#endif
    check(table->ValidateFull());
    // run_uuid may be a plain string or a dictionary of strings with any index width (e.g. from pandas).
    const std::shared_ptr<arrow::DataType> run_uuid = table->schema()->num_fields() == 6 ? table->schema()->field(2)->type() : nullptr;
    const bool string_uuids = run_uuid
                              && (run_uuid->id() == arrow::Type::STRING
                                  || (run_uuid->id() == arrow::Type::DICTIONARY
                                      && static_cast<const arrow::DictionaryType&>(*run_uuid).value_type()->id() == arrow::Type::STRING));
    if (! string_uuids || ! table->schema()->Equals(*linkSchema(run_uuid), false)) invalid("Unexpected identification link table schema");
    // Batches keep the columns aligned without concatenating dictionary chunks.
    arrow::TableBatchReader batches(*table);
    std::shared_ptr<arrow::RecordBatch> batch;
    for (check(batches.ReadNext(&batch)); batch; check(batches.ReadNext(&batch)))
    {
      const auto& feature_ids = static_cast<const arrow::UInt64Array&>(*batch->column(0));
      const auto& links = static_cast<const arrow::StringArray&>(*batch->column(1));
      const RunUuids uuids(batch->column(2));
      const auto& record_ids = static_cast<const arrow::UInt64Array&>(*batch->column(3));
      const auto& encodings = static_cast<const arrow::UInt8Array&>(*batch->column(4));
      const auto& representations = static_cast<const arrow::StringArray&>(*batch->column(5));
      for (int64_t row = 0; row < batch->num_rows(); ++row)
      {
        if (feature_ids.IsNull(row) || links.IsNull(row)) invalid("Identification link without feature or kind");
        auto& entry = by_feature[feature_ids.Value(row)];
        const auto link = links.GetString(row);
        if (link == PRIMARY)
        {
          if (encodings.IsNull(row) || representations.IsNull(row) || encodings.Value(row) > static_cast<uint8_t>(ID::Encoding::DATABASE_ID))
            invalid("Invalid primary molecule link");
          if (entry.primary) invalid("Feature has more than one primary molecule");
          entry.primary = ID::MoleculeIdentity {static_cast<ID::Encoding>(encodings.Value(row)), representations.GetString(row)};
          continue;
        }
        if (uuids.isNull(row) || record_ids.IsNull(row)) invalid("Identification link without run UUID or record ID");
        if (link == QUERY) entry.queries.insert({uuids.at(row), ID::QueryId {record_ids.Value(row)}});
        else if (link == MATCH)
          entry.matches.insert({uuids.at(row), ID::MatchId {record_ids.Value(row)}});
        else
          invalid("Unknown identification link kind: " + link);
      }
    }
  }
  std::map<UInt64, BaseFeature*> by_id;
  std::set<UInt64> duplicates;
  for (auto* feature : features)
    if (! by_id.emplace(feature->getUniqueId(), feature).second) duplicates.insert(feature->getUniqueId());
  for (const auto& [unique_id, entry] : by_feature)
  {
    if (duplicates.contains(unique_id)) invalid("Identification links refer to a duplicated feature unique ID " + std::to_string(unique_id));
    if (! by_id.contains(unique_id)) invalid("Identification links refer to an unknown feature unique ID " + std::to_string(unique_id));
  }
  // Validate against the loaded data with temporary features before touching the destination.
  std::vector<BaseFeature> staged;
  staged.reserve(by_feature.size());
  for (const auto& [unique_id, entry] : by_feature)
  {
    BaseFeature feature;
    feature.setUniqueId(unique_id);
    feature.getIDQueries() = entry.queries;
    feature.getIDMatches() = entry.matches;
    staged.push_back(std::move(feature));
  }
  std::vector<const BaseFeature*> pointers;
  for (const auto& feature : staged)
    pointers.push_back(&feature);
  validateLinks(loaded, pointers);
  for (const auto& [unique_id, entry] : by_feature)
  {
    auto* feature = by_id.at(unique_id);
    if (entry.primary) feature->setPrimaryID(*entry.primary);
    feature->getIDQueries() = entry.queries;
    feature->getIDMatches() = entry.matches;
  }
  data.swap(loaded);
}

StagedBundle::StagedBundle(const std::string& target, const std::string& main_file): target_(target), main_file_(main_file)
{
  if (! target_.empty() && ! target_.has_filename()) target_ = target_.parent_path(); // "out.featureparquet/"
  if (target_.empty()) invalid("Output destination is empty");
  std::error_code ec;
  if (fs::exists(target_, ec) && ! replaceable_())
    fileError(target_, "Output exists and is neither empty nor a bundle containing " + main_file_);
  fs::path parent = target_.parent_path();
  if (parent.empty()) parent = ".";
  if (! fs::is_directory(parent, ec)) fileError(parent, "Output parent directory does not exist");
  for (Size attempt = 0; attempt < 100; ++attempt)
  {
    path_ = sibling_(".tmp-");
    if (fs::create_directory(path_, ec)) return;
    if (ec) fileError(path_, ec.message());
  }
  fileError(target_, "Cannot create a staging directory");
}

StagedBundle::~StagedBundle()
{
  if (! published_)
  {
    std::error_code ec;
    fs::remove_all(path_, ec);
  }
}

void StagedBundle::publish()
{
  std::error_code ec;
  const bool exists = fs::exists(target_, ec);
  fs::path previous;
  if (exists)
  {
    if (! replaceable_()) fileError(target_, "Output appeared while writing and is not a bundle containing " + main_file_);
    // Move the previous output aside, publish, then delete it. If publication fails, restore it.
    previous = sibling_(".old-");
    ec = renameWithRetry(target_, previous);
    if (ec) fileError(target_, "Cannot replace the existing output: " + ec.message());
  }
  ec = renameWithRetry(path_, target_);
  if (ec)
  {
    std::error_code restore;
    if (! previous.empty()) fs::rename(previous, target_, restore);
    fileError(target_, ec.message());
  }
  published_ = true;
  if (! previous.empty()) fs::remove_all(previous, ec); // best effort; the new bundle is already published
}

fs::path StagedBundle::sibling_(const std::string& infix) const
{
  fs::path parent = target_.parent_path();
  if (parent.empty()) parent = ".";
  return parent / (target_.filename().string() + infix + std::to_string(UniqueIdGenerator::getUniqueId()));
}

bool StagedBundle::replaceable_() const
{
  std::error_code ec;
  if (! fs::is_directory(target_, ec)) return false;
  if (fs::is_empty(target_, ec) && ! ec) return true;
  return fs::is_regular_file(target_ / main_file_, ec);
}
} // namespace OpenMS::Internal::MapIdentificationParquet
