// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
#include "IdentificationDataMatchTable.h"

#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/IdentificationDataArrow.h>
#include <arrow/record_batch.h>
#include <arrow/table.h>
#include <algorithm>
#include <cmath>
#include <set>

namespace OpenMS
{
namespace
{
  namespace IO = Internal::IdentificationDataIO;
  using ID = IdentificationData;
  using Arrow = IdentificationDataArrow;

  [[noreturn]] void invalidPatch(const std::string& message)
  { throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Identification data patch: " + message); }

  /// Builds one record batch of a schema; the writer interface of IO::writeMatch().
  class BatchWriter
  {
  public:
    explicit BatchWriter(std::shared_ptr<arrow::Schema> schema): schema_(std::move(schema))
    { builder_ = IO::value(arrow::RecordBatchBuilder::Make(schema_, arrow::default_memory_pool())); }
    arrow::ArrayBuilder& column(Size index)
    { return *builder_->GetField(static_cast<int>(index)); }
    Size columns() const
    { return static_cast<Size>(schema_->num_fields()); }
    void finishRow(Size)
    {}
    std::shared_ptr<arrow::RecordBatch> finish()
    { return IO::value(builder_->Flush()); }

  private:
    std::shared_ptr<arrow::Schema> schema_;
    std::unique_ptr<arrow::RecordBatchBuilder> builder_;
  };

  bool isCatalog(const ID::Run& run)
  {
    const auto& settings = run.getSettings();
    return settings.metaValueExists("identification:catalog") && settings.getMetaValue("identification:catalog").toString() == "true";
  }
  /// Runs of the score schema: those with score definitions or matches, except sequence catalogs.
  bool scored(const ID::Run& run)
  { return ! isCatalog(run) && (! run.getScoreDefinitions().empty() || run.getNumberOfMatches() > 0); }

  /// A run UUID of a key column: a string (large, view) or a dictionary of strings.
  std::string keyText(const arrow::Array& array, int64_t row)
  {
    if (array.IsNull(row)) invalidPatch("run_uuid must not be null");
    switch (array.type_id())
    {
      case arrow::Type::STRING:
        return static_cast<const arrow::StringArray&>(array).GetString(row);
      case arrow::Type::LARGE_STRING:
        return static_cast<const arrow::LargeStringArray&>(array).GetString(row);
      case arrow::Type::STRING_VIEW: // e.g. from polars
        return std::string(static_cast<const arrow::StringViewArray&>(array).GetView(row));
      case arrow::Type::DICTIONARY:
      {
        const auto& dictionary = static_cast<const arrow::DictionaryArray&>(array);
        return keyText(*dictionary.dictionary(), dictionary.GetValueIndex(row));
      }
      default:
        invalidPatch("run_uuid must be a string column");
    }
  }
  /// A non-negative integer of a key or value column.
  UInt64 unsignedValue(const arrow::Array& array, int64_t row, const std::string& column)
  {
    if (array.IsNull(row)) invalidPatch(column + " must not be null");
    Int64 value = 0;
    switch (array.type_id())
    {
      case arrow::Type::UINT64:
        return static_cast<const arrow::UInt64Array&>(array).Value(row);
      case arrow::Type::UINT32:
        return static_cast<const arrow::UInt32Array&>(array).Value(row);
      case arrow::Type::UINT16:
        return static_cast<const arrow::UInt16Array&>(array).Value(row);
      case arrow::Type::UINT8:
        return static_cast<const arrow::UInt8Array&>(array).Value(row);
      case arrow::Type::INT64:
        value = static_cast<const arrow::Int64Array&>(array).Value(row);
        break;
      case arrow::Type::INT32:
        value = static_cast<const arrow::Int32Array&>(array).Value(row);
        break;
      case arrow::Type::INT16:
        value = static_cast<const arrow::Int16Array&>(array).Value(row);
        break;
      case arrow::Type::INT8:
        value = static_cast<const arrow::Int8Array&>(array).Value(row);
        break;
      default:
        invalidPatch(column + " must be an integer column");
    }
    if (value < 0) invalidPatch(column + " must not be negative");
    return static_cast<UInt64>(value);
  }
  /// A score value: null or NaN is a missing value.
  std::optional<double> scoreValue(const arrow::Array& array, int64_t row, const std::string& column)
  {
    if (array.IsNull(row)) return std::nullopt;
    double value = 0;
    switch (array.type_id())
    {
      case arrow::Type::DOUBLE:
        value = static_cast<const arrow::DoubleArray&>(array).Value(row);
        break;
      case arrow::Type::FLOAT:
        value = static_cast<const arrow::FloatArray&>(array).Value(row);
        break;
      default:
        invalidPatch("score column " + column + " must be a floating-point column");
    }
    if (std::isnan(value)) return std::nullopt;
    if (! std::isfinite(value)) invalidPatch("score column " + column + " has a value that is not finite");
    return value;
  }

  /// The edits of one patch row.
  struct Edit
  {
    Size run;
    ID::MatchId match;
    std::vector<std::pair<UInt32, std::optional<double>>> scores; ///< Score index in the schema, value
    std::optional<ID::TargetDecoy> target_decoy;
  };
} // namespace

std::shared_ptr<arrow::Table> IdentificationDataArrow::matchTable(const ID& data)
{ return matchTable(data, IdentificationDataFile::Projection {}); }

std::shared_ptr<arrow::Table> IdentificationDataArrow::matchTable(const ID& data, const IdentificationDataFile::Projection& projection)
{
  data.validate();
  const auto& definitions = data.getScoreDefinitions();
  const auto score_columns = IO::scoreColumns(definitions);
  // One dictionary of metadata descriptors for the whole table, which the metadata column refers to.
  IO::Dictionary dictionary;
  for (const auto& run : data.getRuns())
    for (const auto& source : run.getSources())
      for (const auto& query : source.identifications)
        for (const auto& match : query.getMatches())
          dictionary.collect(match);
  const auto schema = IO::matchSchema(score_columns);
  BatchWriter writer(schema);
  arrow::Int32Builder run_rows;
  arrow::StringBuilder run_uuids;
  IO::Json revisions = IO::Json::object();
  for (Size r = 0; r < data.getRuns().size(); ++r)
  {
    const auto& run = data.getRuns()[r];
    IO::check(run_uuids.Append(run.getUuid()));
    revisions[run.getUuid()] = run.getRevision();
    std::vector<ID::ScoreView> views;
    for (UInt32 i = 0; i < run.getScoreDefinitions().size(); ++i)
      views.push_back(run.bindScore(run.getScoreId(i)));
    std::vector<std::optional<double>> scores(views.size());
    for (const auto& source : run.getSources())
      for (const auto& query : source.identifications)
        for (const auto& match : query.getMatches())
        {
          for (Size i = 0; i < views.size(); ++i)
            scores[i] = views[i](match);
          IO::writeMatch(writer, IO::MatchView {match.getId().value, query.getId().value, match.getData(), scores}, dictionary);
          IO::check(run_rows.Append(static_cast<int32_t>(r)));
        }
  }
  const auto batch = writer.finish();
  const auto run_uuid = IO::value(arrow::DictionaryArray::FromArrays(IO::keyType(), IO::value(run_rows.Finish()), IO::value(run_uuids.Finish())));
  IO::Json scores = IO::Json::array();
  for (const auto& definition : definitions)
    scores.push_back(IO::scoreJson(definition));
  IO::Json primary = nullptr;
  if (const auto definition = data.getPrimaryScoreDefinition())
    primary = static_cast<UInt64>(std::find(definitions.begin(), definitions.end(), *definition) - definitions.begin());
  auto metadata = std::make_shared<arrow::KeyValueMetadata>(
    std::vector<std::string> {"openms:score_definitions", "openms:metadata_descriptors", "openms:primary_score", "openms:revisions"},
    std::vector<std::string> {scores.dump(), dictionary.toJson().dump(), primary.dump(), revisions.dump()});
  auto fields = schema->fields();
  fields.push_back(arrow::field(IO::partitionColumn("matches"), IO::keyType(), false));
  auto columns = batch->columns();
  columns.push_back(run_uuid);
  auto table = arrow::Table::Make(arrow::schema(fields, metadata), columns, batch->num_rows());
  // The projection selects columns as IdentificationDataFile::scan() reads them; the keys always stay.
  std::vector<int> selected;
  for (const auto& name : IO::projectedMatchColumns(projection, score_columns))
    selected.push_back(table->schema()->GetFieldIndex(name));
  selected.push_back(table->num_columns() - 1);
  if (static_cast<int>(selected.size()) == table->num_columns()) return table;
  return IO::value(table->SelectColumns(selected));
}

IdentificationDataArrow::PatchResult IdentificationDataArrow::applyPatch(ID& data, const arrow::Table& patch)
{ return applyPatch(data, patch, PatchOptions {}); }

IdentificationDataArrow::PatchResult IdentificationDataArrow::applyPatch(ID& data, arrow::RecordBatchReader& patch, const PatchOptions& options)
{
  const auto table = IO::value(arrow::Table::FromRecordBatchReader(&patch));
  return applyPatch(data, *table, options);
}

IdentificationDataArrow::PatchResult IdentificationDataArrow::applyPatch(ID& data, const arrow::Table& patch, const PatchOptions& options)
{
  for (const auto& run : data.getRuns())
    run.checkMutation_();
  data.validate();
  PatchResult result;
  // The score schema after the additions, and the names of its columns.
  auto definitions = data.getScoreDefinitions();
  std::vector<ID::ScoreDefinition> added;
  for (const auto& definition : options.add_scores)
  {
    if (definition.name.empty()) invalidPatch("an added score definition needs a name");
    if (std::find(definitions.begin(), definitions.end(), definition) != definitions.end()) continue;
    definitions.push_back(definition);
    added.push_back(definition);
  }
  if (! added.empty() && std::none_of(data.getRuns().begin(), data.getRuns().end(), scored))
    invalidPatch("there is no run to add score definitions to");
  const auto names = IO::scoreColumns(definitions);
  for (Size i = definitions.size() - added.size(); i < definitions.size(); ++i)
    result.added.push_back(names[i]);
  std::optional<UInt32> primary;
  if (const auto definition = data.getPrimaryScoreDefinition())
    primary = static_cast<UInt32>(std::find(definitions.begin(), definitions.end(), *definition) - definitions.begin());

  // Columns: the keys, score columns by name, target_decoy.
  const auto& schema = *patch.schema();
  int uuid_column = -1, id_column = -1, target_decoy_column = -1;
  std::vector<std::pair<int, UInt32>> score_columns;
  std::set<std::string> seen_columns;
  for (int i = 0; i < schema.num_fields(); ++i)
  {
    const auto& name = schema.field(i)->name();
    if (! seen_columns.insert(name).second) invalidPatch("repeated column " + name);
    if (name == IO::partitionColumn("matches")) uuid_column = i;
    else if (name == "match_id")
      id_column = i;
    else if (name == "target_decoy")
      target_decoy_column = i;
    else
    {
      const auto found = std::find(names.begin(), names.end(), name);
      if (found == names.end())
        invalidPatch("unknown column " + name + "; a patch has the key columns run_uuid and match_id, score columns (" + ListUtils::concatenate(names, ", ")
                     + ") and target_decoy");
      score_columns.emplace_back(i, static_cast<UInt32>(found - names.begin()));
    }
    if (name != IO::partitionColumn("matches") && name != "match_id") result.columns.push_back(name);
  }
  if (uuid_column < 0 || id_column < 0) invalidPatch("a patch needs the key columns run_uuid and match_id");

  for (const auto& [uuid, revision] : options.expected_revisions)
  {
    const auto* run = data.findRunByUuid(uuid);
    if (! run) invalidPatch("expected revision of an unknown run " + uuid);
    if (run->getRevision() != revision)
      invalidPatch("run " + run->getIdentifier() + " changed since the patch was made (revision " + std::to_string(run->getRevision()) + ", expected "
                   + std::to_string(revision) + ")");
  }

  // Check every row before anything changes.
  std::map<std::string, Size> runs;
  for (Size r = 0; r < data.getRuns().size(); ++r)
    runs.emplace(data.getRuns()[r].getUuid(), r);
  std::vector<std::set<UInt64>> keys(data.getRuns().size());
  std::vector<Edit> edits;
  edits.reserve(static_cast<Size>(patch.num_rows()));
  const auto batch = IO::value(patch.CombineChunksToBatch());
  for (int64_t row = 0; row < batch->num_rows(); ++row)
  {
    const auto uuid = keyText(*batch->column(uuid_column), row);
    const auto run = runs.find(uuid);
    if (run == runs.end()) invalidPatch("unknown run " + uuid);
    const auto& target = data.getRuns()[run->second];
    const ID::MatchId match {unsignedValue(*batch->column(id_column), row, "match_id")};
    if (! target.findMatch(match)) invalidPatch("run " + target.getIdentifier() + " has no match " + std::to_string(match.value));
    if (! keys[run->second].insert(match.value).second)
      invalidPatch("match " + std::to_string(match.value) + " of run " + target.getIdentifier() + " appears more than once");
    Edit edit {run->second, match, {}, std::nullopt};
    for (const auto& [column, score] : score_columns)
    {
      if (! scored(target)) invalidPatch("run " + target.getIdentifier() + " is a sequence catalog without scores");
      const auto value = scoreValue(*batch->column(column), row, names[score]);
      if (! value && primary == score) invalidPatch("the primary score " + names[score] + " cannot be removed");
      edit.scores.emplace_back(score, value);
    }
    if (target_decoy_column >= 0)
    {
      const auto value = unsignedValue(*batch->column(target_decoy_column), row, "target_decoy");
      if (value > static_cast<UInt64>(ID::TargetDecoy::BOTH)) invalidPatch("target_decoy " + std::to_string(value) + " is out of range (0 unknown, 1 target, 2 decoy, 3 both)");
      edit.target_decoy = static_cast<ID::TargetDecoy>(value);
    }
    edits.push_back(std::move(edit));
  }

  // Edit copies of the runs that change: those with edits and, if scores are added, every run of the score schema.
  // A patch of keys only changes nothing.
  const bool values = ! score_columns.empty() || target_decoy_column >= 0;
  std::map<Size, ID::Run> staged;
  for (Size r = 0; r < data.getRuns().size(); ++r)
    if ((values && ! keys[r].empty()) || (! added.empty() && scored(data.getRuns()[r]))) staged.emplace(r, data.getRuns()[r]);
  for (auto& [index, run] : staged)
    if (scored(run))
      for (const auto& definition : added)
        run.addScore(definition);
  for (const auto& edit : edits)
  {
    if (edit.scores.empty() && ! edit.target_decoy) continue;
    auto& run = staged.at(edit.run);
    for (const auto& [score, value] : edit.scores)
      run.setScore(edit.match, run.getScoreId(score), value);
    if (edit.target_decoy && run.getMatch(edit.match).target_decoy != *edit.target_decoy)
    {
      auto payload = run.getMatch(edit.match).getData();
      payload.target_decoy = *edit.target_decoy;
      run.replaceMatch(edit.match, payload);
    }
  }
  for (const auto& [index, run] : staged)
    run.validate();
  // Commit: swapping run data cannot throw.
  for (auto& [index, run] : staged)
    data.findRunByUuid(run.getUuid())->swapData_(run);
  result.rows = edits.size();
  return result;
}
} // namespace OpenMS
