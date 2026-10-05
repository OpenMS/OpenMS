// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $
#include "IdentificationDataFileSupport.h"

#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/FORMAT/IdentificationDataFile.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>

namespace OpenMS
{
namespace
{
  namespace IO = Internal::IdentificationDataIO;
  namespace fs = std::filesystem;
  using ID = IdentificationData;
  using File = IdentificationDataFile;
  using IO::append;
  using IO::appendText;
  using IO::check;
  using IO::invalid;
  using IO::Json;
  using IO::number;
  using IO::optionalNumber;
  using IO::text;
  const std::string FORMAT = "OpenMS.IdentificationData";

  // All outputs remain private until every table is closed and the manifest is written.
  class StagedDirectory
  {
  public:
    explicit StagedDirectory(const std::string& target): target_(target)
    {
      if (target_.empty() || fs::exists(target_)) invalid("Output destination already exists or is empty");
      fs::path parent = target_.parent_path();
      if (parent.empty()) parent = ".";
      if (! fs::is_directory(parent)) invalid("Output parent directory does not exist");
      path = parent / (target_.filename().string() + ".tmp-" + std::to_string(UniqueIdGenerator::getUniqueId()));
      if (! fs::create_directory(path)) invalid("Cannot create staging directory");
    }
    ~StagedDirectory()
    {
      if (! published_)
      {
        std::error_code ec;
        fs::remove_all(path, ec);
      }
    }
    void publish()
    {
      if (fs::exists(target_)) invalid("Output destination appeared while writing");
      fs::rename(path, target_);
      published_ = true;
    }
    fs::path path;

  private:
    fs::path target_;
    bool published_ = false;
  };
  std::string numberDirectory(Size index)
  {
    std::ostringstream stream;
    stream << std::setw(3) << std::setfill('0') << index;
    return stream.str();
  }
  void validateJsonText(const Json& json)
  {
    if (json.is_string()) IO::validateText(json.get<std::string>());
    else if (json.is_object())
      for (auto i = json.begin(); i != json.end(); ++i)
      {
        IO::validateText(i.key());
        validateJsonText(i.value());
      }
    else if (json.is_array())
      for (const auto& v : json)
        validateJsonText(v);
  }
  Json readManifest(const std::string& path)
  try
  {
    std::ifstream in(fs::path(path) / "manifest.json", std::ios::binary);
    if (! in) invalid("Cannot open manifest.json");
    Json manifest;
    try
    {
      in >> manifest;
    }
    catch (const Json::exception& e)
    {
      invalid(e.what());
    }
    validateJsonText(manifest);
    if (manifest.at("format") != FORMAT || IO::integer<unsigned>(manifest.at("schema_version")) != 1)
      invalid("Unsupported native identification format");
    if (! manifest.at("runs").is_array() || ! manifest.at("inference").is_array()) invalid("Malformed run/result descriptors");
    std::set<std::string> paths;
    const auto claim = [&](const fs::path& base, const Json& entry) {
      const auto relative = entry.get<std::string>();
      const auto name = IO::tablePath(base, relative).lexically_normal().generic_string();
      if (! paths.insert(name).second) invalid("Two tables declare the same file path");
    };
    for (const auto& run : manifest.at("runs"))
    {
      const auto& tables = run.at("tables");
      if (! tables.is_object() || ! tables.contains("queries") || ! tables.contains("matches")
          || tables.size() != (tables.contains("parents") ? 3u : 2u))
        invalid("Invalid run table declarations");
      for (const auto& table : tables)
        claim({}, table);
    }
    for (const auto& result : manifest.at("inference"))
    {
      const auto directory = IO::tablePath({}, result.at("directory"));
      const auto& tables = result.at("tables");
      if (! tables.is_object() || tables.size() != 6) invalid("Inference requires six typed tables");
      for (const auto* name : {"inputs", "input_members", "proteins", "groups", "group_members", "assignments"})
      {
        claim(directory, tables.at(name));
      }
    }
    std::set<std::string> uuids, identifiers;
    for (const auto& run : manifest.at("runs"))
    {
      if (! uuids.insert(run.at("uuid").get<std::string>()).second || ! identifiers.insert(run.at("identifier").get<std::string>()).second)
        invalid("Duplicate run UUID or identifier");
    }
    return manifest;
  }
  catch (const Json::exception& error)
  {
    invalid(error.what());
  }
  void writeManifest(const fs::path& directory, const Json& manifest)
  {
    validateJsonText(manifest);
    std::ofstream out(directory / "manifest.json", std::ios::binary);
    if (! out) invalid("Cannot write manifest");
    out << manifest.dump(2) << '\n';
    out.close();
    if (! out) invalid("Failed to close manifest");
  }
  template<class Builder, class T>
  void appendOptional(arrow::ArrayBuilder& builder, const std::optional<T>& v)
  {
    if (v) append<Builder>(builder, *v);
    else
      check(builder.AppendNull());
  }
  std::shared_ptr<arrow::DataType> identityType()
  { return arrow::struct_({arrow::field("database", arrow::utf8(), false), arrow::field("accession", arrow::utf8(), false)}); }
  std::shared_ptr<arrow::Schema> querySchema()
  {
    return arrow::schema({arrow::field("query_id", arrow::uint64(), false), arrow::field("source_id", arrow::uint32(), false),
                          arrow::field("data_id", arrow::utf8(), false), arrow::field("rt", arrow::float64()), arrow::field("mz", arrow::float64()),
                          arrow::field("selected_match_id", arrow::uint64()), arrow::field("metadata", IO::metadataType(), false)});
  }
  std::shared_ptr<arrow::Schema> matchSchema(Size score_count)
  {
    std::vector<std::shared_ptr<arrow::Field>> fields {
      arrow::field("match_id", arrow::uint64(), false),
      arrow::field("query_id", arrow::uint64(), false),
      arrow::field("representation", arrow::utf8(), false),
      arrow::field("encoding", arrow::uint8(), false),
      arrow::field("charge", arrow::int32(), false),
      arrow::field("calculated_mz", arrow::float64()),
      arrow::field("target_decoy", arrow::uint8(), false),
      arrow::field("name", arrow::utf8(), false),
      arrow::field("formula", arrow::utf8()),
      arrow::field("identifiers", arrow::list(arrow::field("item", identityType(), false)), false),
      arrow::field("adduct", arrow::struct_({arrow::field("name", arrow::utf8(), false), arrow::field("formula", arrow::utf8(), false),
                                             arrow::field("charge", arrow::int32(), false), arrow::field("multiplier", arrow::uint32(), false)})),
      arrow::field(
        "parent_evidence",
        arrow::list(arrow::field("item",
                                 arrow::struct_({arrow::field("database", arrow::utf8(), false), arrow::field("accession", arrow::utf8(), false),
                                                 arrow::field("start", arrow::uint64()), arrow::field("end", arrow::uint64()),
                                                 arrow::field("before", arrow::utf8(), false), arrow::field("after", arrow::utf8(), false)}),
                                 false)),
        false),
      arrow::field(
        "peak_annotations",
        arrow::list(arrow::field("item",
                                 arrow::struct_({arrow::field("annotation", arrow::utf8(), false), arrow::field("charge", arrow::int32(), false),
                                                 arrow::field("mz", arrow::float64(), false), arrow::field("intensity", arrow::float64(), false)}),
                                 false)),
        false),
      arrow::field("metadata", IO::metadataType(), false)};
    for (Size i = 0; i < score_count; ++i)
      fields.push_back(arrow::field("score_" + std::to_string(i), arrow::float64()));
    return arrow::schema(fields);
  }
  Size queryBytes(const File::QueryRecord& q)
  { return 64 + q.data.data_id.size() + IO::metadataBytes(q.data); }
  Size matchBytes(const File::MatchRecord& m)
  {
    const auto& d = m.data;
    Size bytes = 128 + d.representation.size() + d.name.size() + (d.formula ? d.formula->size() : 0) + m.scores.size() * 9 + IO::metadataBytes(d);
    for (const auto& i : d.identifiers)
      bytes += 16 + i.database.size() + i.accession.size();
    if (d.adduct) bytes += 32 + d.adduct->getName().size() + d.adduct->getEmpiricalFormula().toString().size();
    for (const auto& e : d.parent_evidence)
      bytes += 40 + e.parent.database.size() + e.parent.accession.size() + e.before.size() + e.after.size();
    for (const auto& a : d.peak_annotations)
      bytes += 32 + a.annotation.size();
    return bytes;
  }
  void writeQuery(IO::TableWriter& writer, const File::QueryRecord& q, const IO::Dictionary& dictionary)
  {
    append<arrow::UInt64Builder>(writer.column(0), q.query_id);
    append<arrow::UInt32Builder>(writer.column(1), q.source_id);
    appendText(writer.column(2), q.data.data_id);
    appendOptional<arrow::DoubleBuilder>(writer.column(3), q.data.rt);
    appendOptional<arrow::DoubleBuilder>(writer.column(4), q.data.mz);
    appendOptional<arrow::UInt64Builder>(writer.column(5), q.selected_match_id);
    IO::appendMetadata(writer.column(6), q.data, dictionary);
    writer.finishRow(queryBytes(q));
  }
  void writeMatch(IO::TableWriter& writer, const File::MatchRecord& m, const IO::Dictionary& dictionary)
  {
    const auto& d = m.data;
    append<arrow::UInt64Builder>(writer.column(0), m.match_id);
    append<arrow::UInt64Builder>(writer.column(1), m.query_id);
    appendText(writer.column(2), d.representation);
    append<arrow::UInt8Builder>(writer.column(3), static_cast<unsigned>(d.encoding));
    append<arrow::Int32Builder>(writer.column(4), d.charge);
    appendOptional<arrow::DoubleBuilder>(writer.column(5), d.calculated_mz);
    append<arrow::UInt8Builder>(writer.column(6), static_cast<unsigned>(d.target_decoy));
    appendText(writer.column(7), d.name);
    if (d.formula) appendText(writer.column(8), *d.formula);
    else
      check(writer.column(8).AppendNull());
    auto& identifiers = static_cast<arrow::ListBuilder&>(writer.column(9));
    check(identifiers.Append());
    auto& identity = *static_cast<arrow::StructBuilder*>(identifiers.value_builder());
    for (const auto& i : d.identifiers)
    {
      check(identity.Append());
      appendText(*identity.field_builder(0), i.database);
      appendText(*identity.field_builder(1), i.accession);
    }
    auto& adduct = static_cast<arrow::StructBuilder&>(writer.column(10));
    if (! d.adduct) check(adduct.AppendNull());
    else
    {
      check(adduct.Append());
      appendText(*adduct.field_builder(0), d.adduct->getName());
      appendText(*adduct.field_builder(1), d.adduct->getEmpiricalFormula().toString());
      append<arrow::Int32Builder>(*adduct.field_builder(2), d.adduct->getCharge());
      append<arrow::UInt32Builder>(*adduct.field_builder(3), d.adduct->getMolMultiplier());
    }
    auto& evidence = static_cast<arrow::ListBuilder&>(writer.column(11));
    check(evidence.Append());
    auto& entry = *static_cast<arrow::StructBuilder*>(evidence.value_builder());
    for (const auto& e : d.parent_evidence)
    {
      check(entry.Append());
      appendText(*entry.field_builder(0), e.parent.database);
      appendText(*entry.field_builder(1), e.parent.accession);
      appendOptional<arrow::UInt64Builder>(*entry.field_builder(2), e.start);
      appendOptional<arrow::UInt64Builder>(*entry.field_builder(3), e.end);
      appendText(*entry.field_builder(4), e.before);
      appendText(*entry.field_builder(5), e.after);
    }
    auto& annotations = static_cast<arrow::ListBuilder&>(writer.column(12));
    check(annotations.Append());
    auto& annotation = *static_cast<arrow::StructBuilder*>(annotations.value_builder());
    for (const auto& a : d.peak_annotations)
    {
      check(annotation.Append());
      appendText(*annotation.field_builder(0), a.annotation);
      append<arrow::Int32Builder>(*annotation.field_builder(1), a.charge);
      append<arrow::DoubleBuilder>(*annotation.field_builder(2), a.mz);
      append<arrow::DoubleBuilder>(*annotation.field_builder(3), a.intensity);
    }
    IO::appendMetadata(writer.column(13), d, dictionary);
    for (Size i = 0; i < m.scores.size(); ++i)
      appendOptional<arrow::DoubleBuilder>(writer.column(14 + i), m.scores[i]);
    writer.finishRow(matchBytes(m));
  }
  template<class Callback>
  void readList(const arrow::Array& array, int64_t row, const Callback& callback)
  {
    if (array.IsNull(row)) invalid("Null required list");
    const auto& list = static_cast<const arrow::ListArray&>(array);
    const auto& entries = static_cast<const arrow::StructArray&>(*list.values());
    for (int64_t i = list.value_offset(row), end = i + list.value_length(row); i < end; ++i)
    {
      if (entries.IsNull(i)) invalid("Null list item");
      callback(entries, i);
    }
  }
  File::QueryRecord readQuery(const IO::TableReader& reader, const IO::Dictionary& dictionary)
  {
    const auto row = reader.row();
    File::QueryRecord q;
    q.query_id = number<arrow::UInt64Array>(reader.column("query_id"), row);
    q.source_id = number<arrow::UInt32Array>(reader.column("source_id"), row);
    q.data.data_id = text(reader.column("data_id"), row);
    q.data.rt = optionalNumber<arrow::DoubleArray>(reader.column("rt"), row);
    q.data.mz = optionalNumber<arrow::DoubleArray>(reader.column("mz"), row);
    q.selected_match_id = optionalNumber<arrow::UInt64Array>(reader.column("selected_match_id"), row);
    if (reader.hasColumn("metadata")) IO::readMetadata(reader.column("metadata"), row, q.data, dictionary);
    if ((q.data.rt && ! std::isfinite(*q.data.rt)) || (q.data.mz && ! std::isfinite(*q.data.mz))) invalid("Nonfinite observation coordinate");
    return q;
  }
  File::MatchRecord readMatch(const IO::TableReader& reader, const IO::Dictionary& dictionary, Size score_count)
  {
    const auto row = reader.row();
    File::MatchRecord m;
    m.match_id = number<arrow::UInt64Array>(reader.column("match_id"), row);
    m.query_id = number<arrow::UInt64Array>(reader.column("query_id"), row);
    auto& d = m.data;
    if (reader.hasColumn("representation"))
    {
      d.representation = text(reader.column("representation"), row);
      auto encoding = number<arrow::UInt8Array>(reader.column("encoding"), row);
      if (encoding > static_cast<unsigned>(ID::Encoding::DATABASE_ID)) invalid("Unknown molecular encoding");
      d.encoding = static_cast<ID::Encoding>(encoding);
      d.charge = number<arrow::Int32Array>(reader.column("charge"), row);
      d.calculated_mz = optionalNumber<arrow::DoubleArray>(reader.column("calculated_mz"), row);
      if (d.calculated_mz && ! std::isfinite(*d.calculated_mz)) invalid("Nonfinite calculated m/z");
      auto target_decoy = number<arrow::UInt8Array>(reader.column("target_decoy"), row);
      if (target_decoy > static_cast<unsigned>(ID::TargetDecoy::BOTH)) invalid("Unknown target/decoy state");
      d.target_decoy = static_cast<ID::TargetDecoy>(target_decoy);
      d.name = text(reader.column("name"), row);
      if (! reader.column("formula").IsNull(row)) d.formula = text(reader.column("formula"), row);
      readList(reader.column("identifiers"), row,
               [&](const arrow::StructArray& a, int64_t i) { d.identifiers.push_back({text(*a.field(0), i), text(*a.field(1), i)}); });
      const auto& adduct = static_cast<const arrow::StructArray&>(reader.column("adduct"));
      if (! adduct.IsNull(row))
      {
        d.adduct.emplace(text(*adduct.field(0), row), EmpiricalFormula(text(*adduct.field(1), row)), number<arrow::Int32Array>(*adduct.field(2), row),
                         number<arrow::UInt32Array>(*adduct.field(3), row));
        if (d.adduct->getCharge() != d.charge) invalid("Adduct and molecular charge disagree");
      }
    }
    if (reader.hasColumn("parent_evidence"))
      readList(reader.column("parent_evidence"), row, [&](const arrow::StructArray& a, int64_t i) {
        ID::ParentEvidence e;
        e.parent = {text(*a.field(0), i), text(*a.field(1), i)};
        e.start = optionalNumber<arrow::UInt64Array>(*a.field(2), i);
        e.end = optionalNumber<arrow::UInt64Array>(*a.field(3), i);
        if (e.start && e.end && *e.start > *e.end) invalid("Reversed parent evidence positions");
        e.before = text(*a.field(4), i);
        e.after = text(*a.field(5), i);
        d.parent_evidence.push_back(std::move(e));
      });
    if (reader.hasColumn("peak_annotations"))
      readList(reader.column("peak_annotations"), row, [&](const arrow::StructArray& a, int64_t i) {
        PeptideHit::PeakAnnotation item;
        item.annotation = text(*a.field(0), i);
        item.charge = number<arrow::Int32Array>(*a.field(1), i);
        item.mz = number<arrow::DoubleArray>(*a.field(2), i);
        item.intensity = number<arrow::DoubleArray>(*a.field(3), i);
        d.peak_annotations.push_back(std::move(item));
      });
    if (reader.hasColumn("metadata")) IO::readMetadata(reader.column("metadata"), row, d, dictionary);
    m.scores.resize(score_count);
    for (Size i = 0; i < score_count; ++i)
    {
      const std::string name = "score_" + std::to_string(i);
      if (reader.hasColumn(name))
      {
        m.scores[i] = optionalNumber<arrow::DoubleArray>(reader.column(name), row);
        if (m.scores[i] && ! std::isfinite(*m.scores[i])) invalid("Nonfinite score");
      }
    }
    return m;
  }
  Json sourceJson(const ID::SourceFile& source)
  {
    return {
      {"identifier", source.identifier}, {"path", source.path}, {"primary_files", source.primary_files}, {"metadata", IO::metadataJson(source)}};
  }
  ID::SourceFile readSourceJson(const Json& j)
  {
    ID::SourceFile source;
    source.identifier = j.at("identifier");
    source.path = j.at("path");
    source.primary_files = j.at("primary_files").get<std::vector<std::string>>();
    IO::readMetadataJson(j.at("metadata"), source);
    return source;
  }
  Json runJson(const ID::Run& run)
  {
    const auto& processing = run.getProcessingMetadata();
    if (! processing.getHits().empty() || ! processing.getProteinGroups().empty() || ! processing.getIndistinguishableProteins().empty())
      invalid("Run processing metadata must contain configuration only; use parent catalogue and inference results for protein values");
    Json scores = Json::array(), sources = Json::array();
    for (const auto& s : run.getScoreDefinitions())
      scores.push_back(IO::scoreJson(s));
    for (const auto& s : run.getSourceBlocks())
      sources.push_back(sourceJson(s.source));
    return {{"identifier", run.getIdentifier()},
            {"uuid", run.getUuid()},
            {"molecule_kind", run.getMoleculeKind()},
            {"processing", IO::processingJson(run.getProcessingMetadata())},
            {"scores", std::move(scores)},
            {"sources", std::move(sources)},
            {"primary_score", run.getPrimaryScore() ? Json(run.getPrimaryScore()->value) : Json()},
            {"next_query_id", run.getNextQueryId()},
            {"next_match_id", run.getNextMatchId()},
            {"query_count", run.getNumberOfIdentifications()},
            {"match_count", run.getNumberOfMatches()}};
  }
  ID::Run runShell(const Json& j)
  {
    unsigned kind = IO::integer<unsigned>(j.at("molecule_kind"));
    if (kind > static_cast<unsigned>(ID::MoleculeKind::COMPOUND)) invalid("Unknown molecule kind");
    ID::Run run(j.at("identifier").get<std::string>(), static_cast<ID::MoleculeKind>(kind));
    run.setProcessingMetadata(IO::readProcessingJson(j.at("processing")));
    if (! j.at("scores").is_array() || ! j.at("sources").is_array()) invalid("Run descriptors must be arrays");
    for (const auto& s : j.at("scores"))
      run.addScore(IO::readScoreJson(s));
    for (const auto& s : j.at("sources"))
      run.addSource(readSourceJson(s));
    if (! j.at("primary_score").is_null()) run.setPrimaryScore(run.getScoreId(IO::integer<UInt32>(j.at("primary_score"))));
    return run;
  }
  File::RunDescriptor descriptor(const Json& j)
  {
    auto shell = runShell(j);
    shell.restoreIdentity(j.at("uuid"), IO::integer<UInt64>(j.at("next_query_id")), IO::integer<UInt64>(j.at("next_match_id")));
    File::RunDescriptor d;
    d.identifier = shell.getIdentifier();
    d.uuid = shell.getUuid();
    d.molecule_kind = shell.getMoleculeKind();
    d.scores = shell.getScoreDefinitions();
    for (const auto& source : shell.getSourceBlocks())
      d.sources.push_back(source.source);
    if (shell.getPrimaryScore()) d.primary_score = shell.getPrimaryScore()->value;
    d.next_query_id = shell.getNextQueryId();
    d.next_match_id = shell.getNextMatchId();
    d.query_count = IO::integer<UInt64>(j.at("query_count"));
    d.match_count = IO::integer<UInt64>(j.at("match_count"));
    return d;
  }
  std::vector<const Json*> selectRuns(const Json& manifest, const std::vector<std::string>& selection)
  {
    std::set<std::string> selected;
    for (const auto& key : selection)
    {
      Size found = 0;
      for (const auto& j : manifest.at("runs"))
        if (j.at("identifier") == key || j.at("uuid") == key)
        {
          selected.insert(j.at("uuid").get<std::string>());
          ++found;
        }
      if (found != 1) invalid("Unknown or ambiguous run selection: " + key);
    }
    std::vector<const Json*> runs;
    for (const auto& j : manifest.at("runs"))
      if (selection.empty() || selected.contains(j.at("uuid").get<std::string>())) runs.push_back(&j);
    return runs;
  }
  void collectDictionary(const ID::Run& run, IO::Dictionary& dictionary)
  {
    for (const auto& source : run.getSourceBlocks())
      for (const auto& query : source.identifications)
      {
        dictionary.collect(query);
        for (const auto& match : query.getMatches())
          dictionary.collect(match);
      }
    if (run.getParents())
      for (const auto& parent : *run.getParents())
        dictionary.collect(parent);
  }
  std::vector<std::string> projectedMatchColumns(const File::Projection& projection, Size score_count)
  {
    std::vector<std::string> fields {"match_id", "query_id"};
    if (projection.molecule)
      for (const auto* field : {"representation", "encoding", "charge", "calculated_mz", "target_decoy", "name", "formula", "identifiers", "adduct"})
        fields.emplace_back(field);
    if (projection.evidence) fields.emplace_back("parent_evidence");
    if (projection.annotations) fields.emplace_back("peak_annotations");
    if (projection.metadata) fields.emplace_back("metadata");
    std::set<UInt32> scores(projection.score_ids.begin(), projection.score_ids.end());
    for (UInt32 score : scores)
      if (score >= score_count) invalid("Projected score ID does not exist in the run");
    for (Size i = 0; i < score_count; ++i)
      if (projection.all_scores || scores.contains(static_cast<UInt32>(i))) fields.push_back("score_" + std::to_string(i));
    return fields;
  }
  // Both tables retain scientific ordering. A single match row is sufficient lookahead.
  template<class Begin, class Match, class End>
  void walkRun(const fs::path& root, const Json& j, const File::ScanOptions& options, const Begin& begin, const Match& consume, const End& end)
  {
    const auto d = descriptor(j);
    IO::Dictionary dictionary;
    dictionary.load(j.at("metadata_descriptors"));
    const auto& tables = j.at("tables");
    std::vector<std::string> query_columns {"query_id", "source_id", "data_id", "rt", "mz", "selected_match_id"};
    if (options.projection.metadata) query_columns.emplace_back("metadata");
    IO::TableReader queries(IO::tablePath(root, tables.at("queries")), querySchema(), options.buffering, query_columns);
    IO::TableReader matches(IO::tablePath(root, tables.at("matches")), matchSchema(d.scores.size()), options.buffering,
                            projectedMatchColumns(options.projection, d.scores.size()));
    if (queries.rows() != d.query_count || matches.rows() != d.match_count) invalid("Manifest and table row counts disagree");
    std::set<UInt64> query_ids, match_ids;
    bool has_match = matches.next();
    File::MatchRecord pending;
    if (has_match) pending = readMatch(matches, dictionary, d.scores.size());
    UInt64 seen_queries = 0, seen_matches = 0;
    UInt32 prior_source = 0;
    while (queries.next())
    {
      File::QueryRecord q = readQuery(queries, dictionary);
      ++seen_queries;
      if (! q.query_id || q.query_id >= d.next_query_id || q.source_id >= d.sources.size() || q.source_id < prior_source)
        invalid("Invalid query ID, source or source ordering");
      prior_source = q.source_id;
      if (q.selected_match_id && (! *q.selected_match_id || *q.selected_match_id >= d.next_match_id)) invalid("Invalid selected match ID");
      if (queryBytes(q) > options.buffering.max_record_bytes) invalid("Query exceeds max_record_bytes");
      if (options.validate_unique_ids && ! query_ids.insert(q.query_id).second) invalid("Duplicate query ID");
      begin(q);
      bool found_selected = ! q.selected_match_id;
      while (has_match && pending.query_id == q.query_id)
      {
        if (! pending.match_id || pending.match_id >= d.next_match_id) invalid("Invalid match ID or allocation counter");
        if (options.projection.molecule)
        {
          const auto encoding = pending.data.encoding;
          const bool compatible
            = encoding == ID::Encoding::DATABASE_ID || (d.molecule_kind == ID::MoleculeKind::PEPTIDE && encoding == ID::Encoding::AA_SEQUENCE)
              || (d.molecule_kind == ID::MoleculeKind::OLIGONUCLEOTIDE && encoding == ID::Encoding::NA_SEQUENCE)
              || (d.molecule_kind == ID::MoleculeKind::COMPOUND && (encoding == ID::Encoding::SMILES || encoding == ID::Encoding::INCHI));
          if (pending.data.representation.empty() || ! compatible) invalid("Empty or incompatible molecular representation");
        }
        if (d.molecule_kind == ID::MoleculeKind::COMPOUND && ! pending.data.parent_evidence.empty())
          invalid("Compound candidate has sequence evidence");
        for (const auto& evidence : pending.data.parent_evidence)
          if (evidence.parent.accession.empty()) invalid("Parent evidence needs an accession");
        if (matchBytes(pending) > options.buffering.max_record_bytes) invalid("Match exceeds max_record_bytes");
        if (options.validate_unique_ids && ! match_ids.insert(pending.match_id).second) invalid("Duplicate match ID");
        if (d.primary_score)
        {
          const auto name = "score_" + std::to_string(*d.primary_score);
          if (matches.hasColumn(name) && ! pending.scores[*d.primary_score]) invalid("Primary score is missing");
        }
        found_selected = found_selected || q.selected_match_id == pending.match_id;
        consume(pending);
        ++seen_matches;
        has_match = matches.next();
        if (has_match) pending = readMatch(matches, dictionary, d.scores.size());
      }
      if (! found_selected) invalid("Selected match does not belong to its query");
      end(q);
    }
    if (has_match || seen_queries != d.query_count || seen_matches != d.match_count) invalid("Unmatched or out-of-order match rows");
  }
  ID::Run readRun(const fs::path& root, const Json& j, const File::Options& options)
  {
    auto run = runShell(j);
    File::ScanOptions scan;
    scan.buffering = options;
    walkRun(
      root, j, scan, [&](const File::QueryRecord& q) { run.importIdentification(run.getSourceId(q.source_id), ID::QueryId {q.query_id}, q.data); },
      [&](const File::MatchRecord& m) { run.importMatch(ID::QueryId {m.query_id}, ID::MatchId {m.match_id}, m.data, m.scores); },
      [&](const File::QueryRecord& q) {
        if (q.selected_match_id) run.setSelectedMatch(ID::QueryId {q.query_id}, ID::MatchId {*q.selected_match_id});
      });
    const auto& tables = j.at("tables");
    if (tables.contains("parents"))
    {
      IO::Dictionary dictionary;
      dictionary.load(j.at("metadata_descriptors"));
      run.setParents(IO::readParents(IO::tablePath(root, tables.at("parents")), dictionary, options));
    }
    run.restoreIdentity(j.at("uuid"), IO::integer<UInt64>(j.at("next_query_id")), IO::integer<UInt64>(j.at("next_match_id")));
    return run;
  }
  void validateProvenanceCounters(const ID& data, const ID::InferenceResult& result)
  {
    const auto check_id = [&](const std::string& uuid, ID::MatchId match) {
      if (! match.value) invalid("Zero match ID in inference provenance");
      const auto* run = data.findRunByUuid(uuid);
      if (run && match.value >= run->getNextMatchId()) invalid("Run allocation counter does not reserve retained inference IDs");
    };
    for (const auto& input : result.inputs)
      for (auto match : input.matches)
        check_id(input.run_uuid, match);
    for (const auto& assignment : result.assignments)
      check_id(assignment.run_uuid, assignment.match);
  }
} // namespace

bool File::isNativeFile(const std::string& path)
{
  const auto manifest = fs::path(path) / "manifest.json";
  if (! fs::exists(manifest)) return false;
  std::ifstream in(manifest, std::ios::binary);
  if (! in) invalid("Cannot open manifest.json");
  Json j;
  try
  {
    in >> j;
  }
  catch (const Json::exception& e)
  {
    invalid(e.what());
  }
  return j.contains("format") && j.at("format") == FORMAT;
}
std::vector<File::RunDescriptor> File::inspect(const std::string& path)
try
{
  const Json manifest = readManifest(path);
  std::vector<RunDescriptor> result;
  for (const auto& j : manifest.at("runs"))
    result.push_back(descriptor(j));
  return result;
}
catch (const Json::exception& error)
{
  invalid(error.what());
}
void File::store(const std::string& path, const ID& data)
{ store(path, data, Options {}); }
void File::store(const std::string& path, const ID& data, const Options& options)
{
  IO::validateOptions(options);
  data.validate();
  StagedDirectory output(path);
  Json manifest {{"format", FORMAT}, {"schema_version", 1}, {"runs", Json::array()}, {"inference", Json::array()}};
  for (Size index = 0; index < data.getRuns().size(); ++index)
  {
    const auto& run = data.getRuns()[index];
    const std::string directory = "runs/" + numberDirectory(index);
    fs::create_directories(output.path / directory);
    IO::Dictionary dictionary;
    collectDictionary(run, dictionary);
    Json j = runJson(run);
    Json tables {{"queries", directory + "/queries.parquet"}, {"matches", directory + "/matches.parquet"}};
    IO::TableWriter queries(output.path / tables.at("queries").get<std::string>(), querySchema(), options);
    IO::TableWriter matches(output.path / tables.at("matches").get<std::string>(), matchSchema(run.getScoreDefinitions().size()), options);
    for (const auto& source : run.getSourceBlocks())
      for (const auto& query : source.identifications)
      {
        QueryRecord q {query.getId().value, source.id.value, query.getObservation(), std::nullopt};
        if (query.getSelectedMatch()) q.selected_match_id = query.getSelectedMatch()->value;
        writeQuery(queries, q, dictionary);
        for (const auto& match : query.getMatches())
          writeMatch(matches, MatchRecord {match.getId().value, query.getId().value, match.getData(), match.getScores()}, dictionary);
      }
    queries.close();
    matches.close();
    if (run.getParents())
    {
      tables["parents"] = directory + "/parents.parquet";
      IO::writeParents(output.path / tables.at("parents").get<std::string>(), *run.getParents(), dictionary, options);
    }
    j["metadata_descriptors"] = dictionary.toJson();
    j["tables"] = std::move(tables);
    manifest["runs"].push_back(std::move(j));
  }
  for (Size index = 0; index < data.getInferenceResults().size(); ++index)
  {
    const auto& result = data.getInferenceResults()[index];
    validateProvenanceCounters(data, result);
    const std::string directory = "inference/" + numberDirectory(index);
    fs::create_directories(output.path / directory);
    Json j = IO::writeInference(output.path / directory, result, options);
    j["directory"] = directory;
    manifest["inference"].push_back(std::move(j));
  }
  writeManifest(output.path, manifest);
  output.publish();
}
void File::load(const std::string& path, ID& data)
{ load(path, data, Options {}); }
void File::load(const std::string& path, ID& data, const Options& options)
try
{
  IO::validateOptions(options);
  const auto manifest = readManifest(path);
  ID temporary;
  for (const auto& j : manifest.at("runs"))
    temporary.addRun(readRun(path, j, options));
  for (const auto& j : manifest.at("inference"))
  {
    auto result = IO::readInference(IO::tablePath(path, j.at("directory")), j, options);
    validateProvenanceCounters(temporary, result);
    temporary.addInferenceResult(result);
  }
  temporary.validate();
  data.swap(temporary);
}
catch (const Json::exception& error)
{
  invalid(error.what());
}
ID::Run File::loadRun(const std::string& path, const std::string& run)
{ return loadRun(path, run, Options {}); }
ID::Run File::loadRun(const std::string& path, const std::string& run, const Options& options)
try
{
  IO::validateOptions(options);
  const auto manifest = readManifest(path);
  const auto selected = selectRuns(manifest, {run});
  return readRun(path, *selected.front(), options);
}
catch (const Json::exception& error)
{
  invalid(error.what());
}
File::ScanStatistics File::scan(const std::string& path, const ScanOptions& options, const QueryCallback& on_queries, const MatchCallback& on_matches)
try
{
  IO::validateOptions(options.buffering);
  const auto manifest = readManifest(path);
  ScanStatistics statistics;
  statistics.descriptor_bytes = fs::file_size(fs::path(path) / "manifest.json");
  for (const auto* j : selectRuns(manifest, options.runs))
  {
    const std::string uuid = j->at("uuid");
    std::vector<QueryRecord> queries;
    std::vector<MatchRecord> matches;
    Size query_bytes = 0, match_bytes = 0;
    auto flush_queries = [&]() {
      if (! queries.empty() && on_queries) on_queries(uuid, queries);
      queries.clear();
      query_bytes = 0;
    };
    auto flush_matches = [&]() {
      if (! matches.empty() && on_matches) on_matches(uuid, matches);
      matches.clear();
      match_bytes = 0;
    };
    walkRun(
      path, *j, options,
      [&](const QueryRecord& q) {
        ++statistics.queries;
        if (on_queries)
        {
          queries.push_back(q);
          query_bytes += queryBytes(q);
          if (queries.size() >= options.buffering.batch_rows || query_bytes >= options.buffering.batch_bytes) flush_queries();
        }
      },
      [&](const MatchRecord& m) {
        ++statistics.matches;
        if (on_matches)
        {
          matches.push_back(m);
          match_bytes += matchBytes(m);
          if (matches.size() >= options.buffering.batch_rows || match_bytes >= options.buffering.batch_bytes) flush_matches();
        }
      },
      [](const QueryRecord&) {});
    flush_queries();
    flush_matches();
  }
  return statistics;
}
catch (const Json::exception& error)
{
  invalid(error.what());
}
void File::filter(const std::string& input,
                  const std::string& output,
                  const MatchPredicate& keep,
                  ID::InferencePolicy policy,
                  bool keep_empty_queries)
{ filter(input, output, keep, policy, keep_empty_queries, Options {}); }
void File::filter(const std::string& input,
                  const std::string& output,
                  const MatchPredicate& keep,
                  ID::InferencePolicy policy,
                  bool keep_empty_queries,
                  const Options& options)
try
{
  IO::validateOptions(options);
  if (! keep) invalid("Filter requires a predicate");
  if (policy != ID::InferencePolicy::PRESERVE && policy != ID::InferencePolicy::DISCARD) invalid("Unknown inference policy");
  auto manifest = readManifest(input);
  std::map<std::string, UInt64> counters;
  for (const auto& run : manifest.at("runs"))
    counters.emplace(run.at("uuid").get<std::string>(), IO::integer<UInt64>(run.at("next_match_id")));
  StagedDirectory staged(output);
  for (auto& j : manifest["runs"])
  {
    const std::string uuid = j.at("uuid");
    const auto& tables = j.at("tables");
    IO::Dictionary dictionary;
    dictionary.load(j.at("metadata_descriptors"));
    auto query_path = IO::tablePath(staged.path, tables.at("queries"));
    fs::create_directories(query_path.parent_path());
    auto match_path = IO::tablePath(staged.path, tables.at("matches"));
    fs::create_directories(match_path.parent_path());
    IO::TableWriter queries(query_path, querySchema(), options);
    IO::TableWriter matches(match_path, matchSchema(j.at("scores").size()), options);
    File::ScanOptions scan;
    scan.buffering = options;
    bool kept_any = false, kept_selection = false;
    std::optional<UInt64> selected;
    walkRun(
      input, j, scan,
      [&](const QueryRecord& q) {
        kept_any = false;
        kept_selection = false;
        selected = q.selected_match_id;
      },
      [&](const MatchRecord& m) {
        if (keep(uuid, m))
        {
          writeMatch(matches, m, dictionary);
          kept_any = true;
          if (selected == m.match_id) kept_selection = true;
        }
      },
      [&](const QueryRecord& original) {
        // Keep the explicit selection only when its candidate survives.
        QueryRecord q = original;
        if (q.selected_match_id && ! kept_selection) q.selected_match_id.reset();
        if (kept_any || keep_empty_queries) writeQuery(queries, q, dictionary);
      });
    queries.close();
    matches.close();
    j["query_count"] = queries.rows();
    j["match_count"] = matches.rows();
    if (tables.contains("parents"))
    {
      auto dest = IO::tablePath(staged.path, tables.at("parents"));
      fs::create_directories(dest.parent_path());
      const auto source = IO::tablePath(input, tables.at("parents"));
      IO::validateParents(source, dictionary, options);
      fs::copy_file(source, dest);
    }
  }
  if (policy == ID::InferencePolicy::DISCARD) manifest["inference"] = Json::array();
  else
    for (const auto& result : manifest.at("inference"))
    {
      const auto source = IO::tablePath(input, result.at("directory"));
      const auto destination = IO::tablePath(staged.path, result.at("directory"));
      IO::validateInferenceTables(source, result, options, counters);
      fs::create_directories(destination);
      for (auto i = result.at("tables").begin(); i != result.at("tables").end(); ++i)
      {
        auto dest = IO::tablePath(destination, i.value());
        fs::create_directories(dest.parent_path());
        fs::copy_file(IO::tablePath(source, i.value()), dest);
      }
    }
  writeManifest(staged.path, manifest);
  staged.publish();
}
catch (const Json::exception& error)
{
  invalid(error.what());
}
} // namespace OpenMS
