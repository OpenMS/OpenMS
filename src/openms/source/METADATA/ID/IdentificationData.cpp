// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <random>
#include <set>
#include <type_traits>
#include <utility>

namespace OpenMS
{
namespace
{
  using ID = IdentificationData;
  static_assert(std::is_nothrow_swappable_v<std::vector<ID::Match>>);
  static_assert(std::is_nothrow_move_assignable_v<ID::Identification>);
  static_assert(std::is_nothrow_swappable_v<std::unique_ptr<ID::RunSettings>>);
  static_assert(std::is_nothrow_move_assignable_v<ID::Match>);
  [[noreturn]] void invalid(const std::string& message)
  { throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message, "IdentificationData"); }
  std::string describe(const ID::ScoreDefinition& definition)
  {
    std::string text = "'" + definition.name + "' (" + (definition.higher_better ? "higher" : "lower") + " is better, from "
                       + (definition.software.empty() ? std::string("unknown software") : definition.software)
                       + (definition.software_version.empty() ? std::string() : " " + definition.software_version);
    if (! definition.accession.empty()) text += ", " + definition.accession;
    if (! definition.calibration.empty() || ! definition.aggregation.empty() || ! definition.parameters.isMetaEmpty()
        || definition.scope != ID::ScoreScope::MATCH)
      text += ", with scope/parameter/calibration provenance";
    return text + ")";
  }
  std::string describe(const std::vector<ID::ScoreDefinition>& definitions)
  {
    std::string text = "[";
    for (const auto& definition : definitions)
      text += (text.size() > 1 ? "; " : "") + describe(definition);
    return text + "]";
  }
  bool isCatalog(const ID::Run& run)
  {
    const auto& processing = run.getSettings();
    return processing.metaValueExists("identification:catalog") && processing.getMetaValue("identification:catalog").toString() == "true";
  }
  /// Checks the dataset-wide PSM score contract over @p runs; returns the first configured run.
  const ID::Run* checkScoreContract(const std::vector<const ID::Run*>& runs, bool check_primary = true)
  {
    const ID::Run* expected = nullptr;
    for (const auto* run : runs)
    {
      const auto& definitions = run->getScoreDefinitions();
      const auto primary = run->getPrimaryScore();
      const bool catalog = isCatalog(*run);
      if (catalog && (! definitions.empty() || primary)) invalid("A sequence catalog cannot declare PSM scores");
      if (definitions.empty() && ! primary && (run->getNumberOfMatches() == 0 || catalog)) continue;
      if (check_primary && ! primary) invalid("Run '" + run->getIdentifier() + "' must select a primary PSM score");
      if (! expected)
      {
        expected = run;
        continue;
      }
      if (definitions != expected->getScoreDefinitions())
        invalid("PSM score definitions differ between runs '" + expected->getIdentifier() + "' " + describe(expected->getScoreDefinitions())
                + " and '" + run->getIdentifier() + "' " + describe(definitions)
                + ". Runs can only be combined with one ordered score schema, including the producing software and version. "
                  "Normalize the scores first, e.g. rescore each search with PercolatorAdapter or IDPosteriorErrorProbability, "
                  "or select a common score with IDScoreSwitcher");
      if (check_primary && primary->value != expected->getPrimaryScore()->value)
        invalid("Primary PSM score selection differs between runs '" + expected->getIdentifier() + "' and '" + run->getIdentifier() + "'");
    }
    return expected;
  }
  const ID::Run* checkScoreContract(const std::deque<ID::Run>& runs,
                                    const ID::Run* candidate = nullptr,
                                    const std::string* replacing_uuid = nullptr,
                                    bool check_primary = true)
  {
    std::vector<const ID::Run*> selected;
    selected.reserve(runs.size() + 1);
    for (const auto& run : runs)
      if (! replacing_uuid || run.getUuid() != *replacing_uuid) selected.push_back(&run);
    if (candidate) selected.push_back(candidate);
    return checkScoreContract(selected, check_primary);
  }
  UInt64 token()
  {
    static std::atomic<UInt64> next {1};
    UInt64 value = next.fetch_add(1, std::memory_order_relaxed);
    if (! value || value == std::numeric_limits<UInt64>::max()) invalid("Runtime handle space exhausted");
    return value;
  }
  /// A new state of the score columns of a run (never 0, which no run uses).
  UInt32 newTag()
  {
    static std::atomic<UInt32> next {1};
    UInt32 value = next.fetch_add(1, std::memory_order_relaxed);
    while (value == 0)
      value = next.fetch_add(1, std::memory_order_relaxed);
    return value;
  }
  std::optional<double> stored(double value)
  { return std::isnan(value) ? std::nullopt : std::optional<double>(value); }
  std::string uuid()
  {
    thread_local std::mt19937_64 engine([] {
      std::random_device source;
      std::seed_seq seed {source(), source(), source(), source(), source(), source(), source(), source()};
      return std::mt19937_64(seed);
    }());
    std::string result;
    const char* digits = "0123456789abcdef";
    for (Size i = 0; i < 32; ++i)
    {
      if (i == 8 || i == 12 || i == 16 || i == 20) result += '-';
      UInt64 digit = engine() & 15;
      if (i == 12) digit = 4;
      if (i == 16) digit = 8 | (digit & 3);
      result += digits[digit];
    }
    return result;
  }
  bool validUuid(const std::string& value)
  {
    if (value.size() != 36) return false;
    for (Size i = 0; i < value.size(); ++i)
    {
      if (i == 8 || i == 13 || i == 18 || i == 23)
      {
        if (value[i] != '-') return false;
      }
      else if (! ((value[i] >= '0' && value[i] <= '9') || (value[i] >= 'a' && value[i] <= 'f')))
        return false;
    }
    return true;
  }
  bool sameHypothesis(const ID::MatchData& first, const ID::MatchData& second)
  {
    const auto& a = first.details.value_or_default();
    const auto& b = second.details.value_or_default();
    return first.representation == second.representation && first.encoding == second.encoding && first.charge == second.charge
           && a.formula == b.formula && a.adduct == b.adduct;
  }
  void validNumber(const std::optional<double>& number)
  {
    if (number && ! std::isfinite(*number)) invalid("Scores and observation coordinates must be finite");
  }
  void validObservation(const ID::Observation& observation)
  {
    validNumber(observation.rt);
    validNumber(observation.mz);
  }
  UInt64 following(UInt64 id)
  {
    if (! id || id == std::numeric_limits<UInt64>::max()) invalid("ID is zero or its allocation counter would overflow");
    return id + 1;
  }
  struct CallbackGuard
  {
    bool& active;
    explicit CallbackGuard(bool& flag): active(flag)
    {
      if (active) invalid("Reentrant mutation during callback");
      active = true;
    }
    ~CallbackGuard()
    { active = false; }
  };
} // namespace

/// The score values of a run: one column per score definition and one row per match (Match::slot_); NaN is a missing
/// value. The tag names the state of the rows; matches and views record the tag they belong to.
struct ID::ScoreTable
{
  UInt32 tag = newTag();
  std::vector<std::vector<double>> columns;
};

std::optional<double> ID::ScoreView::operator()(const Match& match) const
{
  if (! table_ || match.tag_ != tag_ || table_->tag != tag_ || column_ >= table_->columns.size() || match.slot_ >= table_->columns[column_].size())
    invalid("Score view belongs to a different run, score schema or state of the run; bind it again");
  return stored(table_->columns[column_][match.slot_]);
}

ID::Run::Run(std::string identifier, MoleculeKind kind):
    identifier_(std::move(identifier)), uuid_(uuid()), kind_(kind), score_table_(std::make_shared<ScoreTable>())
{
  if (kind_ != MoleculeKind::PEPTIDE && kind_ != MoleculeKind::OLIGONUCLEOTIDE && kind_ != MoleculeKind::COMPOUND) invalid("Unknown molecule kind");
}

ID::Run::Run(const Run& other):
    identifier_(other.identifier_),
    uuid_(other.uuid_),
    kind_(other.kind_),
    settings_(std::make_unique<RunSettings>(*other.settings_)),
    databases_(other.databases_),
    sequences_(other.sequences_),
    sources_(other.sources_),
    scores_(other.scores_),
    score_owners_(other.score_owners_),
    primary_(other.primary_),
    score_table_(std::make_shared<ScoreTable>(*other.score_table_)),
    next_query_id_(other.next_query_id_),
    next_match_id_(other.next_match_id_),
    revision_(other.revision_),
    import_finalized_(other.import_finalized_),
    query_count_(other.query_count_),
    match_count_(other.match_count_),
    last_query_(other.last_query_),
    last_match_(other.last_match_)
{
  // Persistent identities and existing score/source handles remain valid in copies. The copy's score columns are
  // its own, in a new state, so views of the original cannot read them.
  score_table_->tag = newTag();
  for (auto& source : sources_)
    for (auto& query : source.identifications)
      for (auto& match : query.matches_)
        match.tag_ = score_table_->tag;
}

ID::Run::Run(Run&& other): Run()
{
  other.checkMutation_();
  swapData_(other);
}

void ID::Run::swapData_(Run& other) noexcept
{
  using std::swap;
  swap(identifier_, other.identifier_);
  swap(uuid_, other.uuid_);
  swap(kind_, other.kind_);
  swap(settings_, other.settings_);
  swap(databases_, other.databases_);
  swap(sequences_, other.sequences_);
  swap(sources_, other.sources_);
  swap(scores_, other.scores_);
  swap(score_owners_, other.score_owners_);
  swap(primary_, other.primary_);
  swap(score_table_, other.score_table_);
  swap(next_query_id_, other.next_query_id_);
  swap(next_match_id_, other.next_match_id_);
  swap(revision_, other.revision_);
  swap(import_finalized_, other.import_finalized_);
  swap(query_count_, other.query_count_);
  swap(match_count_, other.match_count_);
  // Atomics are not swappable; swapping is a mutation with exclusive access to both runs.
  query_index_built_.store(other.query_index_built_.exchange(query_index_built_.load()));
  match_index_built_.store(other.match_index_built_.exchange(match_index_built_.load()));
  swap(query_index_, other.query_index_);
  swap(match_index_, other.match_index_);
  swap(last_query_, other.last_query_);
  swap(last_match_, other.last_match_);
}

void ID::Run::checkMutation_() const
{
  if (callback_active_) invalid("Cannot modify a run from its filtering or transformation callback");
}
void ID::Run::checkScore_(ScoreId score) const
{
  if (score.value >= scores_.size() || score_owners_[score.value] != score.owner) invalid("Foreign or invalid score handle");
}
void ID::Run::checkRow_(const Match& match) const
{
  if (match.tag_ != score_table_->tag || (! score_table_->columns.empty() && match.slot_ >= score_table_->columns.front().size()))
    invalid("Match does not belong to this state of the run");
}
void ID::Run::renumberRows_() noexcept
{
  const UInt32 state = newTag();
  UInt32 row = 0;
  for (auto& source : sources_)
    for (auto& query : source.identifications)
      for (auto& match : query.matches_)
      {
        match.slot_ = row++;
        match.tag_ = state;
      }
  score_table_->tag = state;
}
void ID::Run::setSettings(const RunSettings& settings)
{
  checkMutation_();
  if (settings.metaValueExists("spectra_data")
      && ! (settings.getMetaValue("spectra_data").valueType() == DataValue::STRING_LIST && settings.getMetaValue("spectra_data").toStringList().empty()))
    invalid("The files of a run are its sources; its settings must not list them as 'spectra_data'");
  if (! settings.search.db.empty() || ! settings.search.db_version.empty() || ! settings.search.taxonomy.empty())
    invalid("The databases of a run are its own records (Run::addDatabase); the search settings must not name them");
  auto replacement = std::make_unique<RunSettings>(settings);
  settings_.swap(replacement);
  ++revision_;
}
ID::DatabaseId ID::Run::addDatabase(const Database& database)
{
  checkMutation_();
  const auto existing = std::find(databases_.begin(), databases_.end(), database);
  if (existing != databases_.end()) return {static_cast<UInt32>(existing - databases_.begin())};
  if (databases_.size() >= std::numeric_limits<UInt32>::max()) invalid("Too many databases");
  databases_.push_back(database);
  ++revision_;
  return {static_cast<UInt32>(databases_.size() - 1)};
}
ID::DatabaseId ID::Run::getDatabaseId(UInt32 index) const
{
  if (index >= databases_.size()) invalid("Invalid database index");
  return {index};
}
const ID::Database& ID::Run::getDatabase(DatabaseId database) const
{
  if (database.value >= databases_.size()) invalid("Unknown database of the run");
  return databases_[database.value];
}
ID::QualifiedAccession ID::Run::qualify(DatabaseId database, const std::string& accession) const
{ return {getDatabase(database).path, accession}; }
void ID::Run::setDatabaseSequences(std::optional<std::vector<DatabaseSequence>> sequences)
{
  checkMutation_();
  if (sequences)
  {
    std::set<std::pair<UInt32, std::string>> identities;
    for (const auto& sequence : *sequences)
    {
      if (sequence.database.value >= databases_.size()) invalid("Database sequence refers to an unknown database of the run");
      if (sequence.accession.empty() || ! identities.emplace(sequence.database.value, sequence.accession).second)
        invalid("Empty or duplicate database sequence accession");
    }
  }
  sequences_ = std::move(sequences);
  ++revision_;
}
ID::SourceId ID::Run::addSource(const SourceFile& source)
{
  checkMutation_();
  if (sources_.size() >= std::numeric_limits<UInt32>::max()) invalid("Too many sources");
  SourceId id {static_cast<UInt32>(sources_.size()), token()};
  sources_.push_back({id, source, {}});
  ++revision_;
  return id;
}
ID::SourceId ID::Run::getSourceId(UInt32 index) const
{
  if (index >= sources_.size()) invalid("Invalid source index");
  return sources_[index].id;
}
ID::ScoreId ID::Run::addScore(const ScoreDefinition& definition)
{
  checkMutation_();
  if (definition.name.empty()) invalid("A score definition needs a name");
  if (definition.scope < ScoreScope::MATCH || definition.scope > ScoreScope::OTHER) invalid("Invalid score scope");
  for (UInt32 i = 0; i < scores_.size(); ++i)
    if (scores_[i] == definition) return getScoreId(i);
  if (scores_.size() >= std::numeric_limits<UInt32>::max()) invalid("Too many score definitions");
  // Allocate before changing the logical schema, maintaining the strong guarantee.
  ScoreDefinition copy = definition;
  scores_.reserve(scores_.size() + 1);
  score_owners_.reserve(scores_.size() + 1);
  score_table_->columns.reserve(scores_.size() + 1);
  std::vector<double> column(match_count_, std::numeric_limits<double>::quiet_NaN());
  UInt64 owner = token();
  scores_.push_back(std::move(copy));
  score_owners_.push_back(owner);
  score_table_->columns.push_back(std::move(column));
  // A new schema is a new state of the columns: views bound before reject the matches.
  const UInt32 state = newTag();
  for (auto& source : sources_)
    for (auto& query : source.identifications)
      for (auto& match : query.matches_)
        match.tag_ = state;
  score_table_->tag = state;
  ++revision_;
  return getScoreId(static_cast<UInt32>(scores_.size() - 1));
}
ID::ScoreId ID::Run::getScoreId(UInt32 index) const
{
  if (index >= scores_.size()) invalid("Invalid score index");
  return {index, score_owners_[index]};
}
ID::ScoreId ID::Run::findScore(const ScoreDefinition& definition) const
{
  for (UInt32 i = 0; i < scores_.size(); ++i)
    if (scores_[i] == definition) return getScoreId(i);
  invalid("Unknown score definition");
}
const ID::ScoreDefinition& ID::Run::getScoreDefinition(ScoreId score) const
{
  checkScore_(score);
  return scores_[score.value];
}
ID::ScoreView ID::Run::bindScore(ScoreId score) const
{
  checkScore_(score);
  ScoreView view;
  view.table_ = score_table_;
  view.column_ = score.value;
  view.tag_ = score_table_->tag;
  view.definition_ = scores_[score.value];
  return view;
}
void ID::Run::setPrimaryScore(std::optional<ScoreId> score)
{
  checkMutation_();
  if (score)
  {
    checkScore_(*score);
    const auto& column = score_table_->columns[score->value];
    if (std::any_of(column.begin(), column.end(), [](double value) { return std::isnan(value); })) invalid("Primary score is missing on a candidate");
  }
  primary_ = score;
  ++revision_;
}
void ID::Run::validateMatch_(const MatchData& data, const std::vector<std::optional<double>>& values) const
{
  if (values.size() > scores_.size()) invalid("Too many score values");
  for (const auto& value : values)
    validNumber(value);
  if (primary_ && (primary_->value >= values.size() || ! values[primary_->value])) invalid("Primary score is missing on a candidate");
  validateMatchData_(data);
}
void ID::Run::validateMatchData_(const MatchData& data) const
{
  validNumber(data.calculated_mz);
  if (data.representation.empty()) invalid("Molecular representation must not be empty");
  bool compatible = kind_ == MoleculeKind::PEPTIDE ? data.encoding == Encoding::AA_SEQUENCE || data.encoding == Encoding::DATABASE_ID
                    : kind_ == MoleculeKind::OLIGONUCLEOTIDE
                      ? data.encoding == Encoding::NA_SEQUENCE || data.encoding == Encoding::DATABASE_ID
                      : data.encoding == Encoding::SMILES || data.encoding == Encoding::INCHI || data.encoding == Encoding::DATABASE_ID;
  if (! compatible) invalid("Molecular encoding is incompatible with the run kind");
  if (data.target_decoy < TargetDecoy::UNKNOWN || data.target_decoy > TargetDecoy::BOTH) invalid("Invalid target/decoy state");
  if (data.details && data.details->adduct && data.details->adduct->getCharge() != data.charge) invalid("Adduct charge disagrees with ion charge");
  if (kind_ == MoleculeKind::COMPOUND && ! data.sequence_evidence.empty()) invalid("Compound candidates cannot contain sequence evidence");
  for (const auto& evidence : data.sequence_evidence)
  {
    if (evidence.database.value >= databases_.size()) invalid("Sequence evidence refers to an unknown database of the run");
    if (evidence.start && evidence.end && *evidence.start > *evidence.end) invalid("Sequence evidence start exceeds end");
  }
}
void ID::Run::prepareLookupIndexes()
{
  checkMutation_();
  ensureQueryIndex_();
  ensureMatchIndex_();
}
void ID::Run::shrinkToFit()
{
  checkMutation_();
  // Lookup indexes and the cached last positions hold positions, not addresses, so they stay valid.
  databases_.shrink_to_fit();
  if (sequences_) sequences_->shrink_to_fit();
  for (auto& column : score_table_->columns)
    column.shrink_to_fit();
  sources_.shrink_to_fit();
  for (auto& source : sources_)
  {
    source.identifications.shrink_to_fit();
    for (auto& query : source.identifications)
      query.matches_.shrink_to_fit();
  }
}
void ID::Run::ensureQueryIndex_() const
{
  // Double-checked: concurrent const lookups build the index once and then read it lock-free.
  if (query_index_built_.load(std::memory_order_acquire)) return;
  std::lock_guard<std::mutex> lock(index_mutex_);
  if (query_index_built_.load(std::memory_order_relaxed)) return;
  std::unordered_map<UInt64, std::array<Size, 2>> index;
  index.reserve(query_count_);
  for (Size s = 0; s < sources_.size(); ++s)
    for (Size q = 0; q < sources_[s].identifications.size(); ++q)
      index.emplace(sources_[s].identifications[q].id_.value, std::array<Size, 2> {s, q});
  query_index_.swap(index);
  query_index_built_.store(true, std::memory_order_release);
}
void ID::Run::ensureMatchIndex_() const
{
  if (match_index_built_.load(std::memory_order_acquire)) return;
  std::lock_guard<std::mutex> lock(index_mutex_);
  if (match_index_built_.load(std::memory_order_relaxed)) return;
  std::unordered_map<UInt64, std::array<Size, 3>> index;
  index.reserve(match_count_);
  for (Size s = 0; s < sources_.size(); ++s)
    for (Size q = 0; q < sources_[s].identifications.size(); ++q)
      for (Size m = 0; m < sources_[s].identifications[q].matches_.size(); ++m)
        index.emplace(sources_[s].identifications[q].matches_[m].id_.value, std::array<Size, 3> {s, q, m});
  match_index_.swap(index);
  match_index_built_.store(true, std::memory_order_release);
}
void ID::Run::invalidateIndexes_()
{
  query_index_.clear();
  match_index_.clear();
  query_index_built_ = false;
  match_index_built_ = false;
  last_query_.reset();
  last_match_.reset();
}
const ID::Identification* ID::Run::findIdentification(QueryId id) const
{
  if (last_query_)
  {
    const auto& p = *last_query_;
    const auto& query = sources_[p[0]].identifications[p[1]];
    if (query.id_ == id) return &query;
  }
  ensureQueryIndex_();
  auto found = query_index_.find(id.value);
  return found == query_index_.end() ? nullptr : &sources_[found->second[0]].identifications[found->second[1]];
}
const ID::Match* ID::Run::findMatch(MatchId id) const
{
  if (last_match_)
  {
    const auto& p = *last_match_;
    const auto& match = sources_[p[0]].identifications[p[1]].matches_[p[2]];
    if (match.id_ == id) return &match;
  }
  ensureMatchIndex_();
  auto found = match_index_.find(id.value);
  return found == match_index_.end() ? nullptr : &sources_[found->second[0]].identifications[found->second[1]].matches_[found->second[2]];
}
const ID::Identification& ID::Run::getIdentification(QueryId id) const
{
  const auto* query = findIdentification(id);
  if (! query) invalid("Unknown query ID");
  return *query;
}
const ID::Match& ID::Run::getMatch(MatchId id) const
{
  const auto* match = findMatch(id);
  if (! match) invalid("Unknown match ID");
  return *match;
}
const ID::Identification& ID::Run::getIdentificationForMatch(MatchId id) const
{
  ensureMatchIndex_();
  auto found = match_index_.find(id.value);
  if (found == match_index_.end()) invalid("Unknown match ID");
  return sources_[found->second[0]].identifications[found->second[1]];
}
ID::Identification& ID::Run::query_(QueryId id)
{ return const_cast<Identification&>(getIdentification(id)); }
ID::Match& ID::Run::match_(MatchId id)
{ return const_cast<Match&>(getMatch(id)); }

ID::QueryId ID::Run::addIdentification(SourceId source, const Observation& observation)
{ return importIdentification(source, QueryId {next_query_id_}, observation); }
ID::QueryId ID::Run::importIdentification(SourceId source, QueryId id, Observation observation)
{
  checkMutation_();
  if (source.value >= sources_.size() || sources_[source.value].id != source) invalid("Foreign or invalid source handle");
  UInt64 next = following(id.value);
  validObservation(observation);
  if (id.value < next_query_id_ && (import_finalized_ || findIdentification(id))) invalid("Duplicate or historical query ID");
  Identification query;
  static_cast<Observation&>(query) = std::move(observation);
  query.id_ = id;
  auto& queries = sources_[source.value].identifications;
  std::array<Size, 2> position {source.value, queries.size()};
  // Drop an optional index on allocation failure; scientific values remain unchanged.
  if (query_index_built_)
  {
    try
    {
      query_index_.emplace(id.value, position);
    }
    catch (...)
    {
      query_index_.clear();
      query_index_built_ = false;
      throw;
    }
  }
  try
  {
    queries.push_back(std::move(query));
  }
  catch (...)
  {
    if (query_index_built_) query_index_.erase(id.value);
    throw;
  }
  ++query_count_;
  next_query_id_ = std::max(next_query_id_, next);
  last_query_ = position;
  ++revision_;
  return id;
}
ID::MatchId ID::Run::addMatch(QueryId query, const MatchData& data, const std::vector<std::optional<double>>& values)
{ return importMatch(query, MatchId {next_match_id_}, data, values); }
ID::MatchId ID::Run::importMatch(QueryId query_id, MatchId id, MatchData data, const std::vector<std::optional<double>>& values)
{
  checkMutation_();
  UInt64 next = following(id.value);
  validateMatch_(data, values);
  if (id.value < next_match_id_ && (import_finalized_ || findMatch(id))) invalid("Duplicate or historical match ID");
  auto& query = query_(query_id);
  if (match_count_ >= std::numeric_limits<UInt32>::max()) invalid("Too many matches in a run");
  // Details without a value are stored as none.
  if (data.details && data.details->empty()) data.details.reset();
  Match match;
  static_cast<MatchData&>(match) = std::move(data);
  match.id_ = id;
  match.tag_ = score_table_->tag;
  match.slot_ = static_cast<UInt32>(match_count_);
  auto& columns = score_table_->columns;
  for (auto& column : columns)
    column.reserve(match_count_ + 1);
  std::array<Size, 2> qp;
  if (last_query_ && sources_[(*last_query_)[0]].identifications[(*last_query_)[1]].id_ == query_id) qp = *last_query_;
  else
  {
    ensureQueryIndex_();
    qp = query_index_.at(query_id.value);
  }
  std::array<Size, 3> position {qp[0], qp[1], query.matches_.size()};
  if (match_index_built_)
  {
    try
    {
      match_index_.emplace(id.value, position);
    }
    catch (...)
    {
      match_index_.clear();
      match_index_built_ = false;
      throw;
    }
  }
  try
  {
    query.matches_.push_back(std::move(match));
  }
  catch (...)
  {
    if (match_index_built_) match_index_.erase(id.value);
    throw;
  }
  // Capacity was reserved above, so appending the row cannot throw.
  for (Size i = 0; i < columns.size(); ++i)
    columns[i].push_back(i < values.size() && values[i] ? *values[i] : std::numeric_limits<double>::quiet_NaN());
  ++match_count_;
  next_match_id_ = std::max(next_match_id_, next);
  last_match_ = position;
  ++revision_;
  return id;
}
std::optional<double> ID::Run::getScore(MatchId match, ScoreId score) const
{
  checkScore_(score);
  return stored(score_table_->columns[score.value][getMatch(match).slot_]);
}
std::vector<std::optional<double>> ID::Run::getScores(const Match& match) const
{
  checkRow_(match);
  std::vector<std::optional<double>> result;
  result.reserve(scores_.size());
  for (const auto& column : score_table_->columns)
    result.push_back(stored(column[match.slot_]));
  return result;
}
std::vector<std::optional<double>> ID::Run::getScores(MatchId match) const
{ return getScores(getMatch(match)); }
void ID::Run::setScore(MatchId match, ScoreId score, std::optional<double> value)
{
  checkMutation_();
  checkScore_(score);
  validNumber(value);
  if (! value && primary_ == score) invalid("Cannot clear a primary score");
  score_table_->columns[score.value][match_(match).slot_] = value.value_or(std::numeric_limits<double>::quiet_NaN());
  ++revision_;
}
void ID::Run::setSelectedMatch(QueryId query_id, std::optional<MatchId> selected)
{
  checkMutation_();
  auto& query = query_(query_id);
  if (selected && std::none_of(query.matches_.begin(), query.matches_.end(), [&](const Match& match) { return match.id_ == selected; }))
    invalid("Selected match is not owned by the query");
  query.selected_ = selected;
  ++revision_;
}
void ID::Run::replaceObservation(QueryId query, const Observation& observation)
{
  checkMutation_();
  validObservation(observation);
  Observation copy(observation);
  static_cast<Observation&>(query_(query)) = std::move(copy);
  ++revision_;
}
void ID::Run::replaceMatch(MatchId match, const MatchData& data)
{
  checkMutation_();
  auto& target = match_(match);
  if (! sameHypothesis(target, data)) invalid("Changing a molecular/ion hypothesis requires explicitly replacing its scores");
  replaceMatch(match, data, getScores(target));
}
void ID::Run::replaceMatch(MatchId match, const MatchData& data, const std::vector<std::optional<double>>& values)
{
  checkMutation_();
  auto& target = match_(match);
  validateMatch_(data, values);
  Match replacement(target);
  static_cast<MatchData&>(replacement) = data;
  if (replacement.details && replacement.details->empty()) replacement.details.reset();
  // The match keeps its row; its values are replaced after everything that may throw.
  target = std::move(replacement);
  for (Size i = 0; i < score_table_->columns.size(); ++i)
    score_table_->columns[i][target.slot_] = i < values.size() && values[i] ? *values[i] : std::numeric_limits<double>::quiet_NaN();
  ++revision_;
}
Size ID::Run::filterMatches(const std::function<bool(const Match&)>& keep, bool keep_empty_queries)
{
  checkMutation_();
  // Evaluate first: exceptions and attempted reentrant edits cannot leave a partial filter.
  std::vector<bool> decisions;
  decisions.reserve(match_count_);
  {
    CallbackGuard guard(callback_active_);
    for (const auto& source : sources_)
      for (const auto& query : source.identifications)
        for (const auto& match : query.matches_)
          decisions.push_back(keep(match));
  }
  // The score columns of the retained matches, in their order, are built before anything changes.
  const Size retained = static_cast<Size>(std::count(decisions.begin(), decisions.end(), true));
  std::vector<std::vector<double>> columns(score_table_->columns.size());
  for (auto& column : columns)
    column.reserve(retained);
  {
    Size at = 0;
    for (const auto& source : sources_)
      for (const auto& query : source.identifications)
        for (const auto& match : query.matches_)
          if (decisions[at++])
            for (Size i = 0; i < columns.size(); ++i)
              columns[i].push_back(score_table_->columns[i][match.slot_]);
  }
  // Commit: matches move without throwing; erasing queries and swapping vectors cannot throw.
  Size at = 0, removed = 0;
  for (auto& source : sources_)
    for (auto& query : source.identifications)
    {
      auto end = std::remove_if(query.matches_.begin(), query.matches_.end(), [&](const Match&) {
        const bool erase = ! decisions[at++];
        removed += erase;
        return erase;
      });
      query.matches_.erase(end, query.matches_.end());
    }
  for (auto& source : sources_)
    for (auto& query : source.identifications)
    {
      if (query.selected_
          && std::none_of(query.matches_.begin(), query.matches_.end(), [&](const Match& match) { return match.id_ == query.selected_; }))
        query.selected_.reset();
    }
  if (! keep_empty_queries)
    for (auto& source : sources_)
    {
      Size before = source.identifications.size();
      std::erase_if(source.identifications, [](const Identification& query) { return query.matches_.empty(); });
      query_count_ -= before - source.identifications.size();
    }
  match_count_ -= removed;
  score_table_->columns.swap(columns);
  renumberRows_();
  import_finalized_ = true;
  invalidateIndexes_();
  ++revision_;
  return removed;
}
Size ID::Run::eraseMatches(const std::function<bool(const Match&)>& remove, bool keep_empty_queries)
{
  return filterMatches([&](const Match& match) { return ! remove(match); }, keep_empty_queries);
}
Size ID::Run::eraseIdentifications(const std::function<bool(const Identification&)>& remove)
{
  checkMutation_();
  std::vector<bool> decisions;
  decisions.reserve(query_count_);
  {
    CallbackGuard guard(callback_active_);
    for (const auto& source : sources_)
      for (const auto& query : source.identifications)
        decisions.push_back(remove(query));
  }
  const Size erased = static_cast<Size>(std::count(decisions.begin(), decisions.end(), true));
  if (erased == 0) return 0;
  std::set<MatchId> matches;
  {
    Size at = 0;
    for (const auto& source : sources_)
      for (const auto& query : source.identifications)
        if (decisions[at++])
          for (const auto& match : query.matches_)
            matches.insert(match.id_);
  }
  // Removing the matches first rebuilds the score columns; it changes nothing if it throws.
  // Afterwards the identifications to remove are empty, and erasing them cannot throw.
  filterMatches([&](const Match& match) { return ! matches.contains(match.id_); }, true);
  Size offset = 0;
  for (auto& source : sources_)
  {
    const Identification* first = source.identifications.data();
    const auto end = std::remove_if(source.identifications.begin(), source.identifications.end(),
                                    [&](const Identification& query) { return decisions[offset + static_cast<Size>(&query - first)]; });
    offset += source.identifications.size();
    source.identifications.erase(end, source.identifications.end());
  }
  query_count_ -= erased;
  invalidateIndexes_();
  ++revision_;
  return erased;
}
Size ID::Run::retainBest(ScoreId score, bool keep_ties, bool keep_empty_queries)
{
  checkMutation_();
  checkScore_(score);
  std::vector<bool> retained;
  retained.reserve(match_count_);
  for (const auto& source : sources_)
    for (const auto& query : source.identifications)
    {
      const auto& column = score_table_->columns[score.value];
      std::optional<double> best;
      for (const auto& match : query.matches_)
      {
        double value = column[match.slot_];
        if (std::isnan(value)) continue;
        if (! best || (scores_[score.value].higher_better ? value > *best : value < *best)) best = value;
      }
      bool taken = false;
      for (const auto& match : query.matches_)
      {
        bool keep = best && column[match.slot_] == *best && (keep_ties || ! taken);
        retained.push_back(keep);
        taken = taken || keep;
      }
    }
  Size at = 0;
  return filterMatches([&](const Match&) { return retained[at++]; }, keep_empty_queries);
}
void ID::Run::transformMatches(const std::function<void(MatchData&)>& transform)
{
  checkMutation_();
  std::vector<MatchData> replacements;
  replacements.reserve(match_count_);
  {
    CallbackGuard guard(callback_active_);
    for (const auto& source : sources_)
      for (const auto& query : source.identifications)
        for (const auto& match : query.matches_)
        {
          MatchData data(match);
          transform(data);
          if (! sameHypothesis(match, data)) invalid("Transformation cannot change a hypothesis without explicit replacement scores");
          validateMatch_(data, getScores(match));
          if (data.details && data.details->empty()) data.details.reset();
          replacements.push_back(std::move(data));
        }
  }
  // Moving the payloads cannot throw; rows and scores stay.
  Size at = 0;
  for (auto& source : sources_)
    for (auto& query : source.identifications)
      for (auto& match : query.matches_)
        static_cast<MatchData&>(match) = std::move(replacements[at++]);
  ++revision_;
}
Size ID::Run::getNumberOfIdentifications() const
{ return query_count_; }
Size ID::Run::getNumberOfMatches() const
{ return match_count_; }
void ID::Run::restoreIdentity(const std::string& identity, UInt64 next_query, UInt64 next_match)
{
  checkMutation_();
  if (! validUuid(identity)) invalid("Invalid run UUID");
  if (! next_query || ! next_match || next_query < next_query_id_ || next_match < next_match_id_)
    invalid("ID allocation counters cannot move backwards");
  uuid_ = identity;
  next_query_id_ = next_query;
  next_match_id_ = next_match;
  import_finalized_ = true;
  ++revision_;
}
void ID::Run::reserveMatchId(MatchId id)
{
  checkMutation_();
  next_match_id_ = std::max(next_match_id_, following(id.value));
  import_finalized_ = true;
  ++revision_;
}
void ID::Run::validate() const
{
  if (! validUuid(uuid_) || ! next_query_id_ || ! next_match_id_) invalid("Invalid run identity");
  Size nq = 0, nm = 0;
  std::vector<UInt64> queries, matches;
  queries.reserve(query_count_);
  matches.reserve(match_count_);
  const auto& columns = score_table_->columns;
  if (columns.size() != scores_.size()) invalid("Inconsistent score column count");
  for (const auto& column : columns)
  {
    if (column.size() != match_count_) invalid("Inconsistent score column length");
    for (double score : column)
      if (! std::isnan(score) && ! std::isfinite(score)) invalid("Scores must be finite or missing");
  }
  std::vector<bool> rows(match_count_, false);
  for (const auto& source : sources_)
    for (const auto& query : source.identifications)
    {
      ++nq;
      validObservation(query);
      if (! query.id_.value || query.id_.value >= next_query_id_) invalid("Invalid or duplicate query ID");
      queries.push_back(query.id_.value);
      bool selected_found = ! query.selected_;
      for (const auto& match : query.matches_)
      {
        ++nm;
        validateMatchData_(match);
        if (match.tag_ != score_table_->tag || match.slot_ >= rows.size() || rows[match.slot_]) invalid("Inconsistent score rows");
        rows[match.slot_] = true;
        if (primary_ && std::isnan(columns[primary_->value][match.slot_])) invalid("Primary score is missing on a candidate");
        if (! match.id_.value || match.id_.value >= next_match_id_) invalid("Invalid or duplicate match ID");
        matches.push_back(match.id_.value);
        if (query.selected_ == match.id_) selected_found = true;
      }
      if (! selected_found) invalid("Selected match is not in its query");
    }
  if (nq != query_count_ || nm != match_count_) invalid("Inconsistent record counts");
  // Most imports retain increasing IDs. Check these linearly; arbitrary scientific
  // ordering only requires sorting compact ID copies, never moving model records.
  for (auto* ids : {&queries, &matches})
  {
    if (! std::is_sorted(ids->begin(), ids->end())) std::sort(ids->begin(), ids->end());
    if (std::adjacent_find(ids->begin(), ids->end()) != ids->end()) invalid("Duplicate query or match ID");
  }
}

ID::IdentificationData(const IdentificationData& other): runs_(other.runs_), inference_(other.inference_)
{
}
ID::IdentificationData(IdentificationData&& other)
{
  other.checkMutation_();
  runs_ = std::move(other.runs_);
  inference_ = std::move(other.inference_);
}
ID& ID::operator=(const IdentificationData& other)
{
  checkMutation_();
  if (this != &other)
  {
    ID copy(other);
    swap(copy);
  }
  return *this;
}
ID& ID::operator=(IdentificationData&& other)
{
  checkMutation_();
  other.checkMutation_();
  if (this != &other)
  {
    ID moved(std::move(other));
    swap(moved);
  }
  return *this;
}
void ID::checkMutation_() const
{
  if (callback_active_) invalid("Cannot modify a dataset from its transformation callback");
  for (const auto& run : runs_)
    run.checkMutation_();
}
ID::Run& ID::addRun(const std::string& identifier, MoleculeKind kind)
{ return addRun(Run(identifier, kind)); }
ID::Run& ID::addRun(Run run)
{
  checkMutation_();
  for (const auto& existing : runs_)
    if (existing.uuid_ == run.uuid_ || existing.identifier_ == run.identifier_) invalid("Duplicate run UUID or display identifier");
  run.validate();
  checkScoreContract(runs_, &run);
  for (const auto& result : inference_)
  {
    for (const auto& input : result.inputs)
      if (input.run_uuid == run.uuid_) invalid("An independent run cannot reuse an unresolved provenance UUID");
  }
  runs_.push_back(std::move(run));
  return runs_.back();
}
ID::Run& ID::getRun(const std::string& identifier)
{
  auto found = std::find_if(runs_.begin(), runs_.end(), [&](const Run& run) { return run.identifier_ == identifier; });
  if (found == runs_.end()) invalid("Unknown run identifier");
  return *found;
}
const ID::Run& ID::getRun(const std::string& identifier) const
{ return const_cast<ID*>(this)->getRun(identifier); }
ID::Run* ID::findRunByUuid(const std::string& identity)
{
  auto found = std::find_if(runs_.begin(), runs_.end(), [&](const Run& run) { return run.uuid_ == identity; });
  return found == runs_.end() ? nullptr : &*found;
}
const ID::Run* ID::findRunByUuid(const std::string& identity) const
{ return const_cast<ID*>(this)->findRunByUuid(identity); }
void ID::addInferenceResult(InferenceResult result)
{
  checkMutation_();
  for (const auto& existing : inference_)
    if (existing.identifier == result.identifier) invalid("Duplicate inference result identifier");
  for (const auto& input : result.inputs)
  {
    if (! validUuid(input.run_uuid)) invalid("Inference input needs a run UUID");
  }
  inference_.push_back(std::move(result));
}
bool ID::Run::operator==(const Run& other) const
{
  if (uuid_ != other.uuid_ || identifier_ != other.identifier_ || kind_ != other.kind_ || *settings_ != *other.settings_
      || databases_ != other.databases_ || sequences_ != other.sequences_ || scores_ != other.scores_ || next_query_id_ != other.next_query_id_ || next_match_id_ != other.next_match_id_
      || sources_.size() != other.sources_.size())
    return false;
  const auto primary = primary_ ? std::optional<UInt32>(primary_->value) : std::nullopt;
  const auto other_primary = other.primary_ ? std::optional<UInt32>(other.primary_->value) : std::nullopt;
  if (primary != other_primary) return false;
  for (Size i = 0; i < sources_.size(); ++i)
  {
    const auto& source = sources_[i];
    const auto& rhs = other.sources_[i];
    if (source.id.value != rhs.id.value || source.file != rhs.file || source.identifications.size() != rhs.identifications.size()) return false;
    for (Size q = 0; q < source.identifications.size(); ++q)
    {
      const auto& query = source.identifications[q];
      const auto& rq = rhs.identifications[q];
      if (query.getId() != rq.getId() || query.getObservation() != rq.getObservation() || query.getSelectedMatch() != rq.getSelectedMatch()
          || query.getMatches().size() != rq.getMatches().size())
        return false;
      for (Size m = 0; m < query.getMatches().size(); ++m)
      {
        const auto& match = query.getMatches()[m];
        const auto& rm = rq.getMatches()[m];
        if (match.getId() != rm.getId() || match.getData() != rm.getData()) return false;
        for (Size j = 0; j < scores_.size(); ++j)
        {
          double left = score_table_->columns[j][match.slot_], right = other.score_table_->columns[j][rm.slot_];
          if (left != right && ! (std::isnan(left) && std::isnan(right))) return false;
        }
      }
    }
  }
  return true;
}
bool ID::operator==(const IdentificationData& other) const
{ return runs_ == other.runs_ && inference_ == other.inference_; }
void ID::clear()
{
  checkMutation_();
  runs_.clear();
  inference_.clear();
}
std::set<std::string> ID::MatchData::extractProteinAccessionsSet() const
{
  std::set<std::string> accessions;
  for (const auto& evidence : sequence_evidence)
    if (! evidence.accession.empty()) accessions.insert(evidence.accession);
  return accessions;
}

const ID::Match* ID::QueryMatches::getBestMatch() const
{
  const auto primary = run ? run->getPrimaryScore() : std::nullopt;
  if (! primary) return nullptr;
  const bool higher_better = run->getScoreDefinition(*primary).higher_better;
  const Match* best = nullptr;
  double best_score = 0.0;
  for (const auto* match : matches)
  {
    const auto score = run->getScore(match->getId(), *primary);
    if (! score || std::isnan(*score)) continue;
    if (best && ! (higher_better ? *score > best_score : *score < best_score)) continue;
    best = match;
    best_score = *score;
  }
  return best;
}
std::vector<ID::QueryMatches> ID::resolveLinks(const std::set<QueryReference>& queries, const std::set<MatchReference>& matches) const
{
  const auto missing = [](const std::string& what) {
    throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Link refers to a missing " + what);
  };
  std::map<QueryReference, QueryMatches> linked;
  for (const auto& reference : queries)
  {
    const auto* run = findRunByUuid(reference.run_uuid);
    const auto* query = run ? run->findIdentification(reference.query) : nullptr;
    if (! query) missing("identification");
    linked[reference] = {run, query, {}};
  }
  for (const auto& reference : matches)
  {
    const auto* run = findRunByUuid(reference.run_uuid);
    const auto* match = run ? run->findMatch(reference.match) : nullptr;
    if (! match) missing("match");
    const auto& query = run->getIdentificationForMatch(reference.match);
    auto& entry = linked[{reference.run_uuid, query.getId()}];
    entry.run = run;
    entry.query = &query;
    entry.matches.push_back(match);
  }
  std::map<const Run*, Size> positions;
  for (const auto& run : runs_)
    positions.emplace(&run, positions.size());
  std::vector<QueryMatches> result;
  result.reserve(linked.size());
  for (auto& [reference, entry] : linked)
  {
    // The matches of an identification are contiguous, in their order.
    std::sort(entry.matches.begin(), entry.matches.end(), std::less<const Match*>());
    result.push_back(std::move(entry));
  }
  std::stable_sort(result.begin(), result.end(), [&](const QueryMatches& a, const QueryMatches& b) {
    return std::make_pair(a.query->getId(), positions[a.run]) < std::make_pair(b.query->getId(), positions[b.run]);
  });
  return result;
}
std::vector<ID::QueryMatches> ID::getUnlinked(const std::set<QueryReference>& queries, const std::set<MatchReference>& matches) const
{
  std::vector<QueryMatches> order;
  for (const auto& run : runs_)
    for (const auto& source : run.getSources())
      for (const auto& query : source.identifications)
        order.push_back({&run, &query, {}});
  std::stable_sort(order.begin(), order.end(), [](const QueryMatches& a, const QueryMatches& b) { return a.query->getId() < b.query->getId(); });
  // Runs that share query IDs keep their identifications together (the order of an export).
  if (std::adjacent_find(order.begin(), order.end(), [](const auto& a, const auto& b) { return a.query->getId() == b.query->getId(); }) != order.end())
  {
    std::map<const Run*, Size> positions;
    for (const auto& run : runs_)
      positions.emplace(&run, positions.size());
    std::stable_sort(order.begin(), order.end(), [&](const QueryMatches& a, const QueryMatches& b) {
      return std::make_pair(positions[a.run], a.query->getId()) < std::make_pair(positions[b.run], b.query->getId());
    });
  }
  std::vector<QueryMatches> result;
  for (auto& entry : order)
  {
    const auto& uuid = entry.run->getUuid();
    for (const auto& match : entry.query->getMatches())
      if (! matches.contains({uuid, match.getId()})) entry.matches.push_back(&match);
    if (! entry.matches.empty() || (entry.query->getMatches().empty() && ! queries.contains({uuid, entry.query->getId()})))
      result.push_back(std::move(entry));
  }
  return result;
}
void ID::merge(const IdentificationData& other)
{
  checkMutation_();
  validate();
  // Every run and result of a dataset is already present in itself, with equal values.
  if (this == &other) return;
  other.validate();
  // Stage only the incoming runs and results and validate them against this dataset. The
  // commit appends to the run deque, which keeps references to existing runs valid.
  std::set<std::string> run_identifiers;
  for (const auto& run : runs_)
    run_identifiers.insert(run.identifier_);
  std::vector<Run> staged_runs;
  staged_runs.reserve(other.runs_.size());
  for (const auto& run : other.runs_)
  {
    if (const auto* existing = findRunByUuid(run.getUuid()))
    {
      auto comparable = run;
      comparable.identifier_ = existing->identifier_;
      if (*existing != comparable) invalid("Cannot merge conflicting values for the same run UUID");
      continue;
    }
    for (const auto& result : inference_)
      for (const auto& input : result.inputs)
        if (input.run_uuid == run.getUuid()) invalid("An independent run cannot reuse an unresolved provenance UUID");
    auto copy = run;
    const auto original = copy.identifier_;
    Size suffix = 2;
    while (run_identifiers.contains(copy.identifier_))
      copy.identifier_ = original + "#" + std::to_string(suffix++);
    run_identifiers.insert(copy.identifier_);
    staged_runs.push_back(std::move(copy));
  }
  {
    std::vector<const Run*> combined;
    combined.reserve(runs_.size() + staged_runs.size());
    for (const auto& run : runs_)
      combined.push_back(&run);
    for (const auto& run : staged_runs)
      combined.push_back(&run);
    checkScoreContract(combined);
  }
  const auto find_identifier = [&](const std::string& uuid) -> const std::string* {
    if (const auto* run = findRunByUuid(uuid)) return &run->identifier_;
    for (const auto& run : staged_runs)
      if (run.uuid_ == uuid) return &run.identifier_;
    return nullptr;
  };
  std::set<std::string> result_identifiers;
  for (const auto& result : inference_)
    result_identifiers.insert(result.identifier);
  std::vector<InferenceResult> staged_results;
  for (const auto& result : other.inference_)
  {
    auto copy = result;
    // Inputs follow renamed runs; compare after mapping so a repeated merge adds nothing.
    for (auto& input : copy.inputs)
    {
      if (! validUuid(input.run_uuid)) invalid("Inference input needs a run UUID");
      if (const auto* identifier = find_identifier(input.run_uuid)) input.run_identifier = *identifier;
    }
    const auto existing = std::find_if(inference_.begin(), inference_.end(), [&](const auto& r) { return r.identifier == copy.identifier; });
    if (existing != inference_.end() && *existing == copy) continue;
    const auto original = copy.identifier;
    Size suffix = 2;
    while (result_identifiers.contains(copy.identifier))
      copy.identifier = original + "#" + std::to_string(suffix++);
    result_identifiers.insert(copy.identifier);
    staged_results.push_back(std::move(copy));
  }
  // Commit. If an allocation fails part-way, remove what was appended so nothing changes.
  const Size old_runs = runs_.size();
  const Size old_results = inference_.size();
  try
  {
    for (auto& run : staged_runs)
      runs_.push_back(std::move(run));
    inference_.reserve(inference_.size() + staged_results.size());
    for (auto& result : staged_results)
      inference_.push_back(std::move(result));
  }
  catch (...)
  {
    while (runs_.size() > old_runs)
      runs_.pop_back();
    while (inference_.size() > old_results)
      inference_.pop_back();
    throw;
  }
}
void ID::clearInferenceResults()
{
  checkMutation_();
  inference_.clear();
}
Size ID::filterMatches(const std::function<bool(const Match&)>& keep, InferencePolicy policy, bool keep_empty_queries)
{
  checkMutation_();
  if (policy != InferencePolicy::PRESERVE && policy != InferencePolicy::DISCARD) invalid("Invalid inference policy");
  // Copy runs for an atomic operation across runs; streaming filtering avoids this owning cost.
  ID replacement(*this);
  Size removed = 0;
  {
    CallbackGuard guard(callback_active_);
    std::vector<std::unique_ptr<CallbackGuard>> guards;
    for (auto& run : runs_)
      guards.push_back(std::make_unique<CallbackGuard>(run.callback_active_));
    for (auto& run : replacement.runs_)
      removed += run.filterMatches(keep, keep_empty_queries);
  }
  if (policy == InferencePolicy::DISCARD) replacement.inference_.clear();
  for (Size i = 0; i < runs_.size(); ++i)
    runs_[i].swapData_(replacement.runs_[i]);
  inference_.swap(replacement.inference_);
  return removed;
}
const std::vector<ID::ScoreDefinition>& ID::getScoreDefinitions() const
{
  static const std::vector<ScoreDefinition> empty;
  const auto* run = checkScoreContract(runs_);
  return run ? run->getScoreDefinitions() : empty;
}
std::optional<ID::ScoreDefinition> ID::getPrimaryScoreDefinition() const
{
  checkScoreContract(runs_);
  for (const auto& run : runs_)
    if (run.primary_) return run.scores_[run.primary_->value];
  return std::nullopt;
}
void ID::setPrimaryScore(const ScoreDefinition& definition)
{
  checkMutation_();
  checkScoreContract(runs_, nullptr, nullptr, false);
  std::vector<std::pair<Run*, ScoreId>> selections;
  for (auto& run : runs_)
  {
    run.checkMutation_();
    // Same exemptions as the score contract: unconfigured empty runs and scoreless catalogs.
    if (! run.primary_ && run.scores_.empty() && (run.match_count_ == 0 || isCatalog(run))) continue;
    const auto score = run.findScore(definition);
    const auto& column = run.score_table_->columns[score.value];
    if (std::any_of(column.begin(), column.end(), [](double value) { return std::isnan(value); }))
      invalid("Primary score is missing on a candidate in run '" + run.identifier_ + "'");
    selections.emplace_back(&run, score);
  }
  // Every potentially throwing operation precedes the commit.
  for (auto& [run, score] : selections)
  {
    run->primary_ = score;
    ++run->revision_;
  }
}
void ID::validate() const
{
  checkScoreContract(runs_);
  std::set<std::string> identities, identifiers;
  for (const auto& run : runs_)
  {
    run.validate();
    if (! identities.insert(run.uuid_).second || ! identifiers.insert(run.identifier_).second) invalid("Duplicate run identity");
  }
  for (const auto& result : inference_)
  {
    for (const auto& input : result.inputs)
    {
      if (! validUuid(input.run_uuid)) invalid("Invalid inference input");
    }
  }
}
void ID::swap(IdentificationData& other)
{
  checkMutation_();
  other.checkMutation_();
  runs_.swap(other.runs_);
  inference_.swap(other.inference_);
}
} // namespace OpenMS
