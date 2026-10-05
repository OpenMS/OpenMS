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
  static_assert(std::is_nothrow_move_assignable_v<ID::Match>);
  static_assert(std::is_nothrow_move_assignable_v<ID::Identification>);
  static_assert(std::is_nothrow_swappable_v<ProteinIdentification>);
  [[noreturn]] void invalid(const std::string& message)
  { throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message, "IdentificationData"); }
  const ID::Run* checkScoreContract(const std::deque<ID::Run>& runs,
                                    const ID::Run* candidate = nullptr,
                                    const std::string* replacing_uuid = nullptr,
                                    bool check_primary = true)
  {
    const ID::Run* expected = nullptr;
    const auto check = [&](const ID::Run& run) {
      const auto& definitions = run.getScoreDefinitions();
      const auto primary = run.getPrimaryScore();
      if (definitions.empty() && ! primary && run.getNumberOfMatches() == 0) return;
      if (check_primary && ! primary) invalid("Run '" + run.getIdentifier() + "' must select a primary PSM score");
      if (! expected) expected = &run;
      else
      {
        if (definitions != expected->getScoreDefinitions())
          invalid("Ordered PSM score schema differs between runs '" + expected->getIdentifier() + "' and '" + run.getIdentifier()
                  + "'. Normalize score definitions and column order before combining runs");
        if (check_primary && primary->value != expected->getPrimaryScore()->value)
          invalid("Primary PSM score selection differs between runs '" + expected->getIdentifier() + "' and '" + run.getIdentifier() + "'");
      }
    };
    for (const auto& run : runs)
      if (!replacing_uuid || run.getUuid() != *replacing_uuid) check(run);
    if (candidate) check(*candidate);
    return expected;
  }
  UInt64 token()
  {
    static std::atomic<UInt64> next {1};
    UInt64 value = next.fetch_add(1, std::memory_order_relaxed);
    if (! value || value == std::numeric_limits<UInt64>::max()) invalid("Runtime handle space exhausted");
    return value;
  }
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
    return first.representation == second.representation && first.encoding == second.encoding && first.charge == second.charge
           && first.formula == second.formula && first.adduct == second.adduct;
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

std::vector<std::optional<double>> ID::Match::getScores() const
{
  std::vector<std::optional<double>> result;
  result.reserve(scores_.size());
  for (double value : scores_)
    result.push_back(std::isnan(value) ? std::nullopt : std::optional<double>(value));
  return result;
}

std::optional<double> ID::ScoreView::operator()(const Match& match) const
{
  if (match.schema_token_ != schema_token_ || id_.value >= match.scores_.size()) invalid("Score view belongs to a different score schema");
  double value = match.scores_[id_.value];
  return std::isnan(value) ? std::nullopt : std::optional<double>(value);
}

ID::Run::Run(std::string identifier, MoleculeKind kind): identifier_(std::move(identifier)), uuid_(uuid()), kind_(kind), schema_token_(token())
{
  if (kind_ != MoleculeKind::PEPTIDE && kind_ != MoleculeKind::OLIGONUCLEOTIDE && kind_ != MoleculeKind::COMPOUND) invalid("Unknown molecule kind");
}

ID::Run::Run(const Run& other):
    identifier_(other.identifier_),
    uuid_(other.uuid_),
    kind_(other.kind_),
    processing_(other.processing_),
    parents_(other.parents_),
    sources_(other.sources_),
    scores_(other.scores_),
    score_owners_(other.score_owners_),
    primary_(other.primary_),
    schema_token_(other.schema_token_),
    next_query_id_(other.next_query_id_),
    next_match_id_(other.next_match_id_),
    import_finalized_(other.import_finalized_),
    query_count_(other.query_count_),
    match_count_(other.match_count_),
    last_query_(other.last_query_),
    last_match_(other.last_match_)
{
  // Persistent identities and existing score/source handles remain valid in copies.
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
  swap(processing_, other.processing_);
  swap(parents_, other.parents_);
  swap(sources_, other.sources_);
  swap(scores_, other.scores_);
  swap(score_owners_, other.score_owners_);
  swap(primary_, other.primary_);
  swap(schema_token_, other.schema_token_);
  swap(next_query_id_, other.next_query_id_);
  swap(next_match_id_, other.next_match_id_);
  swap(import_finalized_, other.import_finalized_);
  swap(query_count_, other.query_count_);
  swap(match_count_, other.match_count_);
  swap(query_index_built_, other.query_index_built_);
  swap(match_index_built_, other.match_index_built_);
  swap(query_index_, other.query_index_);
  swap(match_index_, other.match_index_);
  swap(last_query_, other.last_query_);
  swap(last_match_, other.last_match_);
}

ID::Run& ID::Run::operator=(const Run& other)
{
  checkMutation_();
  if (this != &other)
  {
    Run copy(other);
    if (uuid_ == copy.uuid_)
    {
      copy.next_query_id_ = std::max(next_query_id_, copy.next_query_id_);
      copy.next_match_id_ = std::max(next_match_id_, copy.next_match_id_);
      copy.import_finalized_ = import_finalized_ || copy.import_finalized_;
    }
    swapData_(copy);
  }
  return *this;
}
ID::Run& ID::Run::operator=(Run&& other)
{
  checkMutation_();
  other.checkMutation_();
  if (this != &other)
  {
    Run moved(std::move(other));
    if (uuid_ == moved.uuid_)
    {
      moved.next_query_id_ = std::max(next_query_id_, moved.next_query_id_);
      moved.next_match_id_ = std::max(next_match_id_, moved.next_match_id_);
      moved.import_finalized_ = import_finalized_ || moved.import_finalized_;
    }
    swapData_(moved);
  }
  return *this;
}
void ID::Run::checkMutation_() const
{
  if (callback_active_) invalid("Cannot modify a run from its filtering or transformation callback");
}
void ID::Run::checkScore_(ScoreId score) const
{
  if (score.value >= scores_.size() || score_owners_[score.value] != score.owner) invalid("Foreign or invalid score handle");
}
void ID::Run::setProcessingMetadata(const ProteinIdentification& metadata)
{
  checkMutation_();
  processing_ = metadata;
}
void ID::Run::setParents(std::optional<std::vector<ParentRecord>> parents)
{
  checkMutation_();
  if (parents)
  {
    std::set<QualifiedAccession> identities;
    for (const auto& parent : *parents)
    {
      if (parent.identity.accession.empty() || ! identities.insert(parent.identity).second) invalid("Empty or duplicate parent identity");
    }
  }
  parents_ = std::move(parents);
}
ID::SourceId ID::Run::addSource(const SourceFile& source)
{
  checkMutation_();
  if (sources_.size() >= std::numeric_limits<UInt32>::max()) invalid("Too many sources");
  SourceId id {static_cast<UInt32>(sources_.size()), token()};
  sources_.push_back({id, source, {}});
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
  // Reserve before changing the logical schema, maintaining the strong guarantee.
  ScoreDefinition copy = definition;
  scores_.reserve(scores_.size() + 1);
  score_owners_.reserve(scores_.size() + 1);
  for (auto& source : sources_)
    for (auto& query : source.identifications)
      for (auto& match : query.matches_)
        match.scores_.reserve(scores_.size() + 1);
  UInt64 owner = token();
  UInt64 new_schema = token();
  scores_.push_back(std::move(copy));
  score_owners_.push_back(owner);
  schema_token_ = new_schema;
  for (auto& source : sources_)
    for (auto& query : source.identifications)
      for (auto& match : query.matches_)
      {
        match.scores_.push_back(std::numeric_limits<double>::quiet_NaN());
        match.schema_token_ = schema_token_;
      }
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
  view.id_ = score;
  view.schema_token_ = schema_token_;
  view.definition_ = scores_[score.value];
  return view;
}
void ID::Run::setPrimaryScore(std::optional<ScoreId> score)
{
  checkMutation_();
  if (score)
  {
    checkScore_(*score);
    for (const auto& source : sources_)
      for (const auto& query : source.identifications)
        for (const auto& match : query.matches_)
          if (std::isnan(match.scores_[score->value])) invalid("Primary score is missing on a candidate");
  }
  primary_ = score;
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
  if (data.adduct && data.adduct->getCharge() != data.charge) invalid("Adduct charge disagrees with ion charge");
  if (kind_ == MoleculeKind::COMPOUND && ! data.parent_evidence.empty()) invalid("Compound candidates cannot contain sequence-to-parent evidence");
  for (const auto& evidence : data.parent_evidence)
  {
    if (evidence.parent.accession.empty()) invalid("Parent evidence needs an accession");
    if (evidence.start && evidence.end && *evidence.start > *evidence.end) invalid("Parent evidence start exceeds end");
  }
}
void ID::Run::prepareLookupIndexes()
{
  checkMutation_();
  ensureQueryIndex_();
  ensureMatchIndex_();
}
void ID::Run::ensureQueryIndex_() const
{
  if (query_index_built_) return;
  std::unordered_map<UInt64, std::array<Size, 2>> index;
  index.reserve(query_count_);
  for (Size s = 0; s < sources_.size(); ++s)
    for (Size q = 0; q < sources_[s].identifications.size(); ++q)
      index.emplace(sources_[s].identifications[q].id_.value, std::array<Size, 2> {s, q});
  query_index_.swap(index);
  query_index_built_ = true;
}
void ID::Run::ensureMatchIndex_() const
{
  if (match_index_built_) return;
  std::unordered_map<UInt64, std::array<Size, 3>> index;
  index.reserve(match_count_);
  for (Size s = 0; s < sources_.size(); ++s)
    for (Size q = 0; q < sources_[s].identifications.size(); ++q)
      for (Size m = 0; m < sources_[s].identifications[q].matches_.size(); ++m)
        index.emplace(sources_[s].identifications[q].matches_[m].id_.value, std::array<Size, 3> {s, q, m});
  match_index_.swap(index);
  match_index_built_ = true;
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
  Match match;
  static_cast<MatchData&>(match) = std::move(data);
  match.id_ = id;
  match.schema_token_ = schema_token_;
  match.scores_.assign(scores_.size(), std::numeric_limits<double>::quiet_NaN());
  for (Size i = 0; i < values.size(); ++i)
    if (values[i]) match.scores_[i] = *values[i];
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
  ++match_count_;
  next_match_id_ = std::max(next_match_id_, next);
  last_match_ = position;
  return id;
}
std::optional<double> ID::Run::getScore(MatchId match, ScoreId score) const
{
  checkScore_(score);
  double value = getMatch(match).scores_[score.value];
  return std::isnan(value) ? std::nullopt : std::optional<double>(value);
}
void ID::Run::setScore(MatchId match, ScoreId score, std::optional<double> value)
{
  checkMutation_();
  checkScore_(score);
  validNumber(value);
  if (! value && primary_ == score) invalid("Cannot clear a primary score");
  match_(match).scores_[score.value] = value.value_or(std::numeric_limits<double>::quiet_NaN());
}
void ID::Run::setSelectedMatch(QueryId query_id, std::optional<MatchId> selected)
{
  checkMutation_();
  auto& query = query_(query_id);
  if (selected && std::none_of(query.matches_.begin(), query.matches_.end(), [&](const Match& match) { return match.id_ == selected; }))
    invalid("Selected match is not owned by the query");
  query.selected_ = selected;
}
void ID::Run::replaceObservation(QueryId query, const Observation& observation)
{
  checkMutation_();
  validObservation(observation);
  Observation copy(observation);
  static_cast<Observation&>(query_(query)) = std::move(copy);
}
void ID::Run::replaceMatch(MatchId match, const MatchData& data)
{
  checkMutation_();
  auto& target = match_(match);
  if (! sameHypothesis(target, data)) invalid("Changing a molecular/ion hypothesis requires explicitly replacing its scores");
  replaceMatch(match, data, target.getScores());
}
void ID::Run::replaceMatch(MatchId match, const MatchData& data, const std::vector<std::optional<double>>& values)
{
  checkMutation_();
  auto& target = match_(match);
  validateMatch_(data, values);
  Match replacement(target);
  static_cast<MatchData&>(replacement) = data;
  replacement.scores_.assign(scores_.size(), std::numeric_limits<double>::quiet_NaN());
  for (Size i = 0; i < values.size(); ++i)
    if (values[i]) replacement.scores_[i] = *values[i];
  target = std::move(replacement);
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
  Size at = 0, removed = 0;
  for (auto& source : sources_)
    for (auto& query : source.identifications)
    {
      auto end = std::remove_if(query.matches_.begin(), query.matches_.end(), [&](const Match&) {
        bool erase = ! decisions[at++];
        removed += erase;
        return erase;
      });
      query.matches_.erase(end, query.matches_.end());
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
  import_finalized_ = true;
  invalidateIndexes_();
  return removed;
}
Size ID::Run::eraseMatches(const std::function<bool(const Match&)>& remove, bool keep_empty_queries)
{
  return filterMatches([&](const Match& match) { return ! remove(match); }, keep_empty_queries);
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
      std::optional<double> best;
      for (const auto& match : query.matches_)
      {
        double value = match.scores_[score.value];
        if (std::isnan(value)) continue;
        if (! best || (scores_[score.value].higher_better ? value > *best : value < *best)) best = value;
      }
      bool taken = false;
      for (const auto& match : query.matches_)
      {
        bool keep = best && match.scores_[score.value] == *best && (keep_ties || ! taken);
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
          validateMatch_(data, match.getScores());
          replacements.push_back(std::move(data));
        }
  }
  Size at = 0;
  for (auto& source : sources_)
    for (auto& query : source.identifications)
      for (auto& match : query.matches_)
        static_cast<MatchData&>(match) = std::move(replacements[at++]);
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
}
void ID::Run::reserveMatchId(MatchId id)
{
  checkMutation_();
  next_match_id_ = std::max(next_match_id_, following(id.value));
  import_finalized_ = true;
}
void ID::Run::validate() const
{
  if (! validUuid(uuid_) || ! next_query_id_ || ! next_match_id_) invalid("Invalid run identity");
  Size nq = 0, nm = 0;
  std::vector<UInt64> queries, matches;
  queries.reserve(query_count_);
  matches.reserve(match_count_);
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
        if (match.scores_.size() != scores_.size()) invalid("Inconsistent score column count");
        for (double score : match.scores_)
          if (! std::isnan(score) && ! std::isfinite(score)) invalid("Scores must be finite or missing");
        if (primary_ && (primary_->value >= match.scores_.size() || std::isnan(match.scores_[primary_->value])))
          invalid("Primary score is missing on a candidate");
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
void ID::replaceRun(const Run& run)
{
  checkMutation_();
  Run copy(run);
  copy.validate();
  auto found = std::find_if(runs_.begin(), runs_.end(), [&](const Run& current) { return current.uuid_ == run.uuid_; });
  if (found == runs_.end() || found->identifier_ != run.identifier_) invalid("Replacement must have the same run UUID and identifier");
  checkScoreContract(runs_, &copy, &run.uuid_);
  copy.next_query_id_ = std::max(copy.next_query_id_, found->next_query_id_);
  copy.next_match_id_ = std::max(copy.next_match_id_, found->next_match_id_);
  *found = std::move(copy);
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
    if (!run.primary_ && run.match_count_ == 0 && run.scores_.empty()) continue;
    const auto score = run.findScore(definition);
    for (const auto& source : run.sources_)
      for (const auto& query : source.identifications)
        for (const auto& match : query.getMatches())
          if (std::isnan(match.getScoreValues()[score.value])) invalid("Primary score is missing on a candidate in run '" + run.identifier_ + "'");
    selections.emplace_back(&run, score);
  }
  // Every potentially throwing operation precedes the commit.
  for (auto& [run, score] : selections) run->primary_ = score;
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
