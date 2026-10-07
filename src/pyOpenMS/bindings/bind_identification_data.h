// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
#pragma once

#include "all_casters.h"

#include <OpenMS/ANALYSIS/ID/IdentificationDataInference.h>
#include <OpenMS/FORMAT/IdentificationDataFile.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/METADATA/ID/IdentificationDataAdapter.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <nanobind/operators.h>
#include <algorithm>
#include <concepts>
#include <functional>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

namespace pyopenms_identification
{
namespace nb = nanobind;
using ID = OpenMS::IdentificationData;
using File = OpenMS::IdentificationDataFile;
using Adapter = OpenMS::IdentificationDataAdapter;
using Inference = OpenMS::IdentificationDataInference;

/// Values compare with their C++ operator==. Match and Identification inherit the comparison of their
/// payload base, which ignores ID and scores; they get complete comparisons in bind().
template<typename T>
constexpr bool comparesAsValue = std::equality_comparable<T> && ! std::is_same_v<T, ID::Match> && ! std::is_same_v<T, ID::Identification>;

/// Value equality without a hash: equal values must have equal hashes, and these values are mutable.
/// Identity types (IDs, references) define a value __hash__ after this.
template<typename Class>
void unhashable(Class& cls)
{ cls.attr("__hash__") = nb::none(); }

// All exposed records are owned Python values. Never retain a pointer into a
// run's vector: appending/filtering its owner may invalidate such a pointer.
template<typename T, typename... Bases>
auto valueClass(nb::handle scope, const char* name)
{
  auto cls = nb::class_<T, Bases...>(scope, name, "Owning value. Nested records and containers are returned as copies; assign them back after editing.")
               .def(nb::init<>())
               .def(nb::init<const T&>())
               .def("__copy__", [](const T& self) { return T(self); })
               .def("__deepcopy__", [](const T& self, nb::dict) { return T(self); }, nb::arg("memo"));
  if constexpr (comparesAsValue<T>)
  {
    cls.def(nb::self == nb::self).def(nb::self != nb::self);
    unhashable(cls);
  }
  return cls;
}

/// Metadata as a Python dict {name: value}, as MetaInfoInterface.getMetaValues() returns it.
inline nb::dict metadataDict(const OpenMS::MetaInfoInterface& metadata)
{
  nb::dict result;
  std::vector<std::string> keys;
  metadata.getKeys(keys);
  for (const auto& key : keys)
    result[nb::str(key.c_str())] = nb::cast(metadata.getMetaValue(key));
  return result;
}
/// Fills metadata from a dict {name: value}; values convert like MetaInfoInterface.setMetaValue.
inline void setMetadata(OpenMS::MetaInfoInterface& metadata, nb::handle values, const std::string& keyword)
{
  if (! nb::isinstance<nb::dict>(values)) throw nb::type_error(("'" + keyword + "' expects a dict {name: value}").c_str());
  for (auto [key, item] : nb::borrow<nb::dict>(values))
  {
    const auto name = nb::cast<std::string>(key);
    try
    {
      metadata.setMetaValue(name, nb::cast<OpenMS::DataValue>(item));
    }
    catch (const nb::cast_error&)
    {
      throw nb::type_error(("invalid value for metadata '" + name + "' in '" + keyword + "'").c_str());
    }
  }
}

/// Python truthiness of a callback result (accepts numpy.bool_, 0/1, None, ...).
inline bool truthy(const nb::object& value)
{
  const int result = PyObject_IsTrue(value.ptr());
  if (result < 0) throw nb::python_error();
  return result != 0;
}

/// Hash for a std::string plus integer identity, consistent with the bound __eq__.
inline size_t hashIdentity(const std::string& text, OpenMS::UInt64 value)
{
  return std::hash<std::string> {}(text) ^ (std::hash<OpenMS::UInt64> {}(value) + 0x9e3779b97f4a7c15ULL + (std::hash<std::string> {}(text) << 6));
}

template<typename Class, typename Same>
void compareWith(Class& cls, Same same)
{
  using Self = typename Class::Type;
  cls.def("__eq__", [same](const Self& a, const Self& b) { return same(a, b); }, nb::is_operator())
    .def("__ne__", [same](const Self& a, const Self& b) { return ! same(a, b); }, nb::is_operator());
  unhashable(cls);
}

/// The fields bound with field(), per bound class, in binding order. They drive the keyword
/// constructor, __repr__ and (for types without a C++ operator==) __eq__ of the class.
template<typename Self>
struct BoundField
{
  std::string name;
  std::function<void(Self&, nb::handle)> set;
  std::function<nb::object(const Self&)> get; ///< owned copy; metadata as a dict
};
template<typename Self>
std::vector<BoundField<Self>>& boundFields()
{
  static std::vector<BoundField<Self>> fields;
  return fields;
}
/// Added once all fields of a class are bound (finishFieldProtocols).
inline std::vector<std::function<void()>>& pendingFieldProtocols()
{
  static std::vector<std::function<void()>> pending;
  return pending;
}

/// Keyword constructor, __repr__ and, unless the class already compares, field-wise __eq__ (like a
/// dataclass) from the bound fields; classes with metadata also take and show `metadata`
/// ({name: value}). The constructor builds the value completely before placing it, so a rejected
/// argument leaves nothing half-constructed. For example:
///   QualifiedAccession(database="uniprot.fasta", accession="P02769")
template<typename Class>
void addFieldProtocols(Class& cls)
{
  using Self = typename Class::Type;
  constexpr bool has_metadata = std::is_base_of_v<OpenMS::MetaInfoInterface, Self>;
  // The signature lists the keywords for help() and the generated stubs; nanobind keeps the pointer.
  static std::vector<std::unique_ptr<std::string>> signatures;
  std::string signature = "def __init__(self, *";
  for (const auto& field : boundFields<Self>())
    signature += ", " + field.name + "=...";
  if (has_metadata) signature += ", metadata=...";
  signatures.push_back(std::make_unique<std::string>(signature + ") -> None"));
  cls.def(
    "__init__",
    [](Self* self, nb::kwargs kwargs) {
      Self value;
      const auto& fields = boundFields<Self>();
      for (auto [key, item] : kwargs)
      {
        const auto name = nb::cast<std::string>(key);
        if constexpr (has_metadata)
        {
          if (name == "metadata")
          {
            setMetadata(value, item, name);
            continue;
          }
        }
        const auto field = std::find_if(fields.begin(), fields.end(), [&](const auto& candidate) { return candidate.name == name; });
        if (field == fields.end()) throw nb::type_error(("unexpected keyword argument '" + name + "'; the keywords are the field names").c_str());
        try
        {
          field->set(value, item);
        }
        catch (const nb::cast_error&)
        {
          throw nb::type_error(("invalid type for field '" + name + "'").c_str());
        }
      }
      new (self) Self(std::move(value));
    },
    nb::sig(signatures.back()->c_str()), "Construct from keyword arguments named after the fields; omitted fields keep their defaults.");
  const auto type_name = nb::cast<std::string>(cls.attr("__qualname__"));
  cls.def("__repr__", [type_name](const Self& self) {
    std::string text = type_name + "(";
    const char* separator = "";
    for (const auto& field : boundFields<Self>())
    {
      text += separator + field.name + "=" + nb::cast<std::string>(nb::repr(field.get(self)));
      separator = ", ";
    }
    if constexpr (has_metadata)
    {
      if (! self.isMetaEmpty()) text += separator + std::string("metadata=") + nb::cast<std::string>(nb::repr(metadataDict(self)));
    }
    return text + ")";
  });
  if (! nb::cast<bool>(cls.attr("__dict__").attr("__contains__")("__eq__")))
  {
    compareWith(cls, [](const Self& a, const Self& b) {
      for (const auto& field : boundFields<Self>())
        if (! field.get(a).equal(field.get(b))) return false;
      return true;
    });
  }
}

inline void finishFieldProtocols()
{
  for (const auto& add : pendingFieldProtocols())
    add();
  pendingFieldProtocols().clear();
}

template<typename Class, typename T, typename Value>
void field(Class& cls, const char* name, Value T::* member)
{
  using Self = typename Class::Type;
  auto& fields = boundFields<Self>();
  if (fields.empty()) pendingFieldProtocols().push_back([cls]() mutable { addFieldProtocols(cls); });
  BoundField<Self> bound;
  bound.name = name;
  bound.set = [member, keyword = std::string(name)](Self& self, nb::handle item) {
    if constexpr (std::is_same_v<Value, OpenMS::MetaInfoInterface>)
    {
      // Metadata-valued fields (e.g. ScoreDefinition.parameters) also take a dict {name: value}.
      if (nb::isinstance<nb::dict>(item))
      {
        OpenMS::MetaInfoInterface metadata;
        setMetadata(metadata, item, keyword);
        static_cast<T&>(self).*member = metadata;
        return;
      }
    }
    static_cast<T&>(self).*member = nb::cast<Value>(item);
  };
  bound.get = [member](const Self& self) -> nb::object {
    if constexpr (std::is_same_v<Value, OpenMS::MetaInfoInterface>) return metadataDict(static_cast<const T&>(self).*member);
    else
      return nb::cast(static_cast<const T&>(self).*member, nb::rv_policy::copy);
  };
  fields.push_back(std::move(bound));
  // Explicit value return also protects nested structs and vector elements
  // when Python subsequently replaces the containing field.
  cls.def_prop_rw(
    name, [member](const T& self) -> Value { return self.*member; }, [member](T& self, const Value& value) { self.*member = value; },
    "Owned value; assign an edited nested record or container back to this property.");
}

/// Complete value comparisons for records that inherit the payload comparison of their base.
inline bool sameMatch(const ID::Match& a, const ID::Match& b)
{ return a.getId() == b.getId() && a.getData() == b.getData() && a.getScores() == b.getScores(); }
inline bool sameIdentification(const ID::Identification& a, const ID::Identification& b)
{
  return a.getId() == b.getId() && a.getObservation() == b.getObservation() && a.getSelectedMatch() == b.getSelectedMatch()
         && std::equal(a.getMatches().begin(), a.getMatches().end(), b.getMatches().begin(), b.getMatches().end(), sameMatch);
}
inline bool sameSource(const ID::Source& a, const ID::Source& b)
{
  return a.id.value == b.id.value && a.file == b.file
         && std::equal(a.identifications.begin(), a.identifications.end(), b.identifications.begin(), b.identifications.end(), sameIdentification);
}

inline std::string moleculeKindName(ID::MoleculeKind kind)
{
  switch (kind)
  {
    case ID::MoleculeKind::PEPTIDE:
      return "PEPTIDE";
    case ID::MoleculeKind::OLIGONUCLEOTIDE:
      return "OLIGONUCLEOTIDE";
    default:
      return "COMPOUND";
  }
}
inline std::string runSummary(const std::string& type, const ID::Run& run)
{
  return type + "('" + run.getIdentifier() + "', uuid='" + run.getUuid() + "', kind=MoleculeKind." + moleculeKindName(run.getMoleculeKind())
         + ", queries=" + std::to_string(run.getNumberOfIdentifications()) + ", matches=" + std::to_string(run.getNumberOfMatches()) + ")";
}

/// A run inside a dataset, addressed by UUID: it cannot dangle when the dataset's runs change,
/// and every edit goes to the live run, so IDs are only ever allocated in one place.
struct RunView
{
  nb::object owner; // the IdentificationData (an owned object or a view into a feature/consensus map)
  std::string uuid;
  ID::Run& run() const
  {
    auto* found = nb::cast<ID&>(owner).findRunByUuid(uuid);
    if (! found) throw nb::key_error(("The run " + uuid + " is no longer part of the dataset").c_str());
    return *found;
  }
};

/// Run methods shared by owned runs (Run) and live runs (RunView); names follow the C++ API.
template<typename Class, typename Get>
void bindRunApi(Class& cls, Get get)
{
  using Self = typename Class::Type;
  cls.def("getIdentifier", [get](Self& self) { return get(self).getIdentifier(); })
    .def("getUuid", [get](Self& self) { return get(self).getUuid(); })
    .def("getMoleculeKind", [get](Self& self) { return get(self).getMoleculeKind(); })
    .def("getSettings", [get](Self& self) { return get(self).getSettings(); })
    .def("setSettings", [get](Self& self, const ID::RunSettings& settings) { get(self).setSettings(settings); }, nb::arg("settings"))
    .def("getParents", [get](Self& self) { return get(self).getParents(); })
    .def("setParents", [get](Self& self, std::optional<std::vector<ID::ParentRecord>> parents) { get(self).setParents(std::move(parents)); },
         nb::arg("parents"))
    .def("getSources", [get](Self& self) { return get(self).getSources(); })
    .def("getScoreDefinitions", [get](Self& self) { return get(self).getScoreDefinitions(); })
    .def("addSource", [get](Self& self, const ID::SourceFile& source) { return get(self).addSource(source); }, nb::arg("source"))
    .def("getSourceId", [get](Self& self, OpenMS::UInt32 index) { return get(self).getSourceId(index); }, nb::arg("index"))
    .def("addScore", [get](Self& self, const ID::ScoreDefinition& definition) { return get(self).addScore(definition); }, nb::arg("definition"))
    .def("getScoreId", [get](Self& self, OpenMS::UInt32 index) { return get(self).getScoreId(index); }, nb::arg("index"))
    .def("findScore", [get](Self& self, const ID::ScoreDefinition& definition) { return get(self).findScore(definition); }, nb::arg("definition"))
    .def("getScoreDefinition", [get](Self& self, ID::ScoreId score) { return ID::ScoreDefinition(get(self).getScoreDefinition(score)); },
         nb::arg("score"))
    .def("bindScore", [get](Self& self, ID::ScoreId score) { return get(self).bindScore(score); }, nb::arg("score"))
    .def("getPrimaryScore", [get](Self& self) { return get(self).getPrimaryScore(); })
    .def("setPrimaryScore", [get](Self& self, std::optional<ID::ScoreId> score) { get(self).setPrimaryScore(score); }, nb::arg("score"))
    .def("addIdentification", [get](Self& self, ID::SourceId source, const ID::Observation& observation) {
           return get(self).addIdentification(source, observation);
         }, nb::arg("source"), nb::arg("observation"))
    .def("addMatch", [get](Self& self, ID::QueryId query, const ID::MatchData& value, const std::vector<std::optional<double>>& scores) {
           return get(self).addMatch(query, value, scores);
         }, nb::arg("query"), nb::arg("data"), nb::arg("scores") = std::vector<std::optional<double>> {})
    .def("findIdentification", [get](Self& self, ID::QueryId id) -> std::optional<ID::Identification> {
           const auto* found = get(self).findIdentification(id);
           return found ? std::optional<ID::Identification>(*found) : std::nullopt;
         }, nb::arg("id"))
    .def("findMatch", [get](Self& self, ID::MatchId id) -> std::optional<ID::Match> {
           const auto* found = get(self).findMatch(id);
           return found ? std::optional<ID::Match>(*found) : std::nullopt;
         }, nb::arg("id"))
    .def("getIdentification", [get](Self& self, ID::QueryId id) { return ID::Identification(get(self).getIdentification(id)); }, nb::arg("id"))
    .def("getMatch", [get](Self& self, ID::MatchId id) { return ID::Match(get(self).getMatch(id)); }, nb::arg("id"))
    .def("getScore", [get](Self& self, ID::MatchId match, ID::ScoreId score) { return get(self).getScore(match, score); }, nb::arg("match"),
         nb::arg("score"))
    .def("setScore", [get](Self& self, ID::MatchId match, ID::ScoreId score, std::optional<double> value) { get(self).setScore(match, score, value); },
         nb::arg("match"), nb::arg("score"), nb::arg("value"))
    .def("setSelectedMatch", [get](Self& self, ID::QueryId query, std::optional<ID::MatchId> selected) { get(self).setSelectedMatch(query, selected); },
         nb::arg("query"), nb::arg("selected"))
    .def("replaceObservation", [get](Self& self, ID::QueryId query, const ID::Observation& observation) {
           get(self).replaceObservation(query, observation);
         }, nb::arg("query"), nb::arg("observation"))
    .def("replaceMatch", [get](Self& self, ID::MatchId id, const ID::MatchData& value) { get(self).replaceMatch(id, value); }, nb::arg("match"),
         nb::arg("data"))
    .def("replaceMatch", [get](Self& self, ID::MatchId id, const ID::MatchData& value, const std::vector<std::optional<double>>& scores) {
           get(self).replaceMatch(id, value, scores);
         }, nb::arg("match"), nb::arg("data"), nb::arg("scores"))
    .def("filterMatches", [get](Self& self, nb::callable keep, bool keep_empty) {
           return get(self).filterMatches([&](const ID::Match& match) { return truthy(keep(nb::cast(match, nb::rv_policy::copy))); }, keep_empty);
         }, nb::arg("keep"), nb::arg("keep_empty_queries") = false)
    .def("eraseMatches", [get](Self& self, nb::callable remove, bool keep_empty) {
           return get(self).eraseMatches([&](const ID::Match& match) { return truthy(remove(nb::cast(match, nb::rv_policy::copy))); }, keep_empty);
         }, nb::arg("remove"), nb::arg("keep_empty_queries") = false)
    .def("retainBest", [get](Self& self, ID::ScoreId score, bool keep_ties, bool keep_empty) { return get(self).retainBest(score, keep_ties, keep_empty); },
         nb::arg("score"), nb::arg("keep_ties") = true, nb::arg("keep_empty_queries") = false)
    .def("transformMatches", [get](Self& self, nb::callable transform) {
           get(self).transformMatches([&](ID::MatchData& payload) {
             auto owned = nb::cast(payload, nb::rv_policy::copy);
             transform(owned);
             payload = nb::cast<ID::MatchData>(owned);
           });
         }, nb::arg("transform"))
    .def("getNumberOfIdentifications", [get](Self& self) { return get(self).getNumberOfIdentifications(); })
    .def("getNumberOfMatches", [get](Self& self) { return get(self).getNumberOfMatches(); })
    .def("prepareLookupIndexes", [get](Self& self) { get(self).prepareLookupIndexes(); })
    .def("getIdentificationForMatch", [get](Self& self, ID::MatchId match) { return ID::Identification(get(self).getIdentificationForMatch(match)); },
         nb::arg("match"))
    .def("getNextQueryId", [get](Self& self) { return get(self).getNextQueryId(); })
    .def("getNextMatchId", [get](Self& self) { return get(self).getNextMatchId(); })
    .def("validate", [get](Self& self) { get(self).validate(); });
}

inline void bind(nb::module_& m)
{
  auto data = valueClass<ID>(m, "IdentificationData");
  nb::enum_<ID::MoleculeKind>(data, "MoleculeKind")
    .value("PEPTIDE", ID::MoleculeKind::PEPTIDE)
    .value("OLIGONUCLEOTIDE", ID::MoleculeKind::OLIGONUCLEOTIDE)
    .value("COMPOUND", ID::MoleculeKind::COMPOUND);
  nb::enum_<ID::Encoding>(data, "Encoding")
    .value("AA_SEQUENCE", ID::Encoding::AA_SEQUENCE)
    .value("NA_SEQUENCE", ID::Encoding::NA_SEQUENCE)
    .value("SMILES", ID::Encoding::SMILES)
    .value("INCHI", ID::Encoding::INCHI)
    .value("DATABASE_ID", ID::Encoding::DATABASE_ID);
  nb::enum_<ID::TargetDecoy>(data, "TargetDecoy")
    .value("UNKNOWN", ID::TargetDecoy::UNKNOWN)
    .value("TARGET", ID::TargetDecoy::TARGET)
    .value("DECOY", ID::TargetDecoy::DECOY)
    .value("BOTH", ID::TargetDecoy::BOTH);
  nb::enum_<ID::ScoreScope>(data, "ScoreScope")
    .value("MATCH", ID::ScoreScope::MATCH)
    .value("PEPTIDE", ID::ScoreScope::PEPTIDE)
    .value("PROTEIN", ID::ScoreScope::PROTEIN)
    .value("PROTEIN_GROUP", ID::ScoreScope::PROTEIN_GROUP)
    .value("OTHER", ID::ScoreScope::OTHER);
  nb::enum_<ID::InferencePolicy>(data, "InferencePolicy")
    .value("PRESERVE", ID::InferencePolicy::PRESERVE)
    .value("DISCARD", ID::InferencePolicy::DISCARD);
  auto queryid = valueClass<ID::QueryId>(data, "QueryId");
  queryid.def_ro("value", &ID::QueryId::value);
  queryid.def(nb::init<OpenMS::UInt64>(), nb::arg("value"));
  queryid.def("__hash__", [](const ID::QueryId& self) { return std::hash<OpenMS::UInt64> {}(self.value); })
    .def("__repr__", [](const ID::QueryId& self) { return "QueryId(" + std::to_string(self.value) + ")"; });
  auto matchid = valueClass<ID::MatchId>(data, "MatchId");
  matchid.def_ro("value", &ID::MatchId::value);
  matchid.def(nb::init<OpenMS::UInt64>(), nb::arg("value"));
  matchid.def("__hash__", [](const ID::MatchId& self) { return std::hash<OpenMS::UInt64> {}(self.value); })
    .def("__repr__", [](const ID::MatchId& self) { return "MatchId(" + std::to_string(self.value) + ")"; });
  auto queryreference = valueClass<ID::QueryReference>(data, "QueryReference");
  field(queryreference, "run_uuid", &ID::QueryReference::run_uuid);
  field(queryreference, "query", &ID::QueryReference::query);
  queryreference.def("__hash__", [](const ID::QueryReference& self) { return hashIdentity(self.run_uuid, self.query.value); });
  auto matchreference = valueClass<ID::MatchReference>(data, "MatchReference");
  field(matchreference, "run_uuid", &ID::MatchReference::run_uuid);
  field(matchreference, "match", &ID::MatchReference::match);
  matchreference.def("__hash__", [](const ID::MatchReference& self) { return hashIdentity(self.run_uuid, self.match.value); });
  auto moleculeidentity = valueClass<ID::MoleculeIdentity>(data, "MoleculeIdentity");
  field(moleculeidentity, "encoding", &ID::MoleculeIdentity::encoding);
  field(moleculeidentity, "representation", &ID::MoleculeIdentity::representation);
  moleculeidentity.def("__hash__", [](const ID::MoleculeIdentity& self) { return hashIdentity(self.representation, static_cast<OpenMS::UInt64>(self.encoding)); });
  auto scoreid = valueClass<ID::ScoreId>(data, "ScoreId");
  scoreid.def_ro("value", &ID::ScoreId::value);
  scoreid.def("__hash__", [](const ID::ScoreId& self) { return std::hash<OpenMS::UInt64> {}(self.owner) ^ (std::hash<OpenMS::UInt64> {}(self.value) << 1); })
    .def("__repr__", [](const ID::ScoreId& self) { return "<ScoreId " + std::to_string(self.value) + ">"; });
  auto sourceid = valueClass<ID::SourceId>(data, "SourceId");
  sourceid.def_ro("value", &ID::SourceId::value);
  sourceid.def("__hash__", [](const ID::SourceId& self) { return std::hash<OpenMS::UInt64> {}(self.owner) ^ (std::hash<OpenMS::UInt64> {}(self.value) << 1); })
    .def("__repr__", [](const ID::SourceId& self) { return "<SourceId " + std::to_string(self.value) + ">"; });
  auto qualifiedaccession = valueClass<ID::QualifiedAccession>(data, "QualifiedAccession");
  field(qualifiedaccession, "database", &ID::QualifiedAccession::database);
  field(qualifiedaccession, "accession", &ID::QualifiedAccession::accession);
  qualifiedaccession.def("__hash__", [](const ID::QualifiedAccession& self) { return hashIdentity(self.database, std::hash<std::string> {}(self.accession)); });
  auto scoredefinition = valueClass<ID::ScoreDefinition>(data, "ScoreDefinition");
  field(scoredefinition, "name", &ID::ScoreDefinition::name);
  field(scoredefinition, "accession", &ID::ScoreDefinition::accession);
  field(scoredefinition, "higher_better", &ID::ScoreDefinition::higher_better);
  field(scoredefinition, "scope", &ID::ScoreDefinition::scope);
  field(scoredefinition, "software", &ID::ScoreDefinition::software);
  field(scoredefinition, "software_version", &ID::ScoreDefinition::software_version);
  field(scoredefinition, "parameters", &ID::ScoreDefinition::parameters);
  field(scoredefinition, "calibration", &ID::ScoreDefinition::calibration);
  field(scoredefinition, "aggregation", &ID::ScoreDefinition::aggregation);
  auto runsettings = valueClass<ID::RunSettings, OpenMS::MetaInfoInterface>(data, "RunSettings");
  field(runsettings, "software", &ID::RunSettings::software);
  field(runsettings, "software_version", &ID::RunSettings::software_version);
  field(runsettings, "date", &ID::RunSettings::date);
  field(runsettings, "search", &ID::RunSettings::search);
  auto sourcefile = valueClass<ID::SourceFile, OpenMS::MetaInfoInterface>(data, "SourceFile");
  field(sourcefile, "identifier", &ID::SourceFile::identifier);
  field(sourcefile, "path", &ID::SourceFile::path);
  auto parentevidence = valueClass<ID::ParentEvidence>(data, "ParentEvidence");
  field(parentevidence, "parent", &ID::ParentEvidence::parent);
  field(parentevidence, "start", &ID::ParentEvidence::start);
  field(parentevidence, "end", &ID::ParentEvidence::end);
  field(parentevidence, "before", &ID::ParentEvidence::before);
  field(parentevidence, "after", &ID::ParentEvidence::after);
  auto parentrecord = valueClass<ID::ParentRecord, OpenMS::MetaInfoInterface>(data, "ParentRecord");
  field(parentrecord, "identity", &ID::ParentRecord::identity);
  field(parentrecord, "sequence", &ID::ParentRecord::sequence);
  field(parentrecord, "description", &ID::ParentRecord::description);
  field(parentrecord, "target_decoy", &ID::ParentRecord::target_decoy);
  auto observation = valueClass<ID::Observation, OpenMS::MetaInfoInterface>(data, "Observation");
  field(observation, "data_id", &ID::Observation::data_id);
  field(observation, "rt", &ID::Observation::rt);
  field(observation, "mz", &ID::Observation::mz);
  auto matchdata = valueClass<ID::MatchData, OpenMS::MetaInfoInterface>(data, "MatchData");
  field(matchdata, "representation", &ID::MatchData::representation);
  field(matchdata, "encoding", &ID::MatchData::encoding);
  field(matchdata, "charge", &ID::MatchData::charge);
  field(matchdata, "calculated_mz", &ID::MatchData::calculated_mz);
  field(matchdata, "target_decoy", &ID::MatchData::target_decoy);
  field(matchdata, "name", &ID::MatchData::name);
  field(matchdata, "formula", &ID::MatchData::formula);
  field(matchdata, "identifiers", &ID::MatchData::identifiers);
  field(matchdata, "adduct", &ID::MatchData::adduct);
  field(matchdata, "parent_evidence", &ID::MatchData::parent_evidence);
  field(matchdata, "peak_annotations", &ID::MatchData::peak_annotations);
  auto inferenceinput = valueClass<ID::InferenceInput>(data, "InferenceInput");
  field(inferenceinput, "run_identifier", &ID::InferenceInput::run_identifier);
  field(inferenceinput, "run_uuid", &ID::InferenceInput::run_uuid);
  field(inferenceinput, "score", &ID::InferenceInput::score);
  field(inferenceinput, "selection", &ID::InferenceInput::selection);
  auto inferenceresult = valueClass<ID::InferenceResult>(data, "InferenceResult");
  field(inferenceresult, "identifier", &ID::InferenceResult::identifier);
  field(inferenceresult, "proteins", &ID::InferenceResult::proteins);
  field(inferenceresult, "parent_score", &ID::InferenceResult::parent_score);
  field(inferenceresult, "group_score", &ID::InferenceResult::group_score);
  field(inferenceresult, "qualified_accessions", &ID::InferenceResult::qualified_accessions);
  field(inferenceresult, "inputs", &ID::InferenceResult::inputs);

  // Match and Identification are records of a run: they are read from a run, not constructed with
  // keywords, and compare by ID, payload, scores and candidates.
  auto match = valueClass<ID::Match, ID::MatchData>(data, "Match");
  match.def("getId", &ID::Match::getId)
    .def("getData", [](const ID::Match& self) { return ID::MatchData(self.getData()); })
    .def("getScores", &ID::Match::getScores)
    .def("__repr__", [](const ID::Match& self) {
      return "Match(id=" + std::to_string(self.getId().value) + ", scores=" + nb::cast<std::string>(nb::repr(nb::cast(self.getScores())))
             + ", data=" + nb::cast<std::string>(nb::repr(nb::cast(ID::MatchData(self.getData())))) + ")";
    });
  compareWith(match, sameMatch);
  auto identification = valueClass<ID::Identification, ID::Observation>(data, "Identification");
  identification.def("getId", &ID::Identification::getId)
    .def("getObservation", [](const ID::Identification& self) { return ID::Observation(self.getObservation()); })
    .def("getMatches", [](const ID::Identification& self) { return self.getMatches(); })
    .def("getSelectedMatch", &ID::Identification::getSelectedMatch)
    .def("__repr__", [](const ID::Identification& self) {
      const auto selected = self.getSelectedMatch();
      return "Identification(id=" + std::to_string(self.getId().value) + ", matches=" + std::to_string(self.getMatches().size())
             + ", selected=" + (selected ? std::to_string(selected->value) : std::string("None"))
             + ", observation=" + nb::cast<std::string>(nb::repr(nb::cast(ID::Observation(self.getObservation())))) + ")";
    });
  compareWith(identification, sameIdentification);
  auto source = valueClass<ID::Source>(data, "Source");
  field(source, "id", &ID::Source::id);
  field(source, "file", &ID::Source::file);
  field(source, "identifications", &ID::Source::identifications);
  compareWith(source, sameSource);
  valueClass<ID::ScoreView>(data, "ScoreView")
    .def("__call__", &ID::ScoreView::operator(), nb::arg("match"))
    .def("getDefinition", [](const ID::ScoreView& self) { return self.getDefinition(); });

  auto run = valueClass<ID::Run>(data, "Run");
  run.def(nb::init<std::string, ID::MoleculeKind>(), nb::arg("identifier"), nb::arg("kind") = ID::MoleculeKind::PEPTIDE);
  bindRunApi(run, [](ID::Run& self) -> ID::Run& { return self; });
  run.def("__repr__", [](const ID::Run& self) { return runSummary("Run", self); });
  // Construction from persisted values; not offered on a RunView, whose identity is fixed.
  run.def("importIdentification", &ID::Run::importIdentification, nb::arg("source"), nb::arg("id"), nb::arg("observation"))
    .def("importMatch", &ID::Run::importMatch, nb::arg("query"), nb::arg("id"), nb::arg("data"),
         nb::arg("scores") = std::vector<std::optional<double>> {})
    .def("restoreIdentity", &ID::Run::restoreIdentity, nb::arg("uuid"), nb::arg("next_query"), nb::arg("next_match"))
    .def("reserveMatchId", &ID::Run::reserveMatchId, nb::arg("id"));

  auto run_view = nb::class_<RunView>(data, "RunView",
                                      "Live access to a run inside an IdentificationData, looked up by UUID on every call. "
                                      "Edits land in the dataset directly; the view keeps the dataset alive and raises if the run is gone.");
  bindRunApi(run_view, [](RunView& self) -> ID::Run& { return self.run(); });
  run_view.def("copy", [](RunView& self) { return ID::Run(self.run()); }, "Return an owned copy of the run (a snapshot; edits to it do not reach the dataset).")
    .def("__repr__", [](const RunView& self) {
      const auto* found = nb::cast<ID&>(self.owner).findRunByUuid(self.uuid);
      return found ? runSummary("RunView", *found) : "RunView(uuid='" + self.uuid + "', removed)";
    })
    // Views are handles: equal when they address the same run of the same dataset object.
    .def("__eq__", [](const RunView& a, const RunView& b) { return a.owner.is(b.owner) && a.uuid == b.uuid; }, nb::is_operator())
    .def("__ne__", [](const RunView& a, const RunView& b) { return ! a.owner.is(b.owner) || a.uuid != b.uuid; }, nb::is_operator())
    .def("__hash__", [](const RunView& self) { return std::hash<std::string> {}(self.uuid); });

  data
    .def(
      "addRun",
      [](nb::object self, const std::string& identifier, ID::MoleculeKind kind) {
        return RunView {self, nb::cast<ID&>(self).addRun(identifier, kind).getUuid()};
      },
      nb::arg("identifier"), nb::arg("kind") = ID::MoleculeKind::PEPTIDE, "Add a run and return a RunView: edits through it land in this dataset.")
    .def(
      "addRun", [](nb::object self, ID::Run value) { return RunView {self, nb::cast<ID&>(self).addRun(std::move(value)).getUuid()}; },
      nb::arg("run"), "Add a copy of a standalone run and return a RunView of the added run.")
    .def(
      "run_view", [](nb::object self, const std::string& identifier) { return RunView {self, nb::cast<ID&>(self).getRun(identifier).getUuid()}; },
      nb::arg("identifier"), "Live view of a run (see RunView). Use this to edit a run of the dataset.")
    .def(
      "run_view_by_uuid",
      [](nb::object self, const std::string& uuid) {
        if (! nb::cast<ID&>(self).findRunByUuid(uuid)) throw nb::key_error(("No run with UUID " + uuid).c_str());
        return RunView {self, uuid};
      },
      nb::arg("uuid"))
    .def(
      "getRun", [](const ID& self, const std::string& identifier) { return ID::Run(self.getRun(identifier)); }, nb::arg("identifier"),
      "Return an owned copy of a run (a snapshot). Edit runs through run_view().")
    .def(
      "findRunByUuid",
      [](const ID& self, const std::string& uuid) -> std::optional<ID::Run> {
        const auto* found = self.findRunByUuid(uuid);
        return found ? std::optional<ID::Run>(*found) : std::nullopt;
      },
      nb::arg("uuid"))
    .def("getRuns", [](const ID& self) { return std::vector<ID::Run>(self.getRuns().begin(), self.getRuns().end()); })
    .def("getInferenceResults", [](const ID& self) { return self.getInferenceResults(); })
    .def("addInferenceResult", &ID::addInferenceResult, nb::arg("result"))
    .def("clearInferenceResults", &ID::clearInferenceResults)
    .def("empty", &ID::empty)
    .def("clear", &ID::clear)
    .def("merge", &ID::merge, nb::arg("other"))
    .def(
      "filterMatches",
      [](ID& self, nb::callable keep, ID::InferencePolicy policy, bool keep_empty) {
        return self.filterMatches([&](const ID::Match& match) { return truthy(keep(nb::cast(match, nb::rv_policy::copy))); }, policy,
                                  keep_empty);
      },
      nb::arg("keep"), nb::arg("inference_policy"), nb::arg("keep_empty_queries") = false)
    .def("getScoreDefinitions", &ID::getScoreDefinitions, nb::rv_policy::copy)
    .def("getPrimaryScoreDefinition", &ID::getPrimaryScoreDefinition)
    .def("setPrimaryScore", &ID::setPrimaryScore, nb::arg("definition"))
    .def("validate", &ID::validate)
    .def("swap", &ID::swap, nb::arg("other"))
    .def("__repr__", [](const ID& self) {
      nb::list identifiers;
      for (const auto& run : self.getRuns())
        identifiers.append(nb::str(run.getIdentifier().c_str()));
      return "IdentificationData(runs=" + nb::cast<std::string>(nb::repr(identifiers))
             + ", inference_results=" + std::to_string(self.getInferenceResults().size()) + ")";
    });

  // --- IdentificationDataFile ---
  auto file = valueClass<File>(m, "IdentificationDataFile");
  auto options = valueClass<File::Options>(file, "Options");
  field(options, "batch_rows", &File::Options::batch_rows);
  field(options, "row_group_rows", &File::Options::row_group_rows);
  field(options, "batch_bytes", &File::Options::batch_bytes);
  field(options, "row_group_bytes", &File::Options::row_group_bytes);
  field(options, "max_record_bytes", &File::Options::max_record_bytes);
  field(options, "threads", &File::Options::threads);
  field(options, "replace_existing", &File::Options::replace_existing);
  auto projection = valueClass<File::Projection>(file, "Projection");
  field(projection, "molecule", &File::Projection::molecule);
  field(projection, "evidence", &File::Projection::evidence);
  field(projection, "annotations", &File::Projection::annotations);
  field(projection, "metadata", &File::Projection::metadata);
  field(projection, "all_scores", &File::Projection::all_scores);
  field(projection, "score_ids", &File::Projection::score_ids);
  auto scanoptions = valueClass<File::ScanOptions>(file, "ScanOptions");
  field(scanoptions, "buffering", &File::ScanOptions::buffering);
  field(scanoptions, "projection", &File::ScanOptions::projection);
  field(scanoptions, "runs", &File::ScanOptions::runs);
  field(scanoptions, "validate_unique_ids", &File::ScanOptions::validate_unique_ids);
  auto rundescriptor = valueClass<File::RunDescriptor>(file, "RunDescriptor");
  field(rundescriptor, "identifier", &File::RunDescriptor::identifier);
  field(rundescriptor, "uuid", &File::RunDescriptor::uuid);
  field(rundescriptor, "molecule_kind", &File::RunDescriptor::molecule_kind);
  field(rundescriptor, "scores", &File::RunDescriptor::scores);
  field(rundescriptor, "score_columns", &File::RunDescriptor::score_columns);
  field(rundescriptor, "sources", &File::RunDescriptor::sources);
  field(rundescriptor, "primary_score", &File::RunDescriptor::primary_score);
  field(rundescriptor, "query_count", &File::RunDescriptor::query_count);
  field(rundescriptor, "match_count", &File::RunDescriptor::match_count);
  field(rundescriptor, "next_query_id", &File::RunDescriptor::next_query_id);
  field(rundescriptor, "next_match_id", &File::RunDescriptor::next_match_id);
  auto queryrecord = valueClass<File::QueryRecord>(file, "QueryRecord");
  field(queryrecord, "query_id", &File::QueryRecord::query_id);
  field(queryrecord, "source_id", &File::QueryRecord::source_id);
  field(queryrecord, "data", &File::QueryRecord::data);
  field(queryrecord, "selected_match_id", &File::QueryRecord::selected_match_id);
  auto matchrecord = valueClass<File::MatchRecord>(file, "MatchRecord");
  field(matchrecord, "match_id", &File::MatchRecord::match_id);
  field(matchrecord, "query_id", &File::MatchRecord::query_id);
  field(matchrecord, "data", &File::MatchRecord::data);
  field(matchrecord, "scores", &File::MatchRecord::scores);
  auto scanstatistics = valueClass<File::ScanStatistics>(file, "ScanStatistics");
  field(scanstatistics, "queries", &File::ScanStatistics::queries);
  field(scanstatistics, "matches", &File::ScanStatistics::matches);
  field(scanstatistics, "descriptor_bytes", &File::ScanStatistics::descriptor_bytes);

  file
    .def_static(
      "store",
      [](const std::string& path, const ID& values, const File::Options& options) {
        nb::gil_scoped_release release;
        File::store(path, values, options);
      },
      nb::arg("path"),
      nb::arg("data"), nb::arg("options") = File::Options {})
    .def_static(
      "load",
      [](const std::string& path, const File::Options& options) {
        ID values;
        {
          nb::gil_scoped_release release;
          File::load(path, values, options);
        }
        return values;
      },
      nb::arg("path"), nb::arg("options") = File::Options {})
    .def_static(
      "load",
      [](const std::string& path, ID& values, const File::Options& options) {
        nb::gil_scoped_release release;
        File::load(path, values, options);
      },
      nb::arg("path"),
      nb::arg("data"), nb::arg("options") = File::Options {})
    .def_static(
      "loadRun",
      [](const std::string& path, const std::string& selected, const File::Options& options) {
        nb::gil_scoped_release release;
        return File::loadRun(path, selected, options);
      },
      nb::arg("path"), nb::arg("run"), nb::arg("options") = File::Options {})
    .def_static("inspect", &File::inspect, nb::arg("path"), nb::call_guard<nb::gil_scoped_release>())
    .def_static(
      "scan",
      [](const std::string& path, const File::ScanOptions& options, nb::object queries, nb::object matches) {
        File::QueryCallback query_callback;
        File::MatchCallback match_callback;
        if (! queries.is_none())
        {
          if (! PyCallable_Check(queries.ptr())) { throw nb::type_error("queries must be callable or None"); }
          query_callback
            = [&](const std::string& uuid, const std::vector<File::QueryRecord>& batch) { queries(uuid, nb::cast(batch, nb::rv_policy::copy)); };
        }
        if (! matches.is_none())
        {
          if (! PyCallable_Check(matches.ptr())) { throw nb::type_error("matches must be callable or None"); }
          match_callback
            = [&](const std::string& uuid, const std::vector<File::MatchRecord>& batch) { matches(uuid, nb::cast(batch, nb::rv_policy::copy)); };
        }
        return File::scan(path, options, query_callback, match_callback);
      },
      nb::arg("path"), nb::arg("options"), nb::arg("queries") = nb::none(), nb::arg("matches") = nb::none(),
      "Deliver owned batches; callbacks may retain records after the scan returns.")
    .def_static(
      "filter",
      [](const std::string& input, const std::string& output, nb::callable keep, ID::InferencePolicy policy, bool keep_empty,
         const File::Options& options) {
        File::filter(
          input, output,
          [&](const std::string& uuid, const File::MatchRecord& record) { return truthy(keep(uuid, nb::cast(record, nb::rv_policy::copy))); },
          policy, keep_empty, options);
      },
      nb::arg("input"), nb::arg("output"), nb::arg("keep"), nb::arg("inference_policy"), nb::arg("keep_empty_queries") = false,
      nb::arg("options") = File::Options {});

  // --- IdentificationDataAdapter ---
  auto adapter = valueClass<Adapter>(m, "IdentificationDataAdapter");
  nb::enum_<Adapter::LossPolicy>(adapter, "LossPolicy").value("STRICT", Adapter::LossPolicy::STRICT).value("ALLOW", Adapter::LossPolicy::ALLOW);
  nb::enum_<Adapter::MissingLinkPolicy>(adapter, "MissingLinkPolicy")
    .value("REJECT", Adapter::MissingLinkPolicy::REJECT)
    .value("PRUNE", Adapter::MissingLinkPolicy::PRUNE);
  auto adapter_exportoptions = valueClass<Adapter::ExportOptions>(adapter, "ExportOptions");
  field(adapter_exportoptions, "loss_policy", &Adapter::ExportOptions::loss_policy);
  field(adapter_exportoptions, "include_inference", &Adapter::ExportOptions::include_inference);
  field(adapter_exportoptions, "inference_result", &Adapter::ExportOptions::inference_result);
  adapter.attr("QueryReference") = data.attr("QueryReference");
  auto adapter_importresult = valueClass<Adapter::ImportResult>(adapter, "ImportResult");
  field(adapter_importresult, "data", &Adapter::ImportResult::data);
  field(adapter_importresult, "queries", &Adapter::ImportResult::queries);
  auto adapter_legacyresult = valueClass<Adapter::LegacyResult>(adapter, "LegacyResult");
  field(adapter_legacyresult, "proteins", &Adapter::LegacyResult::proteins);
  field(adapter_legacyresult, "peptides", &Adapter::LegacyResult::peptides);
  field(adapter_legacyresult, "queries", &Adapter::LegacyResult::queries);
  field(adapter_legacyresult, "losses", &Adapter::LegacyResult::losses);
  auto adapter_featureassociation = valueClass<Adapter::FeatureAssociation>(adapter, "FeatureAssociation");
  field(adapter_featureassociation, "unassigned", &Adapter::FeatureAssociation::unassigned);
  field(adapter_featureassociation, "feature_path", &Adapter::FeatureAssociation::feature_path);
  field(adapter_featureassociation, "query", &Adapter::FeatureAssociation::query);
  field(adapter_featureassociation, "matches", &Adapter::FeatureAssociation::matches);
  auto adapter_featureimportresult = valueClass<Adapter::FeatureImportResult>(adapter, "FeatureImportResult");
  field(adapter_featureimportresult, "data", &Adapter::FeatureImportResult::data);
  field(adapter_featureimportresult, "associations", &Adapter::FeatureImportResult::associations);

  adapter.def_static("importLegacy", &Adapter::importLegacy, nb::arg("proteins"), nb::arg("peptides"))
    .def_static("fromLegacy", &Adapter::fromLegacy, nb::arg("proteins"), nb::arg("peptides"))
    .def_static(
      "toLegacy", [](const ID& values, const Adapter::ExportOptions& options) { return Adapter::toLegacy(values, options); }, nb::arg("data"),
      nb::arg("options") = Adapter::ExportOptions {})
    .def_static("materializePeptide", &Adapter::materializePeptide, nb::arg("run"), nb::arg("match"), nb::arg("score"))
    .def_static("fromFeatureMap", &Adapter::fromFeatureMap, nb::arg("map"))
    .def_static("fromConsensusMap", &Adapter::fromConsensusMap, nb::arg("map"))
    .def_static(
      "reconcileAssociations",
      [](const ID& values, std::vector<Adapter::FeatureAssociation> associations, Adapter::MissingLinkPolicy policy) {
        const auto removed = Adapter::reconcileAssociations(values, associations, policy);
        return std::make_pair(std::move(associations), removed);
      },
      nb::arg("data"), nb::arg("associations"), nb::arg("policy"), "Return (updated associations, removed link count).")
    .def_static("applyToFeatureMap", &Adapter::applyToFeatureMap, nb::arg("data"), nb::arg("associations"), nb::arg("map"), nb::arg("options"),
                nb::arg("policy"))
    .def_static("applyToConsensusMap", &Adapter::applyToConsensusMap, nb::arg("data"), nb::arg("associations"), nb::arg("map"), nb::arg("options"),
                nb::arg("policy"));

  // --- IdentificationDataInference ---
  auto inference = valueClass<Inference>(m, "IdentificationDataInference");
  nb::enum_<Inference::ProbabilityType>(inference, "ProbabilityType")
    .value("POSTERIOR_ERROR_PROBABILITY", Inference::ProbabilityType::POSTERIOR_ERROR_PROBABILITY)
    .value("POSTERIOR_PROBABILITY", Inference::ProbabilityType::POSTERIOR_PROBABILITY);
  auto inference_input = valueClass<Inference::Input>(inference, "Input");
  field(inference_input, "run_uuid", &Inference::Input::run_uuid);
  field(inference_input, "score", &Inference::Input::score);
  field(inference_input, "probability", &Inference::Input::probability);
  inference
    .def_static(
      "infer",
      [](const ID& values, const std::vector<Inference::Input>& inputs, const std::string& identifier) {
        return Inference::infer(values, inputs, identifier);
      },
      nb::arg("data"), nb::arg("inputs"), nb::arg("identifier"))
    .def_static(
      "infer",
      [](const ID& values, const std::vector<Inference::Input>& inputs, const std::string& identifier, const OpenMS::Param& parameters) {
        return Inference::infer(values, inputs, identifier, parameters);
      },
      nb::arg("data"), nb::arg("inputs"), nb::arg("identifier"), nb::arg("parameters"))
    .def_static(
      "retainProteins",
      [](ID::InferenceResult& result, const std::vector<ID::QualifiedAccession>& retained) {
        Inference::retainProteins(result, {retained.begin(), retained.end()});
      },
      nb::arg("result"), nb::arg("retained"));
  // --- IdentificationDataConverter ---
  nb::class_<OpenMS::IdentificationDataConverter>(m, "IdentificationDataConverter",
                                                  "Conversions of IdentificationData that keep the RNA and compound conventions of mzTab")
    .def_static("exportMzTab", &OpenMS::IdentificationDataConverter::exportMzTab, nb::arg("data"),
                "Return an mzTab document of the dataset (store it with MzTabFile().store)");
  finishFieldProtocols();
}
} // namespace pyopenms_identification
