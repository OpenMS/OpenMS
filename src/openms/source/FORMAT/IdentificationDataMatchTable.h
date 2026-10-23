// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $
#pragma once

// The rows of matches.parquet, shared by native bundles (IdentificationDataFile) and in-memory Arrow
// tables (IdentificationDataArrow).

#include "IdentificationDataFileSupport.h"
#include <cmath>
#include <optional>
#include <set>

namespace OpenMS::Internal::IdentificationDataIO
{
template<class Builder, class T>
void appendOptional(arrow::ArrayBuilder& builder, const std::optional<T>& v)
{
  if (v) append<Builder>(builder, *v);
  else
    check(builder.AppendNull());
}

inline std::shared_ptr<arrow::DataType> identityType()
{ return arrow::struct_({arrow::field("database", arrow::utf8(), false), arrow::field("accession", arrow::utf8(), false)}); }

inline std::shared_ptr<arrow::Schema> matchSchema(const std::vector<std::string>& score_columns)
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
      "sequence_evidence",
      arrow::list(arrow::field("item",
                               arrow::struct_({arrow::field("database", arrow::uint32(), false), arrow::field("accession", arrow::utf8(), false),
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
    arrow::field("metadata", metadataType(), false)};
  for (const auto& name : score_columns)
    fields.push_back(arrow::field(name, arrow::float64()));
  return arrow::schema(fields);
}

struct MatchView
{
  UInt64 match_id, query_id;
  const ID::MatchData& data;
  const std::vector<std::optional<double>>& scores;
};

/// A flanking residue as stored: one character, or an empty string if it is not known.
inline std::string flank(char residue)
{ return residue ? std::string(1, residue) : std::string(); }

template<class Match>
Size matchBytes(const Match& m)
{
  const auto& d = m.data;
  const auto& details = d.details.value_or_default();
  Size bytes = 128 + d.representation.size() + details.name.size() + (details.formula ? details.formula->size() : 0) + m.scores.size() * 9
               + metadataBytes(d);
  for (const auto& i : details.identifiers)
    bytes += 16 + i.database.size() + i.accession.size();
  if (details.adduct) bytes += 32 + details.adduct->getName().size() + details.adduct->getEmpiricalFormula().toString().size();
  for (const auto& e : d.sequence_evidence)
    bytes += 42 + e.accession.size();
  for (const auto& a : d.peak_annotations)
    bytes += 32 + a.annotation.size();
  return bytes;
}

/// One row of matches.parquet: @p writer has column(i), columns() (without the run_uuid column) and finishRow(bytes).
template<class Writer, class Match>
void writeMatch(Writer& writer, const Match& m, const Dictionary& dictionary)
{
  const auto& d = m.data;
  const auto& details = d.details.value_or_default();
  append<arrow::UInt64Builder>(writer.column(0), m.match_id);
  append<arrow::UInt64Builder>(writer.column(1), m.query_id);
  appendText(writer.column(2), d.representation);
  append<arrow::UInt8Builder>(writer.column(3), static_cast<unsigned>(d.encoding));
  append<arrow::Int32Builder>(writer.column(4), d.charge);
  appendOptional<arrow::DoubleBuilder>(writer.column(5), d.calculated_mz);
  append<arrow::UInt8Builder>(writer.column(6), static_cast<unsigned>(d.target_decoy));
  appendText(writer.column(7), details.name);
  if (details.formula) appendText(writer.column(8), *details.formula);
  else
    check(writer.column(8).AppendNull());
  auto& identifiers = static_cast<arrow::ListBuilder&>(writer.column(9));
  check(identifiers.Append());
  auto& identity = *static_cast<arrow::StructBuilder*>(identifiers.value_builder());
  for (const auto& i : details.identifiers)
  {
    check(identity.Append());
    appendText(*identity.field_builder(0), i.database);
    appendText(*identity.field_builder(1), i.accession);
  }
  auto& adduct = static_cast<arrow::StructBuilder&>(writer.column(10));
  if (! details.adduct) check(adduct.AppendNull());
  else
  {
    check(adduct.Append());
    appendText(*adduct.field_builder(0), details.adduct->getName());
    appendText(*adduct.field_builder(1), details.adduct->getEmpiricalFormula().toString());
    append<arrow::Int32Builder>(*adduct.field_builder(2), details.adduct->getCharge());
    append<arrow::UInt32Builder>(*adduct.field_builder(3), details.adduct->getMolMultiplier());
  }
  auto& evidence = static_cast<arrow::ListBuilder&>(writer.column(11));
  check(evidence.Append());
  auto& entry = *static_cast<arrow::StructBuilder*>(evidence.value_builder());
  for (const auto& e : d.sequence_evidence)
  {
    check(entry.Append());
    append<arrow::UInt32Builder>(*entry.field_builder(0), e.database.value);
    appendText(*entry.field_builder(1), e.accession);
    appendOptional<arrow::UInt64Builder>(*entry.field_builder(2), e.start ? std::optional<UInt64>(*e.start) : std::nullopt);
    appendOptional<arrow::UInt64Builder>(*entry.field_builder(3), e.end ? std::optional<UInt64>(*e.end) : std::nullopt);
    appendText(*entry.field_builder(4), flank(e.before));
    appendText(*entry.field_builder(5), flank(e.after));
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
  appendMetadata(writer.column(13), d, dictionary);
  // The shared table has the dataset's score columns; scoreless catalog runs leave them null.
  const Size score_columns = writer.columns() - 14;
  if (m.scores.size() > score_columns) invalid("Match has more scores than the dataset score schema");
  for (Size i = 0; i < score_columns; ++i)
  {
    if (i < m.scores.size()) appendOptional<arrow::DoubleBuilder>(writer.column(14 + i), m.scores[i]);
    else
      check(writer.column(14 + i).AppendNull());
  }
  writer.finishRow(matchBytes(m));
}

inline std::vector<std::string> projectedMatchColumns(const IdentificationDataFile::Projection& projection, const std::vector<std::string>& score_columns)
{
  std::vector<std::string> fields {"match_id", "query_id"};
  if (projection.molecule)
    for (const auto* field : {"representation", "encoding", "charge", "calculated_mz", "target_decoy", "name", "formula", "identifiers", "adduct"})
      fields.emplace_back(field);
  if (projection.evidence) fields.emplace_back("sequence_evidence");
  if (projection.annotations) fields.emplace_back("peak_annotations");
  if (projection.metadata) fields.emplace_back("metadata");
  std::set<UInt32> scores(projection.score_ids.begin(), projection.score_ids.end());
  for (UInt32 score : scores)
    if (score >= score_columns.size()) invalid("Projected score ID does not exist in the run");
  for (Size i = 0; i < score_columns.size(); ++i)
    if (projection.all_scores || scores.contains(static_cast<UInt32>(i))) fields.push_back(score_columns[i]);
  return fields;
}
} // namespace OpenMS::Internal::IdentificationDataIO
