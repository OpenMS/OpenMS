// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/QPXDatasetExport.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

#include <optional>
#include <string>
#include <variant>
#include <vector>

using namespace OpenMS;
using namespace std;

namespace
{
  using Column = QPXDatasetExport::Column;

  /// The three elements of the psmColumns() sections: a PSM with scores, meta values and two protein mappings; an
  /// element without hits (no row); a PSM without an isotope error
  PeptideIdentificationList psmIds()
  {
    PeptideIdentificationList ids;

    PeptideIdentification first;
    first.setScoreType("q-value");
    first.setHigherScoreBetter(false);
    first.setSpectrumReference("scan=5");
    PeptideHit hit(0.01, 1, 2, AASequence::fromString("PEPTIDEK"));
    hit.setMetaValue("isotope_error", 1);
    hit.setMetaValue("ln(hyperscore)", 3.5);
    hit.setMetaValue("target_decoy", "target");
    PeptideEvidence p1("P1", 4, 11, PeptideEvidence::N_TERMINAL_AA, 'K');
    PeptideEvidence p2("P2", PeptideEvidence::UNKNOWN_POSITION, PeptideEvidence::UNKNOWN_POSITION,
                       PeptideEvidence::UNKNOWN_AA, PeptideEvidence::C_TERMINAL_AA);
    hit.setPeptideEvidences({p1, p2});
    first.setHits({hit});
    ids.push_back(first);

    PeptideIdentification no_hits;
    no_hits.setScoreType("q-value");
    no_hits.setSpectrumReference("scan=6");
    ids.push_back(no_hits);

    PeptideIdentification third;
    third.setScoreType("q-value");
    third.setHigherScoreBetter(false);
    third.setSpectrumReference("scan=7");
    third.setHits({PeptideHit(0.02, 1, 2, AASequence::fromString("SAMPLER"))});
    ids.push_back(third);
    return ids;
  }

  /// The column @p name of @p columns (fails the test and returns an empty column if there is none)
  Column columnOf(const vector<Column>& columns, const string& name)
  {
    for (const Column& c : columns) { if (c.name == name) { return c; } }
    TEST_TRUE(false) // no column of that name
    return Column();
  }

  /// The entry @p name of the additional_scores row @p scores, if there is one
  optional<Column::Score> scoreOf(const vector<Column::Score>& scores, const string& name)
  {
    for (const Column::Score& s : scores) { if (s.name == name) { return s; } }
    return nullopt;
  }

  bool isNull(const Column::StructValue& v) { return holds_alternative<monostate>(v.value); }
  string text(const Column::StructValue& v) { return holds_alternative<string>(v.value) ? get<string>(v.value) : "(not a string)"; }
  Int32 integer(const Column::StructValue& v) { return holds_alternative<Int32>(v.value) ? get<Int32>(v.value) : -999; }
}

START_TEST(QPXDatasetExport, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static std::string scoreName(const std::string& name)))
{
  TEST_STRING_EQUAL(QPXDatasetExport::scoreName("q-value"), "global_qvalue")
  TEST_STRING_EQUAL(QPXDatasetExport::scoreName("MS:1001492"), "percolator_score")
  TEST_STRING_EQUAL(QPXDatasetExport::scoreName("ln(hyperscore)"), "ln_hyperscore")
  TEST_STRING_EQUAL(QPXDatasetExport::scoreName("Posterior Error Probability"), "")
  TEST_STRING_EQUAL(QPXDatasetExport::scoreName("DeltaMass"), "delta_mass")
  TEST_STRING_EQUAL(QPXDatasetExport::scoreName("COMET:xcorr"), "comet_xcorr")
}
END_SECTION

START_SECTION((static std::vector<Column> psmColumns(const PeptideIdentificationList& peptide_ids, const std::vector<bool>& has_spectrum, bool isotope_observed_minus_theoretical)))
{
  const PeptideIdentificationList ids = psmIds();
  vector<Column> c = QPXDatasetExport::psmColumns(ids, {true, false, false}, false);
  TEST_EQUAL(c.size(), 4)
  ABORT_IF(c.size() != 4)
  TEST_STRING_EQUAL(c[0].name, "additional_scores")
  TEST_STRING_EQUAL(c[1].name, "spectrum_reference")
  TEST_STRING_EQUAL(c[2].name, "isotope_error")
  TEST_STRING_EQUAL(c[3].name, "protein_mappings")
  TEST_EQUAL(c[0].type == Column::Type::SCORE_LIST, true)
  TEST_EQUAL(c[1].type == Column::Type::UTF8, true)
  TEST_EQUAL(c[2].type == Column::Type::INT16, true)
  TEST_EQUAL(c[3].type == Column::Type::STRUCT_LIST, true)
  // the columns build a table of two rows (one per element with hits)
  TEST_EQUAL(ArrowIOHelpers::tableRowCount(ArrowIOHelpers::tableOfColumns(c, 2)), 2)

  // spectrum_reference: the native ID of a PSM whose spectrum is there, null otherwise
  const Column reference = columnOf(c, "spectrum_reference");
  TEST_EQUAL(reference.utf8_values.size(), 2)
  TEST_EQUAL(reference.null_rows.size(), 2)
  ABORT_IF(reference.utf8_values.size() != 2 || reference.null_rows.size() != 2)
  TEST_STRING_EQUAL(reference.utf8_values[0], "scan=5")
  TEST_EQUAL(reference.null_rows[0], false)
  TEST_EQUAL(reference.null_rows[1], true)

  // isotope_error: theoretical minus observed in the meta value, flipped to observed minus theoretical
  const Column isotope = columnOf(c, "isotope_error");
  TEST_EQUAL(isotope.int64_values.size(), 2)
  TEST_EQUAL(isotope.null_rows.size(), 2)
  ABORT_IF(isotope.int64_values.size() != 2 || isotope.null_rows.size() != 2)
  TEST_EQUAL(isotope.int64_values[0], -1)
  TEST_EQUAL(isotope.null_rows[0], false)
  TEST_EQUAL(isotope.null_rows[1], true)

  // additional_scores: the main score, every numeric meta value (isotope error with the same sign; the rank, which
  // PeptideHit keeps as a meta value), no target_decoy
  const Column scores = columnOf(c, "additional_scores");
  TEST_EQUAL(scores.score_lists.size(), 2)
  ABORT_IF(scores.score_lists.size() != 2)
  const vector<Column::Score>& row0 = scores.score_lists[0];
  TEST_EQUAL(row0.size(), 4)
  const optional<Column::Score> qvalue = scoreOf(row0, "global_qvalue");
  TEST_TRUE(qvalue.has_value())
  ABORT_IF(!qvalue)
  TEST_REAL_SIMILAR(qvalue->value, 0.01)
  TEST_EQUAL(qvalue->higher_better.has_value() && *qvalue->higher_better == false, true)
  const optional<Column::Score> hyperscore = scoreOf(row0, "ln_hyperscore");
  TEST_TRUE(hyperscore.has_value())
  ABORT_IF(!hyperscore)
  TEST_REAL_SIMILAR(hyperscore->value, 3.5)
  TEST_EQUAL(hyperscore->higher_better.has_value() && *hyperscore->higher_better == true, true)
  const optional<Column::Score> isotope_score = scoreOf(row0, "isotope_error");
  TEST_TRUE(isotope_score.has_value())
  ABORT_IF(!isotope_score)
  TEST_REAL_SIMILAR(isotope_score->value, -1.0)
  TEST_FALSE(isotope_score->higher_better.has_value())
  const optional<Column::Score> rank = scoreOf(row0, "rank");
  TEST_TRUE(rank.has_value())
  ABORT_IF(!rank)
  TEST_REAL_SIMILAR(rank->value, 1.0)
  TEST_FALSE(rank->higher_better.has_value())
  TEST_FALSE(scoreOf(row0, "target_decoy").has_value())
  TEST_EQUAL(scores.score_lists[1].size(), 2)
  TEST_TRUE(scoreOf(scores.score_lists[1], "global_qvalue").has_value())

  // protein_mappings: 1-based positions, '-' for a terminus, null for unknown
  const Column mappings = columnOf(c, "protein_mappings");
  TEST_EQUAL(mappings.struct_fields.size(), 5)
  ABORT_IF(mappings.struct_fields.size() != 5)
  TEST_STRING_EQUAL(mappings.struct_fields[0].name, "accession")
  TEST_EQUAL(mappings.struct_fields[0].nullable, false)
  TEST_STRING_EQUAL(mappings.struct_fields[1].name, "start")
  TEST_STRING_EQUAL(mappings.struct_fields[2].name, "end")
  TEST_STRING_EQUAL(mappings.struct_fields[3].name, "aa_before")
  TEST_STRING_EQUAL(mappings.struct_fields[4].name, "aa_after")
  TEST_EQUAL(mappings.struct_lists.size(), 2)
  ABORT_IF(mappings.struct_lists.size() != 2)
  TEST_EQUAL(mappings.struct_lists[0].size(), 2)
  ABORT_IF(mappings.struct_lists[0].size() != 2)
  const Column::StructEntry& m1 = mappings.struct_lists[0][0];
  const Column::StructEntry& m2 = mappings.struct_lists[0][1];
  ABORT_IF(m1.size() != 5 || m2.size() != 5)
  TEST_STRING_EQUAL(text(m1[0]), "P1")
  TEST_EQUAL(integer(m1[1]), 5)
  TEST_EQUAL(integer(m1[2]), 12)
  TEST_STRING_EQUAL(text(m1[3]), "-")
  TEST_STRING_EQUAL(text(m1[4]), "K")
  TEST_STRING_EQUAL(text(m2[0]), "P2")
  TEST_TRUE(isNull(m2[1]))
  TEST_TRUE(isNull(m2[2]))
  TEST_TRUE(isNull(m2[3]))
  TEST_STRING_EQUAL(text(m2[4]), "-")
  TEST_EQUAL(mappings.struct_lists[1].size(), 0)

  // empty has_spectrum: every spectrum there; observed minus theoretical: the meta value as it is
  c = QPXDatasetExport::psmColumns(ids, {}, true);
  const Column reference_all = columnOf(c, "spectrum_reference");
  ABORT_IF(reference_all.utf8_values.size() != 2 || reference_all.null_rows.size() != 2)
  TEST_STRING_EQUAL(reference_all.utf8_values[0], "scan=5")
  TEST_STRING_EQUAL(reference_all.utf8_values[1], "scan=7")
  TEST_EQUAL(reference_all.null_rows[0], false)
  TEST_EQUAL(reference_all.null_rows[1], false)
  const Column isotope_all = columnOf(c, "isotope_error");
  ABORT_IF(isotope_all.int64_values.size() != 2)
  TEST_EQUAL(isotope_all.int64_values[0], 1)
  const optional<Column::Score> isotope_score_all = scoreOf(columnOf(c, "additional_scores").score_lists[0], "isotope_error");
  TEST_TRUE(isotope_score_all.has_value())
  ABORT_IF(!isotope_score_all)
  TEST_REAL_SIMILAR(isotope_score_all->value, 1.0)

  // a native ID "" gives a null spectrum_reference
  PeptideIdentificationList no_reference = ids;
  no_reference[0].setSpectrumReference("");
  c = QPXDatasetExport::psmColumns(no_reference, {}, true);
  const Column reference_empty = columnOf(c, "spectrum_reference");
  ABORT_IF(reference_empty.null_rows.size() != 2)
  TEST_EQUAL(reference_empty.null_rows[0], true)
  TEST_EQUAL(reference_empty.null_rows[1], false)

  // no elements: the four columns without rows
  c = QPXDatasetExport::psmColumns(PeptideIdentificationList(), {}, true);
  TEST_EQUAL(c.size(), 4)
  for (const Column& column : c)
  {
    TEST_EQUAL(column.score_lists.size() + column.utf8_values.size() + column.int64_values.size()
               + column.struct_lists.size() + column.null_rows.size(), 0)
  }
  TEST_EQUAL(ArrowIOHelpers::tableOfColumns(c, 0) != nullptr, true)
}
END_SECTION

START_SECTION(psmColumns refuses a has_spectrum of another size)
{
  const PeptideIdentificationList ids = psmIds();
  TEST_EXCEPTION(Exception::IllegalArgument, QPXDatasetExport::psmColumns(ids, {true}, true))
  TEST_EXCEPTION(Exception::IllegalArgument, QPXDatasetExport::psmColumns(ids, {true, true, true, true}, true))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
