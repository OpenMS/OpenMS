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
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>
#include <OpenMS/FORMAT/QPXFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/TempFiles.h>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/reader.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <optional>
#include <sstream>
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

  /// A spectrum of the peaks (m/z, intensity) @p peaks, in their order
  MSSpectrum spectrumOf(const vector<pair<double, float>>& peaks)
  {
    MSSpectrum s;
    for (const auto& p : peaks) { s.push_back(Peak1D(p.first, p.second)); }
    return s;
  }

  /// A fragment annotation at @p mz
  PeptideHit::PeakAnnotation ann(double mz)
  {
    PeptideHit::PeakAnnotation a;
    a.mz = mz;
    a.intensity = 1.0;
    a.charge = 1;
    a.annotation = "y1+";
    return a;
  }

  /// An MS2 spectrum of the mzView() sections: precursor 500 m/z 2+, centroid
  MSSpectrum ms2(const string& native_id, double rt, const vector<pair<double, float>>& peaks)
  {
    MSSpectrum s = spectrumOf(peaks);
    s.setNativeID(native_id);
    s.setRT(rt);
    s.setMSLevel(2);
    s.setType(SpectrumSettings::SpectrumType::CENTROID);
    Precursor p;
    p.setMZ(500.0);
    p.setCharge(2);
    s.setPrecursors({p});
    return s;
  }

  /// The spectra of the mzView() sections: "scan=1" at 60 s, "scan=2" at 120 s (its peaks in reverse m/z order if
  /// @p unsorted), "scan=1" at 180 s
  MSExperiment mzSpectra(bool unsorted)
  {
    vector<pair<double, float>> peaks2 = {{200.0, 7.0f}, {210.0, 3.0f}, {350.0, 9.0f}, {360.0, 1.0f}};
    if (unsorted) { reverse(peaks2.begin(), peaks2.end()); }
    MSExperiment e;
    e.addSpectrum(ms2("scan=1", 60.0, {{100.0, 1.0f}, {150.0, 2.0f}}));
    e.addSpectrum(ms2("scan=2", 120.0, peaks2));
    e.addSpectrum(ms2("scan=1", 180.0, {{300.0, 4.0f}}));
    return e;
  }

  /// The PSMs of the mzView() sections: scan=2 at 120 s, scan=1 at 60 s, scan=1 at 180 s, scan=9 at 1 s, one without hits
  PeptideIdentificationList mzIds()
  {
    PeptideIdentificationList ids;
    auto psm = [&ids](const string& native_id, double rt, bool with_hit)
    {
      PeptideIdentification id;
      id.setSpectrumReference(native_id);
      id.setRT(rt);
      id.setScoreType("q-value");
      if (with_hit) { id.setHits({PeptideHit(0.01, 1, 2, AASequence::fromString("PEPTIDEK"))}); }
      ids.push_back(id);
    };
    psm("scan=2", 120.0, true);
    psm("scan=1", 60.0, true);
    psm("scan=1", 180.0, true);
    psm("scan=9", 1.0, true);
    psm("scan=3", 240.0, false);
    return ids;
  }

  /// The float list of column @p name at @p row of @p table
  vector<float> floatList(const shared_ptr<arrow::Table>& table, const string& name, int64_t row)
  {
    vector<float> out;
    const auto list = static_pointer_cast<arrow::ListArray>(ArrowIOHelpers::getColumn(table, name));
    if (!list) { return out; }
    const auto values = static_pointer_cast<arrow::FloatArray>(list->values());
    for (int64_t i = list->value_offset(row); i < list->value_offset(row + 1); ++i) { out.push_back(values->Value(i)); }
    return out;
  }

  bool isNull(const Column::StructValue& v) { return holds_alternative<monostate>(v.value); }
  string text(const Column::StructValue& v) { return holds_alternative<string>(v.value) ? get<string>(v.value) : "(not a string)"; }
  Int32 integer(const Column::StructValue& v) { return holds_alternative<Int32>(v.value) ? get<Int32>(v.value) : -999; }

  /// The names of the columns of @p table, in their order
  vector<string> columnNames(const shared_ptr<arrow::Table>& table)
  {
    vector<string> names;
    for (const auto& field : table->schema()->fields()) { names.push_back(field->name()); }
    return names;
  }

  /// The structs of the list column @p name at @p row of @p table; null if the row is null
  shared_ptr<arrow::StructArray> structsAt(const shared_ptr<arrow::Table>& table, const string& name, int64_t row)
  {
    const auto list = static_pointer_cast<arrow::ListArray>(ArrowIOHelpers::getColumn(table, name));
    if (!list || list->IsNull(row)) { return nullptr; }
    return static_pointer_cast<arrow::StructArray>(list->value_slice(row));
  }

  /// The strings of the list column @p name at @p row of @p table
  vector<string> stringList(const shared_ptr<arrow::Table>& table, const string& name, int64_t row)
  {
    vector<string> out;
    const auto list = static_pointer_cast<arrow::ListArray>(ArrowIOHelpers::getColumn(table, name));
    if (!list) { return out; }
    const auto values = static_pointer_cast<arrow::StringArray>(list->values());
    for (int64_t i = list->value_offset(row); i < list->value_offset(row + 1); ++i) { out.push_back(values->GetString(i)); }
    return out;
  }

  /// The footer (schema key-value metadata) of the Parquet file @p file; empty if it cannot be read
  map<string, string> footerOf(const string& file)
  {
    map<string, string> out;
    auto opened = arrow::io::ReadableFile::Open(file);
    if (!opened.ok()) { return out; }
    std::unique_ptr<parquet::arrow::FileReader> reader;
    auto reader_result = parquet::arrow::OpenFile(opened.ValueOrDie(), arrow::default_memory_pool());
    if (!reader_result.ok()) { return out; }
    reader = std::move(reader_result).ValueOrDie();
    shared_ptr<arrow::Schema> schema;
    if (!reader->GetSchema(&schema).ok() || !schema->metadata()) { return out; }
    for (int64_t i = 0; i < schema->metadata()->size(); ++i) { out[schema->metadata()->key(i)] = schema->metadata()->value(i); }
    return out;
  }

  /// The rows of the Parquet file @p file; -1 if it cannot be read
  int64_t rowsOf(const string& file)
  {
    auto opened = arrow::io::ReadableFile::Open(file);
    if (!opened.ok()) { return -1; }
    auto reader = parquet::arrow::OpenFile(opened.ValueOrDie(), arrow::default_memory_pool());
    if (!reader.ok()) { return -1; }
    const auto table = reader.ValueOrDie()->ReadTable();
    if (!table.ok()) { return -1; }
    return table.ValueOrDie()->num_rows();
  }

  /// The mz view of a run @p run_file_name with one MS2 spectrum and PSM per native ID of @p native_ids
  shared_ptr<arrow::Table> mzTableOf(const vector<string>& native_ids, const string& run_file_name)
  {
    MSExperiment spectra;
    PeptideIdentificationList ids;
    double rt = 0.0;
    for (const string& native_id : native_ids)
    {
      rt += 60.0;
      spectra.addSpectrum(ms2(native_id, rt, {{100.0, 1.0f}}));
      PeptideIdentification id;
      id.setSpectrumReference(native_id);
      id.setRT(rt);
      id.setScoreType("q-value");
      id.setHits({PeptideHit(0.01, 1, 2, AASequence::fromString("PEPTIDEK"))});
      ids.push_back(id);
    }
    QPXDatasetExport::PeakFilter all;
    all.top = 0;
    return QPXDatasetExport::mzView(spectra, ids, run_file_name, all).table;
  }

  /// The run of the runView() and writeDataset() sections
  QPXDatasetExport::Run runA()
  {
    QPXDatasetExport::Run run;
    run.accession = "a";
    run.run_file_name = "a";
    run.file_name = "a.mzML";
    QPXDatasetExport::RunSample sample;
    sample.accession = "a";
    run.samples = {sample};
    run.enzymes = {"Trypsin"};
    run.fixed_modifications = {"Carbamidomethyl (C)"};
    run.variable_modifications = {"Oxidation (M)", "Foo (X)"};
    return run;
  }

  /// The steps of the provenanceView() and writeDataset() sections
  vector<QPXDatasetExport::ProvenanceStep> provenanceSteps()
  {
    QPXDatasetExport::ProvenanceStep search{"database_search", "peptide_search", "ProSE",
                                            {{"precursor_mass_tolerance", "10"}, {"enzyme", "Trypsin"}}, "{\"a\":1}", {"psm"}};
    QPXDatasetExport::ProvenanceStep fdr{"fdr_filtering", "psm_fdr", "ProSE", {}, "", {}};
    return {search, fdr};
  }
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
  TEST_EQUAL(bool(reference.null_rows[0]), false)
  TEST_EQUAL(bool(reference.null_rows[1]), true)

  // isotope_error: theoretical minus observed in the meta value, flipped to observed minus theoretical
  const Column isotope = columnOf(c, "isotope_error");
  TEST_EQUAL(isotope.int64_values.size(), 2)
  TEST_EQUAL(isotope.null_rows.size(), 2)
  ABORT_IF(isotope.int64_values.size() != 2 || isotope.null_rows.size() != 2)
  TEST_EQUAL(isotope.int64_values[0], -1)
  TEST_EQUAL(bool(isotope.null_rows[0]), false)
  TEST_EQUAL(bool(isotope.null_rows[1]), true)

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
  TEST_EQUAL(bool(reference_all.null_rows[0]), false)
  TEST_EQUAL(bool(reference_all.null_rows[1]), false)
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
  TEST_EQUAL(bool(reference_empty.null_rows[0]), true)
  TEST_EQUAL(bool(reference_empty.null_rows[1]), false)

  // an isotope error that is no int (e.g. read back from a pin file as a double): an integral double as its integer
  // (with the sign flip), anything else null, without a conversion error
  PeptideIdentificationList typed;
  for (const DataValue& error : {DataValue(2.0), DataValue(0.5), DataValue("1"), DataValue(1.0e6)})
  {
    PeptideIdentification id = ids[0];
    vector<PeptideHit> hits = id.getHits();
    hits[0].setMetaValue("isotope_error", error);
    id.setHits(hits);
    typed.push_back(id);
  }
  c = QPXDatasetExport::psmColumns(typed, {}, false);
  const Column isotope_typed = columnOf(c, "isotope_error");
  TEST_EQUAL(isotope_typed.int64_values.size(), 4)
  TEST_EQUAL(isotope_typed.null_rows.size(), 4)
  ABORT_IF(isotope_typed.int64_values.size() != 4 || isotope_typed.null_rows.size() != 4)
  TEST_EQUAL(isotope_typed.int64_values[0], -2)
  TEST_EQUAL(bool(isotope_typed.null_rows[0]), false)
  TEST_EQUAL(bool(isotope_typed.null_rows[1]), true) // 0.5
  TEST_EQUAL(bool(isotope_typed.null_rows[2]), true) // a string
  TEST_EQUAL(bool(isotope_typed.null_rows[3]), true) // out of the INT16 range

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

START_SECTION((static std::vector<Size> keptPeaks(const MSSpectrum& spectrum, const std::vector<PeptideHit::PeakAnnotation>& annotations, int precursor_charge, const PeakFilter& filter)))
{
  auto kept = [](const MSSpectrum& s, const vector<PeptideHit::PeakAnnotation>& a, int charge, const QPXDatasetExport::PeakFilter& f)
  {
    return QPXDatasetExport::keptPeaks(s, a, charge, f);
  };
  QPXDatasetExport::PeakFilter f;
  f.window = 100.0;
  f.top = 1;
  const MSSpectrum s = spectrumOf({{100.0, 5.0f}, {150.0, 10.0f}, {199.0, 1.0f}, {250.0, 3.0f}, {400.0, 2.0f}});
  // jumping windows from each first peak: [100, 199], [250], [400]
  TEST_EQUAL(kept(s, {}, 2, f) == (vector<Size>{1, 3, 4}), true)
  // an annotated m/z keeps the peak nearest to it within 0.01 Th
  TEST_EQUAL(kept(s, {ann(199.0)}, 2, f) == (vector<Size>{1, 2, 3, 4}), true)
  TEST_EQUAL(kept(s, {ann(199.02)}, 2, f) == (vector<Size>{1, 3, 4}), true)

  // top 0: all peaks
  QPXDatasetExport::PeakFilter all;
  all.top = 0;
  TEST_EQUAL(kept(s, {}, 2, all) == (vector<Size>{0, 1, 2, 3, 4}), true)
  TEST_EQUAL(kept(MSSpectrum(), {ann(199.0)}, 2, f).empty(), true)

  // the last window keeps its full quota (WindowMower would scale it by its width)
  QPXDatasetExport::PeakFilter two;
  two.window = 100.0;
  two.top = 2;
  const MSSpectrum last = spectrumOf({{100.0, 1.0f}, {110.0, 2.0f}, {400.0, 3.0f}, {450.0, 4.0f}});
  TEST_EQUAL(kept(last, {}, 2, two) == (vector<Size>{0, 1, 2, 3}), true)

  // deisotoped: an annotated 1+ m/z keeps the raw peak of its fragment at charge 2 and 3, up to the precursor's charge
  const double p = Constants::PROTON_MASS_U;
  const MSSpectrum d = spectrumOf({{(600.0 + 2 * p) / 3, 1.0f}, {(600.0 + p) / 2, 1.0f}, {500.0, 100.0f}});
  QPXDatasetExport::PeakFilter deiso;
  deiso.window = 1000.0;
  deiso.top = 1;
  TEST_EQUAL(kept(d, {ann(600.0)}, 2, deiso) == (vector<Size>{2}), true) // not deisotoped: 1+ only
  deiso.deisotoped = true;
  TEST_EQUAL(kept(d, {ann(600.0)}, 2, deiso) == (vector<Size>{1, 2}), true)
  TEST_EQUAL(kept(d, {ann(600.0)}, 0, deiso) == (vector<Size>{0, 1, 2}), true) // unknown charge: up to 3
  TEST_EQUAL(kept(d, {ann(600.0)}, 1, deiso) == (vector<Size>{2}), true)
  deiso.cap_at_precursor = false;
  TEST_EQUAL(kept(d, {ann(600.0)}, 2, deiso) == (vector<Size>{0, 1, 2}), true)
}
END_SECTION

START_SECTION((static MzView mzView(const MSExperiment& spectra, const PeptideIdentificationList& peptide_ids, const std::string& run_file_name, const PeakFilter& filter)))
{
  QPXDatasetExport::PeakFilter all;
  all.top = 0;
  const QPXDatasetExport::MzView view = QPXDatasetExport::mzView(mzSpectra(false), mzIds(), "run1", all);
  // scan=2 and scan=1 at 60 s have rows; scan=1 at 180 s shares its native ID with the row of the spectrum at 60 s
  // (unresolved); scan=9 has no spectrum; the element without hits has no PSM
  TEST_EQUAL(view.has_spectrum == (vector<bool>{true, true, false, false, false}), true)
  TEST_EQUAL(view.unresolved, 1)
  TEST_EQUAL(view.without_spectrum, 1)
  TEST_EQUAL(view.peaks_read, 6)
  TEST_EQUAL(view.peaks_kept, view.peaks_read)
  TEST_EQUAL(view.table != nullptr, true)
  ABORT_IF(view.table == nullptr)
  TEST_EQUAL(ArrowIOHelpers::tableRowCount(view.table), 2)

  const vector<string> names = {"id", "run_file_name", "scan", "ms_level", "centroid", "scan_start_time",
                                "inverse_ion_mobility", "ion_injection_time", "total_ion_current", "precursors", "mz",
                                "intensity", "cv_params"};
  TEST_EQUAL(view.table->num_columns(), static_cast<int>(names.size()))
  ABORT_IF(view.table->num_columns() != static_cast<int>(names.size()))
  for (Size i = 0; i < names.size(); ++i) { TEST_STRING_EQUAL(view.table->schema()->field(static_cast<int>(i))->name(), names[i]) }

  // the rows in the order of the first PSM of each native ID
  const auto id = ArrowIOHelpers::getColumn(view.table, "id");
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(id, 0), "scan=2")
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(id, 1), "scan=1")
  const auto run = ArrowIOHelpers::getColumn(view.table, "run_file_name");
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(run, 0), "run1")
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(run, 1), "run1")
  const auto ms_level = ArrowIOHelpers::getColumn(view.table, "ms_level");
  TEST_EQUAL(ArrowIOHelpers::getInt32Value(ms_level, 0), 2)
  const auto centroid = ArrowIOHelpers::getColumn(view.table, "centroid");
  TEST_EQUAL(ArrowIOHelpers::getBoolValue(centroid, 0), true)
  const auto start_time = ArrowIOHelpers::getColumn(view.table, "scan_start_time");
  TEST_REAL_SIMILAR(ArrowIOHelpers::getFloatValue(start_time, 0), 2.0)
  TEST_REAL_SIMILAR(ArrowIOHelpers::getFloatValue(start_time, 1), 1.0)
  TEST_TRUE(ArrowIOHelpers::isNull(ArrowIOHelpers::getColumn(view.table, "inverse_ion_mobility"), 0))
  TEST_TRUE(ArrowIOHelpers::isNull(ArrowIOHelpers::getColumn(view.table, "ion_injection_time"), 0))
  const auto tic = ArrowIOHelpers::getColumn(view.table, "total_ion_current");
  TEST_REAL_SIMILAR(ArrowIOHelpers::getFloatValue(tic, 0), 20.0)
  TEST_REAL_SIMILAR(ArrowIOHelpers::getFloatValue(tic, 1), 3.0)
  const auto cv_params = ArrowIOHelpers::getColumn(view.table, "cv_params");
  TEST_TRUE(ArrowIOHelpers::isNull(cv_params, 0))
  TEST_TRUE(ArrowIOHelpers::isNull(cv_params, 1))
  // the peaks of the spectrum at 60 s, not those at 180 s
  TEST_EQUAL(floatList(view.table, "mz", 1) == (vector<float>{100.0f, 150.0f}), true)
  TEST_EQUAL(floatList(view.table, "intensity", 1) == (vector<float>{1.0f, 2.0f}), true)

  // the precursor charge caps the peak protection of a deisotoped spectrum: the spectrum's (2+), not the PSM's
  const double p = Constants::PROTON_MASS_U;
  MSExperiment deisotoped;
  deisotoped.addSpectrum(ms2("scan=1", 60.0, {{(600.0 + 2 * p) / 3, 1.0f}, {(600.0 + p) / 2, 1.0f}, {500.0, 100.0f}}));
  PeptideIdentificationList ids;
  PeptideIdentification psm;
  psm.setSpectrumReference("scan=1");
  psm.setRT(60.0);
  PeptideHit hit(0.01, 1, 3, AASequence::fromString("PEPTIDEK"));
  hit.setPeakAnnotations({ann(600.0)});
  psm.setHits({hit});
  ids.push_back(psm);
  QPXDatasetExport::PeakFilter deiso;
  deiso.window = 1000.0;
  deiso.top = 1;
  deiso.deisotoped = true;
  QPXDatasetExport::MzView capped = QPXDatasetExport::mzView(deisotoped, ids, "run1", deiso);
  ABORT_IF(capped.table == nullptr)
  TEST_EQUAL(capped.peaks_read, 3)
  TEST_EQUAL(capped.peaks_kept, 2)
  TEST_EQUAL(floatList(capped.table, "mz", 0) == (vector<float>{static_cast<float>((600.0 + p) / 2), 500.0f}), true)
  // a spectrum of unknown precursor charge: up to 3+
  deisotoped[0].getPrecursors()[0].setCharge(0);
  capped = QPXDatasetExport::mzView(deisotoped, ids, "run1", deiso);
  TEST_EQUAL(capped.peaks_kept, 3)
  // ... also if its PSM is of charge 2: the search deisotoped it up to 3+, so the fragment's 3+ peak stays
  ids[0].getHits()[0].setCharge(2);
  capped = QPXDatasetExport::mzView(deisotoped, ids, "run1", deiso);
  TEST_EQUAL(capped.peaks_kept, 3)
  TEST_EQUAL(floatList(capped.table, "mz", 0).front(), static_cast<float>((600.0 + 2 * p) / 3))

  // no PSMs: an empty view
  const QPXDatasetExport::MzView empty = QPXDatasetExport::mzView(mzSpectra(false), PeptideIdentificationList(), "run1", all);
  TEST_EQUAL(empty.table != nullptr, true)
  TEST_EQUAL(ArrowIOHelpers::tableRowCount(empty.table), 0)
  TEST_EQUAL(empty.has_spectrum.empty(), true)
}
END_SECTION

START_SECTION(mzView sorts unsorted spectra)
{
  QPXDatasetExport::PeakFilter all;
  all.top = 0;
  QPXDatasetExport::PeakFilter top1;
  top1.top = 1;
  for (const QPXDatasetExport::PeakFilter& f : {all, top1})
  {
    const QPXDatasetExport::MzView sorted = QPXDatasetExport::mzView(mzSpectra(false), mzIds(), "run1", f);
    const QPXDatasetExport::MzView unsorted = QPXDatasetExport::mzView(mzSpectra(true), mzIds(), "run1", f);
    ABORT_IF(sorted.table == nullptr || unsorted.table == nullptr)
    TEST_EQUAL(unsorted.has_spectrum == sorted.has_spectrum, true)
    TEST_EQUAL(unsorted.peaks_kept, sorted.peaks_kept)
    for (int64_t row = 0; row < 2; ++row)
    {
      const vector<float> mz = floatList(unsorted.table, "mz", row);
      TEST_EQUAL(mz == floatList(sorted.table, "mz", row), true)
      TEST_EQUAL(floatList(unsorted.table, "intensity", row) == floatList(sorted.table, "intensity", row), true)
      TEST_EQUAL(is_sorted(mz.begin(), mz.end()), true)
    }
  }
  // the row of "scan=2" (ascending m/z): all peaks; the most intense per 100 Th
  QPXDatasetExport::MzView view = QPXDatasetExport::mzView(mzSpectra(true), mzIds(), "run1", all);
  TEST_EQUAL(floatList(view.table, "mz", 0) == (vector<float>{200.0f, 210.0f, 350.0f, 360.0f}), true)
  TEST_EQUAL(floatList(view.table, "intensity", 0) == (vector<float>{7.0f, 3.0f, 9.0f, 1.0f}), true)
  view = QPXDatasetExport::mzView(mzSpectra(true), mzIds(), "run1", top1);
  TEST_EQUAL(floatList(view.table, "mz", 0) == (vector<float>{200.0f, 350.0f}), true)
}
END_SECTION

START_SECTION((static std::shared_ptr<arrow::Table> runView(const std::vector<Run>& runs)))
{
  QPXDatasetExport::Run b = runA();
  b.accession = "b";
  b.run_file_name = "b";
  b.file_name = "b.mzML";
  b.enzymes.clear();
  b.fixed_modifications.clear();
  b.variable_modifications.clear();
  QPXDatasetExport::RunSample replicate;
  replicate.accession = "b";
  replicate.label = "TMT126";
  replicate.biological_replicate = 2;
  b.samples = {replicate};
  b.fraction = "1";
  b.instrument = ""; // "": null
  const shared_ptr<arrow::Table> table = QPXDatasetExport::runView({runA(), b});
  TEST_EQUAL(table != nullptr, true)
  ABORT_IF(table == nullptr)
  TEST_EQUAL(table->num_rows(), 2)
  const vector<string> names = {"run_accession", "run_file_name", "file_name", "samples", "fraction", "instrument", "enzymes",
                                "dissociation_method", "acquisition_method", "additional_terms", "modification_parameters"};
  TEST_EQUAL(columnNames(table) == names, true)

  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(ArrowIOHelpers::getColumn(table, "run_accession"), 0), "a")
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(ArrowIOHelpers::getColumn(table, "run_file_name"), 0), "a")
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(ArrowIOHelpers::getColumn(table, "file_name"), 0), "a.mzML")
  for (const char* unknown : {"fraction", "instrument", "dissociation_method", "acquisition_method", "additional_terms"})
  {
    TEST_TRUE(ArrowIOHelpers::isNull(ArrowIOHelpers::getColumn(table, unknown), 0))
  }
  // additional_terms: always null
  TEST_TRUE(ArrowIOHelpers::isNull(ArrowIOHelpers::getColumn(table, "additional_terms"), 1))
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(ArrowIOHelpers::getColumn(table, "fraction"), 1), "1")
  TEST_TRUE(ArrowIOHelpers::isNull(ArrowIOHelpers::getColumn(table, "instrument"), 1))
  TEST_EQUAL(stringList(table, "enzymes", 0) == vector<string>{"Trypsin"}, true)
  // no enzymes: null
  TEST_TRUE(ArrowIOHelpers::isNull(ArrowIOHelpers::getColumn(table, "enzymes"), 1))

  // samples: one of label LFQ, replicates null
  auto samples = structsAt(table, "samples", 0);
  ABORT_IF(samples == nullptr)
  TEST_EQUAL(samples->length(), 1)
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(samples->GetFieldByName("sample_accession"), 0), "a")
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(samples->GetFieldByName("label"), 0), "LFQ")
  TEST_TRUE(ArrowIOHelpers::isNull(samples->GetFieldByName("biological_replicate"), 0))
  TEST_TRUE(ArrowIOHelpers::isNull(samples->GetFieldByName("technical_replicate"), 0))
  samples = structsAt(table, "samples", 1);
  ABORT_IF(samples == nullptr)
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(samples->GetFieldByName("label"), 0), "TMT126")
  TEST_EQUAL(ArrowIOHelpers::getInt32Value(samples->GetFieldByName("biological_replicate"), 0), 2)
  TEST_TRUE(ArrowIOHelpers::isNull(samples->GetFieldByName("technical_replicate"), 0))

  // modification_parameters: the fixed, then the variable ones; one unknown to ModificationsDB by its name
  const auto mods = structsAt(table, "modification_parameters", 0);
  ABORT_IF(mods == nullptr)
  TEST_EQUAL(mods->length(), 3)
  ABORT_IF(mods->length() != 3)
  const auto accession = mods->GetFieldByName("accession");
  const auto name = mods->GetFieldByName("name");
  const auto fixed = mods->GetFieldByName("fixed");
  const auto position = mods->GetFieldByName("position");
  const auto target = mods->GetFieldByName("target_amino_acid");
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(accession, 0), "UNIMOD:4")
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(name, 0), "Carbamidomethyl")
  TEST_EQUAL(ArrowIOHelpers::getBoolValue(fixed, 0), true)
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(position, 0), "Anywhere")
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(target, 0), "C")
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(accession, 1), "UNIMOD:35")
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(name, 1), "Oxidation")
  TEST_EQUAL(ArrowIOHelpers::getBoolValue(fixed, 1), false)
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(position, 1), "Anywhere")
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(target, 1), "M")
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(accession, 2), "Foo (X)")
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(name, 2), "Foo (X)")
  TEST_EQUAL(ArrowIOHelpers::getBoolValue(fixed, 2), false)
  TEST_TRUE(ArrowIOHelpers::isNull(position, 2))
  TEST_TRUE(ArrowIOHelpers::isNull(target, 2))
  // no modifications: an empty list
  const auto none = structsAt(table, "modification_parameters", 1);
  ABORT_IF(none == nullptr)
  TEST_EQUAL(none->length(), 0)

  // no runs: an empty view of the same columns
  const shared_ptr<arrow::Table> empty = QPXDatasetExport::runView({});
  ABORT_IF(empty == nullptr)
  TEST_EQUAL(empty->num_rows(), 0)
  TEST_EQUAL(columnNames(empty) == names, true)
}
END_SECTION

START_SECTION((static std::shared_ptr<arrow::Table> sampleView(const std::vector<Sample>& samples)))
{
  QPXDatasetExport::Sample a;
  a.accession = "a";
  QPXDatasetExport::Sample b{"b", "Homo sapiens", "liver"};
  const shared_ptr<arrow::Table> table = QPXDatasetExport::sampleView({a, b});
  ABORT_IF(table == nullptr)
  TEST_EQUAL(table->num_rows(), 2)
  TEST_EQUAL(columnNames(table) == (vector<string>{"sample_accession", "organism", "organism_part"}), true)
  for (int i = 0; i < table->num_columns(); ++i) { TEST_EQUAL(table->schema()->field(i)->nullable(), false) }
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(ArrowIOHelpers::getColumn(table, "sample_accession"), 0), "a")
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(ArrowIOHelpers::getColumn(table, "organism"), 0), "not available")
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(ArrowIOHelpers::getColumn(table, "organism_part"), 0), "not available")
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(ArrowIOHelpers::getColumn(table, "sample_accession"), 1), "b")
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(ArrowIOHelpers::getColumn(table, "organism"), 1), "Homo sapiens")
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(ArrowIOHelpers::getColumn(table, "organism_part"), 1), "liver")
}
END_SECTION

START_SECTION((static std::shared_ptr<arrow::Table> provenanceView(const std::vector<ProvenanceStep>& steps)))
{
  const shared_ptr<arrow::Table> table = QPXDatasetExport::provenanceView(provenanceSteps());
  ABORT_IF(table == nullptr)
  TEST_EQUAL(table->num_rows(), 2)
  TEST_EQUAL(columnNames(table) == (vector<string>{"step_order", "step_category", "step_name", "tool_name", "tool_version",
                                                   "tool_uri", "parameters", "config", "output_views"}), true)
  const auto order = ArrowIOHelpers::getColumn(table, "step_order");
  TEST_EQUAL(order->type_id() == arrow::Type::INT32, true)
  TEST_EQUAL(ArrowIOHelpers::getInt32Value(order, 0), 1)
  TEST_EQUAL(ArrowIOHelpers::getInt32Value(order, 1), 2)
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(ArrowIOHelpers::getColumn(table, "step_category"), 0), "database_search")
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(ArrowIOHelpers::getColumn(table, "step_name"), 1), "psm_fdr")
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(ArrowIOHelpers::getColumn(table, "tool_name"), 0), "ProSE")
  const auto version = ArrowIOHelpers::getColumn(table, "tool_version");
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(version, 0), VersionInfo::getVersion())
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(version, 1), VersionInfo::getVersion())
  TEST_TRUE(ArrowIOHelpers::isNull(ArrowIOHelpers::getColumn(table, "tool_uri"), 0))
  TEST_TRUE(ArrowIOHelpers::isNull(ArrowIOHelpers::getColumn(table, "tool_uri"), 1))
  // the parameters in their order
  const auto parameters = structsAt(table, "parameters", 0);
  ABORT_IF(parameters == nullptr)
  TEST_EQUAL(parameters->length(), 2)
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(parameters->GetFieldByName("key"), 0), "precursor_mass_tolerance")
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(parameters->GetFieldByName("value"), 0), "10")
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(parameters->GetFieldByName("key"), 1), "enzyme")
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(parameters->GetFieldByName("value"), 1), "Trypsin")
  const auto no_parameters = structsAt(table, "parameters", 1);
  ABORT_IF(no_parameters == nullptr)
  TEST_EQUAL(no_parameters->length(), 0)
  // config "" and output_views {}: null
  const auto config = ArrowIOHelpers::getColumn(table, "config");
  TEST_STRING_EQUAL(ArrowIOHelpers::getStringValue(config, 0), "{\"a\":1}")
  TEST_TRUE(ArrowIOHelpers::isNull(config, 1))
  TEST_EQUAL(stringList(table, "output_views", 0) == vector<string>{"psm"}, true)
  TEST_TRUE(ArrowIOHelpers::isNull(ArrowIOHelpers::getColumn(table, "output_views"), 1))
}
END_SECTION

START_SECTION((static std::string formatNumber(double value, const std::string& non_finite = "null")))
{
  TEST_STRING_EQUAL(QPXDatasetExport::formatNumber(0.1), "0.1")
  TEST_STRING_EQUAL(QPXDatasetExport::formatNumber(1.0 / 3.0), "0.3333333333333333")
  TEST_STRING_EQUAL(QPXDatasetExport::formatNumber(100.0), "100")
  TEST_STRING_EQUAL(QPXDatasetExport::formatNumber(-2.5), "-2.5")
  TEST_STRING_EQUAL(QPXDatasetExport::formatNumber(std::nan("")), "null")
  TEST_STRING_EQUAL(QPXDatasetExport::formatNumber(std::nan(""), "0"), "0")
  TEST_STRING_EQUAL(QPXDatasetExport::formatNumber(INFINITY, "0"), "0")
}
END_SECTION

START_SECTION((static std::string paramToJSON(const Param& param)))
{
  Param param;
  param.setValue("a:b", 3);
  param.setValue("c", 0.1);
  param.setValue("d", "x\"y\nz");
  param.setValue("e", std::vector<std::string>{"u", "v"});
  const std::string json = QPXDatasetExport::paramToJSON(param);
  TEST_EQUAL(json.find("\"a:b\":3") != std::string::npos, true)
  TEST_EQUAL(json.find("\"c\":0.1") != std::string::npos, true)
  TEST_EQUAL(json.find("\"d\":\"x\\\"y\\nz\"") != std::string::npos, true)
  TEST_EQUAL(json.find("\"e\":[\"u\",\"v\"]") != std::string::npos, true)
  TEST_EQUAL(json.front(), '{')
  TEST_EQUAL(json.back(), '}')
  TEST_STRING_EQUAL(QPXDatasetExport::paramToJSON(Param()), "{}")
}
END_SECTION

START_SECTION((static bool writeDataset(const std::string& directory, const std::string& prefix, const std::map<std::string, View>& views)))
{
  TempDir temp;
  // the folder does not exist yet
  const std::string dir = temp.getPath() + "/runs";
  TEST_EQUAL(File::exists(dir), false)
  std::map<std::string, QPXDatasetExport::View> views;
  views["psm"].tables = {QPXFile::exportPSMsToQPXArrow(std::vector<ProteinIdentification>(), PeptideIdentificationList(), false)};
  views["run"].tables = {QPXDatasetExport::runView({runA()})};
  views["run"].metadata = {{"k", "v"}};
  views["sample"].tables = {QPXDatasetExport::sampleView({QPXDatasetExport::Sample{"a"}})};
  views["provenance"].tables = {QPXDatasetExport::provenanceView(provenanceSteps())};
  views["mz"]; // no tables: not written
  TEST_EQUAL(QPXDatasetExport::writeDataset(dir, "x", views), true)
  for (const char* view : {"psm", "run", "sample", "provenance"})
  {
    TEST_EQUAL(File::exists(dir + "/x." + view + ".parquet"), true)
  }
  TEST_EQUAL(File::exists(dir + "/x.mz.parquet"), false)
  const map<string, string> footer = footerOf(dir + "/x.run.parquet");
  TEST_EQUAL(footer.count("file_type"), 1)
  TEST_STRING_EQUAL(footer.count("file_type") ? footer.at("file_type") : "", "run_file")
  TEST_STRING_EQUAL(footer.count("k") ? footer.at("k") : "", "v")
  TEST_EQUAL(footer.count("qpx_version"), 1)
  // (the extra keys of a view are its own)
  TEST_EQUAL(footerOf(dir + "/x.provenance.parquet").count("k"), 0)

  // a psm view that cannot be written (no QPX psm table): false, and no mz view
  const std::string failed_dir = temp.getPath() + "/failed";
  Column x;
  x.name = "x";
  x.type = Column::Type::UTF8;
  x.utf8_values = {"y"};
  QPXDatasetExport::PeakFilter all;
  all.top = 0;
  std::map<std::string, QPXDatasetExport::View> broken;
  broken["psm"].tables = {ArrowIOHelpers::tableOfColumns({x}, 1)};
  broken["mz"].tables = {QPXDatasetExport::mzView(mzSpectra(false), mzIds(), "run1", all).table};
  TEST_EQUAL(QPXDatasetExport::writeDataset(failed_dir, "x", broken), false)
  TEST_EQUAL(File::exists(failed_dir + "/x.psm.parquet"), false)
  TEST_EQUAL(File::exists(failed_dir + "/x.mz.parquet"), false)
  // the mz view alone is written
  broken.erase("psm");
  TEST_EQUAL(QPXDatasetExport::writeDataset(failed_dir, "x", broken), true)
  TEST_EQUAL(File::exists(failed_dir + "/x.mz.parquet"), true)
  // a null table: false
  std::map<std::string, QPXDatasetExport::View> null_table;
  null_table["run"].tables = {nullptr};
  TEST_EQUAL(QPXDatasetExport::writeDataset(failed_dir, "y", null_table), false)
  TEST_EQUAL(File::exists(failed_dir + "/y.run.parquet"), false)
}
END_SECTION

START_SECTION(writeDataset refuses an mz view with an id twice (one run per mz view))
{
  TempDir temp;
  const std::string dir = temp.getPath();
  const auto run1 = mzTableOf({"scan=1", "scan=2"}, "run1");
  const auto run2 = mzTableOf({"scan=3", "scan=1"}, "run2");
  const auto run3 = mzTableOf({"scan=3", "scan=4"}, "run3");
  ABORT_IF(run1 == nullptr || run2 == nullptr || run3 == nullptr)
  TEST_EQUAL(ArrowIOHelpers::tableRowCount(run2), 2)

  // the mz views of two runs that both have "scan=1": false, no mz file (the error names the id); the others are written
  std::map<std::string, QPXDatasetExport::View> views;
  views["mz"].tables = {run1, run2};
  views["run"].tables = {QPXDatasetExport::runView({runA()})};
  ostringstream captured;
  OPENMS_LOG_ERROR.insert(captured);
  TEST_EQUAL(QPXDatasetExport::writeDataset(dir, "two", views), false)
  OPENMS_LOG_ERROR.remove(captured);
  TEST_EQUAL(File::exists(dir + "/two.mz.parquet"), false)
  TEST_EQUAL(File::exists(dir + "/two.run.parquet"), true)
  TEST_EQUAL(captured.str().find("'scan=1'") != string::npos, true)

  // one run's mz view is written
  views["mz"].tables = {run1};
  TEST_EQUAL(QPXDatasetExport::writeDataset(dir, "one", views), true)
  TEST_EQUAL(rowsOf(dir + "/one.mz.parquet"), 2)

  // the mz views of two runs of distinct ids are written, concatenated
  views["mz"].tables = {run1, run3};
  TEST_EQUAL(QPXDatasetExport::writeDataset(dir, "distinct", views), true)
  TEST_EQUAL(rowsOf(dir + "/distinct.mz.parquet"), 4)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
