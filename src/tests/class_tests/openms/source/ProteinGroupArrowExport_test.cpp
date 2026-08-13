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
#include <OpenMS/FORMAT/ConsensusMapArrowExport.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/FORMAT/ProteinGroupArrowExport.h>
#include <OpenMS/FORMAT/QPXValueValidation.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/File.h>

#include <arrow/api.h>

#include <optional>
#include <tuple>
#include <utility>

using namespace OpenMS;
using namespace std;

namespace
{
  /// One ProteinIdentification with a single indistinguishable group over @p paths MS runs.
  ProteinIdentification makeIdOnlyRun(const StringList& paths)
  {
    ProteinIdentification prot_id;
    prot_id.setIdentifier("PI_0");
    prot_id.setScoreType("q-value");
    prot_id.setHigherScoreBetter(false);
    prot_id.setPrimaryMSRunPath(paths);

    ProteinHit ph1;
    ph1.setAccession("PROT_A");
    ph1.setScore(0.01);
    ph1.setTargetDecoyType(ProteinHit::TargetDecoyType::TARGET);
    ProteinHit ph2;
    ph2.setAccession("PROT_B");
    ph2.setScore(0.02);
    ph2.setTargetDecoyType(ProteinHit::TargetDecoyType::TARGET);
    prot_id.setHits({ph1, ph2});

    ProteinIdentification::ProteinGroup group;
    group.accessions = {"PROT_A", "PROT_B"};
    group.probability = 0.01;
    prot_id.insertIndistinguishableProteins(group);
    return prot_id;
  }

  PeptideIdentificationList makePeptides()
  {
    PeptideIdentification pep_id;
    pep_id.setIdentifier("PI_0");
    PeptideHit hit;
    hit.setSequence(AASequence::fromString("PEPTIDEK"));
    PeptideEvidence ev;
    ev.setProteinAccession("PROT_A");
    hit.setPeptideEvidences({ev});
    pep_id.setHits({hit});
    return PeptideIdentificationList{pep_id};
  }

  /// One PSM for @p accession, attributed to file @p merge_index of a merged run.
  /// Pass merge_index < 0 to omit the metavalue entirely (an unattributable PSM).
  PeptideIdentification makeMergedPeptide(const std::string& accession, int merge_index,
                                          const std::string& sequence = "PEPTIDEK",
                                          const std::string& identifier = "PI_0")
  {
    PeptideIdentification pep_id;
    pep_id.setIdentifier(identifier);
    if (merge_index >= 0)
    {
      pep_id.setMetaValue(Constants::UserParam::ID_MERGE_INDEX, merge_index);
    }
    PeptideHit hit;
    hit.setSequence(AASequence::fromString(sequence));
    PeptideEvidence ev;
    ev.setProteinAccession(accession);
    hit.setPeptideEvidences({ev});
    pep_id.setHits({hit});
    return pep_id;
  }

  /// The grouped_runs lists of a table, in row order.
  std::vector<std::vector<std::string>> groupedRuns(const std::shared_ptr<arrow::Table>& table)
  {
    std::vector<std::vector<std::string>> out;
    auto col = std::static_pointer_cast<arrow::ListArray>(table->GetColumnByName("grouped_runs")->chunk(0));
    auto values = std::static_pointer_cast<arrow::StringArray>(col->values());
    for (int64_t i = 0; i < col->length(); ++i)
    {
      std::vector<std::string> runs;
      for (int64_t j = 0; j < col->value_length(i); ++j)
      {
        runs.push_back(values->GetString(col->value_offset(i) + j));
      }
      out.push_back(std::move(runs));
    }
    return out;
  }

  /// The single run of every row, for the identification-only overload, which never groups.
  std::vector<std::string> runNames(const std::shared_ptr<arrow::Table>& table)
  {
    std::vector<std::string> out;
    for (const auto& runs : groupedRuns(table))
    {
      out.push_back(runs.size() == 1 ? runs[0] : ListUtils::concatenate(StringList(runs.begin(), runs.end()), "+"));
    }
    return out;
  }

  /// Attach the named parallel arrays emitted by PeptideAndProteinQuant for QPX 1.1.
  void setFractionGroupQuantities(
    ProteinIdentification::ProteinGroup& group,
    const std::vector<std::tuple<UInt, UInt, float>>& quantities)
  {
    auto& abundance = group.getFloatDataArrays().emplace_back();
    abundance.setName("fraction_group_level_abundance");
    auto& integer_arrays = group.getIntegerDataArrays();
    const Size first = integer_arrays.size();
    integer_arrays.resize(first + 2);
    auto& fraction_group = integer_arrays[first];
    fraction_group.setName("fraction_group_level_fraction_group");
    auto& label = integer_arrays[first + 1];
    label.setName("fraction_group_level_label");
    for (const auto& [fraction_group_id, label_id, value] : quantities)
    {
      abundance.push_back(value);
      fraction_group.push_back(static_cast<Int>(fraction_group_id));
      label.push_back(static_cast<Int>(label_id));
    }
  }

  /// The scalar quantity of one row as an empty or single-entry {label -> value} map.
  std::map<std::string, float> intensitiesOf(const std::shared_ptr<arrow::Table>& table, int64_t row)
  {
    std::map<std::string, float> out;
    auto labels = std::static_pointer_cast<arrow::StringArray>(table->GetColumnByName("label")->chunk(0));
    auto values = std::static_pointer_cast<arrow::FloatArray>(table->GetColumnByName("intensity")->chunk(0));
    if (labels->IsNull(row) || values->IsNull(row)) { return out; }
    out[labels->GetString(row)] = values->Value(row);
    return out;
  }

  /// The (label, intensity) pair of one row, keeping the three legal shapes apart:
  ///   {nullopt, nullopt} -- identification-only: the group is not quantified anywhere,
  ///   {label,   nullopt} -- evidence in this unit, but no quantity for this label,
  ///   {label,   value}   -- quantified.
  /// intensitiesOf() collapses the first two into an empty map, so it cannot express the middle
  /// state; that is exactly what the two sections at the end of this file are about.
  using PgCell = std::pair<std::optional<std::string>, std::optional<float>>;
  PgCell cellOf(const std::shared_ptr<arrow::Table>& table, int64_t row)
  {
    auto labels = std::static_pointer_cast<arrow::StringArray>(table->GetColumnByName("label")->chunk(0));
    auto values = std::static_pointer_cast<arrow::FloatArray>(table->GetColumnByName("intensity")->chunk(0));
    PgCell out;
    if (!labels->IsNull(row)) { out.first = labels->GetString(row); }
    if (!values->IsNull(row)) { out.second = values->Value(row); }
    return out;
  }

  /// Collect scalar quantities from all rows (for tests with one grouped_runs unit).
  std::map<std::string, float> allIntensities(const std::shared_ptr<arrow::Table>& table)
  {
    std::map<std::string, float> out;
    for (int64_t row = 0; row < table->num_rows(); ++row)
    {
      const auto quantity = intensitiesOf(table, row);
      out.insert(quantity.begin(), quantity.end());
    }
    return out;
  }

  std::shared_ptr<arrow::Table> replaceColumn(
    const std::shared_ptr<arrow::Table>& table,
    const std::string& name,
    const std::shared_ptr<arrow::Array>& values)
  {
    std::vector<std::shared_ptr<arrow::ChunkedArray>> columns = table->columns();
    columns.at(static_cast<size_t>(table->schema()->GetFieldIndex(name))) =
      std::make_shared<arrow::ChunkedArray>(values);
    return arrow::Table::Make(table->schema(), std::move(columns));
  }
}

START_TEST(ProteinGroupArrowExport, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static std::shared_ptr<arrow::Table> exportToArrow(const std::vector<ProteinIdentification>&, const PeptideIdentificationList&)))
{
  // Single MS run: grouped_runs holds the origin file without path or extension, matching what
  // the psm and feature tables emit for the same run. Identification input has no design to
  // aggregate over, so the list has exactly one element.
  auto prot_id = makeIdOnlyRun({"/data/SimpleSearchEngine_1.mzML"});
  auto table = ProteinGroupArrowExport::exportToArrow({prot_id}, makePeptides());
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 1)
  TEST_TRUE(groupedRuns(table)[0] == std::vector<std::string>({"SimpleSearchEngine_1"}))
  TEST_TRUE(table->GetColumnByName("label")->chunk(0)->IsNull(0))
  TEST_TRUE(table->GetColumnByName("intensity")->chunk(0)->IsNull(0))
}
END_SECTION

START_SECTION(([EXTRA] pg value validation enforces keys, labels, and grouped_runs invariants))
{
  auto prot_id = makeIdOnlyRun({"/data/runA.mzML"});
  auto table = ProteinGroupArrowExport::exportToArrow({prot_id}, makePeptides());
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 1)

  QPXValueValidation validator(QPXValueValidation::View::PROTEIN_GROUP);
  TEST_TRUE(validator.validate(table).valid)

  auto duplicate = arrow::ConcatenateTables({table, table}).ValueOrDie();
  QPXValueValidation duplicate_validator(QPXValueValidation::View::PROTEIN_GROUP);
  TEST_FALSE(duplicate_validator.validate(duplicate).valid)

  // grouped_runs is a non-empty set, not merely a non-null Arrow list.
  auto run_values = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder duplicate_runs_builder(arrow::default_memory_pool(), run_values);
  TEST_TRUE(duplicate_runs_builder.Append().ok())
  TEST_TRUE(run_values->Append("runA").ok())
  TEST_TRUE(run_values->Append("runA").ok())
  auto duplicate_runs = replaceColumn(table, QPXPgSchema::GROUPED_RUNS,
                                      duplicate_runs_builder.Finish().ValueOrDie());
  QPXValueValidation grouped_runs_validator(QPXValueValidation::View::PROTEIN_GROUP);
  TEST_FALSE(grouped_runs_validator.validate(duplicate_runs).valid)

  auto empty_values = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder empty_runs_builder(arrow::default_memory_pool(), empty_values);
  TEST_TRUE(empty_runs_builder.Append().ok())
  auto empty_runs = replaceColumn(table, QPXPgSchema::GROUPED_RUNS,
                                  empty_runs_builder.Finish().ValueOrDie());
  QPXValueValidation empty_runs_validator(QPXValueValidation::View::PROTEIN_GROUP);
  TEST_FALSE(empty_runs_validator.validate(empty_runs).valid)

  auto blank_values = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder blank_run_builder(arrow::default_memory_pool(), blank_values);
  TEST_TRUE(blank_run_builder.Append().ok())
  TEST_TRUE(blank_values->Append(" ").ok())
  auto blank_run = replaceColumn(table, QPXPgSchema::GROUPED_RUNS,
                                 blank_run_builder.Finish().ValueOrDie());
  QPXValueValidation blank_run_validator(QPXValueValidation::View::PROTEIN_GROUP);
  TEST_FALSE(blank_run_validator.validate(blank_run).valid)

  arrow::StringBuilder anchor_builder;
  TEST_TRUE(anchor_builder.Append("").ok())
  auto empty_anchor = replaceColumn(table, QPXPgSchema::ANCHOR_PROTEIN,
                                    anchor_builder.Finish().ValueOrDie());
  QPXValueValidation anchor_validator(QPXValueValidation::View::PROTEIN_GROUP);
  TEST_FALSE(anchor_validator.validate(empty_anchor).valid)

  // A row that carries an intensity must carry a canonical label beside it.
  arrow::StringBuilder bad_label_builder;
  TEST_TRUE(bad_label_builder.Append("Arg10").ok())
  arrow::FloatBuilder intensity_builder;
  TEST_TRUE(intensity_builder.Append(42.0f).ok())
  auto bad_label = replaceColumn(table, QPXPgSchema::LABEL,
                                 bad_label_builder.Finish().ValueOrDie());
  bad_label = replaceColumn(bad_label, QPXPgSchema::INTENSITY,
                            intensity_builder.Finish().ValueOrDie());
  QPXValueValidation label_validator(QPXValueValidation::View::PROTEIN_GROUP);
  auto label_result = label_validator.validate(bad_label);
  TEST_FALSE(label_result.valid)
  TEST_TRUE(label_result.toString().find("non-canonical") != std::string::npos)

  // Scaffold for the additional_intensities check below: an identification-only row with its label
  // overwritten. Deliberately NOT asserted valid either way. pg_id is part of the identity
  // composite (pg_accessions, grouped_runs, label) and is not recomputed by replaceColumn, so this
  // row's id no longer matches its own columns -- QPXValueValidation does not re-derive ids, so
  // asserting on its validity would be asserting on something this fixture cannot speak to. That a
  // label with a null intensity is accepted is pinned instead on real exporter output, in the two
  // ConsensusMap sections at the end of this file.
  arrow::StringBuilder lfq_label_builder;
  TEST_TRUE(lfq_label_builder.Append("LFQ").ok())
  auto unmatched_label = replaceColumn(table, QPXPgSchema::LABEL,
                                       lfq_label_builder.Finish().ValueOrDie());

  // additional_intensities uses the same canonical label as the flattened pg row.
  auto pair_name_builder = std::make_shared<arrow::StringBuilder>();
  auto pair_value_builder = std::make_shared<arrow::FloatBuilder>();
  auto pair_type = arrow::struct_({
    arrow::field("intensity_name", arrow::utf8(), false),
    arrow::field("intensity_value", arrow::float32(), false)});
  auto pair_struct_builder = std::make_shared<arrow::StructBuilder>(
    pair_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{pair_name_builder, pair_value_builder});
  auto pair_list_builder = std::make_shared<arrow::ListBuilder>(
    arrow::default_memory_pool(), pair_struct_builder);
  auto extra_label_builder = std::make_shared<arrow::StringBuilder>();
  auto extra_type = arrow::struct_({
    arrow::field("label", arrow::utf8(), false),
    arrow::field("intensities", arrow::list(pair_type), false)});
  auto extra_struct_builder = std::make_shared<arrow::StructBuilder>(
    extra_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{extra_label_builder, pair_list_builder});
  arrow::ListBuilder extras_builder(arrow::default_memory_pool(), extra_struct_builder);
  TEST_TRUE(extras_builder.Append().ok())
  TEST_TRUE(extra_struct_builder->Append().ok())
  TEST_TRUE(extra_label_builder->Append("Arg10").ok())
  TEST_TRUE(pair_list_builder->Append().ok())
  TEST_TRUE(pair_struct_builder->Append().ok())
  TEST_TRUE(pair_name_builder->Append("normalized").ok())
  TEST_TRUE(pair_value_builder->Append(21.0f).ok())
  auto bad_extras = unmatched_label;
  arrow::FloatBuilder lfq_intensity_builder;
  TEST_TRUE(lfq_intensity_builder.Append(42.0f).ok())
  bad_extras = replaceColumn(bad_extras, QPXPgSchema::INTENSITY,
                             lfq_intensity_builder.Finish().ValueOrDie());
  bad_extras = replaceColumn(bad_extras, QPXPgSchema::ADDITIONAL_INTENSITIES,
                             extras_builder.Finish().ValueOrDie());
  QPXValueValidation extras_validator(QPXValueValidation::View::PROTEIN_GROUP);
  const auto extras_result = extras_validator.validate(bad_extras);
  TEST_FALSE(extras_result.valid)
  TEST_TRUE(extras_result.toString().find("additional_intensities") != std::string::npos)

  // Distinct grouped_runs keys may not overlap for one (anchor_protein, label), or the same raw
  // measurement would contribute to two protein quantities.
  auto two_rows = arrow::ConcatenateTables({table, table}).ValueOrDie();
  auto overlap_values = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder overlap_builder(arrow::default_memory_pool(), overlap_values);
  TEST_TRUE(overlap_builder.Append().ok())
  TEST_TRUE(overlap_values->Append("runA").ok())
  TEST_TRUE(overlap_values->Append("runB").ok())
  TEST_TRUE(overlap_builder.Append().ok())
  TEST_TRUE(overlap_values->Append("runB").ok())
  TEST_TRUE(overlap_values->Append("runC").ok())
  auto overlap = replaceColumn(two_rows, QPXPgSchema::GROUPED_RUNS,
                               overlap_builder.Finish().ValueOrDie());
  QPXValueValidation overlap_validator(QPXValueValidation::View::PROTEIN_GROUP);
  auto overlap_result = overlap_validator.validate(overlap);
  TEST_FALSE(overlap_result.valid)
  TEST_TRUE(overlap_result.toString().find("reuses run") != std::string::npos)

  // The table-taking writer invokes the shared validator before opening its destination.
  std::string output;
  NEW_TMP_FILE(output);
  TEST_FALSE(ProteinGroupArrowExport::exportToParquet(duplicate, output))
  TEST_FALSE(File::exists(output))
}
END_SECTION

START_SECTION(([EXTRA] exportToArrow - a merged run melts onto every run with evidence))
{
  // A group in a merged run is emitted once per run in which one of its members was identified.
  // The previous behaviour -- one row with an empty run name -- put an empty value in a QPX
  // primary-key component, which quantms silently drops.
  auto prot_id = makeIdOnlyRun({"/data/runA.mzML", "/data/runB.mzML"});
  PeptideIdentificationList peps{makeMergedPeptide("PROT_A", 0), makeMergedPeptide("PROT_B", 1)};

  auto table = ProteinGroupArrowExport::exportToArrow({prot_id}, peps);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 2)

  const auto names = runNames(table);
  const std::set<std::string> runs(names.begin(), names.end());
  TEST_TRUE(runs == std::set<std::string>({"runA", "runB"}))
  TEST_TRUE(runs.count("") == 0)
}
END_SECTION

START_SECTION(([EXTRA] exportToArrow - a merged run emits only the runs that have evidence))
{
  // Melting must not fan out across every file in the run: a group identified in one of two
  // files gets one row, not two, or the table asserts an observation that was never made.
  auto prot_id = makeIdOnlyRun({"/data/runA.mzML", "/data/runB.mzML"});
  PeptideIdentificationList peps{makeMergedPeptide("PROT_A", 0)};

  auto table = ProteinGroupArrowExport::exportToArrow({prot_id}, peps);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 1)
  TEST_STRING_EQUAL(runNames(table)[0], "runA")
}
END_SECTION

START_SECTION(([EXTRA] exportToArrow - a merged run without id_merge_index is refused))
{
  // Without the index every PSM of the run resolves to its FIRST file, so the grouped_runs key
  // would be wrong rather than missing. The psm view refuses the same input.
  auto prot_id = makeIdOnlyRun({"/data/runA.mzML", "/data/runB.mzML"});
  PeptideIdentificationList peps{makeMergedPeptide("PROT_A", -1)};

  TEST_EXCEPTION(Exception::MissingInformation,
                 ProteinGroupArrowExport::exportToArrow({prot_id}, peps))
}
END_SECTION

START_SECTION(([EXTRA] exportToArrow - a merged run with no evidence anywhere emits no row))
{
  // Nothing can key the row: several possible origins and no evidence to choose between them.
  // Skipping is the only option that does not put an empty value in the primary key.
  auto prot_id = makeIdOnlyRun({"/data/runA.mzML", "/data/runB.mzML"});
  PeptideIdentificationList peps{makeMergedPeptide("PROT_UNRELATED", 0)};

  auto table = ProteinGroupArrowExport::exportToArrow({prot_id}, peps);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 0)
}
END_SECTION

START_SECTION(([EXTRA] exportToArrow - a run without spectra_data emits no row))
{
  // No path to stem and no id_merge_index target, so the single-file escape hatch does not
  // apply either. Previously this emitted a row with an empty run name.
  auto prot_id = makeIdOnlyRun({});
  auto table = ProteinGroupArrowExport::exportToArrow({prot_id}, makePeptides());
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 0)
}
END_SECTION

START_SECTION(([EXTRA] exportToArrow - a single-file run keeps its group without peptide evidence))
{
  // The escape hatch that makes the melt non-regressive: one possible origin means no primary-key
  // ambiguity, so a group whose members carry no peptide evidence is still emitted.
  auto prot_id = makeIdOnlyRun({"/data/SimpleSearchEngine_1.mzML"});
  auto table = ProteinGroupArrowExport::exportToArrow({prot_id}, PeptideIdentificationList{});
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 1)
  TEST_STRING_EQUAL(runNames(table)[0], "SimpleSearchEngine_1")
}
END_SECTION

START_SECTION(([EXTRA] exportToArrow - peptide_counts counts sequences, not peptidoforms))
{
  // QPX splits the two: peptide_counts.unique_sequences counts SEQUENCES, while distinguishing
  // modified forms is what feature_counts.unique_features is for. The ConsensusMap overload
  // keys on toUnmodifiedString(); this overload must agree, or the same column means two
  // different things depending on which one produced the file.
  auto prot_id = makeIdOnlyRun({"/data/SimpleSearchEngine_1.mzML"});

  PeptideIdentificationList peps;
  for (const auto& seq : {"PEPTIDEK", "PEPT(Phospho)IDEK"})
  {
    PeptideIdentification pid;
    pid.setIdentifier("PI_0");
    PeptideHit hit;
    hit.setSequence(AASequence::fromString(seq));
    PeptideEvidence ev;
    ev.setProteinAccession("PROT_A");
    hit.setPeptideEvidences({ev});
    pid.setHits({hit});
    peps.push_back(pid);
  }

  auto table = ProteinGroupArrowExport::exportToArrow({prot_id}, peps);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 1)

  // Two peptidoforms of one sequence => one unique sequence.
  auto counts = std::static_pointer_cast<arrow::StructArray>(table->GetColumnByName("peptide_counts")->chunk(0));
  auto unique = std::static_pointer_cast<arrow::Int32Array>(counts->field(0));
  TEST_EQUAL(unique->Value(0), 1)
}
END_SECTION

START_SECTION(([EXTRA] exportToArrow - a merged run without groups does not abort the export))
{
  // The id_merge_index refusal is scoped to the PSMs the export actually reads. Input that
  // contributes no protein group must not be able to abort an export it has no say in --
  // e.g. BayesianProteinInference_test.idXML, which has two spectra_data paths, nine PSMs,
  // no id_merge_index and no groups.
  ProteinIdentification empty_groups;
  empty_groups.setIdentifier("PI_0");
  empty_groups.setPrimaryMSRunPath({"/data/runA.mzML", "/data/runB.mzML"});
  PeptideIdentificationList peps{makeMergedPeptide("PROT_A", -1)};

  auto table = ProteinGroupArrowExport::exportToArrow({empty_groups}, peps);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 0)
}
END_SECTION

START_SECTION(([EXTRA] exportToArrow - the merge-index refusal ignores unrelated runs))
{
  // Only the PSMs of a run that actually contributes groups are validated. A second,
  // group-less merged run in the same input must not veto the export of the first.
  auto usable = makeIdOnlyRun({"/data/SimpleSearchEngine_1.mzML"});
  ProteinIdentification unrelated;
  unrelated.setIdentifier("PI_OTHER");
  unrelated.setPrimaryMSRunPath({"/data/x.mzML", "/data/y.mzML"});

  PeptideIdentification stray;               // merged run, no id_merge_index, no groups behind it
  stray.setIdentifier("PI_OTHER");
  PeptideHit sh;
  sh.setSequence(AASequence::fromString("STRAYPEPTIDEK"));
  PeptideEvidence sev;
  sev.setProteinAccession("PROT_A");
  sh.setPeptideEvidences({sev});
  stray.setHits({sh});

  PeptideIdentificationList peps = makePeptides();
  peps.push_back(stray);

  auto table = ProteinGroupArrowExport::exportToArrow({usable, unrelated}, peps);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 1)
  TEST_STRING_EQUAL(runNames(table)[0], "SimpleSearchEngine_1")
}
END_SECTION

START_SECTION(([EXTRA] exportToArrow - two runs sharing a stem keep both rows under one key))
{
  // '/a/run.mzML' and '/b/run.mzML' both stem to 'run', so the two rows collide on the pg identity
  // (pg_accessions, grouped_runs, label) -- same group, same runs, same label, hence the same
  // pg_id. The in-memory builder keeps both so callers can inspect the source conflict; the
  // Parquet writer's shared value validator refuses them.
  auto prot_a = makeIdOnlyRun({"/a/run.mzML"});
  auto prot_b = makeIdOnlyRun({"/b/run.mzML"});
  prot_b.setIdentifier("PI_1");

  auto table = ProteinGroupArrowExport::exportToArrow({prot_a, prot_b}, makePeptides());
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 2)

  const auto names = runNames(table);
  TEST_STRING_EQUAL(names[0], "run")
  TEST_STRING_EQUAL(names[1], "run")
  auto anchors = std::static_pointer_cast<arrow::StringArray>(table->GetColumnByName("anchor_protein")->chunk(0));
  TEST_STRING_EQUAL(anchors->GetString(0), anchors->GetString(1))
}
END_SECTION

START_SECTION(([EXTRA] two groups sharing a leading protein are distinct rows, not a duplicate key))
{
  // pg_id keys on the group's FULL membership, so PROT_A;PROT_B and PROT_A;PROT_C are two
  // different groups even though both lead with PROT_A. Keying on the leader alone -- which the
  // pg view did before -- derived one id for both and made the value validator refuse the export
  // as a repeated primary key. Both must survive, with distinct ids.
  ProteinIdentification prot_id;
  prot_id.setIdentifier("PI_0");
  prot_id.setScoreType("q-value");
  prot_id.setHigherScoreBetter(false);
  prot_id.setPrimaryMSRunPath({"/data/SimpleSearchEngine_1.mzML"});

  std::vector<ProteinHit> hits;
  for (const std::string& acc : {"PROT_A", "PROT_B", "PROT_C"})
  {
    ProteinHit ph;
    ph.setAccession(acc);
    ph.setScore(0.01);
    ph.setTargetDecoyType(ProteinHit::TargetDecoyType::TARGET);
    hits.push_back(ph);
  }
  prot_id.setHits(hits);

  ProteinIdentification::ProteinGroup shared_leader_1;
  shared_leader_1.accessions = {"PROT_A", "PROT_B"};
  shared_leader_1.probability = 0.01;
  prot_id.insertIndistinguishableProteins(shared_leader_1);

  ProteinIdentification::ProteinGroup shared_leader_2;
  shared_leader_2.accessions = {"PROT_A", "PROT_C"};
  shared_leader_2.probability = 0.01;
  prot_id.insertIndistinguishableProteins(shared_leader_2);

  auto table = ProteinGroupArrowExport::exportToArrow({prot_id}, makePeptides());
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 2)

  auto anchors = std::static_pointer_cast<arrow::StringArray>(
    table->GetColumnByName("anchor_protein")->chunk(0));
  TEST_STRING_EQUAL(anchors->GetString(0), anchors->GetString(1)) // same leader ...

  auto ids = std::static_pointer_cast<arrow::Int64Array>(
    table->GetColumnByName("pg_id")->chunk(0));
  TEST_NOT_EQUAL(ids->Value(0), ids->Value(1))                    // ... different identity

  // And the value validator agrees with the identity rather than refusing the pair.
  QPXValueValidation validator(QPXValueValidation::View::PROTEIN_GROUP);
  const auto validation = validator.validate(table);
  TEST_TRUE(validation.valid)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(([EXTRA] features with divergent peptide annotations are excluded from the counts))
{
  // The pg view attributes a feature to ONE peptide. A feature whose identifications disagree
  // on the top peptide has none, so counting it would credit an accession contributed by one
  // peptide with another peptide's sequence. Uses AnnotationState, so several identifications
  // agreeing on the same peptide (FEATURE_ID_MULTIPLE_SAME -- the output of
  // IDConflictResolverAlgorithm::resolve(.., keep_matching=true)) remain acceptable.
  auto make_pid = [](const std::string& seq)
  {
    PeptideIdentification pid;
    pid.setScoreType("q-value");
    pid.setHigherScoreBetter(false);
    PeptideHit hit;
    hit.setSequence(AASequence::fromString(seq));
    hit.setCharge(2);
    hit.setScore(0.01);
    PeptideEvidence ev;
    ev.setProteinAccession("PROT_A");
    hit.setPeptideEvidences({ev});
    pid.setHits({hit});
    return pid;
  };

  auto build = [&](const PeptideIdentificationList& pids)
  {
    ConsensusMap cmap;
    cmap.setExperimentType("label-free");
    ConsensusMap::ColumnHeader ch;
    ch.filename = "/data/run1.mzML";
    ch.label = "label-free";
    cmap.getColumnHeaders()[0] = ch;

    ProteinIdentification prot;
    prot.setIdentifier("run");
    prot.setScoreType("q-value");
    prot.setPrimaryMSRunPath({"/data/run1.mzML"});
    ProteinHit ph; ph.setAccession("PROT_A");
    prot.setHits({ph});
    ProteinIdentification::ProteinGroup g;
    g.accessions = {"PROT_A"};
    g.probability = 0.01;
    prot.insertIndistinguishableProteins(g);
    cmap.setProteinIdentifications({prot});

    ConsensusFeature cf;
    cf.setMZ(500.0); cf.setRT(100.0); cf.setCharge(2);
    BaseFeature bf; bf.setIntensity(10.0f); bf.setMZ(500.0); bf.setRT(100.0); bf.setCharge(2);
    cf.insert(0, bf);
    cf.setPeptideIdentifications(pids);
    cmap.push_back(cf);
    return cmap;
  };

  // Same peptide twice (FEATURE_ID_MULTIPLE_SAME): unambiguous, exported and counted.
  {
    auto t = ProteinGroupArrowExport::exportToArrow(build({make_pid("PEPTIDEK"), make_pid("PEPTIDEK")}));
    TEST_NOT_EQUAL(t, nullptr)
    TEST_TRUE(t->num_rows() > 0)
    TEST_FALSE(t->GetColumnByName("peptide_counts")->chunk(0)->IsNull(0))
  }

  // Divergent peptides cannot be represented by a view that records one peptide per feature,
  // so the export is refused rather than silently publishing one interpretation.
  {
    auto cmap = build({make_pid("PEPTIDEK"), make_pid("OTHERPEPTIDEK")});
    TEST_EXCEPTION(Exception::IllegalArgument, ProteinGroupArrowExport::exportToArrow(cmap))
    TEST_EXCEPTION(Exception::IllegalArgument, ConsensusMapArrowExport::exportToArrow(cmap))
  }
}
END_SECTION

START_SECTION(([EXTRA] ConsensusMap overload - gene lists stay parallel to pg_accessions))
{
  // gg_accessions/gg_names must have one entry per group member, like pg_accessions and
  // pg_names, so consumers can zip them positionally. Appending only for members that carry a
  // gene metavalue shortens the lists and shifts every later gene onto the wrong protein --
  // which is what a partially gene-annotated FASTA produces.
  ConsensusMap cmap;
  cmap.setExperimentType("label-free");
  ConsensusMap::ColumnHeader ch;
  ch.filename = "/data/run1.mzML";
  ch.label = "label-free";
  cmap.getColumnHeaders()[0] = ch;

  ProteinIdentification prot;
  prot.setIdentifier("run");
  prot.setScoreType("q-value");
  prot.setPrimaryMSRunPath({"/data/run1.mzML"});
  ProteinHit pa; pa.setAccession("PROT_A");                       // no gene annotation
  ProteinHit pb; pb.setAccession("PROT_B");
  pb.setMetaValue("gene_name", "GENE_OF_B");
  prot.setHits({pa, pb});
  ProteinIdentification::ProteinGroup g;
  g.accessions = {"PROT_A", "PROT_B"};
  g.probability = 0.01;
  prot.insertIndistinguishableProteins(g);
  cmap.setProteinIdentifications({prot});

  PeptideIdentification pid;
  pid.setIdentifier("run");
  pid.setScoreType("q-value");
  pid.setHigherScoreBetter(false);
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDEK"));
  hit.setCharge(2);
  PeptideEvidence ev; ev.setProteinAccession("PROT_A");
  hit.setPeptideEvidences({ev});
  pid.setHits({hit});

  ConsensusFeature cf;
  cf.setMZ(500.0); cf.setRT(100.0); cf.setCharge(2);
  BaseFeature bf; bf.setIntensity(10.0f); bf.setMZ(500.0); bf.setRT(100.0); bf.setCharge(2);
  cf.insert(0, bf);
  cf.setPeptideIdentifications({pid});
  cmap.push_back(cf);

  auto t = ProteinGroupArrowExport::exportToArrow(cmap);
  TEST_NOT_EQUAL(t, nullptr)
  TEST_EQUAL(t->num_rows(), 1)

  auto pg_acc = std::static_pointer_cast<arrow::ListArray>(t->GetColumnByName("pg_accessions")->chunk(0));
  auto gg_acc = std::static_pointer_cast<arrow::ListArray>(t->GetColumnByName("gg_accessions")->chunk(0));
  auto gg_nam = std::static_pointer_cast<arrow::ListArray>(t->GetColumnByName("gg_names")->chunk(0));
  TEST_EQUAL(gg_acc->value_length(0), pg_acc->value_length(0))
  TEST_EQUAL(gg_nam->value_length(0), pg_acc->value_length(0))

  // The gene belongs to PROT_B, which is member 1 -- not member 0.
  auto names = std::static_pointer_cast<arrow::StringArray>(gg_nam->values());
  const int64_t off = gg_nam->value_offset(0);
  TEST_STRING_EQUAL(names->GetString(off), "")
  TEST_STRING_EQUAL(names->GetString(off + 1), "GENE_OF_B")
}
END_SECTION

START_SECTION((static std::shared_ptr<arrow::Table> exportToArrow(const ConsensusMap&, const ExperimentalDesign&)))
{
  // QPX 1.1 keys the pg view on pg_id, derived from (pg_accessions, grouped_runs, label). Two
  // fraction groups of two fractions each therefore give two LFQ rows of two runs, not four rows
  // of one.
  const std::vector<std::string> paths = {"/data/S1_F1.mzML", "/data/S1_F2.mzML",
                                          "/data/S2_F1.mzML", "/data/S2_F2.mzML"};
  ConsensusMap cmap;
  cmap.setExperimentType("label-free");
  for (Size m = 0; m < paths.size(); ++m)
  {
    ConsensusMap::ColumnHeader ch;
    ch.filename = paths[m];
    ch.label = "label-free";
    ch.setMetaValue("fraction", static_cast<int>(m % 2) + 1);          // F1, F2
    ch.setMetaValue("fraction_group", static_cast<int>(m / 2) + 1);    // S1, S2
    cmap.getColumnHeaders()[m] = ch;
  }

  ProteinIdentification prot;
  prot.setIdentifier("merged");
  prot.setScoreType("q-value");
  prot.setHigherScoreBetter(false);
  prot.setPrimaryMSRunPath(StringList(paths.begin(), paths.end()));
  ProteinHit ph; ph.setAccession("PROT_A");
  prot.setHits({ph});

  ProteinIdentification::ProteinGroup grp;
  grp.accessions = {"PROT_A"};
  grp.probability = 0.01;
  // The fraction-group arrays themselves identify this as a quantified group; no legacy
  // sample-abundance marker is required.
  setFractionGroupQuantities(grp, {{1, 1, 1000.0f}, {2, 1, 2000.0f}});
  prot.insertIndistinguishableProteins(grp);
  cmap.setProteinIdentifications({prot});

  const ExperimentalDesign design = ExperimentalDesign::fromConsensusMap(cmap);
  TEST_EQUAL(design.getNumberOfSamples(), 2)

  auto t = ProteinGroupArrowExport::exportToArrow(cmap, design);
  TEST_NOT_EQUAL(t, nullptr)

  // One row per fraction group, not per file.
  TEST_EQUAL(t->num_rows(), 2)
  const auto runs = groupedRuns(t);
  TEST_TRUE(runs[0] == std::vector<std::string>({"S1_F1", "S1_F2"}))
  TEST_TRUE(runs[1] == std::vector<std::string>({"S2_F1", "S2_F2"}))

  // ... carrying one scalar label/intensity quantity per row.
  const auto i0 = intensitiesOf(t, 0);
  const auto i1 = intensitiesOf(t, 1);
  TEST_EQUAL(i0.size(), 1)
  TEST_EQUAL(i1.size(), 1)
  TEST_REAL_SIMILAR(i0.at("LFQ"), 1000.0f)
  TEST_REAL_SIMILAR(i1.at("LFQ"), 2000.0f)

  // The convenience overload reconstructs the same design from the very same headers.
  auto t_derived = ProteinGroupArrowExport::exportToArrow(cmap);
  TEST_NOT_EQUAL(t_derived, nullptr)
  TEST_EQUAL(t_derived->num_rows(), 2)
  TEST_TRUE(groupedRuns(t_derived) == runs)
}
END_SECTION

START_SECTION(([EXTRA] ConsensusMap overload - fraction groups stay distinct when a sample repeats))
{
  // In active QPX 1.1, ExperimentalDesign::fraction_group maps exactly to grouped_runs. Reusing a
  // sample across technical-replicate fraction groups must neither join the groups nor duplicate
  // one sample-level value; PeptideAndProteinQuant preserves the two upstream quantities.
  const std::vector<std::string> paths = {"/data/S1_A.mzML", "/data/S1_B.mzML"};
  ConsensusMap cmap;
  cmap.setExperimentType("label-free");
  for (Size m = 0; m < paths.size(); ++m)
  {
    ConsensusMap::ColumnHeader ch;
    ch.filename = paths[m];
    ch.label = "label-free";
    ch.setMetaValue("fraction", static_cast<int>(m) + 1);
    ch.setMetaValue("fraction_group", static_cast<int>(m) * 3 + 1);   // 1 and 4: non-consecutive
    ch.setMetaValue("sample_name", "S1");                             // ... but ONE sample
    cmap.getColumnHeaders()[m] = ch;
  }

  ProteinIdentification prot;
  prot.setIdentifier("merged");
  prot.setScoreType("q-value");
  prot.setHigherScoreBetter(false);
  prot.setPrimaryMSRunPath(StringList(paths.begin(), paths.end()));
  ProteinHit ph; ph.setAccession("PROT_A");
  prot.setHits({ph});
  ProteinIdentification::ProteinGroup grp;
  grp.accessions = {"PROT_A"};
  grp.probability = 0.01;
  grp.getFloatDataArrays().resize(1);
  grp.getFloatDataArrays()[0].setName("abundances");
  grp.getFloatDataArrays()[0].push_back(7000.0f);      // legacy sample-level total
  setFractionGroupQuantities(grp, {{1, 1, 3000.0f}, {4, 1, 4000.0f}});
  prot.insertIndistinguishableProteins(grp);
  cmap.setProteinIdentifications({prot});

  const ExperimentalDesign design = ExperimentalDesign::fromConsensusMap(cmap);
  TEST_EQUAL(design.getNumberOfSamples(), 1)
  TEST_EQUAL(design.getNumberOfFractionGroups(), 2)   // the two disagree, which is the point

  auto t = ProteinGroupArrowExport::exportToArrow(cmap, design);
  TEST_NOT_EQUAL(t, nullptr)
  TEST_EQUAL(t->num_rows(), 2)
  TEST_TRUE(groupedRuns(t)[0] == std::vector<std::string>({"S1_A"}))
  TEST_TRUE(groupedRuns(t)[1] == std::vector<std::string>({"S1_B"}))
  TEST_REAL_SIMILAR(intensitiesOf(t, 0).at("LFQ"), 3000.0f)
  TEST_REAL_SIMILAR(intensitiesOf(t, 1).at("LFQ"), 4000.0f)
}
END_SECTION

START_SECTION(([EXTRA] ConsensusMap overload - a design row missing from one fraction is refused))
{
  // ExperimentalDesign_input_2_wrong.tsv is a shipped 4-plex whose second fraction has no label-1
  // row, and the loader accepts it on purpose ("missing fractions and wrong orders should work
  // now"). It still cannot be published: sample 1 was seen only in fractions 1 and 3, so a row
  // naming all three would claim fraction 2 aggregated a sample it never saw, and QPX resolves a
  // sample from ANY file of grouped_runs. Loader-tolerant is not the same as representable.
  //
  // Grouping by fraction_group is not enough on its own: the unit must also be rectangular.
  const std::vector<std::string> paths = {"/d/F1.mzML", "/d/F2.mzML", "/d/F3.mzML"};
  ConsensusMap cmap;
  cmap.setExperimentType("labeled_MS2");
  const char* names[4] = {"126", "127N", "127C", "128N"};
  Size mi = 0;
  for (Size f = 0; f < paths.size(); ++f)
  {
    for (Size c = 0; c < 4; ++c)
    {
      ConsensusMap::ColumnHeader h;
      h.filename = paths[f];
      h.label = std::string("tmt6plex_") + names[c];
      h.setMetaValue("channel_name", names[c]);
      h.setMetaValue("channel_id", static_cast<int>(c));
      cmap.getColumnHeaders()[mi++] = h;
    }
  }

  // The design: 4 labels over 3 fractions of ONE fraction group -- but fraction 2 has no label-1
  // row, exactly the shape of ExperimentalDesign_input_2_wrong.tsv.
  ExperimentalDesign::MSFileSection section;
  for (Size f = 0; f < paths.size(); ++f)
  {
    for (unsigned lab = 1; lab <= 4; ++lab)
    {
      if (f == 1 && lab == 1) { continue; }            // the missing row
      ExperimentalDesign::MSFileSectionEntry e;
      e.path = paths[f];
      e.fraction = static_cast<unsigned>(f) + 1;
      e.fraction_group = 1;
      e.label = lab;
      e.sample = lab - 1;
      e.sample_name = "S" + StringUtils::toStr(lab);
      section.push_back(e);
    }
  }
  ExperimentalDesign::SampleSection samples;
  for (unsigned lab = 1; lab <= 4; ++lab) { samples.addSample("S" + StringUtils::toStr(lab)); }
  ExperimentalDesign design;
  design.setMSFileSection(section);
  design.setSampleSection(samples);
  TEST_EQUAL(design.getNumberOfSamples(), 4)

  ProteinIdentification prot;
  prot.setIdentifier("merged");
  prot.setScoreType("q-value");
  prot.setPrimaryMSRunPath(StringList(paths.begin(), paths.end()));
  ProteinHit ph; ph.setAccession("PROT_A");
  prot.setHits({ph});
  ProteinIdentification::ProteinGroup grp;
  grp.accessions = {"PROT_A"};
  grp.probability = 0.01;
  grp.getFloatDataArrays().resize(1);
  grp.getFloatDataArrays()[0].setName("abundances");
  for (Size i = 0; i < 4; ++i) { grp.getFloatDataArrays()[0].push_back(100.0f * (i + 1)); }
  setFractionGroupQuantities(grp, {{1, 1, 100.0f}, {1, 2, 200.0f},
                                   {1, 3, 300.0f}, {1, 4, 400.0f}});
  prot.insertIndistinguishableProteins(grp);
  cmap.setProteinIdentifications({prot});

  TEST_EQUAL(ProteinGroupArrowExport::exportToArrow(cmap, design), nullptr)

  // Completing the missing row makes the unit rectangular, and it exports as one unit.
  ExperimentalDesign::MSFileSectionEntry missing;
  missing.path = paths[1];
  missing.fraction = 2;
  missing.fraction_group = 1;
  missing.label = 1;
  missing.sample = 0;
  missing.sample_name = "S1";
  section.push_back(missing);
  ExperimentalDesign complete;
  complete.setMSFileSection(section);
  complete.setSampleSection(samples);

  auto t = ProteinGroupArrowExport::exportToArrow(cmap, complete);
  TEST_NOT_EQUAL(t, nullptr)
  TEST_EQUAL(t->num_rows(), 4)
  TEST_TRUE(groupedRuns(t)[0] == std::vector<std::string>({"F1", "F2", "F3"}))
  const auto ints = allIntensities(t);
  TEST_EQUAL(ints.size(), 4)
  TEST_REAL_SIMILAR(ints.at("TMT126"), 100.0f)
  TEST_REAL_SIMILAR(ints.at("TMT128N"), 400.0f)
}
END_SECTION

START_SECTION(([EXTRA] ConsensusMap overload - a design channel with no column header is refused))
{
  // known_cells is per (run, LABEL), not per run. A run-level test admits a design row whose
  // channel the map has no header for; labelFor() could otherwise invent a scalar label -- "LFQ"
  // for channel 1 (silently) or the bare channel number -- and that token joins to nothing.
  //
  // Measured on real data before this was cell-level: a TMT10-shaped design against the shipped
  // 6-plex ProteinQuantifier_16_input.consensusXML wrote quantms.pg.parquet containing the
  // labels "7","8","9","10" and exited 0, while the feature view of the same directory carried
  // only the six real TMT tokens.
  const std::string path = "/d/plex.mzML";
  ConsensusMap cmap;
  cmap.setExperimentType("labeled_MS2");
  const char* names[3] = {"126", "127N", "127C"};
  for (Size c = 0; c < 3; ++c)                       // headers for channels 1..3 only
  {
    ConsensusMap::ColumnHeader h;
    h.filename = path;
    h.label = std::string("tmt6plex_") + names[c];
    h.setMetaValue("channel_name", names[c]);
    h.setMetaValue("channel_id", static_cast<int>(c));
    cmap.getColumnHeaders()[c] = h;
  }

  // The abundance array is sized from the design under test, so the two checks below differ only
  // in the declared channel count -- not in anything the sample-count guard could react to.
  auto quantify = [&](Size n_samples)
  {
    ProteinIdentification prot;
    prot.setIdentifier("merged");
    prot.setScoreType("q-value");
    prot.setPrimaryMSRunPath({path});
    ProteinHit ph; ph.setAccession("PROT_A");
    prot.setHits({ph});
    ProteinIdentification::ProteinGroup grp;
    grp.accessions = {"PROT_A"};
    grp.probability = 0.01;
    grp.getFloatDataArrays().resize(1);
    grp.getFloatDataArrays()[0].setName("abundances");
    for (Size i = 0; i < n_samples; ++i) { grp.getFloatDataArrays()[0].push_back(100.0f * (i + 1)); }
    std::vector<std::tuple<UInt, UInt, float>> quantities;
    for (UInt label = 1; label <= n_samples; ++label)
    {
      quantities.emplace_back(1, label, 100.0f * label);
    }
    setFractionGroupQuantities(grp, quantities);
    prot.insertIndistinguishableProteins(grp);
    cmap.setProteinIdentifications({prot});
  };

  auto designWithLabels = [&](unsigned n_labels)
  {
    ExperimentalDesign::MSFileSection section;
    for (unsigned lab = 1; lab <= n_labels; ++lab)
    {
      ExperimentalDesign::MSFileSectionEntry e;
      e.path = path; e.fraction = 1; e.fraction_group = 1;
      e.label = lab; e.sample = lab - 1;
      e.sample_name = "S" + StringUtils::toStr(lab);
      section.push_back(e);
    }
    ExperimentalDesign::SampleSection samples;
    for (unsigned lab = 1; lab <= n_labels; ++lab) { samples.addSample("S" + StringUtils::toStr(lab)); }
    ExperimentalDesign d;
    d.setMSFileSection(section);
    d.setSampleSection(samples);
    return d;
  };

  // Six declared channels, three headers: refused rather than inventing "4", "5", "6".
  quantify(6);
  TEST_EQUAL(ProteinGroupArrowExport::exportToArrow(cmap, designWithLabels(6)), nullptr)

  // The matching design exports, with exactly the three real tokens.
  quantify(3);
  auto t = ProteinGroupArrowExport::exportToArrow(cmap, designWithLabels(3));
  TEST_NOT_EQUAL(t, nullptr)
  TEST_EQUAL(t->num_rows(), 3)
  const auto ints = allIntensities(t);
  TEST_EQUAL(ints.size(), 3)
  TEST_TRUE(ints.count("TMT126") == 1 && ints.count("TMT127N") == 1 && ints.count("TMT127C") == 1)
}
END_SECTION

START_SECTION(([EXTRA] ConsensusMap overload - a sample reused across fraction groups remains representable))
{
  // A common reference sample may occur in several plexes. Each fraction group has its own
  // pre-aggregated label quantity, so the repeated sample does not join the plexes and is not a
  // reason to refuse the design.
  ConsensusMap cmap;
  cmap.setExperimentType("labeled_MS2");
  const char* names[2][2] = {{"126", "127N"}, {"126", "127N"}};
  for (Size f = 0; f < 2; ++f)          // two plexes
  {
    for (Size c = 0; c < 2; ++c)        // channel 1 = shared reference, channel 2 = the plex sample
    {
      ConsensusMap::ColumnHeader h;
      h.filename = "/data/plex" + StringUtils::toStr(f) + ".mzML";
      h.label = std::string("tmt6plex_") + names[f][c];
      h.setMetaValue("channel_name", names[f][c]);
      h.setMetaValue("channel_id", static_cast<int>(c));
      h.setMetaValue("fraction", 1);
      h.setMetaValue("fraction_group", static_cast<int>(f) + 1);
      // Channel 1 of BOTH plexes is the same sample; channel 2 differs per plex.
      h.setMetaValue("sample_name", c == 0 ? "REF" : ("S" + StringUtils::toStr(f)));
      cmap.getColumnHeaders()[f * 2 + c] = h;
    }
  }

  ProteinIdentification prot;
  prot.setIdentifier("merged");
  prot.setScoreType("q-value");
  prot.setPrimaryMSRunPath({"/data/plex0.mzML", "/data/plex1.mzML"});
  ProteinHit ph; ph.setAccession("PROT_A");
  prot.setHits({ph});
  ProteinIdentification::ProteinGroup grp;
  grp.accessions = {"PROT_A"};
  grp.probability = 0.01;
  grp.getFloatDataArrays().resize(1);
  grp.getFloatDataArrays()[0].setName("abundances");
  for (Size i = 0; i < 3; ++i) { grp.getFloatDataArrays()[0].push_back(100.0f * (i + 1)); }
  setFractionGroupQuantities(grp, {{1, 1, 10.0f}, {1, 2, 20.0f},
                                   {2, 1, 30.0f}, {2, 2, 40.0f}});
  prot.insertIndistinguishableProteins(grp);
  cmap.setProteinIdentifications({prot});

  const ExperimentalDesign design = ExperimentalDesign::fromConsensusMap(cmap);
  TEST_EQUAL(design.getNumberOfSamples(), 3)   // REF, S0, S1 -- REF is reached from both plexes
  auto t = ProteinGroupArrowExport::exportToArrow(cmap, design);
  TEST_NOT_EQUAL(t, nullptr)
  TEST_EQUAL(t->num_rows(), 4)
  TEST_TRUE(groupedRuns(t)[0] == std::vector<std::string>({"plex0"}))
  TEST_TRUE(groupedRuns(t)[2] == std::vector<std::string>({"plex1"}))
  TEST_REAL_SIMILAR(intensitiesOf(t, 0).at("TMT126"), 10.0f)
  TEST_REAL_SIMILAR(intensitiesOf(t, 2).at("TMT126"), 30.0f)
}
END_SECTION

START_SECTION(([EXTRA] ConsensusMap overload - design files absent from the map stay out of grouped_runs))
{
  // grouped_runs names the files that were aggregated. A design row for a file the map has no
  // column header for contributed nothing, so listing it would claim an aggregation that never
  // happened -- and would put a run in the pg view that the psm and feature views cannot join to.
  ConsensusMap cmap;
  cmap.setExperimentType("label-free");
  ConsensusMap::ColumnHeader present;
  present.filename = "/data/S1_A.mzML";
  present.label = "label-free";
  present.setMetaValue("fraction", 1);
  present.setMetaValue("fraction_group", 1);
  present.setMetaValue("sample_name", "S1");
  cmap.getColumnHeaders()[0] = present;

  // The design also lists a second fraction of the same sample that this map does not contain.
  ExperimentalDesign::MSFileSection section;
  for (const auto& [path, fraction] : std::vector<std::pair<std::string, unsigned>>{
         {"/data/S1_A.mzML", 1u}, {"/data/S1_B.mzML", 2u}})
  {
    ExperimentalDesign::MSFileSectionEntry e;
    e.path = path;
    e.fraction = fraction;
    e.fraction_group = 1;
    e.label = 1;
    e.sample = 0;
    e.sample_name = "S1";
    section.push_back(e);
  }
  ExperimentalDesign::SampleSection samples;
  samples.addSample("S1");
  ExperimentalDesign design;
  design.setMSFileSection(section);
  design.setSampleSection(samples);
  TEST_EQUAL(design.getNumberOfSamples(), 1)

  ProteinIdentification prot;
  prot.setIdentifier("merged");
  prot.setScoreType("q-value");
  prot.setPrimaryMSRunPath({"/data/S1_A.mzML"});
  ProteinHit ph; ph.setAccession("PROT_A");
  prot.setHits({ph});
  ProteinIdentification::ProteinGroup grp;
  grp.accessions = {"PROT_A"};
  grp.probability = 0.01;
  grp.getFloatDataArrays().resize(1);
  grp.getFloatDataArrays()[0].setName("abundances");
  grp.getFloatDataArrays()[0].push_back(5000.0f);
  setFractionGroupQuantities(grp, {{1, 1, 5000.0f}});
  prot.insertIndistinguishableProteins(grp);
  cmap.setProteinIdentifications({prot});

  auto t = ProteinGroupArrowExport::exportToArrow(cmap, design);
  TEST_NOT_EQUAL(t, nullptr)
  TEST_EQUAL(t->num_rows(), 1)
  TEST_TRUE(groupedRuns(t)[0] == std::vector<std::string>({"S1_A"}))   // not S1_B
}
END_SECTION

START_SECTION(([EXTRA] ConsensusMap overload - total_features counts feature rows, not features))
{
  // The feature view melts one ConsensusFeature into one row per run it was quantified in, so a
  // feature spanning two fractions of a unit IS two feature rows. Counting bare feature indices
  // made feature_counts.total_features disagree with the number of feature.parquet rows the same
  // (group, unit) has, while unique_features -- keyed on peptidoform+charge -- must NOT grow.
  const std::vector<std::string> paths = {"/data/S1_F1.mzML", "/data/S1_F2.mzML"};
  ConsensusMap cmap;
  cmap.setExperimentType("label-free");
  for (Size m = 0; m < paths.size(); ++m)
  {
    ConsensusMap::ColumnHeader ch;
    ch.filename = paths[m];
    ch.label = "label-free";
    ch.setMetaValue("fraction", static_cast<int>(m) + 1);
    ch.setMetaValue("fraction_group", 1);
    ch.setMetaValue("sample_name", "S1");
    cmap.getColumnHeaders()[m] = ch;
  }

  ProteinIdentification prot;
  prot.setIdentifier("merged");
  prot.setScoreType("q-value");
  prot.setHigherScoreBetter(false);
  prot.setPrimaryMSRunPath(StringList(paths.begin(), paths.end()));
  ProteinHit ph; ph.setAccession("PROT_A");
  prot.setHits({ph});
  ProteinIdentification::ProteinGroup grp;
  grp.accessions = {"PROT_A"};
  grp.probability = 0.01;
  grp.getFloatDataArrays().resize(1);
  grp.getFloatDataArrays()[0].setName("abundances");
  grp.getFloatDataArrays()[0].push_back(1234.0f);
  setFractionGroupQuantities(grp, {{1, 1, 1234.0f}});
  prot.insertIndistinguishableProteins(grp);
  cmap.setProteinIdentifications({prot});

  // ONE consensus feature, quantified in BOTH fractions of the unit.
  PeptideIdentification pid;
  pid.setIdentifier("merged");
  pid.setScoreType("q-value");
  pid.setHigherScoreBetter(false);
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDEK"));
  hit.setCharge(2);
  PeptideEvidence ev; ev.setProteinAccession("PROT_A");
  hit.setPeptideEvidences({ev});
  pid.setHits({hit});

  ConsensusFeature cf;
  cf.setMZ(500.0); cf.setRT(100.0); cf.setCharge(2);
  for (Size m = 0; m < paths.size(); ++m)
  {
    BaseFeature bf; bf.setIntensity(10.0f); bf.setMZ(500.0); bf.setRT(100.0); bf.setCharge(2);
    cf.insert(m, bf);
  }
  cf.setPeptideIdentifications({pid});
  cmap.push_back(cf);

  auto t = ProteinGroupArrowExport::exportToArrow(cmap);
  TEST_NOT_EQUAL(t, nullptr)
  TEST_EQUAL(t->num_rows(), 1)

  auto fc = std::static_pointer_cast<arrow::StructArray>(t->GetColumnByName("feature_counts")->chunk(0));
  auto pc = std::static_pointer_cast<arrow::StructArray>(t->GetColumnByName("peptide_counts")->chunk(0));
  TEST_EQUAL(std::static_pointer_cast<arrow::Int32Array>(fc->field(1))->Value(0), 2)  // two feature rows
  TEST_EQUAL(std::static_pointer_cast<arrow::Int32Array>(fc->field(0))->Value(0), 1)  // one peptidoform+charge
  TEST_EQUAL(std::static_pointer_cast<arrow::Int32Array>(pc->field(0))->Value(0), 1)  // one sequence
}
END_SECTION

START_SECTION(([EXTRA] ConsensusMap overload - a foreign design is refused, not read under its own numbering))
{
  // The abundance array is indexed by the sample numbering of the design that produced it, and
  // nothing in the array says which design that was. A design of a different size is the one
  // symptom that is detectable, and it must fail rather than attach one sample's quantity to
  // another sample's runs.
  ConsensusMap cmap;
  cmap.setExperimentType("label-free");
  for (Size m = 0; m < 2; ++m)
  {
    ConsensusMap::ColumnHeader ch;
    ch.filename = "/data/run" + StringUtils::toStr(m) + ".mzML";
    ch.label = "label-free";
    cmap.getColumnHeaders()[m] = ch;
  }

  ProteinIdentification prot;
  prot.setIdentifier("merged");
  prot.setScoreType("q-value");
  prot.setPrimaryMSRunPath({"/data/run0.mzML", "/data/run1.mzML"});
  ProteinHit ph; ph.setAccession("PROT_A");
  prot.setHits({ph});
  ProteinIdentification::ProteinGroup grp;
  grp.accessions = {"PROT_A"};
  grp.probability = 0.01;
  grp.getFloatDataArrays().resize(1);
  grp.getFloatDataArrays()[0].setName("abundances");
  grp.getFloatDataArrays()[0].push_back(1000.0f);   // ONE sample, but the design describes two
  setFractionGroupQuantities(grp, {{1, 1, 1000.0f}});
  prot.insertIndistinguishableProteins(grp);
  cmap.setProteinIdentifications({prot});

  TEST_EQUAL(ProteinGroupArrowExport::exportToArrow(cmap), nullptr)
}
END_SECTION

START_SECTION(([EXTRA] ConsensusMap overload - quantification with nowhere to put it is refused))
{
  // A group that HAS abundances but no quantification unit used to be written through the
  // identification-evidence branch with null intensities, losing the numbers with only a
  // LOG_WARN. A null intensities list is not "not quantified yet" -- it positively asserts the
  // group was not quantified in those runs, and is indistinguishable from the genuinely
  // unquantified rows the same branch emits.
  const std::string path = "/d/x.mzML";
  ConsensusMap cmap;
  cmap.setExperimentType("label-free");
  ConsensusMap::ColumnHeader ch;
  ch.filename = path;
  ch.label = "label-free";
  cmap.getColumnHeaders()[0] = ch;

  ProteinIdentification prot;
  prot.setIdentifier("merged");
  prot.setScoreType("q-value");
  prot.setHigherScoreBetter(false);
  prot.setPrimaryMSRunPath({path});
  ProteinHit ph; ph.setAccession("PROT_A");
  prot.setHits({ph});
  ProteinIdentification::ProteinGroup grp;
  grp.accessions = {"PROT_A"};
  grp.probability = 0.01;
  grp.getFloatDataArrays().resize(1);
  grp.getFloatDataArrays()[0].setName("abundances");
  grp.getFloatDataArrays()[0].push_back(1234.0f);
  setFractionGroupQuantities(grp, {{1, 1, 1234.0f}});
  prot.insertIndistinguishableProteins(grp);
  cmap.setProteinIdentifications({prot});
  cmap.getUnassignedPeptideIdentifications() = {makeMergedPeptide("PROT_A", 0, "PEPTIDEK", "merged")};

  ExperimentalDesign::SampleSection samples;
  samples.addSample("S0");

  // A design whose sample section is populated but whose MS file section is empty.
  ExperimentalDesign no_files;
  no_files.setMSFileSection({});
  no_files.setSampleSection(samples);
  TEST_EQUAL(no_files.getNumberOfSamples(), 1)
  TEST_EQUAL(ProteinGroupArrowExport::exportToArrow(cmap, no_files), nullptr)

  // ... and one naming a file this map has no column header for: same outcome, the abundance
  // has no unit to be attributed to.
  ExperimentalDesign::MSFileSectionEntry other;
  other.path = "/other/y.mzML"; other.fraction = 1; other.fraction_group = 1;
  other.label = 1; other.sample = 0; other.sample_name = "S0";
  ExperimentalDesign wrong_file;
  wrong_file.setMSFileSection({other});
  wrong_file.setSampleSection(samples);
  TEST_EQUAL(ProteinGroupArrowExport::exportToArrow(cmap, wrong_file), nullptr)

  // The matching design still exports the real intensity.
  ExperimentalDesign::MSFileSectionEntry match = other;
  match.path = path;
  ExperimentalDesign good;
  good.setMSFileSection({match});
  good.setSampleSection(samples);
  auto t = ProteinGroupArrowExport::exportToArrow(cmap, good);
  TEST_NOT_EQUAL(t, nullptr)
  TEST_EQUAL(t->num_rows(), 1)
  TEST_REAL_SIMILAR(intensitiesOf(t, 0).at("LFQ"), 1234.0f)
}
END_SECTION

START_SECTION(([EXTRA] ConsensusMap overload - an unquantified group melts onto units, not files))
{
  // Quantified and unquantified rows must speak one vocabulary: if a group with no abundances
  // were melted onto individual fraction files while its quantified neighbours are keyed on
  // fraction groups, one table would carry two different meanings of grouped_runs.
  const std::vector<std::string> paths = {"/data/S1_F1.mzML", "/data/S1_F2.mzML"};
  ConsensusMap cmap;
  cmap.setExperimentType("label-free");
  for (Size m = 0; m < paths.size(); ++m)
  {
    ConsensusMap::ColumnHeader ch;
    ch.filename = paths[m];
    ch.label = "label-free";
    ch.setMetaValue("fraction", static_cast<int>(m) + 1);
    ch.setMetaValue("fraction_group", 1);
    cmap.getColumnHeaders()[m] = ch;
  }

  ProteinIdentification prot;
  prot.setIdentifier("merged");
  prot.setScoreType("q-value");
  prot.setHigherScoreBetter(false);
  prot.setPrimaryMSRunPath(StringList(paths.begin(), paths.end()));
  ProteinHit ph; ph.setAccession("PROT_A");
  prot.setHits({ph});
  ProteinIdentification::ProteinGroup grp;         // no data arrays at all == unquantified
  grp.accessions = {"PROT_A"};
  grp.probability = 0.01;
  prot.insertIndistinguishableProteins(grp);
  cmap.setProteinIdentifications({prot});

  // Evidence in BOTH fractions of the one unit.
  cmap.getUnassignedPeptideIdentifications() = {makeMergedPeptide("PROT_A", 0, "PEPTIDEK", "merged"),
                                                makeMergedPeptide("PROT_A", 1, "PEPTIDER", "merged")};

  auto t = ProteinGroupArrowExport::exportToArrow(cmap);
  TEST_NOT_EQUAL(t, nullptr)
  TEST_EQUAL(t->num_rows(), 1)                     // one unit, not two files
  TEST_TRUE(groupedRuns(t)[0] == std::vector<std::string>({"S1_F1", "S1_F2"}))
  TEST_TRUE(t->GetColumnByName("label")->chunk(0)->IsNull(0))
  TEST_TRUE(t->GetColumnByName("intensity")->chunk(0)->IsNull(0))
}
END_SECTION

START_SECTION(([EXTRA] ConsensusMap overload - total_sequences counts sequences, not features))
{
  // peptide_counts.total_sequences must agree with what peptides[] sums to. Using the feature
  // count made it identical to feature_counts.total_features and contradicted peptides[] in the
  // same row -- invisible while every fixture had exactly one feature per sequence.
  ConsensusMap cmap;
  cmap.setExperimentType("label-free");
  ConsensusMap::ColumnHeader ch;
  ch.filename = "/data/run1.mzML";
  ch.label = "label-free";
  cmap.getColumnHeaders()[0] = ch;

  ProteinIdentification prot;
  prot.setIdentifier("run");
  prot.setScoreType("q-value");
  prot.setPrimaryMSRunPath({"/data/run1.mzML"});
  ProteinHit ph; ph.setAccession("PROT_A");
  prot.setHits({ph});
  ProteinIdentification::ProteinGroup g;
  g.accessions = {"PROT_A"};
  g.probability = 0.01;
  prot.insertIndistinguishableProteins(g);
  cmap.setProteinIdentifications({prot});

  // One sequence, two charge states => two consensus features, one peptide.
  for (int z : {2, 3})
  {
    PeptideIdentification pid;
    pid.setIdentifier("run");
    pid.setScoreType("q-value");
    pid.setHigherScoreBetter(false);
    PeptideHit hit;
    hit.setSequence(AASequence::fromString("PEPTIDEK"));
    hit.setCharge(z);
    PeptideEvidence ev; ev.setProteinAccession("PROT_A");
    hit.setPeptideEvidences({ev});
    pid.setHits({hit});

    ConsensusFeature cf;
    cf.setMZ(500.0); cf.setRT(100.0 + z); cf.setCharge(z);
    BaseFeature bf; bf.setIntensity(10.0f); bf.setMZ(500.0); bf.setRT(100.0 + z); bf.setCharge(z);
    cf.insert(0, bf);
    cf.setPeptideIdentifications({pid});
    cmap.push_back(cf);
  }

  auto t = ProteinGroupArrowExport::exportToArrow(cmap);
  TEST_NOT_EQUAL(t, nullptr)
  TEST_EQUAL(t->num_rows(), 1)

  auto pc = std::static_pointer_cast<arrow::StructArray>(t->GetColumnByName("peptide_counts")->chunk(0));
  auto fc = std::static_pointer_cast<arrow::StructArray>(t->GetColumnByName("feature_counts")->chunk(0));
  auto pc_unique = std::static_pointer_cast<arrow::Int32Array>(pc->field(0));
  auto pc_total  = std::static_pointer_cast<arrow::Int32Array>(pc->field(1));
  auto fc_total  = std::static_pointer_cast<arrow::Int32Array>(fc->field(1));

  TEST_EQUAL(pc_unique->Value(0), 1)   // one distinct sequence
  TEST_EQUAL(pc_total->Value(0), 1)    // peptides[] sums to 1, so total_sequences must be 1
  TEST_EQUAL(fc_total->Value(0), 2)    // two features: that is what feature_counts is for
}
END_SECTION

START_SECTION(([EXTRA] a group carrying only the legacy sample abundances is refused, not silently nulled))
{
  // A consensusXML written before the assay grain landed annotates protein groups with the
  // sample-level "abundances" array and nothing else. That group WAS quantified, at a grain this
  // version can no longer interpret. It must be refused: emitting it as identification-only would
  // put the protein in the pg view with a null intensity, indistinguishable from "not quantified",
  // and hand back a complete-looking collection with no quantities in it and exit code 0.
  const std::vector<std::string> paths = {"/data/S1.mzML", "/data/S2.mzML"};
  ConsensusMap cmap;
  cmap.setExperimentType("label-free");
  for (Size m = 0; m < paths.size(); ++m)
  {
    ConsensusMap::ColumnHeader ch;
    ch.filename = paths[m];
    ch.label = "label-free";
    ch.setMetaValue("fraction", 1);
    ch.setMetaValue("fraction_group", static_cast<int>(m) + 1);
    cmap.getColumnHeaders()[m] = ch;
  }

  ProteinIdentification prot;
  prot.setIdentifier("merged");
  prot.setScoreType("q-value");
  prot.setHigherScoreBetter(false);
  prot.setPrimaryMSRunPath(StringList(paths.begin(), paths.end()));
  ProteinHit ph; ph.setAccession("PROT_A");
  prot.setHits({ph});

  ProteinIdentification::ProteinGroup grp;
  grp.accessions = {"PROT_A"};
  grp.probability = 0.01;
  grp.getFloatDataArrays().resize(1);
  grp.getFloatDataArrays()[0].setName("abundances");   // legacy sample grain, and ONLY that
  grp.getFloatDataArrays()[0].push_back(1000.0f);
  grp.getFloatDataArrays()[0].push_back(2000.0f);
  prot.insertIndistinguishableProteins(grp);
  cmap.setProteinIdentifications({prot});

  const ExperimentalDesign design = ExperimentalDesign::fromConsensusMap(cmap);
  TEST_EQUAL(design.getNumberOfSamples(), 2)

  TEST_TRUE(ProteinGroupArrowExport::exportToArrow(cmap, design) == nullptr)

  // A group with no quantity annotation at all is a different case: identification-only data is
  // exportable, and stays so.
  ProteinIdentification id_only = prot;
  ProteinIdentification::ProteinGroup plain;
  plain.accessions = {"PROT_A"};
  plain.probability = 0.01;
  id_only.getIndistinguishableProteins().clear();
  id_only.insertIndistinguishableProteins(plain);
  ConsensusMap id_only_map = cmap;
  id_only_map.setProteinIdentifications({id_only});
  TEST_NOT_EQUAL(ProteinGroupArrowExport::exportToArrow(id_only_map, design), nullptr)
}
END_SECTION

START_SECTION(([EXTRA] ConsensusMap overload - label-free: a cell with evidence but no quantity is null, not zero))
{
  // Two fraction groups of one file each, PROT_A quantified in the first only. The second cell has
  // identification evidence but no quantity, so it keeps its label and takes a null intensity.
  // Before this, the quantity arrays had to cover every design cell exactly, so the sparse shape
  // was refused outright and PeptideAndProteinQuant materialised a 0.0 for the missing cell.
  const std::vector<std::string> paths = {"/data/S1_F1.mzML", "/data/S2_F1.mzML", "/data/S3_F1.mzML"};
  ConsensusMap cmap;
  cmap.setExperimentType("label-free");
  for (Size m = 0; m < paths.size(); ++m)
  {
    ConsensusMap::ColumnHeader ch;
    ch.filename = paths[m];
    ch.label = "label-free";
    ch.setMetaValue("fraction", 1);
    ch.setMetaValue("fraction_group", static_cast<int>(m) + 1);
    cmap.getColumnHeaders()[m] = ch;
  }

  ProteinIdentification prot;
  prot.setIdentifier("merged");
  prot.setScoreType("q-value");
  prot.setHigherScoreBetter(false);
  prot.setPrimaryMSRunPath(StringList(paths.begin(), paths.end()));
  ProteinHit ph; ph.setAccession("PROT_A");
  prot.setHits({ph});

  ProteinIdentification::ProteinGroup grp;
  grp.accessions = {"PROT_A"};
  grp.probability = 0.01;
  // Sparse on purpose: PeptideAndProteinQuant's own result map is already sparse this way -- see
  // PeptideAndProteinQuant_test.cpp asserting 'fraction_group_abundances.count(3) == 0'.
  setFractionGroupQuantities(grp, {{1, 1, 1000.0f}});
  prot.insertIndistinguishableProteins(grp);
  cmap.setProteinIdentifications({prot});

  // Evidence in the first two runs only. The third unit has neither a quantity nor evidence, so
  // the group was never seen there at all and must not acquire a row: a null intensity claims
  // "identified here but not quantified", which would be an observation that was never made.
  cmap.getUnassignedPeptideIdentifications() = {makeMergedPeptide("PROT_A", 0, "PEPTIDEK", "merged"),
                                                makeMergedPeptide("PROT_A", 1, "PEPTIDER", "merged")};

  const ExperimentalDesign design = ExperimentalDesign::fromConsensusMap(cmap);
  auto t = ProteinGroupArrowExport::exportToArrow(cmap, design);
  TEST_NOT_EQUAL(t, nullptr)

  if (t != nullptr)
  {
    TEST_EQUAL(t->num_rows(), 2)
    const auto runs = groupedRuns(t);
    std::map<std::string, PgCell> by_unit;
    for (int64_t row = 0; row < t->num_rows(); ++row) { by_unit[runs[row][0]] = cellOf(t, row); }

    // Quantified.
    TEST_EQUAL(by_unit.count("S1_F1"), 1)
    TEST_TRUE(by_unit.count("S1_F1") == 1 && by_unit.at("S1_F1").second.has_value())
    if (by_unit.count("S1_F1") == 1 && by_unit.at("S1_F1").second.has_value())
    {
      TEST_REAL_SIMILAR(*by_unit.at("S1_F1").second, 1000.0f)
    }

    // Evidence, no quantity: the label survives, the intensity is null and specifically not 0.0.
    TEST_EQUAL(by_unit.count("S2_F1"), 1)
    TEST_TRUE(by_unit.count("S2_F1") == 1 && by_unit.at("S2_F1").first.has_value())
    TEST_TRUE(by_unit.count("S2_F1") == 1 && !by_unit.at("S2_F1").second.has_value())

    // Neither quantity nor evidence: no row in any shape. A labelled null here would claim the
    // group was identified in a unit it was never seen in.
    TEST_EQUAL(by_unit.count("S3_F1"), 0)

    QPXValueValidation sparse_validator(QPXValueValidation::View::PROTEIN_GROUP);
    TEST_TRUE(sparse_validator.validate(t).valid)
  }
}
END_SECTION

START_SECTION(([EXTRA] ConsensusMap overload - isobaric: an undetected reporter channel is null, not zero))
{
  // The isobaric shape of the same rule, and the case that makes it more than bookkeeping.
  // IsobaricChannelExtractor resets every channel to 0 before searching the MS2 for its reporter
  // and inserts the channel handle even when no peak is found, so an undetected reporter arrives
  // at quantification as a stored 0.0. PeptideAndProteinQuant filters those out before aggregating
  // (only 'abundance > 0.0' contributes), so a channel whose reporters were all undetected ends up
  // with no key -- and has to be written as null, because aggregating the resulting empty list
  // returns 0.0, which is indistinguishable from a real measurement of zero.
  const std::string path = "/d/plex.mzML";
  ConsensusMap cmap;
  cmap.setExperimentType("labeled_MS2");
  const char* names[3] = {"126", "127N", "127C"};
  for (Size c = 0; c < 3; ++c)
  {
    ConsensusMap::ColumnHeader h;
    h.filename = path;
    h.label = std::string("tmt6plex_") + names[c];
    h.setMetaValue("channel_name", names[c]);
    h.setMetaValue("channel_id", static_cast<int>(c));
    cmap.getColumnHeaders()[c] = h;
  }

  ProteinIdentification prot;
  prot.setIdentifier("merged");
  prot.setScoreType("q-value");
  prot.setHigherScoreBetter(false);
  prot.setPrimaryMSRunPath({path});
  ProteinHit ph; ph.setAccession("PROT_A");
  prot.setHits({ph});
  ProteinIdentification::ProteinGroup grp;
  grp.accessions = {"PROT_A"};
  grp.probability = 0.01;
  setFractionGroupQuantities(grp, {{1, 1, 100.0f}, {1, 3, 300.0f}});   // 127N is the dead channel
  prot.insertIndistinguishableProteins(grp);
  cmap.setProteinIdentifications({prot});
  cmap.getUnassignedPeptideIdentifications() = {makeMergedPeptide("PROT_A", 0, "PEPTIDEK", "merged")};

  ExperimentalDesign::MSFileSection section;
  for (unsigned lab = 1; lab <= 3; ++lab)
  {
    ExperimentalDesign::MSFileSectionEntry e;
    e.path = path; e.fraction = 1; e.fraction_group = 1;
    e.label = lab; e.sample = lab - 1;
    e.sample_name = "S" + StringUtils::toStr(lab);
    section.push_back(e);
  }
  ExperimentalDesign::SampleSection samples;
  for (unsigned lab = 1; lab <= 3; ++lab) { samples.addSample("S" + StringUtils::toStr(lab)); }
  ExperimentalDesign design;
  design.setMSFileSection(section);
  design.setSampleSection(samples);

  auto t = ProteinGroupArrowExport::exportToArrow(cmap, design);
  TEST_NOT_EQUAL(t, nullptr)

  if (t != nullptr)
  {
    // The plex keeps its full label vocabulary: the group has evidence in this run, so all three
    // channels are statements OpenMS is entitled to make.
    TEST_EQUAL(t->num_rows(), 3)
    std::map<std::string, PgCell> by_label;
    for (int64_t row = 0; row < t->num_rows(); ++row)
    {
      const PgCell cell = cellOf(t, row);
      TEST_TRUE(cell.first.has_value())
      if (cell.first.has_value()) { by_label[*cell.first] = cell; }
    }

    TEST_EQUAL(by_label.count("TMT126"), 1)
    TEST_EQUAL(by_label.count("TMT127N"), 1)
    TEST_EQUAL(by_label.count("TMT127C"), 1)

    // The two measured channels must still carry their values. Without these, a regression that
    // nulled every channel would leave the labels and the row count intact and slip through.
    TEST_TRUE(by_label.count("TMT126") == 1 && by_label.at("TMT126").second.has_value())
    TEST_TRUE(by_label.count("TMT127C") == 1 && by_label.at("TMT127C").second.has_value())
    if (by_label.count("TMT126") == 1 && by_label.at("TMT126").second.has_value())
    {
      TEST_REAL_SIMILAR(*by_label.at("TMT126").second, 100.0f)
    }
    if (by_label.count("TMT127C") == 1 && by_label.at("TMT127C").second.has_value())
    {
      TEST_REAL_SIMILAR(*by_label.at("TMT127C").second, 300.0f)
    }

    // The dead channel: null, and specifically not 0.0.
    TEST_TRUE(by_label.count("TMT127N") == 1 && !by_label.at("TMT127N").second.has_value())

    // The sparse table is a legal QPX pg table in its own right -- this is the self-consistent
    // positive check the hand-edited fixture above cannot provide, because here pg_id was derived
    // by the exporter from the label it actually wrote.
    QPXValueValidation sparse_validator(QPXValueValidation::View::PROTEIN_GROUP);
    TEST_TRUE(sparse_validator.validate(t).valid)

    // additional_intensities stays keyed on the label, not on the intensity: an empty list on any
    // labelled row, null only on an identification-only one.
    auto extras = t->GetColumnByName(QPXPgSchema::ADDITIONAL_INTENSITIES)->chunk(0);
    for (int64_t row = 0; row < t->num_rows(); ++row) { TEST_FALSE(extras->IsNull(row)) }
  }
}
END_SECTION

START_SECTION(([EXTRA] ConsensusMap overload - the relaxed annotation check still refuses what it must))
{
  // The subset check gave up exactly one thing: a design that describes more cells than the group
  // was measured in. Everything else it used to catch it must still catch, and neither boundary is
  // covered by the acceptance sections above.
  const std::string path = "/d/one.mzML";
  ConsensusMap cmap;
  cmap.setExperimentType("label-free");
  ConsensusMap::ColumnHeader ch;
  ch.filename = path;
  ch.label = "label-free";
  ch.setMetaValue("fraction", 1);
  ch.setMetaValue("fraction_group", 1);
  cmap.getColumnHeaders()[0] = ch;

  auto mapWith = [&](const std::vector<std::tuple<UInt, UInt, float>>& quantities)
  {
    ProteinIdentification prot;
    prot.setIdentifier("merged");
    prot.setScoreType("q-value");
    prot.setHigherScoreBetter(false);
    prot.setPrimaryMSRunPath({path});
    ProteinHit ph; ph.setAccession("PROT_A");
    prot.setHits({ph});
    ProteinIdentification::ProteinGroup grp;
    grp.accessions = {"PROT_A"};
    grp.probability = 0.01;
    setFractionGroupQuantities(grp, quantities);
    prot.insertIndistinguishableProteins(grp);
    ConsensusMap out = cmap;
    out.setProteinIdentifications({prot});
    return out;
  };

  const ExperimentalDesign design = ExperimentalDesign::fromConsensusMap(cmap);

  // The one design cell, quantified: exportable.
  TEST_NOT_EQUAL(ProteinGroupArrowExport::exportToArrow(mapWith({{1, 1, 500.0f}}), design), nullptr)

  // A key the design cannot name. This is the direction the subset check keeps, and the one
  // PeptideAndProteinQuant relies on when it restricts its cell set to header-backed cells.
  TEST_EQUAL(ProteinGroupArrowExport::exportToArrow(mapWith({{2, 1, 500.0f}}), design), nullptr)

  // Arrays present but carrying nothing: a group claiming to be quantified without saying where.
  TEST_EQUAL(ProteinGroupArrowExport::exportToArrow(mapWith({}), design), nullptr)
}
END_SECTION

START_SECTION(([EXTRA] pg value validation - an intensity without a label stays rejected))
{
  // A pg intensity is only interpretable through its label: the label is what joins the value to
  // a sample (run.samples[].label), and it is a component of the pg_id identity composite. A value
  // with no label therefore cannot be attributed to anything and the row is unusable.
  //
  // This used to be enforced together with the converse, by a single "set both or neither" check
  // in QPXValueValidation. That coupling was split so a group with evidence in a quantification
  // unit but no quantity for a label can be written as a populated label with a null intensity --
  // a shape the QPX schema permits and the reference writer already handles. This section pins the
  // half that had to survive the split, and which nothing else covers on its own.
  auto prot_id = makeIdOnlyRun({"/data/runA.mzML"});
  auto table = ProteinGroupArrowExport::exportToArrow({prot_id}, makePeptides());
  TEST_NOT_EQUAL(table, nullptr)

  if (table != nullptr)
  {
    // The identification-only row has both columns null; give it an intensity and nothing else.
    arrow::FloatBuilder orphan_intensity_builder;
    TEST_TRUE(orphan_intensity_builder.Append(42.0f).ok())
    auto orphan_intensity = replaceColumn(table, QPXPgSchema::INTENSITY,
                                          orphan_intensity_builder.Finish().ValueOrDie());
    TEST_TRUE(orphan_intensity->GetColumnByName(QPXPgSchema::LABEL)->chunk(0)->IsNull(0))

    QPXValueValidation orphan_validator(QPXValueValidation::View::PROTEIN_GROUP);
    TEST_FALSE(orphan_validator.validate(orphan_intensity).valid)
  }
}
END_SECTION

END_TEST
