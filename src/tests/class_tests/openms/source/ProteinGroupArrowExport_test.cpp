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
#include <OpenMS/FORMAT/ProteinGroupArrowExport.h>
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

#include <arrow/api.h>

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

  /// The intensities of one row as {label -> value}.
  std::map<std::string, float> intensitiesOf(const std::shared_ptr<arrow::Table>& table, int64_t row)
  {
    std::map<std::string, float> out;
    auto col = std::static_pointer_cast<arrow::ListArray>(table->GetColumnByName("intensities")->chunk(0));
    if (col->IsNull(row)) { return out; }
    auto st = std::static_pointer_cast<arrow::StructArray>(col->values());
    auto labels = std::static_pointer_cast<arrow::StringArray>(st->field(0));
    auto values = std::static_pointer_cast<arrow::FloatArray>(st->field(1));
    for (int64_t j = 0; j < col->value_length(row); ++j)
    {
      const int64_t k = col->value_offset(row) + j;
      out[labels->GetString(k)] = values->Value(k);
    }
    return out;
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
  // '/a/run.mzML' and '/b/run.mzML' both stem to 'run', so the two rows collide on the
  // (anchor_protein, grouped_runs) primary key. Both are kept -- dropping one loses data and
  // there is no rule for merging two group probabilities -- and the export warns.
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
  // QPX 1.1 keys the pg view on grouped_runs: the set of raw files aggregated into ONE quantity.
  // Two fraction groups of two fractions each therefore give two rows of two runs, not four rows
  // of one -- and each row carries the SAMPLE-level abundance.
  //
  // This is not a rename. Summing the per-file cells does not reproduce the sample value:
  // PeptideAndProteinQuant sums a peptide across a fraction group's fractions before aggregating
  // them into the protein quantity, so the per-file cells are not summands (measured on BSA:
  // 304,924,936 against a true 214,289,808). The row must carry the sample value, which is why
  // the exporter reads the group's "abundances" array rather than the per-file cells.
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
  // Sample-level abundances, as PeptideAndProteinQuant writes them: one entry per design sample,
  // indexed by MSFileSectionEntry::sample. Two label-free fraction groups => two samples.
  grp.getFloatDataArrays().resize(1);
  grp.getFloatDataArrays()[0].setName("abundances");
  grp.getFloatDataArrays()[0].push_back(1000.0f);
  grp.getFloatDataArrays()[0].push_back(2000.0f);
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

  // ... carrying the sample quantity, one intensity per label.
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

START_SECTION(([EXTRA] ConsensusMap overload - the unit is the sample, not the fraction group))
{
  // A design may spread one sample over several fraction groups -- OpenMS ships one that does:
  // share/OpenMS/examples/FRACTIONS/BSA_design_onetable_nonconsec.tsv puts BSA1_F1 in fraction
  // group 1 and BSA1_F2 in fraction group 4, both sample BSA1. PeptideAndProteinQuant sums the
  // abundances per SAMPLE, so keying the row on the fraction group would emit that one quantity
  // twice, under two keys, each claiming to describe a separate measurement.
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
  grp.getFloatDataArrays()[0].push_back(7000.0f);      // one sample => one abundance
  prot.insertIndistinguishableProteins(grp);
  cmap.setProteinIdentifications({prot});

  const ExperimentalDesign design = ExperimentalDesign::fromConsensusMap(cmap);
  TEST_EQUAL(design.getNumberOfSamples(), 1)
  TEST_EQUAL(design.getNumberOfFractionGroups(), 2)   // the two disagree, which is the point

  auto t = ProteinGroupArrowExport::exportToArrow(cmap, design);
  TEST_NOT_EQUAL(t, nullptr)
  TEST_EQUAL(t->num_rows(), 1)                        // one sample => one row, not two
  TEST_TRUE(groupedRuns(t)[0] == std::vector<std::string>({"S1_A", "S1_B"}))
  TEST_REAL_SIMILAR(intensitiesOf(t, 0).at("LFQ"), 7000.0f)
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
  // Connectivity alone does not catch this -- the three runs share samples 2..4, so they form one
  // component whose label->sample relation is still a bijection. The unit must be RECTANGULAR.
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
  TEST_EQUAL(t->num_rows(), 1)
  TEST_TRUE(groupedRuns(t)[0] == std::vector<std::string>({"F1", "F2", "F3"}))
  const auto ints = intensitiesOf(t, 0);
  TEST_EQUAL(ints.size(), 4)
  TEST_REAL_SIMILAR(ints.at("TMT126"), 100.0f)
  TEST_REAL_SIMILAR(ints.at("TMT128N"), 400.0f)
}
END_SECTION

START_SECTION(([EXTRA] ConsensusMap overload - a design channel with no column header is refused))
{
  // known_cells is per (run, LABEL), not per run. A run-level test admits a design row whose
  // channel the map has no header for; labelFor() then invents an intensities[].label -- "LFQ"
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
  const auto ints = intensitiesOf(t, 0);
  TEST_EQUAL(ints.size(), 3)
  TEST_TRUE(ints.count("TMT126") == 1 && ints.count("TMT127N") == 1 && ints.count("TMT127C") == 1)
}
END_SECTION

START_SECTION(([EXTRA] ConsensusMap overload - a sample reused across units is refused))
{
  // A sample carries ONE aggregated quantity. If the design routes two quantification units to
  // the same sample -- a common reference channel shared across plexes is the realistic case --
  // that one number would be written into both rows, and a consumer summing them counts it
  // twice. The row shape cannot express a sample that spans units, so the export is refused
  // rather than publishing a plausible-looking duplicate.
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
  prot.insertIndistinguishableProteins(grp);
  cmap.setProteinIdentifications({prot});

  const ExperimentalDesign design = ExperimentalDesign::fromConsensusMap(cmap);
  TEST_EQUAL(design.getNumberOfSamples(), 3)   // REF, S0, S1 -- REF is reached from both plexes
  TEST_EQUAL(ProteinGroupArrowExport::exportToArrow(cmap, design), nullptr)
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
  TEST_TRUE(t->GetColumnByName("intensities")->chunk(0)->IsNull(0))
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

END_TEST
