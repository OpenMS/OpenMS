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
                                          const std::string& sequence = "PEPTIDEK")
  {
    PeptideIdentification pep_id;
    pep_id.setIdentifier("PI_0");
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

  /// The run_file_name values of a table, in row order.
  std::vector<std::string> runNames(const std::shared_ptr<arrow::Table>& table)
  {
    std::vector<std::string> out;
    auto col = std::static_pointer_cast<arrow::StringArray>(table->GetColumnByName("run_file_name")->chunk(0));
    for (int64_t i = 0; i < col->length(); ++i) { out.push_back(col->GetString(i)); }
    return out;
  }
}

START_TEST(ProteinGroupArrowExport, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static std::shared_ptr<arrow::Table> exportToArrow(const std::vector<ProteinIdentification>&, const PeptideIdentificationList&)))
{
  // Single MS run: run_file_name is the origin file without path or extension, matching what the
  // psm and feature tables emit for the same run.
  auto prot_id = makeIdOnlyRun({"/data/SimpleSearchEngine_1.mzML"});
  auto table = ProteinGroupArrowExport::exportToArrow({prot_id}, makePeptides());
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 1)
  auto rfn = std::static_pointer_cast<arrow::StringArray>(table->GetColumnByName("run_file_name")->chunk(0));
  TEST_STRING_EQUAL(rfn->GetString(0), "SimpleSearchEngine_1")
}
END_SECTION

START_SECTION(([EXTRA] exportToArrow - a merged run melts onto every run with evidence))
{
  // A group in a merged run is emitted once per run in which one of its members was identified.
  // The previous behaviour -- one row with an empty run_file_name -- put an empty value in a QPX
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
  // Without the index every PSM of the run resolves to its FIRST file, so the run_file_name key
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
  // apply either. Previously this emitted a row with an empty run_file_name.
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
  // (anchor_protein, run_file_name) primary key. Both are kept -- dropping one loses data and
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
