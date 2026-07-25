// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/CHEMISTRY/AASequence.h>


using namespace OpenMS;
using namespace std;


START_TEST(PeptideAndProteinQuant, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeptideAndProteinQuant* ptr = nullptr;
PeptideAndProteinQuant* nullPointer = nullptr;
START_SECTION((PeptideAndProteinQuant()))
  ptr = new PeptideAndProteinQuant();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~PeptideAndProteinQuant()))
  delete ptr;
END_SECTION

PeptideAndProteinQuant quantifier_features;
PeptideAndProteinQuant quantifier_consensus;
PeptideAndProteinQuant quantifier_identifications;
Param params;
params.setValue("top:include_all", "true");
quantifier_features.setParameters(params);
quantifier_consensus.setParameters(params);
quantifier_identifications.setParameters(params);


START_SECTION((void readQuantData(FeatureMap& features, ExperimentalDesign& ed)))
{
  FeatureMap features;
  FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("ProteinQuantifier_input.featureXML"), features);
  ExperimentalDesign design = ExperimentalDesign::fromFeatureMap(features);
  TEST_EQUAL(quantifier_features.getPeptideResults().empty(), true);
  quantifier_features.readQuantData(features, design);
  quantifier_features.quantifyPeptides();
  TEST_EQUAL(quantifier_features.getPeptideResults().empty(), false);
}
END_SECTION

START_SECTION((void readQuantData(ConsensusMap& consensus, ExperimentalDesign& ed)))
{
  ConsensusMap consensus;
  ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("ProteinQuantifier_input.consensusXML"), consensus);
  TEST_EQUAL(quantifier_consensus.getPeptideResults().empty(), true);
  ExperimentalDesign design = ExperimentalDesign::fromConsensusMap(consensus);
  quantifier_consensus.readQuantData(consensus, design);
  quantifier_consensus.quantifyPeptides();
  TEST_EQUAL(quantifier_consensus.getPeptideResults().empty(), false);
}
END_SECTION

START_SECTION((void readQuantData(vector<ProteinIdentification>& proteins, PeptideIdentificationList& peptides, ExperimentalDesign& ed)))
{
  vector<ProteinIdentification> proteins;
  PeptideIdentificationList peptides;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("ProteinQuantifier_input.idXML"), proteins, peptides);
  TEST_EQUAL(quantifier_identifications.getPeptideResults().empty(), true);
  ExperimentalDesign design = ExperimentalDesign::fromIdentifications(proteins);
  quantifier_identifications.readQuantData(proteins, peptides, design);
  quantifier_identifications.quantifyPeptides();
  TEST_EQUAL(quantifier_identifications.getPeptideResults().empty(), false);
}
END_SECTION

START_SECTION((void readQuantData(vector<ProteinIdentification>& proteins, PeptideIdentificationList& peptides, ExperimentalDesign& ed) should split a merged ID run by origin file))
{
  // Regression test for https://github.com/OpenMS/OpenMS/issues/5518:
  // a single (merged) ProteinIdentification run with several primary MS run
  // paths must be split into one sample per origin file, based on the
  // Constants::UserParam::ID_MERGE_INDEX ("id_merge_index") annotation that
  // IDMerger writes on each PeptideIdentification. Previously this code read a
  // mistyped meta key ("id_merge_idx") that never exists, so every PSM was
  // attributed to the first origin file (id_merge_index defaulting to 0).
  const std::string run_id = "merged_run";

  ProteinIdentification protein;
  protein.setIdentifier(run_id);
  protein.setPrimaryMSRunPath({"/data/fileA.mzML", "/data/fileB.mzML"});
  vector<ProteinIdentification> proteins{protein};

  // same peptide identified once per origin file
  auto make_pep = [&run_id](Size merge_idx)
  {
    PeptideHit hit;
    hit.setSequence(AASequence::fromString("PEPTIDEK"));
    hit.setCharge(2);
    hit.setScore(1.0);

    PeptideIdentification pep;
    pep.setIdentifier(run_id);
    pep.setHits({hit});
    pep.setMetaValue(Constants::UserParam::ID_MERGE_INDEX, (int)merge_idx);
    return pep;
  };

  PeptideIdentificationList peptides{make_pep(0), make_pep(1)};

  // fromIdentifications expands the two origin files into two samples
  ExperimentalDesign design = ExperimentalDesign::fromIdentifications(proteins);
  TEST_EQUAL(design.getNumberOfSamples(), 2);

  PeptideAndProteinQuant quantifier;
  quantifier.readQuantData(proteins, peptides, design);
  quantifier.quantifyPeptides();

  const auto& pep_quant = quantifier.getPeptideResults();
  const auto seq_it = pep_quant.find(AASequence::fromString("PEPTIDEK"));
  TEST_TRUE(seq_it != pep_quant.end());
  // the two PSMs must be spread over two different origin files (fraction 1)...
  TEST_EQUAL(seq_it->second.abundances.at(1).size(), 2);
  // ...and end up in two distinct samples, each with a single PSM count.
  TEST_EQUAL(seq_it->second.total_abundances.size(), 2);
  TEST_REAL_SIMILAR(seq_it->second.total_abundances.at(0), 1);
  TEST_REAL_SIMILAR(seq_it->second.total_abundances.at(1), 1);
}
END_SECTION

START_SECTION((void readQuantData(ConsensusMap& consensus, ExperimentalDesign& ed) should fail on missing file+label mapping))
{
  ConsensusMap consensus;
  consensus.setExperimentType("label-free");
  ConsensusMap::ColumnHeader ch;
  ch.filename = "not_in_design.mzML";
  ch.label = "label-free";
  consensus.getColumnHeaders()[0] = ch;

  PeptideIdentification pid;
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDE"));
  hit.setCharge(2);
  hit.setScore(1.0);
  pid.setHits({hit});

  ConsensusFeature cf;
  Peak2D p;
  p.setIntensity(42.0);
  p.setRT(1.0);
  p.setMZ(500.0);
  cf.insert(0, p, 0);
  cf.setPeptideIdentifications({pid});
  consensus.push_back(cf);

  ExperimentalDesign::MSFileSection fs;
  ExperimentalDesign::MSFileSectionEntry row;
  row.path = "in_design.mzML";
  row.label = 1;
  row.fraction = 1;
  row.fraction_group = 1;
  row.sample = 0;
  row.sample_name = "S1";
  fs.push_back(row);
  ExperimentalDesign::SampleSection ss;
  ss.addSample("S1");
  ExperimentalDesign design(fs, ss);

  PeptideAndProteinQuant quantifier;
  TEST_EXCEPTION(Exception::MissingInformation, quantifier.readQuantData(consensus, design));
}
END_SECTION

START_SECTION((void readQuantData(ConsensusMap& consensus, ExperimentalDesign& ed) should distinguish all filename+label pairs))
{
  ConsensusMap consensus;
  consensus.setExperimentType("labeled_MS2");

  ConsensusMap::ColumnHeader h0;
  h0.filename = "A1.mzML";
  h0.label = "ch2";
  h0.setMetaValue("channel_id", 1); // label = 2
  consensus.getColumnHeaders()[0] = h0;

  ConsensusMap::ColumnHeader h1;
  h1.filename = "A.mzML";
  h1.label = "ch12";
  h1.setMetaValue("channel_id", 11); // label = 12
  consensus.getColumnHeaders()[1] = h1;

  auto make_cf = [](Size map_idx, double intensity)
  {
    ConsensusFeature cf;
    Peak2D p;
    p.setIntensity(intensity);
    p.setRT(1.0 + map_idx);
    p.setMZ(400.0 + map_idx);
    cf.insert(map_idx, p, map_idx);

    PeptideIdentification pid;
    PeptideHit hit;
    hit.setSequence(AASequence::fromString("PEPTIDE"));
    hit.setCharge(2);
    hit.setScore(1.0);
    pid.setHits({hit});
    cf.setPeptideIdentifications({pid});
    return cf;
  };

  consensus.push_back(make_cf(0, 10.0));
  consensus.push_back(make_cf(1, 20.0));

  ExperimentalDesign::MSFileSection fs;
  ExperimentalDesign::MSFileSectionEntry r1;
  r1.path = "/tmp/A1.mzML";
  r1.label = 2;
  r1.fraction = 1;
  r1.fraction_group = 1;
  r1.sample = 0;
  r1.sample_name = "S1";
  fs.push_back(r1);

  ExperimentalDesign::MSFileSectionEntry r2;
  r2.path = "/tmp/A.mzML";
  r2.label = 12;
  r2.fraction = 2;
  r2.fraction_group = 2;
  r2.sample = 1;
  r2.sample_name = "S2";
  fs.push_back(r2);

  ExperimentalDesign::SampleSection ss;
  ss.addSample("S1");
  ss.addSample("S2");
  ExperimentalDesign design(fs, ss);

  PeptideAndProteinQuant quantifier;
  quantifier.readQuantData(consensus, design);
  quantifier.quantifyPeptides();

  const auto& pep_quant = quantifier.getPeptideResults();
  const auto seq_it = pep_quant.find(AASequence::fromString("PEPTIDE"));
  TEST_TRUE(seq_it != pep_quant.end());
  TEST_EQUAL(seq_it->second.abundances.size(), 2);
}
END_SECTION

START_SECTION((void quantifyPeptides(const PeptideIdentificationList& peptides = PeptideIdentificationList())))
{
  NOT_TESTABLE // tested together with the "readQuantData" methods
}
END_SECTION

START_SECTION((void quantifyProteins(const ProteinIdentification& proteins = ProteinIdentification())))
{
  TEST_EQUAL(quantifier_features.getProteinResults().empty(), true);
  quantifier_features.quantifyProteins();
  TEST_EQUAL(quantifier_features.getProteinResults().empty(), false);

  TEST_EQUAL(quantifier_consensus.getProteinResults().empty(), true);
  quantifier_consensus.quantifyProteins();
  TEST_EQUAL(quantifier_consensus.getProteinResults().empty(), false);

  TEST_EQUAL(quantifier_identifications.getProteinResults().empty(), true);
  quantifier_identifications.quantifyProteins();
  TEST_EQUAL(quantifier_identifications.getProteinResults().empty(), false);
}
END_SECTION

START_SECTION((const Statistics& getStatistics()))
{
  PeptideAndProteinQuant::Statistics stats;

  stats = quantifier_features.getStatistics();
  TEST_EQUAL(stats.n_samples, 1);
  TEST_EQUAL(stats.quant_proteins, 2);
  TEST_EQUAL(stats.too_few_peptides, 1);
  TEST_EQUAL(stats.quant_peptides, 5);
  TEST_EQUAL(stats.total_peptides, 7);
  TEST_EQUAL(stats.quant_features, 7);
  TEST_EQUAL(stats.total_features, 8);
  TEST_EQUAL(stats.blank_features, 0);
  TEST_EQUAL(stats.ambig_features, 1);

  stats = quantifier_consensus.getStatistics();
  TEST_EQUAL(stats.n_samples, 3);
  TEST_EQUAL(stats.quant_proteins, 1);
  TEST_EQUAL(stats.too_few_peptides, 0);
  TEST_EQUAL(stats.quant_peptides, 4);
  TEST_EQUAL(stats.total_peptides, 4);
  TEST_EQUAL(stats.quant_features, 9);
  TEST_EQUAL(stats.total_features, 9);
  TEST_EQUAL(stats.blank_features, 0);
  TEST_EQUAL(stats.ambig_features, 0);

  stats = quantifier_identifications.getStatistics();
  TEST_EQUAL(stats.n_samples, 2);
  TEST_EQUAL(stats.quant_proteins, 10);
  TEST_EQUAL(stats.too_few_peptides, 10);
  TEST_EQUAL(stats.quant_peptides, 13); // one decoy peptide is not quantified
  TEST_EQUAL(stats.total_peptides, 14);
  TEST_EQUAL(stats.quant_features, 17); // feature with a decoy peptide is not quantified
  TEST_EQUAL(stats.total_features, 18);
  TEST_EQUAL(stats.blank_features, 0);
  TEST_EQUAL(stats.ambig_features, 0);
}
END_SECTION

START_SECTION((const PeptideQuant& getPeptideResults()))
{
  PeptideAndProteinQuant::PeptideQuant pep_quant;
  PeptideAndProteinQuant::PeptideData pep_data;

  pep_quant = quantifier_features.getPeptideResults();
  TEST_EQUAL(pep_quant.size(), 7);
  pep_data = pep_quant[AASequence::fromString("AAAAA")];
  TEST_EQUAL(pep_data.abundances.size(), 1);
  TEST_EQUAL(pep_data.abundances[1].size(), 1);
  TEST_EQUAL(pep_data.total_abundances.size(), 1);
  TEST_REAL_SIMILAR(pep_data.total_abundances[0], 3333);
  TEST_EQUAL(pep_data.accessions.size(), 1);
  TEST_EQUAL(pep_data.psm_count, 2);
  pep_data = pep_quant[AASequence::fromString("CCCCC")];
  TEST_EQUAL(pep_data.abundances.size(), 1); // one fraction
  TEST_EQUAL(pep_data.abundances[1].size(), 1); // one file
  auto& map_file_to_charges = *pep_data.abundances[1].begin();
  TEST_EQUAL(map_file_to_charges.second.size(), 2); // two charges

  
  TEST_EQUAL(pep_data.total_abundances.size(), 1);
  TEST_REAL_SIMILAR(pep_data.total_abundances[0], 7777);
  TEST_EQUAL(pep_data.accessions.size(), 1);
  TEST_EQUAL(pep_data.psm_count, 2);
  pep_data = pep_quant[AASequence::fromString("EEEEE")];
  TEST_EQUAL(pep_data.abundances.size(), 0); // it is the second best hit, so it will not be counted
  TEST_EQUAL(pep_data.total_abundances.size(), 0);
  TEST_EQUAL(pep_data.accessions.size(), 1);
  TEST_EQUAL(pep_data.psm_count, 1);
  pep_data = pep_quant[AASequence::fromString("GGGGG")];
  TEST_EQUAL(pep_data.abundances.size(), 1); // one fraction
  TEST_EQUAL(pep_data.abundances[1].size(), 1); // one file
  TEST_EQUAL((*pep_data.abundances[1].begin()).second.size(), 1); // two charges  

  TEST_EQUAL(pep_data.total_abundances.size(), 1);
  TEST_REAL_SIMILAR(pep_data.total_abundances[0], 7777);
  TEST_EQUAL(pep_data.accessions.size(), 2);
  TEST_EQUAL(pep_data.psm_count, 1);

  pep_quant = quantifier_consensus.getPeptideResults();
  TEST_EQUAL(pep_quant.size(), 4);
  pep_data = pep_quant[AASequence::fromString("AAAK")];
  TEST_EQUAL(pep_data.abundances.size(), 1); // one fraction
  TEST_EQUAL(pep_data.abundances[1].size(), 2); // two files
  TEST_EQUAL((*pep_data.abundances[1].begin()).second.size(), 1); // one charge

  TEST_EQUAL(pep_data.total_abundances.size(), 2);
  TEST_REAL_SIMILAR(pep_data.total_abundances[0], 1000);
  TEST_REAL_SIMILAR(pep_data.total_abundances[2], 1000);
  TEST_EQUAL(pep_data.accessions.size(), 1);
  TEST_EQUAL(pep_data.psm_count, 1);
  pep_data = pep_quant[AASequence::fromString("CCCK")];
  TEST_EQUAL(pep_data.abundances.size(), 1); // one fraction
  TEST_EQUAL(pep_data.abundances[1].size(), 2); // two files
  TEST_EQUAL((*pep_data.abundances[1].begin()).second.size(), 1); // one charge

  TEST_EQUAL(pep_data.total_abundances.size(), 2);
  TEST_REAL_SIMILAR(pep_data.total_abundances[0], 200);
  TEST_REAL_SIMILAR(pep_data.total_abundances[1], 200);
  TEST_EQUAL(pep_data.accessions.size(), 1);
  TEST_EQUAL(pep_data.psm_count, 1);
  pep_data = pep_quant[AASequence::fromString("EEEK")];
  TEST_EQUAL(pep_data.abundances.size(), 1); // one fraction
  TEST_EQUAL(pep_data.abundances[1].size(), 3); // three files
  TEST_EQUAL((*pep_data.abundances[1].begin()).second.size(), 1); // one charge

  TEST_EQUAL(pep_data.total_abundances.size(), 3);
  TEST_REAL_SIMILAR(pep_data.total_abundances[0], 30);
  TEST_REAL_SIMILAR(pep_data.total_abundances[1], 30);
  TEST_REAL_SIMILAR(pep_data.total_abundances[2], 30);
  TEST_EQUAL(pep_data.accessions.size(), 1);
  TEST_EQUAL(pep_data.psm_count, 1);
  pep_data = pep_quant[AASequence::fromString("GGG")];
  TEST_EQUAL(pep_data.abundances.size(), 1); // one fraction
  TEST_EQUAL(pep_data.abundances[1].size(), 2); // two files
  TEST_EQUAL((*pep_data.abundances[1].begin()).second.size(), 1); // one charge

  TEST_EQUAL(pep_data.total_abundances.size(), 2);
  TEST_REAL_SIMILAR(pep_data.total_abundances[0], 4);
  TEST_REAL_SIMILAR(pep_data.total_abundances[1], 4);
  TEST_EQUAL(pep_data.accessions.size(), 1);
  TEST_EQUAL(pep_data.psm_count, 1);
}
END_SECTION

START_SECTION((const ProteinQuant& getProteinResults()))
{
  PeptideAndProteinQuant::ProteinQuant prot_quant;
  PeptideAndProteinQuant::ProteinData prot_data;

  prot_quant = quantifier_features.getProteinResults();
  TEST_EQUAL(prot_quant.size(), 2);
  prot_data = prot_quant["Protein0"];
  TEST_EQUAL(prot_data.peptide_abundances.size(), 3);
  TEST_EQUAL(prot_data.total_abundances.size(), 1);
  TEST_REAL_SIMILAR(prot_data.total_abundances[0], 4711);
  TEST_EQUAL(prot_data.psm_count, 6);
  prot_data = prot_quant["Protein1"];
  TEST_EQUAL(prot_data.peptide_abundances.size(), 1);
  TEST_EQUAL(prot_data.total_abundances.size(), 1);
  TEST_REAL_SIMILAR(prot_data.total_abundances[0], 8888);
  TEST_EQUAL(prot_data.psm_count, 2);

  prot_quant = quantifier_consensus.getProteinResults();
  TEST_EQUAL(prot_quant.size(), 1);
  prot_data = prot_quant["Protein"];
  TEST_EQUAL(prot_data.peptide_abundances.size(), 4);
  TEST_EQUAL(prot_data.total_abundances.size(), 3);
  TEST_REAL_SIMILAR(prot_data.total_abundances[0], 200);
  TEST_REAL_SIMILAR(prot_data.total_abundances[1], 30);
  TEST_REAL_SIMILAR(prot_data.total_abundances[2], 515);
  TEST_EQUAL(prot_data.psm_count, 4);
}
END_SECTION

START_SECTION(([PeptideAndProteinQuant::PeptideData] PeptideData()))
{
  PeptideAndProteinQuant::PeptideData data;
  TEST_EQUAL(data.abundances.empty(), true);
  TEST_EQUAL(data.total_abundances.empty(), true);
  TEST_EQUAL(data.accessions.empty(), true);
  TEST_EQUAL(data.psm_count, 0);
}
END_SECTION

START_SECTION(([PeptideAndProteinQuant::ProteinData] ProteinData()))
{
  PeptideAndProteinQuant::ProteinData data;
  TEST_EQUAL(data.peptide_abundances.empty(), true);
  TEST_EQUAL(data.total_abundances.empty(), true);
  TEST_EQUAL(data.psm_count, 0);
}
END_SECTION

START_SECTION(([PeptideAndProteinQuant::Statistics] Statistics()))
{
  PeptideAndProteinQuant::Statistics stats;
  TEST_EQUAL(stats.n_samples, 0);
  TEST_EQUAL(stats.quant_proteins, 0);
  TEST_EQUAL(stats.too_few_peptides, 0);
  TEST_EQUAL(stats.quant_peptides, 0);
  TEST_EQUAL(stats.total_peptides, 0);
  TEST_EQUAL(stats.quant_features, 0);
  TEST_EQUAL(stats.total_features, 0);
  TEST_EQUAL(stats.blank_features, 0);
  TEST_EQUAL(stats.ambig_features, 0);
}
END_SECTION

// testing various averaging strategies
START_SECTION((const ProteinQuant& getProteinResults()))
{
  FeatureMap f;
  FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("ProteinQuantifier_input.featureXML"), f);
  
  PeptideAndProteinQuant quantifier;
  PeptideAndProteinQuant::ProteinQuant quant;
  PeptideAndProteinQuant::ProteinData protein;
  Param parameters;
  parameters.setValue("top:N", 0);
  
  parameters.setValue("top:aggregate", "median");
  quantifier.setParameters(parameters);
  ExperimentalDesign ed = ExperimentalDesign::fromFeatureMap(f);
  quantifier.readQuantData(f, ed);
  quantifier.quantifyPeptides();
  quantifier.quantifyProteins();
  quant = quantifier.getProteinResults();
  protein = quant["Protein0"];
  TEST_REAL_SIMILAR(protein.total_abundances[0], 4711);

  parameters.setValue("top:aggregate", "mean");
  quantifier.setParameters(parameters);
  quantifier.readQuantData(f, ed);
  quantifier.quantifyPeptides();
  quantifier.quantifyProteins();
  quant = quantifier.getProteinResults();
  protein = quant["Protein0"];
  TEST_REAL_SIMILAR(protein.total_abundances[0], 5273.666666);

  parameters.setValue("top:aggregate", "weighted_mean");
  quantifier.setParameters(parameters);
  quantifier.readQuantData(f, ed);
  quantifier.quantifyPeptides();
  quantifier.quantifyProteins();
  quant = quantifier.getProteinResults();
  protein = quant["Protein0"];
  TEST_REAL_SIMILAR(protein.total_abundances[0], 5927.82624360028);

  parameters.setValue("top:aggregate", "sum");
  quantifier.setParameters(parameters);
  quantifier.readQuantData(f, ed);
  quantifier.quantifyPeptides();
  quantifier.quantifyProteins();
  quant = quantifier.getProteinResults();
  protein = quant["Protein0"];
  TEST_REAL_SIMILAR(protein.total_abundances[0], 15821);
}
END_SECTION

// iBAQ test
START_SECTION((const ProteinQuant& getProteinResults()))
{
  PeptideAndProteinQuant quantifier;
  PeptideAndProteinQuant::ProteinQuant quant;
  PeptideAndProteinQuant::ProteinData protein;

  Param parameters = quantifier.getDefaults();
  parameters.setValue("method", "iBAQ");
  quantifier.setParameters(parameters);

  ConsensusMap consensus;
  ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("ProteinQuantifier_input.consensusXML"), consensus);
  ExperimentalDesign ed = ExperimentalDesign::fromConsensusMap(consensus);
  ProteinIdentification proteins_ = consensus.getProteinIdentifications()[0];
  quantifier.readQuantData(consensus, ed);
  quantifier.quantifyPeptides();
  quantifier.quantifyProteins(proteins_);

  quant = quantifier.getProteinResults();
  protein = quant["Protein"];
  TEST_REAL_SIMILAR(protein.total_abundances[0], 308.5);
  TEST_REAL_SIMILAR(protein.total_abundances[1], 58.5);
  TEST_REAL_SIMILAR(protein.total_abundances[2], 257.5);
}
END_SECTION


/////////////////////////////////////////////////////////////
// File+channel level ("detailed") protein abundances.
//
// ProteinData::channel_level_abundances is the source of the
// "-file_and_channel_level_output true" columns of ProteinQuantifier and of the
// per-run intensities in the QPX pg.parquet export. It has to follow the *same*
// peptide-level policy as ProteinData::total_abundances, which is computed from
// the peptidoform-merged ProteinData::peptide_abundances:
//   - all peptidoforms of an unmodified sequence are summed into one peptide,
//   - all charge states of a peptide are summed into that same one peptide,
//   - "top:N" therefore counts peptides, not peptidoforms and not charge states,
//   - a (file, channel) cell that does not reach "top:N" peptides is absent,
//     not filled with a raw sum over everything.
/////////////////////////////////////////////////////////////

// one quantified peptide observation in the synthetic input below
struct DetailedEntry
{
  Size file_index;       ///< index into the file list
  const char* sequence;  ///< (modified) peptide sequence
  Int charge;
  double intensity;
};

// Build a label-free ConsensusMap (one column header and one sample per file)
// with a matching ExperimentalDesign. All peptides are assigned to @p accession.
auto make_detailed_input = [](const vector<std::string>& files,
                              const vector<DetailedEntry>& entries,
                              const std::string& accession,
                              ConsensusMap& consensus,
                              ExperimentalDesign& design)
{
  consensus.setExperimentType("label-free");
  for (Size i = 0; i < files.size(); ++i)
  {
    ConsensusMap::ColumnHeader h;
    h.filename = files[i] + ".mzML";
    h.label = "label-free";
    consensus.getColumnHeaders()[i] = h;
  }

  double rt = 1.0;
  for (const auto& e : entries)
  {
    PeptideEvidence ev;
    ev.setProteinAccession(accession);

    PeptideHit hit;
    hit.setSequence(AASequence::fromString(e.sequence));
    hit.setCharge(e.charge);
    hit.setScore(1.0);
    hit.setPeptideEvidences({ev});

    PeptideIdentification pid;
    pid.setHits({hit});

    ConsensusFeature cf;
    Peak2D p;
    p.setIntensity(e.intensity);
    p.setRT(rt);
    p.setMZ(500.0 + rt);
    rt += 1.0;
    cf.insert(e.file_index, p, e.file_index);
    cf.setPeptideIdentifications({pid});
    consensus.push_back(cf);
  }

  ExperimentalDesign::MSFileSection fs;
  ExperimentalDesign::SampleSection ss;
  for (Size i = 0; i < files.size(); ++i)
  {
    ExperimentalDesign::MSFileSectionEntry row;
    row.path = files[i] + ".mzML";
    row.label = 1;
    row.fraction = 1;
    row.fraction_group = i + 1;
    row.sample = i;
    row.sample_name = files[i];
    fs.push_back(row);
    ss.addSample(files[i]);
  }
  design = ExperimentalDesign(fs, ss);
};

// Build a FRACTIONATED, multi-channel ConsensusMap: a single fraction group whose files are
// consecutive fractions, each measured in @p n_channels channels. A sample is a
// (file, label) pair, so here one sample spans ALL files at a given channel - i.e. one sample
// corresponds to several (file, channel) cells. @p file_index/@p channel index the entry.
struct FractionEntry
{
  Size file_index;
  UInt channel;          ///< 1-based label/channel
  const char* sequence;
  Int charge;
  double intensity;
};

auto make_fractionated_input = [](const vector<std::string>& files,
                                  UInt n_channels,
                                  const vector<FractionEntry>& entries,
                                  const std::string& accession,
                                  ConsensusMap& consensus,
                                  ExperimentalDesign& design)
{
  consensus.setExperimentType("labeled_MS2");
  // map index = file_index * n_channels + (channel - 1)
  for (Size f = 0; f < files.size(); ++f)
  {
    for (UInt c = 1; c <= n_channels; ++c)
    {
      ConsensusMap::ColumnHeader h;
      h.filename = files[f] + ".mzML";
      h.label = "ch" + std::to_string(c);
      h.setMetaValue("channel_id", (int)(c - 1)); // getLabelAsUInt() returns channel_id + 1
      consensus.getColumnHeaders()[f * n_channels + (c - 1)] = h;
    }
  }

  double rt = 1.0;
  for (const auto& e : entries)
  {
    PeptideEvidence ev;
    ev.setProteinAccession(accession);

    PeptideHit hit;
    hit.setSequence(AASequence::fromString(e.sequence));
    hit.setCharge(e.charge);
    hit.setScore(1.0);
    hit.setPeptideEvidences({ev});

    PeptideIdentification pid;
    pid.setHits({hit});

    const Size map_index = e.file_index * n_channels + (e.channel - 1);
    ConsensusFeature cf;
    Peak2D p;
    p.setIntensity(e.intensity);
    p.setRT(rt);
    p.setMZ(500.0 + rt);
    rt += 1.0;
    cf.insert(map_index, p, map_index);
    cf.setPeptideIdentifications({pid});
    consensus.push_back(cf);
  }

  ExperimentalDesign::MSFileSection fs;
  ExperimentalDesign::SampleSection ss;
  for (UInt c = 1; c <= n_channels; ++c)
  {
    ss.addSample("S" + std::to_string(c));
  }
  for (Size f = 0; f < files.size(); ++f)
  {
    for (UInt c = 1; c <= n_channels; ++c)
    {
      ExperimentalDesign::MSFileSectionEntry row;
      row.path = files[f] + ".mzML";
      row.label = c;
      row.fraction = f + 1;   // consecutive fractions ...
      row.fraction_group = 1; // ... of ONE fraction group
      row.sample = c - 1;     // one sample per channel, spanning all fractions
      row.sample_name = "S" + std::to_string(c);
      fs.push_back(row);
    }
  }
  design = ExperimentalDesign(fs, ss);
};

START_SECTION((const ProteinQuant& getProteinResults() fractionated: file+channel cells decompose the sample for sum))
{
  // One sample spans two fraction files, so the (file, channel) cells are a decomposition of
  // the sample - but only exactly so for top:N 0 with "sum", where per-file aggregation
  // commutes with aggregation across fractions.
  //   fileA ch1: PEPTIDEK 100 + AAAAAK 200 = 300     fileB ch1: CCCCCK 300 + PEPTIDEK 50 = 350
  //   fileA ch2: PEPTIDEK  10 + AAAAAK  20 =  30     fileB ch2: CCCCCK  30 + PEPTIDEK  5 =  35
  ConsensusMap consensus;
  ExperimentalDesign design;
  make_fractionated_input({"fileA", "fileB"}, 2,
                          {{0, 1, "PEPTIDEK", 2, 100.0}, {0, 2, "PEPTIDEK", 2, 10.0},
                           {0, 1, "AAAAAK",   2, 200.0}, {0, 2, "AAAAAK",   2, 20.0},
                           {1, 1, "CCCCCK",   2, 300.0}, {1, 2, "CCCCCK",   2, 30.0},
                           {1, 1, "PEPTIDEK", 2,  50.0}, {1, 2, "PEPTIDEK", 2,  5.0}},
                          "Prot", consensus, design);

  TEST_EQUAL(design.getNumberOfSamples(), 2);
  TEST_EQUAL(design.getNumberOfFractions(), 2);

  PeptideAndProteinQuant quantifier;
  Param p = quantifier.getDefaults();
  p.setValue("top:N", 0);
  p.setValue("top:aggregate", "sum");
  quantifier.setParameters(p);

  quantifier.readQuantData(consensus, design);
  quantifier.quantifyPeptides();
  quantifier.quantifyProteins();

  const auto& quant = quantifier.getProteinResults();
  TEST_TRUE(quant.find("Prot") != quant.end());
  const auto& pd = quant.at("Prot");

  TEST_REAL_SIMILAR(pd.channel_level_abundances.at("fileA").at(1), 300.0);
  TEST_REAL_SIMILAR(pd.channel_level_abundances.at("fileB").at(1), 350.0);
  TEST_REAL_SIMILAR(pd.channel_level_abundances.at("fileA").at(2), 30.0);
  TEST_REAL_SIMILAR(pd.channel_level_abundances.at("fileB").at(2), 35.0);

  // the cells of a sample add up to that sample's value
  TEST_REAL_SIMILAR(pd.total_abundances.at(0), 650.0);
  TEST_REAL_SIMILAR(pd.total_abundances.at(1), 65.0);
  TEST_REAL_SIMILAR(pd.channel_level_abundances.at("fileA").at(1) +
                    pd.channel_level_abundances.at("fileB").at(1), pd.total_abundances.at(0));
  TEST_REAL_SIMILAR(pd.channel_level_abundances.at("fileA").at(2) +
                    pd.channel_level_abundances.at("fileB").at(2), pd.total_abundances.at(1));
}
END_SECTION

START_SECTION((const ProteinQuant& getProteinResults() fractionated: "top:N" is enforced per file, not per sample))
{
  // Same input, but with top:N = 3. Each sample has three peptides, so the protein IS
  // quantified at the sample level; no individual fraction file holds three peptides, so
  // every (file, channel) cell is dropped. This is intended - the cells are per-file values -
  // and is documented on the 'file_and_channel_level_output' parameter. It is also strictly
  // stricter than the sample-level rule, which is why it deserves pinning.
  ConsensusMap consensus;
  ExperimentalDesign design;
  make_fractionated_input({"fileA", "fileB"}, 2,
                          {{0, 1, "PEPTIDEK", 2, 100.0}, {0, 2, "PEPTIDEK", 2, 10.0},
                           {0, 1, "AAAAAK",   2, 200.0}, {0, 2, "AAAAAK",   2, 20.0},
                           {1, 1, "CCCCCK",   2, 300.0}, {1, 2, "CCCCCK",   2, 30.0},
                           {1, 1, "PEPTIDEK", 2,  50.0}, {1, 2, "PEPTIDEK", 2,  5.0}},
                          "Prot", consensus, design);

  PeptideAndProteinQuant quantifier;
  Param p = quantifier.getDefaults();
  p.setValue("top:N", 3);
  p.setValue("top:aggregate", "sum");
  p.setValue("top:include_all", "false");
  quantifier.setParameters(p);

  quantifier.readQuantData(consensus, design);
  quantifier.quantifyPeptides();
  quantifier.quantifyProteins();

  const auto& quant = quantifier.getProteinResults();
  TEST_TRUE(quant.find("Prot") != quant.end());
  const auto& pd = quant.at("Prot");
  TEST_EQUAL(pd.peptide_abundances.size(), 3);

  // quantified per sample (three peptides each) ...
  TEST_REAL_SIMILAR(pd.total_abundances.at(0), 650.0);
  TEST_REAL_SIMILAR(pd.total_abundances.at(1), 65.0);
  // ... but no fraction file holds three peptides, so no cell survives
  TEST_TRUE(pd.channel_level_abundances.empty());
}
END_SECTION

START_SECTION((void annotateQuantificationsToProteins(const ProteinQuant& protein_quants, ProteinIdentification& proteins, bool remove_unquantified) isobaric: one file, N channels))
{
  // This is the shape IsobaricAnalyzer produces and IsobaricWorkflow quantifies: ONE file
  // whose column headers differ only in label/channel. A sample is a (file, label) pair, so
  // here every (file, channel) cell IS exactly one sample - the two granularities are
  // degenerate and cannot disagree, which is why the per-file peptide selection and
  // aggregation documented on 'file_and_channel_level_output' are harmless for unfractionated
  // isobaric data (contrast the two fractionated sections above).
  //   fileA ch1: PEPTIDEK 100, AAAAAK 200, CCCCCK 300     fileA ch2: 10, 20, 30
  // This section also covers annotateQuantificationsToProteins(), the route by which
  // channel_level_abundances reaches mzTab/QPX - IsobaricWorkflow's only consumer of it.
  ConsensusMap consensus;
  ExperimentalDesign design;
  make_fractionated_input({"fileA"}, 2,
                          {{0, 1, "PEPTIDEK", 2, 100.0}, {0, 2, "PEPTIDEK", 2, 10.0},
                           {0, 1, "AAAAAK",   2, 200.0}, {0, 2, "AAAAAK",   2, 20.0},
                           {0, 1, "CCCCCK",   2, 300.0}, {0, 2, "CCCCCK",   2, 30.0}},
                          "Prot", consensus, design);

  TEST_EQUAL(design.getNumberOfSamples(), 2);
  TEST_EQUAL(design.getNumberOfFractions(), 1);
  TEST_EQUAL(design.getNumberOfLabels(), 2);

  ProteinIdentification proteins;
  ProteinHit prot_hit;
  prot_hit.setAccession("Prot");
  proteins.setHits({prot_hit});
  ProteinIdentification::ProteinGroup group;
  group.accessions = {"Prot"};
  group.probability = 1.0;
  proteins.getIndistinguishableProteins().push_back(group);

  PeptideAndProteinQuant quantifier;
  // stock defaults on purpose: top:N 3, "median", include_all false - the configuration whose
  // per-file peptide-count gate drops every cell of a fractionated design
  quantifier.setParameters(quantifier.getDefaults());

  quantifier.readQuantData(consensus, design);
  quantifier.quantifyPeptides();
  quantifier.quantifyProteins(proteins);

  const auto& quant = quantifier.getProteinResults();
  TEST_TRUE(quant.find("Prot") != quant.end());
  const auto& pd = quant.at("Prot");

  // each channel holds three peptides in this single file, so the gate passes and the cells
  // are present - median{100,200,300} = 200 and median{10,20,30} = 20
  TEST_REAL_SIMILAR(pd.channel_level_abundances.at("fileA").at(1), 200.0);
  TEST_REAL_SIMILAR(pd.channel_level_abundances.at("fileA").at(2), 20.0);
  // ... and they equal the sample-level values, because cell == sample here
  TEST_REAL_SIMILAR(pd.total_abundances.at(0), 200.0);
  TEST_REAL_SIMILAR(pd.total_abundances.at(1), 20.0);

  quantifier.annotateQuantificationsToProteins(quant, proteins, true);

  TEST_EQUAL(proteins.getIndistinguishableProteins().size(), 1);
  const auto& annotated = proteins.getIndistinguishableProteins()[0];

  TEST_EQUAL(annotated.getFloatDataArrays().size(), 4);
  TEST_EQUAL(annotated.getStringDataArrays().size(), 2);
  TEST_EQUAL(annotated.getIntegerDataArrays().size(), 2);

  const auto& sample_abundances = annotated.getFloatDataArrays()[0];
  const auto& cell_abundances = annotated.getFloatDataArrays()[3];
  const auto& cell_filenames = annotated.getStringDataArrays()[0];
  const auto& cell_channels = annotated.getIntegerDataArrays()[0];

  TEST_STRING_EQUAL(sample_abundances.getName(), "abundances");
  TEST_STRING_EQUAL(cell_abundances.getName(), "file_channel_level_abundance");
  TEST_STRING_EQUAL(cell_filenames.getName(), "file_channel_level_filename");
  TEST_STRING_EQUAL(cell_channels.getName(), "file_channel_level_channel");

  // one entry per (file, channel) of the design, in design order
  TEST_EQUAL(cell_abundances.size(), 2);
  TEST_EQUAL(cell_filenames.size(), 2);
  TEST_EQUAL(cell_channels.size(), 2);

  TEST_STRING_EQUAL(cell_filenames[0], "fileA");
  TEST_STRING_EQUAL(cell_filenames[1], "fileA");
  TEST_EQUAL(cell_channels[0], 1);
  TEST_EQUAL(cell_channels[1], 2);

  // the annotated cells carry the aggregated values, and coincide with the sample-level
  // array element-wise - the property that makes this output well-defined for isobaric data
  TEST_EQUAL(sample_abundances.size(), 2);
  TEST_REAL_SIMILAR(cell_abundances[0], 200.0);
  TEST_REAL_SIMILAR(cell_abundances[1], 20.0);
  TEST_REAL_SIMILAR(cell_abundances[0], sample_abundances[0]);
  TEST_REAL_SIMILAR(cell_abundances[1], sample_abundances[1]);
}
END_SECTION

START_SECTION((const ProteinQuant& getProteinResults() file+channel level must sum all peptidoforms))
{
  // Two peptidoforms of the same unmodified sequence, in one file/channel.
  // transferPeptideDataToProteins_() merges them into a single peptide
  // (10 + 90 = 100), so the file+channel level cell must be 100 as well - not
  // the value of whichever peptidoform happens to come first in pep_quant_.
  ConsensusMap consensus;
  ExperimentalDesign design;
  make_detailed_input({"fileA"},
                      {{0, "PEPTIDEMK", 2, 10.0},
                       {0, "PEPTIDEM(Oxidation)K", 2, 90.0}},
                      "Prot", consensus, design);

  PeptideAndProteinQuant quantifier;
  Param p = quantifier.getDefaults();
  p.setValue("top:N", 0);
  p.setValue("top:aggregate", "sum");
  p.setValue("top:include_all", "true");
  quantifier.setParameters(p);

  quantifier.readQuantData(consensus, design);
  quantifier.quantifyPeptides();
  quantifier.quantifyProteins();

  // peptide level keeps the peptidoforms apart
  TEST_EQUAL(quantifier.getPeptideResults().size(), 2);

  const auto& quant = quantifier.getProteinResults();
  TEST_TRUE(quant.find("Prot") != quant.end());
  const auto& pd = quant.at("Prot");

  // protein level merges them into one peptide
  TEST_EQUAL(pd.peptide_abundances.size(), 1);
  TEST_REAL_SIMILAR(pd.total_abundances.at(0), 100.0);

  // ... and the file+channel level must agree
  TEST_TRUE(pd.channel_level_abundances.find("fileA") != pd.channel_level_abundances.end());
  TEST_REAL_SIMILAR(pd.channel_level_abundances.at("fileA").at(1), 100.0);
}
END_SECTION

START_SECTION((const ProteinQuant& getProteinResults() file+channel level "top:N" must count peptides, not charge states))
{
  // "Prot" has three distinct peptides overall, so it passes the protein-level
  // "at least N peptides" gate. In fileA only ONE of them is observed, in three
  // charge states; in fileB all three are observed once each.
  // With top:N = 3 the fileA cell must be dropped (one peptide < 3), exactly as
  // the sample-level result for fileA is dropped. Charge states must not be
  // counted as separate peptides.
  ConsensusMap consensus;
  ExperimentalDesign design;
  make_detailed_input({"fileA", "fileB"},
                      {{0, "PEPTIDEK", 2, 100.0},
                       {0, "PEPTIDEK", 3, 100.0},
                       {0, "PEPTIDEK", 4, 100.0},
                       {1, "PEPTIDEK", 2, 10.0},
                       {1, "AAAAAK", 2, 10.0},
                       {1, "CCCCCK", 2, 10.0}},
                      "Prot", consensus, design);

  PeptideAndProteinQuant quantifier;
  Param p = quantifier.getDefaults();
  p.setValue("top:N", 3);
  p.setValue("top:aggregate", "sum");
  p.setValue("top:include_all", "false");
  quantifier.setParameters(p);

  quantifier.readQuantData(consensus, design);
  quantifier.quantifyPeptides();
  quantifier.quantifyProteins();

  const auto& quant = quantifier.getProteinResults();
  TEST_TRUE(quant.find("Prot") != quant.end());
  const auto& pd = quant.at("Prot");
  TEST_EQUAL(pd.peptide_abundances.size(), 3); // clears the protein-level gate

  // sample level: fileA (sample 0) has only one peptide -> no abundance
  TEST_TRUE(pd.total_abundances.find(0) == pd.total_abundances.end());
  TEST_REAL_SIMILAR(pd.total_abundances.at(1), 30.0);

  // file+channel level must make the same decision
  TEST_TRUE(pd.channel_level_abundances.find("fileA") == pd.channel_level_abundances.end());
  TEST_TRUE(pd.channel_level_abundances.find("fileB") != pd.channel_level_abundances.end());
  TEST_REAL_SIMILAR(pd.channel_level_abundances.at("fileB").at(1), 30.0);
}
END_SECTION

START_SECTION((const ProteinQuant& getProteinResults() file+channel level must not keep un-aggregated sums))
{
  // fileA sees only two of the protein's three peptides, so with top:N = 3 its
  // cell cannot be quantified. It must then be ABSENT - it must not fall back to
  // a raw sum over all peptides collected while transferring peptide data to
  // proteins (which would be 5 + 7 = 12 here, a value that never passed top:N).
  ConsensusMap consensus;
  ExperimentalDesign design;
  make_detailed_input({"fileA", "fileB"},
                      {{0, "AAAAAK", 2, 5.0},
                       {0, "CCCCCK", 2, 7.0},
                       {1, "AAAAAK", 2, 1.0},
                       {1, "CCCCCK", 2, 1.0},
                       {1, "EEEEEK", 2, 1.0}},
                      "Prot", consensus, design);

  PeptideAndProteinQuant quantifier;
  Param p = quantifier.getDefaults();
  p.setValue("top:N", 3);
  p.setValue("top:aggregate", "sum");
  p.setValue("top:include_all", "false");
  quantifier.setParameters(p);

  quantifier.readQuantData(consensus, design);
  quantifier.quantifyPeptides();
  quantifier.quantifyProteins();

  const auto& quant = quantifier.getProteinResults();
  TEST_TRUE(quant.find("Prot") != quant.end());
  const auto& pd = quant.at("Prot");
  TEST_EQUAL(pd.peptide_abundances.size(), 3);

  // sample level drops fileA (sample 0)
  TEST_TRUE(pd.total_abundances.find(0) == pd.total_abundances.end());
  TEST_REAL_SIMILAR(pd.total_abundances.at(1), 3.0);

  // file+channel level must drop it too, instead of reporting the raw 12.0
  TEST_TRUE(pd.channel_level_abundances.find("fileA") == pd.channel_level_abundances.end());
  TEST_REAL_SIMILAR(pd.channel_level_abundances.at("fileB").at(1), 3.0);
}
END_SECTION

START_SECTION((const ProteinQuant& getProteinResults() file+channel level "top:N" truncation must rank peptides))
{
  // The sections above use "sum", which is insensitive to how the values are grouped as
  // long as nothing is dropped. This one uses "mean" together with an actually truncating
  // top:N, so a wrong grouping changes the result:
  //   peptide AAAAAMK is split over a peptidoform AND a charge state (100 + 90 = 190),
  //   CCCCCK = 80, EEEEEK = 70.
  // Ranking peptides gives mean(190, 80) = 135; ranking raw peptidoform/charge cells
  // instead would rank 100 and 90 into the top 2 and give mean(100, 90) = 95.
  ConsensusMap consensus;
  ExperimentalDesign design;
  make_detailed_input({"fileA"},
                      {{0, "AAAAAMK", 2, 100.0},
                       {0, "AAAAAM(Oxidation)K", 3, 90.0},
                       {0, "CCCCCK", 2, 80.0},
                       {0, "EEEEEK", 2, 70.0}},
                      "Prot", consensus, design);

  PeptideAndProteinQuant quantifier;
  Param p = quantifier.getDefaults();
  p.setValue("top:N", 2);
  p.setValue("top:aggregate", "mean");
  p.setValue("top:include_all", "false");
  quantifier.setParameters(p);

  quantifier.readQuantData(consensus, design);
  quantifier.quantifyPeptides();
  quantifier.quantifyProteins();

  const auto& quant = quantifier.getProteinResults();
  TEST_TRUE(quant.find("Prot") != quant.end());
  const auto& pd = quant.at("Prot");
  TEST_EQUAL(pd.peptide_abundances.size(), 3);

  TEST_REAL_SIMILAR(pd.total_abundances.at(0), 135.0);
  TEST_REAL_SIMILAR(pd.channel_level_abundances.at("fileA").at(1), 135.0);
}
END_SECTION

START_SECTION((const ProteinQuant& getProteinResults() iBAQ must also normalize the file+channel level))
{
  // iBAQ divides the protein abundance by the number of theoretically observable
  // peptides. That normalization has to reach the file+channel level as well,
  // otherwise the detailed output is off by the (protein-dependent) divisor.
  ConsensusMap consensus;
  ExperimentalDesign design;
  make_detailed_input({"fileA"},
                      {{0, "AAAAAK", 2, 100.0},
                       {0, "CCCCCK", 2, 200.0}},
                      "Prot", consensus, design);

  ProteinHit prot_hit;
  prot_hit.setAccession("Prot");
  prot_hit.setSequence("AAAAAKCCCCCKEEEEEK"); // 3 fully tryptic peptides
  ProteinIdentification proteins;
  proteins.setHits({prot_hit});

  PeptideAndProteinQuant quantifier;
  Param p = quantifier.getDefaults();
  p.setValue("method", "iBAQ");
  quantifier.setParameters(p);

  quantifier.readQuantData(consensus, design);
  quantifier.quantifyPeptides();
  quantifier.quantifyProteins(proteins);

  const auto& quant = quantifier.getProteinResults();
  TEST_TRUE(quant.find("Prot") != quant.end());
  const auto& pd = quant.at("Prot");

  // (100 + 200) / 3 theoretical peptides
  TEST_REAL_SIMILAR(pd.total_abundances.at(0), 100.0);
  TEST_REAL_SIMILAR(pd.channel_level_abundances.at("fileA").at(1), 100.0);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

