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
#include <OpenMS/METADATA/ExperimentalDesign.h>
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
/////////////////////////////////////////////////////////////
END_TEST

