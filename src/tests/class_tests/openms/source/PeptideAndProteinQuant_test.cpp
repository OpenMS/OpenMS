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

#include <algorithm>
#include <cmath>


using namespace OpenMS;
using namespace std;


// Builds a one-file isobaric ConsensusMap: one consensus feature per peptide, one sub-feature per
// channel. `intensities[p][c]` is the reporter intensity of peptide p in channel c; a 0 stands for a
// reporter that was not detected, which is exactly what IsobaricChannelExtractor writes out.
namespace
{
  PeptideAndProteinQuant makeNormalizedQuantifier(const std::vector<std::vector<double>>& intensities,
                                                 bool best_charge_and_fraction = false)
  {
    const Size n_channels = intensities.front().size();

    ConsensusMap consensus;
    consensus.setExperimentType("labeled_MS2");
    for (Size c = 0; c < n_channels; ++c)
    {
      ConsensusMap::ColumnHeader h;
      h.filename = "A.mzML";
      h.label = "ch" + std::to_string(c + 1);
      h.setMetaValue("channel_id", (int)c); // label = channel_id + 1
      consensus.getColumnHeaders()[c] = h;
    }

    for (Size p = 0; p < intensities.size(); ++p)
    {
      ConsensusFeature cf;
      for (Size c = 0; c < n_channels; ++c)
      {
        Peak2D peak;
        peak.setIntensity(intensities[p][c]);
        peak.setRT(1.0 + p);
        peak.setMZ(400.0 + p);
        cf.insert(c, peak, c);
      }
      PeptideIdentification pid;
      PeptideHit hit;
      hit.setSequence(AASequence::fromString(std::string("PEPTIDE") + char('A' + p)));
      hit.setCharge(2);
      hit.setScore(1.0);
      pid.setHits({hit});
      cf.setPeptideIdentifications({pid});
      consensus.push_back(cf);
    }

    ExperimentalDesign::MSFileSection fs;
    ExperimentalDesign::SampleSection ss;
    for (Size c = 0; c < n_channels; ++c)
    {
      ExperimentalDesign::MSFileSectionEntry r;
      r.path = "/tmp/A.mzML";
      r.label = c + 1;
      r.fraction = 1;
      r.fraction_group = 1;
      r.sample = c;
      r.sample_name = "S" + std::to_string(c + 1);
      fs.push_back(r);
      ss.addSample("S" + std::to_string(c + 1));
    }
    ExperimentalDesign design(fs, ss);

    PeptideAndProteinQuant quantifier;
    Param p = quantifier.getParameters();
    p.setValue("consensus:normalize", "true");
    if (best_charge_and_fraction) { p.setValue("best_charge_and_fraction", "true"); }
    quantifier.setParameters(p);
    quantifier.readQuantData(consensus, design);
    quantifier.quantifyPeptides();
    return quantifier;
  }

  /// Abundances of one peptide across all assays, sorted so the test does not depend on assay keys.
  std::vector<double> sortedAbundances(PeptideAndProteinQuant& q, const std::string& sequence)
  {
    const auto& pep_quant = q.getPeptideResults();
    const auto it = pep_quant.find(AASequence::fromString(sequence));
    std::vector<double> out;
    if (it == pep_quant.end()) { return out; }
    for (const auto& [fraction_group, label_abundances] : it->second.fraction_group_abundances)
    {
      (void)fraction_group;
      for (const auto& label_abundance : label_abundances) { out.push_back(label_abundance.second); }
    }
    std::sort(out.begin(), out.end());
    return out;
  }

  ExperimentalDesign::MSFileSectionEntry designRow(const std::string& path,
                                                    UInt label,
                                                    UInt fraction_group,
                                                    Size sample)
  {
    ExperimentalDesign::MSFileSectionEntry row;
    row.path = path;
    row.label = label;
    row.fraction = 1;
    row.fraction_group = fraction_group;
    row.sample = sample;
    row.sample_name = "S" + std::to_string(sample + 1);
    return row;
  }

  // Use the setters deliberately: these malformed designs must reach PeptideAndProteinQuant's
  // lookup validation instead of being rejected by ExperimentalDesign's public constructor.
  ExperimentalDesign uncheckedDesign(const ExperimentalDesign::MSFileSection& rows,
                                     Size sample_count)
  {
    ExperimentalDesign::SampleSection samples;
    for (Size i = 0; i < sample_count; ++i)
    {
      samples.addSample("S" + std::to_string(i + 1));
    }
    ExperimentalDesign design;
    design.setMSFileSection(rows);
    design.setSampleSection(samples);
    return design;
  }
}

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
  // ...and end up in two distinct assays, each with a single spectral-count abundance.
  TEST_EQUAL(seq_it->second.fraction_group_abundances.size(), 2);
  // fromIdentifications currently numbers inferred fraction groups from zero.
  TEST_REAL_SIMILAR(seq_it->second.fraction_group_abundances.at(0).at(1), 1);
  TEST_REAL_SIMILAR(seq_it->second.fraction_group_abundances.at(1).at(1), 1);
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

START_SECTION(([EXTRA] experimental-design lookups validate coordinates used by assay quantification))
{
  FeatureMap features;

  // A design cell cannot refer beyond the sample section.
  {
    ExperimentalDesign design = uncheckedDesign({designRow("A.mzML", 1, 1, 1)}, 1);
    PeptideAndProteinQuant quantifier;
    TEST_EXCEPTION(Exception::MissingInformation, quantifier.readQuantData(features, design))
  }

  // One (run, label) cell cannot carry conflicting sample coordinates.
  {
    ExperimentalDesign design = uncheckedDesign({designRow("A.mzML", 1, 1, 0),
                                                 designRow("A.mzML", 1, 1, 1)}, 2);
    PeptideAndProteinQuant quantifier;
    TEST_EXCEPTION(Exception::MissingInformation, quantifier.readQuantData(features, design))
  }

  // Distinct labels of one run cannot place that run in several fraction groups.
  {
    ExperimentalDesign design = uncheckedDesign({designRow("A.mzML", 1, 1, 0),
                                                 designRow("A.mzML", 2, 2, 1)}, 2);
    PeptideAndProteinQuant quantifier;
    TEST_EXCEPTION(Exception::MissingInformation, quantifier.readQuantData(features, design))
  }

  // One assay cannot belong to several SampleSection entries: mzTab still maps each assay to one
  // study variable even though Sample is no longer an abundance grain.
  {
    ExperimentalDesign design = uncheckedDesign({designRow("A.mzML", 1, 1, 0),
                                                 designRow("B.mzML", 1, 1, 1)}, 2);
    PeptideAndProteinQuant quantifier;
    TEST_EXCEPTION(Exception::MissingInformation, quantifier.readQuantData(features, design))
  }
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

START_SECTION(([EXTRA] fraction-group/label quantities remain distinct and omit absent design runs))
{
  // One sample is reused in two fraction groups. They remain two independently aggregated assays;
  // a third design run not represented by a ConsensusMap column header must not invent a zero.
  ConsensusMap consensus;
  consensus.setExperimentType("label-free");
  for (Size i = 0; i < 2; ++i)
  {
    ConsensusMap::ColumnHeader header;
    header.filename = i == 0 ? "/data/technical_a.mzML" : "/data/technical_b.mzML";
    header.label = "label-free";
    consensus.getColumnHeaders()[i] = header;
  }

  PeptideIdentification peptide_id;
  peptide_id.setIdentifier("run");
  PeptideHit peptide_hit;
  peptide_hit.setSequence(AASequence::fromString("PEPTIDEK"));
  peptide_hit.setCharge(2);
  PeptideEvidence evidence;
  evidence.setProteinAccession("PROT_A");
  peptide_hit.setPeptideEvidences({evidence});
  peptide_id.setHits({peptide_hit});

  ConsensusFeature feature;
  for (Size i = 0; i < 2; ++i)
  {
    BaseFeature handle;
    handle.setIntensity(i == 0 ? 30.0 : 40.0);
    feature.insert(i, handle);
  }
  feature.setPeptideIdentifications({peptide_id});
  consensus.push_back(feature);

  ProteinIdentification proteins;
  proteins.setIdentifier("run");
  ProteinHit protein_hit;
  protein_hit.setAccession("PROT_A");
  proteins.setHits({protein_hit});
  ProteinIdentification::ProteinGroup group;
  group.accessions = {"PROT_A"};
  proteins.insertIndistinguishableProteins(group);
  consensus.setProteinIdentifications({proteins});

  ExperimentalDesign::MSFileSection files;
  for (Size i = 0; i < 3; ++i)
  {
    ExperimentalDesign::MSFileSectionEntry entry;
    entry.path = i == 0 ? "/data/technical_a.mzML"
                        : (i == 1 ? "/data/technical_b.mzML" : "/data/not_in_map.mzML");
    entry.label = 1;
    entry.fraction = 1;
    entry.fraction_group = static_cast<UInt>(i) + 1;
    entry.sample = 0;
    entry.sample_name = "S1";
    files.push_back(entry);
  }
  ExperimentalDesign::SampleSection samples;
  samples.addSample("S1");
  ExperimentalDesign design(files, samples);

  PeptideAndProteinQuant quantifier;
  Param p;
  p.setValue("top:include_all", "true");
  quantifier.setParameters(p);
  quantifier.readQuantData(consensus, design);
  quantifier.quantifyPeptides();

  const auto& peptide = quantifier.getPeptideResults().at(AASequence::fromString("PEPTIDEK"));
  TEST_EQUAL(peptide.fraction_group_abundances.size(), 2)
  TEST_REAL_SIMILAR(peptide.fraction_group_abundances.at(1).at(1), 30.0)
  TEST_REAL_SIMILAR(peptide.fraction_group_abundances.at(2).at(1), 40.0)

  quantifier.quantifyProteins(proteins);
  const auto& protein = quantifier.getProteinResults().at("PROT_A");
  TEST_EQUAL(protein.fraction_group_abundances.size(), 2)
  TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(1), 30.0)
  TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(2).at(1), 40.0)
  TEST_EQUAL(protein.fraction_group_abundances.count(3), 0)

  quantifier.annotateQuantificationsToProteins(quantifier.getProteinResults(), proteins, false);
  const auto& annotated = proteins.getIndistinguishableProteins()[0];
  const auto& float_arrays = annotated.getFloatDataArrays();
  const auto& integer_arrays = annotated.getIntegerDataArrays();
  TEST_STRING_EQUAL(float_arrays[3].getName(), "fraction_group_level_abundance")
  TEST_STRING_EQUAL(integer_arrays[2].getName(), "fraction_group_level_fraction_group")
  TEST_STRING_EQUAL(integer_arrays[3].getName(), "fraction_group_level_label")
  TEST_EQUAL(float_arrays[3].size(), 2)
  TEST_REAL_SIMILAR(float_arrays[3][0], 30.0)
  TEST_REAL_SIMILAR(float_arrays[3][1], 40.0)
  TEST_EQUAL(integer_arrays[2][0], 1)
  TEST_EQUAL(integer_arrays[2][1], 2)
  TEST_EQUAL(integer_arrays[3][0], 1)
  TEST_EQUAL(integer_arrays[3][1], 1)
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
  TEST_EQUAL(pep_data.fraction_group_abundances.size(), 1);
  TEST_REAL_SIMILAR(pep_data.fraction_group_abundances.at(1).at(1), 3333);
  TEST_EQUAL(pep_data.accessions.size(), 1);
  TEST_EQUAL(pep_data.psm_count, 2);
  pep_data = pep_quant[AASequence::fromString("CCCCC")];
  TEST_EQUAL(pep_data.abundances.size(), 1); // one fraction
  TEST_EQUAL(pep_data.abundances[1].size(), 1); // one file
  auto& map_file_to_charges = *pep_data.abundances[1].begin();
  TEST_EQUAL(map_file_to_charges.second.size(), 2); // two charges

  
  TEST_EQUAL(pep_data.fraction_group_abundances.size(), 1);
  TEST_REAL_SIMILAR(pep_data.fraction_group_abundances.at(1).at(1), 7777);
  TEST_EQUAL(pep_data.accessions.size(), 1);
  TEST_EQUAL(pep_data.psm_count, 2);
  pep_data = pep_quant[AASequence::fromString("EEEEE")];
  TEST_EQUAL(pep_data.abundances.size(), 0); // it is the second best hit, so it will not be counted
  TEST_TRUE(pep_data.fraction_group_abundances.empty());
  TEST_EQUAL(pep_data.accessions.size(), 1);
  TEST_EQUAL(pep_data.psm_count, 1);
  pep_data = pep_quant[AASequence::fromString("GGGGG")];
  TEST_EQUAL(pep_data.abundances.size(), 1); // one fraction
  TEST_EQUAL(pep_data.abundances[1].size(), 1); // one file
  TEST_EQUAL((*pep_data.abundances[1].begin()).second.size(), 1); // two charges  

  TEST_EQUAL(pep_data.fraction_group_abundances.size(), 1);
  TEST_REAL_SIMILAR(pep_data.fraction_group_abundances.at(1).at(1), 7777);
  TEST_EQUAL(pep_data.accessions.size(), 2);
  TEST_EQUAL(pep_data.psm_count, 1);

  pep_quant = quantifier_consensus.getPeptideResults();
  TEST_EQUAL(pep_quant.size(), 4);
  pep_data = pep_quant[AASequence::fromString("AAAK")];
  TEST_EQUAL(pep_data.abundances.size(), 1); // one fraction
  TEST_EQUAL(pep_data.abundances[1].size(), 2); // two files
  TEST_EQUAL((*pep_data.abundances[1].begin()).second.size(), 1); // one charge

  TEST_EQUAL(pep_data.fraction_group_abundances.size(), 2);
  TEST_REAL_SIMILAR(pep_data.fraction_group_abundances.at(1).at(1), 1000);
  TEST_REAL_SIMILAR(pep_data.fraction_group_abundances.at(3).at(1), 1000);
  TEST_EQUAL(pep_data.accessions.size(), 1);
  TEST_EQUAL(pep_data.psm_count, 1);
  pep_data = pep_quant[AASequence::fromString("CCCK")];
  TEST_EQUAL(pep_data.abundances.size(), 1); // one fraction
  TEST_EQUAL(pep_data.abundances[1].size(), 2); // two files
  TEST_EQUAL((*pep_data.abundances[1].begin()).second.size(), 1); // one charge

  TEST_EQUAL(pep_data.fraction_group_abundances.size(), 2);
  TEST_REAL_SIMILAR(pep_data.fraction_group_abundances.at(1).at(1), 200);
  TEST_REAL_SIMILAR(pep_data.fraction_group_abundances.at(2).at(1), 200);
  TEST_EQUAL(pep_data.accessions.size(), 1);
  TEST_EQUAL(pep_data.psm_count, 1);
  pep_data = pep_quant[AASequence::fromString("EEEK")];
  TEST_EQUAL(pep_data.abundances.size(), 1); // one fraction
  TEST_EQUAL(pep_data.abundances[1].size(), 3); // three files
  TEST_EQUAL((*pep_data.abundances[1].begin()).second.size(), 1); // one charge

  TEST_EQUAL(pep_data.fraction_group_abundances.size(), 3);
  TEST_REAL_SIMILAR(pep_data.fraction_group_abundances.at(1).at(1), 30);
  TEST_REAL_SIMILAR(pep_data.fraction_group_abundances.at(2).at(1), 30);
  TEST_REAL_SIMILAR(pep_data.fraction_group_abundances.at(3).at(1), 30);
  TEST_EQUAL(pep_data.accessions.size(), 1);
  TEST_EQUAL(pep_data.psm_count, 1);
  pep_data = pep_quant[AASequence::fromString("GGG")];
  TEST_EQUAL(pep_data.abundances.size(), 1); // one fraction
  TEST_EQUAL(pep_data.abundances[1].size(), 2); // two files
  TEST_EQUAL((*pep_data.abundances[1].begin()).second.size(), 1); // one charge

  TEST_EQUAL(pep_data.fraction_group_abundances.size(), 2);
  TEST_REAL_SIMILAR(pep_data.fraction_group_abundances.at(1).at(1), 4);
  TEST_REAL_SIMILAR(pep_data.fraction_group_abundances.at(2).at(1), 4);
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
  TEST_EQUAL(prot_data.peptide_fraction_group_abundances.size(), 3);
  TEST_REAL_SIMILAR(prot_data.fraction_group_abundances.at(1).at(1), 4711);
  TEST_EQUAL(prot_data.psm_count, 6);
  prot_data = prot_quant["Protein1"];
  TEST_EQUAL(prot_data.peptide_fraction_group_abundances.size(), 1);
  TEST_REAL_SIMILAR(prot_data.fraction_group_abundances.at(1).at(1), 8888);
  TEST_EQUAL(prot_data.psm_count, 2);

  prot_quant = quantifier_consensus.getProteinResults();
  TEST_EQUAL(prot_quant.size(), 1);
  prot_data = prot_quant["Protein"];
  TEST_EQUAL(prot_data.peptide_fraction_group_abundances.size(), 4);
  TEST_EQUAL(prot_data.fraction_group_abundances.size(), 3);
  TEST_REAL_SIMILAR(prot_data.fraction_group_abundances.at(1).at(1), 200);
  TEST_REAL_SIMILAR(prot_data.fraction_group_abundances.at(2).at(1), 30);
  TEST_REAL_SIMILAR(prot_data.fraction_group_abundances.at(3).at(1), 515);
  TEST_EQUAL(prot_data.psm_count, 4);
}
END_SECTION

START_SECTION(([PeptideAndProteinQuant::PeptideData] PeptideData()))
{
  PeptideAndProteinQuant::PeptideData data;
  TEST_TRUE(data.abundances.empty());
  TEST_TRUE(data.fraction_group_abundances.empty());
  TEST_TRUE(data.accessions.empty());
  TEST_EQUAL(data.psm_count, 0);
}
END_SECTION

START_SECTION(([PeptideAndProteinQuant::ProteinData] ProteinData()))
{
  PeptideAndProteinQuant::ProteinData data;
  TEST_TRUE(data.peptide_fraction_group_abundances.empty());
  TEST_TRUE(data.fraction_group_abundances.empty());
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
  TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(1), 4711);

  parameters.setValue("top:aggregate", "mean");
  quantifier.setParameters(parameters);
  quantifier.readQuantData(f, ed);
  quantifier.quantifyPeptides();
  quantifier.quantifyProteins();
  quant = quantifier.getProteinResults();
  protein = quant["Protein0"];
  TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(1), 5273.666666);

  parameters.setValue("top:aggregate", "weighted_mean");
  quantifier.setParameters(parameters);
  quantifier.readQuantData(f, ed);
  quantifier.quantifyPeptides();
  quantifier.quantifyProteins();
  quant = quantifier.getProteinResults();
  protein = quant["Protein0"];
  TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(1), 5927.82624360028);

  parameters.setValue("top:aggregate", "sum");
  quantifier.setParameters(parameters);
  quantifier.readQuantData(f, ed);
  quantifier.quantifyPeptides();
  quantifier.quantifyProteins();
  quant = quantifier.getProteinResults();
  protein = quant["Protein0"];
  TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(1), 15821);
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
  TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(1), 308.5);
  TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(2).at(1), 58.5);
  TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(3).at(1), 257.5);
}
END_SECTION


START_SECTION(([EXTRA] normalization: a channel with mostly undetected reporters is scaled on the data it has))
{
  // channel 3 has no reporter for 3 of 5 peptides, so its median over ALL values would be 0 and the
  // scale factor overall_median/0 = inf - which used to turn every abundance in the channel into inf
  // and, once protein quantities were written to consensusXML, produced a file that could not be read
  // back at all. The median is taken over the detected values only.
  //   ch1: 10 20 30 40  50   -> median  30 -> factor 60/30  = 2
  //   ch2: 20 40 60 80 100   -> median  60 -> factor 60/60  = 1
  //   ch3:  -  -  - 90 120   -> median 105 -> factor 60/105
  //   overall median of (30, 60, 105) = 60
  PeptideAndProteinQuant q = makeNormalizedQuantifier({{10, 20, 0},
                                                       {20, 40, 0},
                                                       {30, 60, 0},
                                                       {40, 80, 90},
                                                       {50, 100, 120}});

  const std::vector<double> pd = sortedAbundances(q, "PEPTIDED"); // 40, 80, 90
  TEST_EQUAL(pd.size(), 3)
  for (double v : pd) { TEST_TRUE(std::isfinite(v)) }
  TEST_REAL_SIMILAR(pd[0], 90.0 * 60.0 / 105.0) // ~51.43
  TEST_REAL_SIMILAR(pd[1], 80.0)                // 40 * 2
  TEST_REAL_SIMILAR(pd[2], 80.0)                // 80 * 1

  const std::vector<double> pa = sortedAbundances(q, "PEPTIDEA"); // 10, 20, undetected
  TEST_EQUAL(pa.size(), 3)
  TEST_REAL_SIMILAR(pa[0], 0.0)  // stays undetected
  TEST_REAL_SIMILAR(pa[1], 20.0) // 10 * 2
  TEST_REAL_SIMILAR(pa[2], 20.0) // 20 * 1
}
END_SECTION

START_SECTION(([EXTRA] normalization scales technical-replicate assays independently))
{
  // Both assays belong to one SampleSection entry. Their peptide medians are 20 and 200, so the
  // assay-level reference is 110 and the factors are 5.5 and 0.55. Sample-level normalization
  // would see only one sample and leave the four measurements unchanged.
  ConsensusMap consensus;
  consensus.setExperimentType("label-free");
  for (Size map_index = 0; map_index < 2; ++map_index)
  {
    ConsensusMap::ColumnHeader header;
    header.filename = map_index == 0 ? "A.mzML" : "B.mzML";
    consensus.getColumnHeaders()[map_index] = header;
  }

  const std::vector<std::vector<double>> intensities{{10.0, 100.0}, {30.0, 300.0}};
  for (Size peptide_index = 0; peptide_index < intensities.size(); ++peptide_index)
  {
    ConsensusFeature feature;
    for (Size map_index = 0; map_index < 2; ++map_index)
    {
      Peak2D peak;
      peak.setIntensity(intensities[peptide_index][map_index]);
      feature.insert(map_index, peak, map_index);
    }
    PeptideHit hit;
    hit.setSequence(AASequence::fromString(peptide_index == 0 ? "PEPTIDEA" : "PEPTIDEB"));
    hit.setCharge(2);
    hit.setScore(1.0);
    PeptideIdentification identification;
    identification.setHits({hit});
    feature.setPeptideIdentifications({identification});
    consensus.push_back(feature);
  }

  ExperimentalDesign design = uncheckedDesign({designRow("A.mzML", 1, 1, 0),
                                                designRow("B.mzML", 1, 2, 0)}, 1);
  PeptideAndProteinQuant quantifier;
  Param parameters = quantifier.getParameters();
  parameters.setValue("consensus:normalize", "true");
  quantifier.setParameters(parameters);
  quantifier.readQuantData(consensus, design);
  quantifier.quantifyPeptides();

  const auto& peptides = quantifier.getPeptideResults();
  const auto& peptide_a = peptides.at(AASequence::fromString("PEPTIDEA")).fraction_group_abundances;
  const auto& peptide_b = peptides.at(AASequence::fromString("PEPTIDEB")).fraction_group_abundances;
  TEST_REAL_SIMILAR(peptide_a.at(1).at(1), 55.0)
  TEST_REAL_SIMILAR(peptide_a.at(2).at(1), 55.0)
  TEST_REAL_SIMILAR(peptide_b.at(1).at(1), 165.0)
  TEST_REAL_SIMILAR(peptide_b.at(2).at(1), 165.0)
}
END_SECTION

START_SECTION(([EXTRA] normalization: a channel without any detected reporter does not skew the others))
{
  // channel 3 is empty throughout. It must be left alone AND kept out of the reference, otherwise its
  // median of 0 drags the overall median down and mis-scales every other channel - silently, with no
  // infinity to give it away.
  //   ch1: 10 20 30 40  50 -> median 30 |  ch2: 20 40 60 80 100 -> median 60
  //   overall median of (30, 60) = 45   ->  factors 1.5 and 0.75
  //   (were the empty channel counted, the reference would be median(30, 60, 0) = 30 instead)
  PeptideAndProteinQuant q = makeNormalizedQuantifier({{10, 20, 0},
                                                       {20, 40, 0},
                                                       {30, 60, 0},
                                                       {40, 80, 0},
                                                       {50, 100, 0}});

  const std::vector<double> pd = sortedAbundances(q, "PEPTIDED"); // 40, 80, empty
  TEST_EQUAL(pd.size(), 3)
  for (double v : pd) { TEST_TRUE(std::isfinite(v)) }
  TEST_REAL_SIMILAR(pd[0], 0.0)
  TEST_REAL_SIMILAR(pd[1], 60.0) // 40 * 1.5
  TEST_REAL_SIMILAR(pd[2], 60.0) // 80 * 0.75

  // the empty channel keeps its zeros - it is neither scaled by a made-up factor nor by 0
  const std::vector<double> pe = sortedAbundances(q, "PEPTIDEE"); // 50, 100, empty
  TEST_REAL_SIMILAR(pe[0], 0.0)
  TEST_REAL_SIMILAR(pe[1], 75.0)  // 50 * 1.5
  TEST_REAL_SIMILAR(pe[2], 75.0)  // 100 * 0.75
}
END_SECTION

START_SECTION(([EXTRA] normalization: best charge retains every assay before scaling))
{
  // Every reporter channel below is the same charge state. Selecting that charge must retain all
  // four assays, so all four take part in median normalization instead of only the largest cell of
  // each peptide being visible. Per-assay medians are 55, 110, 45, 55; their median is 55.
  PeptideAndProteinQuant q = makeNormalizedQuantifier({{100, 20, 30, 40},
                                                       {10, 200, 60, 70}},
                                                      /*best_charge_and_fraction=*/true);

  const auto& pep_quant = q.getPeptideResults();
  const auto& p1_groups = pep_quant.at(AASequence::fromString("PEPTIDEA")).fraction_group_abundances;
  const auto& p2_groups = pep_quant.at(AASequence::fromString("PEPTIDEB")).fraction_group_abundances;
  TEST_EQUAL(p1_groups.size(), 1)
  TEST_EQUAL(p2_groups.size(), 1)
  ABORT_IF(p1_groups.empty() || p2_groups.empty())
  const auto& p1 = p1_groups.begin()->second;
  const auto& p2 = p2_groups.begin()->second;
  TEST_EQUAL(p1.size(), 4);
  TEST_EQUAL(p2.size(), 4);
  TEST_REAL_SIMILAR(p1.at(1), 100.0);
  TEST_REAL_SIMILAR(p1.at(2), 10.0);
  TEST_REAL_SIMILAR(p1.at(3), 110.0 / 3.0);
  TEST_REAL_SIMILAR(p1.at(4), 40.0);
  TEST_REAL_SIMILAR(p2.at(1), 10.0);
  TEST_REAL_SIMILAR(p2.at(2), 100.0);
  TEST_REAL_SIMILAR(p2.at(3), 220.0 / 3.0);
  TEST_REAL_SIMILAR(p2.at(4), 70.0);

  Size positive = 0, total = 0;
  for (const auto& pq : pep_quant)
  {
    for (const auto& fa : pq.second.abundances)          // fractions
      for (const auto& fna : fa.second)                  // filenames
        for (const auto& ca : fna.second)                // charges
          for (const auto& cha : ca.second)              // channels
          {
            ++total;
            TEST_TRUE(std::isfinite(cha.second))
            if (cha.second > 0.0) { ++positive; }
          }
  }
  TEST_EQUAL(total, 8)     // 2 peptides x 4 channels, all originally positive
  TEST_EQUAL(positive, 8)
}
END_SECTION

/////////////////////////////////////////////////////////////
// File+channel level ("detailed") protein abundances.
//
// ProteinData::channel_level_abundances is the source of the
// "-file_and_channel_level_output true" columns of ProteinQuantifier and the named
// file/channel arrays persisted on protein groups. It has to follow the *same*
// peptide-level policy as ProteinData::fraction_group_abundances, which is computed from
// the peptidoform-merged ProteinData::peptide_fraction_group_abundances:
//   - all peptidoforms of an unmodified sequence are summed into one peptide,
//   - all charge states are summed by default; with 'best_charge_and_fraction',
//     only the globally selected charge of each peptidoform contributes,
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

START_SECTION(([EXTRA] best_charge_and_fraction selects by assay prevalence and retains all selected-charge cells))
{
  // Charge 2 is quantified in two distinct assays (channels 1 and 2), whereas the much
  // more abundant charge 3 occurs in only one (channel 3). The two fraction observations
  // in channel 1 still count as one assay. Charge 2 must therefore win globally, and all
  // of its fraction/assay cells - including the explicit zero in channel 3 - must survive.
  ConsensusMap consensus;
  ExperimentalDesign design;
  make_fractionated_input({"fileA", "fileB"}, 3,
                          {{0, 1, "PEPTIDEK", 2, 10.0},
                           {1, 1, "PEPTIDEK", 2, 5.0},
                           {0, 2, "PEPTIDEK", 2, 20.0},
                           {0, 3, "PEPTIDEK", 2, 0.0},
                           {0, 3, "PEPTIDEK", 3, 1000.0}},
                          "Prot", consensus, design);

  PeptideAndProteinQuant quantifier;
  Param p = quantifier.getDefaults();
  p.setValue("best_charge_and_fraction", "true");
  p.setValue("top:N", 0);
  p.setValue("top:aggregate", "sum");
  p.setValue("top:include_all", "true");
  quantifier.setParameters(p);

  quantifier.readQuantData(consensus, design);
  quantifier.quantifyPeptides();

  const auto& peptide = quantifier.getPeptideResults().at(AASequence::fromString("PEPTIDEK"));
  TEST_EQUAL(peptide.fraction_group_abundances.at(1).size(), 3);
  TEST_REAL_SIMILAR(peptide.fraction_group_abundances.at(1).at(1), 15.0);
  TEST_REAL_SIMILAR(peptide.fraction_group_abundances.at(1).at(2), 20.0);
  TEST_REAL_SIMILAR(peptide.fraction_group_abundances.at(1).at(3), 0.0);

  // Selection affects protein quantification, not the detailed peptide observations.
  TEST_REAL_SIMILAR(peptide.abundances.at(1).at("fileA").at(3).at(3), 1000.0);

  quantifier.quantifyProteins();
  const auto& protein = quantifier.getProteinResults().at("Prot");
  TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(1), 15.0);
  TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(2), 20.0);
  TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(3), 0.0);
  TEST_REAL_SIMILAR(protein.channel_level_abundances.at("fileA").at(1), 10.0);
  TEST_REAL_SIMILAR(protein.channel_level_abundances.at("fileA").at(2), 20.0);
  TEST_REAL_SIMILAR(protein.channel_level_abundances.at("fileA").at(3), 0.0);
  TEST_REAL_SIMILAR(protein.channel_level_abundances.at("fileB").at(1), 5.0);

  // These named arrays are consumed directly by ProteinGroupArrowExport for QPX pg.parquet.
  ProteinIdentification proteins;
  ProteinIdentification::ProteinGroup group;
  group.accessions = {"Prot"};
  group.probability = 1.0;
  proteins.getIndistinguishableProteins().push_back(group);
  quantifier.annotateQuantificationsToProteins(quantifier.getProteinResults(), proteins, true);
  const auto& annotated = proteins.getIndistinguishableProteins().front();
  const auto& qpx_abundances = annotated.getFloatDataArrays().at(3);
  const auto& qpx_fraction_groups = annotated.getIntegerDataArrays().at(2);
  const auto& qpx_labels = annotated.getIntegerDataArrays().at(3);
  TEST_STRING_EQUAL(qpx_abundances.getName(), "fraction_group_level_abundance");
  TEST_EQUAL(qpx_abundances.size(), 3);
  TEST_REAL_SIMILAR(qpx_abundances[0], 15.0);
  TEST_REAL_SIMILAR(qpx_abundances[1], 20.0);
  TEST_REAL_SIMILAR(qpx_abundances[2], 0.0);
  TEST_EQUAL(qpx_fraction_groups[0], 1);
  TEST_EQUAL(qpx_fraction_groups[1], 1);
  TEST_EQUAL(qpx_fraction_groups[2], 1);
  TEST_EQUAL(qpx_labels[0], 1);
  TEST_EQUAL(qpx_labels[1], 2);
  TEST_EQUAL(qpx_labels[2], 3);
}
END_SECTION

START_SECTION(([EXTRA] best_charge_and_fraction breaks equal-prevalence ties by total abundance))
{
  // Both charge states occur in both assays. For PEPTIDEK, charge 3 has the larger total
  // abundance and must be retained in both assays; selection must not collapse to its largest
  // single cell. AAAAAK has an exact prevalence/abundance tie and must retain lower charge 2.
  ConsensusMap consensus;
  ExperimentalDesign design;
  make_detailed_input({"fileA", "fileB"},
                      {{0, "PEPTIDEK", 2, 20.0},
                       {1, "PEPTIDEK", 2, 10.0},
                       {0, "PEPTIDEK", 3, 5.0},
                       {1, "PEPTIDEK", 3, 40.0},
                       {0, "AAAAAK", 2, 10.0},
                       {1, "AAAAAK", 2, 20.0},
                       {0, "AAAAAK", 3, 20.0},
                       {1, "AAAAAK", 3, 10.0}},
                      "Prot", consensus, design);

  PeptideAndProteinQuant quantifier;
  Param p = quantifier.getDefaults();
  p.setValue("best_charge_and_fraction", "true");
  p.setValue("top:N", 0);
  p.setValue("top:aggregate", "sum");
  quantifier.setParameters(p);

  quantifier.readQuantData(consensus, design);
  quantifier.quantifyPeptides();

  const auto& peptide = quantifier.getPeptideResults().at(AASequence::fromString("PEPTIDEK"));
  TEST_EQUAL(peptide.fraction_group_abundances.size(), 2);
  TEST_REAL_SIMILAR(peptide.fraction_group_abundances.at(1).at(1), 5.0);
  TEST_REAL_SIMILAR(peptide.fraction_group_abundances.at(2).at(1), 40.0);
  const auto& exact_tie = quantifier.getPeptideResults().at(AASequence::fromString("AAAAAK"));
  TEST_REAL_SIMILAR(exact_tie.fraction_group_abundances.at(1).at(1), 10.0);
  TEST_REAL_SIMILAR(exact_tie.fraction_group_abundances.at(2).at(1), 20.0);

  quantifier.quantifyProteins();
  const auto& protein = quantifier.getProteinResults().at("Prot");
  TEST_REAL_SIMILAR(protein.channel_level_abundances.at("fileA").at(1), 15.0);
  TEST_REAL_SIMILAR(protein.channel_level_abundances.at("fileB").at(1), 60.0);
}
END_SECTION

START_SECTION(([EXTRA] best_charge_and_fraction selects a charge independently for each peptidoform))
{
  // The unmodified peptidoform selects charge 3 (100 > 10), while its oxidized form selects
  // charge 2 (90 > 1). Protein-level assay and file/channel values must merge the two selected
  // peptidoform quantities, rather than selecting once per unmodified sequence or summing all cells.
  ConsensusMap consensus;
  ExperimentalDesign design;
  make_detailed_input({"fileA"},
                      {{0, "PEPTIDEMK", 2, 10.0},
                       {0, "PEPTIDEMK", 3, 100.0},
                       {0, "PEPTIDEM(Oxidation)K", 2, 90.0},
                       {0, "PEPTIDEM(Oxidation)K", 3, 1.0}},
                      "Prot", consensus, design);

  PeptideAndProteinQuant quantifier;
  Param p = quantifier.getDefaults();
  p.setValue("best_charge_and_fraction", "true");
  p.setValue("top:N", 0);
  p.setValue("top:aggregate", "sum");
  quantifier.setParameters(p);

  quantifier.readQuantData(consensus, design);
  quantifier.quantifyPeptides();
  TEST_REAL_SIMILAR(quantifier.getPeptideResults().at(AASequence::fromString("PEPTIDEMK")).fraction_group_abundances.at(1).at(1), 100.0);
  TEST_REAL_SIMILAR(quantifier.getPeptideResults().at(AASequence::fromString("PEPTIDEM(Oxidation)K")).fraction_group_abundances.at(1).at(1), 90.0);

  quantifier.quantifyProteins();
  const auto& protein = quantifier.getProteinResults().at("Prot");
  TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(1), 190.0);
  TEST_REAL_SIMILAR(protein.channel_level_abundances.at("fileA").at(1), 190.0);
}
END_SECTION

START_SECTION((const ProteinQuant& getProteinResults() fractionated: file+channel cells decompose the assay for sum))
{
  // One assay spans two fraction files, so the (file, channel) cells are a decomposition of
  // the assay - but only exactly so for top:N 0 with "sum", where per-file aggregation
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

  // the cells of an assay add up to that assay's value
  TEST_REAL_SIMILAR(pd.fraction_group_abundances.at(1).at(1), 650.0);
  TEST_REAL_SIMILAR(pd.fraction_group_abundances.at(1).at(2), 65.0);
  TEST_REAL_SIMILAR(pd.channel_level_abundances.at("fileA").at(1) +
                    pd.channel_level_abundances.at("fileB").at(1), pd.fraction_group_abundances.at(1).at(1));
  TEST_REAL_SIMILAR(pd.channel_level_abundances.at("fileA").at(2) +
                    pd.channel_level_abundances.at("fileB").at(2), pd.fraction_group_abundances.at(1).at(2));
}
END_SECTION

START_SECTION((const ProteinQuant& getProteinResults() fractionated: "top:N" is enforced per file, not per assay))
{
  // Same input, but with top:N = 3. Each assay has three peptides, so the protein IS
  // quantified at the assay level; no individual fraction file holds three peptides, so
  // every (file, channel) cell is dropped. This is intended - the cells are per-file values -
  // and is documented on the 'file_and_channel_level_output' parameter. It is also strictly
  // stricter than the assay-level rule, which is why it deserves pinning.
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
  TEST_EQUAL(pd.peptide_fraction_group_abundances.size(), 3);

  // quantified per assay (three peptides each) ...
  TEST_REAL_SIMILAR(pd.fraction_group_abundances.at(1).at(1), 650.0);
  TEST_REAL_SIMILAR(pd.fraction_group_abundances.at(1).at(2), 65.0);
  // ... but no fraction file holds three peptides, so no cell survives
  TEST_TRUE(pd.channel_level_abundances.empty());
}
END_SECTION

START_SECTION((void annotateQuantificationsToProteins(const ProteinQuant& protein_quants, ProteinIdentification& proteins, bool remove_unquantified) isobaric: one file, N channels))
{
  // This is the shape IsobaricAnalyzer produces and IsobaricWorkflow quantifies: ONE file
  // whose column headers differ only in label/channel. Here every (file, channel) cell is exactly
  // one assay - the two granularities are
  // degenerate and cannot disagree, which is why the per-file peptide selection and
  // aggregation documented on 'file_and_channel_level_output' are harmless for unfractionated
  // isobaric data (contrast the two fractionated sections above).
  //   fileA ch1: PEPTIDEK 100, AAAAAK 200, CCCCCK 300     fileA ch2: 10, 20, 30
  // This section also covers annotateQuantificationsToProteins(), which persists the
  // file/channel arrays and the fraction-group arrays consumed by the QPX protein-group
  // export and mzTab.
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
  // ... and they equal the assay-level values, because cell == assay here
  TEST_REAL_SIMILAR(pd.fraction_group_abundances.at(1).at(1), 200.0);
  TEST_REAL_SIMILAR(pd.fraction_group_abundances.at(1).at(2), 20.0);

  quantifier.annotateQuantificationsToProteins(quant, proteins, true);

  TEST_EQUAL(proteins.getIndistinguishableProteins().size(), 1);
  const auto& annotated = proteins.getIndistinguishableProteins()[0];

  TEST_EQUAL(annotated.getFloatDataArrays().size(), 4);
  TEST_EQUAL(annotated.getStringDataArrays().size(), 2);
  TEST_EQUAL(annotated.getIntegerDataArrays().size(), 4);

  const auto& cell_abundances = annotated.getFloatDataArrays()[2];
  const auto& fraction_group_abundances = annotated.getFloatDataArrays()[3];
  const auto& cell_filenames = annotated.getStringDataArrays()[0];
  const auto& cell_channels = annotated.getIntegerDataArrays()[0];
  const auto& fraction_groups = annotated.getIntegerDataArrays()[2];
  const auto& fraction_group_labels = annotated.getIntegerDataArrays()[3];

  TEST_STRING_EQUAL(cell_abundances.getName(), "file_channel_level_abundance");
  TEST_STRING_EQUAL(cell_filenames.getName(), "file_channel_level_filename");
  TEST_STRING_EQUAL(cell_channels.getName(), "file_channel_level_channel");
  TEST_STRING_EQUAL(fraction_group_abundances.getName(), "fraction_group_level_abundance");
  TEST_STRING_EQUAL(fraction_groups.getName(), "fraction_group_level_fraction_group");
  TEST_STRING_EQUAL(fraction_group_labels.getName(), "fraction_group_level_label");

  // one entry per (file, channel) of the design, in design order
  TEST_EQUAL(cell_abundances.size(), 2);
  TEST_EQUAL(cell_filenames.size(), 2);
  TEST_EQUAL(cell_channels.size(), 2);
  TEST_EQUAL(fraction_group_abundances.size(), 2);
  TEST_EQUAL(fraction_groups.size(), 2);
  TEST_EQUAL(fraction_group_labels.size(), 2);

  TEST_STRING_EQUAL(cell_filenames[0], "fileA");
  TEST_STRING_EQUAL(cell_filenames[1], "fileA");
  TEST_EQUAL(cell_channels[0], 1);
  TEST_EQUAL(cell_channels[1], 2);
  TEST_REAL_SIMILAR(fraction_group_abundances[0], 200.0);
  TEST_REAL_SIMILAR(fraction_group_abundances[1], 20.0);
  TEST_EQUAL(fraction_groups[0], 1);
  TEST_EQUAL(fraction_groups[1], 1);
  TEST_EQUAL(fraction_group_labels[0], 1);
  TEST_EQUAL(fraction_group_labels[1], 2);

  // The annotated cells carry the same values as the assays in this unfractionated design.
  TEST_REAL_SIMILAR(cell_abundances[0], 200.0);
  TEST_REAL_SIMILAR(cell_abundances[1], 20.0);
  TEST_REAL_SIMILAR(cell_abundances[0], fraction_group_abundances[0]);
  TEST_REAL_SIMILAR(cell_abundances[1], fraction_group_abundances[1]);
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
  TEST_EQUAL(pd.peptide_fraction_group_abundances.size(), 1);
  TEST_REAL_SIMILAR(pd.fraction_group_abundances.at(1).at(1), 100.0);

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
  // the assay-level result for fileA is dropped. Charge states must not be
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
  TEST_EQUAL(pd.peptide_fraction_group_abundances.size(), 3); // clears the protein-level gate

  // assay level: fileA (fraction group 1) has only one peptide -> no abundance
  TEST_TRUE(pd.fraction_group_abundances.find(1) == pd.fraction_group_abundances.end());
  TEST_REAL_SIMILAR(pd.fraction_group_abundances.at(2).at(1), 30.0);

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
  TEST_EQUAL(pd.peptide_fraction_group_abundances.size(), 3);

  // assay level drops fileA (fraction group 1)
  TEST_TRUE(pd.fraction_group_abundances.find(1) == pd.fraction_group_abundances.end());
  TEST_REAL_SIMILAR(pd.fraction_group_abundances.at(2).at(1), 3.0);

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
  TEST_EQUAL(pd.peptide_fraction_group_abundances.size(), 3);

  TEST_REAL_SIMILAR(pd.fraction_group_abundances.at(1).at(1), 135.0);
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
  TEST_REAL_SIMILAR(pd.fraction_group_abundances.at(1).at(1), 100.0);
  TEST_REAL_SIMILAR(pd.channel_level_abundances.at("fileA").at(1), 100.0);
}
END_SECTION


START_SECTION((const ProteinQuant& getProteinResults() file+channel level aggregates detected peptides, not stored zeros))
{
  // The protein is measured at 1000 in every channel of the one file. In channel 2 only one of
  // its three peptides has a reporter; the other two are stored as 0.0, which is what
  // IsobaricChannelExtractor writes for a reporter it could not find. Aggregating those zeros
  // reports the cell as empty for a channel the protein was measured in - the same defect as at
  // the sample level, in the array that feeds the abundance|<file>|ch<N> columns and the QPX pg
  // view.
  ConsensusMap consensus;
  ExperimentalDesign design;
  make_fractionated_input({"fileA"}, 3,
                          {{0, 1, "PEPTIDEK", 2, 1000.0}, {0, 2, "PEPTIDEK", 2, 1000.0}, {0, 3, "PEPTIDEK", 2, 1000.0},
                           {0, 1, "AAAAAK",   2, 1000.0}, {0, 2, "AAAAAK",   2,    0.0}, {0, 3, "AAAAAK",   2, 1000.0},
                           {0, 1, "CCCCCK",   2, 1000.0}, {0, 2, "CCCCCK",   2,    0.0}, {0, 3, "CCCCCK",   2, 1000.0}},
                          "Prot", consensus, design);

  // 'top:include_all' is what IsobaricWorkflow sets, so this is the shipped isobaric path.
  {
    PeptideAndProteinQuant quantifier;
    Param p = quantifier.getDefaults();
    p.setValue("top:aggregate", "median");
    p.setValue("top:include_all", "true");
    quantifier.setParameters(p);
    quantifier.readQuantData(consensus, design);
    quantifier.quantifyPeptides();
    quantifier.quantifyProteins();

    const auto& protein = quantifier.getProteinResults().at("Prot");
    TEST_REAL_SIMILAR(protein.channel_level_abundances.at("fileA").at(1), 1000.0)
    TEST_REAL_SIMILAR(protein.channel_level_abundances.at("fileA").at(2), 1000.0) // was 0.0
    TEST_REAL_SIMILAR(protein.channel_level_abundances.at("fileA").at(3), 1000.0)
  }

  // 'sum' documents the opposite convention - missing values count as zero - and summing zeros
  // changes nothing, so its cells must not move.
  {
    PeptideAndProteinQuant quantifier;
    Param p = quantifier.getDefaults();
    p.setValue("top:N", 0);
    p.setValue("top:aggregate", "sum");
    quantifier.setParameters(p);
    quantifier.readQuantData(consensus, design);
    quantifier.quantifyPeptides();
    quantifier.quantifyProteins();

    const auto& protein = quantifier.getProteinResults().at("Prot");
    TEST_REAL_SIMILAR(protein.channel_level_abundances.at("fileA").at(1), 3000.0)
    TEST_REAL_SIMILAR(protein.channel_level_abundances.at("fileA").at(2), 1000.0)
    TEST_REAL_SIMILAR(protein.channel_level_abundances.at("fileA").at(3), 3000.0)
  }

  // The per-file "at least N peptides" gate counts detected peptides too: channel 2 has one, so
  // without 'include_all' and with 'top:N' 3 the cell is declined rather than credited with
  // peptides that were never measured. The cell keeps no entry; the writers render that as 0.0.
  {
    PeptideAndProteinQuant quantifier;
    Param p = quantifier.getDefaults();
    p.setValue("top:N", 3);
    p.setValue("top:aggregate", "median");
    quantifier.setParameters(p);
    quantifier.readQuantData(consensus, design);
    quantifier.quantifyPeptides();
    quantifier.quantifyProteins();

    const auto& cells = quantifier.getProteinResults().at("Prot").channel_level_abundances.at("fileA");
    TEST_REAL_SIMILAR(cells.at(1), 1000.0)
    TEST_TRUE(cells.find(2) == cells.end())
    TEST_REAL_SIMILAR(cells.at(3), 1000.0)
  }
}
END_SECTION


START_SECTION((const ProteinQuant& getProteinResults() fraction group level aggregates detected peptides, not stored zeros))
{
  // Third copy of the same aggregation, and it had the same defect: the per-(fraction group,
  // label) cells were built by pushing every stored abundance, zeros included. The fixture is the
  // one used for the file+channel level above - one file, three channels, and in channel 2 only
  // one of the protein's three peptides carries a reporter - so all three levels describe the
  // same measurement and must agree. Before this fix they did not: IsobaricWorkflow_1's reference
  // showed the fraction-group array at half the detected-value result wherever a peptide went
  // undetected, the signature of mean(x, 0).
  ConsensusMap consensus;
  ExperimentalDesign design;
  make_fractionated_input({"fileA"}, 3,
                          {{0, 1, "PEPTIDEK", 2, 1000.0}, {0, 2, "PEPTIDEK", 2, 1000.0}, {0, 3, "PEPTIDEK", 2, 1000.0},
                           {0, 1, "AAAAAK",   2, 1000.0}, {0, 2, "AAAAAK",   2,    0.0}, {0, 3, "AAAAAK",   2, 1000.0},
                           {0, 1, "CCCCCK",   2, 1000.0}, {0, 2, "CCCCCK",   2,    0.0}, {0, 3, "CCCCCK",   2, 1000.0}},
                          "Prot", consensus, design);

  // The helper puts every file in fraction group 1 and labels the channels 1..N, so the cells are
  // (1, channel) and correspond one-to-one to the assays.
  {
    PeptideAndProteinQuant quantifier;
    Param p = quantifier.getDefaults();
    p.setValue("top:aggregate", "median");
    p.setValue("top:include_all", "true");
    quantifier.setParameters(p);
    quantifier.readQuantData(consensus, design);
    quantifier.quantifyPeptides();
    quantifier.quantifyProteins();

    const auto& protein = quantifier.getProteinResults().at("Prot");
    TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(1), 1000.0)
    TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(2), 1000.0) // was 0.0
    TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(3), 1000.0)

    // One file x N channels means a cell is a sample, so the three levels must report the same
    // number. This is the invariant IsobaricWorkflow_1's reference file encodes.
    TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(2),
                      protein.channel_level_abundances.at("fileA").at(2))
  }

  // 'sum' counts missing values as zero by documented convention, and summing zeros changes
  // nothing, so these cells must not move.
  {
    PeptideAndProteinQuant quantifier;
    Param p = quantifier.getDefaults();
    p.setValue("top:N", 0);
    p.setValue("top:aggregate", "sum");
    quantifier.setParameters(p);
    quantifier.readQuantData(consensus, design);
    quantifier.quantifyPeptides();
    quantifier.quantifyProteins();

    const auto& protein = quantifier.getProteinResults().at("Prot");
    TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(1), 3000.0)
    TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(2), 1000.0)
    TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(3), 3000.0)
  }

  // The "at least N peptides" gate counts detected peptides here too: channel 2 has one, so with
  // 'top:N' 3 and no 'include_all' the cell is declined rather than credited with peptides that
  // were never measured. It keeps no entry; the writers render that as 0.0.
  {
    PeptideAndProteinQuant quantifier;
    Param p = quantifier.getDefaults();
    p.setValue("top:N", 3);
    p.setValue("top:aggregate", "median");
    quantifier.setParameters(p);
    quantifier.readQuantData(consensus, design);
    quantifier.quantifyPeptides();
    quantifier.quantifyProteins();

    const auto& cells = quantifier.getProteinResults().at("Prot").fraction_group_abundances.at(1);
    TEST_REAL_SIMILAR(cells.at(1), 1000.0)
    TEST_TRUE(cells.find(2) == cells.end())
    TEST_REAL_SIMILAR(cells.at(3), 1000.0)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// Isobaric map, one file, @p n_channels channels, one consensus feature per peptide with one
// sub-feature per channel. `peptides[i].second[c]` is peptide i's reporter intensity in channel c;
// a 0 is a reporter IsobaricChannelExtractor could not find, which it writes out as 0.0. All
// peptides are assigned to one protein so quantifyProteins() has something to aggregate.
auto makeIsobaricProteinInput = [](
    const std::vector<std::pair<std::string, std::vector<double>>>& peptides,
    ConsensusMap& consensus,
    ExperimentalDesign& design)
  {
    const Size n_channels = peptides.front().second.size();

    consensus.setExperimentType("labeled_MS2");
    for (Size c = 0; c < n_channels; ++c)
    {
      ConsensusMap::ColumnHeader h;
      h.filename = "A.mzML";
      h.label = "ch" + std::to_string(c + 1);
      h.setMetaValue("channel_id", (int)c); // label = channel_id + 1
      consensus.getColumnHeaders()[c] = h;
    }

    for (Size p = 0; p < peptides.size(); ++p)
    {
      ConsensusFeature cf;
      for (Size c = 0; c < n_channels; ++c)
      {
        Peak2D peak;
        peak.setIntensity(peptides[p].second[c]);
        peak.setRT(1.0 + p);
        peak.setMZ(400.0 + p);
        cf.insert(c, peak, c);
      }
      PeptideEvidence ev;
      ev.setProteinAccession("Prot");
      PeptideHit hit;
      hit.setSequence(AASequence::fromString(peptides[p].first));
      hit.setCharge(2);
      hit.setScore(1.0);
      hit.setPeptideEvidences({ev});
      PeptideIdentification pid;
      pid.setHits({hit});
      cf.setPeptideIdentifications({pid});
      consensus.push_back(cf);
    }

    ExperimentalDesign::MSFileSection fs;
    ExperimentalDesign::SampleSection ss;
    for (Size c = 0; c < n_channels; ++c)
    {
      ExperimentalDesign::MSFileSectionEntry r;
      r.path = "/tmp/A.mzML";
      r.label = c + 1;
      r.fraction = 1;
      r.fraction_group = 1;
      r.sample = c;
      r.sample_name = "S" + std::to_string(c + 1);
      fs.push_back(r);
      ss.addSample("S" + std::to_string(c + 1));
    }
    design = ExperimentalDesign(fs, ss);
  };

START_SECTION(([EXTRA] "occurs in every assay" counts detected assays, not stored zeros))
{
  // 'consensus:fix_peptides' with 'top:N' 0 keeps the peptides that occur in EVERY assay.
  // quantifyFeature_ records an undetected reporter as an explicit 0.0, so the assay key of a
  // peptide exists whether or not it was detected there. Counting keys makes every peptide seen
  // anywhere pass, and the filter stops filtering.
  ConsensusMap consensus;
  ExperimentalDesign design;
  makeIsobaricProteinInput({{"PEPTIDEK", {100.0, 100.0, 100.0, 100.0}},   // detected in all four
                            {"AAAAAK",   {500.0, 500.0,   0.0,   0.0}}},  // detected in two
                           consensus, design);

  PeptideAndProteinQuant quantifier;
  Param p = quantifier.getParameters();
  p.setValue("consensus:fix_peptides", "true");
  p.setValue("top:N", 0);
  p.setValue("top:aggregate", "sum");
  quantifier.setParameters(p);

  quantifier.readQuantData(consensus, design);
  quantifier.quantifyPeptides();
  quantifier.quantifyProteins();

  // Only PEPTIDEK qualifies, so the protein is its 100 per assay. Counting the stored zeros
  // would admit AAAAAK too and report 600, 600, 100, 100.
  const auto& protein = quantifier.getProteinResults().at("Prot");
  TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(1), 100.0)
  TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(2), 100.0)
  TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(3), 100.0)
  TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(4), 100.0)
}
END_SECTION

START_SECTION(([EXTRA] top:N ranking prefers the peptide detected in more assays, not the taller one))
{
  // orderBest_ ranks by (number of assays quantified, total abundance). With stored zeros
  // counted, every peptide of a TMT run ties at the channel count and the tie-break by total
  // abundance decides alone - so a peptide with a big signal in 2 of 10 channels takes the
  // top:N slot from one genuinely measured in 6.
  ConsensusMap consensus;
  ExperimentalDesign design;
  makeIsobaricProteinInput({{"BROADK", {100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 0.0, 0.0, 0.0, 0.0}},
                            {"TALLK",  {1000.0, 1000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
                            {"SMALLK", {10.0, 10.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}},
                           consensus, design);

  PeptideAndProteinQuant quantifier;
  Param p = quantifier.getParameters();
  p.setValue("consensus:fix_peptides", "true");
  p.setValue("top:N", 1);
  p.setValue("top:aggregate", "sum");
  quantifier.setParameters(p);

  quantifier.readQuantData(consensus, design);
  quantifier.quantifyPeptides();
  quantifier.quantifyProteins();

  // BROADK is quantified in 6 assays, TALLK in 2, so BROADK takes the single slot and the
  // protein reports 100 in the first six assays. Ranking on the stored zeros would tie all
  // three at 10 assays and hand the slot to TALLK's larger total, giving 1000, 1000, 0, ...
  const auto& protein = quantifier.getProteinResults().at("Prot");
  TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(1), 100.0)
  TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(2), 100.0)
  TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(6), 100.0)
}
END_SECTION

START_SECTION(([EXTRA] median/mean aggregate the detected peptides, not the stored zeros))
{
  // The protein is measured at 1000 in every channel. In channel 2 only one of its three
  // peptides has a reporter; the other two are stored as 0.0, which is what
  // IsobaricChannelExtractor writes for a reporter it could not find. Aggregating those zeros
  // reports the protein as ABSENT from a sample it was measured in - the artefact
  // normalizePeptides_() already excludes zeros to avoid, and the opposite of what this tool
  // documents: "mean and median ignore missing cases, averaging only present values".
  ConsensusMap consensus;
  ExperimentalDesign design;
  makeIsobaricProteinInput({{"PEPTIDEK", {1000.0, 1000.0, 1000.0}},
                            {"AAAAAK",   {1000.0,    0.0, 1000.0}},
                            {"CCCCCK",   {1000.0,    0.0, 1000.0}}},
                           consensus, design);

  // 'top:include_all' is what IsobaricWorkflow sets, so this is the shipped isobaric path.
  for (const std::string& aggregate : {"median", "mean"})
  {
    PeptideAndProteinQuant quantifier;
    Param p = quantifier.getParameters();
    p.setValue("top:aggregate", aggregate);
    p.setValue("top:include_all", "true");
    quantifier.setParameters(p);

    quantifier.readQuantData(consensus, design);
    quantifier.quantifyPeptides();
    quantifier.quantifyProteins();

    const auto& protein = quantifier.getProteinResults().at("Prot");
    TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(1), 1000.0)
    TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(2), 1000.0) // was 0.0 (median) / 333.33 (mean)
    TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(3), 1000.0)
  }

  // 'sum' documents the opposite convention - missing values count as zero - and summing zeros
  // changes nothing, so its numbers must not move.
  {
    PeptideAndProteinQuant quantifier;
    Param p = quantifier.getParameters();
    p.setValue("top:aggregate", "sum");
    p.setValue("top:include_all", "true");
    quantifier.setParameters(p);

    quantifier.readQuantData(consensus, design);
    quantifier.quantifyPeptides();
    quantifier.quantifyProteins();

    const auto& protein = quantifier.getProteinResults().at("Prot");
    TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(1), 3000.0)
    TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(2), 1000.0)
    TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(3), 3000.0)
  }
}
END_SECTION

START_SECTION(([EXTRA] "enough peptides" counts the detected ones, and a dead sample still reports))
{
  ConsensusMap consensus;
  ExperimentalDesign design;
  makeIsobaricProteinInput({{"PEPTIDEK", {1000.0, 1000.0, 0.0}},
                            {"AAAAAK",   {1000.0,    0.0, 0.0}},
                            {"CCCCCK",   {1000.0,    0.0, 0.0}}},
                           consensus, design);

  // Without 'include_all' and with 'top:N' 3, sample 2 has one detected peptide and sample 3 has
  // none, so neither reaches the threshold and both are left unquantified rather than being
  // credited with peptides that were never measured.
  {
    PeptideAndProteinQuant quantifier;
    Param p = quantifier.getParameters();
    p.setValue("top:N", 3);
    p.setValue("top:aggregate", "median");
    quantifier.setParameters(p);
    quantifier.readQuantData(consensus, design);
    quantifier.quantifyPeptides();
    quantifier.quantifyProteins();

    const auto& protein = quantifier.getProteinResults().at("Prot");
    const auto& assays = protein.fraction_group_abundances.at(1);
    TEST_REAL_SIMILAR(assays.at(1), 1000.0)
    TEST_TRUE(assays.find(2) == assays.end())
    TEST_TRUE(assays.find(3) == assays.end())
  }

  // With 'include_all' the gate is bypassed: sample 2 reports its one detected peptide, and
  // sample 3 - where nothing was detected at all - keeps a zero entry rather than vanishing.
  {
    PeptideAndProteinQuant quantifier;
    Param p = quantifier.getParameters();
    p.setValue("top:N", 3);
    p.setValue("top:aggregate", "median");
    p.setValue("top:include_all", "true");
    quantifier.setParameters(p);
    quantifier.readQuantData(consensus, design);
    quantifier.quantifyPeptides();
    quantifier.quantifyProteins();

    const auto& protein = quantifier.getProteinResults().at("Prot");
    TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(1), 1000.0)
    TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(2), 1000.0)
    TEST_REAL_SIMILAR(protein.fraction_group_abundances.at(1).at(3), 0.0)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
