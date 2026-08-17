// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Timo Sachsenberg$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/FORMAT/ExperimentalDesignFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/SYSTEM/File.h>
#include <fstream>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ExperimentalDesign, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ExperimentalDesign* ptr = 0;
ExperimentalDesign* null_ptr = 0;

ExperimentalDesign labelfree_unfractionated_design = ExperimentalDesignFile::load(
  OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_input_1.tsv")
  , false);

ExperimentalDesign fourplex_fractionated_design = ExperimentalDesignFile::load(
  OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_input_2.tsv")
  , false);

ExperimentalDesign labelfree_unfractionated_single_table_design = ExperimentalDesignFile::load(
  OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_input_1_single_table.tsv")
  , false);

ExperimentalDesign fourplex_fractionated_single_table_design = ExperimentalDesignFile::load(
  OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_input_2_single_table.tsv")
  , false);

ExperimentalDesign labelfree_unfractionated_single_table_no_sample_column = ExperimentalDesignFile::load(
  OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_input_3_single_table.tsv")
  , false);

START_SECTION(ExperimentalDesign())
{
  ptr = new ExperimentalDesign();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~ExperimentalDesign())
{
  delete ptr;
}
END_SECTION

START_SECTION((ExperimentalDesign(MSFileSection msfile_section, SampleSection sample_section)))
{
  ExperimentalDesign::MSFileSection fs;
  ExperimentalDesign::SampleSection ss;
  ExperimentalDesign ed(fs, ss);
}
END_SECTION

START_SECTION((const MSFileSection& getMSFileSection() const ))
{
  ExperimentalDesign::MSFileSection fs = labelfree_unfractionated_design.getMSFileSection();
}
END_SECTION

START_SECTION((void setMSFileSection(const MSFileSection &msfile_section)))
{
  ExperimentalDesign labelfree_unfractionated_design2 = labelfree_unfractionated_design;
  ExperimentalDesign::MSFileSection fs;
  labelfree_unfractionated_design2.setMSFileSection(fs);
}
END_SECTION

START_SECTION((const ExperimentalDesign::SampleSection& getSampleSection() const ))
{
  ExperimentalDesign::SampleSection ss = labelfree_unfractionated_design.getSampleSection();
}
END_SECTION

START_SECTION((unsigned getSample(unsigned fraction_group, unsigned label)))
{
  // an unknown (fraction_group, label) combination must throw instead of dereferencing end() (undefined behavior)
  TEST_EXCEPTION(Exception::ElementNotFound, labelfree_unfractionated_design.getSample(99999, 99999))
}
END_SECTION

START_SECTION((void setSampleSection(const ExperimentalDesign::SampleSection &sample_section)))
{
  ExperimentalDesign labelfree_unfractionated_design2 = labelfree_unfractionated_design;
  ExperimentalDesign::SampleSection fs;
  labelfree_unfractionated_design2.setSampleSection(fs);
}
END_SECTION

START_SECTION((std::map<unsigned int, std::vector<std::string> > getFractionToMSFilesMapping() const ))
{
  const auto lf = labelfree_unfractionated_design.getFractionToMSFilesMapping();
  const auto lfst = labelfree_unfractionated_single_table_design.getFractionToMSFilesMapping();
  const auto lfstns = labelfree_unfractionated_single_table_no_sample_column.getFractionToMSFilesMapping();
  // test both two table as well as single table design
  for (const auto& f2ms : { lf, lfst, lfstns})
  { 
    // unfractionated data so only one fraction
    TEST_EQUAL(f2ms.size(), 1);
    // we have unfactionated data so fraction 1 mapps to all 12 files
    TEST_EQUAL(f2ms.at(1).size(), 12);
  }

  const auto fplex = fourplex_fractionated_design.getFractionToMSFilesMapping();
  const auto fplexst = fourplex_fractionated_single_table_design.getFractionToMSFilesMapping();
  // test both two table as well as single table design
  for (const auto& f2ms : { fplex, fplexst})
  { 
    // tripple fractionated data
    TEST_EQUAL(f2ms.size(), 3);
    // three fractions, 24 files, fraction 1..3 map to 8 files each
    TEST_EQUAL(f2ms.at(1).size(), 8);
    TEST_EQUAL(f2ms.at(2).size(), 8);
    TEST_EQUAL(f2ms.at(3).size(), 8);
  }
}
END_SECTION

START_SECTION(( std::map<std::vector<String>, std::set<unsigned> > getConditionToSampleMapping() const ))
{
  ExperimentalDesign ed;
  
  // add 3 samples without factors to simulate a factor-less design
  ed.addSample("sample_1");
  ed.addSample("sample_2");
  ed.addSample("sample_3");
  
  auto map_a = ed.getSampleToConditionMapping();
  auto map_b = ed.getConditionToSampleMapping();
  
  // both mapping functions should return 3 unique conditions (N samples = N conditions)
  TEST_EQUAL(map_a.size(), 3);
  TEST_EQUAL(map_b.size(), 3);

  // stronger assertions: Verify that condition IDs are actually distinct
  std::set<unsigned> condition_ids;
  for (const auto& [sample, condition] : map_a)
  {
    condition_ids.insert(condition);
  }
  TEST_EQUAL(condition_ids.size(), 3);

  // verify that each map_b entry contains exactly one sample row
  for (const auto& [condition, rows] : map_b)
  {
    TEST_EQUAL(condition.size(), 1);
    TEST_EQUAL(rows.size(), 1);
  }
}
END_SECTION

START_SECTION((std::map< std::pair< String, unsigned >, unsigned> getPathLabelToSampleMapping(bool) const ))
{
  const auto lf = labelfree_unfractionated_design.getPathLabelToSampleMapping(true);
  const auto lfst = labelfree_unfractionated_single_table_design.getPathLabelToSampleMapping(true);
  const auto lfstns = labelfree_unfractionated_single_table_no_sample_column.getPathLabelToSampleMapping(true);
  const auto fplex = fourplex_fractionated_design.getPathLabelToSampleMapping(true);
  const auto fplexst = fourplex_fractionated_single_table_design.getPathLabelToSampleMapping(true);

  // 12 quant. values from label-free, unfractionated files map to 12 samples
  for (const auto& pl2s : {lf, lfst, lfstns})
  {
    TEST_EQUAL(pl2s.size(), 12);
  } 

  // 24 quant. values from 4plex, triple fractionated files map to 8 samples
  for (const auto& pl2s : {fplex, fplexst})
  {
    TEST_EQUAL(pl2s.size(), 24);
    for (const auto& i : pl2s) { TEST_EQUAL((i.second <=7), true)}
  } 
}
END_SECTION

START_SECTION((std::map< std::pair< String, unsigned >, unsigned> getPathLabelToSampleMapping(bool) should reject ambiguous basename mappings))
{
  ExperimentalDesign::MSFileSection fs;
  ExperimentalDesign::MSFileSectionEntry r1;
  r1.fraction_group = 1;
  r1.fraction = 1;
  r1.path = "/tmp/run_a/shared_name.mzML";
  r1.label = 1;
  r1.sample = 0;
  r1.sample_name = "S1";
  fs.push_back(r1);

  ExperimentalDesign::MSFileSectionEntry r2 = r1;
  r2.fraction_group = 2;
  r2.path = "/tmp/run_b/shared_name.mzML";
  r2.sample = 1;
  r2.sample_name = "S2";
  fs.push_back(r2);

  ExperimentalDesign::SampleSection ss;
  ss.addSample("S1");
  ss.addSample("S2");

  ExperimentalDesign design(fs, ss);
  TEST_EXCEPTION(Exception::InvalidValue, design.getPathLabelToSampleMapping(true));
}
END_SECTION

START_SECTION((std::map< std::pair< String, unsigned >, unsigned> getPathLabelToFractionMapping(bool) const ))
{
  const auto lf = labelfree_unfractionated_design.getPathLabelToFractionMapping(true);
  const auto lfst = labelfree_unfractionated_single_table_design.getPathLabelToFractionMapping(true);
  const auto lfstns = labelfree_unfractionated_single_table_no_sample_column.getPathLabelToFractionMapping(true);
  const auto fplex = fourplex_fractionated_design.getPathLabelToFractionMapping(true);
  const auto fplexst = fourplex_fractionated_single_table_design.getPathLabelToFractionMapping(true);

  // 12 quant. values from label-free, unfractionated files map to fraction 1 each
  for (const auto& pl2f : { lf, lfst, lfstns })
  {
    TEST_EQUAL(pl2f.size(), 12);
    for (const auto& i : pl2f) { TEST_EQUAL(i.second, 1); }
  }

  // 24 quant. values map to fractions 1..3
  for (const auto& pl2f : { fplex, fplexst})
  {
    TEST_EQUAL(pl2f.size(), 24);
    for (const auto& i : pl2f) { TEST_EQUAL((i.second >=1 && i.second <=3), true)}
  }
}
END_SECTION

START_SECTION((std::map< std::pair< String, unsigned >, unsigned> getPathLabelToFractionGroupMapping(bool) const ))
{
  const auto lf = labelfree_unfractionated_design.getPathLabelToFractionGroupMapping(true);
  const auto lfst = labelfree_unfractionated_single_table_design.getPathLabelToFractionGroupMapping(true);
  const auto lfstns = labelfree_unfractionated_single_table_no_sample_column.getPathLabelToFractionGroupMapping(true);
  const auto fplex = fourplex_fractionated_design.getPathLabelToFractionGroupMapping(true);
  const auto fplexst = fourplex_fractionated_single_table_design.getPathLabelToFractionGroupMapping(true);

  // 12 quant. values from label-free, unfractionated files map to different fraction groups
  for (const auto& pl2fg : { lf, lfst})
  {
    TEST_EQUAL(pl2fg.size(), 12);
    int count = 1;
    for (const auto& i : pl2fg) { TEST_EQUAL(i.second, count); ++count; }
  }

  for (const auto& pl2fg : { fplex, fplexst})
  {
    TEST_EQUAL(pl2fg.size(), 24);
    for (const auto& i : pl2fg) 
    {
      // extract fraction group from file name
      int file(1); 
      if (StringUtils::hasSubstring(i.first.first, "TR2")) { file = 2; }
      TEST_EQUAL(i.second, file); 
    }
  }
}
END_SECTION

START_SECTION((std::set<std::string> ExperimentalDesign::SampleSection::getFactors() const))
  const auto lfac = labelfree_unfractionated_design.getSampleSection().getFactors();
  const auto lfacst = labelfree_unfractionated_single_table_design.getSampleSection().getFactors();
  const auto lfacstns = labelfree_unfractionated_single_table_no_sample_column.getSampleSection().getFactors();
  const auto facplex = fourplex_fractionated_design.getSampleSection().getFactors();
  const auto facplexst = fourplex_fractionated_single_table_design.getSampleSection().getFactors();

  TEST_EQUAL(lfac.size(), 3)
  TEST_EQUAL(lfacst.size(), 3)
  TEST_EQUAL(lfacstns.size(), 3)

  TEST_TRUE(lfac == lfacst)
  TEST_TRUE(lfac == lfacstns)
  TEST_TRUE(facplex == facplexst)

  auto l = lfac.begin();
  TEST_EQUAL(*l++, "MSstats_BioReplicate")
  TEST_EQUAL(*l++, "MSstats_Condition")
  TEST_EQUAL(*l++, "Sample")

  l = lfacst.begin();
  TEST_EQUAL(*l++, "MSstats_BioReplicate")
  TEST_EQUAL(*l++, "MSstats_Condition")
  TEST_EQUAL(*l++, "Sample")

  l = lfacstns.begin();
  TEST_EQUAL(*l++, "MSstats_BioReplicate")
  TEST_EQUAL(*l++, "MSstats_Condition")
  TEST_EQUAL(*l++, "Sample") // dummy sample get's automatically added if not present in ED file
END_SECTION

START_SECTION((unsigned getNumberOfSamples() const ))
{
  const auto lf = labelfree_unfractionated_design.getNumberOfSamples();
  const auto lfst = labelfree_unfractionated_single_table_design.getNumberOfSamples();
  const auto lfstns = labelfree_unfractionated_single_table_no_sample_column.getNumberOfSamples();

  for (const auto& ns : { lf, lfst, lfstns })
  {
    TEST_EQUAL(ns, 12);
  }

  const auto fplex = fourplex_fractionated_design.getNumberOfSamples();
  const auto fplexst = fourplex_fractionated_single_table_design.getNumberOfSamples();
  for (const auto& ns : { fplex, fplexst })
  {
    TEST_EQUAL(ns, 8);
  }
}
END_SECTION


START_SECTION((std::string SampleSection::getFactorValue(const unsigned sample, const std::string &factor) const))
  // Note: Number of samples are the same (correctness tested in ExperimentalDesign::SampleSection::getNumberOfSamples())
  // Note: Factors are the same (correctness tested in ExperimentalDesign::SampleSection::getFactors())
  const auto lns = labelfree_unfractionated_design.getNumberOfSamples();

  auto lss_tt = labelfree_unfractionated_design.getSampleSection();
  auto lss_st = labelfree_unfractionated_single_table_design.getSampleSection();
  auto lss_stns = labelfree_unfractionated_single_table_no_sample_column.getSampleSection();

  // 12 samples (see getNumberOfSamples test)
  for (size_t sample = 0; sample < lns; ++sample)
  {
    for (const auto& factor : lss_tt.getFactors())
    {
      // check if single table and two table design agree
      std::string f1 = lss_st.getFactorValue(sample, factor);
      std::string f2 = lss_tt.getFactorValue(sample, factor);
      std::string f3 = lss_stns.getFactorValue(sample, factor);
      cout << f1 << "\t" << f2 << "\t" << f3 << endl;
      TEST_EQUAL(f1, f2);
      TEST_EQUAL(f1, f3);
    }
  }    

  const auto fns = fourplex_fractionated_design.getNumberOfSamples();
  const auto& fss_tt = fourplex_fractionated_design.getSampleSection();
  const auto& fss_st = fourplex_fractionated_single_table_design.getSampleSection();
  // 8 samples (see getNumberOfSamples test)
  for (size_t sample = 0; sample < fns; ++sample)
  {
    for (const auto& factor : fss_tt.getFactors())
    {
      // check if single table and two table design agree
      std::string f1 = fss_st.getFactorValue(sample, factor);
      std::string f2 = fss_tt.getFactorValue(sample, factor);
      TEST_EQUAL(f1, f2);
    }
  }      

END_SECTION

START_SECTION((unsigned getNumberOfFractions() const ))
{
  const auto lf = labelfree_unfractionated_design.getNumberOfFractions();
  const auto lfst = labelfree_unfractionated_single_table_design.getNumberOfFractions();
  const auto lfstns = labelfree_unfractionated_single_table_no_sample_column.getNumberOfFractions();
  const auto fplex = fourplex_fractionated_design.getNumberOfFractions();
  const auto fplexst = fourplex_fractionated_single_table_design.getNumberOfFractions();

  for (const auto& ns : { lf, lfst, lfstns})
  {
    TEST_EQUAL(ns, 1);
  }

  for (const auto& ns : { fplex, fplexst})
  {
    TEST_EQUAL(ns, 3);
  }
}
END_SECTION

START_SECTION((unsigned getNumberOfLabels() const ))
{
  const auto lf = labelfree_unfractionated_design.getNumberOfLabels();
  const auto lfst = labelfree_unfractionated_single_table_design.getNumberOfLabels();
  const auto lfstns = labelfree_unfractionated_single_table_no_sample_column.getNumberOfLabels();
  const auto fplex = fourplex_fractionated_design.getNumberOfLabels();
  const auto fplexst = fourplex_fractionated_single_table_design.getNumberOfLabels();

  for (const auto& ns : { lf, lfst, lfstns })
  {
    TEST_EQUAL(ns, 1);
  }

  for (const auto& ns : { fplex, fplexst})
  {
    TEST_EQUAL(ns, 4);
  }
}
END_SECTION

START_SECTION((unsigned getNumberOfMSFiles() const ))
{
  const auto lf = labelfree_unfractionated_design.getNumberOfMSFiles();
  const auto lfst = labelfree_unfractionated_single_table_design.getNumberOfMSFiles();
  const auto lfstns = labelfree_unfractionated_single_table_no_sample_column.getNumberOfMSFiles();
  const auto fplex = fourplex_fractionated_design.getNumberOfMSFiles();
  const auto fplexst = fourplex_fractionated_single_table_design.getNumberOfMSFiles();

  for (const auto& ns : {lf, lfst, lfstns} )
  {
    TEST_EQUAL(ns, 12);
  }

  for (const auto& ns : {fplex, fplexst})
  {
    TEST_EQUAL(ns, 6);
  }
}
END_SECTION

START_SECTION((unsigned getNumberOfFractionGroups() const ))
{
  const auto lf = labelfree_unfractionated_design.getNumberOfFractionGroups();
  const auto lfst = labelfree_unfractionated_single_table_design.getNumberOfFractionGroups();
  const auto lfstns = labelfree_unfractionated_single_table_no_sample_column.getNumberOfFractionGroups();
  const auto fplex = fourplex_fractionated_design.getNumberOfFractionGroups();
  const auto fplexst = fourplex_fractionated_single_table_design.getNumberOfFractionGroups();

  for (const auto& ns : {lf, lfst, lfstns} )
  {
    TEST_EQUAL(ns, 12);
  }

  for (const auto& ns : {fplex, fplexst})
  {
    TEST_EQUAL(ns, 2);
  }
}
END_SECTION

START_SECTION((unsigned getSample(unsigned fraction_group, unsigned label=1)))
{
  const auto lf11 = labelfree_unfractionated_design.getSample(1, 1);
  const auto lfst11 = labelfree_unfractionated_single_table_design.getSample(1, 1);
  const auto lfstns11 = labelfree_unfractionated_single_table_no_sample_column.getSample(1, 1);
  const auto fplex11 = fourplex_fractionated_design.getSample(1, 1);
  const auto fplexst11 = fourplex_fractionated_single_table_design.getSample(1, 1);

  for (const auto& s : {lf11, lfst11, lfstns11})
  {
    TEST_EQUAL(s, 0);
  }
  for (const auto& s : {fplex11, fplexst11})
  {
    TEST_EQUAL(s, 0);
  }

  const auto lf12_1 = labelfree_unfractionated_design.getSample(12, 1);
  const auto lfst12_1 = labelfree_unfractionated_single_table_design.getSample(12, 1);
  const auto lfstns11_1 = labelfree_unfractionated_single_table_no_sample_column.getSample(12, 1);
  for (const auto& s : {lf12_1, lfst12_1, lfstns11_1})
  {
    TEST_EQUAL(s, 11);
  }

  const auto fplex24 = fourplex_fractionated_design.getSample(2, 4);
  const auto fplexst24 = fourplex_fractionated_single_table_design.getSample(2, 4);
  for (const auto& s : {fplex24, fplexst24})
  {
    TEST_EQUAL(s, 7);
  }
}
END_SECTION

START_SECTION((bool isFractionated() const ))
{
  bool lf = labelfree_unfractionated_design.isFractionated();
  bool lfst = labelfree_unfractionated_single_table_design.isFractionated();
  bool lfstns = labelfree_unfractionated_single_table_no_sample_column.isFractionated();
  bool fplex = fourplex_fractionated_design.isFractionated();
  bool fplexst = fourplex_fractionated_single_table_design.isFractionated();

  TEST_EQUAL(lf, false);
  TEST_EQUAL(lfst, false);
  TEST_EQUAL(lfstns, false);

  TEST_EQUAL(fplex, true);
  TEST_EQUAL(fplexst, true);
}
END_SECTION

START_SECTION((bool sameNrOfMSFilesPerFraction() const ))
{
  bool lf = labelfree_unfractionated_design.sameNrOfMSFilesPerFraction();
  bool lfst = labelfree_unfractionated_single_table_design.sameNrOfMSFilesPerFraction();
  bool lfstns = labelfree_unfractionated_single_table_no_sample_column.sameNrOfMSFilesPerFraction();
  bool fplex = fourplex_fractionated_design.sameNrOfMSFilesPerFraction();
  bool fplexst = fourplex_fractionated_single_table_design.sameNrOfMSFilesPerFraction();

  TEST_EQUAL(lf, true);
  TEST_EQUAL(lfst, true);
  TEST_EQUAL(lfstns, true);
  TEST_EQUAL(fplex, true);
  TEST_EQUAL(fplexst, true);
}
END_SECTION

START_SECTION((Size filterByBasenames(const set<std::string>&) keeps sample and file section synchronized))
{
  ExperimentalDesign design = ExperimentalDesignFile::load(
    OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_input_1.tsv"), false);

  std::set<std::string> keep{
    "JD_06232014_sample1-A.raw",
    "JD_06232014_sample1_B.raw"
  };
  design.filterByBasenames(keep);

  TEST_EQUAL(design.getMSFileSection().size(), 2);
  TEST_EQUAL(design.getNumberOfSamples(), 2);
}
END_SECTION

START_SECTION((static ExperimentalDesign fromConsensusMap(const ConsensusMap &c)))
{
  ConsensusXMLFile cfile;
  ConsensusMap cmap; 
  cfile.load(OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_input_3.consensusXML"), cmap);
  /* example consensusXML for TMT10Plex
  	<mapList count="10">
		<map id="0" name="C:/dev/OpenMS/src/tests/topp/TMTTenPlexMethod_test.mzML" label="tmt10plex_126" size="6">
			<UserParam type="string" name="channel_name" value="126"/>
			<UserParam type="int" name="channel_id" value="0"/>
			<UserParam type="string" name="channel_description" value=""/>
			<UserParam type="float" name="channel_center" value="126.127726"/>
		</map>
		<map id="1" name="C:/dev/OpenMS/src/tests/topp/TMTTenPlexMethod_test.mzML" label="tmt10plex_127N" size="6">
			<UserParam type="string" name="channel_name" value="127N"/>
			<UserParam type="int" name="channel_id" value="1"/>
			<UserParam type="string" name="channel_description" value=""/>
			<UserParam type="float" name="channel_center" value="127.124761"/>
		</map>
    ...
  */
  ExperimentalDesign ed_tmt10 = ExperimentalDesign::fromConsensusMap(cmap);
  TEST_EQUAL(ed_tmt10.getNumberOfLabels(), 10);
  TEST_EQUAL(ed_tmt10.getNumberOfMSFiles(), 1);
  TEST_EQUAL(ed_tmt10.getMSFileSection().at(0).label, 1); // "channel_id" + 1
  TEST_EQUAL(ed_tmt10.getMSFileSection().at(9).label, 10); // "channel_id" + 1
  TEST_EQUAL(ed_tmt10.getMSFileSection().at(0).fraction_group, 1); // only one fraction
  TEST_EQUAL(ed_tmt10.getMSFileSection().at(9).fraction_group, 1); // only one fraction
  TEST_EQUAL(ed_tmt10.getMSFileSection().at(0).fraction, 1); 
  TEST_EQUAL(ed_tmt10.getMSFileSection().at(9).fraction, 1);
  TEST_EQUAL(ed_tmt10.getMSFileSection().at(0).sample, 0); // default: sample from 0..n-1
  TEST_EQUAL(ed_tmt10.getMSFileSection().at(9).sample, 9);
  TEST_EQUAL(ed_tmt10.getMSFileSection().at(0).path, "C:/dev/OpenMS/src/tests/topp/TMTTenPlexMethod_test.mzML");
  TEST_EQUAL(ed_tmt10.getMSFileSection().at(1).path, "C:/dev/OpenMS/src/tests/topp/TMTTenPlexMethod_test.mzML");    

  /* example consensusXML for dimethyl labeling (FeatureFinderMultiplex) 
    <mapList count="2">
      <map id="0" name="/home/sachsenb/OpenMS/src/tests/topp/FeatureFinderMultiplex_1_input.mzML" label="Dimethyl0" size="2">
        <UserParam type="int" name="channel_id" value="0"/>
      </map>
      <map id="1" name="/home/sachsenb/OpenMS/src/tests/topp/FeatureFinderMultiplex_1_input.mzML" label="Dimethyl8" size="2">
        <UserParam type="int" name="channel_id" value="1"/>
      </map>
    </mapList>
  */
  cmap.clear();
  cfile.load(OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_input_4.consensusXML"), cmap);
  ExperimentalDesign ed_dimethyl = ExperimentalDesign::fromConsensusMap(cmap);
  TEST_EQUAL(ed_dimethyl.getNumberOfLabels(), 2);
  TEST_EQUAL(ed_dimethyl.getNumberOfMSFiles(), 1);
  TEST_EQUAL(ed_dimethyl.getMSFileSection().at(0).label, 1); // "channel_id" + 1
  TEST_EQUAL(ed_dimethyl.getMSFileSection().at(1).label, 2); // "channel_id" + 1
  TEST_EQUAL(ed_dimethyl.getMSFileSection().at(0).fraction_group, 1); // only one fraction
  TEST_EQUAL(ed_dimethyl.getMSFileSection().at(1).fraction_group, 1); // only one fraction
  TEST_EQUAL(ed_dimethyl.getMSFileSection().at(0).fraction, 1); 
  TEST_EQUAL(ed_dimethyl.getMSFileSection().at(1).fraction, 1);
  TEST_EQUAL(ed_dimethyl.getMSFileSection().at(0).sample, 0); // default: sample from 0..n-1
  TEST_EQUAL(ed_dimethyl.getMSFileSection().at(1).sample, 1);
  TEST_EQUAL(ed_dimethyl.getMSFileSection().at(0).path, "/home/sachsenb/OpenMS/src/tests/topp/FeatureFinderMultiplex_1_input.mzML");
  TEST_EQUAL(ed_dimethyl.getMSFileSection().at(1).path, "/home/sachsenb/OpenMS/src/tests/topp/FeatureFinderMultiplex_1_input.mzML");    

  /* example consensusXML for label-free (FeatureLinker*) 
    <mapList count="2">
      <map id="0" name="raw_file1.mzML" unique_id="8706403922746272921" label="" size="470">
      </map>
      <map id="1" name="raw_file2.mzML" unique_id="10253060449047408476" label="" size="423">
      </map>
    </mapList>
  */
  cmap.clear();
  cfile.load(OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_input_5.consensusXML"), cmap);
  ExperimentalDesign ed_labelfree = ExperimentalDesign::fromConsensusMap(cmap);
  TEST_EQUAL(ed_labelfree.getNumberOfLabels(), 1);
  TEST_EQUAL(ed_labelfree.getNumberOfMSFiles(), 2);
  TEST_EQUAL(ed_labelfree.getMSFileSection().at(0).label, 1); // "channel_id" + 1
  TEST_EQUAL(ed_labelfree.getMSFileSection().at(1).label, 1); // "channel_id" + 1
  TEST_EQUAL(ed_labelfree.getMSFileSection().at(0).fraction, 1); // only one fraction
  TEST_EQUAL(ed_labelfree.getMSFileSection().at(1).fraction, 1);
  TEST_EQUAL(ed_labelfree.getMSFileSection().at(0).fraction_group, 1); // each form a different group
  TEST_EQUAL(ed_labelfree.getMSFileSection().at(1).fraction_group, 2); 
  TEST_EQUAL(ed_labelfree.getMSFileSection().at(0).sample, 0); // default: sample from 0..n-1
  TEST_EQUAL(ed_labelfree.getMSFileSection().at(1).sample, 1);
  TEST_EQUAL(ed_labelfree.getMSFileSection().at(0).path, "raw_file1.mzML");
  TEST_EQUAL(ed_labelfree.getMSFileSection().at(1).path, "raw_file2.mzML");    
}
END_SECTION

START_SECTION((static ExperimentalDesign fromConsensusMap(const ConsensusMap &c) keeps sample_name entries in SampleSection))
{
  ConsensusMap cmap;
  cmap.setExperimentType("label-free");

  ConsensusMap::ColumnHeader h0;
  h0.filename = "file_a.mzML";
  h0.label = "label-free";
  h0.setMetaValue("sample_name", "Sample_A");
  cmap.getColumnHeaders()[0] = h0;

  ConsensusMap::ColumnHeader h1;
  h1.filename = "file_b.mzML";
  h1.label = "label-free";
  h1.setMetaValue("sample_name", "Sample_B");
  cmap.getColumnHeaders()[1] = h1;

  ExperimentalDesign design = ExperimentalDesign::fromConsensusMap(cmap);
  TEST_EQUAL(design.getSampleSection().hasSample("Sample_A"), true);
  TEST_EQUAL(design.getSampleSection().hasSample("Sample_B"), true);
}
END_SECTION

START_SECTION((ProteomicsLFQ subset output keeps design fraction_group assignments across fractions))
{
  const std::string design_file = OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_BSA_design_onetable_nonconsec.tsv");
  const std::string subset_output = OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_ProteomicsLFQ_1_subset_out.consensusXML");
  TEST_EQUAL(File::exists(design_file), true);
  TEST_EQUAL(File::exists(subset_output), true);

  ExperimentalDesign design = ExperimentalDesignFile::load(design_file, false);
  const auto pl2fg = design.getPathLabelToFractionGroupMapping(true);

  ConsensusMap cmap;
  ConsensusXMLFile().load(subset_output, cmap);

  for (const auto& [idx, col_header] : cmap.getColumnHeaders())
  {
    (void)idx;
    TEST_EQUAL(col_header.metaValueExists("fraction_group"), true);
    const unsigned observed_fg = static_cast<unsigned>(col_header.getMetaValue("fraction_group"));
    const unsigned label = col_header.getLabelAsUInt(cmap.getExperimentType());
    const unsigned expected_fg = pl2fg.at({File::basename(col_header.filename), label});
    TEST_EQUAL(observed_fg, expected_fg);
  }
}
END_SECTION

START_SECTION((static ExperimentalDesign fromFeatureMap(const FeatureMap &f)))
{
  FeatureXMLFile ffile;
  FeatureMap fmap;
  ffile.load(OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_input_6.featureXML"), fmap);
  ExperimentalDesign ed_labelfree = ExperimentalDesign::fromFeatureMap(fmap);
  TEST_EQUAL(ed_labelfree.getNumberOfLabels(), 1);
  TEST_EQUAL(ed_labelfree.getNumberOfMSFiles(), 1);
  TEST_EQUAL(ed_labelfree.getMSFileSection().at(0).label, 1); // "channel_id" + 1
  TEST_EQUAL(ed_labelfree.getMSFileSection().at(0).fraction_group, 1); // only one fraction
  TEST_EQUAL(ed_labelfree.getMSFileSection().at(0).fraction, 1); 
  TEST_EQUAL(ed_labelfree.getMSFileSection().at(0).sample, 0); // default: sample indices from 0..n-1
  TEST_EQUAL(ed_labelfree.getMSFileSection().at(0).path, "file://C:/raw_file1.mzML");
}
END_SECTION

START_SECTION((isValid_()))
{

  // missing fractions and wrong orders should work now
  std::string foo = OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_input_2_wrong.tsv");
  ExperimentalDesignFile::load(foo,false);
  std::string baz = OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_input_2_wrong_3.tsv");
  ExperimentalDesignFile::load(baz,false);

  // fraction groups still need to be consecutive
  std::string bar = OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_input_2_wrong_2.tsv");
  TEST_EXCEPTION(Exception::InvalidValue, ExperimentalDesignFile::load(bar,false));

}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
START_SECTION((Size annotateColumnHeaders(ConsensusMap& cmap) const))
{
  // The inverse of fromConsensusMap(). Two fraction files of one quantification unit must come
  // out sharing a fraction_group and differing by fraction -- that shared group is the only
  // thing distinguishing "two fractions of one sample" from "two independent runs", and a
  // consumer reading only the map has no other way to tell.
  ConsensusMap cmap;
  cmap.setExperimentType("labeled_MS2");
  for (UInt64 file = 0; file < 2; ++file)
  {
    for (UInt64 ch = 0; ch < 2; ++ch)
    {
      ConsensusMap::ColumnHeader h;
      h.filename = "/data/run" + std::to_string(file + 1) + ".mzML";
      h.label = "tmt10plex_12" + std::to_string(ch + 6);
      h.setMetaValue("channel_id", static_cast<int>(ch));   // getLabelAsUInt() -> channel_id + 1
      cmap.getColumnHeaders()[file * 2 + ch] = h;
    }
  }

  // Both files in fraction group 1, as fractions 1 and 2; one sample per channel.
  const std::string design_path = File::getTemporaryFile();
  {
    std::ofstream os(design_path.c_str());
    os << "Fraction_Group\tFraction\tSpectra_Filepath\tLabel\tSample\n";
    os << "1\t1\trun1.mzML\t1\t1\n";
    os << "1\t1\trun1.mzML\t2\t2\n";
    os << "1\t2\trun2.mzML\t1\t1\n";
    os << "1\t2\trun2.mzML\t2\t2\n";
  }
  ExperimentalDesign design = ExperimentalDesignFile::load(design_path, false);

  TEST_EQUAL(design.annotateColumnHeaders(cmap), 0)   // every header matched

  const auto& headers = cmap.getColumnHeaders();
  for (const auto& [mi, h] : headers)
  {
    TEST_TRUE(h.metaValueExists("fraction_group"))
    TEST_TRUE(h.metaValueExists("fraction"))
    TEST_EQUAL(static_cast<int>(h.getMetaValue("fraction_group")), 1)  // shared across both files
  }
  // ... and the fraction distinguishes the two files, channel by channel.
  TEST_EQUAL(static_cast<int>(headers.at(0).getMetaValue("fraction")), 1)  // run1, ch1
  TEST_EQUAL(static_cast<int>(headers.at(1).getMetaValue("fraction")), 1)  // run1, ch2
  TEST_EQUAL(static_cast<int>(headers.at(2).getMetaValue("fraction")), 2)  // run2, ch1
  TEST_EQUAL(static_cast<int>(headers.at(3).getMetaValue("fraction")), 2)  // run2, ch2

  // A header the design does not describe is counted, not silently skipped: silent
  // non-annotation is the failure mode this return value exists to expose.
  ConsensusMap::ColumnHeader stray;
  stray.filename = "/data/not_in_design.mzML";
  stray.label = "tmt10plex_126";
  stray.setMetaValue("channel_id", 0);
  cmap.getColumnHeaders()[99] = stray;
  TEST_EQUAL(design.annotateColumnHeaders(cmap), 1)
  TEST_FALSE(cmap.getColumnHeaders().at(99).metaValueExists("fraction_group"))
}
END_SECTION

START_SECTION(([EXTRA] annotateColumnHeaders - a design inferred from the map does not throw))
{
  // fromConsensusMap() builds sample ROWS but no "Sample" COLUMN, so getSampleName() throws
  // std::out_of_range on its column lookup. That is the design every caller gets when no
  // design file is supplied, so annotating with it must work: the sample name is optional
  // decoration, while fraction and fraction_group are what callers group on.
  ConsensusMap cmap;
  cmap.setExperimentType("labeled_MS2");
  for (UInt64 ch = 0; ch < 2; ++ch)
  {
    ConsensusMap::ColumnHeader h;
    h.filename = "/data/run1.mzML";
    h.label = "tmt10plex_12" + std::to_string(ch + 6);
    h.setMetaValue("channel_id", static_cast<int>(ch));
    cmap.getColumnHeaders()[ch] = h;
  }

  ExperimentalDesign inferred = ExperimentalDesign::fromConsensusMap(cmap);
  TEST_EQUAL(inferred.annotateColumnHeaders(cmap), 0)
  TEST_TRUE(cmap.getColumnHeaders().at(0).metaValueExists("fraction_group"))
  TEST_TRUE(cmap.getColumnHeaders().at(1).metaValueExists("fraction"))
}
END_SECTION

START_SECTION(([EXTRA] annotateColumnHeaders - headers that collapse onto one design row are refused))
{
  // ColumnHeader::getLabelAsUInt() returns 1 for any header with no 'channel_id' -- what the
  // Multiplex tools emit, carrying the channel in 'label' instead. Two such headers of one file
  // both resolve to label 1, so neither can be attributed to that design row. Annotating either
  // would stamp one channel's sample onto the other; they must be reported instead.
  ConsensusMap cmap;
  cmap.setExperimentType("labeled_MS2");
  for (UInt64 i = 0; i < 2; ++i)
  {
    ConsensusMap::ColumnHeader h;
    h.filename = "/data/run1.mzML";
    h.label = (i == 0) ? "light" : "heavy";   // no channel_id on purpose
    cmap.getColumnHeaders()[i] = h;
  }

  const std::string design_path = File::getTemporaryFile();
  {
    std::ofstream os(design_path.c_str());
    os << "Fraction_Group\tFraction\tSpectra_Filepath\tLabel\tSample\n";
    os << "1\t1\trun1.mzML\t1\t1\n";
    os << "1\t1\trun1.mzML\t2\t2\n";
  }
  ExperimentalDesign design = ExperimentalDesignFile::load(design_path, false);

  TEST_EQUAL(design.annotateColumnHeaders(cmap), 2)   // both refused, neither guessed at
  TEST_FALSE(cmap.getColumnHeaders().at(0).metaValueExists("sample_name"))
  TEST_FALSE(cmap.getColumnHeaders().at(1).metaValueExists("sample_name"))
}
END_SECTION

START_SECTION(([EXTRA] getSampleName works whether the section came from a file or from addSample))
{
  // A parsed design file fills the content table and a column called "Sample"; a section built
  // with addSample() - which is how fromConsensusMap() and fromIdentifications() build theirs -
  // fills only the name->row store and pushes an empty content row. Reading only the column threw
  // std::out_of_range for every inferred design.
  ExperimentalDesign::SampleSection added;
  added.addSample("BSA1");
  added.addSample("BSA2");
  TEST_STRING_EQUAL(added.getSampleName(0), "BSA1")
  TEST_STRING_EQUAL(added.getSampleName(1), "BSA2")
  TEST_EQUAL(added.getSampleRow("BSA2"), 1)

  // The content-backed form keeps answering from the column, as before.
  ExperimentalDesign::SampleSection from_file(
    {{"S_a", "control"}, {"S_b", "treated"}},
    {{"S_a", 0}, {"S_b", 1}},
    {{"Sample", 0}, {"MSstats_Condition", 1}});
  TEST_STRING_EQUAL(from_file.getSampleName(0), "S_a")
  TEST_STRING_EQUAL(from_file.getSampleName(1), "S_b")

  // Both report a missing row / unknown name as an OpenMS exception, not a bare std::out_of_range.
  TEST_EXCEPTION(Exception::ElementNotFound, added.getSampleName(7))
  TEST_EXCEPTION(Exception::ElementNotFound, added.getSampleRow("nope"))
}
END_SECTION

START_SECTION(([EXTRA] an inferred design annotates sample_name onto consensus map headers))
{
  // Consequence of the above: fromConsensusMap() builds its sample section with addSample(), so
  // annotateColumnHeaders() can now supply the name it previously had to drop.
  ConsensusMap cmap;
  cmap.setExperimentType("label-free");
  ConsensusMap::ColumnHeaders headers;
  for (Size i = 0; i < 2; ++i)
  {
    ConsensusMap::ColumnHeader h;
    h.filename = "/data/run_" + StringUtils::toStr(i) + ".mzML";
    h.label = "label-free";
    headers[i] = h;
  }
  cmap.setColumnHeaders(headers);

  const ExperimentalDesign inferred = ExperimentalDesign::fromConsensusMap(cmap);
  TEST_EQUAL(inferred.annotateColumnHeaders(cmap), 0)
  TEST_TRUE(cmap.getColumnHeaders().at(0).metaValueExists("sample_name"))
  TEST_TRUE(cmap.getColumnHeaders().at(1).metaValueExists("sample_name"))
}
END_SECTION

START_SECTION(([EXTRA] sample names that are not stringified row indices must not break the derived mappings))
{
  // A Sample column holds a NAME. sdrf-pipelines writes 1-based numeric names ("1", "2"), and a
  // hand-written design may use anything ("BSA1"). MSFileSectionEntry::sample is the zero-based row
  // index those names intern to, so a name and its index coincide only by accident - and only the
  // designs that OpenMS infers itself (which set sample_name = toStr(sample)) make them coincide.
  auto designWithSampleNames = [](const std::string& name_a, const std::string& name_b)
  {
    ExperimentalDesign::MSFileSection fs;
    for (Size i = 0; i < 2; ++i)
    {
      ExperimentalDesign::MSFileSectionEntry e;
      e.fraction_group = static_cast<unsigned>(i) + 1; // 1-based id
      e.fraction = 1;
      e.label = 1;
      e.path = "/data/run_" + StringUtils::toStr(i) + ".mzML";
      e.sample = static_cast<unsigned>(i);             // 0-based index
      e.sample_name = (i == 0) ? name_a : name_b;
      fs.push_back(e);
    }
    ExperimentalDesign::SampleSection ss;
    ss.addSample(name_a);
    ss.addSample(name_b);
    return ExperimentalDesign(fs, ss);
  };

  // The sdrf-pipelines shape: names are 1-based, indices 0-based, so they never coincide.
  const ExperimentalDesign sdrf_like = designWithSampleNames("1", "2");
  // A hand-written design: names are not numeric at all.
  const ExperimentalDesign named = designWithSampleNames("BSA1", "BSA2");

  for (const auto& design : {sdrf_like, named})
  {
    const std::set<std::string> names = design.getSampleSection().getSamples();

    // Both Sample* mappings are keyed by sample NAME, and cover every sample exactly once.
    const auto s2pf = design.getSampleToPrefractionationMapping();
    const auto s2c = design.getSampleToConditionMapping();
    TEST_EQUAL(s2pf.size(), design.getNumberOfSamples())
    TEST_EQUAL(s2c.size(), design.getNumberOfSamples())
    for (const auto& name : names)
    {
      TEST_TRUE(s2pf.find(name) != s2pf.end())
      TEST_TRUE(s2c.find(name) != s2c.end())
    }

    // Resolving them per (path, label) must go through the name, not the stringified index.
    // Both of these threw a bare std::out_of_range ("map::at") for any design whose sample names
    // are not "0", "1", ... - which is every design a user or sdrf-pipelines writes.
    const auto p2pf = design.getPathLabelToPrefractionationMapping(false);
    const auto p2c = design.getPathLabelToConditionMapping(false);
    TEST_EQUAL(p2pf.size(), 2)
    TEST_EQUAL(p2c.size(), 2)
  }
}
END_SECTION

START_SECTION(([EXTRA] fromIdentifications numbers fraction groups from 1, samples from 0))
{
  // 'sample' is a zero-based index into the sample section, a fraction group is a 1-based id.
  // One counter used to feed both, so inferred designs were rejected by isValid_() - the class
  // produced objects its own validator refuses - and leaked a fraction group of 0 downstream.
  ProteinIdentification protein_id;
  protein_id.setPrimaryMSRunPath({"a.mzML", "b.mzML", "c.mzML"});
  const ExperimentalDesign design = ExperimentalDesign::fromIdentifications({protein_id});

  const ExperimentalDesign::MSFileSection& rows = design.getMSFileSection();
  TEST_EQUAL(rows.size(), 3)
  for (Size i = 0; i < rows.size(); ++i)
  {
    TEST_EQUAL(rows[i].sample, i)             // zero-based index
    TEST_EQUAL(rows[i].fraction_group, i + 1) // 1-based id
  }

  // The decisive invariant: an inferred design must satisfy the same validator as a loaded one.
  // The two-argument constructor runs isValid_(), which requires the fraction-group set to be
  // consecutive and start at 1; before the fix this threw InvalidValue.
  ExperimentalDesign revalidated(rows, design.getSampleSection());
  TEST_EQUAL(revalidated.getNumberOfFractionGroups(), 3)
}
END_SECTION

END_TEST
