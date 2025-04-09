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
#include <OpenMS/FORMAT/ExperimentalDesignFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
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

START_SECTION((void setSampleSection(const ExperimentalDesign::SampleSection &sample_section)))
{
  ExperimentalDesign labelfree_unfractionated_design2 = labelfree_unfractionated_design;
  ExperimentalDesign::SampleSection fs;
  labelfree_unfractionated_design2.setSampleSection(fs);
}
END_SECTION

START_SECTION((std::map<unsigned int, std::vector<String> > getFractionToMSFilesMapping() const ))
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
      if (i.first.first.hasSubstring("TR2")) { file = 2; }
      TEST_EQUAL(i.second, file); 
    }
  }
}
END_SECTION

START_SECTION((std::set< String > ExperimentalDesign::SampleSection::getFactors() const))
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


START_SECTION((String SampleSection::getFactorValue(const unsigned sample, const String &factor) const))
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
      String f1 = lss_st.getFactorValue(sample, factor);
      String f2 = lss_tt.getFactorValue(sample, factor);
      String f3 = lss_stns.getFactorValue(sample, factor);
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
      String f1 = fss_st.getFactorValue(sample, factor);
      String f2 = fss_tt.getFactorValue(sample, factor);
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
  String foo = OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_input_2_wrong.tsv");
  ExperimentalDesignFile::load(foo,false);
  String baz = OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_input_2_wrong_3.tsv");
  ExperimentalDesignFile::load(baz,false);

  // fraction groups still need to be consecutive
  String bar = OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_input_2_wrong_2.tsv");
  TEST_EXCEPTION(Exception::InvalidValue, ExperimentalDesignFile::load(bar,false));

}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
