// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Erhan Kenar, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/AccurateMassSearchEngine.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTabMFile.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(AccurateMassSearchEngine, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

AccurateMassSearchEngine* ptr = nullptr;
AccurateMassSearchEngine* null_ptr = nullptr;
START_SECTION(AccurateMassSearchEngine())
{
    ptr = new AccurateMassSearchEngine();
    TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(virtual ~AccurateMassSearchEngine())
{
    delete ptr;
}
END_SECTION

START_SECTION([EXTRA]AdductInfo)
{
  EmpiricalFormula ef_empty;
  // make sure an empty formula has no weight (we rely on that in AdductInfo's getMZ() and getNeutralMass()
  TEST_EQUAL(ef_empty.getMonoWeight(), 0)

  // now we test if converting from neutral mass to m/z and back recovers the input value using different adducts
  {
  // testing M;-2  // intrinsic doubly negative charge
    AdductInfo ai("TEST_INTRINSIC", ef_empty, -2, 1);
    double neutral_mass=1000; // some mass...
    double mz = ai.getMZ(neutral_mass);
    double neutral_mass_recon = ai.getNeutralMass(mz);
    TEST_REAL_SIMILAR(neutral_mass, neutral_mass_recon);
  }
  { // testing M+Na+H;+2
    EmpiricalFormula simpleAdduct("HNa");
    AdductInfo ai("TEST_WITHADDUCT", simpleAdduct, 2, 1);
    double neutral_mass=1000; // some mass...
    double mz = ai.getMZ(neutral_mass);
    double neutral_mass_recon = ai.getNeutralMass(mz);
    TEST_REAL_SIMILAR(neutral_mass, neutral_mass_recon);
  }

}
END_SECTION

Param ams_param;
ams_param.setValue("db:mapping", std::vector<std::string>{OPENMS_GET_TEST_DATA_PATH("reducedHMDBMapping.tsv")});
ams_param.setValue("db:struct", std::vector<std::string>{OPENMS_GET_TEST_DATA_PATH("reducedHMDB2StructMapping.tsv")});
ams_param.setValue("keep_unidentified_masses", "true");
ams_param.setValue("mzTab:exportIsotopeIntensities", "true");
AccurateMassSearchEngine ams;
ams.setParameters(ams_param);

START_SECTION(void init())
  NOT_TESTABLE // tested below
END_SECTION

START_SECTION((void queryByMZ(const double& observed_mz, const Int& observed_charge, const String& ion_mode, std::vector<AccurateMassSearchResult>& results) const))
{
  std::vector<AccurateMassSearchResult> hmdb_results_pos;

  // test 'ams' not initialized
  TEST_EXCEPTION(Exception::IllegalArgument, ams.queryByMZ(1234, 1, "positive", hmdb_results_pos));
  ams.init();

  // test invalid scan polarity
  TEST_EXCEPTION(Exception::InvalidParameter, ams.queryByMZ(1234, 1, "this_is_an_invalid_ionmode", hmdb_results_pos));

  // test the actual query
  {
    Param ams_param_tmp = ams_param;
    ams_param_tmp.setValue("mass_error_value", 17.0);
    ams.setParameters(ams_param_tmp);
    ams.init();
    // -- positive mode
    // expected hit: C17H11N5 with neutral mass ~285.101445377
    double m = EmpiricalFormula("C17H11N5").getMonoWeight(); 
    double mz = m / 1 + EmpiricalFormula("Na").getMonoWeight() - Constants::ELECTRON_MASS_U; // assume M+Na;+1 as charge
    std::cout << "mz query mass:" << mz << "\n\n";
    // we'll get some other hits as well...
    String id_list_pos[] = {"C10H17N3O6S", "C15H16O7", "C14H14N2OS2", "C16H15NO4",
                            "C17H11N5" /* this one we want! */,
                            "C10H14NO6P", "C14H12O4", "C7H6O2"};
                         //{"C10H17N3O6S", "C15H16O7", "C14H14N2OS2", "C16H15NO4", "C17H11N5", "C10H14NO6P", "C14H12O4", "C7H6O2"};

                         // 290.05475446	C14H14N2OS2	HMDB:HMDB38641 missing

    Size id_list_pos_length(sizeof(id_list_pos)/sizeof(id_list_pos[0]));
    ams.queryByMZ(mz, 1, "positive", hmdb_results_pos);
    ams.setParameters(ams_param); // reset to default 5ppm
    ams.init();
    TEST_EQUAL(hmdb_results_pos.size(), id_list_pos_length)
    ABORT_IF(hmdb_results_pos.size() != id_list_pos_length)
    for (Size i = 0; i < id_list_pos_length; ++i)
    {
      TEST_STRING_EQUAL(hmdb_results_pos[i].getFormulaString(), id_list_pos[i])
      std::cout << hmdb_results_pos[i] << std::endl;
    }
    TEST_EQUAL(hmdb_results_pos[4].getFormulaString(), "C17H11N5"); // correct hit?
    TEST_REAL_SIMILAR(hmdb_results_pos[4].getQueryMass(), m); // was the mass correctly reconstructed internally?
    TEST_REAL_SIMILAR(abs(hmdb_results_pos[4].getMZErrorPPM()), 0.0); // ppm error within float precision? 

  }
  
  // -- negative mode 
  // expected hit: C17H20N2S with neutral mass ~284.13472	
  {
    std::vector<AccurateMassSearchResult> hmdb_results_neg;
    double m = EmpiricalFormula("C17H20N2S").getMonoWeight(); 
    double mz = m / 3 - Constants::PROTON_MASS_U; // assume M-3H;-3 as charge
    // manual check:
    // double mass_recovered = mz * 3 - EmpiricalFormula("H-3").getMonoWeight() - Constants::ELECTRON_MASS_U*3;
    ams.queryByMZ(mz, 3, "negative", hmdb_results_neg);
    ABORT_IF(hmdb_results_neg.size() != 1)
    std::cout << hmdb_results_neg[0] << std::endl;
    TEST_EQUAL(hmdb_results_neg[0].getFormulaString(), "C17H20N2S"); // correct hit?
    TEST_REAL_SIMILAR(hmdb_results_neg[0].getQueryMass(), m); // was the mass correctly reconstructed internally?
    TEST_EQUAL(abs(hmdb_results_neg[0].getMZErrorPPM()) < 0.0002, true); // ppm error within float precision? .. should be ~0.0001576..
  }
}
END_SECTION

AccurateMassSearchEngine ams_feat_test;
ams_feat_test.setParameters(ams_param);
ams_feat_test.init();
String feat_query_pos[] = {"C23H45NO4", "C20H37NO3", "C22H41NO"};

START_SECTION((void queryByFeature(const Feature& feature, const Size& feature_index, const String& ion_mode, std::vector<AccurateMassSearchResult>& results) const))
{
  Feature test_feat;
  test_feat.setRT(300.0);
  test_feat.setMZ(399.33486);
  test_feat.setIntensity(100.0);
  test_feat.setMetaValue(Constants::UserParam::NUM_OF_MASSTRACES, 3);
  test_feat.setCharge(1.0);

  vector<double> masstrace_intenstiy = {100.0, 26.1, 4.0};
  test_feat.setMetaValue("masstrace_intensity", masstrace_intenstiy);

  //test_feat.setMetaValue("masstrace_intensity_0", 100.0);
  //test_feat.setMetaValue("masstrace_intensity_1", 26.1);
  //test_feat.setMetaValue("masstrace_intensity_2", 4.0);

  std::vector<AccurateMassSearchResult> results;
  
  // invalid scan_polarity
  TEST_EXCEPTION(Exception::InvalidParameter, ams_feat_test.queryByFeature(test_feat, 0, "invalid_scan_polatority", results));
  
  // actual test
  ams_feat_test.queryByFeature(test_feat, 0, "positive", results);

  TEST_EQUAL(results.size(), 3)

  for (Size i = 0; i < results.size(); ++i)
  {
    TEST_REAL_SIMILAR(results[i].getObservedRT(), 300.0)
    TEST_REAL_SIMILAR(results[i].getObservedIntensity(), 100.0)
  }

  Size feat_query_size(sizeof(feat_query_pos)/sizeof(feat_query_pos[0]));

  ABORT_IF(results.size() != feat_query_size)
  for (Size i = 0; i < feat_query_size; ++i)
  {
    TEST_STRING_EQUAL(results[i].getFormulaString(), feat_query_pos[i])
  }
}
END_SECTION


START_SECTION((void queryByConsensusFeature(const ConsensusFeature& cfeat, const Size& cf_index, const Size& number_of_maps, const String& ion_mode, std::vector<AccurateMassSearchResult>& results) const))
{
  ConsensusFeature cons_feat;
  cons_feat.setRT(300.0);
  cons_feat.setMZ(399.33486);
  cons_feat.setIntensity(100.0);
  cons_feat.setCharge(1.0);

  FeatureHandle fh1, fh2, fh3;
  fh1.setRT(300.0);
  fh1.setMZ(399.33485);
  fh1.setIntensity(100.0);
  fh1.setCharge(1.0);
  fh1.setMapIndex(0);

  fh2.setRT(310.0);
  fh2.setMZ(399.33486);
  fh2.setIntensity(300.0);
  fh2.setCharge(1.0);
  fh2.setMapIndex(1);

  fh3.setRT(290.0);
  fh3.setMZ(399.33487);
  fh3.setIntensity(500.0);
  fh3.setCharge(1.0);
  fh3.setMapIndex(2);

  cons_feat.insert(fh1);
  cons_feat.insert(fh2);
  cons_feat.insert(fh3);
  cons_feat.computeConsensus();
  
  std::vector<AccurateMassSearchResult> results;

  TEST_EXCEPTION(Exception::InvalidParameter, ams_feat_test.queryByConsensusFeature(cons_feat, 0, 3, "blabla", results)); // invalid scan_polarity
  ams_feat_test.queryByConsensusFeature(cons_feat, 0, 3, "positive", results);

  TEST_EQUAL(results.size(), 3)

  for (Size i = 0; i < results.size(); ++i)
  {
      TEST_REAL_SIMILAR(results[i].getObservedRT(), 300.0)
      TEST_REAL_SIMILAR(results[i].getObservedIntensity(), 0.0)
  }

  // std::cout << cons_feat.getMZ() << " " << results.size() << std::endl;

  for (Size i = 0; i < results.size(); ++i)
  {
    std::vector<double> indiv_ints = results[i].getIndividualIntensities();
    TEST_EQUAL(indiv_ints.size(), 3)

    ABORT_IF(indiv_ints.size() != 3)
    TEST_REAL_SIMILAR(indiv_ints[0], fh1.getIntensity());
    TEST_REAL_SIMILAR(indiv_ints[1], fh2.getIntensity());
    TEST_REAL_SIMILAR(indiv_ints[2], fh3.getIntensity());
  }

  Size feat_query_size(sizeof(feat_query_pos)/sizeof(feat_query_pos[0]));

  ABORT_IF(results.size() != feat_query_size)
  for (Size i = 0; i < feat_query_size; ++i)
  {
    TEST_STRING_EQUAL(results[i].getFormulaString(), feat_query_pos[i])
  }
}
END_SECTION

FuzzyStringComparator fsc;
// fsc.setAcceptableAbsolute((3.04011223650013 - 3.04011223637974)*1.1); // 1.3242891228060217e-10
// also Linux may give slightly different results depending on optimization level (O0 vs O1) 
// note that the default value for TEST_REAL_SIMILAR is 1e-5, see ./source/CONCEPT/ClassTest.cpp
fsc.setAcceptableAbsolute(1e-8);
StringList sl;
sl.push_back("xml-stylesheet");
sl.push_back("IdentificationRun id=\"PI_0\" date=");
sl.push_back("software[1]");
sl.push_back("database[1]-uri");
fsc.setWhitelist(sl);

START_SECTION((void run(FeatureMap&, MzTab&) const))
{
  FeatureMap exp_fm;
  FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("AccurateMassSearchEngine_input1.featureXML"), exp_fm);
  {
    MzTab test_mztab;
    ams_feat_test.run(exp_fm, test_mztab);

    // test annotation of input
    String tmp_file;
    NEW_TMP_FILE(tmp_file);
    FeatureXMLFile ff;
    ff.store(tmp_file, exp_fm);
    TEST_EQUAL(fsc.compareFiles(tmp_file, OPENMS_GET_TEST_DATA_PATH("AccurateMassSearchEngine_output1.featureXML")), true);

    String tmp_mztab_file;
    NEW_TMP_FILE(tmp_mztab_file);
    MzTabFile().store(tmp_mztab_file, test_mztab);
    TEST_EQUAL(fsc.compareFiles(tmp_mztab_file, OPENMS_GET_TEST_DATA_PATH("AccurateMassSearchEngine_output1_featureXML.mzTab")), true);
    
    // test use of adduct information
    Param ams_param_tmp = ams_param;
    ams_param_tmp.setValue("use_feature_adducts", "true");
      
    AccurateMassSearchEngine ams_feat_test2;
    ams_feat_test2.setParameters(ams_param_tmp);
    ams_feat_test2.init();

    FeatureMap exp_fm2;
    FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("AccurateMassSearchEngine_input1.featureXML"), exp_fm2);
    MzTab test_mztab2;
    ams_feat_test2.run(exp_fm2, test_mztab2);

    String tmp_mztab_file2;
    NEW_TMP_FILE(tmp_mztab_file2);
    MzTabFile().store(tmp_mztab_file2, test_mztab2);
    TEST_EQUAL(fsc.compareFiles(tmp_mztab_file2, OPENMS_GET_TEST_DATA_PATH("AccurateMassSearchEngine_output2_featureXML.mzTab")), true);
  }
}
END_SECTION

START_SECTION((void run(FeatureMap&, MzTabM&) const))
{
  FeatureMap exp_fm;
  FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("AccurateMassSearchEngine_input1.featureXML"), exp_fm);
  {
    MzTabM test_mztabm;
    ams_feat_test.run(exp_fm, test_mztabm);

    String tmp_mztabm_file;
    NEW_TMP_FILE(tmp_mztabm_file);
    MzTabMFile().store(tmp_mztabm_file, test_mztabm);
    TEST_EQUAL(fsc.compareFiles(tmp_mztabm_file, OPENMS_GET_TEST_DATA_PATH("AccurateMassSearchEngine_output1_mztabm_featureXML.mzTab")), true);

    // test use of adduct information
    Param ams_param_tmp = ams_param;
    ams_param_tmp.setValue("use_feature_adducts", "true");

    AccurateMassSearchEngine ams_feat_test2;
    ams_feat_test2.setParameters(ams_param_tmp);
    ams_feat_test2.init();

    FeatureMap exp_fm2;
    FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("AccurateMassSearchEngine_input1.featureXML"), exp_fm2);
    MzTabM test_mztabm2;
    ams_feat_test2.run(exp_fm2, test_mztabm2);

    String tmp_mztabm_file2;
    NEW_TMP_FILE(tmp_mztabm_file2);
    MzTabMFile().store(tmp_mztabm_file2, test_mztabm2);
    TEST_EQUAL(fsc.compareFiles(tmp_mztabm_file2, OPENMS_GET_TEST_DATA_PATH("AccurateMassSearchEngine_output2_mztabm_featureXML.mzTab")), true);
  }
}
END_SECTION

START_SECTION((void run(ConsensusMap&, MzTab&) const))
  ConsensusMap exp_cm;
  ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("AccurateMassSearchEngine_input1.consensusXML"), exp_cm);
  MzTab test_mztab2;
  ams_feat_test.run(exp_cm, test_mztab2);

  // test annotation of input
  String tmp_file;
  NEW_TMP_FILE(tmp_file);
  ConsensusXMLFile ff;
  ff.store(tmp_file, exp_cm);
  TEST_EQUAL(fsc.compareFiles(tmp_file, OPENMS_GET_TEST_DATA_PATH("AccurateMassSearchEngine_output1.consensusXML")), true);

  String tmp_mztab_file;
  NEW_TMP_FILE(tmp_mztab_file);
  MzTabFile().store(tmp_mztab_file, test_mztab2);
  TEST_EQUAL(fsc.compareFiles(tmp_mztab_file, OPENMS_GET_TEST_DATA_PATH("AccurateMassSearchEngine_output1_consensusXML.mzTab")), true);
END_SECTION

START_SECTION([EXTRA] template <typename MAPTYPE> void resolveAutoMode_(const MAPTYPE& map))
  FeatureMap exp_fm;
  FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("AccurateMassSearchEngine_input1.featureXML"), exp_fm);
  FeatureMap fm_p = exp_fm;
  AccurateMassSearchEngine ams;
  MzTab mzt;
  Param p;
  p.setValue("ionization_mode","auto");
  p.setValue("db:mapping", std::vector<std::string>{OPENMS_GET_TEST_DATA_PATH("reducedHMDBMapping.tsv")});
  p.setValue("db:struct", std::vector<std::string>{OPENMS_GET_TEST_DATA_PATH("reducedHMDB2StructMapping.tsv")});
  ams.setParameters(p);
  ams.init();

  TEST_EXCEPTION(Exception::InvalidParameter, ams.run(fm_p, mzt)); // 'fm_p' has no scan_polarity meta value
  for (auto& f : fm_p)
  {
    f.setMetaValue("scan_polarity", "something;somethingelse");
  }
  TEST_EXCEPTION(Exception::InvalidParameter, ams.run(fm_p, mzt)); // 'fm_p' scan_polarity meta value wrong

  for (auto& f : fm_p)
  {
    f.setMetaValue("scan_polarity", "positive");
  }
  ams.run(fm_p, mzt);

  for (auto& f : fm_p)
  {
    f.setMetaValue("scan_polarity", "negative");
  }
  ams.run(fm_p, mzt);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
