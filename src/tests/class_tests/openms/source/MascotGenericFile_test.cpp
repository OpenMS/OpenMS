// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/MascotGenericFile.h>
#include <sstream>
///////////////////////////

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

using namespace OpenMS;
using namespace std;

START_TEST(MascotGenericFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MascotGenericFile* ptr = nullptr;
MascotGenericFile* nullPointer = nullptr;
START_SECTION(MascotGenericFile())
{
  ptr = new MascotGenericFile();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~MascotGenericFile())
{
  delete ptr;
}
END_SECTION

ptr = new MascotGenericFile();

START_SECTION((template < typename MapType > void load(const String &filename, MapType &exp)))
{
  PeakMap exp;
  ptr->load(OPENMS_GET_TEST_DATA_PATH("MascotInfile_test.mascot_in"), exp);
  TEST_EQUAL(exp.size(), 1)

  TEST_EQUAL(exp.begin()->size(), 9)
}
END_SECTION

START_SECTION((void store(std::ostream &os, const String &filename, const PeakMap &experiment, bool compact = false)))
{
  PeakMap exp;
  ptr->load(OPENMS_GET_TEST_DATA_PATH("MascotInfile_test.mascot_in"), exp);

  // handling of modifications:
  Param params = ptr->getParameters();
  params.setValue("fixed_modifications", std::vector<std::string>{"Carbamidomethyl (C)","Phospho (S)"});
  params.setValue("variable_modifications", std::vector<std::string>{"Oxidation (M)","Deamidated (N)","Deamidated (Q)"});
  ptr->setParameters(params);

  stringstream ss;
  ptr->store(ss, "test", exp);

  vector<String> strings;
  strings.push_back("BEGIN IONS\n"
                    "TITLE=Testtitle_index=0\n" // different from input!
                    "PEPMASS=1998.0\n"
                    "RTINSECONDS=25.379000000000001\n"
                    "SCANS=0");
  strings.push_back("1.0 1.0\n"
                    "2.0 4.0\n"
                    "3.0 9.0\n"
                    "4.0 16.0\n"
                    "5.0 25.0\n"
                    "6.0 36.0\n"
                    "7.0 49.0\n"
                    "8.0 64.0\n"
                    "9.0 81.0\n"
                    "END IONS\n");
  strings.push_back("MODS=Carbamidomethyl (C)\n");
  strings.push_back("MODS=Phospho (ST)\n");
  strings.push_back("IT_MODS=Deamidated (NQ)");
  strings.push_back("IT_MODS=Oxidation (M)");

  String mgf_file(ss.str());
  for (Size i = 0; i < strings.size(); ++i)
  {
    TEST_EQUAL(mgf_file.hasSubstring(strings[i]), true)
  }

  // test of making default TITLE
  exp[0].removeMetaValue("TITLE");
  stringstream ss2;
  ptr->store(ss2, "test", exp);
  vector<String> strings2;
  strings2.push_back("BEGIN IONS\n"
                    "TITLE=1998.0_25.379000000000001_index=0_test\n" // different from input!
                    "PEPMASS=1998.0\n"
                    "RTINSECONDS=25.379000000000001\n"
                    "SCANS=0");
  strings2.push_back("1.0 1.0\n"
                    "2.0 4.0\n"
                    "3.0 9.0\n"
                    "4.0 16.0\n"
                    "5.0 25.0\n"
                    "6.0 36.0\n"
                    "7.0 49.0\n"
                    "8.0 64.0\n"
                    "9.0 81.0\n"
                    "END IONS\n");
  strings2.push_back("MODS=Carbamidomethyl (C)\n");
  strings2.push_back("MODS=Phospho (ST)\n");
  strings2.push_back("IT_MODS=Deamidated (NQ)");
  strings2.push_back("IT_MODS=Oxidation (M)");
  String mgf_file2(ss2.str());
  for (Size i = 0; i < strings2.size(); ++i)
  {
    TEST_EQUAL(mgf_file2.hasSubstring(strings2[i]), true)
  }

  ptr->setParameters(ptr->getDefaults()); // reset parameters

  // test compact format:
  MSSpectrum spec;
  spec.setNativeID("index=250");
  spec.setMSLevel(2);
  spec.setRT(234.5678901);
  Precursor prec;
  prec.setMZ(901.2345678);
  spec.getPrecursors().push_back(prec);
  Peak1D peak;
  peak.setMZ(567.8901234);
  peak.setIntensity(0.0);
  spec.push_back(peak); // intensity zero -> not present in output
  peak.setMZ(890.1234567);
  peak.setIntensity(2345.678901);
  spec.push_back(peak);
  exp.clear(true);
  exp.addSpectrum(spec);

  ss.str("");
  ptr->store(ss, "test", exp, true);
  mgf_file = ss.str();
  String content = ("BEGIN IONS\n"
                    "TITLE=901.23457_234.568_index=250_test\n"
                    "PEPMASS=901.23457\n"
                    "RTINSECONDS=234.568\n"
                    "SCANS=250\n"
                    "890.12346 2345.679\n"
                    "END IONS");
  TEST_EQUAL(mgf_file.hasSubstring(content), true);
}
END_SECTION

START_SECTION((void store(const String &filename, const PeakMap &experiment, bool compact = false)))
{
  String tmp_name("MascotGenericFile_1.tmp");
  NEW_TMP_FILE(tmp_name)
  PeakMap exp;
  ptr->load(OPENMS_GET_TEST_DATA_PATH("MascotInfile_test.mascot_in"), exp);


  ptr->store(tmp_name, exp);

  PeakMap exp2;
  ptr->load(tmp_name, exp2);
  TEST_EQUAL(exp.size() == exp2.size(), true)
  TEST_EQUAL(exp.begin()->size() == exp2.begin()->size(), true)
  TEST_REAL_SIMILAR(exp.begin()->getRT(), exp2.begin()->getRT())
  TEST_REAL_SIMILAR(exp.begin()->getPrecursors().begin()->getMZ(), exp2.begin()->getPrecursors().begin()->getMZ())
}
END_SECTION

START_SECTION((COMPOUND_NAME to Metabolite_Name mapping))
{
  // Test that COMPOUND_NAME in MGF is correctly mapped to Constants::UserParam::MSM_METABOLITE_NAME
  String mgf_content = "BEGIN IONS\n"
                       "TITLE=Test spectrum\n"
                       "PEPMASS=500.0\n"
                       "CHARGE=1\n"
                       "COMPOUND_NAME=Caffeine\n"
                       "100.0 1000.0\n"
                       "200.0 2000.0\n"
                       "END IONS\n";

  stringstream mgf_stream(mgf_content);
  
  PeakMap exp;
  MascotGenericFile mgf_file;
  
  // Create a temporary file to test the loading functionality
  String tmp_name("test_compound_name.mgf");
  NEW_TMP_FILE(tmp_name)
  
  // Write MGF content to temporary file
  std::ofstream ofs(tmp_name.c_str());
  ofs << mgf_content;
  ofs.close();
  
  // Load the MGF file
  mgf_file.load(tmp_name, exp);
  
  // Test that we have one spectrum
  TEST_EQUAL(exp.size(), 1)
  
  // Test that the spectrum has the correct metabolite name metadata
  TEST_EQUAL(exp[0].metaValueExists(Constants::UserParam::MSM_METABOLITE_NAME), true)
  TEST_EQUAL(String(exp[0].getMetaValue(Constants::UserParam::MSM_METABOLITE_NAME)), "Caffeine")
  
  // Test that other expected properties are also parsed correctly
  TEST_EQUAL(exp[0].size(), 2) // Two peaks
  TEST_REAL_SIMILAR(exp[0].getPrecursors()[0].getMZ(), 500.0)
  TEST_EQUAL(exp[0].getPrecursors()[0].getCharge(), 1)
}
END_SECTION

START_SECTION((GNPS MGF file - 3-Des-Microcystein_LR))
{
  // Test loading a real GNPS library spectrum with COMPOUND_NAME and SPECTRUMID metadata
  PeakMap exp;
  MascotGenericFile mgf_file;
  
  // Load the GNPS MGF file
  mgf_file.load(OPENMS_GET_TEST_DATA_PATH("MascotGenericFile_GNPS.mgf"), exp);
  
  // Test that we have one spectrum
  TEST_EQUAL(exp.size(), 1)
  
  // Test that SPECTRUMID was correctly parsed and stored as GNPS_Spectrum_ID
  TEST_EQUAL(exp[0].metaValueExists("GNPS_Spectrum_ID"), true)
  TEST_EQUAL(String(exp[0].getMetaValue("GNPS_Spectrum_ID")), "CCMSLIB00000001547")
  
  // Test that COMPOUND_NAME was correctly mapped to MSM_METABOLITE_NAME
  TEST_EQUAL(exp[0].metaValueExists(Constants::UserParam::MSM_METABOLITE_NAME), true)
  TEST_EQUAL(String(exp[0].getMetaValue(Constants::UserParam::MSM_METABOLITE_NAME)), "3-Des-Microcystein_LR")
  
  // Test precursor m/z
  TEST_REAL_SIMILAR(exp[0].getPrecursors()[0].getMZ(), 981.54)
  
  // Test charge state
  TEST_EQUAL(exp[0].getPrecursors()[0].getCharge(), 1)
  
  // Test that we have the expected number of peaks (43 peaks in the file)
  TEST_EQUAL(exp[0].size(), 43)
  
  // Test the base peak (m/z 599.352783 with intensity 764523.0)
  bool found_base_peak = false;
  for (Size i = 0; i < exp[0].size(); ++i)
  {
    if (std::abs(exp[0][i].getMZ() - 599.352783) < 0.001)
    {
      found_base_peak = true;
      TEST_REAL_SIMILAR(exp[0][i].getIntensity(), 764523.0)
      break;
    }
  }
  TEST_EQUAL(found_base_peak, true)
}
END_SECTION

START_SECTION((SEQ sequence query field - single and multiple))
{
  // Single SEQ line: parsed, stored as StringList (always), round-tripped on write.
  {
    String mgf_content = "BEGIN IONS\n"
                         "TITLE=seq_single\n"
                         "PEPMASS=500.0\n"
                         "CHARGE=2+\n"
                         "SEQ=PEPTIDER\n"
                         "100.0 1000.0\n"
                         "200.0 2000.0\n"
                         "END IONS\n";

    String tmp_in("MascotGenericFile_SEQ_single_in.mgf");
    NEW_TMP_FILE(tmp_in)
    std::ofstream ofs(tmp_in.c_str());
    ofs << mgf_content;
    ofs.close();

    PeakMap exp;
    MascotGenericFile mgf_file;
    mgf_file.load(tmp_in, exp);

    TEST_EQUAL(exp.size(), 1)
    TEST_TRUE(exp[0].metaValueExists("SEQ"))
    StringList seqs = exp[0].getMetaValue("SEQ").toStringList();
    TEST_EQUAL(seqs.size(), 1)
    TEST_EQUAL(seqs[0], "PEPTIDER")

    // Round-trip: write and re-load, SEQ must survive unchanged.
    String tmp_out("MascotGenericFile_SEQ_single_out.mgf");
    NEW_TMP_FILE(tmp_out)
    mgf_file.store(tmp_out, exp);

    PeakMap exp2;
    mgf_file.load(tmp_out, exp2);
    TEST_EQUAL(exp2.size(), 1)
    TEST_TRUE(exp2[0].metaValueExists("SEQ"))
    StringList seqs_rt = exp2[0].getMetaValue("SEQ").toStringList();
    TEST_EQUAL(seqs_rt.size(), 1)
    TEST_EQUAL(seqs_rt[0], "PEPTIDER")
  }

  // Multiple SEQ lines in one query: accumulated into a StringList.
  {
    String mgf_content = "BEGIN IONS\n"
                         "TITLE=seq_multi\n"
                         "PEPMASS=600.0\n"
                         "CHARGE=2+\n"
                         "SEQ=PEPTIDEA\n"
                         "SEQ=PEPTIDEB\n"
                         "SEQ=PEPTIDEC\n"
                         "100.0 1000.0\n"
                         "END IONS\n";

    String tmp_in("MascotGenericFile_SEQ_multi_in.mgf");
    NEW_TMP_FILE(tmp_in)
    std::ofstream ofs(tmp_in.c_str());
    ofs << mgf_content;
    ofs.close();

    PeakMap exp;
    MascotGenericFile mgf_file;
    mgf_file.load(tmp_in, exp);

    TEST_EQUAL(exp.size(), 1)
    TEST_TRUE(exp[0].metaValueExists("SEQ"))
    StringList seqs = exp[0].getMetaValue("SEQ").toStringList();
    TEST_EQUAL(seqs.size(), 3)
    TEST_EQUAL(seqs[0], "PEPTIDEA")
    TEST_EQUAL(seqs[1], "PEPTIDEB")
    TEST_EQUAL(seqs[2], "PEPTIDEC")

    // Round-trip preserves all three SEQ lines.
    String tmp_out("MascotGenericFile_SEQ_multi_out.mgf");
    NEW_TMP_FILE(tmp_out)
    mgf_file.store(tmp_out, exp);

    PeakMap exp2;
    mgf_file.load(tmp_out, exp2);
    TEST_EQUAL(exp2.size(), 1)
    StringList seqs2 = exp2[0].getMetaValue("SEQ").toStringList();
    TEST_EQUAL(seqs2.size(), 3)
    TEST_EQUAL(seqs2[0], "PEPTIDEA")
    TEST_EQUAL(seqs2[1], "PEPTIDEB")
    TEST_EQUAL(seqs2[2], "PEPTIDEC")
  }

  // Writing: a user sets SEQ programmatically before export.
  // Must work in both default and compact store modes.
  {
    MSSpectrum spec;
    spec.setNativeID("index=0");
    spec.setMSLevel(2);
    spec.setRT(100.0);
    Precursor prec;
    prec.setMZ(500.0);
    prec.setCharge(2);
    spec.getPrecursors().push_back(prec);
    Peak1D peak;
    peak.setMZ(100.0);
    peak.setIntensity(1000.0);
    spec.push_back(peak);
    spec.setMetaValue("SEQ", StringList{"PEPTIDER"});

    PeakMap exp;
    exp.addSpectrum(spec);

    MascotGenericFile mgf_file;
    stringstream ss;
    mgf_file.store(ss, "test", exp);
    TEST_TRUE(String(ss.str()).hasSubstring("SEQ=PEPTIDER"))

    stringstream compact_ss;
    mgf_file.store(compact_ss, "test", exp, true);
    TEST_TRUE(String(compact_ss.str()).hasSubstring("SEQ=PEPTIDER"))
  }

  // SEQ must not bleed across spectra during sequential load.
  {
    String mgf_content = "BEGIN IONS\n"
                         "TITLE=first\n"
                         "PEPMASS=500.0\n"
                         "SEQ=FIRSTONE\n"
                         "100.0 1000.0\n"
                         "END IONS\n"
                         "BEGIN IONS\n"
                         "TITLE=second\n"
                         "PEPMASS=600.0\n"
                         "200.0 2000.0\n"
                         "END IONS\n";

    String tmp_in("MascotGenericFile_SEQ_bleed.mgf");
    NEW_TMP_FILE(tmp_in)
    std::ofstream ofs(tmp_in.c_str());
    ofs << mgf_content;
    ofs.close();

    PeakMap exp;
    MascotGenericFile mgf_file;
    mgf_file.load(tmp_in, exp);

    TEST_EQUAL(exp.size(), 2)
    TEST_TRUE(exp[0].metaValueExists("SEQ"))
    StringList seqs_first = exp[0].getMetaValue("SEQ").toStringList();
    TEST_EQUAL(seqs_first.size(), 1)
    TEST_EQUAL(seqs_first[0], "FIRSTONE")
    // second spectrum had no SEQ - must not inherit it from the first
    TEST_FALSE(exp[1].metaValueExists("SEQ"))
  }
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
