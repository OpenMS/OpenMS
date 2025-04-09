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

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
