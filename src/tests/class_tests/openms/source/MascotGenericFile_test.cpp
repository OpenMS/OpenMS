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

START_SECTION((template < typename MapType > void load(const std::string &filename, MapType &exp)))
{
  PeakMap exp;
  ptr->load(OPENMS_GET_TEST_DATA_PATH("MascotInfile_test.mascot_in"), exp);
  TEST_EQUAL(exp.size(), 1)

  TEST_EQUAL(exp.begin()->size(), 9)
}
END_SECTION

START_SECTION((void store(std::ostream &os, const std::string &filename, const PeakMap &experiment, bool compact = false)))
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

  vector<std::string> strings;
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

  std::string mgf_file(ss.str());
  for (Size i = 0; i < strings.size(); ++i)
  {
    TEST_EQUAL(StringUtils::hasSubstring(mgf_file, strings[i]), true)
  }

  // test of making default TITLE
  exp[0].removeMetaValue("TITLE");
  stringstream ss2;
  ptr->store(ss2, "test", exp);
  vector<std::string> strings2;
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
  std::string mgf_file2(ss2.str());
  for (Size i = 0; i < strings2.size(); ++i)
  {
    TEST_EQUAL(StringUtils::hasSubstring(mgf_file2, strings2[i]), true)
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
  std::string content = ("BEGIN IONS\n"
                    "TITLE=901.23457_234.568_index=250_test\n"
                    "PEPMASS=901.23457\n"
                    "RTINSECONDS=234.568\n"
                    "SCANS=250\n"
                    "890.12346 2345.679\n"
                    "END IONS");
  TEST_EQUAL(StringUtils::hasSubstring(mgf_file, content), true);
}
END_SECTION

START_SECTION((void store(const std::string &filename, const PeakMap &experiment, bool compact = false)))
{
  std::string tmp_name("MascotGenericFile_1.tmp");
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
  std::string mgf_content = "BEGIN IONS\n"
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
  std::string tmp_name("test_compound_name.mgf");
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
  TEST_EQUAL(StringUtils::toStr(exp[0].getMetaValue(Constants::UserParam::MSM_METABOLITE_NAME)), "Caffeine")
  
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
  TEST_EQUAL(StringUtils::toStr(exp[0].getMetaValue("GNPS_Spectrum_ID")), "CCMSLIB00000001547")
  
  // Test that COMPOUND_NAME was correctly mapped to MSM_METABOLITE_NAME
  TEST_EQUAL(exp[0].metaValueExists(Constants::UserParam::MSM_METABOLITE_NAME), true)
  TEST_EQUAL(StringUtils::toStr(exp[0].getMetaValue(Constants::UserParam::MSM_METABOLITE_NAME)), "3-Des-Microcystein_LR")
  
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
    std::string mgf_content = "BEGIN IONS\n"
                         "TITLE=seq_single\n"
                         "PEPMASS=500.0\n"
                         "CHARGE=2+\n"
                         "SEQ=PEPTIDER\n"
                         "100.0 1000.0\n"
                         "200.0 2000.0\n"
                         "END IONS\n";

    std::string tmp_in("MascotGenericFile_SEQ_single_in.mgf");
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
    std::string tmp_out("MascotGenericFile_SEQ_single_out.mgf");
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
    std::string mgf_content = "BEGIN IONS\n"
                         "TITLE=seq_multi\n"
                         "PEPMASS=600.0\n"
                         "CHARGE=2+\n"
                         "SEQ=PEPTIDEA\n"
                         "SEQ=PEPTIDEB\n"
                         "SEQ=PEPTIDEC\n"
                         "100.0 1000.0\n"
                         "END IONS\n";

    std::string tmp_in("MascotGenericFile_SEQ_multi_in.mgf");
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
    std::string tmp_out("MascotGenericFile_SEQ_multi_out.mgf");
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
    TEST_TRUE(StringUtils::hasSubstring(ss.str(), "SEQ=PEPTIDER"))

    stringstream compact_ss;
    mgf_file.store(compact_ss, "test", exp, true);
    TEST_TRUE(StringUtils::hasSubstring(compact_ss.str(), "SEQ=PEPTIDER"))
  }

  // SEQ must not bleed across spectra during sequential load.
  {
    std::string mgf_content = "BEGIN IONS\n"
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

    std::string tmp_in("MascotGenericFile_SEQ_bleed.mgf");
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

START_SECTION((precursor charge, precursor intensity and RT are not inherited from the previous block))
{
  // Regression test: the reader reused one spectrum object for all blocks and only reset the
  // peaks, TITLE and SEQ. A block without CHARGE, RTINSECONDS or a PEPMASS intensity therefore
  // silently received the values of the previous block (wrong precursor charge/mass, RT, ...).
  std::string mgf_content = "BEGIN IONS\n"
                            "TITLE=s1\n"
                            "PEPMASS=500.25 12000\n"
                            "CHARGE=2+\n"
                            "RTINSECONDS=1200\n"
                            "100 10\n"
                            "END IONS\n"
                            "BEGIN IONS\n"
                            "TITLE=s2\n"
                            "PEPMASS=612.8\n"
                            "150 20\n"
                            "END IONS\n";

  std::string tmp_in("MascotGenericFile_no_carry_over.mgf");
  NEW_TMP_FILE(tmp_in)
  std::ofstream ofs(tmp_in.c_str());
  ofs << mgf_content;
  ofs.close();

  PeakMap exp;
  MascotGenericFile mgf_file;
  mgf_file.load(tmp_in, exp);

  TEST_EQUAL(exp.size(), 2)
  ABORT_IF(exp.size() != 2)

  // s1: all fields are given in its block
  TEST_EQUAL(exp[0].getPrecursors().size(), 1)
  TEST_REAL_SIMILAR(exp[0].getPrecursors()[0].getMZ(), 500.25)
  TEST_REAL_SIMILAR(exp[0].getPrecursors()[0].getIntensity(), 12000.0)
  TEST_EQUAL(exp[0].getPrecursors()[0].getCharge(), 2)
  TEST_REAL_SIMILAR(exp[0].getRT(), 1200.0)
  TEST_EQUAL(exp[0].getMSLevel(), 2)
  TEST_STRING_EQUAL(exp[0].getNativeID(), "index=0")
  TEST_STRING_EQUAL(StringUtils::toStr(exp[0].getMetaValue("TITLE")), "s1_index=0")
  TEST_EQUAL(exp[0].size(), 1)
  TEST_REAL_SIMILAR(exp[0][0].getMZ(), 100.0)

  // s2: only the PEPMASS m/z is given -> everything else must be at its default value
  // (charge 0, intensity 0, RT unset i.e. -1) and not at the values of s1
  TEST_EQUAL(exp[1].getPrecursors().size(), 1)
  TEST_REAL_SIMILAR(exp[1].getPrecursors()[0].getMZ(), 612.8)
  TEST_REAL_SIMILAR(exp[1].getPrecursors()[0].getIntensity(), 0.0)
  TEST_EQUAL(exp[1].getPrecursors()[0].getCharge(), 0)
  TEST_REAL_SIMILAR(exp[1].getRT(), -1.0) // MSSpectrum default; getRT() < 0 means 'no RT'
  TEST_EQUAL(exp[1].getMSLevel(), 2)
  TEST_TRUE(exp[1].getType() == SpectrumSettings::SpectrumType::CENTROID)
  TEST_STRING_EQUAL(exp[1].getNativeID(), "index=1")
  TEST_STRING_EQUAL(StringUtils::toStr(exp[1].getMetaValue("TITLE")), "s2_index=1")
  TEST_EQUAL(exp[1].size(), 1)
  TEST_REAL_SIMILAR(exp[1][0].getMZ(), 150.0)
  TEST_REAL_SIMILAR(exp[1][0].getIntensity(), 20.0)
}
END_SECTION

START_SECTION((MS level and library fields (NAME, SMILES, INCHI, ...) are not inherited from the previous block))
{
  // Spectral-library entries differ in which fields they list; a block without NAME/SMILES/...
  // must not report the compound of the previous block.
  std::string mgf_content = "BEGIN IONS\n"
                            "TITLE=lib1\n"
                            "PEPMASS=195.0877\n"
                            "CHARGE=1\n"
                            "MSLEVEL=3\n"
                            "NAME=Caffeine\n"
                            "SMILES=CC\n"
                            "INCHI=InChI=1S/C8H10N4O2\n"
                            "IONMODE=positive\n"
                            "SPECTRUMID=CCMSLIB00000000001\n"
                            "SCANS=1\n"
                            "100 10\n"
                            "END IONS\n"
                            "BEGIN IONS\n"
                            "TITLE=lib2\n"
                            "PEPMASS=181.0720\n"
                            "150 20\n"
                            "END IONS\n";

  std::string tmp_in("MascotGenericFile_no_carry_over_library.mgf");
  NEW_TMP_FILE(tmp_in)
  std::ofstream ofs(tmp_in.c_str());
  ofs << mgf_content;
  ofs.close();

  PeakMap exp;
  MascotGenericFile mgf_file;
  mgf_file.load(tmp_in, exp);

  TEST_EQUAL(exp.size(), 2)
  ABORT_IF(exp.size() != 2)

  // block 1: everything present
  TEST_EQUAL(exp[0].getMSLevel(), 3)
  TEST_EQUAL(exp[0].getPrecursors()[0].getCharge(), 1)
  TEST_TRUE(exp[0].metaValueExists(Constants::UserParam::MSM_METABOLITE_NAME))
  TEST_STRING_EQUAL(StringUtils::toStr(exp[0].getMetaValue(Constants::UserParam::MSM_METABOLITE_NAME)), "Caffeine")
  TEST_TRUE(exp[0].metaValueExists(Constants::UserParam::MSM_SMILES_STRING))
  TEST_STRING_EQUAL(StringUtils::toStr(exp[0].getMetaValue(Constants::UserParam::MSM_SMILES_STRING)), "CC")
  TEST_TRUE(exp[0].metaValueExists(Constants::UserParam::MSM_INCHI_STRING))
  TEST_STRING_EQUAL(StringUtils::toStr(exp[0].getMetaValue(Constants::UserParam::MSM_INCHI_STRING)), "InChI=1S/C8H10N4O2")
  TEST_TRUE(exp[0].metaValueExists("IONMODE"))
  TEST_TRUE(exp[0].metaValueExists("GNPS_Spectrum_ID"))
  TEST_TRUE(exp[0].metaValueExists("Scan_ID"))

  // block 2: none of these fields is given -> none of them may be present
  TEST_REAL_SIMILAR(exp[1].getPrecursors()[0].getMZ(), 181.0720)
  TEST_EQUAL(exp[1].getPrecursors()[0].getCharge(), 0)
  TEST_EQUAL(exp[1].getMSLevel(), 2) // MGF default, not the MSLEVEL=3 of block 1
  TEST_FALSE(exp[1].metaValueExists(Constants::UserParam::MSM_METABOLITE_NAME))
  TEST_FALSE(exp[1].metaValueExists(Constants::UserParam::MSM_SMILES_STRING))
  TEST_FALSE(exp[1].metaValueExists(Constants::UserParam::MSM_INCHI_STRING))
  TEST_FALSE(exp[1].metaValueExists("IONMODE"))
  TEST_FALSE(exp[1].metaValueExists("GNPS_Spectrum_ID"))
  TEST_FALSE(exp[1].metaValueExists("Scan_ID"))
  TEST_STRING_EQUAL(StringUtils::toStr(exp[1].getMetaValue("TITLE")), "lib2_index=1")
}
END_SECTION

START_SECTION((a block without peak lines is read as an empty spectrum and not merged with the next block))
{
  // Regression test: 'END IONS' was only detected inside the peak loop, so a block without
  // peak lines ran on into the following block (merging both, and shifting all native IDs)
  // and was dropped when it was the last block of the file.
  {
    std::string mgf_content = "BEGIN IONS\n"
                              "TITLE=first\n"
                              "PEPMASS=500.0\n"
                              "CHARGE=2+\n"
                              "100 10\n"
                              "END IONS\n"
                              "BEGIN IONS\n"
                              "TITLE=empty\n"
                              "PEPMASS=600.0\n"
                              "CHARGE=3+\n"
                              "END IONS\n"
                              "BEGIN IONS\n"
                              "TITLE=third\n"
                              "PEPMASS=700.0\n"
                              "200 20\n"
                              "END IONS\n";

    std::string tmp_in("MascotGenericFile_empty_block.mgf");
    NEW_TMP_FILE(tmp_in)
    std::ofstream ofs(tmp_in.c_str());
    ofs << mgf_content;
    ofs.close();

    PeakMap exp;
    MascotGenericFile mgf_file;
    mgf_file.load(tmp_in, exp);

    TEST_EQUAL(exp.size(), 3)
    ABORT_IF(exp.size() != 3)

    TEST_STRING_EQUAL(StringUtils::toStr(exp[0].getMetaValue("TITLE")), "first_index=0")
    TEST_REAL_SIMILAR(exp[0].getPrecursors()[0].getMZ(), 500.0)
    TEST_EQUAL(exp[0].getPrecursors()[0].getCharge(), 2)
    TEST_EQUAL(exp[0].size(), 1)

    // the empty block keeps its own header fields, but has no peaks
    TEST_STRING_EQUAL(exp[1].getNativeID(), "index=1")
    TEST_STRING_EQUAL(StringUtils::toStr(exp[1].getMetaValue("TITLE")), "empty_index=1")
    TEST_REAL_SIMILAR(exp[1].getPrecursors()[0].getMZ(), 600.0)
    TEST_EQUAL(exp[1].getPrecursors()[0].getCharge(), 3)
    TEST_EQUAL(exp[1].size(), 0)

    // the third block is a spectrum of its own: its TITLE and PEPMASS are on the third
    // spectrum, and it does not inherit the CHARGE of the empty block
    TEST_STRING_EQUAL(exp[2].getNativeID(), "index=2")
    TEST_STRING_EQUAL(StringUtils::toStr(exp[2].getMetaValue("TITLE")), "third_index=2")
    TEST_REAL_SIMILAR(exp[2].getPrecursors()[0].getMZ(), 700.0)
    TEST_EQUAL(exp[2].getPrecursors()[0].getCharge(), 0)
    TEST_EQUAL(exp[2].size(), 1)
    TEST_REAL_SIMILAR(exp[2][0].getMZ(), 200.0)
  }

  // an empty block at the end of the file is kept as well
  {
    std::string mgf_content = "BEGIN IONS\n"
                              "TITLE=first\n"
                              "PEPMASS=500.0\n"
                              "100 10\n"
                              "END IONS\n"
                              "BEGIN IONS\n"
                              "TITLE=last_empty\n"
                              "PEPMASS=600.0\n"
                              "END IONS\n";

    std::string tmp_in("MascotGenericFile_empty_last_block.mgf");
    NEW_TMP_FILE(tmp_in)
    std::ofstream ofs(tmp_in.c_str());
    ofs << mgf_content;
    ofs.close();

    PeakMap exp;
    MascotGenericFile mgf_file;
    mgf_file.load(tmp_in, exp);

    TEST_EQUAL(exp.size(), 2)
    ABORT_IF(exp.size() != 2)
    TEST_STRING_EQUAL(StringUtils::toStr(exp[1].getMetaValue("TITLE")), "last_empty_index=1")
    TEST_REAL_SIMILAR(exp[1].getPrecursors()[0].getMZ(), 600.0)
    TEST_EQUAL(exp[1].size(), 0)
  }
}
END_SECTION

START_SECTION((parameters in the file header are not applied to spectra))
{
  // positive control for the documented behaviour: lines outside BEGIN IONS/END IONS
  // (Mascot search parameters) are skipped and do not act as defaults for the blocks
  std::string mgf_content = "CHARGE=1,2,3\n"
                            "MASS=monoisotopic\n"
                            "\n"
                            "BEGIN IONS\n"
                            "TITLE=no_charge\n"
                            "PEPMASS=500.0\n"
                            "100 10\n"
                            "END IONS\n";

  std::string tmp_in("MascotGenericFile_header_params.mgf");
  NEW_TMP_FILE(tmp_in)
  std::ofstream ofs(tmp_in.c_str());
  ofs << mgf_content;
  ofs.close();

  PeakMap exp;
  MascotGenericFile mgf_file;
  mgf_file.load(tmp_in, exp);

  TEST_EQUAL(exp.size(), 1)
  TEST_REAL_SIMILAR(exp[0].getPrecursors()[0].getMZ(), 500.0)
  TEST_EQUAL(exp[0].getPrecursors()[0].getCharge(), 0)
  TEST_EQUAL(exp[0].size(), 1)
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
