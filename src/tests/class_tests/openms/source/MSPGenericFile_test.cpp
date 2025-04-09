// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/MSPGenericFile.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MSPGenericFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MSPGenericFile* ptr = nullptr;
MSPGenericFile* null_ptr = nullptr;
const String input_filepath = OPENMS_GET_TEST_DATA_PATH("MSPGenericFile_input.msp");

START_SECTION(MSPGenericFile())
{
  ptr = new MSPGenericFile();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MSPGenericFile())
{
  delete ptr;
}
END_SECTION

START_SECTION(void load(const String& filename, MSExperiment& experiment) const)
{
  MSPGenericFile msp;
  MSExperiment experiment;
  msp.load(input_filepath, experiment);
  const vector<MSSpectrum>& spectra = experiment.getSpectra();
  TEST_EQUAL(spectra.size(), 3)

  const MSSpectrum& s1 = spectra[0];
  TEST_EQUAL(s1.size(), 14)
  TEST_EQUAL(s1.getName(), "name1 of first")

  TEST_EQUAL(s1.metaValueExists("Synon"), true)
  TEST_STRING_EQUAL(s1.getMetaValue("Synon"), "name2 of 1st|name3 of firsttt")

  TEST_EQUAL(s1.metaValueExists("Formula"), true)
  TEST_STRING_EQUAL(s1.getMetaValue("Formula"), "A11B22C333")

  TEST_EQUAL(s1.metaValueExists("MW"), true)
  TEST_STRING_EQUAL(s1.getMetaValue("MW"), "156")

  TEST_EQUAL(s1.metaValueExists("CAS#"), true)
  TEST_STRING_EQUAL(s1.getMetaValue("CAS#"), "0123-45-6")

  TEST_EQUAL(s1.metaValueExists("NIST#"), true)
  TEST_STRING_EQUAL(s1.getMetaValue("NIST#"), "654321")

  TEST_EQUAL(s1.metaValueExists("DB#"), true)
  TEST_STRING_EQUAL(s1.getMetaValue("DB#"), "1")

  TEST_EQUAL(s1.metaValueExists("Comments"), true)
  TEST_STRING_EQUAL(s1.getMetaValue("Comments"), "Some comment")

  TEST_EQUAL(s1.metaValueExists("Num Peaks"), true)
  TEST_STRING_EQUAL(s1.getMetaValue("Num Peaks"), "14")

  TEST_EQUAL(s1[0].getPos(), 27)
  TEST_EQUAL(s1[0].getIntensity(), 29)
  TEST_EQUAL(s1[5].getPos(), 60)
  TEST_EQUAL(s1[5].getIntensity(), 41)
  TEST_EQUAL(s1[10].getPos(), 90)
  TEST_EQUAL(s1[10].getIntensity(), 168)
  TEST_EQUAL(s1[13].getPos(), 105)
  TEST_EQUAL(s1[13].getIntensity(), 36)

  const MSSpectrum& s2 = spectra[1];
  TEST_EQUAL(s2.size(), 15)
  TEST_EQUAL(s2.getName(), "name1 of second")

  TEST_EQUAL(s2.metaValueExists("Synon"), true)
  TEST_STRING_EQUAL(s2.getMetaValue("Synon"), "name2 of 2nd|name3 of seconddd")

  TEST_EQUAL(s2.metaValueExists("Formula"), true)
  TEST_STRING_EQUAL(s2.getMetaValue("Formula"), "A44B55C666")

  TEST_EQUAL(s2.metaValueExists("MW"), true)
  TEST_STRING_EQUAL(s2.getMetaValue("MW"), "589")

  TEST_EQUAL(s2.metaValueExists("CAS#"), true)
  TEST_STRING_EQUAL(s2.getMetaValue("CAS#"), "3210-45-6")

  TEST_EQUAL(s2.metaValueExists("NIST#"), true)
  TEST_STRING_EQUAL(s2.getMetaValue("NIST#"), "789564")

  TEST_EQUAL(s2.metaValueExists("DB#"), true)
  TEST_STRING_EQUAL(s2.getMetaValue("DB#"), "2")

  TEST_EQUAL(s2.metaValueExists("Comments"), true)
  TEST_STRING_EQUAL(s2.getMetaValue("Comments"), "Some other comment")

  TEST_EQUAL(s2.metaValueExists("Num Peaks"), true)
  TEST_STRING_EQUAL(s2.getMetaValue("Num Peaks"), "15")

  TEST_EQUAL(s2[0].getPos(), 27)
  TEST_EQUAL(s2[0].getIntensity(), 29)
  TEST_EQUAL(s2[5].getPos(), 260)
  TEST_EQUAL(s2[5].getIntensity(), 41)
  TEST_EQUAL(s2[10].getPos(), 290)
  TEST_EQUAL(s2[10].getIntensity(), 168)
  TEST_EQUAL(s2[14].getPos(), 310)
  TEST_EQUAL(s2[14].getIntensity(), 20)

  const MSSpectrum& s3 = spectra[2];
  TEST_EQUAL(s3.size(), 16)
  TEST_EQUAL(s3.getName(), "name1 of third")

  TEST_EQUAL(s3.metaValueExists("Synon"), true)
  TEST_STRING_EQUAL(s3.getMetaValue("Synon"), "name2 of 3rd|name3 of thirddd")

  TEST_EQUAL(s3.metaValueExists("Formula"), true)
  TEST_STRING_EQUAL(s3.getMetaValue("Formula"), "A12B12C123")

  TEST_EQUAL(s3.metaValueExists("MW"), true)
  TEST_STRING_EQUAL(s3.getMetaValue("MW"), "562")

  TEST_EQUAL(s3.metaValueExists("CAS#"), true)
  TEST_STRING_EQUAL(s3.getMetaValue("CAS#"), "4210-47-4")

  TEST_EQUAL(s3.metaValueExists("NIST#"), true)
  TEST_STRING_EQUAL(s3.getMetaValue("NIST#"), "749514")

  TEST_EQUAL(s3.metaValueExists("DB#"), true)
  TEST_STRING_EQUAL(s3.getMetaValue("DB#"), "3")

  TEST_EQUAL(s3.metaValueExists("Comments"), false) // this spectrum doesn't have a comment

  TEST_EQUAL(s3.metaValueExists("Num Peaks"), true)
  TEST_STRING_EQUAL(s3.getMetaValue("Num Peaks"), "16")

  TEST_EQUAL(s3[0].getPos(), 27)
  TEST_EQUAL(s3[0].getIntensity(), 29)
  TEST_EQUAL(s3[5].getPos(), 260)
  TEST_EQUAL(s3[5].getIntensity(), 41)
  TEST_EQUAL(s3[10].getPos(), 290)
  TEST_EQUAL(s3[10].getIntensity(), 168)
  TEST_EQUAL(s3[14].getPos(), 310)
  TEST_EQUAL(s3[14].getIntensity(), 20)
  TEST_EQUAL(s3[15].getPos(), 111)
  TEST_EQUAL(s3[15].getIntensity(), 44)
}
END_SECTION

START_SECTION(void store(const String& filename, const MSExperiment& library) const)
{
  MSPGenericFile msp;

  MSExperiment exp;
  PeakMap::SpectrumType spec;
  PeakMap::PeakType peak;

  spec.setName("first spectrum");
  spec.setMetaValue("Synon", "first1|first2|first3");
  spec.setMetaValue("CAS#", "0123-45-6");
  spec.setMetaValue("NIST#", "654321");
  spec.setRT(11.1);
  spec.setMSLevel(1);
  peak.getPosition()[0] = 1;
  peak.setIntensity(1.50f);
  spec.push_back(peak);
  peak.getPosition()[0] = 2;
  peak.setIntensity(2.5f);
  spec.push_back(peak);
  peak.getPosition()[0] = 3;
  peak.setIntensity(3.5f);
  spec.push_back(peak);
  exp.addSpectrum(spec);

  spec.clear(true);
  spec.setName("second spectrum");
  spec.setMetaValue("Synon", "second1");
  spec.setMetaValue("CAS#", "0123-45-2");
  spec.setMetaValue("NIST#", "654322");
  spec.setMetaValue("other_metadata1", "value1");
  spec.setMetaValue("other_metadata2", "value2");
  spec.setRT(22.2);
  spec.setMSLevel(1);
  peak.getPosition()[0] = 11;
  peak.setIntensity(11.50f);
  spec.push_back(peak);
  peak.getPosition()[0] = 12;
  peak.setIntensity(12.5f);
  spec.push_back(peak);
  peak.getPosition()[0] = 13;
  peak.setIntensity(13.5f);
  spec.push_back(peak);
  peak.getPosition()[0] = 14;
  peak.setIntensity(14.5f);
  spec.push_back(peak);
  peak.getPosition()[0] = 15;
  peak.setIntensity(15.5f);
  spec.push_back(peak);
  peak.getPosition()[0] = 16;
  peak.setIntensity(16.5f);
  spec.push_back(peak);
  exp.addSpectrum(spec);

  spec.clear(true);
  spec.setName("third spectrum");
  spec.setMetaValue("CAS#", "0123-45-3");
  spec.setMetaValue("NIST#", "654323");
  spec.setRT(33.3);
  spec.setMSLevel(1);
  peak.getPosition()[0] = 101;
  peak.setIntensity(101.50f);
  spec.push_back(peak);
  peak.getPosition()[0] = 102;
  peak.setIntensity(102.5f);
  spec.push_back(peak);
  peak.getPosition()[0] = 103;
  peak.setIntensity(103.5f);
  spec.push_back(peak);
  exp.addSpectrum(spec);
  
  String output_filepath;
  NEW_TMP_FILE(output_filepath)
  msp.store(output_filepath, exp);

  // read back created file
  MSExperiment exp_test;
  msp.load(output_filepath, exp_test);

  const vector<MSSpectrum>& spectra = exp_test.getSpectra();
  TEST_EQUAL(spectra.size(), 3)

  const MSSpectrum& s1 = spectra[0];
  TEST_EQUAL(s1.getName(), "first spectrum")
  TEST_EQUAL(s1.metaValueExists("Synon"), true)
  TEST_STRING_EQUAL(s1.getMetaValue("Synon"), "first1|first2|first3")
  TEST_EQUAL(s1.metaValueExists("CAS#"), true)
  TEST_STRING_EQUAL(s1.getMetaValue("CAS#"), "0123-45-6")
  TEST_EQUAL(s1.metaValueExists("NIST#"), true)
  TEST_STRING_EQUAL(s1.getMetaValue("NIST#"), "654321")
  TEST_EQUAL(s1.size(), 3)
  TEST_REAL_SIMILAR(s1[0].getMZ(), 1)
  TEST_REAL_SIMILAR(s1[0].getIntensity(), 1.5)
  TEST_REAL_SIMILAR(s1[1].getMZ(), 2)
  TEST_REAL_SIMILAR(s1[1].getIntensity(), 2.5)
  TEST_REAL_SIMILAR(s1[2].getMZ(), 3)
  TEST_REAL_SIMILAR(s1[2].getIntensity(), 3.5)

  const MSSpectrum& s2 = spectra[1];
  TEST_EQUAL(s2.getName(), "second spectrum")
  TEST_EQUAL(s2.metaValueExists("Synon"), true)
  TEST_STRING_EQUAL(s2.getMetaValue("Synon"), "second1")
  TEST_EQUAL(s2.metaValueExists("CAS#"), true)
  TEST_STRING_EQUAL(s2.getMetaValue("CAS#"), "0123-45-2")
  TEST_EQUAL(s2.metaValueExists("NIST#"), true)
  TEST_STRING_EQUAL(s2.getMetaValue("NIST#"), "654322")
  TEST_EQUAL(s2.size(), 6)
  TEST_REAL_SIMILAR(s2[0].getMZ(), 11)
  TEST_REAL_SIMILAR(s2[0].getIntensity(), 11.5)
  TEST_REAL_SIMILAR(s2[1].getMZ(), 12)
  TEST_REAL_SIMILAR(s2[1].getIntensity(), 12.5)
  TEST_REAL_SIMILAR(s2[2].getMZ(), 13)
  TEST_REAL_SIMILAR(s2[2].getIntensity(), 13.5)
  TEST_REAL_SIMILAR(s2[3].getMZ(), 14)
  TEST_REAL_SIMILAR(s2[3].getIntensity(), 14.5)
  TEST_REAL_SIMILAR(s2[4].getMZ(), 15)
  TEST_REAL_SIMILAR(s2[4].getIntensity(), 15.5)
  TEST_REAL_SIMILAR(s2[5].getMZ(), 16)
  TEST_REAL_SIMILAR(s2[5].getIntensity(), 16.5)

  const MSSpectrum& s3 = spectra[2];
  TEST_EQUAL(s3.getName(), "third spectrum")
  TEST_EQUAL(s3.metaValueExists("Synon"), false)
  TEST_EQUAL(s3.metaValueExists("CAS#"), true)
  TEST_STRING_EQUAL(s3.getMetaValue("CAS#"), "0123-45-3")
  TEST_EQUAL(s3.metaValueExists("NIST#"), true)
  TEST_STRING_EQUAL(s3.getMetaValue("NIST#"), "654323")
  TEST_EQUAL(s3.size(), 3)
  TEST_REAL_SIMILAR(s3[0].getMZ(), 101)
  TEST_REAL_SIMILAR(s3[0].getIntensity(), 101.5)
  TEST_REAL_SIMILAR(s3[1].getMZ(), 102)
  TEST_REAL_SIMILAR(s3[1].getIntensity(), 102.5)
  TEST_REAL_SIMILAR(s3[2].getMZ(), 103)
  TEST_REAL_SIMILAR(s3[2].getIntensity(), 103.5)

  // test invalid spectrum (no name).
  MSExperiment invalid_exp;
  PeakMap::SpectrumType invalid_spec;
  PeakMap::PeakType invalid_peak;
  invalid_spec.setMetaValue("Synon", "first1|first2|first3");
  invalid_spec.setMetaValue("CAS#", "0123-45-6");
  invalid_spec.setMetaValue("NIST#", "654321");
  invalid_spec.setRT(11.1);
  invalid_spec.setMSLevel(1);
  invalid_peak.getPosition()[0] = 1;
  invalid_peak.setIntensity(1.50f);
  invalid_spec.push_back(invalid_peak);
  invalid_peak.getPosition()[0] = 2;
  invalid_peak.setIntensity(2.5f);
  invalid_spec.push_back(invalid_peak);
  invalid_peak.getPosition()[0] = 3;
  invalid_peak.setIntensity(3.5f);
  invalid_spec.push_back(invalid_peak);
  invalid_exp.addSpectrum(invalid_spec);

  NEW_TMP_FILE(output_filepath)
  TEST_EXCEPTION(Exception::MissingInformation, msp.store(output_filepath, invalid_exp))
}
END_SECTION

START_SECTION(void addSpectrumToLibrary(
  MSSpectrum& spectrum,
  MSExperiment& library
))
{
  MSPGenericFile_friend msp_f;
  MSExperiment lib;

  MSSpectrum spec;
  spec.setName(""); // empty name
  spec.setMetaValue("is_valid", 1);

  TEST_EXCEPTION(Exception::MissingInformation, msp_f.addSpectrumToLibrary(spec, lib))
  TEST_EQUAL(lib.size(), 0)

  spec.setName("foo"); // Num Peaks still absent!
  TEST_EXCEPTION(Exception::MissingInformation, msp_f.addSpectrumToLibrary(spec, lib))
  TEST_EQUAL(lib.size(), 0)

  spec.setMetaValue("Num Peaks", "2");
  // Num Peaks is set but raw data poins have not been added
  TEST_EXCEPTION(Exception::ParseError, msp_f.addSpectrumToLibrary(spec, lib))
  TEST_EQUAL(lib.size(), 0)

  spec.push_back(Peak1D(1.0, 2.0));
  spec.push_back(Peak1D(3.0, 4.0)); // now the spectrum is valid
  msp_f.addSpectrumToLibrary(spec, lib);
  TEST_EQUAL(lib.size(), 1)

  spec.setName("bar");
  spec.setMetaValue("is_valid", 1);
  msp_f.addSpectrumToLibrary(spec, lib);
  TEST_EQUAL(lib.size(), 2)

  spec.setMetaValue("is_valid", 1);
  msp_f.addSpectrumToLibrary(spec, lib); // duplicate, won't be added
  TEST_EQUAL(lib.size(), 2)

  spec.setMetaValue("is_valid", 0);
  spec.setName("not a duplicate");
  msp_f.addSpectrumToLibrary(spec, lib);
  TEST_EQUAL(lib.size(), 2)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
