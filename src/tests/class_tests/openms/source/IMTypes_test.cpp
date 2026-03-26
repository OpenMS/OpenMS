// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(MSRunIMSplitter, "$Id$")

/////////////////////////////////////////////////////////////

IMTypes* e_ptr = nullptr;
IMTypes* e_nullPointer = nullptr;

START_SECTION((IMTypes()))
	e_ptr = new IMTypes;
  TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~IMTypes()))
	delete e_ptr;
END_SECTION

START_SECTION((DriftTimeUnit toDriftTimeUnit(const String& dtu_string)))
  TEST_EQUAL(toDriftTimeUnit("<NONE>") == DriftTimeUnit::NONE, true)
  for (size_t i = 0; i < (size_t)DriftTimeUnit::SIZE_OF_DRIFTTIMEUNIT; ++i)
  {
    TEST_EQUAL((size_t)toDriftTimeUnit(NamesOfDriftTimeUnit[i]), i)
  }
  TEST_EXCEPTION(Exception::InvalidValue, toDriftTimeUnit("haha"));
END_SECTION

START_SECTION(const String& driftTimeUnitToString(const DriftTimeUnit value))
  TEST_EQUAL(driftTimeUnitToString(DriftTimeUnit::NONE), "<NONE>")
  for (size_t i = 0; i < (size_t)DriftTimeUnit::SIZE_OF_DRIFTTIMEUNIT; ++i)
  {
    TEST_EQUAL(driftTimeUnitToString(DriftTimeUnit(i)), NamesOfDriftTimeUnit[i])
  }
  TEST_EXCEPTION(Exception::InvalidValue, driftTimeUnitToString(DriftTimeUnit::SIZE_OF_DRIFTTIMEUNIT));
END_SECTION


START_SECTION((IMFormat toIMFormat(const String& IM_format)))
  TEST_EQUAL(toIMFormat("none") == IMFormat::NONE, true)
  for (size_t i = 0; i < (size_t) IMFormat::SIZE_OF_IMFORMAT; ++i)
  {
    TEST_EQUAL((size_t) toIMFormat(NamesOfIMFormat[i]), i)
  }
  TEST_EXCEPTION(Exception::InvalidValue, toIMFormat("haha"));
END_SECTION

START_SECTION(const String& imFormatToString(const IMFormat value))
  TEST_EQUAL(imFormatToString(IMFormat::NONE), "none")
  for (size_t i = 0; i < (size_t)IMFormat::SIZE_OF_IMFORMAT; ++i)
  {
    TEST_EQUAL(imFormatToString(IMFormat(i)), NamesOfIMFormat[i])
  }
  TEST_EXCEPTION(Exception::InvalidValue, imFormatToString(IMFormat::SIZE_OF_IMFORMAT));

END_SECTION

START_SECTION(IMPeakType string conversions)
{
  TEST_EQUAL(toIMPeakType("im_profile"), IMPeakType::IM_PROFILE)
  TEST_EQUAL(toIMPeakType("im_centroided"), IMPeakType::IM_CENTROIDED)
  TEST_EQUAL(toIMPeakType("unknown"), IMPeakType::UNKNOWN)
  TEST_EQUAL(imPeakTypeToString(IMPeakType::IM_PROFILE), "im_profile")
  TEST_EQUAL(imPeakTypeToString(IMPeakType::IM_CENTROIDED), "im_centroided")
  TEST_EQUAL(imPeakTypeToString(IMPeakType::UNKNOWN), "unknown")
  TEST_EXCEPTION(Exception::InvalidValue, toIMPeakType("garbage"))
}
END_SECTION


// single IM value for whole spec
const MSSpectrum IMwithDrift = [&]() { 
  MSSpectrum spec; 
  spec.setDriftTime(123.4);
  spec.setDriftTimeUnit(DriftTimeUnit::VSSC);
  return spec;
}();

// convert to IM-Frame with float meta-data array
const MSSpectrum IMwithFDA = [&]() {
  MSExperiment exp;
  exp.addSpectrum(IMwithDrift);
  auto single = IMDataConverter::reshapeIMFrameToSingle(exp);
  return single[0];
}();

START_SECTION(static IMFormat determineIMFormat(const MSExperiment& exp, int ms_level))

  // empty experiment
  TEST_EQUAL(IMTypes::determineIMFormat(MSExperiment(), 1) == IMFormat::NONE, true)

  {
    MSExperiment exp;
    exp.addSpectrum(MSSpectrum()); // default MS level = 1
    exp.addSpectrum(MSSpectrum());
    TEST_EQUAL(IMTypes::determineIMFormat(exp, 1) == IMFormat::NONE, true)
  }

  {
    auto ms1_drift = IMwithDrift;
    ms1_drift.setMSLevel(1);
    MSExperiment exp;
    exp.addSpectrum(MSSpectrum());
    exp.addSpectrum(ms1_drift);
    TEST_EQUAL(IMTypes::determineIMFormat(exp, 1) == IMFormat::IM_SPECTRUM, true)
  }

  {
    auto ms1_peak = IMwithFDA;
    ms1_peak.setMSLevel(1);
    MSExperiment exp;
    exp.addSpectrum(MSSpectrum());
    exp.addSpectrum(ms1_peak);
    TEST_EQUAL(IMTypes::determineIMFormat(exp, 1) == IMFormat::IM_PEAK, true)
  }

  {
    // MS1 = IM_PEAK, MS2 = IM_SPECTRUM — per-level queries work independently
    auto ms1_peak = IMwithFDA;
    ms1_peak.setMSLevel(1);
    auto ms2_drift = IMwithDrift;
    ms2_drift.setMSLevel(2);
    MSExperiment exp;
    exp.addSpectrum(ms1_peak);
    exp.addSpectrum(ms2_drift);
    TEST_EQUAL(IMTypes::determineIMFormat(exp, 1) == IMFormat::IM_PEAK, true)
    TEST_EQUAL(IMTypes::determineIMFormat(exp, 2) == IMFormat::IM_SPECTRUM, true)
  }

  {
    // no spectra of requested level
    MSExperiment exp;
    auto ms1 = IMwithDrift;
    ms1.setMSLevel(1);
    exp.addSpectrum(ms1);
    TEST_EQUAL(IMTypes::determineIMFormat(exp, 2) == IMFormat::NONE, true)
  }

  {
    // mixed formats within same MS level throws
    auto ms1_drift = IMwithDrift;
    ms1_drift.setMSLevel(1);
    auto ms1_peak = IMwithFDA;
    ms1_peak.setMSLevel(1);
    MSExperiment exp;
    exp.addSpectrum(ms1_drift);
    exp.addSpectrum(ms1_peak);
    TEST_EXCEPTION(Exception::InvalidValue, IMTypes::determineIMFormat(exp, 1))
  }

END_SECTION

START_SECTION(static IMFormat determineIMFormat(const MSSpectrum& spec))
   TEST_EQUAL(IMTypes::determineIMFormat(MSSpectrum()) == IMFormat::NONE, true)
   
   // single IM value for whole spec
   TEST_EQUAL(IMTypes::determineIMFormat(IMwithDrift) == IMFormat::IM_SPECTRUM, true)

   // convert to IM-Frame with float meta-data array
   TEST_EQUAL(IMTypes::determineIMFormat(IMwithFDA) == IMFormat::IM_PEAK, true)

   // set both ... is valid (typically concatenated + some average value)
   auto IMwithFDA2 = IMwithFDA;
   IMwithFDA2.setDriftTime(123.4);
   TEST_EQUAL(IMTypes::determineIMFormat(IMwithFDA2) == IMFormat::IM_PEAK, true)
END_SECTION

START_SECTION(determineIMFormat returns IM_PEAK for centroided IM data)
{
  MSSpectrum s;
  MSSpectrum::FloatDataArray fda;
  fda.setName("Ion Mobility");
  s.getFloatDataArrays().push_back(fda);
  s.setIMPeakType(IMPeakType::IM_CENTROIDED);
  TEST_EQUAL(IMTypes::determineIMFormat(s), IMFormat::IM_PEAK)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
