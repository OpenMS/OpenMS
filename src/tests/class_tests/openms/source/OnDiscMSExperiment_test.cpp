// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
///////////////////////////

START_TEST(OnDiscMSExperiment, "$Id$");

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

OnDiscPeakMap* ptr = nullptr;
OnDiscPeakMap* nullPointer = nullptr;
START_SECTION((OnDiscMSExperiment()))
{
  ptr = new OnDiscPeakMap();
  TEST_NOT_EQUAL(ptr, nullPointer);
}
END_SECTION

START_SECTION((~OnDiscMSExperiment()))
{
  delete ptr;
}
END_SECTION

START_SECTION((OnDiscMSExperiment(const OnDiscMSExperiment& filename)))
{
  OnDiscPeakMap tmp; 
  tmp.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  OnDiscPeakMap tmp2(tmp);
  TEST_EQUAL(tmp2.getExperimentalSettings()->getInstrument().getName(), tmp.getExperimentalSettings()->getInstrument().getName() )
  TEST_EQUAL(tmp2.getExperimentalSettings()->getInstrument().getVendor(), tmp.getExperimentalSettings()->getInstrument().getVendor() )
  TEST_EQUAL(tmp2.getExperimentalSettings()->getInstrument().getModel(), tmp.getExperimentalSettings()->getInstrument().getModel() )
  TEST_EQUAL(tmp2.getExperimentalSettings()->getInstrument().getMassAnalyzers().size(), tmp.getExperimentalSettings()->getInstrument().getMassAnalyzers().size() )
  TEST_EQUAL(tmp2.size(),tmp.size());
}
END_SECTION

// START_SECTION((OnDiscMSExperiment(const String& filename)))
// {
//   OnDiscPeakMap tmp(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
//   TEST_EQUAL(tmp.size(), 2);
// }
// END_SECTION

START_SECTION((bool operator== (const OnDiscMSExperiment& rhs) const))
{
  OnDiscPeakMap tmp; tmp.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  OnDiscPeakMap tmp2; tmp2.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"), true);
  OnDiscPeakMap same; same.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  OnDiscPeakMap failed; failed.openFile(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"));

  TEST_TRUE(tmp == same);
  TEST_EQUAL(tmp2==same, false);
  TEST_TRUE(tmp2 == tmp2);
  TEST_EQUAL((*tmp.getExperimentalSettings())==(*same.getExperimentalSettings()), true);
  TEST_EQUAL(tmp==failed, false);
}
END_SECTION

START_SECTION((bool operator!= (const OnDiscMSExperiment& rhs) const))
{
  OnDiscPeakMap tmp; tmp.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  OnDiscPeakMap tmp2; tmp2.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"), true);
  OnDiscPeakMap same; same.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  OnDiscPeakMap failed; failed.openFile(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"));

  TEST_EQUAL(tmp!=same, false);
  TEST_FALSE(tmp2 == same);
  TEST_FALSE(tmp == failed);
}
END_SECTION

START_SECTION(( bool openFile(const String& filename, bool skipMetaData = false) ))
{
  OnDiscPeakMap tmp;
  OnDiscPeakMap same;
  OnDiscPeakMap failed;

  bool res;
  res = tmp.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_EQUAL(res, true)

  res = tmp.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"), true);
  TEST_EQUAL(res, true)

  res = same.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_EQUAL(res, true)

  res = same.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"), true);
  TEST_EQUAL(res, true)

  res = failed.openFile(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"));
  TEST_EQUAL(res, false)

  res = failed.openFile(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), true);
  TEST_EQUAL(res, false)
}
END_SECTION

START_SECTION((bool isSortedByRT() const))
{
  OnDiscPeakMap tmp; tmp.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_EQUAL(tmp.isSortedByRT(), true);
}
END_SECTION

START_SECTION((Size size() const))
{
  OnDiscPeakMap tmp; tmp.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  OnDiscPeakMap tmp2; tmp2.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"), true);
  OnDiscPeakMap failed; failed.openFile(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"));
  TEST_EQUAL(tmp.size(), 2);
  TEST_EQUAL(tmp2.size(), 2);
  TEST_EQUAL(failed.size(), 0);
}
END_SECTION

START_SECTION((bool empty() const))
{
  OnDiscPeakMap tmp; tmp.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  OnDiscPeakMap tmp2; tmp2.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"), true);
  OnDiscPeakMap failed; failed.openFile(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"));
  TEST_EQUAL(tmp.empty(), false);
  TEST_EQUAL(tmp2.empty(), false);
  TEST_EQUAL(failed.empty(), true);
}
END_SECTION

START_SECTION((Size getNrSpectra() const))
{
  OnDiscPeakMap tmp; tmp.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  OnDiscPeakMap tmp2; tmp2.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"), true);
  OnDiscPeakMap failed; failed.openFile(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"));
  TEST_EQUAL(tmp.getNrSpectra(), 2);
  TEST_EQUAL(tmp2.getNrSpectra(), 2);
  TEST_EQUAL(failed.getNrSpectra(), 0);
}
END_SECTION

START_SECTION((Size getNrChromatograms() const))
{
  OnDiscPeakMap tmp; tmp.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  OnDiscPeakMap tmp2; tmp2.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"), true);
  OnDiscPeakMap failed; failed.openFile(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"));
  TEST_EQUAL(tmp.getNrChromatograms(), 1);
  TEST_EQUAL(tmp2.getNrChromatograms(), 1);
  TEST_EQUAL(failed.getNrChromatograms(), 0);
}
END_SECTION

START_SECTION((std::shared_ptr<const ExperimentalSettings> getExperimentalSettings() const))
{
  OnDiscPeakMap tmp; tmp.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  OnDiscPeakMap tmp2; tmp2.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"), true);
  std::shared_ptr<const ExperimentalSettings> settings = tmp.getExperimentalSettings();

  TEST_EQUAL(settings->getInstrument().getName(), "LTQ FT")
  TEST_EQUAL(settings->getInstrument().getMassAnalyzers().size(), 1)

  settings = tmp2.getExperimentalSettings();
  TEST_TRUE(settings == nullptr)
}
END_SECTION

START_SECTION((MSSpectrum operator[] (Size n) const))
{
  OnDiscPeakMap tmp; tmp.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_EQUAL(tmp.empty(), false);
  MSSpectrum s = tmp[0];
  TEST_EQUAL(s.empty(), false);
  TEST_EQUAL(s.size(), 19914);

  OnDiscPeakMap tmp2; tmp2.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"), true);
  TEST_EQUAL(tmp2.empty(), false);
  s = tmp2[0];
  TEST_EQUAL(s.empty(), false);
  TEST_EQUAL(s.size(), 19914);
}
END_SECTION

START_SECTION((MSSpectrum getSpectrum(Size id)))
{
  OnDiscPeakMap tmp; tmp.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_EQUAL(tmp.empty(), false);
  MSSpectrum s = tmp.getSpectrum(0);
  TEST_EQUAL(s.empty(), false);
  TEST_EQUAL(s.size(), 19914);

  OnDiscPeakMap tmp2; tmp2.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"), true);
  TEST_EQUAL(tmp2.empty(), false);
  MSSpectrum s2 = tmp2.getSpectrum(0);
  TEST_EQUAL(s2.empty(), false);
  TEST_EQUAL(s2.size(), 19914);
  MSSpectrum s3 = tmp2.getSpectrum(1);
  TEST_EQUAL(s3.empty(), false);
  TEST_EQUAL(s3.size(), 19800);
}
END_SECTION

START_SECTION(OpenMS::Interfaces::SpectrumPtr getSpectrumById(Size id))
{
  OnDiscPeakMap tmp; tmp.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_EQUAL(tmp.empty(), false);
  OpenMS::Interfaces::SpectrumPtr s = tmp.getSpectrumById(0);
  TEST_EQUAL(s->getMZArray()->data.empty(), false);
  TEST_EQUAL(s->getMZArray()->data.size(), 19914);
  TEST_EQUAL(s->getIntensityArray()->data.empty(), false);
  TEST_EQUAL(s->getIntensityArray()->data.size(), 19914);

  OnDiscPeakMap tmp2; tmp2.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"), true);
  TEST_EQUAL(tmp2.empty(), false);
  s = tmp2.getSpectrumById(0);
  TEST_EQUAL(s->getMZArray()->data.empty(), false);
  TEST_EQUAL(s->getMZArray()->data.size(), 19914);
  TEST_EQUAL(s->getIntensityArray()->data.empty(), false);
  TEST_EQUAL(s->getIntensityArray()->data.size(), 19914);
}
END_SECTION

START_SECTION((MSChromatogram getChromatogram(Size id)))
{
  OnDiscPeakMap tmp; tmp.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_EQUAL(tmp.getNrChromatograms(), 1);
  TEST_EQUAL(tmp.empty(), false);
  MSChromatogram c = tmp.getChromatogram(0);
  TEST_EQUAL(c.empty(), false);
  TEST_EQUAL(c.size(), 48);

  OnDiscPeakMap tmp2; tmp2.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"), true);
  TEST_EQUAL(tmp2.getNrChromatograms(), 1);
  TEST_EQUAL(tmp2.empty(), false);
  c = tmp2.getChromatogram(0);
  TEST_EQUAL(c.empty(), false);
  TEST_EQUAL(c.size(), 48);
}
END_SECTION

START_SECTION(OpenMS::Interfaces::ChromatogramPtr getChromatogramById(Size id))
{
  OnDiscPeakMap tmp; tmp.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_EQUAL(tmp.empty(), false);
  OpenMS::Interfaces::ChromatogramPtr s = tmp.getChromatogramById(0);
  TEST_EQUAL(s->getTimeArray()->data.empty(), false);
  TEST_EQUAL(s->getTimeArray()->data.size(), 48);
  TEST_EQUAL(s->getIntensityArray()->data.empty(), false);
  TEST_EQUAL(s->getIntensityArray()->data.size(), 48);

  OnDiscPeakMap tmp2; tmp2.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"), true);
  TEST_EQUAL(tmp2.empty(), false);
  s = tmp2.getChromatogramById(0);
  TEST_EQUAL(s->getTimeArray()->data.empty(), false);
  TEST_EQUAL(s->getTimeArray()->data.size(), 48);
  TEST_EQUAL(s->getIntensityArray()->data.empty(), false);
  TEST_EQUAL(s->getIntensityArray()->data.size(), 48);
}
END_SECTION

START_SECTION(MSChromatogram getChromatogramByNativeId(const std::string& id))
{
  OnDiscPeakMap tmp; tmp.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_EQUAL(tmp.empty(), false);
  OpenMS::MSChromatogram s = tmp.getChromatogramByNativeId("TIC");
  TEST_EQUAL(s.empty(), false);
  TEST_EQUAL(s.size(), 48);
  TEST_EXCEPTION(Exception::IllegalArgument, tmp.getChromatogramByNativeId("TIK"))

  OnDiscPeakMap tmp2; tmp2.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"), true);
  TEST_EQUAL(tmp2.empty(), false);
  s = tmp2.getChromatogramByNativeId("TIC");
  TEST_EQUAL(s.empty(), false);
  TEST_EQUAL(s.size(), 48);
  TEST_EXCEPTION(Exception::IllegalArgument, tmp2.getChromatogramByNativeId("TIK"))
}
END_SECTION

START_SECTION(MSMSSpectrum getSpectrumByNativeId(const std::string& id))
{
  OnDiscPeakMap tmp; tmp.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_EQUAL(tmp.empty(), false);
  OpenMS::MSSpectrum s = tmp.getSpectrumByNativeId("controllerType=0 controllerNumber=1 scan=1");
  TEST_EQUAL(s.empty(), false);
  TEST_EQUAL(s.size(), 19914);
  s = tmp.getSpectrumByNativeId("controllerType=0 controllerNumber=1 scan=2");
  TEST_EQUAL(s.empty(), false);
  TEST_EQUAL(s.size(), 19800);
  TEST_EXCEPTION(Exception::IllegalArgument, tmp.getSpectrumByNativeId("TIK"))

  OnDiscPeakMap tmp2; tmp2.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"), true);
  TEST_EQUAL(tmp2.empty(), false);
  s = tmp2.getSpectrumByNativeId("controllerType=0 controllerNumber=1 scan=1");
  TEST_EQUAL(s.empty(), false);
  TEST_EQUAL(s.size(), 19914);
  TEST_EXCEPTION(Exception::IllegalArgument, tmp2.getSpectrumByNativeId("TIK"))
}
END_SECTION

START_SECTION((PeakFileOptions& getOptions()))
{
  OnDiscPeakMap tmp;
  PeakFileOptions& options = tmp.getOptions();
  TEST_EQUAL(options.hasMZRange(), false);
  TEST_EQUAL(options.hasRTRange(), false);
  TEST_EQUAL(options.hasIntensityRange(), false);

  // Test that modifications persist
  options.setMZRange(DRange<1>(400.0, 600.0));
  TEST_EQUAL(tmp.getOptions().hasMZRange(), true);
  TEST_REAL_SIMILAR(tmp.getOptions().getMZRange().minPosition()[0], 400.0);
  TEST_REAL_SIMILAR(tmp.getOptions().getMZRange().maxPosition()[0], 600.0);
}
END_SECTION

START_SECTION((const PeakFileOptions& getOptions() const))
{
  OnDiscPeakMap tmp;
  tmp.getOptions().setMZRange(DRange<1>(400.0, 600.0));

  const OnDiscPeakMap& const_tmp = tmp;
  const PeakFileOptions& options = const_tmp.getOptions();
  TEST_EQUAL(options.hasMZRange(), true);
  TEST_REAL_SIMILAR(options.getMZRange().minPosition()[0], 400.0);
}
END_SECTION

START_SECTION((void setOptions(const PeakFileOptions& options)))
{
  OnDiscPeakMap tmp;
  PeakFileOptions options;
  options.setMZRange(DRange<1>(400.0, 600.0));
  options.setIntensityRange(DRange<1>(100.0, 1000.0));

  tmp.setOptions(options);
  TEST_EQUAL(tmp.getOptions().hasMZRange(), true);
  TEST_EQUAL(tmp.getOptions().hasIntensityRange(), true);
  TEST_REAL_SIMILAR(tmp.getOptions().getMZRange().minPosition()[0], 400.0);
  TEST_REAL_SIMILAR(tmp.getOptions().getIntensityRange().minPosition()[0], 100.0);
}
END_SECTION

START_SECTION(([EXTRA] Test m/z range filtering on getSpectrum))
{
  OnDiscPeakMap tmp;
  tmp.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));

  // Get unfiltered spectrum first
  MSSpectrum s_unfiltered = tmp.getSpectrum(0);
  TEST_EQUAL(s_unfiltered.size(), 19914);

  // Now set m/z range filter and verify filtering works
  tmp.getOptions().setMZRange(DRange<1>(400.0, 600.0));
  MSSpectrum s_filtered = tmp.getSpectrum(0);

  // Filtered spectrum should be smaller
  TEST_EQUAL(s_filtered.size() < s_unfiltered.size(), true);
  TEST_EQUAL(s_filtered.size() > 0, true);

  // All peaks should be within range
  for (const auto& peak : s_filtered)
  {
    TEST_EQUAL(peak.getMZ() >= 400.0, true);
    TEST_EQUAL(peak.getMZ() <= 600.0, true);
  }
}
END_SECTION

START_SECTION(([EXTRA] Test intensity range filtering on getSpectrum))
{
  OnDiscPeakMap tmp;
  tmp.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));

  // Get unfiltered spectrum first
  MSSpectrum s_unfiltered = tmp.getSpectrum(0);
  Size unfiltered_size = s_unfiltered.size();

  // Set intensity range filter
  tmp.getOptions().setIntensityRange(DRange<1>(1000.0, 1000000.0));
  MSSpectrum s_filtered = tmp.getSpectrum(0);

  // Filtered spectrum should be smaller (fewer peaks above 1000 intensity)
  TEST_EQUAL(s_filtered.size() < unfiltered_size, true);

  // All peaks should be within intensity range
  for (const auto& peak : s_filtered)
  {
    TEST_EQUAL(peak.getIntensity() >= 1000.0, true);
    TEST_EQUAL(peak.getIntensity() <= 1000000.0, true);
  }
}
END_SECTION

START_SECTION(([EXTRA] Test copy constructor copies options))
{
  OnDiscPeakMap tmp;
  tmp.getOptions().setMZRange(DRange<1>(400.0, 600.0));
  tmp.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));

  OnDiscPeakMap tmp2(tmp);
  TEST_EQUAL(tmp2.getOptions().hasMZRange(), true);
  TEST_REAL_SIMILAR(tmp2.getOptions().getMZRange().minPosition()[0], 400.0);
  TEST_REAL_SIMILAR(tmp2.getOptions().getMZRange().maxPosition()[0], 600.0);
}
END_SECTION

START_SECTION(([EXTRA] Test RT range filter skips loading peak data))
{
  OnDiscPeakMap tmp;
  tmp.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));

  // Get first spectrum to know its RT
  MSSpectrum s1 = tmp.getSpectrum(0);
  double rt = s1.getRT();
  TEST_EQUAL(s1.size(), 19914);  // has peaks

  // Now set RT range that excludes this spectrum
  tmp.getOptions().setRTRange(DRange<1>(rt + 1000, rt + 2000));
  MSSpectrum s2 = tmp.getSpectrum(0);

  // Spectrum should have metadata but no peaks (filtered out, peaks not loaded)
  TEST_EQUAL(s2.empty(), true);  // no peaks loaded
  TEST_REAL_SIMILAR(s2.getRT(), rt);  // but metadata is preserved
}
END_SECTION

START_SECTION(([EXTRA] Test MS level filter skips loading peak data))
{
  OnDiscPeakMap tmp;
  tmp.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));

  // Get first spectrum info
  MSSpectrum s1 = tmp.getSpectrum(0);
  int ms_level = s1.getMSLevel();
  TEST_EQUAL(s1.size(), 19914);  // has peaks

  // Set MS level filter to exclude this spectrum's level
  std::vector<Int> levels;
  levels.push_back(ms_level + 1);  // filter for a different MS level
  tmp.getOptions().setMSLevels(levels);
  MSSpectrum s2 = tmp.getSpectrum(0);

  // Spectrum should have metadata but no peaks (filtered out, peaks not loaded)
  TEST_EQUAL(s2.empty(), true);  // no peaks loaded
  TEST_EQUAL(s2.getMSLevel(), ms_level);  // but metadata is preserved
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

