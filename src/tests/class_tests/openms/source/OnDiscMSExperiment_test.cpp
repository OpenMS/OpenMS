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
#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>
#include <OpenMS/KERNEL/MSExperiment.h>
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

// START_SECTION((OnDiscMSExperiment(const std::string& filename)))
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

START_SECTION(( bool openFile(const std::string& filename, bool skipMetaData = false) ))
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

START_SECTION(([EXTRA] Test precursor m/z range filter skips loading peak data for MS2))
{
  OnDiscPeakMap tmp;
  tmp.openFile(OPENMS_GET_TEST_DATA_PATH("MzMLFile_4_indexed.mzML"));

  // File has 4 spectra: MS1 (index 0), MS2 (index 1 with precursor m/z 5.55), MS1 (index 2), MS1 (index 3)
  // Get MS2 spectrum (index 1) without filter
  MSSpectrum s1 = tmp.getSpectrum(1);
  TEST_EQUAL(s1.getMSLevel(), 2);
  Size unfiltered_size = s1.size();
  TEST_EQUAL(unfiltered_size > 0, true);  // Should have peaks

  // Get precursor m/z for verification
  TEST_EQUAL(s1.getPrecursors().empty(), false);
  double precursor_mz = s1.getPrecursors()[0].getMZ();
  TEST_REAL_SIMILAR(precursor_mz, 5.55);  // Verify expected precursor m/z

  // Set precursor m/z range that excludes this spectrum (range: 10.0 - 20.0)
  tmp.getOptions().setPrecursorMZRange(DRange<1>(10.0, 20.0));
  MSSpectrum s2 = tmp.getSpectrum(1);

  // Spectrum should have metadata but no peaks (filtered out by precursor m/z)
  TEST_EQUAL(s2.empty(), true);  // no peaks loaded
  TEST_EQUAL(s2.getMSLevel(), 2);  // but metadata is preserved
  TEST_EQUAL(s2.getPrecursors().empty(), false);  // precursor info is preserved

  // Set precursor m/z range that includes this spectrum (range: 5.0 - 6.0)
  tmp.getOptions().setPrecursorMZRange(DRange<1>(5.0, 6.0));
  MSSpectrum s3 = tmp.getSpectrum(1);

  // Spectrum should have peaks (passes precursor m/z filter)
  TEST_EQUAL(s3.size(), unfiltered_size);  // peaks loaded
  TEST_EQUAL(s3.getMSLevel(), 2);

  // MS1 spectra (index 0, 2, 3) should not be affected by precursor m/z filter
  MSSpectrum s_ms1 = tmp.getSpectrum(0);
  TEST_EQUAL(s_ms1.getMSLevel(), 1);
  TEST_EQUAL(s_ms1.size() > 0, true);  // MS1 should still have peaks despite precursor filter
}
END_SECTION

START_SECTION(([EXTRA] Range filters keep spectrum metadata and filter data arrays in step with the peaks))
{
  // indexed mzML with an MS2 spectrum carrying float, string and integer data arrays and a chromatogram with arrays
  PeakMap exp;
  MSSpectrum spec;
  spec.setRT(123.5);
  spec.setMSLevel(2);
  spec.setName("my_spectrum");
  spec.setNativeID("scan=1");
  Precursor prec;
  prec.setMZ(500.0);
  spec.getPrecursors().push_back(prec);
  MSSpectrum::FloatDataArray fda;
  fda.setName("Ion Mobility");
  MSSpectrum::StringDataArray sda;
  sda.setName("labels");
  MSSpectrum::IntegerDataArray ida;
  ida.setName("indices");
  for (Size i = 0; i < 10; ++i)
  {
    spec.emplace_back(100.0 + 10.0 * i, 10.0 * (i + 1));   // mz 100..190, intensity 10..100
    fda.push_back(0.5 + i);
    sda.push_back("p" + std::to_string(i));
    ida.push_back(Int64(1000 + i));
  }
  spec.getFloatDataArrays().push_back(fda);
  spec.getStringDataArrays().push_back(sda);
  spec.getIntegerDataArrays().push_back(ida);
  exp.addSpectrum(spec);

  MSChromatogram chrom;
  chrom.setName("my_chrom");
  chrom.setNativeID("chrom=1");
  MSChromatogram::FloatDataArray cfda;
  cfda.setName("chrom_array");
  for (Size i = 0; i < 10; ++i)
  {
    chrom.emplace_back(1.0 + i, 10.0 * (i + 1));   // rt 1..10, intensity 10..100
    cfda.push_back(0.25 + i);
  }
  chrom.getFloatDataArrays().push_back(cfda);
  exp.addChromatogram(chrom);

  std::string filename;
  NEW_TMP_FILE(filename);
  IndexedMzMLFileLoader().store(filename, exp);

  OnDiscPeakMap od;
  TEST_EQUAL(od.openFile(filename), true);
  MSSpectrum ref = od.getSpectrum(0);
  TEST_EQUAL(ref.size(), 10);
  // The unfiltered indexed read is the reference. It does not round-trip the spectrum name and returns
  // string data arrays without content (the on-disc decoder does not fill them), so those two are compared
  // against the reference instead of the originally stored values.
  const bool ref_has_strings = !ref.getStringDataArrays().empty() && !ref.getStringDataArrays()[0].empty();

  // checks shared by all range types; kept = indices of the peaks that must survive
  auto check = [&](const MSSpectrum& s, const std::vector<Size>& kept)
  {
    TEST_REAL_SIMILAR(s.getRT(), 123.5);
    TEST_EQUAL(s.getMSLevel(), 2);
    TEST_EQUAL(s.getName(), ref.getName());
    TEST_EQUAL(s.getNativeID(), "scan=1");
    TEST_EQUAL(s.getPrecursors().size(), 1);
    TEST_EQUAL(s.size(), kept.size());
    TEST_EQUAL(s.getFloatDataArrays().size(), 1);
    TEST_EQUAL(s.getStringDataArrays().size(), 1);
    TEST_EQUAL(s.getIntegerDataArrays().size(), 1);
    if (s.getFloatDataArrays().size() != 1 || s.getStringDataArrays().size() != 1 || s.getIntegerDataArrays().size() != 1) return;
    TEST_EQUAL(s.getFloatDataArrays()[0].getName(), "Ion Mobility");
    TEST_EQUAL(s.getStringDataArrays()[0].getName(), "labels");
    TEST_EQUAL(s.getIntegerDataArrays()[0].getName(), "indices");
    TEST_EQUAL(s.getFloatDataArrays()[0].size(), kept.size());
    TEST_EQUAL(s.getStringDataArrays()[0].size(), ref_has_strings ? kept.size() : 0);
    TEST_EQUAL(s.getIntegerDataArrays()[0].size(), kept.size());
    if (s.getFloatDataArrays()[0].size() != kept.size()
        || s.getIntegerDataArrays()[0].size() != kept.size()
        || s.size() != kept.size()
        || (ref_has_strings && s.getStringDataArrays()[0].size() != kept.size())) return;
    for (Size j = 0; j < kept.size(); ++j)
    {
      TEST_REAL_SIMILAR(s[j].getMZ(), ref[kept[j]].getMZ());
      TEST_REAL_SIMILAR(s[j].getIntensity(), ref[kept[j]].getIntensity());
      TEST_REAL_SIMILAR(s.getFloatDataArrays()[0][j], ref.getFloatDataArrays()[0][kept[j]]);
      if (ref_has_strings) TEST_EQUAL(s.getStringDataArrays()[0][j], ref.getStringDataArrays()[0][kept[j]]);
      TEST_EQUAL(s.getIntegerDataArrays()[0][j], ref.getIntegerDataArrays()[0][kept[j]]);
    }
  };

  od.getOptions().setMZRange(DRange<1>(125.0, 165.0));   // keeps mz 130..160 -> indices 3..6
  MSSpectrum s_mz = od.getSpectrum(0);
  check(s_mz, {3, 4, 5, 6});
  // the cached ranges describe the kept peaks, not the unfiltered spectrum
  TEST_REAL_SIMILAR(s_mz.getMinMZ(), 130.0);
  TEST_REAL_SIMILAR(s_mz.getMaxMZ(), 160.0);
  TEST_REAL_SIMILAR(s_mz.getMinIntensity(), 40.0);
  TEST_REAL_SIMILAR(s_mz.getMaxIntensity(), 70.0);

  od.getOptions() = PeakFileOptions();
  od.getOptions().setIntensityRange(DRange<1>(25.0, 65.0));   // keeps intensity 30..60 -> indices 2..5
  check(od.getSpectrum(0), {2, 3, 4, 5});

  // both filters at once; only their intersection survives (mz alone keeps 3, 4; intensity alone keeps 4, 5)
  od.getOptions().setIntensityRange(DRange<1>(45.0, 65.0));   // intensity 50, 60 -> indices 4, 5
  od.getOptions().setMZRange(DRange<1>(125.0, 145.0));   // mz 130, 140 -> indices 3, 4
  check(od.getSpectrum(0), {4});

  // chromatogram: RT range and intensity range keep name, native ID and data arrays in step with the peaks
  od.getOptions() = PeakFileOptions();
  MSChromatogram cref = od.getChromatogram(0);
  TEST_EQUAL(cref.size(), 10);
  auto checkChrom = [&](const MSChromatogram& c, const std::vector<Size>& kept)
  {
    TEST_EQUAL(c.getName(), "my_chrom");
    TEST_EQUAL(c.getNativeID(), "chrom=1");
    TEST_EQUAL(c.size(), kept.size());
    TEST_EQUAL(c.getFloatDataArrays().size(), 1);
    if (c.getFloatDataArrays().size() != 1 || c.size() != kept.size()) return;
    TEST_EQUAL(c.getFloatDataArrays()[0].getName(), "chrom_array");
    TEST_EQUAL(c.getFloatDataArrays()[0].size(), kept.size());
    if (c.getFloatDataArrays()[0].size() != kept.size()) return;
    for (Size j = 0; j < kept.size(); ++j)
    {
      TEST_REAL_SIMILAR(c[j].getRT(), cref[kept[j]].getRT());
      TEST_REAL_SIMILAR(c[j].getIntensity(), cref[kept[j]].getIntensity());
      TEST_REAL_SIMILAR(c.getFloatDataArrays()[0][j], cref.getFloatDataArrays()[0][kept[j]]);
    }
  };
  od.getOptions().setRTRange(DRange<1>(3.5, 6.5));   // rt 4..6 -> indices 3..5
  MSChromatogram c_rt = od.getChromatogram(0);
  checkChrom(c_rt, {3, 4, 5});
  TEST_REAL_SIMILAR(c_rt.getMinRT(), 4.0);
  TEST_REAL_SIMILAR(c_rt.getMaxRT(), 6.0);
  od.getOptions() = PeakFileOptions();
  od.getOptions().setIntensityRange(DRange<1>(25.0, 55.0));   // intensity 30..50 -> indices 2..4
  checkChrom(od.getChromatogram(0), {2, 3, 4});
}
END_SECTION

START_SECTION(([EXTRA] Range filter on a spectrum with a data array of the wrong length does not throw))
{
  // a data array shorter than the peak list cannot be aligned with the kept peaks; it is emptied instead of throwing
  PeakMap exp;
  MSSpectrum spec;
  spec.setRT(10.0);
  for (Size i = 0; i < 10; ++i) spec.emplace_back(100.0 + 10.0 * i, 10.0 * (i + 1));
  MSSpectrum::FloatDataArray fda;
  fda.setName("short array");
  fda.push_back(1.0);
  fda.push_back(2.0);
  spec.getFloatDataArrays().push_back(fda);
  exp.addSpectrum(spec);

  std::string filename;
  NEW_TMP_FILE(filename);
  IndexedMzMLFileLoader().store(filename, exp);

  OnDiscPeakMap od;
  TEST_EQUAL(od.openFile(filename), true);
  MSSpectrum ref = od.getSpectrum(0);
  TEST_EQUAL(ref.size(), 10);
  od.getOptions().setMZRange(DRange<1>(125.0, 165.0));
  MSSpectrum s = od.getSpectrum(0);
  TEST_EQUAL(s.size(), 4);
  TEST_REAL_SIMILAR(s.getRT(), 10.0);
  if (!ref.getFloatDataArrays().empty() && ref.getFloatDataArrays()[0].size() != ref.size())
  {
    // the reader kept the mis-sized array: the filtered spectrum must keep the array (and its name) but empty
    TEST_EQUAL(s.getFloatDataArrays().size(), 1);
    if (!s.getFloatDataArrays().empty()) TEST_EQUAL(s.getFloatDataArrays()[0].empty(), true);
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

