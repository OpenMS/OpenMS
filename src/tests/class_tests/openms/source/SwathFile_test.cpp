// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/SystemSettings.h>

///////////////////////////
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/CommonEnums.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>

#include <cmath>
#include <utility>


using namespace OpenMS;

bool sortSwathMaps(const OpenSwath::SwathMap& left, const OpenSwath::SwathMap& right)
{
  // true if left is smaller
  if (left.ms1) return true;
  if (right.ms1) return false;
  return left.lower < right.lower;
}


void storeSwathFile(std::string filename, int nr_swathes=32)
{
  PeakMap exp;
  {
    MSSpectrum s;
    s.setMSLevel(1);
    Peak1D p; p.setMZ(101); p.setIntensity(201);
    s.push_back(p);
    exp.addSpectrum(s);
  }
  for (int i = 0; i< nr_swathes; i++)
  {
    MSSpectrum s;
    s.setMSLevel(2);
    std::vector<Precursor> prec(1);
    prec[0].setIsolationWindowLowerOffset(12.5);
    prec[0].setIsolationWindowUpperOffset(12.5);
    prec[0].setMZ(400 + i*25 + 12.5);
    s.setPrecursors(prec);
    Peak1D p; p.setMZ(101 + i); p.setIntensity(201 + i);
    s.push_back(p);
    exp.addSpectrum(s);
  }
  MzMLFile().store(filename, exp);
}

void storeSplitSwathFile(std::vector<std::string> filenames)
{
  {
    PeakMap exp;
    MSSpectrum s;
    s.setMSLevel(1);
    Peak1D p; p.setMZ(101); p.setIntensity(201);
    s.push_back(p);
    exp.addSpectrum(s);
    MzMLFile().store(filenames[0], exp);
  }
  for (Size i = 0; i< filenames.size() -1; i++)
  {
    PeakMap exp;
    MSSpectrum s;
    s.setMSLevel(2);
    std::vector<Precursor> prec(1);
    prec[0].setIsolationWindowLowerOffset(12.5);
    prec[0].setIsolationWindowUpperOffset(12.5);
    prec[0].setMZ(400 + i*25 + 12.5);
    s.setPrecursors(prec);
    Peak1D p; p.setMZ(101 + i); p.setIntensity(201 + i);
    s.push_back(p);
    exp.addSpectrum(s);
    MzMLFile().store(filenames[i+1], exp);
  }
}

namespace
{
  MSSpectrum makeFAIMSSpectrum(int ms_level,
                               double rt,
                               double peak_mz,
                               float peak_intensity,
                               double faims_cv,
                               bool annotate_cv,
                               double isolation_center = 0.0,
                               double isolation_lower = 0.0,
                               double isolation_upper = 0.0)
  {
    MSSpectrum spectrum;
    spectrum.setMSLevel(ms_level);
    spectrum.setRT(rt);
    Peak1D peak;
    peak.setMZ(peak_mz);
    peak.setIntensity(peak_intensity);
    spectrum.push_back(peak);

    if (annotate_cv)
    {
      spectrum.setDriftTime(faims_cv);
      spectrum.setDriftTimeUnit(DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
    }

    if (ms_level > 1)
    {
      Precursor precursor;
      precursor.setMZ(isolation_center);
      precursor.setIsolationWindowLowerOffset(isolation_lower);
      precursor.setIsolationWindowUpperOffset(isolation_upper);
      spectrum.setPrecursors({precursor});
    }

    return spectrum;
  }

  PeakMap makeTwoWindowDIAExperiment(bool faims,
                                    double faims_cv,
                                    double intensity_offset = 0.0,
                                    bool omit_ms2_cv = false,
                                    double rt_offset = 0.0)
  {
    PeakMap exp;

    exp.addSpectrum(makeFAIMSSpectrum(
      1,
      rt_offset,
      100.0,
      static_cast<float>(1000.0 + intensity_offset),
      faims_cv,
      faims));

    exp.addSpectrum(makeFAIMSSpectrum(
      2,
      rt_offset + 1.0,
      110.0,
      static_cast<float>(1100.0 + intensity_offset),
      faims_cv,
      faims && !omit_ms2_cv,
      412.5,
      12.5,
      12.5));

    exp.addSpectrum(makeFAIMSSpectrum(
      2,
      rt_offset + 2.0,
      120.0,
      static_cast<float>(1200.0 + intensity_offset),
      faims_cv,
      faims,
      437.5,
      12.5,
      12.5));

    return exp;
  }
}

START_TEST(SwathFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SwathFile* swath_file_ptr = nullptr;
SwathFile* swath_file_nullPointer = nullptr;

START_SECTION((SwathFile()))
  swath_file_ptr = new SwathFile;
  TEST_NOT_EQUAL(swath_file_ptr, swath_file_nullPointer)
END_SECTION

START_SECTION(([EXTRA]virtual ~SwathFile()))
    delete swath_file_ptr;
END_SECTION

// fast
START_SECTION(std::vector< OpenSwath::SwathMap > loadMzML(std::string file, std::string tmp, std::shared_ptr<ExperimentalSettings>& exp_meta, std::string readoptions="normal") )
{
  Size nr_swathes = 6;
  storeSwathFile("swathFile_1.tmp", nr_swathes);
  std::shared_ptr<ExperimentalSettings> meta = std::shared_ptr<ExperimentalSettings>(new ExperimentalSettings());
  std::vector< OpenSwath::SwathMap > maps = SwathFile().loadMzML("swathFile_1.tmp", SystemSettings::getTempDirectory() + "/", meta);

  TEST_EQUAL(maps.size(), nr_swathes+1)
  TEST_EQUAL(maps[0].ms1, true)
  for (Size i = 0; i< nr_swathes; i++)
  {
    TEST_EQUAL(maps[i+1].ms1, false)
    TEST_EQUAL(maps[i+1].sptr->getNrSpectra(), 1)
    TEST_EQUAL(maps[i+1].sptr->getSpectrumById(0)->getMZArray()->data.size(), 1)
    TEST_REAL_SIMILAR(maps[i+1].sptr->getSpectrumById(0)->getMZArray()->data[0], 101.0+i)
    TEST_REAL_SIMILAR(maps[i+1].sptr->getSpectrumById(0)->getIntensityArray()->data[0], 201.0+i)
    TEST_REAL_SIMILAR(maps[i+1].lower, 400+i*25.0)
    TEST_REAL_SIMILAR(maps[i+1].upper, 425+i*25.0)
  }
}
END_SECTION

// medium (2x slower than normal mzML)
START_SECTION([EXTRA]std::vector< OpenSwath::SwathMap > loadMzML(std::string file, std::string tmp, std::shared_ptr<ExperimentalSettings>& exp_meta, std::string readoptions="cache") )
{
  Size nr_swathes = 2;
  storeSwathFile("swathFile_1.tmp", nr_swathes);
  std::shared_ptr<ExperimentalSettings> meta = std::shared_ptr<ExperimentalSettings>(new ExperimentalSettings());
  std::vector< OpenSwath::SwathMap > maps = SwathFile().loadMzML("swathFile_1.tmp", SystemSettings::getTempDirectory() + "/", meta, "cache");

  TEST_EQUAL(maps.size(), nr_swathes+1)
  TEST_EQUAL(maps[0].ms1, true)
  for (Size i = 0; i< nr_swathes; i++)
  {
    TEST_EQUAL(maps[i+1].ms1, false)
    TEST_EQUAL(maps[i+1].sptr->getNrSpectra(), 1)
    TEST_EQUAL(maps[i+1].sptr->getSpectrumById(0)->getMZArray()->data.size(), 1)
    TEST_REAL_SIMILAR(maps[i+1].sptr->getSpectrumById(0)->getMZArray()->data[0], 101.0+i)
    TEST_REAL_SIMILAR(maps[i+1].sptr->getSpectrumById(0)->getIntensityArray()->data[0], 201.0+i)
    TEST_REAL_SIMILAR(maps[i+1].lower, 400+i*25.0)
    TEST_REAL_SIMILAR(maps[i+1].upper, 425+i*25.0)
  }
}
END_SECTION

// medium (2x slower than normal mzML)
START_SECTION(std::vector< OpenSwath::SwathMap > loadSplit(StringList file_list, std::string tmp, std::shared_ptr<ExperimentalSettings>& exp_meta, std::string readoptions="normal"))
{
  std::vector<std::string> swath_filenames;
  Size nr_swathes = 3;
  swath_filenames.push_back("swathFile_2_ms1.tmp");
  for (Size i = 0; i < nr_swathes; i++)
  {
    swath_filenames.push_back("swathFile_2_sw"  + StringUtils::toStr(i) + ".tmp");
  }
  storeSplitSwathFile(swath_filenames);
  std::shared_ptr<ExperimentalSettings> meta = std::shared_ptr<ExperimentalSettings>(new ExperimentalSettings());
  std::vector< OpenSwath::SwathMap > maps = SwathFile().loadSplit(swath_filenames, SystemSettings::getTempDirectory() + "/", meta);

  // ensure they are sorted ... 
  std::sort(maps.begin(), maps.end(), sortSwathMaps);

  TEST_EQUAL(maps.size(), nr_swathes + 1)
  TEST_EQUAL(maps[0].ms1, true)
  for (Size i = 0; i< maps.size() -1; i++)
  {
    TEST_EQUAL(maps[i+1].ms1, false)
    TEST_EQUAL(maps[i+1].sptr->getNrSpectra(), 1)
    TEST_EQUAL(maps[i+1].sptr->getSpectrumById(0)->getMZArray()->data.size(), 1)
    TEST_REAL_SIMILAR(maps[i+1].sptr->getSpectrumById(0)->getMZArray()->data[0], 101.0+i)
    TEST_REAL_SIMILAR(maps[i+1].sptr->getSpectrumById(0)->getIntensityArray()->data[0], 201.0+i)
    TEST_REAL_SIMILAR(maps[i+1].lower, 400+i*25.0)
    TEST_REAL_SIMILAR(maps[i+1].upper, 425+i*25.0)
  }

}
END_SECTION

// slow (7x slower than normal mzML)
START_SECTION([EXTRA]std::vector< OpenSwath::SwathMap > loadSplit(StringList file_list, std::string tmp, std::shared_ptr<ExperimentalSettings>& exp_meta, std::string readoptions="cache"))
{
  std::vector<std::string> swath_filenames;
  Size nr_swathes = 2;
  swath_filenames.push_back("swathFile_3_ms1.tmp");
  for (Size i = 0; i < nr_swathes; i++)
  {
    swath_filenames.push_back("swathFile_3_sw"  + StringUtils::toStr(i) + ".tmp");
  }
  storeSplitSwathFile(swath_filenames);
  std::shared_ptr<ExperimentalSettings> meta = std::shared_ptr<ExperimentalSettings>(new ExperimentalSettings());
  std::vector< OpenSwath::SwathMap > maps = SwathFile().loadSplit(swath_filenames, SystemSettings::getTempDirectory() + "/", meta, "cache");
  // ensure they are sorted ... 
  std::sort(maps.begin(), maps.end(), sortSwathMaps);

  TEST_EQUAL(maps.size(), nr_swathes + 1)
  TEST_EQUAL(maps[0].ms1, true)
  for (Size i = 0; i< maps.size() -1; i++)
  {
    TEST_EQUAL(maps[i+1].ms1, false)
    TEST_EQUAL(maps[i+1].sptr->getNrSpectra(), 1)
    TEST_EQUAL(maps[i+1].sptr->getSpectrumById(0)->getMZArray()->data.size(), 1)
    TEST_REAL_SIMILAR(maps[i+1].sptr->getSpectrumById(0)->getMZArray()->data[0], 101.0+i)
    TEST_REAL_SIMILAR(maps[i+1].sptr->getSpectrumById(0)->getIntensityArray()->data[0], 201.0+i)
    TEST_REAL_SIMILAR(maps[i+1].lower, 400+i*25.0)
    TEST_REAL_SIMILAR(maps[i+1].upper, 425+i*25.0)
  }

}
END_SECTION

START_SECTION((std::vector<FAIMSSwathMapGroup> loadFromMSExperimentByFAIMSCV(PeakMap&& exp, const std::string& tmp, std::shared_ptr<ExperimentalSettings>& exp_meta, const std::string& readoptions) - non-FAIMS))
{
  PeakMap exp = makeTwoWindowDIAExperiment(false, 0.0);
  std::shared_ptr<ExperimentalSettings> meta;

  auto groups = SwathFile().loadFromMSExperimentByFAIMSCV(
    std::move(exp), SystemSettings::getTempDirectory() + "/", meta, "normal");

  TEST_EQUAL(groups.size(), 1)
  TEST_EQUAL(std::isnan(groups[0].faims_cv), true)
  TEST_EQUAL(groups[0].swath_maps.size(), 3)
  TEST_EQUAL(groups[0].swath_maps[0].ms1, true)
  TEST_REAL_SIMILAR(groups[0].swath_maps[1].lower, 400.0)
  TEST_REAL_SIMILAR(groups[0].swath_maps[1].upper, 425.0)
  TEST_REAL_SIMILAR(groups[0].swath_maps[2].lower, 425.0)
  TEST_REAL_SIMILAR(groups[0].swath_maps[2].upper, 450.0)
}
END_SECTION

START_SECTION((std::vector<FAIMSSwathMapGroup> loadFromMSExperimentByFAIMSCV(PeakMap&& exp, const std::string& tmp, std::shared_ptr<ExperimentalSettings>& exp_meta, const std::string& readoptions) - single CV))
{
  // The first MS2 spectrum deliberately lacks explicit FAIMS metadata. It should inherit -45 V
  // from the preceding MS1 spectrum through the existing splitByFAIMSCV() policy.
  PeakMap exp = makeTwoWindowDIAExperiment(true, -45.0, 0.0, true);
  std::shared_ptr<ExperimentalSettings> meta;

  auto groups = SwathFile().loadFromMSExperimentByFAIMSCV(
    std::move(exp), SystemSettings::getTempDirectory() + "/", meta, "normal");

  TEST_EQUAL(groups.size(), 1)
  TEST_REAL_SIMILAR(groups[0].faims_cv, -45.0)
  TEST_EQUAL(groups[0].swath_maps.size(), 3)
  TEST_EQUAL(groups[0].swath_maps[0].sptr->getNrSpectra(), 1)
  TEST_EQUAL(groups[0].swath_maps[1].sptr->getNrSpectra(), 1)
  TEST_EQUAL(groups[0].swath_maps[2].sptr->getNrSpectra(), 1)
  TEST_REAL_SIMILAR(groups[0].swath_maps[1].lower, 400.0)
  TEST_REAL_SIMILAR(groups[0].swath_maps[1].upper, 425.0)
  TEST_REAL_SIMILAR(groups[0].swath_maps[2].lower, 425.0)
  TEST_REAL_SIMILAR(groups[0].swath_maps[2].upper, 450.0)
}
END_SECTION

START_SECTION((std::vector<FAIMSSwathMapGroup> loadFromMSExperimentByFAIMSCV(PeakMap&& exp, const std::string& tmp, std::shared_ptr<ExperimentalSettings>& exp_meta, const std::string& readoptions) - multiple CVs))
{
  // First complete DIA cycle.
  PeakMap minus55 = makeTwoWindowDIAExperiment(true, -55.0, 5500.0, true, 0.0);

  // Following DIA cycle at the second CV.
  PeakMap minus45 = makeTwoWindowDIAExperiment(true, -45.0, 4500.0, true, 3.0);

  // Preserve acquisition order: one complete DIA cycle at -55 V followed by one at -45 V.
  PeakMap exp;
  for (auto& spectrum : minus55.getSpectra()) exp.addSpectrum(std::move(spectrum));
  for (auto& spectrum : minus45.getSpectra()) exp.addSpectrum(std::move(spectrum));

  std::shared_ptr<ExperimentalSettings> meta;
  auto groups = SwathFile().loadFromMSExperimentByFAIMSCV(
    std::move(exp), SystemSettings::getTempDirectory() + "/", meta, "normal");

  TEST_EQUAL(groups.size(), 2)
  TEST_REAL_SIMILAR(groups[0].faims_cv, -55.0)
  TEST_REAL_SIMILAR(groups[1].faims_cv, -45.0)

  for (const auto& group : groups)
  {
    TEST_EQUAL(group.swath_maps.size(), 3)
    TEST_EQUAL(group.swath_maps[0].ms1, true)
    TEST_REAL_SIMILAR(group.swath_maps[1].lower, 400.0)
    TEST_REAL_SIMILAR(group.swath_maps[1].upper, 425.0)
    TEST_REAL_SIMILAR(group.swath_maps[2].lower, 425.0)
    TEST_REAL_SIMILAR(group.swath_maps[2].upper, 450.0)
  }

  // Peak intensities encode the source CV group, so this verifies that spectra were not mixed
  // between CV buckets during SWATH conversion.
  TEST_REAL_SIMILAR(groups[0].swath_maps[1].sptr->getSpectrumById(0)->getIntensityArray()->data[0], 6600.0)
  TEST_REAL_SIMILAR(groups[1].swath_maps[1].sptr->getSpectrumById(0)->getIntensityArray()->data[0], 5600.0)
}
END_SECTION

START_SECTION((std::vector<FAIMSSwathMapGroup> loadOnDiscByFAIMSCV(const std::shared_ptr<IMSExperimentRandomAccess>& experiment, std::shared_ptr<ExperimentalSettings>& exp_meta, Interfaces::IMSDataConsumer* plugin_consumer)))
{
  // Two CVs, MS1 and MS2 at each; the first MS2 of each cycle has no CV and goes to the CV before it
  PeakMap minus55 = makeTwoWindowDIAExperiment(true, -55.0, 5500.0, true, 0.0);
  PeakMap minus45 = makeTwoWindowDIAExperiment(true, -45.0, 4500.0, true, 3.0);
  PeakMap exp;
  for (auto& spectrum : minus55.getSpectra()) exp.addSpectrum(spectrum);
  for (auto& spectrum : minus45.getSpectra()) exp.addSpectrum(spectrum);
  for (auto& spectrum : minus55.getSpectra()) // a third cycle back at -55 V
  {
    MSSpectrum s = spectrum;
    s.setRT(s.getRT() + 6.0);
    exp.addSpectrum(s);
  }
  // mzML stores the drift time only with a scan (acquisition) entry
  for (auto& spectrum : exp)
  {
    spectrum.getAcquisitionInfo().push_back(Acquisition());
  }
  std::string file;
  NEW_TMP_FILE(file)
  MzMLFile().store(file, exp);

  auto on_disc = std::make_shared<OnDiscMSExperiment>();
  TEST_TRUE(on_disc->openFile(file))
  std::shared_ptr<ExperimentalSettings> meta;
  auto groups = SwathFile().loadOnDiscByFAIMSCV(on_disc, meta);

  TEST_EQUAL(groups.size(), 2)
  ABORT_IF(groups.size() != 2)
  TEST_REAL_SIMILAR(groups[0].faims_cv, -55.0)
  TEST_REAL_SIMILAR(groups[1].faims_cv, -45.0)
  const std::vector<Size> nr_spectra_per_cv = {2, 1};
  for (Size g = 0; g < groups.size(); ++g)
  {
    const auto& maps = groups[g].swath_maps;
    TEST_EQUAL(maps.size(), 3)
    ABORT_IF(maps.size() != 3)
    TEST_EQUAL(maps[0].ms1, true)
    TEST_EQUAL(maps[1].ms1, false)
    TEST_REAL_SIMILAR(maps[1].lower, 400.0)
    TEST_REAL_SIMILAR(maps[1].upper, 425.0)
    TEST_REAL_SIMILAR(maps[2].lower, 425.0)
    TEST_REAL_SIMILAR(maps[2].upper, 450.0)
    for (const auto& map : maps)
    {
      TEST_EQUAL(map.sptr->getNrSpectra(), nr_spectra_per_cv[g])
    }
  }
  // Intensities encode the CV of the cycle: no spectra are mixed between the CVs
  for (Size i = 0; i < 2; ++i)
  {
    TEST_REAL_SIMILAR(groups[0].swath_maps[0].sptr->getSpectrumById(int(i))->getIntensityArray()->data[0], 6500.0)
    TEST_REAL_SIMILAR(groups[0].swath_maps[1].sptr->getSpectrumById(int(i))->getIntensityArray()->data[0], 6600.0)
    TEST_REAL_SIMILAR(groups[0].swath_maps[2].sptr->getSpectrumById(int(i))->getIntensityArray()->data[0], 6700.0)
  }
  TEST_REAL_SIMILAR(groups[1].swath_maps[0].sptr->getSpectrumById(0)->getIntensityArray()->data[0], 5500.0)
  TEST_REAL_SIMILAR(groups[1].swath_maps[1].sptr->getSpectrumById(0)->getIntensityArray()->data[0], 5600.0)
  TEST_REAL_SIMILAR(groups[1].swath_maps[2].sptr->getSpectrumById(0)->getIntensityArray()->data[0], 5700.0)

  // Same groups as the in-memory FAIMS adapter
  std::shared_ptr<ExperimentalSettings> meta_in_memory;
  auto in_memory = SwathFile().loadFromMSExperimentByFAIMSCV(PeakMap(exp), SystemSettings::getTempDirectory() + "/", meta_in_memory, "normal");
  TEST_EQUAL(in_memory.size(), groups.size())
  for (Size g = 0; g < std::min(in_memory.size(), groups.size()); ++g)
  {
    TEST_REAL_SIMILAR(in_memory[g].faims_cv, groups[g].faims_cv)
    TEST_EQUAL(in_memory[g].swath_maps.size(), groups[g].swath_maps.size())
    for (Size m = 0; m < std::min(in_memory[g].swath_maps.size(), groups[g].swath_maps.size()); ++m)
    {
      TEST_EQUAL(in_memory[g].swath_maps[m].sptr->getNrSpectra(), groups[g].swath_maps[m].sptr->getNrSpectra())
    }
  }

  // The flat map list has no CV: multi-CV input is rejected instead of merging the CVs
  std::shared_ptr<ExperimentalSettings> meta_flat;
  TEST_EXCEPTION(Exception::IllegalArgument, SwathFile().loadOnDisc(on_disc, meta_flat))

  // Non-FAIMS input: one group with CV NaN, same maps as loadOnDisc()
  std::string file_no_faims;
  NEW_TMP_FILE(file_no_faims)
  MzMLFile().store(file_no_faims, makeTwoWindowDIAExperiment(false, 0.0));
  auto on_disc_no_faims = std::make_shared<OnDiscMSExperiment>();
  TEST_TRUE(on_disc_no_faims->openFile(file_no_faims))
  auto groups_no_faims = SwathFile().loadOnDiscByFAIMSCV(on_disc_no_faims, meta);
  TEST_EQUAL(groups_no_faims.size(), 1)
  TEST_TRUE(std::isnan(groups_no_faims[0].faims_cv))
  TEST_EQUAL(groups_no_faims[0].swath_maps.size(), 3)
  TEST_EQUAL(SwathFile().loadOnDisc(on_disc_no_faims, meta).size(), 3)
}
END_SECTION

START_SECTION(([EXTRA] plugin consumers with readoptions normal, cache and onDisc))
{
  Size nr_swathes = 2;
  std::string file;
  NEW_TMP_FILE(file)
  storeSwathFile(file, nr_swathes);
  const std::string tmp = SystemSettings::getTempDirectory() + "/";

  MSDataTransformingConsumer scaling;
  scaling.setSpectraProcessingFunc([](MSSpectrum& s) { for (auto& p : s) { p.setIntensity(p.getIntensity() * 2); } });

  // normal and cache keep the changes of a transforming plugin (MS2 intensity 201 -> 402)
  for (const std::string readoptions : {"normal", "cache"})
  {
    std::shared_ptr<ExperimentalSettings> meta;
    auto maps = SwathFile().loadMzML(file, tmp, meta, readoptions, &scaling);
    TEST_EQUAL(maps.size(), nr_swathes + 1)
    TEST_REAL_SIMILAR(maps[1].sptr->getSpectrumById(0)->getIntensityArray()->data[0], 402.0)
  }

  // onDisc reads the spectra from the file again: a transforming plugin is rejected, not ignored
  auto on_disc = std::make_shared<OnDiscMSExperiment>();
  TEST_TRUE(on_disc->openFile(file))
  std::shared_ptr<ExperimentalSettings> meta;
  TEST_EXCEPTION(Exception::IllegalArgument, SwathFile().loadOnDisc(on_disc, meta, &scaling))
  TEST_EXCEPTION(Exception::IllegalArgument, SwathFile().loadOnDiscByFAIMSCV(on_disc, meta, &scaling))

  // An observing plugin (like SwathQC) sees every spectrum once and the maps equal those of normal
  Size nr_seen = 0;
  MSDataTransformingConsumer observer;
  observer.setSpectraProcessingFunc([&nr_seen](MSSpectrum&) { ++nr_seen; });
  auto maps_on_disc = SwathFile().loadOnDisc(on_disc, meta, &observer);
  TEST_EQUAL(nr_seen, on_disc->getNrSpectra())
  auto maps_normal = SwathFile().loadMzML(file, tmp, meta, "normal");
  TEST_EQUAL(maps_on_disc.size(), maps_normal.size())
  for (Size m = 0; m < std::min(maps_on_disc.size(), maps_normal.size()); ++m)
  {
    TEST_EQUAL(maps_on_disc[m].ms1, maps_normal[m].ms1)
    TEST_REAL_SIMILAR(maps_on_disc[m].lower, maps_normal[m].lower)
    TEST_REAL_SIMILAR(maps_on_disc[m].upper, maps_normal[m].upper)
    TEST_EQUAL(maps_on_disc[m].sptr->getNrSpectra(), maps_normal[m].sptr->getNrSpectra())
    TEST_EQUAL(maps_on_disc[m].sptr->getSpectrumById(0)->getIntensityArray()->data == maps_normal[m].sptr->getSpectrumById(0)->getIntensityArray()->data, true)
    TEST_EQUAL(maps_on_disc[m].sptr->getSpectrumById(0)->getMZArray()->data == maps_normal[m].sptr->getSpectrumById(0)->getMZArray()->data, true)
  }
}
END_SECTION

START_SECTION((std::vector< OpenSwath::SwathMap > loadMzXML(std::string file, std::string tmp, std::shared_ptr<ExperimentalSettings>& exp_meta, std::string readoptions="normal") ) )
{
  NOT_TESTABLE // mzXML is not supported
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
