// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/config.h>

#ifdef WITH_THERMO_RAW

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/FORMAT/ThermoRawFile.h>
#include <OpenMS/FORMAT/ZipArchiveFile.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/METADATA/IonSource.h>
#include <OpenMS/METADATA/MassAnalyzer.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/SystemSettings.h>
#include <OpenMS/SYSTEM/TempFiles.h>

#include <cmath>
#include <fstream>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <utility>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace OpenMS;
using namespace std;

namespace
{
  /// Collects a streamed run and checks the order of the consumer calls.
  struct RecordingConsumer : Interfaces::IMSDataConsumer
  {
    MSExperiment exp;
    Size expected_spectra = 0, expected_chromatograms = 0;
    int stage = 0; // 1: expected size, 2: settings, 3: spectra, 4: chromatograms
    bool order_ok = true;

    void advance(int to)
    {
      if (to < stage || (stage < 2 && to > 2) || (stage == 0 && to == 2)) { order_ok = false; }
      stage = to;
    }
    void setExpectedSize(Size spectra, Size chromatograms) override
    {
      advance(1);
      expected_spectra = spectra;
      expected_chromatograms = chromatograms;
    }
    void setExperimentalSettings(const ExperimentalSettings& settings) override
    {
      advance(2);
      exp.getExperimentalSettings() = settings;
    }
    void consumeSpectrum(SpectrumType& s) override
    {
      advance(3);
      exp.addSpectrum(s);
    }
    void consumeChromatogram(ChromatogramType& c) override
    {
      advance(4);
      exp.addChromatogram(c);
    }
  };

  /// The message of the Exception::ParseError that @p f throws (empty if it throws none)
  std::string parseErrorOf(const std::function<void()>& f)
  {
    try
    {
      f();
    }
    catch (const Exception::ParseError& e)
    {
      return e.what();
    }
    return std::string();
  }
}

START_TEST(ThermoRawFile, "$Id$")

START_SECTION(void load(const std::string& path, MSExperiment& exp))
{
  ThermoRawFile file;
  MSExperiment exp;

  TEST_EXCEPTION(Exception::FileNotFound, file.load("/definitely/not/present.raw", exp))
  TEST_EXCEPTION(Exception::FileNotFound, file.load("/definitely/not/present.raw.zip", exp))
}
END_SECTION

START_SECTION([EXTRA] compressed RAW files that cannot be read)
{
  // gzip- and bzip2-compressed RAW files, and ZIP archives without exactly one .raw file, are reported as such
  // (not as the vendor library's failure to open the file), by load(), transform() and FileHandler alike
  TempDir tmp;
  const std::string dummy = tmp.getPath() + "/dummy.raw";
  {
    std::ofstream ofs(dummy.c_str(), std::ios::binary);
    ofs << "not a RAW file";
  }
  const std::string gzipped = dummy + ".gz", bzipped = dummy + ".bz2";
  File::copy(dummy, gzipped);
  File::copy(dummy, bzipped);
  const std::string none = tmp.getPath() + "/none.raw.zip";
  ZipArchiveFile::addOrReplaceFromFile(none, "dummy.txt", dummy);
  ZipArchiveFile::addOrReplaceFromFile(none, "__MACOSX/._dummy.raw", dummy);
  const std::string two = tmp.getPath() + "/two.raw.zip";
  ZipArchiveFile::addOrReplaceFromFile(two, "a.raw", dummy);
  ZipArchiveFile::addOrReplaceFromFile(two, "b.RAW", dummy);

  const std::vector<std::pair<std::string, std::string>> cases = {
    {gzipped, "A gzip-compressed Thermo RAW file cannot be read: decompress it first"},
    {bzipped, "A bzip2-compressed Thermo RAW file cannot be read: decompress it first"},
    {none, "ZIP archive holds no file with the extension '.raw'"},
    {two, "ZIP archive holds 2 files with the extension '.raw' (a.raw, b.RAW)"}};
  for (const auto& [path, expected] : cases)
  {
    TEST_EQUAL(FileHandler::getType(path), FileTypes::RAW)
    MSExperiment exp;
    std::string message = parseErrorOf([&] { ThermoRawFile().load(path, exp); });
    TEST_EQUAL(message.find(expected) != std::string::npos, true)
    // reported once, not wrapped as a failure of the vendor library
    TEST_EQUAL(message.find("Thermo RAW reader") == std::string::npos, true)
    RecordingConsumer collected;
    message = parseErrorOf([&] { ThermoRawFile().transform(path, &collected); });
    TEST_EQUAL(message.find(expected) != std::string::npos, true)
    message = parseErrorOf([&] { FileHandler().loadExperiment(path, exp); });
    TEST_EQUAL(message.find(expected) != std::string::npos, true)
  }
}
END_SECTION

#ifdef THERMO_RAW_TEST_DATA
START_SECTION(round-trip load raw -> mzML -> reload MSExperiment)
{
  // Test data is PNNL's Angiotensin_AllScans.raw (BSD-2-Clause, Battelle 2018):
  // 87 MS1 (FTMS Full ms) + hundreds of MS2 scans with HCD (@hcd30) and
  // supplemental-activation filters (ETD+HCD, ETD+CID) at precursor m/z 432.9.
  ThermoRawFile file;
  MSExperiment original;
  file.load(THERMO_RAW_TEST_DATA, original);

  TEST_NOT_EQUAL(original.size(), 0)
  TEST_EQUAL(original.getSourceFiles().empty(), false)
  // Broad shape rather than exact totals — keeps the test robust to
  // bridge-side scan-filter handling.
  TEST_EQUAL(original.size() > 1000, true)

  // MS-level breakdown: file contains both MS1 and MSn data.
  Size ms1_count = 0;
  Size msn_count = 0;
  for (const MSSpectrum& s : original.getSpectra())
  {
    if (s.getMSLevel() == 1) { ++ms1_count; }
    else if (s.getMSLevel() >= 2) { ++msn_count; }
  }
  TEST_EQUAL(ms1_count > 0, true)
  TEST_EQUAL(msn_count > 0, true)

  // Polarity is positive on at least one spectrum.
  bool found_positive = false;
  for (const MSSpectrum& s : original.getSpectra())
  {
    if (s.getInstrumentSettings().getPolarity() == IonSource::Polarity::POSITIVE)
    {
      found_positive = true;
      break;
    }
  }
  TEST_EQUAL(found_positive, true)

  std::string tmp_mzml = SystemSettings::getTempDirectory() + "/" + File::getUniqueName() + "_thermo_roundtrip.mzML";
  MzMLFile().store(tmp_mzml, original);

  MSExperiment reloaded;
  MzMLFile().load(tmp_mzml, reloaded);
  File::remove(tmp_mzml);

  TEST_EQUAL(original.getSourceFiles()[0].getChecksum().size(), 40)
  TEST_EQUAL(reloaded.getSourceFiles()[0].getChecksum(), original.getSourceFiles()[0].getChecksum())
  TEST_EQUAL(reloaded.getInstrument().getMetaValue("instrument serial number"), original.getInstrument().getMetaValue("instrument serial number"))
  TEST_EQUAL(reloaded.getMetaValue("Thermo instrument methods"), original.getMetaValue("Thermo instrument methods"))
  TEST_TRUE(original.getSample().metaValueExists("Thermo injection volume"))
  TEST_EQUAL(reloaded.getSample().getMetaValue("Thermo injection volume"), original.getSample().getMetaValue("Thermo injection volume"))
  bool saw_supplemental = false;
  for (Size i = 0; i < original.size(); ++i)
  {
    TEST_EQUAL(original[i].getNativeID(), reloaded[i].getNativeID())
    TEST_EQUAL(original[i].getPrecursors().size(), reloaded[i].getPrecursors().size())
    TEST_EQUAL(original[i].getAcquisitionInfo()[0].getMetaValue("Thermo trailer extra"), reloaded[i].getAcquisitionInfo()[0].getMetaValue("Thermo trailer extra"))
    for (Size j = 0; j < original[i].getPrecursors().size(); ++j)
    {
      const auto& before = original[i].getPrecursors()[j];
      const auto& after = reloaded[i].getPrecursors()[j];
      TEST_REAL_SIMILAR(before.getMZ(), after.getMZ())
      TEST_REAL_SIMILAR(before.getMetaValue("selected ion m/z"), after.getMetaValue("selected ion m/z", after.getMZ()))
      TEST_REAL_SIMILAR(before.getIsolationWindowLowerOffset(), after.getIsolationWindowLowerOffset())
      TEST_REAL_SIMILAR(before.getIsolationWindowUpperOffset(), after.getIsolationWindowUpperOffset())
      TEST_EQUAL(before.getMetaValue("spectrum_ref"), after.getMetaValue("spectrum_ref"))
      if (before.metaValueExists("supplemental collision energy"))
      {
        saw_supplemental = true;
        TEST_EQUAL(before.getMetaValue("supplemental collision energy"), after.getMetaValue("supplemental collision energy"))
      }
    }
  }
  TEST_TRUE(saw_supplemental)

  TEST_EQUAL(original.size(), reloaded.size())
  TEST_EQUAL(original.getSourceFiles().size(), reloaded.getSourceFiles().size())

  if (!original.empty())
  {
    TEST_EQUAL(original[0].getMSLevel(), reloaded[0].getMSLevel())
    TEST_REAL_SIMILAR(original[0].getRT(), reloaded[0].getRT())
    TEST_EQUAL(original[0].size(), reloaded[0].size())
  }

  // Per-spectrum MS-level round-trip and precursor / activation-method survey.
  // AllScans has MSn so the precursor assertion must fire.
  bool found_msn_precursor = false;
  std::set<Precursor::ActivationMethod> activation_methods_seen;
  for (Size i = 0; i < reloaded.size(); ++i)
  {
    TEST_EQUAL(original[i].getMSLevel(), reloaded[i].getMSLevel())
    if (reloaded[i].getMSLevel() > 1 && !reloaded[i].getPrecursors().empty())
    {
      const Precursor& prec = reloaded[i].getPrecursors()[0];
      if (prec.getMZ() > 0.0) { found_msn_precursor = true; }
      for (auto am : prec.getActivationMethods()) { activation_methods_seen.insert(am); }
    }
  }
  TEST_EQUAL(found_msn_precursor, true)
  // Bridge reports at least one activation method for the dependent MS2 scans.
  TEST_EQUAL(activation_methods_seen.empty(), false)

  // --- Extended instrument metadata ---
  // All Angiotensin scans are FTMS on an Orbitrap, acquired by ESI.
  TEST_EQUAL(original.getInstrument().getMassAnalyzers().empty(), false)
  TEST_EQUAL(original.getInstrument().getIonSources().empty(), false)
  TEST_EQUAL(original.getInstrument().getIonDetectors().empty(), false)
  TEST_EQUAL(original.getInstrument().getSoftware().getVersion().empty(), false)
  bool has_orbitrap = false;
  for (const MassAnalyzer& ma : original.getInstrument().getMassAnalyzers())
  {
    if (ma.getType() == MassAnalyzer::AnalyzerType::ORBITRAP) { has_orbitrap = true; }
  }
  TEST_EQUAL(has_orbitrap, true)

  // --- Extended per-spectrum metadata ---
  Size with_filter = 0, with_tic = 0, with_basepeak = 0, with_window = 0, with_mzrange = 0, with_injection = 0,
       with_resolving_power = 0;
  for (const MSSpectrum& s : original.getSpectra())
  {
    if (s.metaValueExists("filter string")) { ++with_filter; }
    // the resolution trailer of the FTMS scans, as msconvert reports it (MS:1000800)
    if (s.metaValueExists("mass resolving power")) { ++with_resolving_power; }
    if (s.metaValueExists("total ion current")) { ++with_tic; }
    if (s.metaValueExists("base peak m/z")) { ++with_basepeak; }
    if (s.metaValueExists("lowest observed m/z") && s.metaValueExists("highest observed m/z")) { ++with_mzrange; }
    if (!s.getInstrumentSettings().getScanWindows().empty()) { ++with_window; }
    for (Size ai = 0; ai < s.getAcquisitionInfo().size(); ++ai)
    {
      if (s.getAcquisitionInfo()[ai].metaValueExists("ion injection time")) { ++with_injection; break; }
    }
  }
  TEST_EQUAL(with_filter > 0, true)
  TEST_EQUAL(with_tic > 0, true)
  TEST_EQUAL(with_basepeak > 0, true)
  TEST_EQUAL(with_mzrange > 0, true)
  TEST_EQUAL(with_window > 0, true)
  TEST_EQUAL(with_injection > 0, true)
  TEST_EQUAL(with_resolving_power > 0, true)

  // --- TIC chromatogram extracted and surviving the mzML round-trip ---
  TEST_EQUAL(original.getChromatograms().empty(), false)
  TEST_EQUAL(reloaded.getChromatograms().empty(), false)
  bool reloaded_has_tic = false;
  for (const MSChromatogram& c : reloaded.getChromatograms())
  {
    if (c.getChromatogramType() == ChromatogramSettings::ChromatogramType::TOTAL_ION_CURRENT_CHROMATOGRAM)
    {
      reloaded_has_tic = true;
    }
  }
  TEST_EQUAL(reloaded_has_tic, true)

  // --- Spectrum metadata survives the mzML round-trip ---
  Size rt_filter = 0, rt_tic = 0, rt_window = 0, rt_resolving_power = 0;
  for (const MSSpectrum& s : reloaded.getSpectra())
  {
    if (s.metaValueExists("filter string")) { ++rt_filter; }
    if (s.metaValueExists("total ion current")) { ++rt_tic; }
    if (!s.getInstrumentSettings().getScanWindows().empty()) { ++rt_window; }
    if (s.metaValueExists("mass resolving power")) { ++rt_resolving_power; }
  }
  TEST_EQUAL(rt_filter > 0, true)
  TEST_EQUAL(rt_tic > 0, true)
  TEST_EQUAL(rt_window > 0, true)
  // written as MS:1000800 in the scan and read back as the same value
  TEST_EQUAL(rt_resolving_power, with_resolving_power)
  for (Size i = 0; i < original.size() && i < reloaded.size(); ++i)
  {
    if (!original[i].metaValueExists("mass resolving power")) { continue; }
    TEST_EQUAL(reloaded[i].getMetaValue("mass resolving power").toString(), original[i].getMetaValue("mass resolving power").toString())
  }
}
END_SECTION

START_SECTION(FileHandler::loadExperiment applies the vendor peak picking)
{
  // FileHandler reads .raw as FileConverter does: centroided by Thermo's peak picking.
  // ThermoRawFile keeps the scans as acquired; the MS1 scans of this file are profile.
  TEST_EQUAL(ThermoRawFile().getOptions().centroid, false)
  MSExperiment acquired, centroided;
  ThermoRawFile().load(THERMO_RAW_TEST_DATA, acquired);
  FileHandler().loadExperiment(THERMO_RAW_TEST_DATA, centroided, {FileTypes::RAW});
  TEST_EQUAL(centroided.size(), acquired.size())

  Size same_id = 0, centroid_type = 0, ms1_profile = 0, profile_peaks = 0, picked_peaks = 0;
  for (Size i = 0; i < acquired.size() && i < centroided.size(); ++i)
  {
    if (centroided[i].getNativeID() == acquired[i].getNativeID()) { ++same_id; }
    if (centroided[i].getType() == SpectrumSettings::SpectrumType::CENTROID) { ++centroid_type; }
    if (acquired[i].getMSLevel() == 1 && acquired[i].getType() == SpectrumSettings::SpectrumType::PROFILE)
    {
      ++ms1_profile;
      profile_peaks += acquired[i].size();
      picked_peaks += centroided[i].size();
    }
  }
  TEST_EQUAL(same_id, acquired.size())
  TEST_EQUAL(centroid_type, centroided.size())
  TEST_EQUAL(ms1_profile > 0, true)
  TEST_EQUAL(picked_peaks < profile_peaks, true)
}
END_SECTION

START_SECTION([EXTRA] sampled noise survives RAW to mzML independently of spectrum peaks)
{
  const vector<string> noise_names = {"sampled noise m/z array", "sampled noise intensity array", "sampled noise baseline array"};
  for (bool centroid : {false, true})
  {
    ThermoRawFile file;
    auto options = file.getOptions();
    options.centroid = centroid;
    options.checksum = false;
    options.last_scan = 100;
    options.noise_data = true;
    file.setOptions(options);
    MSExperiment original;
    file.load(THERMO_RAW_TEST_DATA, original);

    Size with_noise = 0;
    for (const auto& spectrum : original)
    {
      const bool present = spectrum.metaValueExists(noise_names[0]);
      for (const auto& name : noise_names)
      {
        TEST_EQUAL(spectrum.metaValueExists(name), present)
        if (present)
        {
          TEST_EQUAL(spectrum.getMetaValue(name).toDoubleList().size(),
                     spectrum.getMetaValue(noise_names[0]).toDoubleList().size())
        }
      }
      if (present)
      {
        ++with_noise;
        TEST_FALSE(spectrum.getMetaValue(noise_names[0]).toDoubleList().empty())
      }
    }
    TEST_TRUE(with_noise > 0)

    MzMLFile mzml;
    string encoded;
    mzml.storeBuffer(encoded, original);
    for (const auto& accession : {"MS:1002743", "MS:1002744", "MS:1002745"})
    {
      TEST_TRUE(encoded.find(accession) != string::npos)
    }
    MSExperiment reloaded;
    mzml.loadBuffer(encoded, reloaded);
    TEST_EQUAL(reloaded.size(), original.size())
    for (Size i = 0; i < original.size() && i < reloaded.size(); ++i)
    {
      for (const auto& name : noise_names)
      {
        TEST_EQUAL(reloaded[i].metaValueExists(name), original[i].metaValueExists(name))
        if (original[i].metaValueExists(name))
        {
          TEST_EQUAL(reloaded[i].getMetaValue(name).toDoubleList(), original[i].getMetaValue(name).toDoubleList())
        }
      }
    }

    options.noise_data = false;
    file.setOptions(options);
    MSExperiment without_noise;
    file.load(THERMO_RAW_TEST_DATA, without_noise);
    TEST_EQUAL(without_noise.size(), original.size())
    for (Size i = 0; i < original.size() && i < without_noise.size(); ++i)
    {
      TEST_EQUAL(without_noise[i].size(), original[i].size())
      for (Size j = 0; j < original[i].size() && j < without_noise[i].size(); ++j)
      {
        TEST_EQUAL(without_noise[i][j], original[i][j])
      }
      for (const auto& name : noise_names)
      {
        TEST_FALSE(without_noise[i].metaValueExists(name))
      }
    }
  }
}
END_SECTION

START_SECTION([EXTRA] ThermoRawFile::Options scan selection - metadata_only - skip_chromatograms)
{
  // Every selection must give the matching subset of a full load, with identical spectra:
  // precursors (and their intensity, estimated from the parent scan) included.
  auto scanNumber = [](const MSSpectrum& s) { return std::stoi(s.getNativeID().substr(s.getNativeID().rfind('=') + 1)); };
  ThermoRawFile::Options base;
  base.centroid = true;
  base.checksum = false;
  auto load = [&](const ThermoRawFile::Options& options) {
    ThermoRawFile file;
    file.setOptions(options);
    MSExperiment exp;
    file.load(THERMO_RAW_TEST_DATA, exp);
    return exp;
  };
  auto testSubset = [](const MSExperiment& selected, const MSExperiment& full, const std::function<bool(const MSSpectrum&)>& keep) {
    Size expected = 0, same = 0;
    for (const MSSpectrum& s : full.getSpectra())
    {
      if (! keep(s)) { continue; }
      if (expected < selected.size() && selected[expected] == s) { ++same; }
      ++expected;
    }
    TEST_EQUAL(selected.size(), expected)
    TEST_EQUAL(same, expected)
  };

  const MSExperiment full = load(base);
  TEST_EQUAL(full.size() > 200, true)
  TEST_EQUAL(full.getChromatograms().size(), 1)

  ThermoRawFile::Options ms2 = base;
  ms2.ms_levels = {2};
  const MSExperiment only_ms2 = load(ms2);
  TEST_EQUAL(only_ms2.empty(), false)
  testSubset(only_ms2, full, [](const MSSpectrum& s) { return s.getMSLevel() == 2; });
  // the TIC still covers all scans
  TEST_EQUAL(only_ms2.getChromatograms().size(), 1)
  TEST_EQUAL(only_ms2.getChromatograms()[0] == full.getChromatograms()[0], true)

  const int first = scanNumber(full[50]), last = scanNumber(full[150]);
  ThermoRawFile::Options scans = base;
  scans.first_scan = first;
  scans.last_scan = last;
  testSubset(load(scans), full, [&](const MSSpectrum& s) { return scanNumber(s) >= first && scanNumber(s) <= last; });

  const double rt_min = full[100].getRT(), rt_max = full[200].getRT();
  ThermoRawFile::Options rt = base;
  rt.rt_min = rt_min;
  rt.rt_max = rt_max;
  const MSExperiment rt_range = load(rt);
  testSubset(rt_range, full, [&](const MSSpectrum& s) { return s.getRT() >= rt_min && s.getRT() < rt_max; });
  // the TIC is read in full from the vendor chromatogram when the load stops early
  TEST_EQUAL(rt_range.getChromatograms().size(), 1)
  TEST_EQUAL(rt_range.getChromatograms()[0].size(), full.getChromatograms()[0].size())

  // MS1 spectra and spectra without precursor are kept, as by PeakFileOptions::setPrecursorMZRange
  for (const auto& [lo, hi] : {std::pair{400.0, 500.0}, std::pair{500.0, 600.0}, std::pair{659.0, 660.0}})
  {
    ThermoRawFile::Options precursor = base;
    precursor.precursor_mz_min = lo;
    precursor.precursor_mz_max = hi;
    const MSExperiment selected = load(precursor);
    testSubset(selected, full, [lo, hi](const MSSpectrum& s) {
      return s.getMSLevel() <= 1 || s.getPrecursors().empty() || (s.getPrecursors()[0].getMZ() >= lo && s.getPrecursors()[0].getMZ() < hi);
    });
    TEST_EQUAL(selected.size() < full.size(), true)
  }

  ThermoRawFile::Options metadata = base;
  metadata.metadata_only = true;
  metadata.skip_chromatograms = true;
  const MSExperiment headers = load(metadata);
  TEST_EQUAL(headers.size(), full.size())
  TEST_EQUAL(headers.getChromatograms().size(), 0)
  Size empty = 0, same_id = 0, same_precursor = 0;
  for (Size i = 0; i < headers.size() && i < full.size(); ++i)
  {
    if (headers[i].empty()) { ++empty; }
    if (headers[i].getNativeID() == full[i].getNativeID() && headers[i].getRT() == full[i].getRT()
        && headers[i].getMSLevel() == full[i].getMSLevel()) { ++same_id; }
    if (headers[i].getPrecursors().size() == full[i].getPrecursors().size()
        && (full[i].getPrecursors().empty() || headers[i].getPrecursors()[0].getMZ() == full[i].getPrecursors()[0].getMZ())) { ++same_precursor; }
  }
  TEST_EQUAL(empty, headers.size())
  TEST_EQUAL(same_id, headers.size())
  TEST_EQUAL(same_precursor, headers.size())
}
END_SECTION

START_SECTION([EXTRA] read-ahead on worker threads gives the single-threaded result)
{
#ifdef _OPENMP
  const int threads = omp_get_max_threads();
  auto load = [](int n, const ThermoRawFile::Options& options) {
    omp_set_num_threads(n);
    ThermoRawFile file;
    file.setOptions(options);
    MSExperiment exp;
    file.load(THERMO_RAW_TEST_DATA, exp);
    return exp;
  };
  ThermoRawFile::Options options;
  options.checksum = false;
  const MSExperiment single = load(1, options);
  const MSExperiment ahead = load(4, options);
  TEST_TRUE(single.size() > 0)
  TEST_TRUE(single == ahead)
  // Stopping early (last_scan) and header-only scans (MS-level filter) on the read-ahead path.
  options.last_scan = 150;
  TEST_TRUE(load(1, options) == load(4, options))
  options.last_scan = 0;
  options.ms_levels = {2};
  TEST_TRUE(load(1, options) == load(4, options))
  omp_set_num_threads(threads);
#endif
}
END_SECTION

START_SECTION(FileHandler::loadExperiment passes PeakFileOptions to the reader)
{
  ThermoRawFile::Options base;
  base.centroid = true;
  base.checksum = false;
  ThermoRawFile file;
  file.setOptions(base);
  MSExperiment full;
  file.load(THERMO_RAW_TEST_DATA, full);

  const double rt_min = full[100].getRT(), rt_max = full[300].getRT();
  FileHandler handler;
  handler.getOptions().addMSLevel(2);
  handler.getOptions().setRTRange(DRange<1>(DPosition<1>(rt_min), DPosition<1>(rt_max)));
  MSExperiment loaded;
  handler.loadExperiment(THERMO_RAW_TEST_DATA, loaded, {FileTypes::RAW});
  Size expected = 0, same = 0;
  for (const MSSpectrum& s : full.getSpectra())
  {
    if (s.getMSLevel() != 2 || s.getRT() < rt_min || s.getRT() >= rt_max) { continue; }
    if (expected < loaded.size() && loaded[expected] == s) { ++same; }
    ++expected;
  }
  TEST_EQUAL(expected > 0, true)
  TEST_EQUAL(loaded.size(), expected)
  TEST_EQUAL(same, expected)
  // no checksum is computed by default
  TEST_EQUAL(loaded.getSourceFiles().size(), 1)
  TEST_EQUAL(loaded.getSourceFiles()[0].getChecksum(), "")

  // fill data off: spectra without peaks; skip chromatograms is honoured
  FileHandler headers;
  headers.getOptions().setFillData(false);
  headers.getOptions().setSkipChromatograms(true);
  MSExperiment no_peaks;
  headers.loadExperiment(THERMO_RAW_TEST_DATA, no_peaks, {FileTypes::RAW});
  TEST_EQUAL(no_peaks.size(), full.size())
  TEST_EQUAL(no_peaks.getChromatograms().size(), 0)
  Size empty = 0;
  for (const MSSpectrum& s : no_peaks.getSpectra())
  {
    if (s.empty()) { ++empty; }
  }
  TEST_EQUAL(empty, no_peaks.size())
}
END_SECTION

START_SECTION(void transform(const std::string& path, Interfaces::IMSDataConsumer* consumer))
{
  TEST_EXCEPTION(Exception::NullPointer, ThermoRawFile().transform(THERMO_RAW_TEST_DATA, nullptr))

  ThermoRawFile::Options base;
  base.centroid = true;
  base.checksum = false;
  auto load = [](const ThermoRawFile::Options& options) {
    ThermoRawFile file;
    file.setOptions(options);
    MSExperiment exp;
    file.load(THERMO_RAW_TEST_DATA, exp);
    return exp;
  };
  // transform() into a collecting consumer must give load(), announced in advance and in RT order
  auto testTransform = [&](const ThermoRawFile::Options& options) {
    RecordingConsumer collected;
    ThermoRawFile file;
    file.setOptions(options);
    file.transform(THERMO_RAW_TEST_DATA, &collected);
    const MSExperiment loaded = load(options);
    TEST_TRUE(collected.order_ok)
    TEST_EQUAL(collected.expected_spectra, collected.exp.size())
    TEST_EQUAL(collected.expected_chromatograms, collected.exp.getChromatograms().size())
    TEST_TRUE(collected.exp.isSorted(false))
    TEST_EQUAL(collected.exp.size(), loaded.size())
    TEST_TRUE(collected.exp == loaded)
  };

  testTransform(base);
  ThermoRawFile::Options acquired = base;
  acquired.centroid = false;
  testTransform(acquired);
  ThermoRawFile::Options ms2 = base;
  ms2.ms_levels = {2};
  testTransform(ms2);
  ThermoRawFile::Options ms1 = base;
  ms1.ms_levels = {1};
  testTransform(ms1);
  const MSExperiment full = load(base);
  ThermoRawFile::Options rt = base;
  rt.rt_min = full[100].getRT();
  rt.rt_max = full[200].getRT();
  testTransform(rt);
  ThermoRawFile::Options scans = base;
  scans.first_scan = 51;
  scans.last_scan = 151;
  testTransform(scans);
  ThermoRawFile::Options precursor = base;
  precursor.precursor_mz_min = 500.0;
  precursor.precursor_mz_max = 600.0;
  testTransform(precursor);
  ThermoRawFile::Options metadata = base;
  metadata.metadata_only = true;
  metadata.skip_chromatograms = true;
  testTransform(metadata);
  ThermoRawFile::Options extras = base;
  extras.charge_data = true;
  extras.noise_data = true;
  extras.checksum = true;
  testTransform(extras);
}
END_SECTION

START_SECTION(void FileHandler::transformExperiment(const std::string& filename, Interfaces::IMSDataConsumer* consumer, ...))
{
  // The same PeakFileOptions give the same spectra as loadExperiment, the peak-level filters included
  FileHandler handler;
  handler.getOptions().addMSLevel(2);
  handler.getOptions().setMZRange(DRange<1>(DPosition<1>(200.0), DPosition<1>(800.0)));
  handler.getOptions().setSkipChromatograms(true);
  MSExperiment loaded;
  handler.loadExperiment(THERMO_RAW_TEST_DATA, loaded, {FileTypes::RAW});
  RecordingConsumer collected;
  handler.transformExperiment(THERMO_RAW_TEST_DATA, &collected, {FileTypes::RAW});
  TEST_TRUE(collected.order_ok)
  TEST_TRUE(loaded.size() > 0)
  TEST_EQUAL(collected.exp.size(), loaded.size())
  TEST_EQUAL(collected.exp.getChromatograms().size(), 0)
  TEST_TRUE(collected.exp == loaded)

  // with the source file rewritten (absolute path, file name)
  MSExperiment loaded_rewritten;
  FileHandler().loadExperiment(THERMO_RAW_TEST_DATA, loaded_rewritten, {FileTypes::RAW}, ProgressLogger::NONE, true);
  RecordingConsumer rewritten;
  FileHandler().transformExperiment(THERMO_RAW_TEST_DATA, &rewritten, {FileTypes::RAW}, ProgressLogger::NONE, true);
  TEST_EQUAL(rewritten.exp.getSourceFiles().size(), 1)
  TEST_TRUE(rewritten.exp == loaded_rewritten)

  TEST_EXCEPTION(Exception::ParseError, FileHandler().transformExperiment(THERMO_RAW_TEST_DATA, &collected, {FileTypes::MZML}))
}
END_SECTION

START_SECTION([EXTRA] a ZIP archive that holds the RAW file is read like the RAW file)
{
  // Compressed (addOrReplaceFromFile) at the archive root, and stored (zipDirectory) below a directory next to the
  // AppleDouble file that macOS adds, with the upper-case names that Windows acquisition PCs often write
  TempDir tmp;
  const std::string deflated = tmp.getPath() + "/deflated.raw.zip";
  ZipArchiveFile::addOrReplaceFromFile(deflated, "Angiotensin_AllScans.raw", THERMO_RAW_TEST_DATA);
  const std::string tree = tmp.getPath() + "/tree";
  File::makeDir(tree + "/run");
  File::makeDir(tree + "/__MACOSX/run");
  File::copy(THERMO_RAW_TEST_DATA, tree + "/run/ANGIOTENSIN.RAW");
  {
    std::ofstream ofs((tree + "/__MACOSX/run/._ANGIOTENSIN.RAW").c_str(), std::ios::binary);
    ofs << "AppleDouble";
  }
  const std::string stored = tmp.getPath() + "/STORED.RAW.ZIP";
  ZipArchiveFile::zipDirectory(tree, stored);

  ThermoRawFile file; // with the checksum
  MSExperiment expected;
  file.load(THERMO_RAW_TEST_DATA, expected);
  MSExperiment handled_raw;
  FileHandler().loadExperiment(THERMO_RAW_TEST_DATA, handled_raw);
  for (const std::string& archive : {deflated, stored})
  {
    MSExperiment exp;
    file.load(archive, exp);
    TEST_EQUAL(exp.size(), expected.size())
    TEST_TRUE(exp.getSpectra() == expected.getSpectra())
    TEST_TRUE(exp.getChromatograms() == expected.getChromatograms())
    // The source file is the archive, as for a Bruker .d.zip; its checksum is the one of the RAW file
    TEST_EQUAL(exp.getSourceFiles().size(), 1)
    TEST_STRING_EQUAL(exp.getSourceFiles()[0].getNameOfFile(), File::basename(archive))
    TEST_STRING_EQUAL(exp.getSourceFiles()[0].getChecksum(), expected.getSourceFiles()[0].getChecksum())

    RecordingConsumer collected;
    file.transform(archive, &collected);
    TEST_TRUE(collected.order_ok)
    TEST_TRUE(collected.exp == exp)

    TEST_EQUAL(FileHandler::getType(archive), FileTypes::RAW)
    MSExperiment handled;
    FileHandler().loadExperiment(archive, handled);
    TEST_EQUAL(handled.size(), handled_raw.size())
    TEST_TRUE(handled.getSpectra() == handled_raw.getSpectra())
  }
}
END_SECTION
#endif

#ifdef THERMO_FAIMS_DIA_TEST_DATA
START_SECTION(real Thermo FAIMS-DIA RAW -> FAIMS-aware SWATH maps)
{
  ThermoRawFile file;
  MSExperiment exp;
  file.load(THERMO_FAIMS_DIA_TEST_DATA, exp);

  TEST_EQUAL(exp.empty(), false)
  TEST_EQUAL(exp.getSourceFiles().empty(), false)

  Size ms1_count = 0;
  Size ms2_count = 0;
  Size explicit_faims_count = 0;
  Size dia_ms2_count = 0;
  std::set<double> faims_cvs;
  std::map<std::pair<UInt, std::string>, Size> resolving_powers; // (MS level, value) -> spectra

  for (const MSSpectrum& spectrum : exp.getSpectra())
  {
    ++resolving_powers[{spectrum.getMSLevel(), spectrum.getMetaValue("mass resolving power", "none").toString()}];
    if (spectrum.getMSLevel() == 1)
    {
      ++ms1_count;
    }
    else if (spectrum.getMSLevel() == 2)
    {
      ++ms2_count;

      if (!spectrum.getPrecursors().empty())
      {
        const Precursor& precursor = spectrum.getPrecursors()[0];
        if (precursor.getMZ() > 0.0 &&
            precursor.getIsolationWindowLowerOffset() > 0.0 &&
            precursor.getIsolationWindowUpperOffset() > 0.0)
        {
          ++dia_ms2_count;
        }
      }
    }

    if (spectrum.getDriftTimeUnit() == DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE)
    {
      ++explicit_faims_count;
      faims_cvs.insert(spectrum.getDriftTime());
    }
  }

  // Pin the native-reader side of the integration boundary to this public fixture.
  TEST_EQUAL(exp.size(), 92778)
  TEST_EQUAL(ms1_count, 1974)
  TEST_EQUAL(ms2_count, 90804)
  TEST_EQUAL(dia_ms2_count, ms2_count)
  TEST_EQUAL(explicit_faims_count, exp.size())
  TEST_EQUAL(faims_cvs.size(), 1)
  if (!faims_cvs.empty())
  {
    TEST_REAL_SIMILAR(*faims_cvs.begin(), -5.0)
  }
  // resolving power from the resolution entry of the scan trailer: 120000 for MS1, 30000 for MS2 (as msconvert reports it)
  TEST_EQUAL(resolving_powers.size(), 2)
  TEST_EQUAL((resolving_powers[{1, "120000"}]), 1974)
  TEST_EQUAL((resolving_powers[{2, "30000"}]), 90804)

  cout << "Thermo FAIMS-DIA validation: spectra=" << exp.size()
       << ", MS1=" << ms1_count
       << ", MS2=" << ms2_count
       << ", DIA-like MS2=" << dia_ms2_count
       << ", explicit FAIMS spectra=" << explicit_faims_count
       << ", unique CVs=" << faims_cvs.size() << '\n';
  cout << "FAIMS CVs:";
  for (double cv : faims_cvs) { cout << ' ' << cv; }
  cout << '\n';

  std::shared_ptr<ExperimentalSettings> exp_meta;
  auto groups = SwathFile().loadFromMSExperimentByFAIMSCV(
    std::move(exp), SystemSettings::getTempDirectory() + "/", exp_meta, "normal");

  TEST_EQUAL(groups.size(), 1)
  TEST_EQUAL(exp_meta != nullptr, true)
  if (!groups.empty())
  {
    TEST_REAL_SIMILAR(groups[0].faims_cv, -5.0)
  }

  Size groups_with_ms1 = 0;
  Size total_swath_windows = 0;
  for (const auto& group : groups)
  {
    TEST_EQUAL(std::isfinite(group.faims_cv), true)
    TEST_EQUAL(faims_cvs.count(group.faims_cv), 1)
    TEST_EQUAL(group.swath_maps.empty(), false)

    Size group_ms1_maps = 0;
    Size group_swath_windows = 0;
    Size group_ms2_spectra = 0;

    for (const auto& map : group.swath_maps)
    {
      TEST_EQUAL(map.sptr != nullptr, true)
      if (map.sptr == nullptr) { continue; }

      if (map.ms1)
      {
        ++group_ms1_maps;
      }
      else
      {
        ++group_swath_windows;
        group_ms2_spectra += map.sptr->getNrSpectra();
        TEST_EQUAL(map.lower < map.upper, true)
        TEST_EQUAL(map.sptr->getNrSpectra() > 0, true)
      }
    }

    // SwathFile should emit at most one MS1 map and at least one DIA window per CV.
    TEST_EQUAL(group_ms1_maps <= 1, true)
    TEST_EQUAL(group_swath_windows > 0, true)
    TEST_EQUAL(group_ms2_spectra > 0, true)

    if (group_ms1_maps == 1) { ++groups_with_ms1; }
    total_swath_windows += group_swath_windows;

    cout << "  CV " << group.faims_cv
         << ": maps=" << group.swath_maps.size()
         << ", MS1 maps=" << group_ms1_maps
         << ", DIA windows=" << group_swath_windows
         << ", MS2 spectra=" << group_ms2_spectra << '\n';
  }

  TEST_EQUAL(groups_with_ms1, 1)
  TEST_EQUAL(total_swath_windows, 75)

  Size total_ms2_spectra = 0;
  if (!groups.empty())
  {
    TEST_EQUAL(groups[0].swath_maps.size(), 76)
    for (const auto& map : groups[0].swath_maps)
    {
      if (!map.ms1 && map.sptr != nullptr)
      {
        total_ms2_spectra += map.sptr->getNrSpectra();
      }
    }
  }
  TEST_EQUAL(total_ms2_spectra, ms2_count)

  cout << "FAIMS SWATH groups=" << groups.size()
       << ", groups with MS1=" << groups_with_ms1
       << ", total DIA windows=" << total_swath_windows << '\n';
}
END_SECTION

START_SECTION([EXTRA] SwathFile::loadThermoRaw with readoptions onDisc equals normal)
{
  // The onDisc maps read the spectra from the file; they must expose what the in-memory maps hold.
  std::shared_ptr<ExperimentalSettings> meta_normal, meta_on_disc;
  const std::string tmp = SystemSettings::getTempDirectory() + "/";
  std::vector<OpenSwath::SwathMap> normal = SwathFile().loadThermoRaw(THERMO_FAIMS_DIA_TEST_DATA, tmp, meta_normal, "normal");
  std::vector<OpenSwath::SwathMap> on_disc = SwathFile().loadThermoRaw(THERMO_FAIMS_DIA_TEST_DATA, tmp, meta_on_disc, "onDisc");
  TEST_EQUAL(on_disc.size(), normal.size())
  TEST_TRUE(*meta_on_disc == *meta_normal)
  Size mismatches = 0, compared = 0;
  for (Size m = 0; m < std::min(normal.size(), on_disc.size()); ++m)
  {
    const OpenSwath::SwathMap& a = normal[m];
    const OpenSwath::SwathMap& b = on_disc[m];
    TEST_EQUAL(b.ms1, a.ms1)
    TEST_EQUAL(b.lower, a.lower)
    TEST_EQUAL(b.upper, a.upper)
    TEST_EQUAL(b.center, a.center)
    TEST_EQUAL(b.imLower, a.imLower)
    TEST_EQUAL(b.imUpper, a.imUpper)
    TEST_EQUAL(b.sptr->getNrSpectra(), a.sptr->getNrSpectra())
    if (b.sptr->getNrSpectra() != a.sptr->getNrSpectra()) { continue; }
    const int n = static_cast<int>(a.sptr->getNrSpectra());
    for (int i : {0, n / 2, n - 1})
    {
      if (i < 0) { continue; }
      ++compared;
      const OpenSwath::SpectrumMeta ma = a.sptr->getSpectrumMetaById(i), mb = b.sptr->getSpectrumMetaById(i);
      const OpenSwath::SpectrumPtr sa = a.sptr->getSpectrumById(i), sb = b.sptr->getSpectrumById(i);
      if (ma.RT != mb.RT || ma.ms_level != mb.ms_level ||
          sa->getMZArray()->data != sb->getMZArray()->data ||
          sa->getIntensityArray()->data != sb->getIntensityArray()->data)
      {
        ++mismatches;
      }
    }
    const double rt = a.sptr->getSpectrumMetaById(n / 2).RT;
    TEST_TRUE(b.sptr->getSpectraByRT(rt, 10.0) == a.sptr->getSpectraByRT(rt, 10.0))
  }
  TEST_EQUAL(compared > 0, true)
  TEST_EQUAL(mismatches, 0)
}
END_SECTION
#endif

END_TEST

#else  // !WITH_THERMO_RAW

// Placeholder main so the test target still links when the bridge is disabled
// (e.g. Linux aarch64 default). Returns 0 (treated as a passing no-op test).
int main(int /*argc*/, const char** /*argv*/)
{
  return 0;
}

#endif // WITH_THERMO_RAW
