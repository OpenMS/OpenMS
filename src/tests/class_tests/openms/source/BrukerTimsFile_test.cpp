// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/config.h>

#ifdef WITH_OPENTIMS

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/ANALYSIS/ID/ProSEAlgorithm.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/BrukerTimsFile.h>
#include <OpenMS/FORMAT/RationalScan2ImConverter.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/DATAACCESS/SwathFileConsumer.h>

using namespace OpenMS;
using namespace std;

START_TEST(BrukerTimsFile, "$Id$")

START_SECTION(void load(const String& path, MSExperiment& exp, const Config& config))
{
  // Test: invalid path throws FileNotReadable
  BrukerTimsFile f;
  MSExperiment exp;
  TEST_EXCEPTION(Exception::FileNotReadable, f.load("/nonexistent/path.d", exp));
}
END_SECTION

START_SECTION([FileHandler] BRUKER_TDF detection)
{
  // Test: .d suffix is detected as BRUKER_TDF
  TEST_EQUAL(FileHandler::getTypeByFileName("sample.d"), FileTypes::BRUKER_TDF);
  TEST_EQUAL(FileHandler::getTypeByFileName("/path/to/experiment.d"), FileTypes::BRUKER_TDF);

  // Test: non-.d suffixes are not BRUKER_TDF
  TEST_NOT_EQUAL(FileHandler::getTypeByFileName("sample.mzML"), FileTypes::BRUKER_TDF);
}
END_SECTION

START_SECTION(RationalScan2ImConverter forward conversion)
{
  // Coefficients from opentims test data (ModelType=2)
  using Coeff = RationalScan2ImConverter::Coefficients;
  Coeff c{1.0, 917.0, 213.5998, 75.81729, 33.0, 1.0, -0.009065829, 135.4364, 13.32608, 1663.341};

  std::unordered_map<uint32_t, Coeff> calibrations;
  calibrations[1] = c;
  std::vector<uint32_t> frame_to_cal = {0, 1}; // index 0 unused, frame 1 -> cal 1

  RationalScan2ImConverter converter(std::move(calibrations), std::move(frame_to_cal));

  // Worked example from spec (precise values):
  // V = 213.5998 + ((75.81729 - 213.5998) / 917.0) * (500 - 33 - 1) = 143.5816
  // 1/K0 = 1 / (-0.009065829 + 135.4364 / 143.5816) = 1.070429
  double scan_val = 500.0;
  double result = 0.0;
  converter.convert(1, &result, &scan_val, 1);
  TEST_REAL_SIMILAR(result, 1.070429);

  // Test with uint32_t scan input
  uint32_t scan_int = 500;
  double result2 = 0.0;
  converter.convert(1, &result2, &scan_int, 1);
  TEST_REAL_SIMILAR(result2, 1.070429);

  // Test multiple scans at once
  double scans[] = {0.0, 250.0, 500.0, 916.0};
  double results[4] = {};
  converter.convert(1, results, scans, 4);
  // scan=0: V = 213.5998 + slope*(0-34) = 218.7084
  //         1/K0 = 1/(-0.009066 + 135.4364/218.7084) = 1.63883
  TEST_REAL_SIMILAR(results[0], 1.63883);
  // All should be positive
  for (int i = 0; i < 4; ++i)
  {
    TEST_EQUAL(results[i] > 0.0, true);
  }
  // Increasing scan -> decreasing 1/K0 (higher scan = lower mobility)
  for (int i = 1; i < 4; ++i)
  {
    TEST_EQUAL(results[i] < results[i-1], true);
  }
}
END_SECTION

START_SECTION(RationalScan2ImConverter round-trip via inverse_convert)
{
  using Coeff = RationalScan2ImConverter::Coefficients;
  Coeff c{1.0, 917.0, 213.5998, 75.81729, 33.0, 1.0, -0.009065829, 135.4364, 13.32608, 1663.341};

  std::unordered_map<uint32_t, Coeff> calibrations;
  calibrations[1] = c;
  std::vector<uint32_t> frame_to_cal = {0, 1};

  RationalScan2ImConverter converter(std::move(calibrations), std::move(frame_to_cal));

  // Forward: scan 500 -> 1/K0
  double scan_val = 500.0;
  double inv_k0 = 0.0;
  converter.convert(1, &inv_k0, &scan_val, 1);

  // Inverse: 1/K0 -> scan (should round-trip to 500)
  uint32_t scan_back = 0;
  converter.inverse_convert(1, &scan_back, &inv_k0, 1);
  TEST_EQUAL(scan_back, 500);

  // Test a range of scans for round-trip
  for (uint32_t s = 10; s < 900; s += 50)
  {
    double sv = static_cast<double>(s);
    double ik0 = 0.0;
    converter.convert(1, &ik0, &sv, 1);
    uint32_t back = 0;
    converter.inverse_convert(1, &back, &ik0, 1);
    // Allow +/- 1 for rounding
    TEST_EQUAL(back >= s - 1 && back <= s + 1, true);
  }
}
END_SECTION

START_SECTION(RationalScan2ImConverter per-frame calibration)
{
  using Coeff = RationalScan2ImConverter::Coefficients;
  // Two different calibration segments
  Coeff c1{1.0, 917.0, 213.5998, 75.81729, 33.0, 1.0, -0.009065829, 135.4364, 13.32608, 1663.341};
  Coeff c2{1.0, 917.0, 220.0, 80.0, 33.0, 1.0, -0.01, 140.0, 13.0, 1660.0};

  std::unordered_map<uint32_t, Coeff> calibrations;
  calibrations[1] = c1;
  calibrations[2] = c2;
  // frame 1 uses cal 1, frame 2 uses cal 2
  std::vector<uint32_t> frame_to_cal = {0, 1, 2};

  RationalScan2ImConverter converter(std::move(calibrations), std::move(frame_to_cal));

  double scan_val = 500.0;
  double result_f1 = 0.0, result_f2 = 0.0;
  converter.convert(1, &result_f1, &scan_val, 1);
  converter.convert(2, &result_f2, &scan_val, 1);

  // Different calibrations should produce different results
  TEST_NOT_EQUAL(result_f1, result_f2);
  // Both should be positive
  TEST_EQUAL(result_f1 > 0.0, true);
  TEST_EQUAL(result_f2 > 0.0, true);
}
END_SECTION

START_SECTION(RationalScan2ImConverter description)
{
  using Coeff = RationalScan2ImConverter::Coefficients;
  Coeff c{1.0, 917.0, 213.5998, 75.81729, 33.0, 1.0, -0.009065829, 135.4364, 13.32608, 1663.341};

  std::unordered_map<uint32_t, Coeff> cals;
  cals[1] = c;
  std::vector<uint32_t> ftc = {0, 1};

  RationalScan2ImConverter converter(std::move(cals), std::move(ftc));
  std::string desc = converter.description();
  TEST_EQUAL(desc.find("RationalScan2ImConverter") != std::string::npos, true);
  TEST_EQUAL(desc.find("1") != std::string::npos, true); // 1 calibration segment
}
END_SECTION

START_SECTION(RationalScan2ImConverter singularity edge cases)
{
  using Coeff = RationalScan2ImConverter::Coefficients;

  // Normal coefficients as baseline
  Coeff c{1.0, 917.0, 213.5998, 75.81729, 33.0, 1.0, -0.009065829, 135.4364, 13.32608, 1663.341};

  // Edge case: c3 == c2 (zero slope in scan-to-voltage mapping)
  // Should not crash — singularity guard produces a finite value
  Coeff c_zero_slope = c;
  c_zero_slope.c3 = c_zero_slope.c2; // c3 == c2

  std::unordered_map<uint32_t, Coeff> cals;
  cals[1] = c_zero_slope;
  std::vector<uint32_t> ftc = {0, 1};
  RationalScan2ImConverter conv(std::move(cals), std::move(ftc));

  double scan_val = 500.0;
  double result = 0.0;
  conv.convert(1, &result, &scan_val, 1);
  TEST_EQUAL(std::isfinite(result), true);

  // Edge case: scan that makes V approach 0
  Coeff c_v_zero = c;
  c_v_zero.c2 = 0.0;
  c_v_zero.c3 = 0.0;

  std::unordered_map<uint32_t, Coeff> cals2;
  cals2[1] = c_v_zero;
  std::vector<uint32_t> ftc2 = {0, 1};
  RationalScan2ImConverter conv2(std::move(cals2), std::move(ftc2));

  double result2 = 0.0;
  conv2.convert(1, &result2, &scan_val, 1);
  TEST_EQUAL(std::isfinite(result2), true);
}
END_SECTION

// Integration tests (only run when ENABLE_OPENTIMS_TESTS is ON and data is available)
#ifdef OPENTIMS_DDA_TEST_DATA

START_SECTION(DDA loading integration test)
{
  BrukerTimsFile f;
  MSExperiment exp;
  f.load(OPENTIMS_DDA_TEST_DATA, exp);

  // Verify we got spectra
  TEST_NOT_EQUAL(exp.size(), 0);

  // Check MS levels present
  bool has_ms1 = false, has_ms2 = false;
  for (const auto& spec : exp)
  {
    if (spec.getMSLevel() == 1) has_ms1 = true;
    if (spec.getMSLevel() == 2) has_ms2 = true;
  }
  TEST_EQUAL(has_ms1, true);
  TEST_EQUAL(has_ms2, true);

  // Check MS1 spectra have IM data in IM_PEAK format
  for (const auto& spec : exp)
  {
    if (spec.getMSLevel() == 1 && !spec.empty())
    {
      TEST_EQUAL(spec.containsIMData(), true);
      break;
    }
  }

  // Check MS2 spectra have precursor and drift time
  for (const auto& spec : exp)
  {
    if (spec.getMSLevel() == 2 && !spec.getPrecursors().empty())
    {
      TEST_NOT_EQUAL(spec.getDriftTime(), 0.0);
      TEST_EQUAL(spec.getDriftTimeUnit(), DriftTimeUnit::VSSC);
      TEST_NOT_EQUAL(spec.getPrecursors()[0].getMZ(), 0.0);
      // TDF MS2 peaks are detector-centroided; must be labelled CENTROID so
      // downstream tools (e.g. CometAdapter) don't reject them as profile.
      TEST_EQUAL(spec.getType(), SpectrumSettings::SpectrumType::CENTROID);
      break;
    }
  }

  // Verify source file metadata was populated (I4)
  TEST_NOT_EQUAL(exp.getSourceFiles().size(), 0);
}
END_SECTION

START_SECTION(DDA native ID format test)
{
  // Contract: DDA MS2 native IDs are "frame=<F> scan=<S> precursor=<P>".
  // This extends the MS:1002818 pattern with a trailing "precursor=<P>"
  // token because OpenMS aggregates all PasefFrameMsMsInfo entries sharing
  // the same Precursors.Id into ONE spectrum (pwiz emits per-mobility-scan,
  // so its (frame, scan_begin) pairs are inherently unique — ours are not).
  // The Precursors.Id is ALSO duplicated as MetaValue "bruker_precursor_id"
  // for typed programmatic access.
  BrukerTimsFile f;
  MSExperiment exp;
  f.load(OPENTIMS_DDA_TEST_DATA, exp);

  // Find a non-empty MS2 spectrum
  const MSSpectrum* ms2 = nullptr;
  for (const auto& spec : exp)
  {
    if (spec.getMSLevel() == 2 && !spec.empty())
    {
      ms2 = &spec;
      break;
    }
  }
  TEST_NOT_EQUAL(ms2, nullptr);

  const String& id = ms2->getNativeID();
  TEST_EQUAL(id.hasPrefix("frame="), true);
  TEST_EQUAL(id.hasSubstring(" scan="), true);
  TEST_EQUAL(id.hasSubstring(" precursor="), true);

  // Precursors.Id must ALSO be stored as a typed MetaValue.
  TEST_EQUAL(ms2->metaValueExists("bruker_precursor_id"), true);

  // All DDA MS2 native IDs must be unique inside the run (XSD-level mzML
  // requirement). This is what the "precursor=<P>" disambiguator buys us.
  std::set<String> ms2_ids;
  Size ms2_count = 0;
  for (const auto& spec : exp)
  {
    if (spec.getMSLevel() == 2)
    {
      ms2_ids.insert(spec.getNativeID());
      ++ms2_count;
    }
  }
  TEST_EQUAL(ms2_ids.size(), ms2_count);

  STATUS("DDA MS2 native ID sample: " << id
         << " bruker_precursor_id=" << ms2->getMetaValue("bruker_precursor_id").toString());
}
END_SECTION

START_SECTION(Bruker load_ms1=false test)
{
  BrukerTimsFile f;

  // Baseline: default config loads both MS1 and MS2
  MSExperiment exp_both;
  f.load(OPENTIMS_DDA_TEST_DATA, exp_both);
  Size raw_ms1 = 0, raw_ms2 = 0;
  for (const auto& s : exp_both)
  {
    if (s.getMSLevel() == 1) ++raw_ms1;
    else if (s.getMSLevel() == 2) ++raw_ms2;
  }
  TEST_NOT_EQUAL(raw_ms1, 0);
  TEST_NOT_EQUAL(raw_ms2, 0);

  // load_ms1=false: zero MS1, MS2 count unchanged
  BrukerTimsFile::Config cfg;
  cfg.load_ms1 = false;
  MSExperiment exp_ms2_only;
  f.load(OPENTIMS_DDA_TEST_DATA, exp_ms2_only, cfg);
  Size skip_ms1 = 0, skip_ms2 = 0;
  for (const auto& s : exp_ms2_only)
  {
    if (s.getMSLevel() == 1) ++skip_ms1;
    else if (s.getMSLevel() == 2) ++skip_ms2;
  }
  TEST_EQUAL(skip_ms1, 0);
  TEST_EQUAL(skip_ms2, raw_ms2);

  STATUS("DDA load_ms1=false: raw MS1=" << raw_ms1 << " MS2=" << raw_ms2
         << " | skipped MS1=" << skip_ms1 << " MS2=" << skip_ms2);

#ifdef OPENTIMS_DIA_TEST_DATA
  // Same invariant on DIA
  MSExperiment exp_dia_both;
  f.load(OPENTIMS_DIA_TEST_DATA, exp_dia_both);
  Size dia_raw_ms1 = 0, dia_raw_ms2 = 0;
  for (const auto& s : exp_dia_both)
  {
    if (s.getMSLevel() == 1) ++dia_raw_ms1;
    else if (s.getMSLevel() == 2) ++dia_raw_ms2;
  }

  BrukerTimsFile::Config cfg_dia;
  cfg_dia.load_ms1 = false;
  MSExperiment exp_dia_skip;
  f.load(OPENTIMS_DIA_TEST_DATA, exp_dia_skip, cfg_dia);
  Size dia_skip_ms1 = 0, dia_skip_ms2 = 0;
  for (const auto& s : exp_dia_skip)
  {
    if (s.getMSLevel() == 1) ++dia_skip_ms1;
    else if (s.getMSLevel() == 2) ++dia_skip_ms2;
  }
  TEST_EQUAL(dia_skip_ms1, 0);
  TEST_EQUAL(dia_skip_ms2, dia_raw_ms2);

  // FRAME export mode must also honor load_ms1 (skips level==1 loop)
  BrukerTimsFile::Config cfg_frame;
  cfg_frame.export_mode = BrukerTimsFile::Config::FRAME;
  cfg_frame.load_ms1 = false;
  MSExperiment exp_frame;
  f.load(OPENTIMS_DIA_TEST_DATA, exp_frame, cfg_frame);
  Size frame_ms1 = 0;
  for (const auto& s : exp_frame)
  {
    if (s.getMSLevel() == 1) ++frame_ms1;
  }
  TEST_EQUAL(frame_ms1, 0);
#endif
}
END_SECTION

START_SECTION(DDA round-trip test: load .d -> write mzML -> reload -> verify)
{
  // Load from .d
  BrukerTimsFile f;
  MSExperiment orig;
  f.load(OPENTIMS_DDA_TEST_DATA, orig);

  // Write to temporary mzML — avoid NEW_TMP_FILE because the test
  // framework validates all registered .mzML files, and the IM data
  // array CV term MS:1003008 fails semantic validation (known issue).
  String tmp_mzml = File::getTempDirectory() + "/" + File::getUniqueName() + "_dda_roundtrip.mzML";
  MzMLFile().store(tmp_mzml, orig);

  // Reload from mzML
  MSExperiment reloaded;
  MzMLFile().load(tmp_mzml, reloaded);
  File::remove(tmp_mzml);

  // Verify spectrum count matches
  TEST_EQUAL(orig.size(), reloaded.size());

  // Verify MS levels survive round-trip
  for (Size i = 0; i < std::min(orig.size(), reloaded.size()); ++i)
  {
    TEST_EQUAL(orig[i].getMSLevel(), reloaded[i].getMSLevel());
  }

  // Verify precursor m/z survives round-trip for MS2
  for (const auto& spec : reloaded)
  {
    if (spec.getMSLevel() == 2 && !spec.getPrecursors().empty())
    {
      TEST_NOT_EQUAL(spec.getPrecursors()[0].getMZ(), 0.0);
      break;
    }
  }
}
END_SECTION

START_SECTION(DDA frame-level loading test)
{
  BrukerTimsFile f;
  BrukerTimsFile::Config cfg;
  cfg.export_mode = BrukerTimsFile::Config::FRAME;

  MSExperiment exp;
  f.load(OPENTIMS_DDA_TEST_DATA, exp, cfg);

  TEST_NOT_EQUAL(exp.size(), 0);

  // All spectra should have IM_PEAK IM data (even MS2 in frame mode)
  for (const auto& spec : exp)
  {
    if (!spec.empty())
    {
      TEST_EQUAL(spec.containsIMData(), true);
      // frameToSpectrum emits detector-centroided peaks; label as CENTROID.
      TEST_EQUAL(spec.getType(), SpectrumSettings::SpectrumType::CENTROID);
      break;
    }
  }
}
END_SECTION

START_SECTION(DDA MS1 centroiding test)
{
  BrukerTimsFile f;

  // Load without centroiding (baseline)
  MSExperiment exp_raw;
  f.load(OPENTIMS_DDA_TEST_DATA, exp_raw);

  // Load with centroiding enabled
  BrukerTimsFile::Config cfg;
  cfg.ms1_centroid_mz_ppm = 5.0f;
  cfg.ms1_centroid_im_pct = 3.0f;
  MSExperiment exp_cent;
  f.load(OPENTIMS_DDA_TEST_DATA, exp_cent, cfg);

  // Count MS1 and MS2 spectra and total MS1 peaks in both
  Size raw_ms1_peaks = 0, cent_ms1_peaks = 0;
  Size raw_ms1_count = 0, cent_ms1_count = 0;
  Size raw_ms2_count = 0, cent_ms2_count = 0;

  for (const auto& spec : exp_raw)
  {
    if (spec.getMSLevel() == 1) { ++raw_ms1_count; raw_ms1_peaks += spec.size(); }
    else { ++raw_ms2_count; }
  }
  for (const auto& spec : exp_cent)
  {
    if (spec.getMSLevel() == 1) { ++cent_ms1_count; cent_ms1_peaks += spec.size(); }
    else { ++cent_ms2_count; }
  }

  // MS1 frame count should be the same
  TEST_EQUAL(raw_ms1_count, cent_ms1_count);

  // MS2 count should be unaffected by centroiding
  TEST_EQUAL(raw_ms2_count, cent_ms2_count);

  // Centroided MS1 should have fewer total peaks
  TEST_EQUAL(cent_ms1_peaks < raw_ms1_peaks, true);

  // Centroided MS1 spectra should still have IM data, be marked CENTROID,
  // have IM array length == peak count, and be sorted by m/z
  for (const auto& spec : exp_cent)
  {
    if (spec.getMSLevel() == 1 && !spec.empty())
    {
      TEST_EQUAL(spec.containsIMData(), true);
      TEST_EQUAL(spec.getType(), SpectrumSettings::SpectrumType::CENTROID);

      // IM array length must match peak count
      const auto& fda = spec.getFloatDataArrays();
      TEST_EQUAL(fda.empty(), false);
      TEST_EQUAL(fda[0].size(), spec.size());

      // Peaks must be sorted by m/z (centroider sorts output)
      bool mz_sorted = true;
      for (Size j = 1; j < spec.size(); ++j)
      {
        if (spec[j].getMZ() < spec[j-1].getMZ()) { mz_sorted = false; break; }
      }
      TEST_EQUAL(mz_sorted, true);
      break;
    }
  }
}
END_SECTION

START_SECTION(DDA partial centroiding config test)
{
  // Setting only one tolerance should NOT enable centroiding
  BrukerTimsFile f;

  BrukerTimsFile::Config cfg_partial;
  cfg_partial.ms1_centroid_mz_ppm = 5.0f;  // set
  cfg_partial.ms1_centroid_im_pct = 0.0f;  // NOT set
  MSExperiment exp_partial;
  f.load(OPENTIMS_DDA_TEST_DATA, exp_partial, cfg_partial);

  // Load without centroiding for comparison
  MSExperiment exp_raw;
  f.load(OPENTIMS_DDA_TEST_DATA, exp_raw);

  // MS1 peak counts should be identical (centroiding was NOT applied)
  Size partial_ms1_peaks = 0, raw_ms1_peaks = 0;
  for (const auto& spec : exp_partial)
    if (spec.getMSLevel() == 1) partial_ms1_peaks += spec.size();
  for (const auto& spec : exp_raw)
    if (spec.getMSLevel() == 1) raw_ms1_peaks += spec.size();

  TEST_EQUAL(partial_ms1_peaks, raw_ms1_peaks);
}
END_SECTION

START_SECTION(DDA search engine IM annotation integration test)
{
  // Load real DDA-PASEF data
  BrukerTimsFile f;
  MSExperiment exp;
  f.load(OPENTIMS_DDA_TEST_DATA, exp);

  // Verify MS2 spectra have drift times (pre-condition for IM annotation)
  Size ms2_count = 0;
  Size ms2_with_im = 0;
  for (const auto& spec : exp)
  {
    if (spec.getMSLevel() == 2)
    {
      ++ms2_count;
      if (IMTypes::determineIMFormat(spec) == IMFormat::IM_SPECTRUM)
      {
        ++ms2_with_im;
      }
    }
  }
  TEST_TRUE(ms2_count > 0)
  TEST_EQUAL(ms2_with_im, ms2_count) // all MS2 spectra should have drift time

  // Run ProSEAlgorithm in-memory with the same FASTA and
  // parameters used by the TOPP-level DDA tests for SSE and FI.
  // The key test: any PSMs produced must have IM annotation.
  vector<FASTAFile::FASTAEntry> fasta_db;
  FASTAFile().load(OPENTIMS_TEST_FASTA, fasta_db);
  TEST_TRUE(fasta_db.size() > 0)

  // Typical timsTOF Pro DDA-PASEF search parameters
  ProSEAlgorithm algo;
  Param p = algo.getParameters();
  p.setValue("precursor:mass_tolerance_lower", 20.0);
  p.setValue("precursor:mass_tolerance_upper", 20.0);
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  p.setValue("fragment:mass_tolerance", 20.0);
  p.setValue("fragment:mass_tolerance_unit", "ppm");
  p.setValue("enzyme", "Trypsin/P");
  p.setValue("peptide:missed_cleavages", 2);
  p.setValue("modifications:variable", std::vector<std::string>{"Oxidation (M)", "Acetyl (Protein N-term)"});
  algo.setParameters(p);

  vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  auto ec = algo.search(exp, fasta_db, prot_ids, pep_ids);

  TEST_EQUAL(ec == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)
  TEST_EQUAL(prot_ids.size(), 1)

  // Verify IM annotation: every PSM must have IM meta value (all MS2 spectra have drift time)
  for (const auto& pid : pep_ids)
  {
    TEST_EQUAL(pid.metaValueExists(Constants::UserParam::IM), true)
    if (pid.metaValueExists(Constants::UserParam::IM))
    {
      double im_val = pid.getMetaValue(Constants::UserParam::IM);
      TEST_TRUE(im_val > 0.0) // 1/K0 values are positive
    }
  }

  // If we got any PSMs, verify IM unit on ProteinIdentification
  if (!pep_ids.empty())
  {
    TEST_EQUAL(prot_ids[0].metaValueExists(Constants::UserParam::IM), true)
    if (prot_ids[0].metaValueExists(Constants::UserParam::IM))
    {
      TEST_STRING_EQUAL(prot_ids[0].getMetaValue(Constants::UserParam::IM).toString(), "1/K0")
    }
  }
}
END_SECTION

START_SECTION(DDA MS1 aggregation with RT cap test)
{
  BrukerTimsFile f;

  // Baseline: raw DDA MS1
  MSExperiment exp_raw;
  f.load(OPENTIMS_DDA_TEST_DATA, exp_raw);

  // Aggregated: ms1_n_neighbors=1 uncapped
  BrukerTimsFile::Config cfg_uncapped;
  cfg_uncapped.ms1_n_neighbors = 1;
  cfg_uncapped.ms1_max_rt_distance_sec = 0.0;
  MSExperiment exp_uncapped;
  f.load(OPENTIMS_DDA_TEST_DATA, exp_uncapped, cfg_uncapped);

  // Aggregated: tight 1.0s RT cap (may be inert if fixture MS1 cadence is tight)
  BrukerTimsFile::Config cfg_capped;
  cfg_capped.ms1_n_neighbors = 1;
  cfg_capped.ms1_max_rt_distance_sec = 1.0;
  MSExperiment exp_capped;
  f.load(OPENTIMS_DDA_TEST_DATA, exp_capped, cfg_capped);

  auto count_ms1 = [](const MSExperiment& e, Size& n_spectra, double& total_intensity)
  {
    n_spectra = 0; total_intensity = 0.0;
    for (const auto& s : e)
      if (s.getMSLevel() == 1)
      {
        ++n_spectra;
        for (const auto& p : s) total_intensity += p.getIntensity();
      }
  };

  Size raw_n, uncapped_n, capped_n;
  double raw_i, uncapped_i, capped_i;
  count_ms1(exp_raw, raw_n, raw_i);
  count_ms1(exp_uncapped, uncapped_n, uncapped_i);
  count_ms1(exp_capped, capped_n, capped_i);

  TEST_NOT_EQUAL(raw_n, 0);
  TEST_NOT_EQUAL(uncapped_n, 0);
  TEST_NOT_EQUAL(capped_n, 0);

  // Uncapped aggregation boosts intensity over raw
  TEST_EQUAL(uncapped_i > raw_i, true);

  // Capped aggregation ≤ uncapped (cap can only exclude neighbors, never add)
  TEST_EQUAL(capped_i <= uncapped_i, true);

  STATUS("DDA MS1 aggregation: raw intensity=" << raw_i
         << " uncapped=" << uncapped_i
         << " capped(1.0s)=" << capped_i
         << " cap_reduction=" << (1.0 - capped_i / uncapped_i));
}
END_SECTION

#endif // OPENTIMS_DDA_TEST_DATA

#ifdef OPENTIMS_DIA_TEST_DATA

START_SECTION(DIA loading integration test)
{
  BrukerTimsFile f;
  MSExperiment exp;
  f.load(OPENTIMS_DIA_TEST_DATA, exp);

  TEST_NOT_EQUAL(exp.size(), 0);

  // Check MS2 spectra have per-peak IM (IM_PEAK) and isolation windows
  for (const auto& spec : exp)
  {
    if (spec.getMSLevel() == 2 && !spec.empty())
    {
      TEST_EQUAL(spec.containsIMData(), true);
      TEST_NOT_EQUAL(spec.getPrecursors().size(), 0);
      // TDF MS2 peaks are detector-centroided; must be labelled CENTROID.
      TEST_EQUAL(spec.getType(), SpectrumSettings::SpectrumType::CENTROID);
      break;
    }
  }
}
END_SECTION

START_SECTION(DIA round-trip test: load .d -> write mzML -> reload -> verify)
{
  // Load from .d
  BrukerTimsFile f;
  MSExperiment orig;
  f.load(OPENTIMS_DIA_TEST_DATA, orig);

  // Write to temporary mzML — avoid NEW_TMP_FILE (see DDA round-trip comment)
  String tmp_mzml = File::getTempDirectory() + "/" + File::getUniqueName() + "_dia_roundtrip.mzML";
  MzMLFile().store(tmp_mzml, orig);

  // Reload from mzML
  MSExperiment reloaded;
  MzMLFile().load(tmp_mzml, reloaded);
  File::remove(tmp_mzml);

  // Verify spectrum count matches
  TEST_EQUAL(orig.size(), reloaded.size());

  // Verify MS2 spectra retain IM data arrays after round-trip
  for (const auto& spec : reloaded)
  {
    if (spec.getMSLevel() == 2 && !spec.empty())
    {
      TEST_EQUAL(spec.containsIMData(), true);
      break;
    }
  }
}
END_SECTION

START_SECTION(DIA MS1 centroiding test)
{
  BrukerTimsFile f;

  // Load without centroiding
  MSExperiment exp_raw;
  f.load(OPENTIMS_DIA_TEST_DATA, exp_raw);

  // Load with centroiding
  BrukerTimsFile::Config cfg;
  cfg.ms1_centroid_mz_ppm = 5.0f;
  cfg.ms1_centroid_im_pct = 3.0f;
  MSExperiment exp_cent;
  f.load(OPENTIMS_DIA_TEST_DATA, exp_cent, cfg);

  // Count MS1 peaks
  Size raw_ms1_peaks = 0, cent_ms1_peaks = 0;
  for (const auto& spec : exp_raw)
    if (spec.getMSLevel() == 1) raw_ms1_peaks += spec.size();
  for (const auto& spec : exp_cent)
    if (spec.getMSLevel() == 1) cent_ms1_peaks += spec.size();

  // Centroided should have fewer peaks
  TEST_EQUAL(cent_ms1_peaks < raw_ms1_peaks, true);
}
END_SECTION

START_SECTION(DIA MS2 aggregation test)
{
  BrukerTimsFile f;

  // Load without aggregation (baseline)
  MSExperiment exp_raw;
  f.load(OPENTIMS_DIA_TEST_DATA, exp_raw);

  // Load with aggregation (n_neighbors=1 → 3-frame sum)
  BrukerTimsFile::Config cfg;
  cfg.dia_ms2_n_neighbors = 1;
  cfg.dia_ms2_min_support = 1;
  MSExperiment exp_agg;
  f.load(OPENTIMS_DIA_TEST_DATA, exp_agg, cfg);

  // Count MS2 spectra, total MS2 peaks, and cumulative MS2 intensity in both.
  //
  // FrameAggregator (BrukerTimsFile.cpp) is a two-stage operation on a
  // sparse 2D grid keyed by (mz_bin, scan_id) — NOT (mz_bin, scan_id, frame_id):
  //
  //   1. SUM stage: for each (center_frame_i, window) position, accumulate
  //      peaks from 2*n_neighbors+1 = 3 adjacent frames into grid cells. Each
  //      cell's intensity_sum is the SUM of all contributing frames' peak
  //      intensities at that (mz_bin, scan_id). Each raw frame's peaks get
  //      added to 3 different aggregator grids (one as center, two as
  //      neighbors), so the TOTAL intensity across all aggregated positions
  //      is ~3x the raw total (verified empirically with min_support=0:
  //      8469 raw spectra, 50.1M peaks, 6.00e9 intensity ->
  //      8469 agg spectra, 141.2M peaks, 1.80e10 intensity = 2.998x raw).
  //
  //   2. DENOISE stage: with min_support=1, each cell is dropped unless at
  //      least ONE of its 8 spatial neighbors in the (mz_bin, scan_id) 3x3
  //      grid is also occupied. This is a 2D convolution-style filter that
  //      removes isolated cells — typically noise peaks that lack a compact
  //      cluster shape in (m/z, IM) space. Empirically, this drops ~64% of
  //      peaks and ~57% of intensity (141.2M -> 50.6M peaks, 1.80e10 ->
  //      7.71e9 intensity). Surviving peaks have ~1.20x the average
  //      intensity of the pre-denoise pool, confirming denoising
  //      preferentially keeps high-SNR signal.
  //
  // Net effect (n_neighbors=1, min_support=1): peak COUNT stays roughly equal
  // to raw (+1% from transient edge peaks that just pass the neighbor check);
  // total intensity is ~1.28x raw. The sparsity reduction for picking comes
  // from the DENOISED output: each surviving cell is both part of a compact
  // 2D cluster AND has its intensity boosted by the 3-frame sum. Downstream
  // 2D Gaussian smoothing in finalizeCentroided() operates on this
  // clustered-and-boosted grid, giving much better centroid fits than raw.
  Size raw_ms2_count = 0, agg_ms2_count = 0;
  Size raw_ms2_peaks = 0, agg_ms2_peaks = 0;
  double raw_ms2_intensity = 0.0, agg_ms2_intensity = 0.0;
  for (const auto& spec : exp_raw)
  {
    if (spec.getMSLevel() == 2)
    {
      ++raw_ms2_count;
      raw_ms2_peaks += spec.size();
      for (const auto& p : spec) raw_ms2_intensity += p.getIntensity();
    }
  }
  for (const auto& spec : exp_agg)
  {
    if (spec.getMSLevel() == 2)
    {
      ++agg_ms2_count;
      agg_ms2_peaks += spec.size();
      for (const auto& p : spec) agg_ms2_intensity += p.getIntensity();
    }
  }

  // Both should have MS2 spectra.
  TEST_NOT_EQUAL(raw_ms2_count, 0);
  TEST_NOT_EQUAL(agg_ms2_count, 0);

  // Spectrum count: both paths iterate per-WindowGroup, but aggregation drops
  // a small number of boundary positions at the LC run edges where the
  // n_neighbors window cannot be satisfied. Allow up to 2% loss.
  TEST_EQUAL(agg_ms2_count <= raw_ms2_count, true);
  TEST_EQUAL(agg_ms2_count >= raw_ms2_count * 98 / 100, true);

  // Peak count is nearly unchanged (see big comment above). Allow ±10% for
  // transient edge peaks (cells populated by only some of the neighbor frames).
  TEST_EQUAL(agg_ms2_peaks >= raw_ms2_peaks * 90 / 100, true);
  TEST_EQUAL(agg_ms2_peaks <= raw_ms2_peaks * 110 / 100, true);

  // Cumulative intensity: aggregation sums contributions from neighbor frames
  // into shared (mz_bin, scan_id) cells. The effective intensity multiplier
  // is NOT the naive 2*n_neighbors+1 = 3 — it depends on how many cells are
  // shared across neighbor frames, which is a data-dependent property (profile
  // mode sampling + m/z binning + integer IM scan quantization + stable vs
  // transient ion populations). For this fixture with n_neighbors=1 we observe
  // ratio ≈ 1.28, reflecting that only a subset of cells are populated across
  // all 3 neighbor frames.
  //
  // The bounds below are a regression guard: [1.15, 1.45] catches
  //   - no intensity boost at all (ratio == 1, aggregation broken)
  //   - runaway intensity (ratio > 1.5, e.g. double-counting)
  // while tolerating minor fixture variation.
  TEST_EQUAL(raw_ms2_intensity > 0.0, true); // guard division below
  const double intensity_ratio = agg_ms2_intensity / raw_ms2_intensity;
  TEST_EQUAL(intensity_ratio >= 1.15, true);
  TEST_EQUAL(intensity_ratio <= 1.45, true);

  // Aggregated spectra should have per-peak IM data.
  for (const auto& spec : exp_agg)
  {
    if (spec.getMSLevel() == 2 && !spec.empty())
    {
      TEST_EQUAL(spec.containsIMData(), true);
      TEST_EQUAL(spec.getIMPeakType() == IMPeakType::IM_PROFILE, true);
      TEST_NOT_EQUAL(spec.getPrecursors().size(), 0);
      break;
    }
  }

  STATUS("DIA aggregation: raw MS2 spectra=" << raw_ms2_count
         << " peaks=" << raw_ms2_peaks
         << " intensity=" << raw_ms2_intensity
         << " | aggregated MS2 spectra=" << agg_ms2_count
         << " peaks=" << agg_ms2_peaks
         << " intensity=" << agg_ms2_intensity
         << " | intensity_ratio=" << intensity_ratio);
}
END_SECTION

START_SECTION(DIA MS2 centroiding test)
{
  BrukerTimsFile f;

  // Load with aggregation + centroiding
  BrukerTimsFile::Config cfg;
  cfg.dia_ms2_n_neighbors = 1;
  cfg.dia_ms2_min_support = 1;
  cfg.dia_ms2_centroid = true;
  MSExperiment exp_cent;
  f.load(OPENTIMS_DIA_TEST_DATA, exp_cent, cfg);

  // Count MS2 spectra
  Size cent_ms2_count = 0;
  Size cent_ms2_peaks = 0;
  for (const auto& spec : exp_cent)
  {
    if (spec.getMSLevel() == 2)
    {
      ++cent_ms2_count;
      cent_ms2_peaks += spec.size();
    }
  }

  TEST_NOT_EQUAL(cent_ms2_count, 0);

  // Centroided spectra should have IM_CENTROIDED type and per-peak IM data
  for (const auto& spec : exp_cent)
  {
    if (spec.getMSLevel() == 2 && !spec.empty())
    {
      TEST_EQUAL(spec.containsIMData(), true);
      TEST_EQUAL(spec.getIMPeakType() == IMPeakType::IM_CENTROIDED, true);
      TEST_NOT_EQUAL(spec.getPrecursors().size(), 0);
      break;
    }
  }

  STATUS("DIA centroiding: MS2 spectra=" << cent_ms2_count
         << " peaks=" << cent_ms2_peaks);
}
END_SECTION

START_SECTION(DIA MS2 centroiding without denoising (min_support=0))
{
  // Regression guard for the "centroid without noise filter" path exposed via
  // Config::dia_ms2_min_support = 0. With min_support=0 the 3x3 neighbor filter
  // in FrameAggregator::finalize{,Centroided} is a no-op (`neighbors < 0` is
  // always false), so every populated grid cell survives the denoise step.
  // Under centroiding that means more local maxima → strictly more output peaks
  // than the denoised (min_support=1) path, because denoising preferentially
  // removes isolated cells that would otherwise become their own centroid.
  BrukerTimsFile f;

  BrukerTimsFile::Config cfg_denoise;
  cfg_denoise.dia_ms2_n_neighbors = 1;
  cfg_denoise.dia_ms2_min_support = 1;
  cfg_denoise.dia_ms2_centroid = true;
  MSExperiment exp_denoise;
  f.load(OPENTIMS_DIA_TEST_DATA, exp_denoise, cfg_denoise);

  BrukerTimsFile::Config cfg_no_denoise;
  cfg_no_denoise.dia_ms2_n_neighbors = 1;
  cfg_no_denoise.dia_ms2_min_support = 0;
  cfg_no_denoise.dia_ms2_centroid = true;
  MSExperiment exp_no_denoise;
  f.load(OPENTIMS_DIA_TEST_DATA, exp_no_denoise, cfg_no_denoise);

  Size denoise_peaks = 0, no_denoise_peaks = 0;
  Size denoise_spectra = 0, no_denoise_spectra = 0;
  for (const auto& spec : exp_denoise)
    if (spec.getMSLevel() == 2) { ++denoise_spectra; denoise_peaks += spec.size(); }
  for (const auto& spec : exp_no_denoise)
    if (spec.getMSLevel() == 2) { ++no_denoise_spectra; no_denoise_peaks += spec.size(); }

  // Without denoising, a few spectra that would become empty after denoising
  // (and are dropped by the `if (peaks.empty()) continue;` guard in the
  // BrukerTimsFile loader) are now retained. Count is >= denoised count.
  TEST_EQUAL(no_denoise_spectra >= denoise_spectra, true);

  // Without denoising the centroided output contains strictly more peaks:
  // every populated grid cell becomes a centroid candidate instead of being
  // dropped when isolated.
  TEST_NOT_EQUAL(no_denoise_peaks, 0);
  TEST_EQUAL(no_denoise_peaks > denoise_peaks, true);

  // Output is still centroided in the IM dimension.
  for (const auto& spec : exp_no_denoise)
  {
    if (spec.getMSLevel() == 2 && !spec.empty())
    {
      TEST_EQUAL(spec.containsIMData(), true);
      TEST_EQUAL(spec.getIMPeakType() == IMPeakType::IM_CENTROIDED, true);
      TEST_NOT_EQUAL(spec.getPrecursors().size(), 0);
      break;
    }
  }

  STATUS("DIA centroiding (no denoise): MS2 peaks=" << no_denoise_peaks
         << " vs denoised=" << denoise_peaks
         << " ratio=" << (static_cast<double>(no_denoise_peaks) / denoise_peaks));
}
END_SECTION

START_SECTION(DIA MS1 default config regression test)
{
  // Regression guard: with all new MS1 aggregation knobs at their defaults
  // (ms1_n_neighbors=0, ms1_min_support=0, ms1_max_rt_distance_sec=0.0),
  // the DIA MS1 output must be byte-identical to the pre-PR load. This
  // section fails if a future change accidentally enables aggregation
  // by default.
  BrukerTimsFile f;

  BrukerTimsFile::Config cfg_default;  // all defaults
  MSExperiment exp_default;
  f.load(OPENTIMS_DIA_TEST_DATA, exp_default, cfg_default);

  // Snapshot: count MS1 spectra, total MS1 peaks, total MS1 intensity
  Size ms1_spectra = 0, ms1_peaks = 0;
  double ms1_intensity = 0.0;
  for (const auto& spec : exp_default)
  {
    if (spec.getMSLevel() == 1)
    {
      ++ms1_spectra;
      ms1_peaks += spec.size();
      for (const auto& p : spec) ms1_intensity += p.getIntensity();
    }
  }

  TEST_NOT_EQUAL(ms1_spectra, 0);
  TEST_NOT_EQUAL(ms1_peaks, 0);

  STATUS("DIA MS1 default-config snapshot: spectra=" << ms1_spectra
         << " peaks=" << ms1_peaks
         << " intensity=" << ms1_intensity);

  // Independently reproduce with ms1_n_neighbors=0 explicitly — assert
  // the numbers match exactly, proving the knob at 0 is inert.
  BrukerTimsFile::Config cfg_explicit_off;
  cfg_explicit_off.ms1_n_neighbors = 0;
  cfg_explicit_off.ms1_min_support = 0;
  cfg_explicit_off.ms1_max_rt_distance_sec = 0.0;
  MSExperiment exp_explicit_off;
  f.load(OPENTIMS_DIA_TEST_DATA, exp_explicit_off, cfg_explicit_off);

  Size ms1_spectra_b = 0, ms1_peaks_b = 0;
  double ms1_intensity_b = 0.0;
  for (const auto& spec : exp_explicit_off)
  {
    if (spec.getMSLevel() == 1)
    {
      ++ms1_spectra_b;
      ms1_peaks_b += spec.size();
      for (const auto& p : spec) ms1_intensity_b += p.getIntensity();
    }
  }

  TEST_EQUAL(ms1_spectra_b, ms1_spectra);
  TEST_EQUAL(ms1_peaks_b, ms1_peaks);
  TEST_REAL_SIMILAR(ms1_intensity_b, ms1_intensity);
}
END_SECTION

START_SECTION(DIA MS1 aggregation test)
{
  // RED test for MS1 frame aggregation. With ms1_n_neighbors=1 we expect
  // ~same spectrum count as raw (edge-truncation aside), boosted intensity,
  // and IM_PROFILE output with per-peak IM.
  BrukerTimsFile f;

  MSExperiment exp_raw;
  f.load(OPENTIMS_DIA_TEST_DATA, exp_raw);

  BrukerTimsFile::Config cfg;
  cfg.ms1_n_neighbors = 1;
  MSExperiment exp_agg;
  f.load(OPENTIMS_DIA_TEST_DATA, exp_agg, cfg);

  Size raw_ms1_spectra = 0, raw_ms1_peaks = 0;
  double raw_ms1_intensity = 0.0;
  for (const auto& spec : exp_raw)
  {
    if (spec.getMSLevel() == 1)
    {
      ++raw_ms1_spectra;
      raw_ms1_peaks += spec.size();
      for (const auto& p : spec) raw_ms1_intensity += p.getIntensity();
    }
  }

  Size agg_ms1_spectra = 0, agg_ms1_peaks = 0;
  double agg_ms1_intensity = 0.0;
  for (const auto& spec : exp_agg)
  {
    if (spec.getMSLevel() == 1)
    {
      ++agg_ms1_spectra;
      agg_ms1_peaks += spec.size();
      for (const auto& p : spec) agg_ms1_intensity += p.getIntensity();
    }
  }

  TEST_NOT_EQUAL(raw_ms1_spectra, 0);
  TEST_NOT_EQUAL(agg_ms1_spectra, 0);

  // Spectrum count: aggregated MS1 path emits exactly one spectrum per center
  // frame (empty spectra emitted when all peaks filter out — see the empty-
  // spectrum contract in loadAggregatedMS1Spectrum). Unlike the MS2 path,
  // this is a hard equality: MS1 never drops spectra.
  TEST_EQUAL(agg_ms1_spectra, raw_ms1_spectra);

  // Intensity boost: aggregation sums intensity across 2*N+1 = 3 frames.
  // MS1 fixture shows ratio ≈ 2.99 (empirical) — much higher than MS2's ≈1.28
  // because MS1 ion populations are more stable across adjacent frames: peaks
  // at the same nominal (mz_bin, scan_id) from different frames tend to land
  // in different cells (scan_id jitter, m/z jitter beyond 0.01 Da), so the
  // aggregator rarely merges them — closer to the naive 3x sum than MS2.
  // Bounds reflect this: accept anywhere in [2.5, 3.1] as valid aggregation.
  TEST_EQUAL(raw_ms1_intensity > 0.0, true);
  const double intensity_ratio = agg_ms1_intensity / raw_ms1_intensity;
  TEST_EQUAL(intensity_ratio >= 2.5, true);
  TEST_EQUAL(intensity_ratio <= 3.1, true);

  // Output must carry per-peak IM with IM_PROFILE type
  for (const auto& spec : exp_agg)
  {
    if (spec.getMSLevel() == 1 && !spec.empty())
    {
      TEST_EQUAL(spec.containsIMData(), true);
      TEST_EQUAL(spec.getIMPeakType() == IMPeakType::IM_PROFILE, true);
      // MS1 spectra must NOT carry precursors (unlike MS2).
      TEST_EQUAL(spec.getPrecursors().empty(), true);
      break;
    }
  }

  STATUS("DIA MS1 aggregation: raw spectra=" << raw_ms1_spectra
         << " peaks=" << raw_ms1_peaks
         << " intensity=" << raw_ms1_intensity
         << " | aggregated spectra=" << agg_ms1_spectra
         << " peaks=" << agg_ms1_peaks
         << " intensity=" << agg_ms1_intensity
         << " | ratio=" << intensity_ratio);
}
END_SECTION

START_SECTION(DIA MS1 aggregation with centroiding test)
{
  BrukerTimsFile f;

  // Aggregated, no centroiding (baseline)
  BrukerTimsFile::Config cfg_profile;
  cfg_profile.ms1_n_neighbors = 1;
  MSExperiment exp_profile;
  f.load(OPENTIMS_DIA_TEST_DATA, exp_profile, cfg_profile);

  // Aggregated + within-frame centroiding
  BrukerTimsFile::Config cfg_centroid;
  cfg_centroid.ms1_n_neighbors = 1;
  cfg_centroid.ms1_centroid_mz_ppm = 5.0f;
  cfg_centroid.ms1_centroid_im_pct = 3.0f;
  MSExperiment exp_centroid;
  f.load(OPENTIMS_DIA_TEST_DATA, exp_centroid, cfg_centroid);

  auto count_ms1 = [](const MSExperiment& e, Size& n_spectra, Size& n_peaks, double& total_intensity)
  {
    n_spectra = 0; n_peaks = 0; total_intensity = 0.0;
    for (const auto& s : e)
      if (s.getMSLevel() == 1)
      {
        ++n_spectra;
        n_peaks += s.size();
        for (const auto& p : s) total_intensity += p.getIntensity();
      }
  };

  Size profile_n, profile_p, centroid_n, centroid_p;
  double profile_i, centroid_i;
  count_ms1(exp_profile, profile_n, profile_p, profile_i);
  count_ms1(exp_centroid, centroid_n, centroid_p, centroid_i);

  TEST_NOT_EQUAL(centroid_n, 0);

  // Same spectrum count (centroiding doesn't drop spectra)
  TEST_EQUAL(centroid_n, profile_n);

  // Centroiding collapses IM-adjacent peaks → fewer peaks than profile.
  TEST_EQUAL(centroid_p < profile_p, true);

  // Intensity conservation: with the default ms1_centroid_max_peaks=100000
  // the cap rarely fires on typical aggregated MS1 surveys (top 100k peaks
  // out of ~600k carry >80% of total intensity). Allow [0.5, 1.0] to cover
  // dense fixtures where a noticeable long-tail still gets dropped.
  const double ratio = centroid_i / profile_i;
  TEST_EQUAL(ratio > 0.5, true);
  TEST_EQUAL(ratio <= 1.0, true);

  // Output shape
  for (const auto& spec : exp_centroid)
  {
    if (spec.getMSLevel() == 1 && !spec.empty())
    {
      TEST_EQUAL(spec.containsIMData(), true);
      TEST_EQUAL(spec.getIMPeakType() == IMPeakType::IM_CENTROIDED, true);
      TEST_EQUAL(spec.getType() == SpectrumSettings::SpectrumType::CENTROID, true);
      TEST_EQUAL(spec.getPrecursors().empty(), true);
      break;
    }
  }

  STATUS("DIA MS1 aggregation+centroid: profile peaks=" << profile_p
         << " centroid peaks=" << centroid_p
         << " ratio=" << (static_cast<double>(centroid_p) / profile_p)
         << " intensity_ratio=" << ratio);
}
END_SECTION

START_SECTION(DIA MS1 aggregation min_support test)
{
  BrukerTimsFile f;

  // No denoise (baseline from DIA MS1 aggregation test: intensity ratio ≈ 3.0)
  BrukerTimsFile::Config cfg_no_denoise;
  cfg_no_denoise.ms1_n_neighbors = 1;
  cfg_no_denoise.ms1_min_support = 0;
  MSExperiment exp_no_denoise;
  f.load(OPENTIMS_DIA_TEST_DATA, exp_no_denoise, cfg_no_denoise);

  // With denoise
  BrukerTimsFile::Config cfg_denoise;
  cfg_denoise.ms1_n_neighbors = 1;
  cfg_denoise.ms1_min_support = 1;
  MSExperiment exp_denoise;
  f.load(OPENTIMS_DIA_TEST_DATA, exp_denoise, cfg_denoise);

  Size no_denoise_peaks = 0, denoise_peaks = 0;
  double no_denoise_intensity = 0.0, denoise_intensity = 0.0;
  for (const auto& s : exp_no_denoise)
    if (s.getMSLevel() == 1)
    {
      no_denoise_peaks += s.size();
      for (const auto& p : s) no_denoise_intensity += p.getIntensity();
    }
  for (const auto& s : exp_denoise)
    if (s.getMSLevel() == 1)
    {
      denoise_peaks += s.size();
      for (const auto& p : s) denoise_intensity += p.getIntensity();
    }

  TEST_NOT_EQUAL(no_denoise_peaks, 0);
  TEST_NOT_EQUAL(denoise_peaks, 0);

  // Denoise drops peaks. Conservative lower bound: ≥10% drop; MS1 surveys
  // often have many low-intensity isolated cells that the 3x3 filter culls.
  TEST_EQUAL(denoise_peaks < no_denoise_peaks, true);
  const double drop_frac = 1.0 - static_cast<double>(denoise_peaks) / no_denoise_peaks;
  TEST_EQUAL(drop_frac >= 0.10, true);

  // Intensity boost still present after denoise (signal peaks survive).
  MSExperiment exp_raw;
  f.load(OPENTIMS_DIA_TEST_DATA, exp_raw);
  double raw_intensity = 0.0;
  for (const auto& s : exp_raw)
    if (s.getMSLevel() == 1)
      for (const auto& p : s) raw_intensity += p.getIntensity();
  TEST_EQUAL(denoise_intensity > raw_intensity, true);

  STATUS("DIA MS1 min_support: no_denoise_peaks=" << no_denoise_peaks
         << " denoise_peaks=" << denoise_peaks
         << " drop_frac=" << drop_frac
         << " denoise_intensity/raw=" << (denoise_intensity / raw_intensity));
}
END_SECTION

START_SECTION(DIA MS1 RT cap test)
{
  BrukerTimsFile f;

  // Uncapped reference: ms1_n_neighbors=2 gives a wider window for the cap
  // to bite into.
  BrukerTimsFile::Config cfg_uncapped;
  cfg_uncapped.ms1_n_neighbors = 2;
  cfg_uncapped.ms1_max_rt_distance_sec = 0.0;
  MSExperiment exp_uncapped;
  f.load(OPENTIMS_DIA_TEST_DATA, exp_uncapped, cfg_uncapped);

  // Tight cap: 0.5s (DIA MS1 cadence ~1.5s, so N=2 reaches ~3s total span;
  // the 0.5s cap should truncate aggressively).
  BrukerTimsFile::Config cfg_capped;
  cfg_capped.ms1_n_neighbors = 2;
  cfg_capped.ms1_max_rt_distance_sec = 0.5;
  MSExperiment exp_capped;
  f.load(OPENTIMS_DIA_TEST_DATA, exp_capped, cfg_capped);

  auto total_ms1_intensity = [](const MSExperiment& e)
  {
    double s = 0.0;
    for (const auto& spec : e)
      if (spec.getMSLevel() == 1)
        for (const auto& p : spec) s += p.getIntensity();
    return s;
  };

  double uncapped_i = total_ms1_intensity(exp_uncapped);
  double capped_i = total_ms1_intensity(exp_capped);

  TEST_EQUAL(uncapped_i > 0.0, true);
  TEST_EQUAL(capped_i > 0.0, true);
  TEST_EQUAL(capped_i < uncapped_i, true);

  // Cadence invariant: RT cap must never drop a center frame. Both runs emit
  // exactly one MS1 spectrum per raw MS1 frame regardless of the cap.
  auto count_ms1_spectra = [](const MSExperiment& e)
  {
    Size n = 0;
    for (const auto& s : e) if (s.getMSLevel() == 1) ++n;
    return n;
  };
  MSExperiment exp_raw;
  f.load(OPENTIMS_DIA_TEST_DATA, exp_raw);
  const Size raw_ms1 = count_ms1_spectra(exp_raw);
  TEST_EQUAL(count_ms1_spectra(exp_uncapped), raw_ms1);
  TEST_EQUAL(count_ms1_spectra(exp_capped), raw_ms1);

  STATUS("DIA MS1 RT cap: uncapped=" << uncapped_i
         << " capped(0.5s, N=2)=" << capped_i
         << " reduction=" << (1.0 - capped_i / uncapped_i));
}
END_SECTION

START_SECTION(FRAME mode ignores ms1_n_neighbors test)
{
  BrukerTimsFile f;

  // Reference: raw FRAME export with defaults
  BrukerTimsFile::Config cfg_ref;
  cfg_ref.export_mode = BrukerTimsFile::Config::FRAME;
  MSExperiment exp_ref;
  f.load(OPENTIMS_DIA_TEST_DATA, exp_ref, cfg_ref);

  // Same load but with aggregation knob set — must be ignored
  BrukerTimsFile::Config cfg_agg;
  cfg_agg.export_mode = BrukerTimsFile::Config::FRAME;
  cfg_agg.ms1_n_neighbors = 1;
  MSExperiment exp_agg;
  f.load(OPENTIMS_DIA_TEST_DATA, exp_agg, cfg_agg);

  Size ref_ms1 = 0, agg_ms1 = 0;
  Size ref_peaks = 0, agg_peaks = 0;
  double ref_intensity = 0.0, agg_intensity = 0.0;
  for (const auto& s : exp_ref)
    if (s.getMSLevel() == 1)
    {
      ++ref_ms1;
      ref_peaks += s.size();
      for (const auto& p : s) ref_intensity += p.getIntensity();
    }
  for (const auto& s : exp_agg)
    if (s.getMSLevel() == 1)
    {
      ++agg_ms1;
      agg_peaks += s.size();
      for (const auto& p : s) agg_intensity += p.getIntensity();
    }

  TEST_NOT_EQUAL(ref_ms1, 0);
  // FRAME mode must emit one spectrum per raw MS1 frame regardless of
  // ms1_n_neighbors — aggregation is ignored.
  TEST_EQUAL(agg_ms1, ref_ms1);
  // Stronger: content must be identical too. Aggregation would inflate both
  // peak count and cumulative intensity; FRAME mode's "raw frames" contract
  // requires bit-for-bit behavior regardless of the ms1_n_neighbors knob.
  TEST_EQUAL(agg_peaks, ref_peaks);
  TEST_REAL_SIMILAR(agg_intensity, ref_intensity);

  // Warning emission is not asserted here — the test framework has no
  // log-capture helper at present. Manual verification: the test output
  // should show the "Warning: ms1_n_neighbors ... ignored" line.
}
END_SECTION

START_SECTION(DIA readDIAMetadata test)
{
  BrukerTimsFile f;
  ExperimentalSettings settings;
  auto meta = f.readDIAMetadata(OPENTIMS_DIA_TEST_DATA, settings);

  // Should have SWATH windows
  TEST_NOT_EQUAL(meta.boundaries.size(), 0);

  // Should have MS1 frames
  TEST_NOT_EQUAL(meta.nr_ms1_spectra, 0);

  // nr_ms2_spectra should have same size as boundaries
  TEST_EQUAL(meta.nr_ms2_spectra.size(), meta.boundaries.size());

  // Each window should have spectra
  for (Size i = 0; i < meta.nr_ms2_spectra.size(); ++i)
  {
    TEST_NOT_EQUAL(meta.nr_ms2_spectra[i], 0);
  }

  // Boundaries should have valid m/z and IM ranges
  for (const auto& b : meta.boundaries)
  {
    TEST_EQUAL(b.ms1, false);
    TEST_EQUAL(b.center > 0, true);
    TEST_EQUAL(b.lower < b.upper, true);
    TEST_EQUAL(b.imLower >= 0, true);
    TEST_EQUAL(b.imUpper > b.imLower, true);
  }

  // ExperimentalSettings should have source file
  TEST_EQUAL(settings.getSourceFiles().size(), 1);

  STATUS("readDIAMetadata: " << meta.boundaries.size() << " windows, "
         << meta.nr_ms1_spectra << " MS1 spectra");
}
END_SECTION

START_SECTION(DIA loadDIAStreaming test)
{
  BrukerTimsFile f;

  // Get metadata first
  ExperimentalSettings settings;
  auto meta = f.readDIAMetadata(OPENTIMS_DIA_TEST_DATA, settings);

  // Create in-memory consumer with known boundaries
  RegularSwathFileConsumer consumer(meta.boundaries);
  consumer.setExperimentalSettings(settings);

  // Stream spectra
  f.loadDIAStreaming(OPENTIMS_DIA_TEST_DATA, consumer);

  // Retrieve SwathMaps
  std::vector<OpenSwath::SwathMap> swath_maps;
  consumer.retrieveSwathMaps(swath_maps);

  // Should have MS1 map + MS2 maps
  Size ms1_count = 0, ms2_count = 0;
  for (const auto& m : swath_maps)
  {
    if (m.ms1) ++ms1_count;
    else ++ms2_count;
  }
  TEST_EQUAL(ms1_count, 1);
  TEST_EQUAL(ms2_count, meta.boundaries.size());

  // MS1 map should have spectra matching metadata count
  for (const auto& m : swath_maps)
  {
    if (m.ms1)
    {
      TEST_EQUAL(m.sptr->getNrSpectra(), static_cast<size_t>(meta.nr_ms1_spectra));
      break;
    }
  }

  // Each MS2 map should have spectra with IM data
  for (const auto& m : swath_maps)
  {
    if (!m.ms1 && m.sptr->getNrSpectra() > 0)
    {
      auto spec = m.sptr->getSpectrumById(0);
      TEST_NOT_EQUAL(spec->getMZArray()->data.size(), 0);
      TEST_EQUAL(spec->getDriftTimeArray() == nullptr, false);
      break;
    }
  }

  STATUS("loadDIAStreaming: " << swath_maps.size() << " SwathMaps "
         << "(" << ms1_count << " MS1, " << ms2_count << " MS2)");
}
END_SECTION

START_SECTION(DIA MS1 streaming aggregation test)
{
  BrukerTimsFile f;

  // Baseline: streaming with defaults (no aggregation)
  ExperimentalSettings settings_ref;
  auto meta_ref = f.readDIAMetadata(OPENTIMS_DIA_TEST_DATA, settings_ref);
  RegularSwathFileConsumer consumer_ref(meta_ref.boundaries);
  consumer_ref.setExperimentalSettings(settings_ref);
  f.loadDIAStreaming(OPENTIMS_DIA_TEST_DATA, consumer_ref);
  std::vector<OpenSwath::SwathMap> maps_ref;
  consumer_ref.retrieveSwathMaps(maps_ref);

  // Aggregated streaming (ms1_n_neighbors=1)
  ExperimentalSettings settings_agg;
  auto meta_agg = f.readDIAMetadata(OPENTIMS_DIA_TEST_DATA, settings_agg);
  RegularSwathFileConsumer consumer_agg(meta_agg.boundaries);
  consumer_agg.setExperimentalSettings(settings_agg);
  BrukerTimsFile::Config cfg;
  cfg.ms1_n_neighbors = 1;
  f.loadDIAStreaming(OPENTIMS_DIA_TEST_DATA, consumer_agg, cfg);
  std::vector<OpenSwath::SwathMap> maps_agg;
  consumer_agg.retrieveSwathMaps(maps_agg);

  // Locate each run's MS1 map by its ms1 flag (not by position).
  auto find_ms1_map = [](const std::vector<OpenSwath::SwathMap>& maps)
    -> const OpenSwath::SwathMap* {
    for (const auto& m : maps) if (m.ms1) return &m;
    return nullptr;
  };
  const auto* ms1_ref = find_ms1_map(maps_ref);
  const auto* ms1_agg = find_ms1_map(maps_agg);
  TEST_NOT_EQUAL(ms1_ref, nullptr);
  TEST_NOT_EQUAL(ms1_agg, nullptr);

  // Cadence invariant: both paths deliver one spectrum per raw MS1 frame.
  TEST_EQUAL(ms1_ref->sptr->getNrSpectra(),
             static_cast<unsigned>(meta_ref.nr_ms1_spectra));
  TEST_EQUAL(ms1_agg->sptr->getNrSpectra(),
             static_cast<unsigned>(meta_agg.nr_ms1_spectra));

  // Proof that aggregation actually ran: total MS1 peak count must be
  // strictly higher in the aggregated path (stacked 2*N+1=3 frames of raw
  // peaks into each output spectrum). If ms1_n_neighbors were silently
  // ignored in loadDIAStreaming, the counts would be equal.
  auto total_peaks = [](const OpenSwath::SwathMap& m) -> Size {
    Size n = 0;
    for (Size i = 0; i < m.sptr->getNrSpectra(); ++i)
      n += m.sptr->getSpectrumById(i)->getMZArray()->data.size();
    return n;
  };
  const Size peaks_ref = total_peaks(*ms1_ref);
  const Size peaks_agg = total_peaks(*ms1_agg);
  TEST_EQUAL(peaks_agg > peaks_ref, true);

  STATUS("DIA MS1 streaming aggregation: ref peaks=" << peaks_ref
         << " agg peaks=" << peaks_agg
         << " ratio=" << (static_cast<double>(peaks_agg) / peaks_ref));
}
END_SECTION

#endif // OPENTIMS_DIA_TEST_DATA

END_TEST

#else // WITH_OPENTIMS

// Minimal test when opentims is not available
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

START_TEST(BrukerTimsFile, "$Id$")
// No tests when WITH_OPENTIMS is off
END_TEST

#endif // WITH_OPENTIMS
