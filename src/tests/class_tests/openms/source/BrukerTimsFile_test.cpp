// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/config.h>

#ifdef WITH_OPENTIMS

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/ANALYSIS/ID/PeptideSearchEngineFIAlgorithm.h>
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
      break;
    }
  }

  // Verify source file metadata was populated (I4)
  TEST_NOT_EQUAL(exp.getSourceFiles().size(), 0);
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

  // Run PeptideSearchEngineFIAlgorithm in-memory with the same FASTA and
  // parameters used by the TOPP-level DDA tests for SSE and FI.
  // The key test: any PSMs produced must have IM annotation.
  vector<FASTAFile::FASTAEntry> fasta_db;
  FASTAFile().load(OPENTIMS_TEST_FASTA, fasta_db);
  TEST_TRUE(fasta_db.size() > 0)

  // Typical timsTOF Pro DDA-PASEF search parameters
  PeptideSearchEngineFIAlgorithm algo;
  Param p = algo.getParameters();
  p.setValue("precursor:mass_tolerance", 20.0);
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

  TEST_EQUAL(ec == PeptideSearchEngineFIAlgorithm::ExitCodes::EXECUTION_OK, true)
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

  // Count MS2 spectra and total MS2 peaks in both
  Size raw_ms2_count = 0, agg_ms2_count = 0;
  Size raw_ms2_peaks = 0, agg_ms2_peaks = 0;
  for (const auto& spec : exp_raw)
  {
    if (spec.getMSLevel() == 2)
    {
      ++raw_ms2_count;
      raw_ms2_peaks += spec.size();
    }
  }
  for (const auto& spec : exp_agg)
  {
    if (spec.getMSLevel() == 2)
    {
      ++agg_ms2_count;
      agg_ms2_peaks += spec.size();
    }
  }

  // Both should have MS2 spectra
  TEST_NOT_EQUAL(raw_ms2_count, 0);
  TEST_NOT_EQUAL(agg_ms2_count, 0);

  // Same spectrum count (both use per-WindowGroup iteration)
  TEST_EQUAL(raw_ms2_count, agg_ms2_count);

  // Aggregation + denoising should reduce total peak count
  TEST_EQUAL(agg_ms2_peaks < raw_ms2_peaks, true);

  // Aggregated spectra should have per-peak IM data
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
         << " | aggregated MS2 spectra=" << agg_ms2_count
         << " peaks=" << agg_ms2_peaks);
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
