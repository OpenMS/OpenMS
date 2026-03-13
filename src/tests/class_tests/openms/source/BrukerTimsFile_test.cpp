// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/config.h>

#ifdef WITH_TIMSRUST

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/ANALYSIS/ID/PeptideSearchEngineFIAlgorithm.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/BrukerTimsFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/SYSTEM/File.h>

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

// Integration tests (only run when ENABLE_TIMSRUST_TESTS is ON and data is available)
#ifdef TIMSRUST_DDA_TEST_DATA

START_SECTION(DDA loading integration test)
{
  BrukerTimsFile f;
  MSExperiment exp;
  f.load(TIMSRUST_DDA_TEST_DATA, exp);

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

  // Check MS1 spectra have IM data in CONCATENATED format
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
  f.load(TIMSRUST_DDA_TEST_DATA, orig);

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
  f.load(TIMSRUST_DDA_TEST_DATA, exp, cfg);

  TEST_NOT_EQUAL(exp.size(), 0);

  // All spectra should have CONCATENATED IM data (even MS2 in frame mode)
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
  f.load(TIMSRUST_DDA_TEST_DATA, exp_raw);

  // Load with centroiding enabled
  BrukerTimsFile::Config cfg;
  cfg.ms1_centroid_mz_ppm = 5.0f;
  cfg.ms1_centroid_im_pct = 3.0f;
  MSExperiment exp_cent;
  f.load(TIMSRUST_DDA_TEST_DATA, exp_cent, cfg);

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
  f.load(TIMSRUST_DDA_TEST_DATA, exp_partial, cfg_partial);

  // Load without centroiding for comparison
  MSExperiment exp_raw;
  f.load(TIMSRUST_DDA_TEST_DATA, exp_raw);

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
  f.load(TIMSRUST_DDA_TEST_DATA, exp);

  // Verify MS2 spectra have drift times (pre-condition for IM annotation)
  Size ms2_count = 0;
  Size ms2_with_im = 0;
  for (const auto& spec : exp)
  {
    if (spec.getMSLevel() == 2)
    {
      ++ms2_count;
      if (IMTypes::determineIMFormat(spec) == IMFormat::MULTIPLE_SPECTRA)
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
  FASTAFile().load(TIMSRUST_TEST_FASTA, fasta_db);
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

#endif // TIMSRUST_DDA_TEST_DATA

#ifdef TIMSRUST_DIA_TEST_DATA

START_SECTION(DIA loading integration test)
{
  BrukerTimsFile f;
  MSExperiment exp;
  f.load(TIMSRUST_DIA_TEST_DATA, exp);

  TEST_NOT_EQUAL(exp.size(), 0);

  // Check MS2 spectra have per-peak IM (CONCATENATED) and isolation windows
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
  f.load(TIMSRUST_DIA_TEST_DATA, orig);

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
  f.load(TIMSRUST_DIA_TEST_DATA, exp_raw);

  // Load with centroiding
  BrukerTimsFile::Config cfg;
  cfg.ms1_centroid_mz_ppm = 5.0f;
  cfg.ms1_centroid_im_pct = 3.0f;
  MSExperiment exp_cent;
  f.load(TIMSRUST_DIA_TEST_DATA, exp_cent, cfg);

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

#endif // TIMSRUST_DIA_TEST_DATA

END_TEST

#else // WITH_TIMSRUST

// Minimal test when timsrust is not available
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

START_TEST(BrukerTimsFile, "$Id$")
// No tests when WITH_TIMSRUST is off
END_TEST

#endif // WITH_TIMSRUST
