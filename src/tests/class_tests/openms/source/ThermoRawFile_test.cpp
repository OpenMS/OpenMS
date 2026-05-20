// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/config.h>

#ifdef WITH_THERMO_RAW

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/ThermoRawFile.h>
#include <OpenMS/METADATA/IonSource.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/SYSTEM/File.h>

#include <set>

using namespace OpenMS;
using namespace std;

START_TEST(ThermoRawFile, "$Id$")

START_SECTION(void load(const String& path, MSExperiment& exp))
{
  ThermoRawFile file;
  MSExperiment exp;

  TEST_EXCEPTION(Exception::FileNotFound, file.load("/definitely/not/present.raw", exp))
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

  String tmp_mzml = File::getTempDirectory() + "/" + File::getUniqueName() + "_thermo_roundtrip.mzML";
  MzMLFile().store(tmp_mzml, original);

  MSExperiment reloaded;
  MzMLFile().load(tmp_mzml, reloaded);
  File::remove(tmp_mzml);

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
