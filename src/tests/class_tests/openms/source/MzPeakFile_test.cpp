// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Kohlbacher $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzPeakFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/METADATA/SpectrumSettings.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(MzPeakFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MzPeakFile* ptr = nullptr;
MzPeakFile* nullPointer = nullptr;
START_SECTION((MzPeakFile()))
ptr = new MzPeakFile;
TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~MzPeakFile()))
delete ptr;
END_SECTION

START_SECTION(void load(const String& filename, MapType& map))
{
  // ----------------------------------------------------------------------
  // Missing file must raise FileNotFound.
  // ----------------------------------------------------------------------
  {
    MSExperiment exp;
    TEST_EXCEPTION(Exception::FileNotFound, MzPeakFile().load("does_not_exist.mzpeak", exp))
  }

  // ----------------------------------------------------------------------
  // Load the bundled point-layout fixture (14 profile MS1 + 34 centroid MS2).
  // ----------------------------------------------------------------------
  MSExperiment exp;
  MzPeakFile().load(OPENMS_GET_TEST_DATA_PATH("small.mzpeak"), exp);

  // 48 spectra total.
  TEST_EQUAL(exp.size(), 48)

  // ms_level split: 14 MS1 + 34 MS2.
  Size n_ms1 = 0;
  Size n_ms2 = 0;
  for (Size i = 0; i < exp.size(); ++i)
  {
    if (exp[i].getMSLevel() == 1) ++n_ms1;
    else if (exp[i].getMSLevel() == 2)
      ++n_ms2;
  }
  TEST_EQUAL(n_ms1, 14)
  TEST_EQUAL(n_ms2, 34)

  // updateRanges() must have populated sane, non-default ranges.
  TEST_EQUAL(exp.getMinMZ() < exp.getMaxMZ(), true)
  TEST_EQUAL(exp.getMaxMZ() > 1000.0, true)

  // ----------------------------------------------------------------------
  // Profile spectrum (mzpeak index 0): native id "index=0", MS1, RT 0.004935 s.
  // After null-mz reconstruction it has 13589 points, strictly ascending m/z,
  // NO interior zero m/z; first/last m/z ~ 202.6066 / 1999.8404.
  // ----------------------------------------------------------------------
  const MSSpectrum* prof = nullptr;
  for (Size i = 0; i < exp.size(); ++i)
  {
    if (exp[i].getNativeID() == "index=0")
    {
      prof = &exp[i];
      break;
    }
  }
  TEST_NOT_EQUAL(prof, nullptr)
  if (prof)
  {
    TEST_EQUAL(prof->getMSLevel(), 1)
    TEST_REAL_SIMILAR(prof->getRT(), 0.004935)
    TEST_EQUAL(prof->size(), 13589)
    TEST_EQUAL(prof->empty(), false)

    // Strictly ascending m/z and no interior zero m/z (reconstruction worked).
    bool ascending = true;
    bool interior_zero = false;
    for (Size k = 1; k < prof->size(); ++k)
    {
      if ((*prof)[k].getMZ() <= (*prof)[k - 1].getMZ()) ascending = false;
      if ((*prof)[k].getMZ() == 0.0) interior_zero = true;
    }
    if ((*prof)[0].getMZ() == 0.0) interior_zero = true;
    TEST_EQUAL(ascending, true)
    TEST_EQUAL(interior_zero, false)

    TEST_REAL_SIMILAR((*prof)[0].getMZ(), 202.60657495520474)
    TEST_REAL_SIMILAR((*prof)[prof->size() - 1].getMZ(), 1999.8404377599534)

    // Checked intensity: point 1 has intensity 1938.1174 (point 0 intensity 0).
    TOLERANCE_ABSOLUTE(1e-2)
    TEST_REAL_SIMILAR((*prof)[1].getIntensity(), 1938.117431640625)
    TOLERANCE_ABSOLUTE(1e-4)
  }

  // ----------------------------------------------------------------------
  // Centroid spectrum (mzpeak index 2): native id "index=2", MS2, RT 0.011218 s.
  // 485 peaks; first peak (231.38884, 26.5451), last (1560.71985, 22.9731).
  // ----------------------------------------------------------------------
  const MSSpectrum* cent = nullptr;
  for (Size i = 0; i < exp.size(); ++i)
  {
    if (exp[i].getNativeID() == "index=2")
    {
      cent = &exp[i];
      break;
    }
  }
  TEST_NOT_EQUAL(cent, nullptr)
  if (cent)
  {
    TEST_EQUAL(cent->getMSLevel(), 2)
    TEST_REAL_SIMILAR(cent->getRT(), 0.011218333333)
    TEST_EQUAL(cent->size(), 485)
    TEST_EQUAL(cent->empty(), false)

    // Peaks are sorted by m/z.
    bool sorted = true;
    for (Size k = 1; k < cent->size(); ++k)
    {
      if ((*cent)[k].getMZ() < (*cent)[k - 1].getMZ()) sorted = false;
    }
    TEST_EQUAL(sorted, true)

    TEST_REAL_SIMILAR((*cent)[0].getMZ(), 231.3888397216797)
    TEST_REAL_SIMILAR((*cent)[cent->size() - 1].getMZ(), 1560.7198486328125)

    TOLERANCE_ABSOLUTE(1e-3)
    TEST_REAL_SIMILAR((*cent)[0].getIntensity(), 26.54511260986328)
    TEST_REAL_SIMILAR((*cent)[cent->size() - 1].getIntensity(), 22.973094940185547)
    TOLERANCE_ABSOLUTE(1e-4)
  }

  // ======================================================================
  // METADATA layer (02-02).
  // ======================================================================

  // ----------------------------------------------------------------------
  // Spectrum type: profile (index 0) is PROFILE, centroid (index 2) is CENTROID.
  // Polarity is positive on both; native id is non-empty.
  // ----------------------------------------------------------------------
  if (prof)
  {
    TEST_EQUAL(prof->getType() == SpectrumSettings::SpectrumType::PROFILE, true)
    TEST_EQUAL(prof->getInstrumentSettings().getPolarity() == IonSource::Polarity::POSITIVE, true)
    TEST_EQUAL(prof->getNativeID().empty(), false)
  }
  if (cent)
  {
    TEST_EQUAL(cent->getType() == SpectrumSettings::SpectrumType::CENTROID, true)
    TEST_EQUAL(cent->getInstrumentSettings().getPolarity() == IonSource::Polarity::POSITIVE, true)
    TEST_EQUAL(cent->getNativeID().empty(), false)

    // The mzPeak native id (controllerType=...) is preserved as a meta value.
    TEST_EQUAL(cent->metaValueExists("mzpeak_native_id"), true)

    // --------------------------------------------------------------------
    // Precursor (index 2): isolation target ~810.789, offsets 1.0/1.0;
    // activation CID + collision energy 35; selected-ion m/z ~810.789428,
    // intensity ~1994039.125.
    // --------------------------------------------------------------------
    TEST_EQUAL(cent->getPrecursors().size() >= 1, true)
    if (! cent->getPrecursors().empty())
    {
      const Precursor& prec = cent->getPrecursors()[0];

      // m/z follows the selected ion (~810.789428); within isolation tolerance.
      TOLERANCE_ABSOLUTE(1e-2)
      TEST_REAL_SIMILAR(prec.getMZ(), 810.789428)
      TEST_REAL_SIMILAR(prec.getIsolationWindowLowerOffset(), 1.0)
      TEST_REAL_SIMILAR(prec.getIsolationWindowUpperOffset(), 1.0)

      // Activation: CID method present, collision energy 35.
      TEST_EQUAL(prec.getActivationMethods().count(Precursor::ActivationMethod::CID) == 1, true)
      TEST_REAL_SIMILAR(prec.getActivationEnergy(), 35.0)

      // Selected-ion intensity ~1994039.125.
      TOLERANCE_ABSOLUTE(1.0)
      TEST_REAL_SIMILAR(prec.getIntensity(), 1994039.125)
      TOLERANCE_ABSOLUTE(1e-4)

      // Isolation target m/z kept as a meta value (selected ion overrode m/z).
      TEST_EQUAL(prec.metaValueExists("isolation window target m/z"), true)
    }
  }

  // ----------------------------------------------------------------------
  // Run-level metadata (ExperimentalSettings): instrument, source file, sample.
  // small.mzpeak carries an instrument_configuration_list (LTQ FT), a
  // source_files block (small.RAW), and a sample_list.
  // ----------------------------------------------------------------------
  TEST_EQUAL(exp.getInstrument().getName().empty(), false)
  TEST_EQUAL(exp.getInstrument().getIonSources().empty(), false)
  TEST_EQUAL(exp.getSourceFiles().empty(), false)
  if (! exp.getSourceFiles().empty())
  {
    TEST_EQUAL(exp.getSourceFiles()[0].getNameOfFile(), "small.RAW")
    // SHA-1 checksum from MS:1000569 mapped via the source-file CV dispatch.
    TEST_EQUAL(exp.getSourceFiles()[0].getChecksum().empty(), false)
  }
  // An unrecognized CvParam surfaces as a meta value: the instrument serial
  // number (MS:1000529) is kept on the Instrument as a meta value.
  TEST_EQUAL(exp.getInstrument().metaValueExists("instrument serial number"), true)

  // Run date parsed from the ISO start_time (2005-07-20T19:44:22Z).
  TEST_EQUAL(exp.getDateTime().isNull(), false)
}
END_SECTION

START_SECTION([EXTRA] load via generic FileHandler API)
{
  // .mzpeak must resolve by extension and by content, and route to MzPeakFile.
  TEST_EQUAL(FileHandler::getTypeByFileName("x.mzpeak") == FileTypes::MZPEAK, true)
  TEST_EQUAL(FileHandler::getTypeByContent(OPENMS_GET_TEST_DATA_PATH("small.mzpeak")) == FileTypes::MZPEAK, true)

  FileHandler fh;
  MSExperiment exp;
  fh.loadExperiment(OPENMS_GET_TEST_DATA_PATH("small.mzpeak"), exp);
  TEST_EQUAL(exp.size(), 48)
}
END_SECTION

START_SECTION(void store(const String& filename, const MapType& map))
{
  // ----------------------------------------------------------------------
  // Round trip: load the bundled fixture, store it to a fresh .mzpeak, reload,
  // and assert the reloaded experiment is EQUIVALENT to the source (count,
  // ms_level, type, RT, peak count, and first/last m/z + intensity per
  // spectrum within tolerance). Spectra are sorted ascending by m/z.
  // ----------------------------------------------------------------------
  MSExperiment src;
  MzPeakFile().load(OPENMS_GET_TEST_DATA_PATH("small.mzpeak"), src);
  TEST_EQUAL(src.size(), 48)

  String tmp;
  NEW_TMP_FILE_EXT(tmp, ".mzpeak")
  MzPeakFile().store(tmp, src);

  MSExperiment rt;
  MzPeakFile().load(tmp, rt);

  TEST_EQUAL(rt.size(), src.size())

  // load() sorts spectra by RT; map source spectra by native id so the compare
  // is robust to ordering differences between the two loads.
  for (Size i = 0; i < src.size(); ++i)
  {
    const MSSpectrum& s = src[i];
    const MSSpectrum* r = nullptr;
    for (Size j = 0; j < rt.size(); ++j)
    {
      if (rt[j].getNativeID() == s.getNativeID())
      {
        r = &rt[j];
        break;
      }
    }
    TEST_NOT_EQUAL(r, nullptr)
    if (! r) continue;

    TEST_EQUAL(r->getMSLevel(), s.getMSLevel())
    TEST_EQUAL(r->getType() == s.getType(), true)
    TOLERANCE_ABSOLUTE(1e-6)
    TEST_REAL_SIMILAR(r->getRT(), s.getRT())
    TEST_EQUAL(r->size(), s.size())

    if (! s.empty() && r->size() == s.size())
    {
      // First and last peak m/z (tight absolute) and intensity (relative).
      TOLERANCE_ABSOLUTE(1e-6)
      TEST_REAL_SIMILAR((*r)[0].getMZ(), s[0].getMZ())
      TEST_REAL_SIMILAR((*r)[r->size() - 1].getMZ(), s[s.size() - 1].getMZ())

      TOLERANCE_RELATIVE(1.0 + 1e-3)
      TEST_REAL_SIMILAR((*r)[0].getIntensity(), s[0].getIntensity())
      TEST_REAL_SIMILAR((*r)[r->size() - 1].getIntensity(), s[s.size() - 1].getIntensity())
      TOLERANCE_RELATIVE(1.0)
      TOLERANCE_ABSOLUTE(1e-4)
    }
  }

  // Counts of MS1/MS2 must be preserved exactly.
  Size src_ms1 = 0, rt_ms1 = 0;
  for (Size i = 0; i < src.size(); ++i)
    if (src[i].getMSLevel() == 1) ++src_ms1;
  for (Size i = 0; i < rt.size(); ++i)
    if (rt[i].getMSLevel() == 1) ++rt_ms1;
  TEST_EQUAL(rt_ms1, src_ms1)

  // --- Precursor round-trip: MS2 precursor m/z, isolation, activation, selected ion ---
  {
    // Find the first MS2 spectrum in src that has a precursor
    const MSSpectrum* src_ms2 = nullptr;
    for (Size i = 0; i < src.size(); ++i)
      if (src[i].getMSLevel() == 2 && ! src[i].getPrecursors().empty())
      {
        src_ms2 = &src[i];
        break;
      }
    TEST_NOT_EQUAL(src_ms2, nullptr)
    if (src_ms2)
    {
      // Find matching spectrum in rt
      const MSSpectrum* rt_ms2 = nullptr;
      for (Size i = 0; i < rt.size(); ++i)
        if (rt[i].getNativeID() == src_ms2->getNativeID())
        {
          rt_ms2 = &rt[i];
          break;
        }
      TEST_NOT_EQUAL(rt_ms2, nullptr)
      if (rt_ms2 && ! src_ms2->getPrecursors().empty() && ! rt_ms2->getPrecursors().empty())
      {
        const Precursor& sp = src_ms2->getPrecursors()[0];
        const Precursor& rp = rt_ms2->getPrecursors()[0];
        TOLERANCE_ABSOLUTE(1e-3)
        TEST_REAL_SIMILAR(rp.getMZ(), sp.getMZ())
        TEST_REAL_SIMILAR(rp.getIsolationWindowLowerOffset(), sp.getIsolationWindowLowerOffset())
        TEST_REAL_SIMILAR(rp.getIsolationWindowUpperOffset(), sp.getIsolationWindowUpperOffset())
        // Activation: CID method preserved
        TEST_EQUAL(rp.getActivationMethods().count(Precursor::ActivationMethod::CID)
                     == sp.getActivationMethods().count(Precursor::ActivationMethod::CID),
                   true)
        if (! sp.getActivationMethods().empty()) TEST_REAL_SIMILAR(rp.getActivationEnergy(), sp.getActivationEnergy())
        TOLERANCE_ABSOLUTE(1.0)
        TEST_REAL_SIMILAR(rp.getIntensity(), sp.getIntensity())
        TOLERANCE_ABSOLUTE(1e-4)
      }
    }
  }

  // --- Run-level metadata round-trip ---
  {
    // Re-load into fresh experiment to avoid side effects
    MSExperiment rt2;
    MzPeakFile().load(tmp, rt2);
    TEST_EQUAL(rt2.getInstrument().getName().empty(), false)
    TEST_EQUAL(rt2.getSourceFiles().empty(), false)
    if (! rt2.getSourceFiles().empty()) TEST_EQUAL(rt2.getSourceFiles()[0].getNameOfFile().empty(), false)
  }
}
END_SECTION

START_SECTION([EXTRA] store via generic FileHandler API)
{
  // FileHandler().storeExperiment must route .mzpeak to MzPeakFile::store.
  MSExperiment src;
  MzPeakFile().load(OPENMS_GET_TEST_DATA_PATH("small.mzpeak"), src);

  String tmp;
  NEW_TMP_FILE_EXT(tmp, ".mzpeak")
  FileHandler fh;
  fh.storeExperiment(tmp, src);

  MSExperiment rt;
  fh.loadExperiment(tmp, rt);
  TEST_EQUAL(rt.size(), src.size())
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
