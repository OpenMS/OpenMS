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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
