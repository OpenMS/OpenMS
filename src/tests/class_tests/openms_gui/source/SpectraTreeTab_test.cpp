// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/VISUAL/SpectraTreeTab.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;
using namespace std;

namespace
{
  // Appends a spectrum of MS level @p ms_level. MSn spectra get a precursor, which refers to
  // the spectrum with native ID @p ref (if not empty), like the mzML attribute 'spectrumRef'.
  void addSpectrum(MSExperiment& exp, UInt ms_level, const std::string& native_id = "", const std::string& ref = "")
  {
    MSSpectrum spec;
    spec.setMSLevel(ms_level);
    spec.setNativeID(native_id);
    if (ms_level > 1)
    {
      Precursor precursor;
      if (!ref.empty())
      {
        precursor.setMetaValue("spectrum_ref", ref);
      }
      spec.getPrecursors().push_back(precursor);
    }
    exp.addSpectrum(spec);
  }

  // spectra of the given MS levels, without native IDs or precursor references
  MSExperiment levelsOnly(const std::vector<UInt>& ms_levels)
  {
    MSExperiment exp;
    for (UInt ms_level : ms_levels)
    {
      addSpectrum(exp, ms_level);
    }
    return exp;
  }
}

START_TEST(SpectraTreeTab, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static std::vector<int> getParentIndices(const MSExperiment& exp)))
{
  TEST_EQUAL(SpectraTreeTab::getParentIndices(MSExperiment()).empty(), true)

  { // no precursor references: the order of the spectra defines the hierarchy
    // (the second MS2 belongs to the MS1, not to the MS2 before the MS3 in between)
    std::vector<int> expected = {-1, 0, 1, 0, 3, 3, -1, 6};
    TEST_EQUAL(SpectraTreeTab::getParentIndices(levelsOnly({1, 2, 3, 2, 3, 3, 1, 2})), expected)
  }

  { // SPS-MS3: the MS3 scans are interleaved with later MS2 scans, but refer to their own MS2 scan
    MSExperiment exp;
    addSpectrum(exp, 1, "scan=1");
    addSpectrum(exp, 2, "scan=2", "scan=1");
    addSpectrum(exp, 2, "scan=3", "scan=1");
    addSpectrum(exp, 2, "scan=4", "scan=1");
    addSpectrum(exp, 3, "scan=5", "scan=3");
    addSpectrum(exp, 2, "scan=6", "scan=1");
    addSpectrum(exp, 3, "scan=7", "scan=4");
    addSpectrum(exp, 3, "scan=8", "scan=6");
    addSpectrum(exp, 3, "scan=9", "scan=2");
    addSpectrum(exp, 1, "scan=10");
    addSpectrum(exp, 2, "scan=11", "scan=10");
    addSpectrum(exp, 3, "scan=12", "scan=2"); // acquired after the next MS1 scan
    std::vector<int> expected = {-1, 0, 0, 0, 2, 0, 3, 5, 1, -1, 9, 1};
    TEST_EQUAL(SpectraTreeTab::getParentIndices(exp), expected)
  }

  { // references that cannot be used fall back to the order of the spectra
    MSExperiment exp;
    addSpectrum(exp, 1, "a");
    addSpectrum(exp, 2, "b", "a");
    addSpectrum(exp, 3, "c", "unknown"); // no such spectrum
    addSpectrum(exp, 3, "d", "a");       // refers to an MS1 spectrum
    addSpectrum(exp, 2, "e", "f");       // refers to a later spectrum
    addSpectrum(exp, 2, "f");            // no reference
    MSSpectrum no_precursor;
    no_precursor.setMSLevel(3);
    no_precursor.setNativeID("g");
    exp.addSpectrum(no_precursor);
    std::vector<int> expected = {-1, 0, 1, 1, 0, 0, 5};
    TEST_EQUAL(SpectraTreeTab::getParentIndices(exp), expected)
  }

  { // a reference wins over the order of the spectra, even across MS1 spectra;
    // if several spectra share a native ID, the last one before the referring spectrum counts
    MSExperiment exp;
    addSpectrum(exp, 1, "x");
    addSpectrum(exp, 2, "y", "x");
    addSpectrum(exp, 1, "w");
    addSpectrum(exp, 2, "z", "x");
    addSpectrum(exp, 1, "x");
    addSpectrum(exp, 2, "v", "x");
    std::vector<int> expected = {-1, 0, -1, 0, -1, 4};
    TEST_EQUAL(SpectraTreeTab::getParentIndices(exp), expected)
  }

  { // no spectrum one MS level lower since the last spectrum of an even lower level: top-level entry
    std::vector<int> expected = {-1, -1, 1, -1};
    TEST_EQUAL(SpectraTreeTab::getParentIndices(levelsOnly({2, 2, 3, 2})), expected)
    expected = {-1, -1, 0, 2};
    TEST_EQUAL(SpectraTreeTab::getParentIndices(levelsOnly({1, 3, 2, 3})), expected)
    expected = {-1, 0, -1, -1};
    TEST_EQUAL(SpectraTreeTab::getParentIndices(levelsOnly({1, 2, 1, 3})), expected)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
