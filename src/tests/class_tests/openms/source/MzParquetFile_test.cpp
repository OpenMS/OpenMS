// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/MzParquetFile.h>
///////////////////////////

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/METADATA/Precursor.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(MzParquetFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MzParquetFile* ptr = nullptr;

START_SECTION((MzParquetFile()))
  ptr = new MzParquetFile();
  TEST_EQUAL(ptr == nullptr, false)
END_SECTION

START_SECTION((~MzParquetFile()))
  delete ptr;
END_SECTION

START_SECTION((void load(const String& filename, PeakMap& map)))
  MzParquetFile file;
  PeakMap exp;

  // Test file not found exception
  TEST_EXCEPTION(Exception::FileNotFound, file.load("MzParquetFile_test_this_file_does_not_exist.parquet", exp))

  // Load the test parquet file
  file.load(OPENMS_GET_TEST_DATA_PATH("MzParquetFile_test.parquet"), exp);

  // Verify basic structure
  TEST_EQUAL(exp.size(), 4)  // 4 spectra as defined in minimal test

  // Check spectrum 1 (MS1, scan 1)
  MSSpectrum spec1 = exp[0];
  TEST_EQUAL(spec1.getMSLevel(), 1)
  TEST_EQUAL(spec1.size(), 3)  // 3 peaks
  TEST_REAL_SIMILAR(spec1.getRT(), 10.5)
  TEST_EQUAL(spec1.getNativeID(), "scan=1")
  
  // Check some peaks in spectrum 1
  TEST_REAL_SIMILAR(spec1[0].getMZ(), 445.2)
  TEST_REAL_SIMILAR(spec1[0].getIntensity(), 50000.0)
  TEST_REAL_SIMILAR(spec1[1].getMZ(), 523.8)
  TEST_REAL_SIMILAR(spec1[1].getIntensity(), 75000.0)
  TEST_REAL_SIMILAR(spec1[2].getMZ(), 612.3)
  TEST_REAL_SIMILAR(spec1[2].getIntensity(), 120000.0)

  // Check spectrum 2 (MS2, scan 2)
  MSSpectrum spec2 = exp[1];
  TEST_EQUAL(spec2.getMSLevel(), 2)
  TEST_EQUAL(spec2.size(), 3)  // 3 peaks
  TEST_REAL_SIMILAR(spec2.getRT(), 10.6)
  TEST_EQUAL(spec2.getNativeID(), "scan=2")
  
  // Check precursor information
  TEST_EQUAL(spec2.getPrecursors().size(), 1)
  Precursor prec2 = spec2.getPrecursors()[0];
  TEST_REAL_SIMILAR(prec2.getMZ(), 612.3)
  TEST_EQUAL(prec2.getCharge(), 2)
  TEST_REAL_SIMILAR(prec2.getActivationEnergy(), 30.0)
  
  // Check isolation window
  TEST_REAL_SIMILAR(prec2.getIsolationWindowLowerOffset(), 1.0)  // 612.3 - 611.3
  TEST_REAL_SIMILAR(prec2.getIsolationWindowUpperOffset(), 1.0)  // 613.3 - 612.3
  
  // Check ion mobility
  TEST_REAL_SIMILAR(spec2.getDriftTime(), 1.2)
  
  // Check some peaks in spectrum 2
  TEST_REAL_SIMILAR(spec2[0].getMZ(), 175.1)
  TEST_REAL_SIMILAR(spec2[0].getIntensity(), 25000.0)
  TEST_REAL_SIMILAR(spec2[1].getMZ(), 288.2)
  TEST_REAL_SIMILAR(spec2[1].getIntensity(), 45000.0)
  TEST_REAL_SIMILAR(spec2[2].getMZ(), 401.3)
  TEST_REAL_SIMILAR(spec2[2].getIntensity(), 35000.0)

  // Check spectrum 3 (MS1, scan 3)
  MSSpectrum spec3 = exp[2];
  TEST_EQUAL(spec3.getMSLevel(), 1)
  TEST_EQUAL(spec3.size(), 2)  // 2 peaks
  TEST_REAL_SIMILAR(spec3.getRT(), 15.2)
  TEST_EQUAL(spec3.getNativeID(), "scan=3")
  
  // Check peaks in spectrum 3
  TEST_REAL_SIMILAR(spec3[0].getMZ(), 389.7)
  TEST_REAL_SIMILAR(spec3[0].getIntensity(), 80000.0)
  TEST_REAL_SIMILAR(spec3[1].getMZ(), 456.8)
  TEST_REAL_SIMILAR(spec3[1].getIntensity(), 95000.0)

  // Check spectrum 4 (MS2, scan 4)
  MSSpectrum spec4 = exp[3];
  TEST_EQUAL(spec4.getMSLevel(), 2)
  TEST_EQUAL(spec4.size(), 3)  // 3 peaks
  TEST_REAL_SIMILAR(spec4.getRT(), 15.3)
  TEST_EQUAL(spec4.getNativeID(), "scan=4")
  
  // Check precursor information for spectrum 4
  TEST_EQUAL(spec4.getPrecursors().size(), 1)
  Precursor prec4 = spec4.getPrecursors()[0];
  TEST_REAL_SIMILAR(prec4.getMZ(), 456.8)
  TEST_EQUAL(prec4.getCharge(), 3)
  TEST_REAL_SIMILAR(prec4.getActivationEnergy(), 28.0)
  
  // Check isolation window for spectrum 4
  TEST_REAL_SIMILAR(prec4.getIsolationWindowLowerOffset(), 1.0)  // 456.8 - 455.8
  TEST_REAL_SIMILAR(prec4.getIsolationWindowUpperOffset(), 1.0)  // 457.8 - 456.8
  
  // Check that spectrum 4 has no ion mobility (should be -1 for missing value)
  TEST_EQUAL(spec4.getDriftTime(), -1.0)
  
  // Check some peaks in spectrum 4
  TEST_REAL_SIMILAR(spec4[0].getMZ(), 147.1)
  TEST_REAL_SIMILAR(spec4[0].getIntensity(), 15000.0)
  TEST_REAL_SIMILAR(spec4[1].getMZ(), 261.2)
  TEST_REAL_SIMILAR(spec4[1].getIntensity(), 28000.0)
  TEST_REAL_SIMILAR(spec4[2].getMZ(), 374.3)
  TEST_REAL_SIMILAR(spec4[2].getIntensity(), 22000.0)

  // Verify that spectra are sorted by RT
  TEST_EQUAL(exp[0].getRT() <= exp[1].getRT(), true)
  TEST_EQUAL(exp[1].getRT() <= exp[2].getRT(), true)
  TEST_EQUAL(exp[2].getRT() <= exp[3].getRT(), true)
  
  // Verify that peaks within each spectrum are sorted by m/z
  for (Size i = 0; i < exp.size(); ++i)
  {
    const MSSpectrum& spec = exp[i];
    for (Size j = 1; j < spec.size(); ++j)
    {
      TEST_EQUAL(spec[j-1].getMZ() <= spec[j].getMZ(), true)
    }
  }
END_SECTION

START_SECTION((void store(const String& filename, const PeakMap& map) const))
  MzParquetFile file;
  PeakMap exp;
  
  // Store is not implemented, should throw NotImplemented exception
  TEST_EXCEPTION(Exception::NotImplemented, file.store("test_output.parquet", exp))
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
