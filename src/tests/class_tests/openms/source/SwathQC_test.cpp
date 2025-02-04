// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
///////////////////////////

#include <OpenMS/ANALYSIS/OPENSWATH/SwathQC.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/SYSTEM/File.h>


using namespace OpenMS;
using namespace std;
using namespace OpenSwath;

class SwathQCTest : public SwathQC
{
  public:
    static bool isSubsampledSpectrum_(const size_t total_spec_count, const size_t subsample_count, const size_t idx)
    {
      return SwathQC::isSubsampledSpectrum_(total_spec_count, subsample_count, idx);
    }
};

START_TEST(SwathQC, "$Id$")
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SwathQC* nullPointer = nullptr;
SwathQC* ptr = nullptr;

START_SECTION(SwathQC())
{
  ptr = new SwathQC(10, 0.04);
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~SwathQC())
{
  delete ptr;
}
END_SECTION

// Create a mock spectrum fitting to the transition group
boost::shared_ptr<MSExperiment> exp(new MSExperiment);
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_orbitrap_sn1_out.mzML"), *exp);
OpenSwath::SpectrumAccessPtr sptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);

std::vector< OpenSwath::SwathMap > swath_maps(1);
swath_maps.back().sptr = sptr;
swath_maps.back().ms1 = true;

START_SECTION((static ChargeDistribution getChargeDistribution(const std::vector<SwathMap>& swath_maps, const size_t nr_samples, const double mz_tol)))
{
  auto cd = SwathQC::getChargeDistribution(swath_maps, 10, 0.04);
  SwathQC::ChargeDistribution cde = { {1,17}, {2,4}, {5,1}, {6,2}, {8,2}, {9,1}, {10,5} };
  TEST_EQUAL(cd.size(), cde.size());
  TEST_TRUE(cd == cde)
}
END_SECTION

START_SECTION((static bool isSubsampledSpectrum_(const size_t total_spec_count, const size_t subsample_count, const size_t idx)))
{
  TEST_EQUAL(SwathQCTest::isSubsampledSpectrum_(0, 100, 4), true); // always true (unknown number of total spectra)
  TEST_EQUAL(SwathQCTest::isSubsampledSpectrum_(10, 100, 4), true); // always true (not enough samples)
  TEST_EQUAL(SwathQCTest::isSubsampledSpectrum_(10, 4, 10), false); // always false (index beyond # of total spectra)
  TEST_EQUAL(SwathQCTest::isSubsampledSpectrum_(10, 4, 11), false); // always false (index beyond # of total spectra)

  int r[] = {1, 0, 0, 1, 0, 1, 0, 0, 1, 0};
  int c = 10;
  for (int i = 0; i < c; ++i)
  {
    //std::cout << i << ": " << SwathQCTest::isSubsampledSpectrum_(c, 4, i) << "\n";
    TEST_EQUAL(SwathQCTest::isSubsampledSpectrum_(c, 4, i), r[i]);
  }

  // sample none
  c = 10;
  for (int i = 0; i < c; ++i)
  {
    //std::cout << i << ": " << SwathQCTest::isSubsampledSpectrum_(c, 0, i) << "\n";
    TEST_EQUAL(SwathQCTest::isSubsampledSpectrum_(c, 0, i), false);
  }

  // sample all
  c = 4;
  for (int i = 0; i < c; ++i)
  {
    //std::cout << i << ": " << SwathQCTest::isSubsampledSpectrum_(c, c, i) << "\n";
    TEST_EQUAL(SwathQCTest::isSubsampledSpectrum_(c, c, i), true);
  }

  // sample 2 of 5
  c = 5;
  int r5[] = {1,0,0,1,0};
  for (int i = 0; i < c; ++i)
  {
    //std::cout << i << ": " << SwathQCTest::isSubsampledSpectrum_(5, 2, i) << "\n";
    TEST_EQUAL(SwathQCTest::isSubsampledSpectrum_(c, 2, i), r5[i]);
  }

}
END_SECTION

START_SECTION((static void storeJSON(const OpenMS::String& filename)))
{
  SwathQC qc(10, 0.04);
  int count{};
  for (auto& s : *exp)
  {
    if (s.getMSLevel()==1) ++count;
  }
  qc.setNrMS1Spectra(count);
  auto f = qc.getSpectraProcessingFunc();
  for (auto& s : *exp)
  {
    if (s.getMSLevel()==1) f(s);
  }

  // getChargeDistribution(swath_maps, 10, 0.04);
  String tmp_json = File::getTemporaryFile();
  qc.storeJSON(tmp_json);
  String tmp_expected = File::getTemporaryFile();
  TextFile tf;
  tf.addLine(R"({
  "ChargeDistributionMS1": [
    [
      1,
      17
    ],
    [
      2,
      4
    ],
    [
      5,
      1
    ],
    [
      6,
      2
    ],
    [
      8,
      2
    ],
    [
      9,
      1
    ],
    [
      10,
      5
    ]
  ]
})");
  tf.store(tmp_expected);
  TEST_EQUAL(FuzzyStringComparator().compareFiles(tmp_json, tmp_expected), true);

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
