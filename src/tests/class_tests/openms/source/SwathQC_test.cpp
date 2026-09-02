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
std::shared_ptr<MSExperiment> exp(new MSExperiment);
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_orbitrap_sn1_out.mzML"), *exp);
OpenSwath::SpectrumAccessPtr sptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);

std::vector< OpenSwath::SwathMap > swath_maps(1);
swath_maps.back().sptr = sptr;
swath_maps.back().ms1 = true;

START_SECTION((static ChargeDistribution getChargeDistribution(const std::vector<SwathMap>& swath_maps, const size_t nr_samples, const double mz_tol)))
{
  auto cd = SwathQC::getChargeDistribution(swath_maps, 10, 0.04);
  // expected values derive from PeakPickerHiRes_orbitrap_sn1_out.mzML, which is a
  // peak-picked reference (S/N threshold 1) and changes when the picker's noise
  // estimation changes
  SwathQC::ChargeDistribution cde = { {1,12}, {6,1}, {10,4} };
  TEST_EQUAL(cd.size(), cde.size());
  if (cd != cde)
  {
    std::cout << "Expected:\n";
    for (auto& c : cde)
    {
      std::cout << c.first << " " << c.second << "\n";
    }
    std::cout << "Got:\n";
    for (auto& c : cd)
    {
      std::cout << c.first << " " << c.second << "\n";
    }
  }
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

START_SECTION((static void storeJSON(const std::string& filename)))
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
  std::string tmp_json = File::getTemporaryFile();
  qc.storeJSON(tmp_json);
  std::string tmp_expected = File::getTemporaryFile();
  TextFile tf;
  tf.addLine(R"({
  "ChargeDistributionMS1": [
    [
      1,
      12
    ],
    [
      6,
      1
    ],
    [
      10,
      4
    ]
  ]
})");
  tf.store(tmp_expected);
  TEST_EQUAL(FuzzyStringComparator().compareFiles(tmp_json, tmp_expected), true);
}
END_SECTION

START_SECTION([EXTRA] getSpectraProcessingFunc uniform subsampling over ALL MS1 spectra (issue #9488, ANSW-5))
{
  // The test data has 3 MS1 spectra. Subsampling 2-of-3 must select the UNIFORMLY
  // spaced indices 0 and 2 (confirmed by isSubsampledSpectrum_), NOT the front-loaded
  // indices 0 and 1 that the pre-fix lambda produced (it only advanced its index for
  // ACCEPTED spectra, so spectrum 2 was checked with idx==1 and rejected).
  TEST_EQUAL(SwathQCTest::isSubsampledSpectrum_(3, 2, 0), true)
  TEST_EQUAL(SwathQCTest::isSubsampledSpectrum_(3, 2, 1), false)
  TEST_EQUAL(SwathQCTest::isSubsampledSpectrum_(3, 2, 2), true)

  // Member path: stream all 3 MS1 spectra through the processing function.
  SwathQC qc_member(2, 0.04);   // cd_spectra = 2
  qc_member.setNrMS1Spectra(3); // nr_ms1_spectra = 3
  auto f = qc_member.getSpectraProcessingFunc();
  for (auto& s : *exp)
  {
    if (s.getMSLevel() == 1) f(s);
  }
  const auto cd_member = qc_member.getChargeDistribution();

  // Static path: subsamples 2-of-3 using its own true loop index (always indices 0 and 2).
  const auto cd_static = SwathQC::getChargeDistribution(swath_maps, 2, 0.04);

  // Both paths must process the SAME two spectra (0 and 2) and yield identical distributions.
  // Pre-fix, the member path processed spectra 0 and 1 -> a different distribution.
  TEST_EQUAL(cd_member.size(), cd_static.size())
  TEST_TRUE(cd_member == cd_static)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
