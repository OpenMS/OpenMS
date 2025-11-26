// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz, Chris Bielow $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(IMDataConverter, "$Id$")

/////////////////////////////////////////////////////////////

IMDataConverter* e_ptr = nullptr;
IMDataConverter* e_nullPointer = nullptr;

START_SECTION((IMDataConverter()))
	e_ptr = new IMDataConverter;
  TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~IMDataConverter()))
	delete e_ptr;
END_SECTION


START_SECTION((std::vector<PeakMap> splitByFAIMSCV(PeakMap& exp)))
  MzMLFile IM_file;
  PeakMap exp;
  IM_file.load(OPENMS_GET_TEST_DATA_PATH("IM_FAIMS_test.mzML"), exp);

  TEST_EQUAL(exp.getSpectra().size(), 19)

  vector<PeakMap> splitPeakMap = IMDataConverter::splitByFAIMSCV(std::move(exp));
  TEST_EQUAL(exp.empty(), true) // moved out
  TEST_EQUAL(splitPeakMap.size(), 3)

	TEST_EQUAL(splitPeakMap[0].size(), 4)
	TEST_EQUAL(splitPeakMap[1].size(), 9)
	TEST_EQUAL(splitPeakMap[2].size(), 6)

	for (PeakMap::Iterator it = splitPeakMap[0].begin(); it != splitPeakMap[0].end(); ++it)
	{
		TEST_EQUAL(it->getDriftTime(), -65.0)
	}
	for (PeakMap::Iterator it = splitPeakMap[1].begin(); it != splitPeakMap[1].end(); ++it)
	{
		TEST_EQUAL(it->getDriftTime(), -55.0)
	}
	for (PeakMap::Iterator it = splitPeakMap[2].begin(); it != splitPeakMap[2].end(); ++it)
	{
		TEST_EQUAL(it->getDriftTime(), -45.0)
	}

	TEST_EQUAL(splitPeakMap[1].getExperimentalSettings().getDateTime().toString(), "2019-09-07T09:40:04")

END_SECTION


START_SECTION(static void setIMUnit(DataArrays::FloatDataArray& fda, const DriftTimeUnit unit))
	MSSpectrum::FloatDataArray fda;
  TEST_EXCEPTION(Exception::InvalidValue, IMDataConverter::setIMUnit(fda, DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE))
  TEST_EXCEPTION(Exception::InvalidValue, IMDataConverter::setIMUnit(fda, DriftTimeUnit::NONE))
  
	DriftTimeUnit unit;
  IMDataConverter::setIMUnit(fda, DriftTimeUnit::MILLISECOND);
  TEST_EQUAL(IMDataConverter::getIMUnit(fda, unit), true)
  TEST_EQUAL(DriftTimeUnit::MILLISECOND == unit, true)

	IMDataConverter::setIMUnit(fda, DriftTimeUnit::VSSC);
  TEST_EQUAL(IMDataConverter::getIMUnit(fda, unit), true)
  TEST_EQUAL(DriftTimeUnit::VSSC == unit, true)
END_SECTION

START_SECTION(static bool getIMUnit(const DataArrays::FloatDataArray& fda, DriftTimeUnit& unit))
	NOT_TESTABLE // tested above
END_SECTION


MSSpectrum frame;
frame.push_back({1.0, 11.0f});
frame.push_back({1.0 + Math::ppmToMass(4.0, 1.0), 12.0f}); // should merge with the one above
frame.push_back({1.2, 13.0f});
frame.push_back({2.0, 20.0f});
frame.push_back({3.0 - Math::ppmToMass(3.0, 3.0), 32.0f}); // should merge with the one below
frame.push_back({3.0, 31.0f});

frame.push_back({4.0, 40.0f});

frame.push_back({5.0, 50.0f});
frame.push_back({6.0, 60.0f});
frame.push_back({7.0, 70.0f});
frame.setRT(1);
MSSpectrum::FloatDataArray& afa = frame.getFloatDataArrays().emplace_back();
//           <---------- bin 1 ----------->   < bin 2 >  < -- bin 3 --->
afa.assign({1.1, 1.11, 1.11, 2.2, 3.2, 3.22,    4.4,      5.6, 6.6, 7.7});
IMDataConverter::setIMUnit(afa, DriftTimeUnit::MILLISECOND);

MSSpectrum spec;
spec.push_back({111.0, -1.0f});
spec.push_back({222.0, -2.0f});
spec.push_back({333.0, -3.0f});
spec.setRT(2); // just a spectrum with RT = 2


START_SECTION(static MSExperiment reshapeIMFrameToMany(MSSpectrum im_frame))
{
  // not am IM frame:
  TEST_EXCEPTION(Exception::MissingInformation, IMDataConverter::reshapeIMFrameToMany(spec))
	
  {
		auto exp = IMDataConverter::reshapeIMFrameToMany(frame);
		TEST_EQUAL(exp.size(), 9); // nine different IM-values
		TEST_EQUAL(exp[0].size(), 1);
    TEST_EQUAL(exp[1].size(), 2);
    TEST_EQUAL(exp[1][0].getIntensity(), 12.0f);
    TEST_EQUAL(exp[1][1].getIntensity(), 13.0f);

		TEST_EQUAL(exp[0].getDriftTime(), 1.1f);
		TEST_TRUE(exp[0].getDriftTimeUnit() == DriftTimeUnit::MILLISECOND);
		TEST_EQUAL(exp[0].getRT(), 1);
    TEST_EQUAL(exp[1].getDriftTime(), 1.11f);
    TEST_EQUAL(exp[8].getDriftTime(), 7.7f);
    TEST_TRUE(exp.isIMFrame());

		auto frame_reconstruct = IMDataConverter::reshapeIMFrameToSingle(exp);
		TEST_EQUAL(frame_reconstruct.size(), 1)
		TEST_EQUAL(frame_reconstruct[0], frame);
	}
}
END_SECTION

START_SECTION((static std::tuple<std::vector<MSExperiment>, Math::BinContainer> splitExperimentByIonMobility(MSExperiment&& in, UInt number_of_IM_bins, double bin_extension_abs, double mz_binning_width, MZ_UNITS mz_binning_width_unit)))
{
	MSExperiment e_in;
	e_in.addSpectrum(frame);
  auto frame2 = frame; // a second frame so we can test if two RT's show up in the result
  frame2.setRT(3);
  e_in.addSpectrum(frame2);
  
  // IM-range is 7.7-1.1 = 6.6
  // --> each bin is 2.2 wide
	const auto [exp_slices, bin_values] = IMDataConverter::splitExperimentByIonMobility(std::move(e_in), 3, 0.0, 5.0, MZ_UNITS::PPM);
  const auto ranges = Math::BinContainer { {1.1, 3.3},  {3.3, 5.5},  {5.5, 7.7}};
  for (int i = 0; i < 3; ++i)
  {
    TEST_REAL_SIMILAR(bin_values[i].getMin(), ranges[i].getMin());
    TEST_REAL_SIMILAR(bin_values[i].getMax(), ranges[i].getMax());
  }
  TEST_EQUAL(exp_slices.size(), 3);
  const auto& exp11 = exp_slices[0];
  const auto& exp33 = exp_slices[1];
  const auto& exp55 = exp_slices[2];
  TEST_EQUAL(exp11[0].size(), 4);
  TEST_EQUAL(exp33[0].size(), 1);
  TEST_EQUAL(exp55[0].size(), 3);
  TEST_EQUAL(exp11[1].size(), 4); // second frame. Identical to first frame
  TEST_EQUAL(exp33[1].size(), 1);
  TEST_EQUAL(exp55[1].size(), 3);
  TEST_EQUAL(exp11[0][0].getIntensity(), 11+12);
  TEST_REAL_SIMILAR(exp11[0].getDriftTime(), 2.2f); // center of bin 1.1-3.3
  TEST_TRUE(exp11[0].getDriftTimeUnit() == DriftTimeUnit::MILLISECOND);
  TEST_EQUAL(exp11[0].getRT(), 1);
  TEST_EQUAL(exp11[1].getRT(), 3);
  TEST_EQUAL(exp33[0].getRT(), 1);
  TEST_EQUAL(exp33[1].getRT(), 3);
  TEST_EQUAL(exp55[0].getRT(), 1);
  TEST_EQUAL(exp55[1].getRT(), 3);
}
END_SECTION

START_SECTION(static MSExperiment reshapeIMFrameToSingle(const MSExperiment& in))
	NOT_TESTABLE // tested_above
END_SECTION

START_SECTION((std::vector<PeakMap> splitByFAIMSCV assigns MS2 without explicit CV to last seen FAIMS CV))
{
  PeakMap exp_synth;

  MSSpectrum ms1a;
  ms1a.setMSLevel(1);
  ms1a.setDriftTimeUnit(DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
  ms1a.setDriftTime(-55.0);

  MSSpectrum ms2a;
  ms2a.setMSLevel(2);
  ms2a.setDriftTimeUnit(DriftTimeUnit::NONE);

  MSSpectrum ms2b;
  ms2b.setMSLevel(2);
  ms2b.setDriftTimeUnit(DriftTimeUnit::NONE);

  MSSpectrum ms1b;
  ms1b.setMSLevel(1);
  ms1b.setDriftTimeUnit(DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
  ms1b.setDriftTime(-45.0);

  MSSpectrum ms2c;
  ms2c.setMSLevel(2);
  ms2c.setDriftTimeUnit(DriftTimeUnit::NONE);

  exp_synth.addSpectrum(ms1a);
  exp_synth.addSpectrum(ms2a);
  exp_synth.addSpectrum(ms2b);
  exp_synth.addSpectrum(ms1b);
  exp_synth.addSpectrum(ms2c);

  std::vector<PeakMap> bins = IMDataConverter::splitByFAIMSCV(std::move(exp_synth));
  TEST_EQUAL(bins.size(), 2)

  // bins ordered by ascending CV: -55 first, -45 second
  TEST_EQUAL(bins[0].size(), 3) // ms1a + ms2a + ms2b
  TEST_EQUAL(bins[1].size(), 2) // ms1b + ms2c

  Size ms2_bin0 = 0;
  for (const auto& s : bins[0]) if (s.getMSLevel() > 1) ++ms2_bin0;
  TEST_EQUAL(ms2_bin0, 2)

  Size ms2_bin1 = 0;
  for (const auto& s : bins[1]) if (s.getMSLevel() > 1) ++ms2_bin1;
  TEST_EQUAL(ms2_bin1, 1)
}
END_SECTION

START_SECTION((std::vector<PeakMap> splitByFAIMSCV throws for non-FAIMS dataset))
{
  PeakMap exp_nonfaims;
  MSSpectrum s;
  s.setMSLevel(1);
  s.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
  s.setDriftTime(10.0);
  exp_nonfaims.addSpectrum(s);

  TEST_EXCEPTION(Exception::MissingInformation, IMDataConverter::splitByFAIMSCV(std::move(exp_nonfaims)))
}
END_SECTION

START_SECTION((std::vector<PeakMap> splitByFAIMSCV with include_non_faims_cv0=true includes non-FAIMS spectra))
{
  PeakMap exp_mixed;

  // Non-FAIMS spectrum
  MSSpectrum non_faims;
  non_faims.setMSLevel(1);
  non_faims.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
  non_faims.setDriftTime(10.0);

  // FAIMS spectrum with CV=-55
  MSSpectrum faims1;
  faims1.setMSLevel(1);
  faims1.setDriftTimeUnit(DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
  faims1.setDriftTime(-55.0);

  // FAIMS spectrum with CV=-45
  MSSpectrum faims2;
  faims2.setMSLevel(1);
  faims2.setDriftTimeUnit(DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
  faims2.setDriftTime(-45.0);

  exp_mixed.addSpectrum(non_faims);
  exp_mixed.addSpectrum(faims1);
  exp_mixed.addSpectrum(faims2);

  std::vector<PeakMap> bins = IMDataConverter::splitByFAIMSCV(std::move(exp_mixed), true);
  // Expect 3 bins: CV=0 for non-FAIMS, CV=-55, CV=-45
  TEST_EQUAL(bins.size(), 3)

  // First bin should be the non-FAIMS spectra (CV=0)
  TEST_EQUAL(bins[0].size(), 1)
  TEST_EQUAL(bins[0][0].getDriftTimeUnit(), DriftTimeUnit::MILLISECOND)
  
  // Second bin should be CV=-55
  TEST_EQUAL(bins[1].size(), 1)
  TEST_EQUAL(bins[1][0].getDriftTime(), -55.0)
  
  // Third bin should be CV=-45
  TEST_EQUAL(bins[2].size(), 1)
  TEST_EQUAL(bins[2][0].getDriftTime(), -45.0)
}
END_SECTION

START_SECTION((std::vector<PeakMap> splitByFAIMSCV with include_non_faims_cv0=true handles only non-FAIMS data))
{
  PeakMap exp_nonfaims;
  MSSpectrum s;
  s.setMSLevel(1);
  s.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
  s.setDriftTime(10.0);
  exp_nonfaims.addSpectrum(s);

  // Should NOT throw when include_non_faims_cv0=true
  std::vector<PeakMap> bins = IMDataConverter::splitByFAIMSCV(std::move(exp_nonfaims), true);
  TEST_EQUAL(bins.size(), 1) // Just CV=0 group
  TEST_EQUAL(bins[0].size(), 1)
}
END_SECTION

START_SECTION((std::vector<double> getFAIMSCVs returns correct CV values))
{
  PeakMap exp_mixed;

  // Non-FAIMS spectrum
  MSSpectrum non_faims;
  non_faims.setMSLevel(1);
  non_faims.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);

  // FAIMS spectra
  MSSpectrum faims1;
  faims1.setMSLevel(1);
  faims1.setDriftTimeUnit(DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
  faims1.setDriftTime(-55.0);

  MSSpectrum faims2;
  faims2.setMSLevel(1);
  faims2.setDriftTimeUnit(DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
  faims2.setDriftTime(-45.0);

  exp_mixed.addSpectrum(non_faims);
  exp_mixed.addSpectrum(faims1);
  exp_mixed.addSpectrum(faims2);

  // Without non-FAIMS CV=0
  std::vector<double> cvs = IMDataConverter::getFAIMSCVs(exp_mixed, false);
  TEST_EQUAL(cvs.size(), 2)
  TEST_REAL_SIMILAR(cvs[0], -55.0)
  TEST_REAL_SIMILAR(cvs[1], -45.0)

  // With non-FAIMS CV=0
  std::vector<double> cvs_with_zero = IMDataConverter::getFAIMSCVs(exp_mixed, true);
  TEST_EQUAL(cvs_with_zero.size(), 3)
  TEST_REAL_SIMILAR(cvs_with_zero[0], 0.0)
  TEST_REAL_SIMILAR(cvs_with_zero[1], -55.0)
  TEST_REAL_SIMILAR(cvs_with_zero[2], -45.0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
