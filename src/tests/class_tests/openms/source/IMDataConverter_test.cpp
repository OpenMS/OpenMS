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


START_SECTION((std::vector<std::pair<double, MSExperiment>> splitByFAIMSCV(PeakMap& exp)))
  MzMLFile IM_file;
  PeakMap exp;
  IM_file.load(OPENMS_GET_TEST_DATA_PATH("IM_FAIMS_test.mzML"), exp);

  TEST_EQUAL(exp.getSpectra().size(), 19)

  auto splitPeakMap = IMDataConverter::splitByFAIMSCV(std::move(exp));
  TEST_EQUAL(exp.empty(), true) // moved out
  TEST_EQUAL(splitPeakMap.size(), 3)

  // expect keys -65, -55, -45 in ascending order
  TEST_EQUAL(splitPeakMap[0].first, -65.0)
  TEST_EQUAL(splitPeakMap[1].first, -55.0)
  TEST_EQUAL(splitPeakMap[2].first, -45.0)

  TEST_EQUAL(splitPeakMap[0].second.size(), 4)
  TEST_EQUAL(splitPeakMap[1].second.size(), 9)
  TEST_EQUAL(splitPeakMap[2].second.size(), 6)

  for (const auto& spec : splitPeakMap[0].second)
  {
    TEST_EQUAL(spec.getDriftTime(), -65.0)
  }
  for (const auto& spec : splitPeakMap[1].second)
  {
    TEST_EQUAL(spec.getDriftTime(), -55.0)
  }
  for (const auto& spec : splitPeakMap[2].second)
  {
    TEST_EQUAL(spec.getDriftTime(), -45.0)
  }

  TEST_EQUAL(splitPeakMap[1].second.getExperimentalSettings().getDateTime().toString(), "2019-09-07T09:40:04")

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

START_SECTION((std::vector<std::pair<double, MSExperiment>> splitByFAIMSCV assigns MS2 without explicit CV to last seen FAIMS CV))
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

  auto bins = IMDataConverter::splitByFAIMSCV(std::move(exp_synth));
  TEST_EQUAL(bins.size(), 2)

  // bins ordered by ascending CV: -55 first, -45 second
  TEST_EQUAL(bins[0].first, -55.0)
  TEST_EQUAL(bins[1].first, -45.0)

  const auto& bin_minus55 = bins[0].second;
  const auto& bin_minus45 = bins[1].second;

  TEST_EQUAL(bin_minus55.size(), 3) // ms1a + ms2a + ms2b
  TEST_EQUAL(bin_minus45.size(), 2) // ms1b + ms2c

  Size ms2_bin0 = 0;
  for (const auto& s : bin_minus55) if (s.getMSLevel() > 1) ++ms2_bin0;
  TEST_EQUAL(ms2_bin0, 2)

  Size ms2_bin1 = 0;
  for (const auto& s : bin_minus45) if (s.getMSLevel() > 1) ++ms2_bin1;
  TEST_EQUAL(ms2_bin1, 1)
}
END_SECTION

START_SECTION((std::vector<std::pair<double, MSExperiment>> splitByFAIMSCV returns single-element group for non-FAIMS dataset))
{
  PeakMap exp_nonfaims;
  MSSpectrum s;
  s.setMSLevel(1);
  s.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
  s.setDriftTime(10.0);
  exp_nonfaims.addSpectrum(s);

  auto bins = IMDataConverter::splitByFAIMSCV(std::move(exp_nonfaims));
  TEST_EQUAL(bins.size(), 1)
  TEST_EQUAL(bins[0].second.size(), 1)
}
END_SECTION

START_SECTION(static void convertVSSCToCCS(MSExperiment& spectra))
{
  // Create a simple MSExperiment with spectra that have VSSC drift times
  MSExperiment exp;
  
  // Spectrum 1: m/z = 500, charge = 2, drift time (1/k0) = 1.0
  MSSpectrum s1;
  s1.setMSLevel(2);
  s1.setDriftTime(1.0);
  s1.setDriftTimeUnit(DriftTimeUnit::VSSC);
  Precursor p1;
  p1.setMZ(500.0);
  p1.setCharge(2);
  s1.setPrecursors({p1});
  exp.addSpectrum(s1);
  
  // Spectrum 2: m/z = 750, charge = 3, drift time (1/k0) = 1.5
  MSSpectrum s2;
  s2.setMSLevel(2);
  s2.setDriftTime(1.5);
  s2.setDriftTimeUnit(DriftTimeUnit::VSSC);
  Precursor p2;
  p2.setMZ(750.0);
  p2.setCharge(3);
  s2.setPrecursors({p2});
  exp.addSpectrum(s2);
  
  // Spectrum 3: m/z = 400, charge = 0 (should be skipped)
  MSSpectrum s3;
  s3.setMSLevel(2);
  s3.setDriftTime(0.8);
  s3.setDriftTimeUnit(DriftTimeUnit::VSSC);
  Precursor p3;
  p3.setMZ(400.0);
  p3.setCharge(0); // Unknown charge
  s3.setPrecursors({p3});
  exp.addSpectrum(s3);
  
  // Spectrum 4: no precursor (should be skipped)
  MSSpectrum s4;
  s4.setMSLevel(1);
  s4.setDriftTime(0.5);
  s4.setDriftTimeUnit(DriftTimeUnit::VSSC);
  exp.addSpectrum(s4);
  
  // Spectrum 5: negative drift time (should be skipped)
  MSSpectrum s5;
  s5.setMSLevel(2);
  s5.setDriftTime(-1.0);
  s5.setDriftTimeUnit(DriftTimeUnit::VSSC);
  Precursor p5;
  p5.setMZ(500.0);
  p5.setCharge(2);
  s5.setPrecursors({p5});
  exp.addSpectrum(s5);
  
  // Spectrum 6: zero drift time (should be skipped)
  MSSpectrum s6;
  s6.setMSLevel(2);
  s6.setDriftTime(0.0);
  s6.setDriftTimeUnit(DriftTimeUnit::VSSC);
  Precursor p6;
  p6.setMZ(500.0);
  p6.setCharge(2);
  s6.setPrecursors({p6});
  exp.addSpectrum(s6);
  
  // Run conversion
  IMDataConverter::convertVSSCToCCS(exp);
  
  // Verify first spectrum: manually calculate expected CCS
  // mass = 500 * 2 = 1000
  // reduced_mass = (1000 * 28) / (1000 + 28) = 27.237...
  // CCS = 1.0 * 2 * 1059.62245 / sqrt(27.237...) ≈ 406.09
  TEST_REAL_SIMILAR(exp[0].getDriftTime(), 406.09) // Expected CCS for spectrum 1
  
  // Verify second spectrum
  // mass = 750 * 3 = 2250
  // reduced_mass = (2250 * 28) / (2250 + 28) = 27.656...
  // CCS = 1.5 * 3 * 1059.62245 / sqrt(27.656...) ≈ 905.22
  TEST_REAL_SIMILAR(exp[1].getDriftTime(), 905.22) // Expected CCS for spectrum 2
  
  // Third spectrum (charge 0) should not be converted - stays the same
  TEST_REAL_SIMILAR(exp[2].getDriftTime(), 0.8)
  
  // Fourth spectrum (no precursor) should not be converted - stays the same
  TEST_REAL_SIMILAR(exp[3].getDriftTime(), 0.5)
  
  // Fifth spectrum (negative drift time) should not be converted - stays the same
  TEST_REAL_SIMILAR(exp[4].getDriftTime(), -1.0)
  
  // Sixth spectrum (zero drift time) should not be converted - stays the same
  TEST_REAL_SIMILAR(exp[5].getDriftTime(), 0.0)
  
  // Test IM float data array conversion
  MSExperiment exp_fda;
  MSSpectrum s_fda;
  s_fda.setMSLevel(2);
  s_fda.setDriftTime(1.0);
  s_fda.setDriftTimeUnit(DriftTimeUnit::VSSC);
  Precursor p_fda;
  p_fda.setMZ(500.0);
  p_fda.setCharge(2);
  s_fda.setPrecursors({p_fda});
  
  // Add IM float data array with VSSC values
  MSSpectrum::FloatDataArray& im_fda = s_fda.getFloatDataArrays().emplace_back();
  IMDataConverter::setIMUnit(im_fda, DriftTimeUnit::VSSC);
  im_fda.push_back(1.0f); // same as spectrum drift time
  im_fda.push_back(1.2f);
  im_fda.push_back(0.0f); // invalid - should stay same
  
  exp_fda.addSpectrum(s_fda);
  IMDataConverter::convertVSSCToCCS(exp_fda);
  
  // Spectrum-level drift time should be converted
  TEST_REAL_SIMILAR(exp_fda[0].getDriftTime(), 406.09)
  
  // Float data array values should also be converted
  // CCS = 1.0 * 2 * 1059.62245 / sqrt(27.237...) ≈ 406.09
  TEST_REAL_SIMILAR(exp_fda[0].getFloatDataArrays()[0][0], 406.09)
  // CCS = 1.2 * 2 * 1059.62245 / sqrt(27.237...) ≈ 487.31
  TEST_REAL_SIMILAR(exp_fda[0].getFloatDataArrays()[0][1], 487.31)
  // Invalid value (0.0) should stay the same
  TEST_REAL_SIMILAR(exp_fda[0].getFloatDataArrays()[0][2], 0.0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
