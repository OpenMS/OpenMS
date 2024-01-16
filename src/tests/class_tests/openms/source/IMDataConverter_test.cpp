// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
frame.push_back({1.0, 29.0f});
frame.push_back({2.0, 60.0f});
frame.push_back({3.0, 34.0f});
frame.push_back({4.0, 29.0f});
frame.push_back({5.0, 37.0f});
frame.push_back({6.0, 31.0f});
frame.setRT(1);
MSSpectrum::FloatDataArray& afa = frame.getFloatDataArrays().emplace_back();
afa.assign({1.1, 2.2, 3.3, 3.3, 5.5, 6.6});
IMDataConverter::setIMUnit(afa, DriftTimeUnit::MILLISECOND);

MSSpectrum spec;
spec.push_back({111.0, -1.0f});
spec.push_back({222.0, -2.0f});
spec.push_back({333.0, -3.0f});
spec.setRT(2); // just a spectrum with RT = 2


START_SECTION(static MSExperiment splitByIonMobility(MSSpectrum im_frame, UInt number_of_bins = -1))
	
TEST_EXCEPTION(Exception::MissingInformation, IMDataConverter::splitByIonMobility(spec))
	{
		auto exp = IMDataConverter::splitByIonMobility(frame);
		TEST_EQUAL(exp.size(), 5);
		TEST_EQUAL(exp[0].size(), 1);
		TEST_EQUAL(exp[2].size(), 2);
		TEST_EQUAL(exp[0][0].getIntensity(), 29.0f);

		TEST_EQUAL(exp[0].getDriftTime(), 1.1f);
		TEST_EQUAL(exp[0].getDriftTimeUnit() == DriftTimeUnit::MILLISECOND, true);
		TEST_EQUAL(exp[0].getRT(), 1);

		auto frame_reconstruct = IMDataConverter::collapseFramesToSingle(exp);
		TEST_EQUAL(frame_reconstruct.size(), 1)
		TEST_EQUAL(frame_reconstruct[0], frame);
	}
  {
		auto exp_binned = IMDataConverter::splitByIonMobility(frame, 1);
		TEST_EQUAL(exp_binned.size(), 1);
		TEST_EQUAL(exp_binned[0].size(), frame.size());
		TEST_EQUAL(exp_binned[0][0].getIntensity(), 29.0f);
		TEST_REAL_SIMILAR(exp_binned[0].getDriftTime(), (6.6-1.1)/2 + 1.1);
		TEST_EQUAL(exp_binned[0].getDriftTimeUnit() == DriftTimeUnit::MILLISECOND, true);
		TEST_EQUAL(exp_binned[0].getRT(), 1);
	}
END_SECTION

START_SECTION(static MSExperiment splitByIonMobility(MSExperiment&& in, UInt number_of_bins = -1))
	MSExperiment e_in;
	e_in.addSpectrum(frame);
  e_in.addSpectrum(spec); // just copy it...
  auto frame3 = frame;
  frame3.setRT(3);
  e_in.addSpectrum(frame3);
  
	auto exp = IMDataConverter::splitByIonMobility(std::move(e_in));
	TEST_EQUAL(exp.size(), 5+1+5);
	TEST_EQUAL(exp[0].size(), 1);
	TEST_EQUAL(exp[2].size(), 2);
	TEST_EQUAL(exp[0][0].getIntensity(), 29.0f);
  TEST_EQUAL(exp[0].getDriftTime(), 1.1f);
  TEST_EQUAL(exp[0].getDriftTimeUnit() == DriftTimeUnit::MILLISECOND, true);
  TEST_EQUAL(exp[0].getRT(), 1);
  
	TEST_EQUAL(exp[5], spec); // copied

	TEST_EQUAL(exp[6].getRT(), 3);					  

	auto frame_reconstruct = IMDataConverter::collapseFramesToSingle(exp);
	TEST_EQUAL(frame_reconstruct.size(), 3)
	TEST_EQUAL(frame_reconstruct[0] == frame, true);
  TEST_EQUAL(frame_reconstruct[1] == spec, true);
  TEST_EQUAL(frame_reconstruct[2] == frame3, true);
END_SECTION

START_SECTION(static MSExperiment collapseFramesToSingle(const MSExperiment& in))
	NOT_TESTABLE // tested_above
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
