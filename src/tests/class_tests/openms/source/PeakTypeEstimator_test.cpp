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

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(String, "$Id$")

/////////////////////////////////////////////////////////////

PeakTypeEstimator* ptr = nullptr;
PeakTypeEstimator* nullPointer = nullptr;

START_SECTION(([EXTRA]PeakTypeEstimator()))
	ptr = new PeakTypeEstimator();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(([EXTRA] ~PeakTypeEstimator()))
	delete ptr;
END_SECTION

START_SECTION((template<typename PeakConstIterator> SpectrumSettings::SpectrumType estimateType(const PeakConstIterator& begin, const PeakConstIterator& end) const))
	DTAFile file;
	MSSpectrum spec;
	// raw data (with zeros)
	file.load(OPENMS_GET_TEST_DATA_PATH("PeakTypeEstimator_raw.dta"), spec);
  TEST_EQUAL(PeakTypeEstimator::estimateType(spec.begin(), spec.end()), SpectrumSettings::SpectrumType::PROFILE);
  // TOF raw data (without zeros)
	file.load(OPENMS_GET_TEST_DATA_PATH("PeakTypeEstimator_rawTOF.dta"), spec);
  TEST_EQUAL(PeakTypeEstimator::estimateType(spec.begin(), spec.end()), SpectrumSettings::SpectrumType::PROFILE);
  // peak data
	file.load(OPENMS_GET_TEST_DATA_PATH("PeakTypeEstimator_peak.dta"), spec);
  TEST_EQUAL(PeakTypeEstimator::estimateType(spec.begin(), spec.end()), SpectrumSettings::SpectrumType::CENTROID);
  // too few data points
  spec.resize(4);
	TEST_EQUAL(PeakTypeEstimator::estimateType(spec.begin(), spec.end()), SpectrumSettings::SpectrumType::UNKNOWN);
END_SECTION

START_SECTION([EXTRA] estimateType() of peaks of very few sampling points (ion mobility frames))
  // m/z 600 to 608 of a timsTOF MS1 frame (first spectrum of PeakPickerIM_input_1.mzML of the TOPP tests): all mobility
  // scans are stored in one spectrum, i.e. many peaks have the same m/z. An ion covers about three bins of the
  // time-of-flight digitizer, and the intensity drops steeply to the bins next to them.
  DTAFile file;
  MSSpectrum spec;
  file.load(OPENMS_GET_TEST_DATA_PATH("PeakTypeEstimator_rawTIMS.dta"), spec);
  TEST_EQUAL(spec.size(), 4360)
  TEST_EQUAL(PeakTypeEstimator::estimateType(spec.begin(), spec.end()), SpectrumSettings::SpectrumType::PROFILE);

  // the same frame with the mobility scans added up: one peak per bin of the digitizer
  MSSpectrum summed;
  for (const Peak1D& p : spec)
  {
    if (!summed.empty() && summed.back().getMZ() == p.getMZ())
    {
      summed.back().setIntensity(summed.back().getIntensity() + p.getIntensity());
    }
    else
    {
      summed.push_back(p);
    }
  }
  TEST_EQUAL(summed.size(), 1200)
  TEST_EQUAL(PeakTypeEstimator::estimateType(summed.begin(), summed.end()), SpectrumSettings::SpectrumType::PROFILE);

  // five ions with five sampling points each: @p spacing Th apart, with intensities which drop steeply.
  // With @p skip, the leftmost sampling point of each ion is one grid position further away (a point without intensity,
  // which is not stored, is in between).
  auto steepPeaks = [](const double spacing, const bool skip)
  {
    const float intensities[] = {30, 2500, 40000, 3000, 50};
    MSSpectrum s;
    for (Size ion = 0; ion < 5; ++ion)
    {
      for (int i = 0; i < 5; ++i)
      {
        const int grid_position = (i == 0 && skip) ? -3 : i - 2;
        s.emplace_back(500.0 + ion * 10.0 + grid_position * spacing, intensities[i] * float(ion + 1));
      }
    }
    return s;
  };
  // spaced like the bins of a digitizer
  MSSpectrum steep = steepPeaks(0.004, false);
  TEST_EQUAL(PeakTypeEstimator::estimateType(steep.begin(), steep.end()), SpectrumSettings::SpectrumType::PROFILE);
  steep = steepPeaks(0.004, true);
  TEST_EQUAL(PeakTypeEstimator::estimateType(steep.begin(), steep.end()), SpectrumSettings::SpectrumType::PROFILE);
  // ... but not with three points missing: no sampling grid
  steep = steepPeaks(0.004, false);
  for (Size i = 0; i < steep.size(); i += 5)
  {
    steep[i].setMZ(steep[i].getMZ() - 3 * 0.004);
  }
  TEST_EQUAL(PeakTypeEstimator::estimateType(steep.begin(), steep.end()), SpectrumSettings::SpectrumType::CENTROID);
  // ... and spaced like the isotopic peaks of ions with charge 5: centroided
  steep = steepPeaks(0.2, false);
  TEST_EQUAL(PeakTypeEstimator::estimateType(steep.begin(), steep.end()), SpectrumSettings::SpectrumType::CENTROID);

  // centroided isotopic peaks of a triply charged ion are equally spaced as well
  // (a spectrum of HighResPrecursorMassCorrector_3070_1191_3.mzML of the TOPP tests)
  MSSpectrum isotopes;
  isotopes.emplace_back(1191.5837, 1672920.0f);
  isotopes.emplace_back(1191.9164, 2968130.0f);
  isotopes.emplace_back(1192.2468, 2974250.0f);
  isotopes.emplace_back(1192.5872, 2070300.0f);
  isotopes.emplace_back(1192.9198, 1612660.0f);
  TEST_EQUAL(PeakTypeEstimator::estimateType(isotopes.begin(), isotopes.end()), SpectrumSettings::SpectrumType::CENTROID);

  // peaks at the same m/z are one sampling point: five of them are not five peaks of a profile
  MSSpectrum same_mz;
  for (int i = 0; i < 6; ++i)
  {
    same_mz.emplace_back(500.0, 100.0f + i);
  }
  TEST_EQUAL(PeakTypeEstimator::estimateType(same_mz.begin(), same_mz.end()), SpectrumSettings::SpectrumType::CENTROID);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
