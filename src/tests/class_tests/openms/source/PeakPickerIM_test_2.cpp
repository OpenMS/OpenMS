// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Author: Timo Sachsenberg, Mohammed Alhigaylan $
// $Maintainer: Timo Sachsenberg $
// -------------------------------------------------------------------------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

///////////////////////////
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerIM.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PeakPickerIM, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakPickerIM* ptr = nullptr;
PeakPickerIM* nullPointer = nullptr;

START_SECTION((PeakPickerIM()))
  ptr = new PeakPickerIM();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~PeakPickerIM()))
  delete ptr;
END_SECTION


// create dummy ion mobility spectrum
MSSpectrum input;
{
  double start_mz = 0.8;
  double end_mz = 1.2;
  double step = 0.01;
  double max_intensity = 1000.0; // Maximum intensity at m/z = 1.0

  for (double mz = start_mz; mz <= end_mz; mz += step)
  {
    double intensity = 0.0;
    if (mz <= 1.0)
    {
      // Linearly increase intensity from start_mz to 1.0
      intensity = ((mz - start_mz) / (1.0 - start_mz)) * max_intensity;
    }
    else
    {
      // Linearly decrease intensity from 1.0 to end_mz
      intensity = ((end_mz - mz) / (end_mz - 1.0)) * max_intensity;
    }
    input.emplace_back(mz, intensity);
  }
  // Set the IM format to CONCATENATED to force the peak picking branch
  input.setIMFormat(IMFormat::CONCATENATED);
//  input.getFloatDataArrays().resize(1);
//  input.getFloatDataArrays()[0].setName("Ion Mobility");
}


START_SECTION(void pickIMTraces(MSSpectrum& spectrum))
{
    PeakPickerIM pp_im;
    // print all peaks in our current input.
    std::cout << "start printing dummy spectrum BEFORE picking! " << std::endl;
    for (const auto& peak : input)
    {
        std::cout << "m/z: " << peak.getMZ()
                  << ", intensity: " << peak.getIntensity() << std::endl;
    }

    pp_im.pickIMTraces(input);
    std::cout << "start printing dummy spectrum AFTER picking! " << std::endl;
    for (const auto& peak : input)
    {
        std::cout << "m/z: " << peak.getMZ()
                  << ", intensity: " << peak.getIntensity() << std::endl;
    }

    // TODO adapt
    TEST_EQUAL(input.size(), 10)
    TEST_REAL_SIMILAR(input[0].getIntensity(), 450)
    TEST_REAL_SIMILAR(input[0].getMZ(), 100.02)

//    TEST_EQUAL(input.getFloatDataArrays().size(), 1)
 //   TEST_EQUAL(input.getFloatDataArrays()[0].getName(), "Ion Mobility")
 //   TEST_REAL_SIMILAR(input.getFloatDataArrays()[0][0], 150.0)
}
END_SECTION

END_TEST