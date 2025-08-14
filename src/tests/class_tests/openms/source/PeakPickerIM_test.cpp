// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Author: Timo Sachsenberg $
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
    const size_t points_per_im = 10;

    // For each ion mobility value (1.0 to 10.0)
    for (double im = 1.0; im <= 10.0; im += 1.0)
    {
        double baseIntensity = (im <= 5.0) ? 
                            200.0 + ((im - 1.0) * 200.0) :    // Increasing intensity from IM 1.0 to 5.0
                            1000.0 - ((im - 5.0) * 200.0);    // Decreasing intensity from IM 6.0 to 10.0
        
        // Create 10 data points for each ion mobility value
        for (size_t i = 0; i < points_per_im; i++)
        {
            double mz = 100.0 + (i * 0.01);
            double intensity = baseIntensity + (i * 20.0);
            input.emplace_back(mz, intensity);
        }
    }

    // Set up the ion mobility array
    input.getFloatDataArrays().resize(1);
    input.getFloatDataArrays()[0].setName("Ion Mobility");

    // Add ion mobility values (each repeated 10 times)
    for (double im = 1.0; im <= 10.0; im += 1.0)
    {
        for (size_t i = 0; i < points_per_im; i++)
        {
            input.getFloatDataArrays()[0].push_back(im);
        }
    }
}

START_SECTION(void pickIMTraces(MSSpectrum& spectrum))
{
    PeakPickerIM pp_im;
    pp_im.pickIMTraces(input);
    // TODO adapt
    TEST_EQUAL(input.size(), 10)
    TEST_REAL_SIMILAR(input[0].getIntensity(), 450)
    TEST_REAL_SIMILAR(input[0].getMZ(), 100.02)

    TEST_EQUAL(input.getFloatDataArrays().size(), 1)
    TEST_EQUAL(input.getFloatDataArrays()[0].getName(), "Ion Mobility")
    TEST_REAL_SIMILAR(input.getFloatDataArrays()[0][0], 150.0)
}
END_SECTION

END_TEST