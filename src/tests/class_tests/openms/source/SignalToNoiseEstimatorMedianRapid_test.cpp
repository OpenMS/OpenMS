// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedianRapid.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SignalToNoiseEstimatorMedianRapid, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SignalToNoiseEstimatorMedianRapid* ptr = nullptr;
SignalToNoiseEstimatorMedianRapid* nullPointer = nullptr;
START_SECTION((SignalToNoiseEstimatorMedianRapid()))
	ptr = new SignalToNoiseEstimatorMedianRapid(200);
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~SignalToNoiseEstimatorMedianRapid()))
	delete ptr;
END_SECTION

/*
 * Python code:
 *
 
import random
mz = [200 + 10*i for i in range(40)]
int = [ random.random() * 10 for i in range(40)]
int = [5.4332, 5.6189, 4.3025, 4.5705, 5.4538, 9.7202, 8.805, 8.5391, 6.6257, 5.809, 6.5518, 7.9273, 5.3875, 9.826, 5.139, 5.8588, 0.7806, 4.2054, 9.9171, 4.0198, 1.1462, 5.1042, 7.8318, 4.8553, 6.691, 4.2377, 7.2344, 4.0124, 3.8565, 6.2867, 1.0817, 8.2412, 5.0589, 7.0478, 5.9388, 1.2747, 2.4228, 4.909, 6.856, 1.9665]

import numpy
numpy.median( int[0:21] )
numpy.median( int[21:40] )

*/

START_SECTION( (NoiseEstimator estimateNoise(std::vector<double>& mz_array, std::vector<double>& int_array)))
{
  static const double arr1[] = 
  {
    200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340,
    350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490,
    500, 510, 520, 530, 540, 550, 560, 570, 580, 590
  };
  std::vector<double> mz (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
  static const double arr2[] = 
  {
    5.4332, 5.6189, 4.3025, 4.5705, 5.4538, 9.7202, 8.805, 8.5391, 6.6257,
    5.809, 6.5518, 7.9273, 5.3875, 9.826, 5.139, 5.8588, 0.7806, 4.2054,
    9.9171, 4.0198, 1.1462, 5.1042, 7.8318, 4.8553, 6.691, 4.2377, 7.2344,
    4.0124, 3.8565, 6.2867, 1.0817, 8.2412, 5.0589, 7.0478, 5.9388, 1.2747,
    2.4228, 4.909, 6.856, 1.9665
  };
  std::vector<double> intensity (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );

  // Large window (200)
  {
    SignalToNoiseEstimatorMedianRapid sne(200);
    SignalToNoiseEstimatorMedianRapid::NoiseEstimator e = sne.estimateNoise(mz, intensity);
    TEST_REAL_SIMILAR(e.get_noise_even(200), 5.71395) // numpy.median(int[:20])
    TEST_REAL_SIMILAR(e.get_noise_even(500), 4.98395) // numpy.median(int[20:])

    TEST_REAL_SIMILAR(e.get_noise_odd(200), 5.71395)  // numpy.median(int[:10])
    TEST_REAL_SIMILAR(e.get_noise_odd(400), 5.26325)  // numpy.median(int[10:30])
    TEST_REAL_SIMILAR(e.get_noise_odd(500), 4.98395)  // numpy.median(int[30:])

    TEST_REAL_SIMILAR(e.get_noise_value(200), 5.71395) // numpy.median(int[:20])
    TEST_REAL_SIMILAR(e.get_noise_value(410), (5.26325+4.98395)/2 ) // (numpy.median(int[10:30])+numpy.median(int[30:])) /2
    TEST_REAL_SIMILAR(e.get_noise_value(500), 4.98395)  // numpy.median(int[30:])

  }

  // Smaller window (100)
  {
    SignalToNoiseEstimatorMedianRapid sne(100);
    SignalToNoiseEstimatorMedianRapid::NoiseEstimator e = sne.estimateNoise(mz, intensity);
    TEST_REAL_SIMILAR(e.get_noise_even(250), 5.71395)  // numpy.median( int[:10] )
    TEST_REAL_SIMILAR(e.get_noise_even(350), 5.62315 ) // numpy.median( int[10:20] )
    TEST_REAL_SIMILAR(e.get_noise_even(450), 4.97975 ) // numpy.median( int[20:30] )
    TEST_REAL_SIMILAR(e.get_noise_even(550), 4.98395)  // numpy.median( int[30:] )

    TEST_REAL_SIMILAR(e.get_noise_odd(200), 5.4332)  // numpy.median( int[:5] )
    TEST_REAL_SIMILAR(e.get_noise_odd(300), 7.2765)  // numpy.median( int[5:15] )
    TEST_REAL_SIMILAR(e.get_noise_odd(400), 4.97975)  // numpy.median( int[15:25] )
    TEST_REAL_SIMILAR(e.get_noise_odd(500), 5.49885)  // numpy.median( int[25:35] )

    TEST_REAL_SIMILAR(e.get_noise_value(510), (5.49885+4.98395)/2)  // (numpy.median( int[25:35] ) + numpy.median( int[30:] ) )/2
  }

  // Uneven window size (50)
  {
    SignalToNoiseEstimatorMedianRapid sne(50);
    SignalToNoiseEstimatorMedianRapid::NoiseEstimator e = sne.estimateNoise(mz, intensity);
    TEST_REAL_SIMILAR(e.get_noise_even(220), 5.4332 ) // numpy.median( int[:5] )
    TEST_REAL_SIMILAR(e.get_noise_even(420), 5.1042 ) // numpy.median(int[20:25])
    TEST_REAL_SIMILAR(e.get_noise_even(460), 4.2377 ) // numpy.median(int[25:30])
  }

  // Uneven window size (110)
  {
    SignalToNoiseEstimatorMedianRapid sne(110);
    SignalToNoiseEstimatorMedianRapid::NoiseEstimator e = sne.estimateNoise(mz, intensity);
    TEST_REAL_SIMILAR(e.get_noise_even(250), 5.809  ) // numpy.median( int[:11] )
    TEST_REAL_SIMILAR(e.get_noise_even(350), 5.139  ) // numpy.median( int[11:22] )
    TEST_REAL_SIMILAR(e.get_noise_even(450), 5.05890) // numpy.median( int[22:33] )
    TEST_REAL_SIMILAR(e.get_noise_even(550), 4.909  ) // numpy.median( int[33:] )
  }

}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


