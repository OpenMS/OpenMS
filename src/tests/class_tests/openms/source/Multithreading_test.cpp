// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------


/**

  Most of the tests, generously provided by the BALL people, taken from version 1.2

*/

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

// OpenMP support
#ifdef _OPENMP
    #include <omp.h>
#endif
#include <climits>

///////////////////////////

using namespace OpenMS;
//using namespace Logger;
using namespace std;


START_TEST(Multithreading, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((OpenMP test))
{
  int wanted_threads = 1;
  int threads = 1;
  #ifdef _OPENMP
  wanted_threads = 2;
  omp_set_dynamic(0);  // Explicitly disable dynamic teams
  // Use 2 threads for all consecutive parallel regions
  omp_set_num_threads(wanted_threads);
  threads = omp_get_max_threads();
  #endif

  int max = INT_MIN;
  int i = 0;
  #pragma omp parallel private(i)
  {
    int maxi = INT_MIN;

    #pragma omp for
    for (i = 0; i < 10; i++)
    {
      #ifdef _OPENMP
      int threadnum = omp_get_thread_num() + 1;
      #else
      int threadnum = 1;
      #endif

      if (threadnum > maxi)
        maxi = threadnum;
    }

    #pragma omp critical
    {
      if (maxi > max)
        max = maxi;
    }
  }

  TEST_EQUAL(threads, wanted_threads)
  TEST_EQUAL(max, wanted_threads)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
