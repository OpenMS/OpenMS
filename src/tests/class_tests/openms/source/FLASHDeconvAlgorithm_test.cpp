// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Jaekwan Kim$
// $Authors: Jaekwan Kim$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/FORMAT/MzMLFile.h>
//#include <OpenMS/FORMAT/MzMLFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FLASHDeconvAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FLASHDeconvAlgorithm* ptr = 0;
FLASHDeconvAlgorithm* null_ptr = 0;
START_SECTION(FLASHDeconvAlgorithm())
{
  ptr = new FLASHDeconvAlgorithm();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~FLASHDeconvAlgorithm())
{
  delete ptr;
}
END_SECTION

FLASHDeconvAlgorithm fd_algo = FLASHDeconvAlgorithm();


START_SECTION(double getNoiseDecoyWeight())
{
  TEST_EQUAL(fd_algo.getNoiseDecoyWeight(),1);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// load test data
PeakMap input;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("thermo.mzML"), input);

START_SECTION(int getScanNumber())
{
  TEST_EQUAL(fd_algo.getScanNumber(input, 0), 1);
}
END_SECTION
std::vector<DeconvolvedSpectrum> deconvolved_spectra;
std::vector<FLASHDeconvHelperStructs::MassFeature> deconvolved_features;

START_SECTION(void run())
{
  fd_algo.run(input, deconvolved_spectra, deconvolved_features);
  TEST_EQUAL(deconvolved_spectra.size(), 2432);
  TEST_EQUAL(deconvolved_features.size(), 6);
}
END_SECTION

START_SECTION(FLASHDeconvHelperStructs::PrecalculatedAveragine& getAveragine())
{
  TEST_EQUAL(fd_algo.getAveragine().getMaxIsotopeIndex(), 199);
  TEST_EQUAL(fd_algo.getAveragine().getLeftCountFromApex(50), 2);
  TOLERANCE_ABSOLUTE(0.001)
  TEST_REAL_SIMILAR(fd_algo.getAveragine().getAverageMassDelta(50), 0.0251458);
  TEST_REAL_SIMILAR(fd_algo.getAveragine().getSNRMultiplicationFactor(50), 1.051538);

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
