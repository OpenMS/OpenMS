// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Müller, Jaekwan Kim$
// $Authors: Tom David Müller, Jaekwan Kim$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>

///////////////////////////

START_TEST(FLASHDeconvAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

FLASHDeconvAlgorithm* ptr = nullptr;
FLASHDeconvAlgorithm* null_ptr = nullptr;

START_SECTION(FLASHDeconvAlgorithm())
  ptr = new FLASHDeconvAlgorithm();
  TEST_NOT_EQUAL(ptr, null_ptr)
END_SECTION

START_SECTION(FLASHDeconvAlgorithm(const FLASHDeconvAlgorithm& source))
FLASHDeconvAlgorithm copy(*ptr);
  TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(FLASHDeconvAlgorithm& operator=(const FLASHDeconvAlgorithm& source))
FLASHDeconvAlgorithm copy;
  copy = *ptr;
  TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(~FLASHDeconvAlgorithm())
  delete ptr;
END_SECTION

ptr = new FLASHDeconvAlgorithm();
START_SECTION(FLASHDeconvAlgorithm(FLASHDeconvAlgorithm&& source))
  FLASHDeconvAlgorithm temp;
  temp.setParameters(ptr->getParameters());
  FLASHDeconvAlgorithm moved(std::move(temp));
  TEST_EQUAL(moved.getParameters(), ptr->getParameters());
END_SECTION

START_SECTION(FLASHDeconvAlgorithm& operator=(FLASHDeconvAlgorithm&& source))
  FLASHDeconvAlgorithm temp;
  temp.setParameters(ptr->getParameters());
  FLASHDeconvAlgorithm moved;
  moved = std::move(temp);
  TEST_EQUAL(moved.getParameters(), ptr->getParameters());
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ptr = new FLASHDeconvAlgorithm();
Param params;
START_SECTION(std::vector<double> getTolerances())
  params.setValue("SD:tol", ListUtils::create<double>("10.0,5.0"));
  ptr->setParameters(params);
  auto tolerances = ptr->getTolerances();
  TEST_EQUAL(tolerances.size(), 2)
  TEST_REAL_SIMILAR(tolerances[0], 10.0);
  TEST_REAL_SIMILAR(tolerances[1], 5.0);
END_SECTION

// load test data
// TODO: minimize
PeakMap input;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("FLASHDeconv_1_input.mzML"), input);

// Store FD outputs
std::vector<DeconvolvedSpectrum> deconvolved_spectra;
std::vector<FLASHHelperClasses::MassFeature> deconvolved_features;

START_SECTION(void run())
  ptr->run(input, deconvolved_spectra, deconvolved_features);
  TEST_EQUAL(deconvolved_spectra.size(), input.size());
  TEST_FALSE(deconvolved_features.empty());
END_SECTION

START_SECTION(int getScanNumber())
  TEST_EQUAL(ptr->getScanNumber(input, 0), 1);
END_SECTION

// Check if the precalculated averagine was initialized during runtime
START_SECTION(FLASHHelperClasses::PrecalculatedAveragine getAveragine())
  TEST_NOT_EQUAL(ptr->getAveragine().getMaxIsotopeIndex(), 0);
  TEST_FALSE(ptr->getAveragine().get(500.0).getContainer().empty());
END_SECTION

START_SECTION(FLASHHelperClasses::PrecalculatedAveragine getDecoyAveragine())
  TEST_NOT_EQUAL(ptr->getDecoyAveragine().getMaxIsotopeIndex(), 0);
  TEST_FALSE(ptr->getDecoyAveragine().get(500.0).getContainer().empty());
END_SECTION

// Decoy Averagine needs to be handled differently if FDR is reported
START_SECTION(FLASHHelperClasses::PrecalculatedAveragine getDecoyAveragine())
  params.setValue("FD:report_FD", true);
  ptr = new FLASHDeconvAlgorithm();  
  ptr->setParameters(params);
  ptr->run(input, deconvolved_spectra, deconvolved_features);
  TEST_NOT_EQUAL(ptr->getDecoyAveragine().getMaxIsotopeIndex(), 0);
  TEST_FALSE(ptr->getDecoyAveragine().get(500.0).getContainer().empty());
END_SECTION


// // TODO: Move
// START_SECTION(FLASHHelperClasses::PrecalculatedAveragine& getAveragine())
//   TEST_EQUAL(ptr->getAveragine().getMaxIsotopeIndex(), 199);
//   TEST_EQUAL(ptr->getAveragine().getLeftCountFromApex(50), 2);
//   TOLERANCE_ABSOLUTE(0.001)
//   TEST_REAL_SIMILAR(fd_algo->getAveragine().getAverageMassDelta(50), 0.0251458);
//   TEST_REAL_SIMILAR(fd_algo->getAveragine().getSNRMultiplicationFactor(50), 1.051538);
// END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
