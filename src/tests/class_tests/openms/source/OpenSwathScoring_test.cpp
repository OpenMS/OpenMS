// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------
#include <boost/shared_ptr.hpp>

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathScoring.h>
///////////////////////////

// we don't want any inclusion of OpenMS Kernel classes here ...
#ifdef OPENMS_KERNEL_MSSPECTRUM_H
ThisShouldFailAtCompileTime = 0
#endif
#ifdef OPENMS_KERNEL_MRMFEATURE_H
ThisShouldFailAtCompileTime = 0
#endif

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>

using namespace OpenMS;
using namespace std;

MSSpectrum generateImSpec(int k_min, int k_max, double rt)
{
  MSSpectrum imSpec;
  DataArrays::FloatDataArray fda;
  imSpec.setRT(rt);

  for (int k = k_min; k < k_max; k++)
  {
    Peak1D p;
    p.setMZ(100. + k);
    p.setIntensity(1.);
    imSpec.push_back(p);
    fda.push_back((double) k);
    fda.setName("Ion Mobility");
  }
  imSpec.getFloatDataArrays().push_back(fda);

  return imSpec;
}



START_TEST(OpenSwathScoring, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

OpenSwathScoring* ptr = nullptr;
OpenSwathScoring* nullPointer = nullptr;

START_SECTION(OpenSwathScoring())
{
	ptr = new OpenSwathScoring();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~OpenSwathScoring())
{
  delete ptr;
}
END_SECTION

START_SECTION(( void initialize(double rt_normalization_factor, int add_up_spectra, double spacing_for_spectra_resampling, const OpenSwath_Scores_Usage & su, const std::string& spectrum_addition_method, bool use_ms1_ion_mobility) ))
{
	ptr = new OpenSwathScoring();
  OpenSwath_Scores_Usage su;
	TEST_NOT_EQUAL(ptr, nullPointer)
  ptr->initialize(100.0, 1, 0.01, 0.0, su, "simple", true);
  delete ptr;
}
END_SECTION

START_SECTION((void calculateChromatographicScores( OpenSwath::IMRMFeature* imrmfeature, const std::vector<std::string>& native_ids, const std::vector<double>& normalized_library_intensity, std::vector<OpenSwath::ISignalToNoisePtr>& signal_noise_estimators, OpenSwath_Scores & scores) ))
{
  NOT_TESTABLE // see MRMFeatureFinderScoring_test.cpp
  // - the OpenSwathScoring is a facade and thus does not need testing on its own
}
END_SECTION

START_SECTION((void calculateChromatographicIdScores( OpenSwath::IMRMFeature* imrmfeature, const std::vector<std::string>& native_ids_identification,, const std::vector<std::string>& native_ids_detection, std::vector<OpenSwath::ISignalToNoisePtr>& signal_noise_estimators, OpenSwath_Scores & idscores) ))
{
  NOT_TESTABLE // see MRMFeatureFinderScoring_test.cpp
  // - the OpenSwathScoring is a facade and thus does not need testing on its own
}
END_SECTION

START_SECTION((void calculateLibraryScores( OpenSwath::IMRMFeature* imrmfeature, const std::vector<TransitionType> & transitions, const PeptideType& pep, const double normalized_feature_rt, OpenSwath_Scores & scores)))
{
  NOT_TESTABLE // see MRMFeatureFinderScoring_test.cpp
  // - the OpenSwathScoring is a facade and thus does not need testing on its own
}
END_SECTION

START_SECTION((void calculateDIAScores(OpenSwath::IMRMFeature* imrmfeature, const std::vector<TransitionType> & transitions, OpenSwath::SpectrumAccessPtr swath_map, OpenMS::DIAScoring & diascoring, const PeptideType& pep, OpenSwath_Scores & scores)))
{
  NOT_TESTABLE // see MRMFeatureFinderScoring_test.cpp
  // - the OpenSwathScoring is a facade and thus does not need testing on its own
}
END_SECTION

START_SECTION((void getNormalized_library_intensities_(const std::vector<TransitionType> & transitions, std::vector<double>& normalized_library_intensity)))
{
  NOT_TESTABLE // see MRMFeatureFinderScoring_test.cpp
  // - the OpenSwathScoring is a facade and thus does not need testing on its own
}
END_SECTION

START_SECTION((OpenSwath::SpectrumPtr OpenSwathScoring::fetchSpectrumSwath(std::vector<OpenSwath::SwathMap> swath_maps,
                                                              double RT, int nr_spectra_to_add, RangeMobility im_range)))
{

  OpenMS::RangeMobility im_range_empty; // use this empty im range as input for all examples
  // test result for empty map
  {
    boost::shared_ptr<PeakMap > swath_map (new PeakMap);
    OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);

    OpenSwathScoring sc;
    std::vector<OpenSwath::SpectrumPtr> sp = sc.fetchSpectrumSwath(swath_ptr, 20.0, 1, im_range_empty);

    TEST_EQUAL(sp.empty(), true);
  }


  // test result for map with single spectrum
  {
    PeakMap* eptr = new PeakMap;
    MSSpectrum s;
    Peak1D p;
    p.setMZ(20.0);
    p.setIntensity(200.0);
    s.push_back(p);
    s.setRT(20.0);
    eptr->addSpectrum(s);
    boost::shared_ptr<PeakMap > swath_map (eptr);
    OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);

    TEST_EQUAL(swath_ptr->getNrSpectra(), 1)
    OpenSwathScoring sc;
    OpenSwath_Scores_Usage su;
    sc.initialize(1.0, 1, 0.005, 0.0, su, "resample", true);

    std::vector<OpenSwath::SpectrumPtr> sp = sc.fetchSpectrumSwath(swath_ptr, 20.0, 1, im_range_empty);

    TEST_EQUAL(sp.size(), 1);
    TEST_EQUAL(sp[0]->getMZArray()->data.size(), 1);
    TEST_EQUAL(sp[0]->getIntensityArray()->data.size(), 1);

    TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[0], 20.0);
    TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[0], 200.0);

    sc.initialize(1.0, 1, 0.005, 0.0, su, "simple", true);
    sp = sc.fetchSpectrumSwath(swath_ptr, 20.0, 1, im_range_empty);

    TEST_EQUAL(sp.size(), 1);
    TEST_EQUAL(sp[0]->getMZArray()->data.size(), 1);
    TEST_EQUAL(sp[0]->getIntensityArray()->data.size(), 1);

    TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[0], 20.0);
    TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[0], 200.0);
    // delete eptr;
  }

  // test result for map with three spectra
  {
    PeakMap* eptr = new PeakMap;
    MSSpectrum s;
    Peak1D p;
    p.setMZ(20.0);
    p.setIntensity(200.0);
    s.push_back(p);
    s.setRT(10.0);
    eptr->addSpectrum(s);
    s.setRT(20.0);
    eptr->addSpectrum(s);
    s.setRT(30.0);
    eptr->addSpectrum(s);
    boost::shared_ptr<PeakMap > swath_map (eptr);
    OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);

    TEST_EQUAL(swath_ptr->getNrSpectra(), 3)
    OpenSwathScoring sc;
    OpenSwath_Scores_Usage su;
    sc.initialize(1.0, 1, 0.005, 0.0, su, "resample", true);
    std::vector<OpenSwath::SpectrumPtr> sp = sc.fetchSpectrumSwath(swath_ptr, 20.0, 3, im_range_empty);

    TEST_EQUAL(sp.size(), 1);
    TEST_EQUAL(sp[0]->getMZArray()->data.size(), 1);
    TEST_EQUAL(sp[0]->getIntensityArray()->data.size(), 1);

    TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[0], 20.0);
    TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[0], 600.0);

    sc.initialize(1.0, 1, 0.005, 0.0, su, "simple", true);
    sp = sc.fetchSpectrumSwath(swath_ptr, 20.0, 3, im_range_empty);
    TEST_EQUAL(sp.size(), 3);

    TEST_EQUAL(sp[0]->getMZArray()->data.size(), 1);
    TEST_EQUAL(sp[0]->getIntensityArray()->data.size(), 1);


    TEST_EQUAL(sp[1]->getMZArray()->data.size(), 1);
    TEST_EQUAL(sp[1]->getIntensityArray()->data.size(), 1);


    TEST_EQUAL(sp[2]->getMZArray()->data.size(), 1);
    TEST_EQUAL(sp[2]->getIntensityArray()->data.size(), 1);

    TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[0], 20.0);
    TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[0], 200.0);
    TEST_REAL_SIMILAR(sp[1]->getMZArray()->data[0], 20.0);
    TEST_REAL_SIMILAR(sp[1]->getIntensityArray()->data[0], 200.0);
    TEST_REAL_SIMILAR(sp[2]->getMZArray()->data[0], 20.0);
    TEST_REAL_SIMILAR(sp[2]->getIntensityArray()->data[0], 200.0);
    // delete eptr;
  }

  // test result for map with uneven number of spectra
  {
    PeakMap* eptr = new PeakMap;
    {
      MSSpectrum s;
      s.emplace_back(20.0, 200.0);
      s.setRT(10.0);
      eptr->addSpectrum(s);
    }
    {
      MSSpectrum s;
      s.emplace_back(20.001, 200.0);
      s.setRT(20.0);
      eptr->addSpectrum(s);
    }
    {
      MSSpectrum s;
      s.emplace_back(250.001, 300.0);
      s.setRT(50.0);
      eptr->addSpectrum(s);
    }
    {
      MSSpectrum s;
      s.emplace_back(250.002, 500.0);
      s.setRT(60.0);
      eptr->addSpectrum(s);
    }
    boost::shared_ptr<PeakMap > swath_map (eptr);
    OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);

    TEST_EQUAL(swath_ptr->getNrSpectra(), 4)
    OpenSwathScoring sc;
    OpenSwath_Scores_Usage su;
    sc.initialize(1.0, 1, 0.005, 0.0, su, "resample", true);
    std::vector<OpenSwath::SpectrumPtr> sp = sc.fetchSpectrumSwath(swath_ptr, 20.0, 3, im_range_empty);

    TEST_EQUAL(sp.size(), 1);
    TEST_EQUAL(sp[0]->getMZArray()->data.size(), 3);
    TEST_EQUAL(sp[0]->getIntensityArray()->data.size(), 3);

    TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[0], 20.0);
    TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[0], 360.0);
    TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[1], 20.005);
    TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[1], 40.0);
    TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[2], 250.0);
    TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[2], 300.0);

    // in simple method all 3 spectra should be returned
    sc.initialize(1.0, 1, 0.005, 0.0, su, "simple", true);
    sp = sc.fetchSpectrumSwath(swath_ptr, 20.0, 3, im_range_empty);
    TEST_EQUAL(sp.size(), 3);
    TEST_EQUAL(sp[0]->getMZArray()->data.size(), 1);
    TEST_EQUAL(sp[1]->getMZArray()->data.size(), 1);
    TEST_EQUAL(sp[2]->getMZArray()->data.size(), 1);
    TEST_EQUAL(sp[0]->getIntensityArray()->data.size(), 1);
    TEST_EQUAL(sp[1]->getIntensityArray()->data.size(), 1);
    TEST_EQUAL(sp[2]->getIntensityArray()->data.size(), 1);
    TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[0], 20.001);
    TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[0], 200.0);
    TEST_REAL_SIMILAR(sp[1]->getMZArray()->data[0], 20.0);
    TEST_REAL_SIMILAR(sp[1]->getIntensityArray()->data[0], 200.0);
    TEST_REAL_SIMILAR(sp[2]->getMZArray()->data[0], 250.0);
    TEST_REAL_SIMILAR(sp[2]->getIntensityArray()->data[0], 300.0);
  }
}
END_SECTION

START_SECTION((OpenSwath::SpectrumPtr OpenSwathScoring::fetchSpectrumSwath(std::vector<OpenSwath::SwathMap> swath_maps,
                                                              double RT, int nr_spectra_to_add, RangeMobility im_range))- extra)
{

  // im range from 2-4
  OpenMS::RangeMobility im_range(3); // use this empty im range as input for all examples
  im_range.minSpanIfSingular(2); //


    // test result for map with single spectrum, should filter by IM because resampling is set
  {
    PeakMap* eptr = new PeakMap;
    eptr->addSpectrum(generateImSpec(1,6,20.0));
    boost::shared_ptr<PeakMap > swath_map (eptr);
    OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);
    TEST_EQUAL(swath_ptr->getNrSpectra(), 1);

    OpenSwathScoring sc;
    OpenSwath_Scores_Usage su;

    // test resample - IM filtering should occur
    {
      sc.initialize(1.0, 1, 0.005, 0.0, su, "resample", true);

      SpectrumSequence sp = sc.fetchSpectrumSwath(swath_ptr, 20.0, 1, im_range);

      TEST_EQUAL(sp.size(), 1);
      TEST_EQUAL(sp[0]->getMZArray()->data.size(), 3);
      TEST_EQUAL(sp[0]->getIntensityArray()->data.size(), 3);
      TEST_EQUAL(sp[0]->getDriftTimeArray()->data.size(), 3);

      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[0], 102.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[1], 103.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[2], 104.);

      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[0], 1.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[1], 1.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[2], 1.);

      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[0], 2.);
      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[1], 3.);
      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[2], 4.);
    }
    // test simple, since downstream functions are IM aware no filtering needs to occur.
    {
      sc.initialize(1.0, 1, 0.005, 0.0, su, "simple", true);
      SpectrumSequence sp = sc.fetchSpectrumSwath(swath_ptr, 20.0, 1, im_range);

      TEST_EQUAL(sp.size(), 1);
      TEST_EQUAL(sp[0]->getMZArray()->data.size(), 5);
      TEST_EQUAL(sp[0]->getIntensityArray()->data.size(), 5);

      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[0], 101.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[1], 102.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[2], 103.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[3], 104.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[4], 105.);

      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[0], 1.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[1], 1.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[2], 1.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[3], 1.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[4], 1.);

      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[0], 1.);
      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[1], 2.);
      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[2], 3.);
      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[3], 4.);
      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[4], 5.);
    }
    // delete eptr;
  }

  // Test result for 3 spectra
  {
    PeakMap* eptr = new PeakMap;
    eptr->addSpectrum(generateImSpec(1,3,19.0));
    eptr->addSpectrum(generateImSpec(1,6,20.0));
    eptr->addSpectrum(generateImSpec(3,6,21.0));
    boost::shared_ptr<PeakMap > swath_map (eptr);
    OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);
    TEST_EQUAL(swath_ptr->getNrSpectra(), 3);

    OpenSwathScoring sc;
    OpenSwath_Scores_Usage su;

    // test resample - IM filtering should occur, also IM information is not needed so is cleared
    {
      sc.initialize(1.0, 1, 0.005, 0.0, su, "resample", true);
      SpectrumSequence sp = sc.fetchSpectrumSwath(swath_ptr, 20.0, 3, im_range);

      TEST_EQUAL(sp.size(), 1);

      TEST_EQUAL(sp[0]->getMZArray()->data.size(), 3);
      TEST_EQUAL(sp[0]->getIntensityArray()->data.size(), 3);
      TEST_TRUE( (sp[0]->getDriftTimeArray() == nullptr ) ); // for resampling we do not use IM array

      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[0], 102.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[1], 103.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[2], 104.);

      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[0], 2.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[1], 2.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[2], 2.);
    }
    // test simple, since downstream functions are IM aware no filtering needs to occur. Should just return all the original spectra
    {
      sc.initialize(1.0, 1, 0.005, 0.0, su, "simple", true);
      SpectrumSequence sp = sc.fetchSpectrumSwath(swath_ptr, 20.0, 3, im_range);

      //test sizing
      TEST_EQUAL(sp.size(), 3);
      TEST_EQUAL(sp[0]->getMZArray()->data.size(), 5);
      TEST_EQUAL(sp[0]->getIntensityArray()->data.size(), 5);
      TEST_EQUAL(sp[0]->getDriftTimeArray()->data.size(), 5);

      TEST_EQUAL(sp[1]->getMZArray()->data.size(), 2);
      TEST_EQUAL(sp[1]->getIntensityArray()->data.size(), 2);
      TEST_EQUAL(sp[1]->getDriftTimeArray()->data.size(), 2);

      TEST_EQUAL(sp[2]->getMZArray()->data.size(), 3);
      TEST_EQUAL(sp[2]->getIntensityArray()->data.size(), 3);
      TEST_EQUAL(sp[2]->getDriftTimeArray()->data.size(), 3);

      // Spectrum #1
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[0], 101.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[1], 102.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[2], 103.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[3], 104.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[4], 105.);

      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[0], 1.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[1], 1.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[2], 1.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[3], 1.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[4], 1.);

      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[0], 1.);
      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[1], 2.);
      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[2], 3.);
      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[3], 4.);
      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[4], 5.);

      // Spectrum #2
      TEST_REAL_SIMILAR(sp[1]->getMZArray()->data[0], 101.);
      TEST_REAL_SIMILAR(sp[1]->getMZArray()->data[1], 102.);

      TEST_REAL_SIMILAR(sp[1]->getIntensityArray()->data[0], 1.);
      TEST_REAL_SIMILAR(sp[1]->getIntensityArray()->data[1], 1.);

      TEST_REAL_SIMILAR(sp[1]->getDriftTimeArray()->data[0], 1.);
      TEST_REAL_SIMILAR(sp[1]->getDriftTimeArray()->data[1], 2.);

      // Spectrum #3
      TEST_REAL_SIMILAR(sp[2]->getMZArray()->data[0], 103.);
      TEST_REAL_SIMILAR(sp[2]->getMZArray()->data[1], 104.);

      TEST_REAL_SIMILAR(sp[2]->getIntensityArray()->data[0], 1.);
      TEST_REAL_SIMILAR(sp[2]->getIntensityArray()->data[1], 1.);

      TEST_REAL_SIMILAR(sp[2]->getDriftTimeArray()->data[0], 3.);
      TEST_REAL_SIMILAR(sp[2]->getDriftTimeArray()->data[1], 4.);
    }
    // delete eptr;
  }

  // test result for map with 4 spectra (select 3)
  {
    PeakMap* eptr = new PeakMap;
    eptr->addSpectrum(generateImSpec(1,3,19.0));
    eptr->addSpectrum(generateImSpec(1,6,20.0));
    eptr->addSpectrum(generateImSpec(3,6,21.0));
    eptr->addSpectrum(generateImSpec(1,6,250));
    boost::shared_ptr<PeakMap > swath_map (eptr);
    OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);
    TEST_EQUAL(swath_ptr->getNrSpectra(), 4);

    OpenSwathScoring sc;
    OpenSwath_Scores_Usage su;

    // Test resampling, IM filtering should occur and the 4th spectrum should not be selected
    {
      sc.initialize(1.0, 1, 0.005, 0.0, su, "resample", true);

      SpectrumSequence sp = sc.fetchSpectrumSwath(swath_ptr, 20.0, 3, im_range);

      TEST_EQUAL(sp.size(), 1);
      TEST_EQUAL(sp[0]->getMZArray()->data.size(), 3);
      TEST_EQUAL(sp[0]->getIntensityArray()->data.size(), 3);
      TEST_TRUE( (sp[0]->getDriftTimeArray() == nullptr ) ); // for resampling we do not use IM array

      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[0], 102.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[1], 103.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[2], 104.);

      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[0], 2.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[1], 2.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[2], 2.);
    }

    // test simple, since downstream functions are IM aware no filtering needs to occur. Should just return all the original spectra, but the 4th spectrum should not be selected
    {
      sc.initialize(1.0, 1, 0.005, 0.0, su, "simple", true);
      SpectrumSequence sp = sc.fetchSpectrumSwath(swath_ptr, 20.0, 3, im_range);

      //test sizing
      TEST_EQUAL(sp.size(), 3);
      TEST_EQUAL(sp[0]->getMZArray()->data.size(), 5);
      TEST_EQUAL(sp[0]->getIntensityArray()->data.size(), 5);
      TEST_EQUAL(sp[0]->getDriftTimeArray()->data.size(), 5);

      TEST_EQUAL(sp[1]->getMZArray()->data.size(), 2);
      TEST_EQUAL(sp[1]->getIntensityArray()->data.size(), 2);
      TEST_EQUAL(sp[1]->getDriftTimeArray()->data.size(), 2);

      TEST_EQUAL(sp[2]->getMZArray()->data.size(), 3);
      TEST_EQUAL(sp[2]->getIntensityArray()->data.size(), 3);
      TEST_EQUAL(sp[2]->getDriftTimeArray()->data.size(), 3);

      // Spectrum #1
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[0], 101.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[1], 102.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[2], 103.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[3], 104.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[4], 105.);

      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[0], 1.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[1], 1.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[2], 1.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[3], 1.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[4], 1.);

      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[0], 1.);
      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[1], 2.);
      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[2], 3.);
      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[3], 4.);
      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[4], 5.);

      // Spectrum #2
      TEST_REAL_SIMILAR(sp[1]->getMZArray()->data[0], 101.);
      TEST_REAL_SIMILAR(sp[1]->getMZArray()->data[1], 102.);

      TEST_REAL_SIMILAR(sp[1]->getIntensityArray()->data[0], 1.);
      TEST_REAL_SIMILAR(sp[1]->getIntensityArray()->data[1], 1.);

      TEST_REAL_SIMILAR(sp[1]->getDriftTimeArray()->data[0], 1.);
      TEST_REAL_SIMILAR(sp[1]->getDriftTimeArray()->data[1], 2.);

      // Spectrum #3
      TEST_REAL_SIMILAR(sp[2]->getMZArray()->data[0], 103.);
      TEST_REAL_SIMILAR(sp[2]->getMZArray()->data[1], 104.);

      TEST_REAL_SIMILAR(sp[2]->getIntensityArray()->data[0], 1.);
      TEST_REAL_SIMILAR(sp[2]->getIntensityArray()->data[1], 1.);

      TEST_REAL_SIMILAR(sp[2]->getDriftTimeArray()->data[0], 3.);
      TEST_REAL_SIMILAR(sp[2]->getDriftTimeArray()->data[1], 4.);
    }
    // delete eptr;
  }
}
END_SECTION
END_TEST
