// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------
#include <memory>

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
#include <OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/MockObjects.h>

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



// names of the scores that differ between @p a and @p b
std::string differingScores(const OpenSwath_Scores& a, const OpenSwath_Scores& b)
{
  std::string differing;
#define CHECK_SCORE(field) if (!(a.field == b.field)) differing += std::string(#field) + " ";
  CHECK_SCORE(elution_model_fit_score);
  CHECK_SCORE(library_corr);
  CHECK_SCORE(library_norm_manhattan);
  CHECK_SCORE(library_rootmeansquare);
  CHECK_SCORE(library_sangle);
  CHECK_SCORE(norm_rt_score);
  CHECK_SCORE(isotope_correlation);
  CHECK_SCORE(isotope_overlap);
  CHECK_SCORE(massdev_score);
  CHECK_SCORE(xcorr_coelution_score);
  CHECK_SCORE(xcorr_shape_score);
  CHECK_SCORE(yseries_score);
  CHECK_SCORE(bseries_score);
  CHECK_SCORE(log_sn_score);
  CHECK_SCORE(weighted_coelution_score);
  CHECK_SCORE(weighted_xcorr_shape);
  CHECK_SCORE(weighted_massdev_score);
  CHECK_SCORE(ms1_xcorr_coelution_score);
  CHECK_SCORE(ms1_xcorr_coelution_contrast_score);
  CHECK_SCORE(ms1_xcorr_coelution_combined_score);
  CHECK_SCORE(ms1_xcorr_shape_score);
  CHECK_SCORE(ms1_xcorr_shape_contrast_score);
  CHECK_SCORE(ms1_xcorr_shape_combined_score);
  CHECK_SCORE(ms1_ppm_score);
  CHECK_SCORE(ms1_isotope_correlation);
  CHECK_SCORE(ms1_isotope_overlap);
  CHECK_SCORE(ms1_mi_score);
  CHECK_SCORE(ms1_mi_contrast_score);
  CHECK_SCORE(ms1_mi_combined_score);
  CHECK_SCORE(im_xcorr_coelution_score);
  CHECK_SCORE(im_xcorr_shape_score);
  CHECK_SCORE(im_delta_score);
  CHECK_SCORE(im_ms1_delta_score);
  CHECK_SCORE(im_drift);
  CHECK_SCORE(im_drift_left);
  CHECK_SCORE(im_drift_right);
  CHECK_SCORE(im_drift_weighted);
  CHECK_SCORE(im_delta);
  CHECK_SCORE(im_log_intensity);
  CHECK_SCORE(im_ms1_contrast_coelution);
  CHECK_SCORE(im_ms1_contrast_shape);
  CHECK_SCORE(im_ms1_sum_contrast_coelution);
  CHECK_SCORE(im_ms1_sum_contrast_shape);
  CHECK_SCORE(im_ms1_drift);
  CHECK_SCORE(im_ms1_delta);
  CHECK_SCORE(im_ind_contrast_coelution);
  CHECK_SCORE(im_ind_contrast_shape);
  CHECK_SCORE(im_ind_sum_contrast_coelution);
  CHECK_SCORE(im_ind_sum_contrast_shape);
  CHECK_SCORE(library_manhattan);
  CHECK_SCORE(library_dotprod);
  CHECK_SCORE(intensity);
  CHECK_SCORE(total_xic);
  CHECK_SCORE(nr_peaks);
  CHECK_SCORE(sn_ratio);
  CHECK_SCORE(mi_score);
  CHECK_SCORE(weighted_mi_score);
  CHECK_SCORE(rt_difference);
  CHECK_SCORE(normalized_experimental_rt);
  CHECK_SCORE(raw_rt_score);
  CHECK_SCORE(dotprod_score_dia);
  CHECK_SCORE(manhatt_score_dia);
#undef CHECK_SCORE
  return differing;
}

MSSpectrum peaksSpectrum(double rt, const std::vector<std::pair<double, double>>& peaks)
{
  MSSpectrum spectrum;
  spectrum.setRT(rt);
  for (const auto& peak : peaks)
  {
    spectrum.push_back(Peak1D(peak.first, peak.second));
  }
  return spectrum;
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
  ptr->initialize(100.0, 1, 0.01, 0.20, 0.0, su, "simple", "fixed", true, true);
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

START_SECTION((OpenSwath::SpectrumPtr OpenSwathScoring::fetchSpectrumSwath(std::vector<OpenSwath::SwathMap> swath_maps,
                                                              double RT, int nr_spectra_to_add, RangeMobility im_range)))
{

  OpenMS::RangeMobility im_range_empty; // use this empty im range as input for all examples
  // test result for empty map
  {
    std::shared_ptr<PeakMap > swath_map (new PeakMap);
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
    std::shared_ptr<PeakMap > swath_map (eptr);
    OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);

    TEST_EQUAL(swath_ptr->getNrSpectra(), 1)
    OpenSwathScoring sc;
    OpenSwath_Scores_Usage su;
    sc.initialize(1.0, 1, 0.005, 0.20, 0.0, su, "resample", "fixed", true, true);

    std::vector<OpenSwath::SpectrumPtr> sp = sc.fetchSpectrumSwath(swath_ptr, 20.0, 1, im_range_empty);

    TEST_EQUAL(sp.size(), 1);
    TEST_EQUAL(sp[0]->getMZArray()->data.size(), 1);
    TEST_EQUAL(sp[0]->getIntensityArray()->data.size(), 1);

    TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[0], 20.0);
    TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[0], 200.0);

    sc.initialize(1.0, 1, 0.005, 0.20, 0.0, su, "simple", "fixed", true, true);
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
    std::shared_ptr<PeakMap > swath_map (eptr);
    OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);

    TEST_EQUAL(swath_ptr->getNrSpectra(), 3)
    OpenSwathScoring sc;
    OpenSwath_Scores_Usage su;
    sc.initialize(1.0, 1, 0.005, 0.20, 0.0, su, "resample", "fixed", true, true);
    std::vector<OpenSwath::SpectrumPtr> sp = sc.fetchSpectrumSwath(swath_ptr, 20.0, 3, im_range_empty);

    TEST_EQUAL(sp.size(), 1);
    TEST_EQUAL(sp[0]->getMZArray()->data.size(), 1);
    TEST_EQUAL(sp[0]->getIntensityArray()->data.size(), 1);

    TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[0], 20.0);
    TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[0], 600.0);

    sc.initialize(1.0, 1, 0.005, 0.20, 0.0, su, "simple", "fixed", true, true);
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
    std::shared_ptr<PeakMap > swath_map (eptr);
    OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);

    TEST_EQUAL(swath_ptr->getNrSpectra(), 4)
    OpenSwathScoring sc;
    OpenSwath_Scores_Usage su;
    sc.initialize(1.0, 1, 0.005, 0.20, 0.0, su, "resample", "fixed", true, true);
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
    sc.initialize(1.0, 1, 0.005, 0.20, 0.0, su, "simple", "fixed", true, true);
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
    std::shared_ptr<PeakMap > swath_map (eptr);
    OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);
    TEST_EQUAL(swath_ptr->getNrSpectra(), 1);

    OpenSwathScoring sc;
    OpenSwath_Scores_Usage su;

    // test resample - IM filtering should occur
    {
      sc.initialize(1.0, 1, 0.005, 0.20, 0.0, su, "resample", "fixed", true, true);

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
      sc.initialize(1.0, 1, 0.005, 0.20, 0.0, su, "simple", "fixed", true, true);
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
    std::shared_ptr<PeakMap > swath_map (eptr);
    OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);
    TEST_EQUAL(swath_ptr->getNrSpectra(), 3);

    OpenSwathScoring sc;
    OpenSwath_Scores_Usage su;

    // test resample - IM filtering should occur and IM information is preserved
    {
      sc.initialize(1.0, 1, 0.005, 0.20, 0.0, su, "resample", "fixed", true, true);
      SpectrumSequence sp = sc.fetchSpectrumSwath(swath_ptr, 20.0, 3, im_range);

      TEST_EQUAL(sp.size(), 1);

      TEST_EQUAL(sp[0]->getMZArray()->data.size(), 6);
      TEST_EQUAL(sp[0]->getIntensityArray()->data.size(), 6);
      TEST_EQUAL(sp[0]->getDriftTimeArray()->data.size(), 6);

      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[0], 102.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[1], 102.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[2], 103.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[3], 103.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[4], 104.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[5], 104.);

      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[0], 1.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[1], 1.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[2], 1.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[3], 1.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[4], 1.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[5], 1.);

      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[0], 2.);
      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[1], 2.);
      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[2], 3.);
      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[3], 3.);
      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[4], 4.);
      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[5], 4.);
    }
    // test simple, since downstream functions are IM aware no filtering needs to occur. Should just return all the original spectra
    {
      sc.initialize(1.0, 1, 0.005, 0.20, 0.0, su, "simple", "fixed", true, true);
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

  // Test IM-enhanced input without an extraction range. Resampling is not
  // applied to IM data; spectra are concatenated to keep drift-time arrays.
  {
    OpenMS::RangeMobility im_range_empty;
    PeakMap* eptr = new PeakMap;
    eptr->addSpectrum(generateImSpec(1,3,19.0));
    eptr->addSpectrum(generateImSpec(1,6,20.0));
    eptr->addSpectrum(generateImSpec(3,6,21.0));
    std::shared_ptr<PeakMap > swath_map (eptr);
    OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);
    TEST_EQUAL(swath_ptr->getNrSpectra(), 3);

    OpenSwathScoring sc;
    OpenSwath_Scores_Usage su;
    sc.initialize(1.0, 1, 0.005, 0.20, 0.0, su, "resample", "fixed", true, true);

    SpectrumSequence sp = sc.fetchSpectrumSwath(swath_ptr, 20.0, 3, im_range_empty);

    TEST_EQUAL(sp.size(), 1);
    TEST_EQUAL(sp[0]->getMZArray()->data.size(), 10);
    TEST_EQUAL(sp[0]->getIntensityArray()->data.size(), 10);
    TEST_EQUAL(sp[0]->getDriftTimeArray()->data.size(), 10);

    TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[0], 101.);
    TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[1], 101.);
    TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[8], 105.);
    TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[9], 105.);
  }

  // test result for map with 4 spectra (select 3)
  {
    PeakMap* eptr = new PeakMap;
    eptr->addSpectrum(generateImSpec(1,3,19.0));
    eptr->addSpectrum(generateImSpec(1,6,20.0));
    eptr->addSpectrum(generateImSpec(3,6,21.0));
    eptr->addSpectrum(generateImSpec(1,6,250));
    std::shared_ptr<PeakMap > swath_map (eptr);
    OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);
    TEST_EQUAL(swath_ptr->getNrSpectra(), 4);

    OpenSwathScoring sc;
    OpenSwath_Scores_Usage su;

    // Test resampling, IM filtering should occur and the 4th spectrum should not be selected
    {
      sc.initialize(1.0, 1, 0.005, 0.20, 0.0, su, "resample", "fixed", true, true);

      SpectrumSequence sp = sc.fetchSpectrumSwath(swath_ptr, 20.0, 3, im_range);

      TEST_EQUAL(sp.size(), 1);
      TEST_EQUAL(sp[0]->getMZArray()->data.size(), 6);
      TEST_EQUAL(sp[0]->getIntensityArray()->data.size(), 6);
      TEST_EQUAL(sp[0]->getDriftTimeArray()->data.size(), 6);

      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[0], 102.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[1], 102.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[2], 103.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[3], 103.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[4], 104.);
      TEST_REAL_SIMILAR(sp[0]->getMZArray()->data[5], 104.);

      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[0], 1.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[1], 1.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[2], 1.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[3], 1.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[4], 1.);
      TEST_REAL_SIMILAR(sp[0]->getIntensityArray()->data[5], 1.);

      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[0], 2.);
      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[1], 2.);
      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[2], 3.);
      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[3], 3.);
      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[4], 4.);
      TEST_REAL_SIMILAR(sp[0]->getDriftTimeArray()->data[5], 4.);
    }

    // test simple, since downstream functions are IM aware no filtering needs to occur. Should just return all the original spectra, but the 4th spectrum should not be selected
    {
      sc.initialize(1.0, 1, 0.005, 0.20, 0.0, su, "simple", "fixed", true, true);
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

START_SECTION([EXTRA] ScoreSubset: RANKING followed by REMAINDER gives the scores of ALL)
{
  // a peak group of two transitions (precursor m/z 500, charge 2) with two precursor isotope traces
  std::vector<double> intensity1 {5.97543668746948, 4.2749171257019, 3.3301842212677, 4.08597040176392, 5.50307035446167, 5.24326848983765,
                                  8.40812492370605, 2.83419919013977, 6.94378805160522, 7.69957494735718, 4.08597040176392};
  std::vector<double> intensity2 {15.8951349258423, 41.5446395874023, 76.0746307373047, 109.069435119629, 111.90364074707, 169.79216003418,
                                  121.043930053711, 63.0136985778809, 44.6150207519531, 21.4926776885986, 7.93575811386108};
  std::vector<double> ms1intensity1 {0.0, 110.0, 200.0, 270.0, 320.0, 350.0, 360.0, 350.0, 320.0, 270.0, 200.0};
  std::vector<double> ms1intensity2 {10.0, 115.0, 180.0, 280.0, 330.0, 340.0, 390.0, 320.0, 300.0, 250.0, 100.0};
  auto trace = [](const std::vector<double>& intensity)
  {
    std::shared_ptr<OpenSwath::MockFeature> feature(new OpenSwath::MockFeature());
    feature->m_intensity_vec = intensity;
    return feature;
  };
  OpenSwath::MockMRMFeature feature;
  feature.m_features["tr1"] = trace(intensity1);
  feature.m_features["tr2"] = trace(intensity2);
  feature.m_precursor_features["prec_0"] = trace(ms1intensity1);
  feature.m_precursor_features["prec_1"] = trace(ms1intensity2);
  feature.m_rt = 20.0;
  feature.m_intensity = 1000.0;
  feature.m_metavalue = 0.0;
  const std::vector<std::string> native_ids {"tr1", "tr2"};
  const std::vector<std::string> precursor_ids {"prec_0", "prec_1"};
  const std::vector<double> library_intensity {0.3, 0.7};

  std::vector<OpenSwath::LightTransition> transitions(2);
  transitions[0].transition_name = "tr1";
  transitions[0].product_mz = 400.2;
  transitions[0].library_intensity = 30;
  transitions[1].transition_name = "tr2";
  transitions[1].product_mz = 600.3;
  transitions[1].library_intensity = 70;
  for (OpenSwath::LightTransition& transition : transitions)
  {
    transition.peptide_ref = "pep";
    transition.precursor_mz = 500.0;
    transition.fragment_charge = 1;
  }
  OpenSwath::LightCompound compound;
  compound.id = "pep";
  compound.sequence = "PEPTIDEK";
  compound.charge = 2;

  // fragments and precursor with their isotopes
  std::shared_ptr<PeakMap> swath_map(new PeakMap);
  std::shared_ptr<PeakMap> ms1_map(new PeakMap);
  for (double rt : {19.0, 20.0, 21.0})
  {
    swath_map->addSpectrum(peaksSpectrum(rt, {{400.2, 1000.0}, {401.2034, 450.0}, {402.2067, 120.0},
                                              {600.3, 2000.0}, {601.3034, 900.0}, {602.3067, 300.0}}));
    ms1_map->addSpectrum(peaksSpectrum(rt, {{500.0, 1000.0}, {500.5017, 600.0}, {501.0034, 250.0}}));
  }
  std::vector<OpenSwath::SwathMap> swath_maps(1);
  swath_maps[0].sptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);
  OpenSwath::SpectrumAccessPtr ms1_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(ms1_map);

  std::vector<OpenSwath::ISignalToNoisePtr> signal_noise;
  for (Size i = 0; i < native_ids.size(); ++i)
  {
    std::shared_ptr<OpenSwath::MockSignalToNoise> sn(new OpenSwath::MockSignalToNoise());
    sn->m_sn_value = 10.0 + i;
    signal_noise.push_back(sn);
  }

  // the scores of MRMFeatureFinderScoring::scorePeakgroups() when it defers the unranked scores
  OpenSwath_Scores_Usage su;
  su.use_im_scores = false;
  su.use_uis_scores = false;
  OpenSwathScoring sc;
  sc.initialize(1.0, 1, 0.005, 0.20, 0.0, su, "simple", "fixed", false, false);
  DIAScoring diascoring;
  const OpenMS::RangeMobility im_range;

  OpenSwath_Scores all;
  std::vector<double> masserror_all;
  sc.calculateChromatographicScores(&feature, native_ids, precursor_ids, library_intensity, signal_noise, all);
  sc.calculateDIAScores(&feature, transitions, library_intensity, swath_maps, ms1_ptr, diascoring, compound, all, masserror_all,
                        -1, im_range, nullptr, -1);

  OpenSwath_Scores split;
  std::vector<double> masserror_ranking;
  sc.calculateChromatographicScores(&feature, native_ids, precursor_ids, library_intensity, signal_noise, split,
                                    OpenSwathScoring::ScoreSubset::RANKING);
  sc.calculateDIAScores(&feature, transitions, library_intensity, swath_maps, ms1_ptr, diascoring, compound, split, masserror_ranking,
                        -1, im_range, nullptr, -1, OpenSwathScoring::ScoreSubset::RANKING);
  const OpenSwath_Scores ranking = split;
  std::vector<double> masserror_remainder;
  sc.calculateChromatographicScores(&feature, native_ids, precursor_ids, library_intensity, signal_noise, split,
                                    OpenSwathScoring::ScoreSubset::REMAINDER);
  sc.calculateDIAScores(&feature, transitions, library_intensity, swath_maps, ms1_ptr, diascoring, compound, split, masserror_remainder,
                        -1, im_range, nullptr, -1, OpenSwathScoring::ScoreSubset::REMAINDER);

  TEST_EQUAL(differingScores(split, all), "")
  TEST_EQUAL(masserror_ranking == masserror_all, true)
  // the ranking scores already give the prescore of all scores
  TEST_EQUAL(ranking.calculate_swath_lda_prescore(ranking), all.calculate_swath_lda_prescore(all))

  // the test covers the scores of both subsets
  TEST_NOT_EQUAL(all.xcorr_shape_score, 0)
  TEST_NOT_EQUAL(all.isotope_correlation, 0)
  TEST_NOT_EQUAL(all.mi_score, 0)
  TEST_NOT_EQUAL(all.ms1_xcorr_shape_score, -1)
  TEST_NOT_EQUAL(all.ms1_mi_score, -1)
  TEST_NOT_EQUAL(all.dotprod_score_dia, 0)
  TEST_NOT_EQUAL(all.ms1_isotope_correlation, 0)
  TEST_EQUAL(ranking.mi_score, 0)
  TEST_EQUAL(ranking.ms1_xcorr_shape_score, -1)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
