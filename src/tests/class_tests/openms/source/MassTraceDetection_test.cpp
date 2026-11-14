// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Erhan Kenar$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/test_config.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

///////////////////////////
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

/// provide access to private functions
class MassTraceDetectionAccess : public OpenMS::MassTraceDetection
{
public:
  using OpenMS::MassTraceDetection::updateIterativeWeightedMean_;
};

START_TEST(MassTraceDetection, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MassTraceDetection* ptr = nullptr;
MassTraceDetection* null_ptr = nullptr;
START_SECTION(MassTraceDetection())
{
    ptr = new MassTraceDetection();
    TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MassTraceDetection())
{
    delete ptr;
}
END_SECTION

MassTraceDetection test_mtd;

START_SECTION((void updateIterativeWeightedMean_(const double &, const double &, double &, double &, double &)))
{
    double centroid_mz(150.22), centroid_int(25000000);
    double new_mz1(150.34), new_int1(23043030);
    double new_mz2(150.11), new_int2(1932392);

    std::vector<double> mzs, ints;
    mzs.push_back(centroid_mz);
    mzs.push_back(new_mz1);
    mzs.push_back(new_mz2);
    ints.push_back(centroid_int);
    ints.push_back(new_int1);
    ints.push_back(new_int2);

    double total_weight1(centroid_int + new_int1);
    double total_weight2(centroid_int + new_int1 + new_int2);

    double wmean1((centroid_mz * centroid_int + new_mz1 * new_int1)/total_weight1);
    double wmean2((centroid_mz * centroid_int + new_mz1 * new_int1 + new_mz2 * new_int2)/total_weight2);

    double prev_count(centroid_mz * centroid_int);
    double prev_denom(centroid_int);

    MassTraceDetectionAccess::updateIterativeWeightedMean_(new_mz1, new_int1, centroid_mz, prev_count, prev_denom);


    TEST_REAL_SIMILAR(centroid_mz, wmean1);

    MassTraceDetectionAccess::updateIterativeWeightedMean_(new_mz2, new_int2, centroid_mz, prev_count, prev_denom);

    TEST_REAL_SIMILAR(centroid_mz, wmean2);

}
END_SECTION

// load a mzML file for testing the algorithm
PeakMap input;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MassTraceDetection_input1.mzML"),input);

// The seed cap excludes weak-isotope observations at +20.84 and +26.68 ppm.
Size exp_mt_lengths[3] = {86, 31, 14};
double exp_mt_rts[3] = {348.618754, 347.086557, 346.812732}; // centroid RTs should be reasonably similar (isotopic traces)
// note: these shifted slightly earlier when updateWeightedMeanRT() started to weight the first
// peak of a trace as well (it used to be skipped, biasing the centroid late -- see issue #2777)
double exp_mt_mzs[3] = {437.26675, 438.27241, 439.275271947714};
double exp_mt_ints[3] = {3381.72226139326, 664.763828332733, 111.367062259816};

std::vector<MassTrace> output_mt;

Param p_mtd = MassTraceDetection().getDefaults();
p_mtd.setValue("min_trace_length", 3.0);

START_SECTION((void run(const PeakMap &, std::vector< MassTrace > &)))
{
    test_mtd.run(input, output_mt);

    // with default parameters, only 2 of 3 traces will be found
    TEST_EQUAL(output_mt.size(), 2);

    // if min_trace_length is set to 3 seconds, another mass trace is detected
    test_mtd.setParameters(p_mtd);
    output_mt.clear();

    test_mtd.run(input, output_mt);

    TEST_EQUAL(output_mt.size(), 3);

    for (Size i = 0; i < output_mt.size(); ++i)
    {
        TEST_EQUAL(output_mt[i].getSize(), exp_mt_lengths[i]);
        TEST_REAL_SIMILAR(output_mt[i].getCentroidRT(), exp_mt_rts[i]);
        TEST_REAL_SIMILAR(output_mt[i].getCentroidMZ(), exp_mt_mzs[i]);
        TEST_REAL_SIMILAR(output_mt[i].computePeakArea(), exp_mt_ints[i]);
    }

    // Regression test for bug #1633
    // Test by adding MS2 spectra to the input
    {
      PeakMap input_new;
      MSSpectrum s;
      s.setMSLevel(2);
      {
        Peak1D p;
        p.setMZ( 500 );
        p.setIntensity( 6000 );
        s.push_back(p);
      }

      // add a few additional MS2 spectra in front
      for (Size i = 0; i < input.size(); ++i)
      {
        input_new.addSpectrum(s);
      }
      // now add the "real" spectra at the end
      for (Size i = 0; i < input.size(); ++i)
      {
        input_new.addSpectrum(input[i]);
      }
      output_mt.clear();
      test_mtd.run(input_new, output_mt);
      TEST_EQUAL(output_mt.size(), 3);

      for (Size i = 0; i < output_mt.size(); ++i)
      {
          TEST_EQUAL(output_mt[i].getSize(), exp_mt_lengths[i]);
          TEST_REAL_SIMILAR(output_mt[i].getCentroidRT(), exp_mt_rts[i]);
          TEST_REAL_SIMILAR(output_mt[i].getCentroidMZ(), exp_mt_mzs[i]);
          TEST_REAL_SIMILAR(output_mt[i].computePeakArea(), exp_mt_ints[i]);
      }

    }
}
END_SECTION

std::vector<MassTrace> filt;

//START_SECTION((void filterByPeakWidth(std::vector< MassTrace > &, std::vector< MassTrace > &)))
//{
//    test_mtd.filterByPeakWidth(output_mt, filt);

//    TEST_EQUAL(output_mt.size(), filt.size());

////    for (Size i = 0; i < output_mt.size(); ++i)
////    {
////        TEST_EQUAL(output_mt[i].getFWHMScansNum(), filt[i].getFWHMScansNum());
////    }
//}
//END_SECTION

PeakMap::ConstAreaIterator mt_it1 = input.areaBeginConst(335.0, 385.0, 437.1, 437.4);
PeakMap::ConstAreaIterator mt_it2 = input.areaBeginConst(335.0, 385.0, 438.2, 438.4);
PeakMap::ConstAreaIterator mt_it3 = input.areaBeginConst(335.0, 385.0, 439.2, 439.4);

std::vector<MassTrace> found_mtraces;

PeakMap::ConstAreaIterator mt_end = input.areaEndConst();

START_SECTION((void run(PeakMap::ConstAreaIterator &begin, PeakMap::ConstAreaIterator &end, std::vector< MassTrace > &found_masstraces)))
{

    NOT_TESTABLE
//    test_mtd.run(mt_it1, mt_end, found_mtraces);
//    TEST_EQUAL(found_mtraces.size(), 1);
//    TEST_EQUAL(found_mtraces[0].getSize(), exp_mt_lengths[0]);

//    TEST_REAL_SIMILAR(found_mtraces[0].getCentroidRT(), exp_mt_rts[0]);
//    TEST_REAL_SIMILAR(found_mtraces[0].getCentroidMZ(), exp_mt_mzs[0]);
//    TEST_REAL_SIMILAR(found_mtraces[0].computePeakArea(), exp_mt_ints[0]);

//    found_mtraces.clear();


//    test_mtd.run(mt_it2, mt_end, found_mtraces);
//    TEST_EQUAL(found_mtraces.size(), 1);
//    TEST_EQUAL(found_mtraces[0].getSize(), exp_mt_lengths[1]);

//    TEST_REAL_SIMILAR(found_mtraces[0].getCentroidRT(), exp_mt_rts[1]);
//    TEST_REAL_SIMILAR(found_mtraces[0].getCentroidMZ(), exp_mt_mzs[1]);
//    TEST_REAL_SIMILAR(found_mtraces[0].computePeakArea(), exp_mt_ints[1]);

//    found_mtraces.clear();


//    test_mtd.run(mt_it3, mt_end, found_mtraces);
//    TEST_EQUAL(found_mtraces.size(), 1);
//    TEST_EQUAL(found_mtraces[0].getSize(), exp_mt_lengths[0]);

//    TEST_REAL_SIMILAR(found_mtraces[0].getCentroidRT(), exp_mt_rts[2]);
//    TEST_REAL_SIMILAR(found_mtraces[0].getCentroidMZ(), exp_mt_mzs[2]);
//    TEST_REAL_SIMILAR(found_mtraces[0].computePeakArea(), exp_mt_ints[2]);

//    found_mtraces.clear();
}
END_SECTION

START_SECTION((fixed apex cap prevents adaptive mass drift))
{
  PeakMap drifting;
  for (Size scan = 0; scan < 16; ++scan)
  {
    MSSpectrum spectrum;
    spectrum.setRT(scan);
    spectrum.emplace_back(500.0 + 0.008 * scan, 1000.0 - 10.0 * scan);
    drifting.addSpectrum(spectrum);
  }
  MassTraceDetection detector;
  detector.setLogType(ProgressLogger::NONE);
  Param parameters = detector.getParameters();
  parameters.setValue("min_trace_length", 1.0);
  detector.setParameters(parameters);
  vector<MassTrace> capped, legacy;
  detector.run(drifting, capped);
  TEST_FALSE(capped.empty())
  if (! capped.empty())
  {
    TEST_EQUAL(capped.front().getSize(), 2)
    for (const auto& peak : capped.front())
    {
      TEST_TRUE(peak.getMZ() >= 499.99 && peak.getMZ() <= 500.01)
    }
  }
  parameters.setValue("max_apex_deviation_ppm", 0.0);
  detector.setParameters(parameters);
  detector.run(drifting, legacy);
  TEST_FALSE(legacy.empty())
  if (! legacy.empty())
  {
    TEST_TRUE(legacy.front().getSize() > 2)
    TEST_TRUE(any_of(legacy.front().begin(), legacy.front().end(), [](const auto& peak) { return peak.getMZ() > 500.01; }))
  }
}
END_SECTION

START_SECTION((disabled apex cap reproduces legacy fixture))
{
  MassTraceDetection legacy;
  legacy.setLogType(ProgressLogger::NONE);
  Param parameters = legacy.getParameters();
  parameters.setValue("min_trace_length", 3.0);
  parameters.setValue("max_apex_deviation_ppm", 0.0);
  legacy.setParameters(parameters);
  vector<MassTrace> traces;
  legacy.run(input, traces);
  TEST_EQUAL(traces.size(), 3)
  if (traces.size() == 3)
  {
    TEST_EQUAL(traces[2].getSize(), 16)
    TEST_REAL_SIMILAR(traces[2].getCentroidRT(), 346.81803158056221)
    TEST_REAL_SIMILAR(traces[2].computePeakArea(), 109.49010862067628)
  }
}
END_SECTION

START_SECTION((candidate search excludes peaks outside cap before nearest selection))
{
  PeakMap experiment;
  for (Size scan = 0; scan < 5; ++scan)
  {
    MSSpectrum spectrum;
    spectrum.setRT(scan);
    if (scan == 2) { spectrum.emplace_back(499.9899, 20.0); }
    spectrum.emplace_back(scan == 0 ? 500.0 : scan == 1 ? 499.99 : 500.004, 1000.0 - 100.0 * scan);
    experiment.addSpectrum(spectrum);
  }
  MassTraceDetection detector;
  detector.setLogType(ProgressLogger::NONE);
  Param parameters = detector.getParameters();
  parameters.setValue("min_trace_length", 3.0);
  detector.setParameters(parameters);
  vector<MassTrace> traces;
  detector.run(experiment, traces);
  TEST_FALSE(traces.empty())
  if (! traces.empty())
  {
    TEST_EQUAL(traces.front().getSize(), 5)
    if (traces.front().getSize() == 5) { TEST_EQUAL(traces.front()[2].getMZ(), 500.004) }
  }

  PeakMap endpoints;
  for (Size scan = 0; scan < 5; ++scan)
  {
    MSSpectrum spectrum;
    spectrum.setRT(scan);
    spectrum.emplace_back(scan == 2 ? 500.0 : scan % 2 == 0 ? 499.99 : 500.01, scan == 2 ? 1000.0 : 900.0);
    endpoints.addSpectrum(spectrum);
  }
  detector.run(endpoints, traces);
  TEST_EQUAL(traces.size(), 1)
  if (! traces.empty()) { TEST_EQUAL(traces.front().getSize(), 5) }
}
END_SECTION

START_SECTION((component extraction preserves serial peak ownership, order and metadata))
{
  PeakMap experiment;
  for (Size scan = 0; scan <= 16; ++scan)
  {
    MSSpectrum spectrum;
    spectrum.setRT(scan);
    MSSpectrum::FloatDataArray mz_width, im, im_width;
    mz_width.setName(Constants::UserParam::FWHM_MZ_ppm);
    im.setName(Constants::UserParam::ION_MOBILITY);
    im_width.setName(Constants::UserParam::FWHM_IM);
    if (scan % 4 == 0)
    {
      for (Size component = 0; component < 32; ++component)
      {
        spectrum.emplace_back(100.0 + 5.0 * component, 1000.0);
        mz_width.push_back(2.0 + component);
        im.push_back(1.0);
        im_width.push_back(0.005);
      }
      spectrum.emplace_back(scan < 8 ? 500.0 : scan == 8 ? 500.009 : 500.018, scan == 8 ? 20.0 : scan % 12 == 0 ? 1000.0 : 950.0);
      mz_width.push_back(3.0);
      im.push_back(1.0);
      im_width.push_back(0.006);
    }
    spectrum.getFloatDataArrays() = {mz_width, im, im_width};
    experiment.addSpectrum(spectrum);
  }

  auto compare = [&](const vector<MassTrace>& serial, const vector<MassTrace>& parallel) {
    TEST_EQUAL(serial.size(), parallel.size())
    for (Size trace = 0; trace < min(serial.size(), parallel.size()); ++trace)
    {
      TEST_EQUAL(serial[trace].getLabel(), parallel[trace].getLabel())
      TEST_EQUAL(serial[trace].getSize(), parallel[trace].getSize())
      TEST_EQUAL(serial[trace].getCentroidMZ(), parallel[trace].getCentroidMZ())
      TEST_EQUAL(serial[trace].getCentroidRT(), parallel[trace].getCentroidRT())
      TEST_EQUAL(serial[trace].getCentroidIM(), parallel[trace].getCentroidIM())
      TEST_EQUAL(serial[trace].getCentroidSD(), parallel[trace].getCentroidSD())
      TEST_EQUAL(serial[trace].computePeakArea(), parallel[trace].computePeakArea())
      TEST_EQUAL(serial[trace].fwhm_mz_avg, parallel[trace].fwhm_mz_avg)
      TEST_EQUAL(serial[trace].fwhm_im_avg, parallel[trace].fwhm_im_avg)
      for (Size peak = 0; peak < min(serial[trace].getSize(), parallel[trace].getSize()); ++peak)
      {
        TEST_TRUE(serial[trace][peak] == parallel[trace][peak])
      }
    }
  };
#ifdef _OPENMP
  const int previous_threads = omp_get_max_threads();
#endif
  for (const string& criterion : {string("outlier"), string("sample_rate")})
  {
    for (double cap : {0.0, 20.0})
    {
      for (int outliers : {0, 1, 5})
      {
        MassTraceDetection detector;
        detector.setLogType(ProgressLogger::NONE);
        Param parameters = detector.getParameters();
        parameters.setValue("max_apex_deviation_ppm", cap);
        parameters.setValue("trace_termination_criterion", criterion);
        parameters.setValue("trace_termination_outliers", outliers);
        parameters.setValue("min_trace_length", 3.0);
        parameters.setValue("min_sample_rate", 0.1);
        detector.setParameters(parameters);
        for (Size maximum : {Size(0), Size(1), Size(2)})
        {
          vector<MassTrace> serial, parallel;
#ifdef _OPENMP
          omp_set_num_threads(1);
#endif
          detector.run(experiment, serial, maximum);
#ifdef _OPENMP
          omp_set_num_threads(2);
#endif
          detector.run(experiment, parallel, maximum);
          compare(serial, parallel);
          TEST_TRUE(detector.hasCentroidIm())
          TEST_TRUE(detector.hasFwhmMz())
          TEST_TRUE(detector.hasFwhmIm())
        }
      }
    }
  }
#ifdef _OPENMP
  omp_set_num_threads(previous_threads);
#endif

  MassTraceDetection reused;
  reused.setLogType(ProgressLogger::NONE);
  vector<MassTrace> traces;
  reused.run(experiment, traces);
  for (auto& spectrum : experiment)
  {
    spectrum.getFloatDataArrays().clear();
  }
  reused.run(experiment, traces);
  TEST_FALSE(reused.hasCentroidIm())
  TEST_FALSE(reused.hasFwhmMz())
  TEST_FALSE(reused.hasFwhmIm())
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
