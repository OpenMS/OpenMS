// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Erhan Kenar, Prachi Agrawal$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/MzMLFile.h>

///////////////////////////
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
///////////////////////////

#include <algorithm>
#include <cmath>

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

Size exp_mt_lengths[3] = {86, 31, 16};
double exp_mt_rts[3] = {348.618754, 347.086557, 346.818032}; // centroid RTs should be reasonably similar (isotopic traces)
// note: these shifted slightly earlier when updateWeightedMeanRT() started to weight the first
// peak of a trace as well (it used to be skipped, biasing the centroid late -- see issue #2777)
double exp_mt_mzs[3] = {437.26675, 438.27241, 439.27594};
double exp_mt_ints[3] = {3381.72226139326, 664.763828332733, 109.490108620676};

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

START_SECTION((run() characterizes the mass trace drift reported in issue 10905))
{
  // This records known buggy behavior, not a chemically validated ion channel.
  // A future drift fix should replace these wide-span expectations.
  PeakMap drift_input;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MassTraceDetection_drift_10905.mzML.gz"), drift_input);

  MassTraceDetection detector;
  Param parameters = detector.getDefaults();
  parameters.setValue("mass_error_ppm", 20.0);
  parameters.setValue("noise_threshold_int", 10.0);
  parameters.setValue("chrom_peak_snr", 3.0);
  parameters.setValue("reestimate_mt_sd", "true");
  parameters.setValue("quant_method", "area");
  parameters.setValue("trace_termination_criterion", "outlier");
  parameters.setValue("trace_termination_outliers", 5);
  parameters.setValue("min_sample_rate", 0.5);
  parameters.setValue("min_trace_length", 5.0);
  parameters.setValue("max_trace_length", -1.0);
  detector.setParameters(parameters);

  std::vector<MassTrace> drift_traces;
  detector.run(drift_input, drift_traces);

  // Cropping changes extraction order and labels (T1930 becomes T14).
  // Identify the trace by its most intense original observation instead.
  const MassTrace* drift_trace = nullptr;
  Size matching_traces = 0;
  for (const auto& trace : drift_traces)
  {
    if (trace.getSize() == 0) { continue; }
    const auto apex = std::max_element(trace.begin(), trace.end(), [](const auto& a, const auto& b) { return a.getIntensity() < b.getIntensity(); });
    if (std::abs(apex->getRT() - 440.696651312) < 1e-6 && std::abs(apex->getMZ() - 694.59539794921875) < 1e-6)
    {
      drift_trace = &trace;
      ++matching_traces;
    }
  }
  TEST_EQUAL(matching_traces, 1)

  if (drift_trace != nullptr)
  {
    TOLERANCE_RELATIVE(1.0)
    TOLERANCE_ABSOLUTE(1e-6)
    TEST_EQUAL(drift_trace->getSize(), 119)
    TEST_REAL_SIMILAR(drift_trace->getTraceLength(), 211.062323823)
    TEST_REAL_SIMILAR(drift_trace->getCentroidMZ(), 694.8168204114992)
    const auto bounds
      = std::minmax_element(drift_trace->begin(), drift_trace->end(), [](const auto& a, const auto& b) { return a.getMZ() < b.getMZ(); });
    TEST_REAL_SIMILAR(bounds.first->getMZ(), 694.24554443359375)
    TEST_REAL_SIMILAR(bounds.second->getMZ(), 695.4852294921875)
    TEST_REAL_SIMILAR(bounds.second->getMZ() - bounds.first->getMZ(), 1.23968505859375)

    // Check the endpoints, seed, intensity maximum, an intermediate excursion,
    // and both mass extrema.
    struct ExpectedPeak
    {
      Size index;
      double rt;
      double mz;
      double intensity;
    };
    const ExpectedPeak expected_peaks[] = {
      {0, 388.858178160, 694.7713623046875, 5607.41455078125},     {3, 393.578101856, 694.822998046875, 7240.923828125},
      {18, 440.696651312, 694.59539794921875, 12166.330078125},    {65, 516.666194319, 695.2978515625, 3250.697265625},
      {104, 577.907320015, 695.4852294921875, 1828.6632080078125}, {117, 598.670037503, 694.24554443359375, 2687.224609375},
      {118, 599.920501983, 694.89813232421875, 3135.842529296875},
    };
    for (const auto& expected : expected_peaks)
    {
      // Keep a changed trace length from causing an out-of-bounds test access.
      if (expected.index < drift_trace->getSize())
      {
        const auto& peak = (*drift_trace)[expected.index];
        TOLERANCE_ABSOLUTE(1e-6)
        TEST_REAL_SIMILAR(peak.getRT(), expected.rt)
        TEST_REAL_SIMILAR(peak.getMZ(), expected.mz)
        TOLERANCE_ABSOLUTE(1e-3)
        TEST_REAL_SIMILAR(peak.getIntensity(), expected.intensity)
      }
    }
    TOLERANCE_RELATIVE(1.00001)
    TOLERANCE_ABSOLUTE(1e-5)
  }
}
END_SECTION

std::vector<MassTrace> filt;

// START_SECTION((void filterByPeakWidth(std::vector< MassTrace > &, std::vector< MassTrace > &)))
//{
//     test_mtd.filterByPeakWidth(output_mt, filt);

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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
