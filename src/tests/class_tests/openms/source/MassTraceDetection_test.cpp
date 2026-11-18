// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Erhan Kenar$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/MzMLFile.h>

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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
// Edit
START_SECTION((max_apex_deviation_ppm parameter))
{
  // Test that max_apex_deviation_ppm prevents runaway drift from the apex
  PeakMap exp;
  // Construct a synthetic run where peaks progressively drift away from apex
  for (Size s = 0; s < 10; ++s)
  {
    MSSpectrum spec;
    spec.setRT(s * 1.0);
    spec.setMSLevel(1);
    // Apex is at scan 0 with m/z 500.0; subsequent scans drift by 0.05 Th each
    Peak1D p;
    p.setMZ(500.0 + s * 0.015);
    p.setIntensity(10000.0 / (1.0 + s));
    spec.push_back(p);
    exp.addSpectrum(spec);
  }

  std::vector<MassTrace> traces;
  MassTraceDetection mtd;
  Param p = mtd.getDefaults();
  p.setValue("min_trace_length", 0.0);
  p.setValue("noise_threshold_int", 1.0);
  p.setValue("chrom_peak_snr", 1.0);

  // 1. With guard disabled (0.0): trace keeps collecting drifting peaks
  p.setValue("max_apex_deviation_ppm", 0.0);
  mtd.setParameters(p);
  mtd.run(exp, traces);
  TEST_EQUAL(traces.empty(), false);
  Size unconstrained_len = traces[0].getSize();

  // 2. With default guard (-1.0, automatically derived as 3 * mass_error_ppm = 60.0 ppm):
  // drifters beyond apex_mz + 0.030 Th are rejected!
  p.setValue("max_apex_deviation_ppm", -1.0);
  mtd.setParameters(p);
  traces.clear();
  mtd.run(exp, traces);
  TEST_EQUAL(traces.empty(), false);
  TEST_TRUE(traces[0].getSize() < unconstrained_len);

  // 3. Verify all collected peaks in the guarded trace are within 60 ppm of the apex
  for (const auto& pt : traces[0])
  {
    double ppm_dev = std::abs(pt.getMZ() - 500.0) / 500.0 * 1e6;
    TEST_TRUE(ppm_dev <= 60.0);
  }
}
    // Regression test for PR #10929:
  // Verify that findBestPeak_ does not pick an invalid peak outside the apex window
  // when an admissible alternative exists in the same scan.
  {
    PeakMap exp_two_peaks;

    // Scan 0: Apex at 500.000 (highest intensity)
    MSSpectrum s0;
    s0.setRT(0.0);
    s0.setMSLevel(1);
    s0.push_back(Peak1D(500.000, 20000.0));
    exp_two_peaks.addSpectrum(s0);

    // Scans 1-4: Subsequent scans at 500.024 shift the centroid upward (within 50 ppm)
    for (Size i = 1; i <= 4; ++i)
    {
      MSSpectrum s;
      s.setRT(i * 1.0);
      s.setMSLevel(1);
      s.push_back(Peak1D(500.024, 15000.0));
      exp_two_peaks.addSpectrum(s);
    }

    // Scan 5 contains two peaks:
    // Peak A: 500.000 (valid, 0 ppm from apex)
    // Peak B: 500.026 (outside 50 ppm apex envelope [limit is 500.025], but closer to centroid ~500.018)
    MSSpectrum s5;
    s5.setRT(5.0);
    s5.setMSLevel(1);
    s5.push_back(Peak1D(500.000, 15000.0));
    s5.push_back(Peak1D(500.026, 15000.0));
    exp_two_peaks.addSpectrum(s5);

    MassTraceDetection mtd_two_peaks;
    Param p_tp = mtd_two_peaks.getDefaults();
    p_tp.setValue("min_trace_length", 0.0);
    p_tp.setValue("noise_threshold_int", 1.0);
    p_tp.setValue("chrom_peak_snr", 1.0);
    p_tp.setValue("mass_error_ppm", 20.0);
    p_tp.setValue("max_apex_deviation_ppm", 50.0);
    mtd_two_peaks.setParameters(p_tp);

    std::vector<MassTrace> traces_tp;
    mtd_two_peaks.run(exp_two_peaks, traces_tp);

    TEST_EQUAL(traces_tp.empty(), false);
    TEST_EQUAL(traces_tp[0].getSize(), 6);
    TEST_REAL_SIMILAR(traces_tp[0].rbegin()->getMZ(), 500.000);
  }
END_SECTION

END_TEST
