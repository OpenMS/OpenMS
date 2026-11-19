// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Erhan Kenar$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#include <algorithm>
#include <cmath>
#include <vector>

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

    // The spectra must be sorted by RT, MS2 spectra included (checked in debug builds). A file that lists all MS1 spectra
    // before its MS2 spectra (like the share/OpenMS/examples/FRACTIONS files ProteomicsLFQ is tested with) is not, even
    // though its MS1 spectra are in order; after MSExperiment::sortSpectra() the same traces are found.
    {
      PeakMap ms1_then_ms2 = input;
      MSSpectrum ms2;
      ms2.setMSLevel(2);
      ms2.setRT((input[0].getRT() + input[1].getRT()) / 2.0);
      ms2.push_back(Peak1D(500.0, 6000.0f));
      ms1_then_ms2.addSpectrum(ms2);
      TEST_EQUAL(ms1_then_ms2.isSorted(false), false)
      TEST_PRECONDITION_VIOLATED(test_mtd.run(ms1_then_ms2, output_mt))

      ms1_then_ms2.sortSpectra(false);
      output_mt.clear();
      test_mtd.run(ms1_then_ms2, output_mt);
      TEST_EQUAL(output_mt.size(), 3)
      for (Size i = 0; i < std::min(output_mt.size(), Size(3)); ++i) // the expected values are for 3 traces
      {
        TEST_EQUAL(output_mt[i].getSize(), exp_mt_lengths[i]);
        TEST_REAL_SIMILAR(output_mt[i].getCentroidRT(), exp_mt_rts[i]);
      }
    }
}
END_SECTION

START_SECTION([EXTRA] the re-estimated m/z standard deviation only narrows the m/z window (issue #10905))
{
  // A peak extends a trace if it lies within three standard deviations of the trace's centroid m/z.
  // The standard deviation starts at 'mass_error_ppm' (20 ppm, i.e. a window of +/- 60 ppm) and, with
  // 'reestimate_mt_sd', is re-estimated from the collected peaks. That estimate comes from the peaks the
  // window let in, so it may narrow the window, but must not widen it: every peak at the edge of the window
  // inflates the estimate, the larger window lets in more distant peaks, and the trace runs away (#10905).
  const double max_window_ppm = 3.0 * 20.0;

  // m/z range of a trace relative to its centroid m/z, in ppm
  auto span_ppm = [](const MassTrace& mt)
  {
    double mz_min = mt.begin()->getMZ();
    double mz_max = mz_min;
    for (const Peak2D& p : mt)
    {
      mz_min = std::min(mz_min, p.getMZ());
      mz_max = std::max(mz_max, p.getMZ());
    }
    return (mz_max - mz_min) / mt.getCentroidMZ() * 1e6;
  };

  // 1. A weak trace (8 scans at m/z 500), followed by peaks of the same intensity whose distance to it
  //    alternates in sign and grows by 17% per scan: +15, -17.6, ..., -45.0, +52.7, -61.6, +72.1, -84.4 ppm.
  //    Each of them is inside the window that its predecessors have inflated. The same holds for the peak
  //    search of ion mobility data, here with one ion mobility for all peaks.
  for (const bool with_ion_mobility : {false, true})
  {
    std::vector<double> offsets_ppm(8, 0.0);
    for (Size k = 0; k < 12; ++k)
    {
      offsets_ppm.push_back((k % 2 == 0 ? 1.0 : -1.0) * 15.0 * std::pow(1.17, static_cast<double>(k)));
    }
    PeakMap ladder;
    for (Size i = 0; i < offsets_ppm.size(); ++i)
    {
      MSSpectrum s;
      s.setRT(static_cast<double>(i));
      s.setMSLevel(1);
      Peak1D p;
      p.setMZ(500.0 * (1.0 + offsets_ppm[i] * 1e-6));
      p.setIntensity(i == 0 ? 1001.0 : 1000.0); // the trace is seeded in the first scan
      s.push_back(p);
      if (with_ion_mobility)
      {
        s.getFloatDataArrays().resize(1);
        s.getFloatDataArrays()[0].setName(Constants::UserParam::ION_MOBILITY);
        s.getFloatDataArrays()[0].push_back(1.0);
      }
      ladder.addSpectrum(s);
    }

    std::vector<MassTrace> traces;
    MassTraceDetection mtd;
    mtd.run(ladder, traces);
    TEST_EQUAL(mtd.hasCentroidIm(), with_ion_mobility)

    // the trace ends with the last peak within 60 ppm (+52.7 ppm); it used to take all 20 peaks (156 ppm)
    TEST_EQUAL(traces.size(), 1)
    ABORT_IF(traces.size() != 1)
    TEST_EQUAL(traces[0].getSize(), 17)
    for (const Peak2D& p : traces[0])
    {
      TEST_TRUE(std::fabs(p.getMZ() - traces[0].getCentroidMZ()) / traces[0].getCentroidMZ() * 1e6 <= max_window_ppm)
    }
  }

  // 2. A noise floor above 'noise_threshold_int', as in the data of #10905: peaks of similar intensity at
  //    random m/z (1.5 peaks per Th and scan). The mass traces through it must not become wider than the
  //    window to both sides of the centroid (the widest one has 97 ppm); they used to span up to 1.0 Th
  //    (1700 ppm) here.
  {
    UInt32 state = 2; // linear congruential generator, so that the data is the same on every platform
    auto next_random = [&state]()
    {
      state = state * 1664525u + 1013904223u;
      return state / 4294967296.0;
    };
    PeakMap noise;
    for (Size i = 0; i < 150; ++i)
    {
      MSSpectrum s;
      s.setRT(static_cast<double>(i));
      s.setMSLevel(1);
      for (Size j = 0; j < 30; ++j)
      {
        Peak1D p;
        p.setMZ(600.0 + 20.0 * next_random());
        p.setIntensity(1500.0 + 1500.0 * next_random());
        s.push_back(p);
      }
      s.sortByPosition();
      noise.addSpectrum(s);
    }

    std::vector<MassTrace> traces;
    MassTraceDetection().run(noise, traces);

    TEST_TRUE(traces.size() > 50)
    double max_span_ppm = 0.0;
    for (const MassTrace& mt : traces)
    {
      max_span_ppm = std::max(max_span_ppm, span_ppm(mt));
    }
    TEST_TRUE(max_span_ppm <= 2.0 * max_window_ppm)
  }

  // 3. The re-estimation still narrows the window. Ion A (m/z 500, a chromatographic peak with its apex in
  //    scan 12) is followed by the weaker ion B, 40 ppm above it. The 25 peaks of A narrow the window of its
  //    trace to less than 40 ppm, so that B is a trace of its own. Without re-estimation, the window stays at
  //    60 ppm and the trace of A takes the peaks of B as well.
  {
    PeakMap two_ions;
    for (Size i = 0; i < 40; ++i)
    {
      MSSpectrum s;
      s.setRT(static_cast<double>(i));
      s.setMSLevel(1);
      Peak1D p;
      if (i < 25)
      {
        p.setMZ(500.0);
        p.setIntensity(1e5 * std::exp(-0.5 * std::pow((static_cast<double>(i) - 12.0) / 4.0, 2)));
      }
      else
      {
        p.setMZ(500.0 * (1.0 + 40.0 * 1e-6));
        p.setIntensity(2000.0);
      }
      s.push_back(p);
      two_ions.addSpectrum(s);
    }

    std::vector<MassTrace> traces;
    MassTraceDetection mtd;
    mtd.run(two_ions, traces);
    TEST_EQUAL(traces.size(), 2)
    ABORT_IF(traces.size() != 2)
    TEST_EQUAL(traces[0].getSize(), 25)
    TEST_REAL_SIMILAR(traces[0].getCentroidMZ(), 500.0)
    TEST_EQUAL(traces[1].getSize(), 15)
    TEST_REAL_SIMILAR(traces[1].getCentroidMZ(), 500.02)

    Param p = mtd.getDefaults();
    p.setValue("reestimate_mt_sd", "false");
    mtd.setParameters(p);
    mtd.run(two_ions, traces);
    TEST_EQUAL(traces.size(), 1)
    ABORT_IF(traces.size() != 1)
    TEST_EQUAL(traces[0].getSize(), 40)
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
END_TEST
