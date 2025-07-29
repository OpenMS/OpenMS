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

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/weighted_mean.hpp>
#include <boost/accumulators/statistics/stats.hpp>

///////////////////////////
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

// Test access class to access protected static methods
class MassTraceDetectionAccess : public MassTraceDetection
{
public:
    using MassTraceDetection::updateIterativeWeightedMean_;
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



// load a mzML file for testing the algorithm
PeakMap input;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MassTraceDetection_input1.mzML"),input);

Size exp_mt_lengths[3] = {85, 31, 16};
double exp_mt_rts[3] = {348.673712773767, 347.107, 346.888}; // centroid RTs should be reasonably similar (isotopic traces)
double exp_mt_mzs[3] = {437.26675, 438.27241, 439.27594};
double exp_mt_ints[3] = {3382.20934460708, 664.763828332733, 109.490108620676};

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

// Test for updateIterativeWeightedMean_ function with boost accumulators
START_SECTION((template<typename Accumulator> void updateIterativeWeightedMean_(double&, Accumulator&, double, double)))
{
    // Define accumulator type for weighted mean
    namespace ba = boost::accumulators;
    typedef ba::accumulator_set<double, ba::stats<ba::tag::weighted_mean>, double> accumulator_t;

    double apex_mz(150.22), apex_int(25000000);
    double new_mz1(150.34), new_int1(23043030);
    double new_mz2(150.11), new_int2(1932392);

    // Test the iterative weighted mean calculation using boost accumulators
    double centroid_mz(0.0);
    accumulator_t mz_accumulator;

    // First call: add the apex peak
    MassTraceDetectionAccess::updateIterativeWeightedMean_(centroid_mz, mz_accumulator, apex_int, apex_mz);
    
    // After first call, centroid should be the apex m/z since we only have one peak
    TEST_REAL_SIMILAR(centroid_mz, apex_mz);

    // Second call: add the first new peak
    MassTraceDetectionAccess::updateIterativeWeightedMean_(centroid_mz, mz_accumulator, new_int1, new_mz1);
    
    // Expected weighted mean after two peaks
    double wmean1 = (apex_mz * apex_int + new_mz1 * new_int1) / (apex_int + new_int1);
    TEST_REAL_SIMILAR(centroid_mz, wmean1);

    // Third call: add the second new peak
    MassTraceDetectionAccess::updateIterativeWeightedMean_(centroid_mz, mz_accumulator, new_int2, new_mz2);
    
    // Expected weighted mean after three peaks
    double wmean2 = (apex_mz * apex_int + new_mz1 * new_int1 + new_mz2 * new_int2) / (apex_int + new_int1 + new_int2);
    TEST_REAL_SIMILAR(centroid_mz, wmean2);

}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
