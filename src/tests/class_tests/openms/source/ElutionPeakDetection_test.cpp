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
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <OpenMS/PROCESSING/SMOOTHING/LowessSmoothing.h>

///////////////////////////
#include <OpenMS/FEATUREFINDER/ElutionPeakDetection.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ElutionPeakDetection, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ElutionPeakDetection* ptr = nullptr;
ElutionPeakDetection* null_ptr = nullptr;
START_SECTION(ElutionPeakDetection())
{
  ptr = new ElutionPeakDetection();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~ElutionPeakDetection())
{
  delete ptr;
}
END_SECTION

PeakMap input;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("ElutionPeakDetection_input1.mzML"), input);

std::vector<MassTrace> output_mt, splitted_mt, filtered_mt;

MassTraceDetection test_mtd;
Param mtd_def = MassTraceDetection().getDefaults();
test_mtd.setParameters(mtd_def);

test_mtd.run(input, output_mt);

ElutionPeakDetection test_epd;
Param epd_def = ElutionPeakDetection().getDefaults();
epd_def.setValue("width_filtering", "off");
epd_def.setValue("masstrace_snr_filtering", "false");
test_epd.setParameters(epd_def);

//String mt_labels[] = {"T6", "T7", "T9", "T3", "T4", "T8", "T5", "T2", "T1"};
//String split_labels[] = {"T6.1", "T6.2", "T6.3", "T6.4", "T7.1", "T7.2", "T7.3", "T7.4", "T7.5", "T9.1", "T9.2", "T9.3", "T3.1", "T3.2", "T3.3", "T3.4", "T3.5", "T3.6", "T3.7", "T3.8", "T3.9", "T3.10", "T3.11", "T3.12", "T3.13", "T4.1", "T4.2", "T4.3", "T8.1", "T8.2", "T5", "T2", "T1.1", "T1.2", "T1.3"};
//String filt_labels[] = {"T6.1", "T6.2", "T7.4", "T3.2", "T3.3", "T3.4", "T3.6", "T3.8", "T3.9", "T3.10", "T4.2", "T1.1", "T6.4", "T3.1", "T3.5", "T3.12", "T4.1", "T4.3", "T8.2", "T6.3", "T7.2", "T7.3", "T9.1", "T3.11", "T8.1", "T7.1", "T5", "T9.2", "T2", "T1.2", "T1.3", "T3.7"};


/* NOTE: The lowess smoothing was changed from using the GSL to a direct regression.
 * The smoothing test work fine for the modification. The ElutionPeakPicker shows
 * large differences and I (and Stefan) have no explanation for the reason.
 *
 * The test was adapted to use a SavitzkyGolay of polynomial-order 2 instead of a lowess smoothing
 * using the GSL because the performance seems to be much more robust (at least when going with
 * our unit tests). In addition, I tested how a lowess smoothing with regression performs. The
 * differences of the test are given as comments. The maintainer must decide on how
 * to handle the situation and which smoothing to keep. I had problems getting all unit
 * test running when using lowess smoothing with regression.
 *
 */
TOLERANCE_RELATIVE(1.01)
START_SECTION((void detectPeaks(std::vector< MassTrace > &, std::vector< MassTrace > &)))
{
    TEST_EQUAL(output_mt.size(), 1);

    if (!output_mt.empty())
    {
        TEST_EQUAL(output_mt[0].getLabel(), "T1");

        test_epd.detectPeaks(output_mt, splitted_mt);

        // mass traces split to local peaks
        //TEST_EQUAL(splitted_mt.size(), 2); // lowess and GSL
        TEST_EQUAL(splitted_mt.size(), 3); // SavitzkyGolay
        //TEST_EQUAL(splitted_mt.size(), 6); // lowess with regression

        // correct labeling if subtraces?
        TEST_EQUAL(splitted_mt[0].getLabel(), "T1.1");//lowess and GSL / SavitzkyGolay / lowess with regression
        TEST_EQUAL(splitted_mt[1].getLabel(), "T1.2");//lowess and GSL / SavitzkyGolay / lowess with regression
        //        TEST_EQUAL(splitted_mt[2].getLabel(), "T1.3");//lowess with regression
        //        TEST_EQUAL(splitted_mt[3].getLabel(), "T1.4");//lowess with regression
        //        TEST_EQUAL(splitted_mt[4].getLabel(), "T1.5");//lowess with regression
        //        TEST_EQUAL(splitted_mt[5].getLabel(), "T1.6");//lowess with regression
    }
}
END_SECTION

START_SECTION((void detectPeaks(MassTrace &, std::vector< MassTrace > &)))
{
  NOT_TESTABLE; // see above
}
END_SECTION

START_SECTION((void filterByPeakWidth(std::vector< MassTrace > &, std::vector< MassTrace > &)))
{
  NOT_TESTABLE;
  //    test_epd.filterByPeakWidth(splitted_mt, filtered_mt);

  //    TEST_EQUAL(filtered_mt.size(), 32);

  //    for (Size i = 0; i < filtered_mt.size(); ++i)
  //    {
  //        TEST_EQUAL(filtered_mt[i].getLabel(), filt_labels[i]);
  //    }
}
END_SECTION

START_SECTION((void findLocalExtrema(const MassTrace &, const Size &, std::vector< Size > &, std::vector< Size > &)))
{
  std::vector<Size> maxes, mins;

  if (!output_mt.empty())
  {
    MassTrace mt(output_mt[0]);

    std::vector<double> rts, ints;

    for (MassTrace::const_iterator c_it = mt.begin(); c_it != mt.end(); ++c_it)
    {
        rts.push_back(c_it->getRT());
        ints.push_back(c_it->getIntensity());
    }

    std::vector<double> smoothed_data;

    double win_size = 20;
    test_epd.smoothData(mt, win_size);
    test_epd.findLocalExtrema(mt, win_size/2, maxes, mins);

    // lowess and GSL
    //TEST_EQUAL(maxes.size(), 5);
    //TEST_EQUAL(mins.size(), 1);
    // SavitzkyGolay
    TEST_EQUAL(maxes.size(), 4);
    TEST_EQUAL(mins.size(), 2);

    // test window overlap
    mt = output_mt[0];
    test_epd.smoothData(mt, win_size);

    // The two largest peaks in the elution profile are about 90 spectra appart
    double distance_between_peaks = 90 - 20; // don't include other maximum but induce overlap 
    test_epd.findLocalExtrema(mt, distance_between_peaks, maxes, mins);
    TEST_EQUAL(maxes.size(), 2);
    TEST_EQUAL(mins.size(), 1);

    // lowess with regression
    //TEST_EQUAL(maxes.size(), 10);
    //TEST_EQUAL(mins.size(), 6);
  }
  // std::cerr << maxes.size() << " " << mins.size() << std::endl;
}
END_SECTION

splitted_mt.clear();
test_epd.detectPeaks(output_mt, splitted_mt);

START_SECTION((double computeMassTraceNoise(const MassTrace &)))
{
    TEST_EQUAL(output_mt.size(), 1);
    
    ABORT_IF(output_mt.empty())
    double est_noise(test_epd.computeMassTraceNoise(output_mt[0]));

    //TEST_REAL_SIMILAR(est_noise, 515.297);//using lowess and GSL
    TEST_REAL_SIMILAR(est_noise, 573.8585);//using SavitzkyGolay
    //TEST_REAL_SIMILAR(est_noise, 49027.69);//using lowess with regression
}
END_SECTION

START_SECTION((double computeMassTraceSNR(const MassTrace &)))
{
    ABORT_IF(splitted_mt.size() != 3);

    double snr1(test_epd.computeMassTraceSNR(splitted_mt[0]));
    double snr2(test_epd.computeMassTraceSNR(splitted_mt[1]));
    double snr3(test_epd.computeMassTraceSNR(splitted_mt[2]));

    // using lowess and GSL
    //TEST_REAL_SIMILAR(snr1, 8.6058);
    //TEST_REAL_SIMILAR(snr2, 8.946);
    // using SavitzkyGolay
    TEST_REAL_SIMILAR(snr1, 0.1907);
    TEST_REAL_SIMILAR(snr2, 9.8855);
    TEST_REAL_SIMILAR(snr3, 7.6432);
    // using lowess with regression
    //TEST_REAL_SIMILAR(snr1, 0.0497);
    //TEST_REAL_SIMILAR(snr2, 0.1450);
}
END_SECTION

START_SECTION((double computeApexSNR(const MassTrace &)))
{
    ABORT_IF(splitted_mt.size() != 3);

    double snr1(test_epd.computeApexSNR(splitted_mt[0]));
    double snr2(test_epd.computeApexSNR(splitted_mt[1]));
    double snr3(test_epd.computeApexSNR(splitted_mt[2]));

    // using lowess and GSL
    //TEST_REAL_SIMILAR(snr1, 40.0159);
    //TEST_REAL_SIMILAR(snr2, 58.5950);
    // using SavitzkyGolay
    TEST_REAL_SIMILAR(snr1,  2.0427);
    TEST_REAL_SIMILAR(snr2, 37.7893);
    TEST_REAL_SIMILAR(snr3, 52.9933);
    // using lowess with regression
    //TEST_REAL_SIMILAR(snr1, 6.5177);
    //TEST_REAL_SIMILAR(snr2, 7.3813);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



