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

#include <algorithm>
#include <vector>

using namespace OpenMS;
using namespace std;

// helper: build a mass trace (first RT 100 s) from an intensity profile with the given scan spacing
static MassTrace make_trace(const std::vector<double>& ints, double spacing)
{
  std::vector<Peak2D> pts;
  for (Size i = 0; i < ints.size(); ++i)
  {
    Peak2D p;
    p.setRT(100.0 + i * spacing);
    p.setMZ(230.1);
    p.setIntensity((Peak2D::IntensityType)ints[i]);
    pts.push_back(p);
  }
  MassTrace mt(pts);
  mt.setLabel("T1");
  return mt;
}

// Mass traces of an untargeted metabolomics run (Orbitrap, MS1 every 1.2 s), found by MassTraceDetection.
// Each one shows one way in which the smoothing used to misrepresent a mass trace (issue #2777).

// A peak with a long tail. Savitzky-Golay undershoots in the tail and was clamped to zero on the last but
// one point, which then counted as a valley and cut the last point off as a trace of its own.
static const std::vector<double> tail_ints = {3880189.0, 35452464.0, 278961792.0, 299236192.0, 98084320.0, 61078160.0,
                                              15377399.0, 5181145.5, 1230206.0, 352417.9, 170774.1, 82065.6, 97492.0};

// A peak that rises within one scan: the apex is the second point. The fit of the first frame, extrapolated
// to the first point, exceeded every other smoothed intensity.
static const std::vector<double> steep_onset_ints = {558573.1, 2166397.5, 1540729.8, 947377.1, 638997.1, 505671.9,
                                                     276991.7, 112347.6, 231729.8, 58555.4, 35435.2};

// A low plateau followed by a peak two orders of magnitude higher. The valley between them is 0.6 of the
// plateau, i.e. not a split point, but the fit undershot to zero at the foot of the peak.
static const std::vector<double> plateau_ints = {13315.2, 16256.3, 15380.5, 18509.5, 15250.1, 16484.5, 13546.0, 11016.4, 11527.2,
                                                 40689.8, 714812.6, 1197706.9, 1474887.0, 1806559.9, 1810322.6, 1240083.0,
                                                 1326237.9, 1018520.6, 1565404.0, 2018843.5, 1614379.2, 609284.1, 1460009.0,
                                                 633715.2, 680450.1, 1003226.9, 861034.8, 435865.6, 785215.8, 468304.4, 470155.0,
                                                 464454.9, 300043.4, 266869.3, 153281.4, 151645.1};

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

//std::string mt_labels[] = {"T6", "T7", "T9", "T3", "T4", "T8", "T5", "T2", "T1"};
//std::string split_labels[] = {"T6.1", "T6.2", "T6.3", "T6.4", "T7.1", "T7.2", "T7.3", "T7.4", "T7.5", "T9.1", "T9.2", "T9.3", "T3.1", "T3.2", "T3.3", "T3.4", "T3.5", "T3.6", "T3.7", "T3.8", "T3.9", "T3.10", "T3.11", "T3.12", "T3.13", "T4.1", "T4.2", "T4.3", "T8.1", "T8.2", "T5", "T2", "T1.1", "T1.2", "T1.3"};
//std::string filt_labels[] = {"T6.1", "T6.2", "T7.4", "T3.2", "T3.3", "T3.4", "T3.6", "T3.8", "T3.9", "T3.10", "T4.2", "T1.1", "T6.4", "T3.1", "T3.5", "T3.12", "T4.1", "T4.3", "T8.2", "T6.3", "T7.2", "T7.3", "T9.1", "T3.11", "T8.1", "T7.1", "T5", "T9.2", "T2", "T1.2", "T1.3", "T3.7"};


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
        // SavitzkyGolay: the two elution peaks (apex 8.6e3 at 963.9 s, apex 3.6e4 at 986.3 s). This used to be 3:
        // the first four points of the trace, baseline of intensity 40-60, were cut off as a sub-trace of their
        // own, because the Savitzky-Golay fit extrapolated to the first point gave 469 and the five points after
        // it were negative and clamped to zero -- a 'maximum' and a 'valley' that are not in the data (#2777).
        TEST_EQUAL(splitted_mt.size(), 2);
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
  // These pin the behaviour issue #2777 is about: a peak next to a split boundary, or one
  // sampled too sparsely to put its maximum in the interior, must come out with its apex as
  // the centroid RT rather than being dropped or positioned at the edge of its RT range.
  // Rejecting fragments by where their apex sits was tried and does not work -- the split loop
  // keeps the minimum as the last point of the left fragment, so the right fragment starts one
  // sample past it and its first sample is normally a genuine peak.

  // Two equally tall peaks separated by a minimum. The second fragment begins one sample past
  // the minimum, so its apex is its own first point -- yet it is higher than both neighbours in
  // the parent trace. Both peaks must survive, each centred on its apex.
  MassTrace two_peak_mt(make_trace({100, 500, 1000, 500, 100, 1000, 500, 100}, 2.0));
  std::vector<MassTrace> two_peak_out;
  test_epd.detectPeaks(two_peak_mt, two_peak_out); // test_epd has width_filtering "off"
  TEST_EQUAL(two_peak_out.size(), 2);
  if (two_peak_out.size() == 2)
  {
    TEST_REAL_SIMILAR(two_peak_out[0].getCentroidRT(), 104.0);
    TEST_REAL_SIMILAR(two_peak_out[1].getCentroidRT(), 110.0);
    TEST_REAL_SIMILAR(two_peak_out[0].getMaxIntensity(false), 1000.0);
    TEST_REAL_SIMILAR(two_peak_out[1].getMaxIntensity(false), 1000.0);
  }

  // Same, with a flat top starting on the first sample after the cut, so the apex equals its
  // neighbour instead of exceeding it.
  MassTrace flat_mt(make_trace({100, 500, 1000, 500, 100, 1000, 1000, 500, 100}, 2.0));
  std::vector<MassTrace> flat_out;
  test_epd.detectPeaks(flat_mt, flat_out);
  TEST_EQUAL(flat_out.size(), 2);
  if (flat_out.size() == 2)
  {
    TEST_REAL_SIMILAR(flat_out[1].getCentroidRT(), 110.0);
    TEST_REAL_SIMILAR(flat_out[1].getMaxIntensity(false), 1000.0);
  }

  // A peak sampled too sparsely to put its maximum in the interior is still a peak. Nothing is
  // split off it, and estimateFWHM() cannot bracket a half maximum for it, but it survives --
  // when the scan rate is close to the peak width most of a run's traces look like this.
  MassTrace sparse_mt(make_trace({10000, 5000, 2500, 1200, 600, 300, 150, 80, 40, 20, 10, 10}, 1.0));
  std::vector<MassTrace> sparse_out;
  test_epd.detectPeaks(sparse_mt, sparse_out);
  TEST_EQUAL(sparse_out.size(), 1);
  if (!sparse_out.empty())
  {
    TEST_EQUAL(sparse_out[0].findMaxByIntPeak(true), 0); // apex on the first point, and kept
  }

  // The remaining cases are real mass traces (see the top of this file) that the smoothing used to distort.

  // A tail does not end in a trace of its own. test_epd uses chrom_fwhm 5, i.e. a frame of 5 points here.
  MassTrace tail_mt(make_trace(tail_ints, 1.1));
  std::vector<MassTrace> tail_out;
  test_epd.detectPeaks(tail_mt, tail_out);
  TEST_EQUAL(tail_out.size(), 1); // was 2: T1.1 and the last point as T1.2
  if (!tail_out.empty())
  {
    TEST_EQUAL(tail_out[0].getSize(), tail_ints.size());
    TEST_REAL_SIMILAR(tail_out[0].getCentroidRT(), 100.0 + 3 * 1.1);
  }

  // A low plateau in front of a peak is not cut off at the foot of the peak.
  MassTrace plateau_mt(make_trace(plateau_ints, 1.1));
  std::vector<MassTrace> plateau_out;
  test_epd.detectPeaks(plateau_mt, plateau_out);
  TEST_EQUAL(plateau_out.size(), 1); // was 2: the plateau as T1.1, the peak as T1.2
  if (!plateau_out.empty())
  {
    TEST_EQUAL(plateau_out[0].getSize(), plateau_ints.size());
    TEST_REAL_SIMILAR(plateau_out[0].getCentroidRT(), 100.0 + 19 * 1.1);
  }

  // A peak that rises within one scan is not positioned at the first point of the trace (the symptom in the
  // title of the issue) ...
  ElutionPeakDetection wide_epd;
  Param wide_param = ElutionPeakDetection().getDefaults();
  wide_param.setValue("chrom_fwhm", 10.0); // a frame of 11 points, as long as the trace
  wide_param.setValue("width_filtering", "off");
  wide_epd.setParameters(wide_param);

  MassTrace steep_mt(make_trace(steep_onset_ints, 1.0));
  std::vector<MassTrace> steep_out;
  wide_epd.detectPeaks(steep_mt, steep_out);
  TEST_EQUAL(steep_out.size(), 1);
  if (!steep_out.empty())
  {
    // was 100.0, the first point. The raw apex is at 101.0; smoothing a peak this narrow with a frame as long
    // as the whole trace moves the smoothed apex one scan towards its tail.
    TEST_REAL_SIMILAR(steep_out[0].getCentroidRT(), 102.0);
    TEST_EQUAL(steep_out[0].findMaxByIntPeak(true), 2);
    TEST_TRUE(steep_out[0].getFWHM() > 0.0); // was 0, i.e. no area
  }

  // ... and, having a width now, is no longer removed by the default width filter.
  wide_param.setValue("width_filtering", "fixed");
  wide_epd.setParameters(wide_param);
  steep_mt = make_trace(steep_onset_ints, 1.0);
  steep_out.clear();
  wide_epd.detectPeaks(steep_mt, steep_out);
  TEST_EQUAL(steep_out.size(), 1); // was 0
  if (!steep_out.empty())
  {
    TEST_REAL_SIMILAR(steep_out[0].getCentroidRT(), 102.0);
  }
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
    TEST_EQUAL(mins.size(), 1); // was 2: the baseline in front of the first peak was split off (see detectPeaks)

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

START_SECTION((void smoothData(MassTrace& mt, int win_size) const))
{
  // A trace that is shorter than the frame is not smoothed.
  for (Size length = 0; length < 9; ++length)
  {
    std::vector<double> ints;
    for (Size i = 0; i < length; ++i) { ints.push_back(100.0 * (i + 1) * (i % 2 + 1)); }
    MassTrace short_mt(make_trace(ints, 1.0));
    test_epd.smoothData(short_mt, 9);
    TEST_EQUAL(short_mt.getSmoothedIntensities().size(), length);
    for (Size i = 0; i < short_mt.getSmoothedIntensities().size(); ++i) { TEST_REAL_SIMILAR(short_mt.getSmoothedIntensities()[i], ints[i]); }
  }

  // Where a frame fits, this is Savitzky-Golay smoothing with a polynomial of order 2. A frame of five points
  // has the coefficients (-3, 12, 17, 12, -3) / 35.
  std::vector<double> ints = {10.0, 20.0, 40.0, 90.0, 70.0, 30.0, 20.0, 10.0};
  MassTrace mt(make_trace(ints, 1.0));
  test_epd.smoothData(mt, 5);
  std::vector<double> smoothed = mt.getSmoothedIntensities();
  TEST_EQUAL(smoothed.size(), ints.size());
  ABORT_IF(smoothed.size() != ints.size());
  for (Size i = 2; i + 2 < ints.size(); ++i)
  {
    TEST_REAL_SIMILAR(smoothed[i], (-3 * ints[i - 2] + 12 * ints[i - 1] + 17 * ints[i] + 12 * ints[i + 1] - 3 * ints[i + 2]) / 35.0);
  }
  // No frame of five points fits around the two outermost points at each end: they are not extrapolated from
  // the fit of the neighbouring frame. The second and the last but one point are the mean of three points, the
  // first and the last point keep their intensity.
  TEST_REAL_SIMILAR(smoothed[0], ints[0]);
  TEST_REAL_SIMILAR(smoothed[1], (10.0 + 20.0 + 40.0) / 3.0);
  TEST_REAL_SIMILAR(smoothed[6], (30.0 + 20.0 + 10.0) / 3.0);
  TEST_REAL_SIMILAR(smoothed[7], ints[7]);

  // an even frame length is the next odd one, and three points (or less) do not smooth
  test_epd.smoothData(mt, 4);
  for (Size i = 0; i < ints.size(); ++i) { TEST_REAL_SIMILAR(mt.getSmoothedIntensities()[i], smoothed[i]); }
  for (int win_size : {-1, 0, 1, 2, 3})
  {
    test_epd.smoothData(mt, win_size);
    for (Size i = 0; i < ints.size(); ++i) { TEST_REAL_SIMILAR(mt.getSmoothedIntensities()[i], ints[i]); }
  }

  // Towards the ends the frame narrows to what fits around the point: with a frame of seven points,
  // (-2, 3, 6, 7, 6, 3, -2) / 21, the third point from either end is smoothed with five.
  test_epd.smoothData(mt, 7);
  smoothed = mt.getSmoothedIntensities();
  TEST_REAL_SIMILAR(smoothed[3], (-2 * 10.0 + 3 * 20.0 + 6 * 40.0 + 7 * 90.0 + 6 * 70.0 + 3 * 30.0 - 2 * 20.0) / 21.0);
  TEST_REAL_SIMILAR(smoothed[4], (-2 * 20.0 + 3 * 40.0 + 6 * 90.0 + 7 * 70.0 + 6 * 30.0 + 3 * 20.0 - 2 * 10.0) / 21.0);
  TEST_REAL_SIMILAR(smoothed[2], (-3 * 10.0 + 12 * 20.0 + 17 * 40.0 + 12 * 90.0 - 3 * 70.0) / 35.0);
  TEST_REAL_SIMILAR(smoothed[5], (-3 * 90.0 + 12 * 70.0 + 17 * 30.0 + 12 * 20.0 - 3 * 10.0) / 35.0);
  TEST_REAL_SIMILAR(smoothed[1], (10.0 + 20.0 + 40.0) / 3.0);
  TEST_REAL_SIMILAR(smoothed[0], ints[0]);

  // A noise spike on the second point does not keep its height (left unsmoothed, it would compete with the
  // smoothed, i.e. lowered, apex of the peak).
  std::vector<double> spike = {100.0, 1000.0, 300.0, 500.0, 800.0, 900.0, 800.0, 500.0, 300.0};
  MassTrace spike_mt(make_trace(spike, 1.0));
  test_epd.smoothData(spike_mt, 5);
  TEST_REAL_SIMILAR(spike_mt.getSmoothedIntensities()[1], (100.0 + 1000.0 + 300.0) / 3.0);
  TEST_EQUAL(spike_mt.findMaxByIntPeak(true), 5);

  // As before, smoothed intensities are never negative (here: a baseline-subtracted trace).
  std::vector<double> negative = {-50.0, -40.0, -30.0, 200.0, 600.0, 200.0, -30.0, -40.0, -50.0};
  MassTrace negative_mt(make_trace(negative, 1.0));
  test_epd.smoothData(negative_mt, 5);
  for (Size i = 1; i + 1 < negative.size(); ++i) { TEST_TRUE(negative_mt.getSmoothedIntensities()[i] >= 0.0); }

  // The smoothed intensity stays within the raw intensities of its frame. At the foot of a steep flank the
  // fit is negative (here (-3 * 100 + 12 * 100 + 17 * 100 + 12 * 110 - 3 * 100000) / 35 < 0) and used to be
  // clamped to zero, which made the foot a valley to split the trace at.
  std::vector<double> flank = {100.0, 100.0, 100.0, 100.0, 100.0, 110.0, 100000.0, 250000.0, 180000.0, 60000.0, 20000.0};
  MassTrace flank_mt(make_trace(flank, 1.0));
  test_epd.smoothData(flank_mt, 5);
  smoothed = flank_mt.getSmoothedIntensities();
  TEST_REAL_SIMILAR(smoothed[4], 100.0);
  for (Size i = 0; i < flank.size(); ++i)
  {
    const Size half_frame = std::min({Size(2), i, flank.size() - 1 - i});
    TEST_TRUE(smoothed[i] >= *std::min_element(flank.begin() + (i - half_frame), flank.begin() + (i + half_frame + 1)));
    TEST_TRUE(smoothed[i] <= *std::max_element(flank.begin() + (i - half_frame), flank.begin() + (i + half_frame + 1)));
  }

  // The smoothed apex of a peak that rises within one scan is not the first point of the trace (extrapolating
  // the fit of the first frame used to give 1430471 for it, more than anywhere else).
  MassTrace steep_mt(make_trace(steep_onset_ints, 1.0));
  test_epd.smoothData(steep_mt, 10);
  TEST_NOT_EQUAL(steep_mt.findMaxByIntPeak(true), 0);
  TEST_REAL_SIMILAR(steep_mt.getSmoothedIntensities()[0], steep_onset_ints[0]);
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
    ABORT_IF(splitted_mt.size() != 2);

    double snr1(test_epd.computeMassTraceSNR(splitted_mt[0]));
    double snr2(test_epd.computeMassTraceSNR(splitted_mt[1]));

    // using lowess and GSL
    //TEST_REAL_SIMILAR(snr1, 8.6058);
    //TEST_REAL_SIMILAR(snr2, 8.946);
    // using SavitzkyGolay
    // (the two elution peaks had 9.8855 and 7.6432 before the baseline in front of the first one stayed with it)
    TEST_REAL_SIMILAR(snr1, 9.8288);
    TEST_REAL_SIMILAR(snr2, 7.6175);
    // using lowess with regression
    //TEST_REAL_SIMILAR(snr1, 0.0497);
    //TEST_REAL_SIMILAR(snr2, 0.1450);
}
END_SECTION

START_SECTION((double computeApexSNR(const MassTrace &)))
{
    ABORT_IF(splitted_mt.size() != 2);

    double snr1(test_epd.computeApexSNR(splitted_mt[0]));
    double snr2(test_epd.computeApexSNR(splitted_mt[1]));

    // using lowess and GSL
    //TEST_REAL_SIMILAR(snr1, 40.0159);
    //TEST_REAL_SIMILAR(snr2, 58.5950);
    // using SavitzkyGolay
    // (37.7893 and 52.9933 before, see computeMassTraceSNR)
    TEST_REAL_SIMILAR(snr1, 38.5923);
    TEST_REAL_SIMILAR(snr2, 53.2332);
    // using lowess with regression
    //TEST_REAL_SIMILAR(snr1, 6.5177);
    //TEST_REAL_SIMILAR(snr2, 7.3813);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



