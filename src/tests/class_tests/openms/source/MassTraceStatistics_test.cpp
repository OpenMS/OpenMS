// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/KERNEL/MassTraceStatistics.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/MassTrace.h>

#include <vector>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MassTraceStatistics, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// a symmetric elution peak: RT 1..5, m/z 100.0..100.4, intensities 0, 10, 20, 10, 0
const double its[] = {0.0, 10.0, 20.0, 10.0, 0.0};
vector<Peak2D> peaks;
for (Size i = 0; i < 5; ++i)
{
  Peak2D p;
  p.setRT(1.0 + i);
  p.setMZ(100.0 + 0.1 * i);
  p.setIntensity(static_cast<Peak2D::IntensityType>(its[i]));
  peaks.push_back(p);
}
const vector<double> smoothed = {1.0, 8.0, 16.0, 12.0, 2.0};

const vector<Peak2D> no_peaks;
const vector<double> no_smoothed;
const vector<Peak2D> one_peak(1, peaks[2]);
const vector<double> one_smoothed(1, 16.0);

// the three central peaks (intensities 10, 20, 10) as a part of the trace
const MassTraceStatistics::PeakRange center = MassTraceStatistics::PeakRange(peaks).subspan(1, 3);
const MassTraceStatistics::IntensityRange center_smoothed = MassTraceStatistics::IntensityRange(smoothed).subspan(1, 3);

START_SECTION((double meanMZ(PeakRange peaks)))
{
  TEST_REAL_SIMILAR(MassTraceStatistics::meanMZ(peaks), 100.2)
  TEST_REAL_SIMILAR(MassTraceStatistics::meanMZ(MassTraceStatistics::PeakRange(peaks).first(2)), 100.05)
  TEST_REAL_SIMILAR(MassTraceStatistics::meanMZ(one_peak), 100.2)
  TEST_EXCEPTION(Exception::InvalidValue, MassTraceStatistics::meanMZ(no_peaks))
}
END_SECTION

START_SECTION((double medianMZ(PeakRange peaks)))
{
  TEST_REAL_SIMILAR(MassTraceStatistics::medianMZ(peaks), 100.2)
  TEST_REAL_SIMILAR(MassTraceStatistics::medianMZ(MassTraceStatistics::PeakRange(peaks).first(4)), 100.15)
  TEST_REAL_SIMILAR(MassTraceStatistics::medianMZ(one_peak), 100.2)
  TEST_EXCEPTION(Exception::InvalidValue, MassTraceStatistics::medianMZ(no_peaks))
}
END_SECTION

START_SECTION((double weightedMeanMZ(PeakRange peaks)))
{
  TEST_REAL_SIMILAR(MassTraceStatistics::weightedMeanMZ(peaks), 100.2)
  // (0 * 100.0 + 10 * 100.1 + 20 * 100.2) / 30
  TEST_REAL_SIMILAR(MassTraceStatistics::weightedMeanMZ(MassTraceStatistics::PeakRange(peaks).first(3)), 100.1666666667)
  TEST_REAL_SIMILAR(MassTraceStatistics::weightedMeanMZ(one_peak), 100.2)
  TEST_EXCEPTION(Exception::InvalidValue, MassTraceStatistics::weightedMeanMZ(no_peaks))

  // several peaks without any intensity
  vector<Peak2D> zero_int(peaks.begin(), peaks.begin() + 2);
  zero_int[1].setIntensity(0.0f);
  TEST_EXCEPTION(Exception::InvalidValue, MassTraceStatistics::weightedMeanMZ(zero_int))
}
END_SECTION

START_SECTION((double weightedMZsd(PeakRange peaks, double centroid_mz)))
{
  // sqrt((10 * 0.1^2 + 20 * 0 + 10 * 0.1^2) / 40)
  TEST_REAL_SIMILAR(MassTraceStatistics::weightedMZsd(peaks, 100.2), 0.0707106781)
  TEST_REAL_SIMILAR(MassTraceStatistics::weightedMZsd(center, 100.2), 0.0707106781)
  TEST_EXCEPTION(Exception::InvalidValue, MassTraceStatistics::weightedMZsd(no_peaks, 100.2))
  TEST_EXCEPTION(Exception::InvalidValue, MassTraceStatistics::weightedMZsd(MassTraceStatistics::PeakRange(peaks).first(1), 100.2))
}
END_SECTION

START_SECTION((double weightedMeanRT(PeakRange peaks)))
{
  TEST_REAL_SIMILAR(MassTraceStatistics::weightedMeanRT(peaks), 3.0)
  // weights (trapezoidal rule): 0 * 0.5, 10 * 1, 20 * 0.5 -> (10 * 2 + 10 * 3) / 20
  TEST_REAL_SIMILAR(MassTraceStatistics::weightedMeanRT(MassTraceStatistics::PeakRange(peaks).first(3)), 2.5)
  TEST_REAL_SIMILAR(MassTraceStatistics::weightedMeanRT(one_peak), 3.0)
  TEST_EXCEPTION(Exception::InvalidValue, MassTraceStatistics::weightedMeanRT(no_peaks))
}
END_SECTION

START_SECTION((double smoothedWeightedMeanRT(PeakRange peaks, IntensityRange smoothed)))
{
  // (1 * 1 + 8 * 2 + 16 * 3 + 12 * 4 + 2 * 5) / 39
  TEST_REAL_SIMILAR(MassTraceStatistics::smoothedWeightedMeanRT(peaks, smoothed), 3.1538461538)
  // (8 * 2 + 16 * 3 + 12 * 4) / 36
  TEST_REAL_SIMILAR(MassTraceStatistics::smoothedWeightedMeanRT(center, center_smoothed), 3.1111111111)
  TEST_REAL_SIMILAR(MassTraceStatistics::smoothedWeightedMeanRT(one_peak, one_smoothed), 3.0)
  TEST_EXCEPTION(Exception::InvalidValue, MassTraceStatistics::smoothedWeightedMeanRT(peaks, no_smoothed))
  TEST_EXCEPTION(Exception::InvalidValue, MassTraceStatistics::smoothedWeightedMeanRT(peaks, vector<double>(5, 0.0)))
}
END_SECTION

START_SECTION((double smoothedMaxRT(PeakRange peaks, IntensityRange smoothed)))
{
  TEST_REAL_SIMILAR(MassTraceStatistics::smoothedMaxRT(peaks, smoothed), 3.0)
  TEST_REAL_SIMILAR(MassTraceStatistics::smoothedMaxRT(MassTraceStatistics::PeakRange(peaks).last(2), MassTraceStatistics::IntensityRange(smoothed).last(2)), 4.0)
  TEST_REAL_SIMILAR(MassTraceStatistics::smoothedMaxRT(one_peak, one_smoothed), 3.0)
  TEST_EXCEPTION(Exception::InvalidValue, MassTraceStatistics::smoothedMaxRT(peaks, no_smoothed))
  TEST_EXCEPTION(Exception::InvalidValue, MassTraceStatistics::smoothedMaxRT(peaks, vector<double>(5, 0.0)))
}
END_SECTION

START_SECTION((double medianRT(PeakRange peaks)))
{
  TEST_REAL_SIMILAR(MassTraceStatistics::medianRT(peaks), 3.0)
  TEST_REAL_SIMILAR(MassTraceStatistics::medianRT(MassTraceStatistics::PeakRange(peaks).first(4)), 2.5)
  TEST_REAL_SIMILAR(MassTraceStatistics::medianRT(one_peak), 3.0)
  TEST_EXCEPTION(Exception::InvalidValue, MassTraceStatistics::medianRT(no_peaks))
}
END_SECTION

START_SECTION((Size apexIndex(PeakRange peaks)))
{
  TEST_EQUAL(MassTraceStatistics::apexIndex(peaks), 2)
  TEST_EQUAL(MassTraceStatistics::apexIndex(center), 1)
  TEST_EQUAL(MassTraceStatistics::apexIndex(MassTraceStatistics::PeakRange(peaks).subspan(3, 2)), 0)
  // the first of several equally intense peaks
  TEST_EQUAL(MassTraceStatistics::apexIndex(vector<Peak2D>{peaks[1], peaks[3]}), 0)
  TEST_EXCEPTION(Exception::InvalidValue, MassTraceStatistics::apexIndex(no_peaks))
}
END_SECTION

START_SECTION((Size apexIndex(PeakRange peaks, IntensityRange smoothed)))
{
  TEST_EQUAL(MassTraceStatistics::apexIndex(peaks, smoothed), 2)
  TEST_EQUAL(MassTraceStatistics::apexIndex(MassTraceStatistics::PeakRange(peaks).last(2), MassTraceStatistics::IntensityRange(smoothed).last(2)), 0)
  TEST_EXCEPTION(Exception::InvalidValue, MassTraceStatistics::apexIndex(peaks, no_smoothed))
  TEST_EXCEPTION(Exception::InvalidValue, MassTraceStatistics::apexIndex(no_peaks, smoothed))
}
END_SECTION

START_SECTION((FWHM estimateFWHM(PeakRange peaks)))
{
  // half maximum (10) is reached exactly at RT 2 and RT 4
  MassTraceStatistics::FWHM fwhm = MassTraceStatistics::estimateFWHM(peaks);
  TEST_REAL_SIMILAR(fwhm.fwhm, 2.0)
  TEST_EQUAL(fwhm.start_idx, 0)
  TEST_EQUAL(fwhm.end_idx, 4)

  // a part of the trace: both flanks end at half maximum
  fwhm = MassTraceStatistics::estimateFWHM(center);
  TEST_REAL_SIMILAR(fwhm.fwhm, 2.0)
  TEST_EQUAL(fwhm.start_idx, 0)
  TEST_EQUAL(fwhm.end_idx, 2)

  // interpolation between the peaks enclosing the half maximum (7.5): RT 1.75 and RT 3.5
  vector<Peak2D> skewed(peaks.begin(), peaks.begin() + 4);
  skewed[2].setIntensity(15.0f);
  skewed[3].setIntensity(0.0f);
  fwhm = MassTraceStatistics::estimateFWHM(skewed);
  TEST_REAL_SIMILAR(fwhm.fwhm, 1.75)
  TEST_EQUAL(fwhm.start_idx, 0)
  TEST_EQUAL(fwhm.end_idx, 3)

  // no peaks left or right of the apex
  fwhm = MassTraceStatistics::estimateFWHM(MassTraceStatistics::PeakRange(peaks).first(3));
  TEST_REAL_SIMILAR(fwhm.fwhm, 0.0)
  TEST_EQUAL(fwhm.start_idx, 0)
  TEST_EQUAL(fwhm.end_idx, 0)
  fwhm = MassTraceStatistics::estimateFWHM(MassTraceStatistics::PeakRange(peaks).last(3));
  TEST_REAL_SIMILAR(fwhm.fwhm, 0.0)
  TEST_EQUAL(fwhm.start_idx, 0)
  TEST_EQUAL(fwhm.end_idx, 0)
  fwhm = MassTraceStatistics::estimateFWHM(one_peak);
  TEST_REAL_SIMILAR(fwhm.fwhm, 0.0)

  TEST_EXCEPTION(Exception::InvalidValue, MassTraceStatistics::estimateFWHM(no_peaks))

  // same result as the caching member of MassTrace
  MassTrace mt(peaks);
  TEST_REAL_SIMILAR(mt.estimateFWHM(false), 2.0)
  TEST_REAL_SIMILAR(mt.getFWHM(), 2.0)
  TEST_EQUAL(mt.getFWHMborders().first, 0)
  TEST_EQUAL(mt.getFWHMborders().second, 4)
}
END_SECTION

START_SECTION((FWHM estimateFWHM(PeakRange peaks, IntensityRange smoothed)))
{
  // half maximum (8) is reached at RT 2 and, between 12 (RT 4) and 2 (RT 5), at RT 4.4
  MassTraceStatistics::FWHM fwhm = MassTraceStatistics::estimateFWHM(peaks, smoothed);
  TEST_REAL_SIMILAR(fwhm.fwhm, 2.4)
  TEST_EQUAL(fwhm.start_idx, 0)
  TEST_EQUAL(fwhm.end_idx, 4)

  // a part of the trace: the right flank stays above the half maximum -> RT of its last peak
  fwhm = MassTraceStatistics::estimateFWHM(center, center_smoothed);
  TEST_REAL_SIMILAR(fwhm.fwhm, 2.0)
  TEST_EQUAL(fwhm.start_idx, 0)
  TEST_EQUAL(fwhm.end_idx, 2)

  TEST_EXCEPTION(Exception::InvalidValue, MassTraceStatistics::estimateFWHM(peaks, no_smoothed))
  TEST_EXCEPTION(Exception::InvalidValue, MassTraceStatistics::estimateFWHM(no_peaks, smoothed))

  // same result as the caching member of MassTrace
  MassTrace mt(peaks);
  mt.setSmoothedIntensities(smoothed);
  TEST_REAL_SIMILAR(mt.estimateFWHM(true), 2.4)
  TEST_EQUAL(mt.getFWHMborders().first, 0)
  TEST_EQUAL(mt.getFWHMborders().second, 4)
}
END_SECTION

START_SECTION((double intensitySum(PeakRange peaks)))
{
  TEST_REAL_SIMILAR(MassTraceStatistics::intensitySum(peaks), 40.0)
  TEST_REAL_SIMILAR(MassTraceStatistics::intensitySum(MassTraceStatistics::PeakRange(peaks).first(3)), 30.0)
  TEST_REAL_SIMILAR(MassTraceStatistics::intensitySum(no_peaks), 0.0)
}
END_SECTION

START_SECTION((double maxIntensity(PeakRange peaks)))
{
  TEST_REAL_SIMILAR(MassTraceStatistics::maxIntensity(peaks), 20.0)
  TEST_REAL_SIMILAR(MassTraceStatistics::maxIntensity(MassTraceStatistics::PeakRange(peaks).last(2)), 10.0)
  TEST_REAL_SIMILAR(MassTraceStatistics::maxIntensity(no_peaks), 0.0)
}
END_SECTION

START_SECTION((double maxIntensity(IntensityRange intensities)))
{
  TEST_REAL_SIMILAR(MassTraceStatistics::maxIntensity(smoothed), 16.0)
  TEST_REAL_SIMILAR(MassTraceStatistics::maxIntensity(MassTraceStatistics::IntensityRange(smoothed).last(2)), 12.0)
  TEST_REAL_SIMILAR(MassTraceStatistics::maxIntensity(no_smoothed), 0.0)
}
END_SECTION

START_SECTION((double medianIntensity(PeakRange peaks)))
{
  TEST_REAL_SIMILAR(MassTraceStatistics::medianIntensity(peaks), 10.0)
  TEST_REAL_SIMILAR(MassTraceStatistics::medianIntensity(MassTraceStatistics::PeakRange(peaks).first(4)), 10.0)
  TEST_REAL_SIMILAR(MassTraceStatistics::medianIntensity(MassTraceStatistics::PeakRange(peaks).first(2)), 5.0)
}
END_SECTION

START_SECTION((double peakArea(PeakRange peaks)))
{
  // 5 + 15 + 15 + 5
  TEST_REAL_SIMILAR(MassTraceStatistics::peakArea(peaks), 40.0)
  TEST_REAL_SIMILAR(MassTraceStatistics::peakArea(center), 30.0)
  TEST_REAL_SIMILAR(MassTraceStatistics::peakArea(one_peak), 0.0)
  TEST_REAL_SIMILAR(MassTraceStatistics::peakArea(no_peaks), 0.0)

  // the RT distance of the peaks is taken into account
  vector<Peak2D> sparse = {peaks[1], peaks[3]};
  TEST_REAL_SIMILAR(MassTraceStatistics::peakArea(sparse), 20.0)
}
END_SECTION

START_SECTION((double peakArea(PeakRange peaks, IntensityRange smoothed)))
{
  // 4.5 + 12 + 14 + 7
  TEST_REAL_SIMILAR(MassTraceStatistics::peakArea(peaks, smoothed), 37.5)
  TEST_REAL_SIMILAR(MassTraceStatistics::peakArea(center, center_smoothed), 26.0)
  TEST_REAL_SIMILAR(MassTraceStatistics::peakArea(one_peak, one_smoothed), 0.0)
  TEST_REAL_SIMILAR(MassTraceStatistics::peakArea(no_peaks, no_smoothed), 0.0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
