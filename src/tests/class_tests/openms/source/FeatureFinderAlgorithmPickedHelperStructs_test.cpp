// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
///////////////////////////

#include <OpenMS/KERNEL/Peak1D.h>

using namespace OpenMS;
using namespace std;

START_TEST(FeatureFinderAlgorithmPickedHelperStructs, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(([FeatureFinderAlgorithmPickedHelperStructs::IsotopePattern] IsotopePattern(Size size)))
{
  const Size expected_size = 10;
  FeatureFinderAlgorithmPickedHelperStructs::IsotopePattern pattern(expected_size);

  TEST_EQUAL(pattern.intensity.size(), expected_size)
  TEST_EQUAL(pattern.mz_score.size(), expected_size)
  TEST_EQUAL(pattern.peak.size(), expected_size)
  TEST_EQUAL(pattern.spectrum.size(), expected_size)
  TEST_EQUAL(pattern.theoretical_mz.size(), expected_size)
}
END_SECTION

// MassTrace for testing
FeatureFinderAlgorithmPickedHelperStructs::MassTrace mt1;
mt1.theoretical_int = 0.8;

/////////////////////////////////////////////////////////////
Peak1D p1_1;
p1_1.setIntensity(1.08268226589f);
p1_1.setMZ(1000);
mt1.peaks.push_back(std::make_pair(677.1 , &p1_1));
Peak1D p1_2;
p1_2.setIntensity(1.58318959267f);
p1_2.setMZ(1000);
mt1.peaks.push_back(std::make_pair(677.4 , &p1_2));
Peak1D p1_3;
p1_3.setIntensity(2.22429840363f);
p1_3.setMZ(1000);
mt1.peaks.push_back(std::make_pair(677.7 , &p1_3));
Peak1D p1_4;
p1_4.setIntensity(3.00248879081f);
p1_4.setMZ(1000);
mt1.peaks.push_back(std::make_pair(678 , &p1_4));
Peak1D p1_5;
p1_5.setIntensity(3.89401804768f);
p1_5.setMZ(1000);
mt1.peaks.push_back(std::make_pair(678.3 , &p1_5));
Peak1D p1_6;
p1_6.setIntensity(4.8522452777f);
p1_6.setMZ(1000);
mt1.peaks.push_back(std::make_pair(678.6 , &p1_6));
Peak1D p1_7;
p1_7.setIntensity(5.80919229659f);
p1_7.setMZ(1000);
mt1.peaks.push_back(std::make_pair(678.9 , &p1_7));
Peak1D p1_8;
p1_8.setIntensity(6.68216169129f);
p1_8.setMZ(1000);
mt1.peaks.push_back(std::make_pair(679.2 , &p1_8));
Peak1D p1_9;
p1_9.setIntensity(7.38493077109f);
p1_9.setMZ(1000);
mt1.peaks.push_back(std::make_pair(679.5 , &p1_9));
Peak1D p1_10;
p1_10.setIntensity(7.84158938645f);
p1_10.setMZ(1000);
mt1.peaks.push_back(std::make_pair(679.8 , &p1_10));

START_SECTION(([FeatureFinderAlgorithmPickedHelperStructs::MassTrace] ConvexHull2D getConvexhull() const ))
{
  ConvexHull2D ch = mt1.getConvexhull();

  DPosition<2> point;
  point[0] = 679.8;
  point[1] = p1_10.getMZ();

  TEST_EQUAL(ch.encloses(point),true);

  point[1] = p1_10.getMZ() + 1.0;
  TEST_EQUAL(ch.encloses(point),false);

  point[1] = p1_10.getMZ();
  point[0] = 679.9;
  TEST_EQUAL(ch.encloses(point),false);
}
END_SECTION

START_SECTION(([FeatureFinderAlgorithmPickedHelperStructs::MassTrace] void updateMaximum()))
{
  mt1.updateMaximum();
  TEST_EQUAL(mt1.max_peak, &p1_10)
  TEST_EQUAL(mt1.max_rt, 679.8)
}
END_SECTION

START_SECTION(([FeatureFinderAlgorithmPickedHelperStructs::MassTrace] double getAvgMZ() const ))
{
  // getAvgMZ computes intensity weighted avg of the mass trace
  TEST_EQUAL(mt1.getAvgMZ(), 1000)

  FeatureFinderAlgorithmPickedHelperStructs::MassTrace mt_avg;

  Peak1D pAvg1;
  pAvg1.setMZ(10.5);
  pAvg1.setIntensity(1000);
  mt_avg.peaks.push_back(std::make_pair(100.0, &pAvg1));

  Peak1D pAvg2;
  pAvg2.setMZ(10.0);
  pAvg2.setIntensity(100);
  mt_avg.peaks.push_back(std::make_pair(100.0, &pAvg2));

  Peak1D pAvg3;
  pAvg3.setMZ(9.5);
  pAvg3.setIntensity(10);
  mt_avg.peaks.push_back(std::make_pair(100.0, &pAvg3));

  TEST_REAL_SIMILAR(mt_avg.getAvgMZ(), 10.4459)
}
END_SECTION

START_SECTION(([FeatureFinderAlgorithmPickedHelperStructs::MassTrace] bool isValid() const ))
{
  TEST_EQUAL(mt1.isValid(), true)
  FeatureFinderAlgorithmPickedHelperStructs::MassTrace mt_non_valid;

  mt_non_valid.peaks.push_back(std::make_pair(679.8 , &p1_10));
  TEST_EQUAL(mt_non_valid.isValid(), false)

  mt_non_valid.peaks.push_back(std::make_pair(679.5 , &p1_9));
  TEST_EQUAL(mt_non_valid.isValid(), false)

  mt_non_valid.peaks.push_back(std::make_pair(679.2 , &p1_8));
  TEST_EQUAL(mt_non_valid.isValid(), true)
}
END_SECTION

// testing mass trace
FeatureFinderAlgorithmPickedHelperStructs::MassTraces mt;
FeatureFinderAlgorithmPickedHelperStructs::MassTraces empty_traces;

// add a mass trace
mt.push_back(mt1);

START_SECTION(([FeatureFinderAlgorithmPickedHelperStructs::MassTraces] MassTraces()))
{
  TEST_EQUAL(mt.max_trace, 0)
}
END_SECTION

START_SECTION(([FeatureFinderAlgorithmPickedHelperStructs::MassTraces] Size getPeakCount() const ))
{
  TEST_EQUAL(mt.getPeakCount(), 10)
  TEST_EQUAL(empty_traces.getPeakCount(), 0)
}
END_SECTION

FeatureFinderAlgorithmPickedHelperStructs::MassTrace mt2;
mt2.theoretical_int = 0.2;

Peak1D p2_4;
p2_4.setIntensity(0.750622197703f);
p2_4.setMZ(1001);
mt2.peaks.push_back(std::make_pair(678, &p2_4));
Peak1D p2_5;
p2_5.setIntensity(0.97350451192f);
p2_5.setMZ(1001);
mt2.peaks.push_back(std::make_pair(678.3, &p2_5));
Peak1D p2_6;
p2_6.setIntensity(1.21306131943f);
p2_6.setMZ(1001);
mt2.peaks.push_back(std::make_pair(678.6, &p2_6));

mt.push_back(mt2);

START_SECTION(([FeatureFinderAlgorithmPickedHelperStructs::MassTraces] bool isValid(double seed_mz, double trace_tolerance)))
{
  // isValid checks if if we have enough traces
  FeatureFinderAlgorithmPickedHelperStructs::MassTraces invalid_traces;
  invalid_traces.push_back(mt1);

  TEST_EQUAL(invalid_traces.isValid(600.0, 0.03), false) // contains only one mass trace

  // and if the given seed is inside one of the mass traces
  TEST_EQUAL(mt.isValid(1000.0, 0.00), true)
  TEST_EQUAL(mt.isValid(1001.003, 0.03), true)
  TEST_EQUAL(mt.isValid(1002, 0.003), false)
}
END_SECTION

START_SECTION(([FeatureFinderAlgorithmPickedHelperStructs::MassTraces] Size getTheoreticalmaxPosition() const ))
{
  TEST_EXCEPTION(Exception::Precondition, empty_traces.getTheoreticalmaxPosition())

  TEST_EQUAL(mt.getTheoreticalmaxPosition(), 0)
}
END_SECTION

START_SECTION(([FeatureFinderAlgorithmPickedHelperStructs::MassTraces] void updateBaseline()))
{
  empty_traces.updateBaseline();
  TEST_EQUAL(empty_traces.baseline, 0.0)

  mt.updateBaseline();
  TEST_EQUAL(mt.baseline, p2_4.getIntensity())
}
END_SECTION

START_SECTION(([FeatureFinderAlgorithmPickedHelperStructs::MassTraces] std::pair<double,double> getRTBounds() const ))
{
  TEST_EXCEPTION(Exception::Precondition, empty_traces.getRTBounds())

  std::pair<double, double> bounds = mt.getRTBounds();
  TEST_EQUAL(bounds.first, 677.1)
  TEST_EQUAL(bounds.second, 679.8)
}
END_SECTION

// add some border cases to the traces that should be checked in computeIntensityProfile()

// add a leading peak to the second trace
Peak1D p2_0;
p2_0.setIntensity(0.286529652f);
p2_0.setMZ(1001);
mt[1].peaks.insert(mt[1].peaks.begin(), std::make_pair(676.8, &p2_0));

// .. add a peak after a gap
Peak1D p2_7;
p2_7.setIntensity(0.72952935f);
p2_7.setMZ(1001);
mt[1].peaks.push_back(std::make_pair(679.2, &p2_7));

// .. and a trailing peak
Peak1D p2_8;
p2_8.setIntensity(0.672624672f);
p2_8.setMZ(1001);
mt[1].peaks.push_back(std::make_pair(680.1, &p2_8));

START_SECTION(([FeatureFinderAlgorithmPickedHelperStructs::MassTraces] void computeIntensityProfile(std::list< std::pair<double, double> > intensity_profile) const))
{
  std::list< std::pair<double, double> > intensity_profile;
  mt.computeIntensityProfile(intensity_profile);

  TEST_EQUAL(intensity_profile.size(), 12)
  ABORT_IF(intensity_profile.size() != 12)

  std::list< std::pair<double, double> >::iterator profile = intensity_profile.begin();

  // the leading peak
  // 676.8 -> 0.286529652f
  TEST_REAL_SIMILAR(profile->first, 676.8)
  TEST_REAL_SIMILAR(profile->second, 0.286529652f)
  ++profile;

  // 677.1 -> 1.08268226589f
  TEST_REAL_SIMILAR(profile->first, 677.1)
  TEST_REAL_SIMILAR(profile->second, 1.08268226589f)
  ++profile;

  // 677.4 -> 1.58318959267f
  TEST_REAL_SIMILAR(profile->first, 677.4)
  TEST_REAL_SIMILAR(profile->second, 1.58318959267f)
  ++profile;

  // 677.7 -> 2.22429840363f
  TEST_REAL_SIMILAR(profile->first, 677.7)
  TEST_REAL_SIMILAR(profile->second, 2.22429840363f)
  ++profile;

  // 678.0 -> 3.00248879081f + 0.750622197703f
  TEST_REAL_SIMILAR(profile->first, 678.0)
  TEST_REAL_SIMILAR(profile->second, (3.00248879081f + 0.750622197703f))
  ++profile;

  // 678.3 -> 3.89401804768f + 0.97350451192f
  TEST_REAL_SIMILAR(profile->first, 678.3)
  TEST_REAL_SIMILAR(profile->second, (3.89401804768f + 0.97350451192f))
  ++profile;

  // 678.6 -> 4.8522452777f + 1.21306131943f
  TEST_REAL_SIMILAR(profile->first, 678.6)
  TEST_REAL_SIMILAR(profile->second, (4.8522452777f + 1.21306131943f))
  ++profile;

  // 678.9 -> 5.80919229659f
  TEST_REAL_SIMILAR(profile->first, 678.9)
  TEST_REAL_SIMILAR(profile->second, 5.80919229659f)
  ++profile;

  // 679.2 -> 6.68216169129f + 0.72952935f
  TEST_REAL_SIMILAR(profile->first, 679.2)
  TEST_REAL_SIMILAR(profile->second, (6.68216169129f + 0.72952935f))
  ++profile;

  // 679.5 -> 7.38493077109f
  TEST_REAL_SIMILAR(profile->first, 679.5)
  TEST_REAL_SIMILAR(profile->second, 7.38493077109f)
  ++profile;

  // 679.8 -> 7.84158938645f
  TEST_REAL_SIMILAR(profile->first, 679.8)
  TEST_REAL_SIMILAR(profile->second, 7.84158938645f)
  ++profile;

  // 680.1 -> 0.672624672f
  TEST_REAL_SIMILAR(profile->first, 680.1)
  TEST_REAL_SIMILAR(profile->second, 0.672624672f)
  ++profile;

  TEST_EQUAL(profile == intensity_profile.end(), true)

}
END_SECTION

START_SECTION(([FeatureFinderAlgorithmPickedHelperStructs::Seed] bool operator<(const Seed &rhs) const ))
{
  FeatureFinderAlgorithmPickedHelperStructs::Seed s1,s2,s3;
  s1.intensity = 100.0;
  s2.intensity = 200.0;
  s3.intensity = 300.0;

  TEST_EQUAL(s1.operator <(s2), true)
  TEST_EQUAL(s1.operator <(s3), true)
  TEST_EQUAL(s2.operator <(s3), true)

  TEST_EQUAL(s2.operator <(s1), false)
  TEST_EQUAL(s3.operator <(s1), false)
  TEST_EQUAL(s3.operator <(s2), false)
}
END_SECTION

START_SECTION(([FeatureFinderAlgorithmPickedHelperStructs::TheoreticalIsotopePattern] Size size() const ))
{
  FeatureFinderAlgorithmPickedHelperStructs::TheoreticalIsotopePattern theo_pattern;
  TEST_EQUAL(theo_pattern.size(), 0)

  theo_pattern.intensity.push_back(0.7);
  theo_pattern.intensity.push_back(0.2);
  theo_pattern.intensity.push_back(0.1);

  TEST_EQUAL(theo_pattern.size(), 3)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



