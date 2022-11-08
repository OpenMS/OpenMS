// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
///////////////////////////

#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

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
double intensities[] = { 1.08268226589f, 0.270670566473f, 1.58318959267f, 0.395797398167f, 2.22429840363f, 0.556074600906f, 3.00248879081f, 0.750622197703f, 3.89401804768f, 0.97350451192f, 4.8522452777f, 1.21306131943f, 5.80919229659f, 1.45229807415f, 6.68216169129f, 1.67054042282f, 7.38493077109f, 1.84623269277f, 7.84158938645f, 1.96039734661f, 8.0f, 2.0f, 7.84158938645f, 1.96039734661f, 7.38493077109f, 1.84623269277f, 6.68216169129f, 1.67054042282f, 5.80919229659f, 1.45229807415f, 4.8522452777f, 1.21306131943f, 3.89401804768f, 0.97350451192f, 3.00248879081f, 0.750622197703f, 2.22429840363f, 0.556074600906f, 1.58318959267f, 0.395797398167f, 1.08268226589f, 0.270670566473f};
double rts[] = { 677.1, 677.1, 677.4, 677.4, 677.7, 677.7, 678, 678, 678.3, 678.3, 678.6, 678.6, 678.9, 678.9, 679.2, 679.2, 679.5, 679.5, 679.8, 679.8, 680.1, 680.1, 680.4, 680.4, 680.7, 680.7, 681, 681, 681.3, 681.3, 681.6, 681.6, 681.9, 681.9, 682.2, 682.2, 682.5, 682.5, 682.8, 682.8, 683.1, 683.1};

std::vector<Peak1D> all_peaks;
std::vector<MSSpectrum> all_spectra;
all_spectra.reserve(10);
all_peaks.reserve(10);

// Generate 10 peaks for mt1 (skip all mt2 peaks for now)
for (int k = 0; k < 20; k++)
{
  Peak1D p1; MSSpectrum s1;
  p1.setIntensity( intensities[k] );
  p1.setMZ(1000);
  s1.setRT( rts[k] );
  all_peaks.push_back(p1); all_spectra.push_back(s1);
  mt1.peaks.push_back(make_pair(&all_spectra.back(), &all_peaks.back()));
  std::cout << "1_" << k /2 +1 << " :: " << mt1.peaks.back().first->getRT() << " : " << mt1.peaks.back().second->getIntensity() << std::endl;

  k++;
}

START_SECTION(([FeatureFinderAlgorithmPickedHelperStructs::MassTrace] ConvexHull2D getConvexhull() const ))
{
  ConvexHull2D ch = mt1.getConvexhull();

  DPosition<2> point;
  point[0] = 679.8;
  point[1] = all_peaks[10-1].getMZ();

  TEST_EQUAL(ch.encloses(point),true);

  point[1] = all_peaks[10-1].getMZ() + 1.0;
  TEST_EQUAL(ch.encloses(point),false);

  point[1] = all_peaks[10-1].getMZ();
  point[0] = 679.9;
  TEST_EQUAL(ch.encloses(point),false);
}
END_SECTION

START_SECTION(([FeatureFinderAlgorithmPickedHelperStructs::MassTrace] void updateMaximum()))
{
  mt1.updateMaximum();
  TEST_EQUAL(mt1.max_peak, &all_peaks[9])
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
  MSSpectrum sAvg1;
  sAvg1.setRT(100.0);
  mt_avg.peaks.push_back(std::make_pair(&sAvg1, &pAvg1));

  Peak1D pAvg2;
  pAvg2.setMZ(10.0);
  pAvg2.setIntensity(100);
  MSSpectrum sAvg2;
  sAvg2.setRT(100.0);
  mt_avg.peaks.push_back(std::make_pair(&sAvg2, &pAvg2));

  Peak1D pAvg3;
  pAvg3.setMZ(9.5);
  pAvg3.setIntensity(10);
  MSSpectrum sAvg3;
  sAvg3.setRT(100.0);
  mt_avg.peaks.push_back(std::make_pair(&sAvg3, &pAvg3));

  TEST_REAL_SIMILAR(mt_avg.getAvgMZ(), 10.4459)
}
END_SECTION

START_SECTION(([FeatureFinderAlgorithmPickedHelperStructs::MassTrace] bool isValid() const ))
{
  TEST_EQUAL(mt1.isValid(), true)
  FeatureFinderAlgorithmPickedHelperStructs::MassTrace mt_non_valid;

  MSSpectrum s;
  s.setRT(679.8);
  mt_non_valid.peaks.push_back(std::make_pair(&s, &all_peaks[9]));
  TEST_EQUAL(mt_non_valid.isValid(), false)

  s.setRT(679.5);
  mt_non_valid.peaks.push_back(std::make_pair(&s, &all_peaks[8]));
  TEST_EQUAL(mt_non_valid.isValid(), false)

  s.setRT(679.2);
  mt_non_valid.peaks.push_back(std::make_pair(&s, &all_peaks[7]));
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
MSSpectrum s2_4;
s2_4.setRT(678);
mt2.peaks.push_back(std::make_pair(&s2_4, &p2_4));
Peak1D p2_5;
p2_5.setIntensity(0.97350451192f);
p2_5.setMZ(1001);
MSSpectrum s2_5;
s2_5.setRT(678.3);
mt2.peaks.push_back(std::make_pair(&s2_5, &p2_5));
Peak1D p2_6;
p2_6.setIntensity(1.21306131943f);
p2_6.setMZ(1001);
MSSpectrum s2_6;
s2_6.setRT(678.6);
mt2.peaks.push_back(std::make_pair(&s2_6, &p2_6));

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
MSSpectrum s2_0;
s2_0.setRT(676.8);
mt[1].peaks.insert(mt[1].peaks.begin(), std::make_pair(&s2_0, &p2_0));

// .. add a peak after a gap
Peak1D p2_7;
p2_7.setIntensity(0.72952935f);
p2_7.setMZ(1001);
MSSpectrum s2_7;
s2_7.setRT(679.2);
mt[1].peaks.push_back(std::make_pair(&s2_7, &p2_7));

// .. and a trailing peak
Peak1D p2_8;
p2_8.setIntensity(0.672624672f);
p2_8.setMZ(1001);
MSSpectrum s2_8;
s2_8.setRT(680.1);
mt[1].peaks.push_back(std::make_pair(&s2_8, &p2_8));

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



