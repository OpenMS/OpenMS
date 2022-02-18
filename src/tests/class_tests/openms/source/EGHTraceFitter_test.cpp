// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EGHTraceFitter.h>
///////////////////////////

#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <cmath>

using namespace OpenMS;
using namespace std;

#define PI 3.14159265358979323846

// TODO: include a more asymmetric trace in the test

START_TEST(EGHTraceFitter, "$Id$")

FeatureFinderAlgorithmPickedHelperStructs::MassTraces mts;

FeatureFinderAlgorithmPickedHelperStructs::MassTrace mt1;
mt1.theoretical_int = 0.8;

FeatureFinderAlgorithmPickedHelperStructs::MassTrace mt2;
mt2.theoretical_int = 0.2;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
// set up mass traces to fit


double intensities[] = { 1.08268226589f, 0.270670566473f, 1.58318959267f, 0.395797398167f, 2.22429840363f, 0.556074600906f, 3.00248879081f, 0.750622197703f, 3.89401804768f, 0.97350451192f, 4.8522452777f, 1.21306131943f, 5.80919229659f, 1.45229807415f, 6.68216169129f, 1.67054042282f, 7.38493077109f, 1.84623269277f, 7.84158938645f, 1.96039734661f, 8.0f, 2.0f, 7.84158938645f, 1.96039734661f, 7.38493077109f, 1.84623269277f, 6.68216169129f, 1.67054042282f, 5.80919229659f, 1.45229807415f, 4.8522452777f, 1.21306131943f, 3.89401804768f, 0.97350451192f, 3.00248879081f, 0.750622197703f, 2.22429840363f, 0.556074600906f, 1.58318959267f, 0.395797398167f, 1.08268226589f, 0.270670566473f};
double rts[] = { 677.1, 677.1, 677.4, 677.4, 677.7, 677.7, 678, 678, 678.3, 678.3, 678.6, 678.6, 678.9, 678.9, 679.2, 679.2, 679.5, 679.5, 679.8, 679.8, 680.1, 680.1, 680.4, 680.4, 680.7, 680.7, 681, 681, 681.3, 681.3, 681.6, 681.6, 681.9, 681.9, 682.2, 682.2, 682.5, 682.5, 682.8, 682.8, 683.1, 683.1};

std::vector<Peak1D> all_peaks;
std::vector<MSSpectrum> all_spectra;
all_spectra.reserve(42);
all_peaks.reserve(42);

for (int k = 0; k < 42; k++)
{
  Peak1D p1; MSSpectrum s1;
  p1.setIntensity( intensities[k] );
  p1.setMZ(1000);
  s1.setRT( rts[k] );
  all_peaks.push_back(p1); all_spectra.push_back(s1);
  mt1.peaks.push_back(make_pair(&all_spectra.back(), &all_peaks.back()));
  // std::cout << "1_" << k /2 +1 << " :: " << mt1.peaks.back().first->getRT() << " : " << mt1.peaks.back().second->getIntensity() << std::endl;

  k++;
  Peak1D p2; MSSpectrum s2;
  p2.setIntensity(intensities[k]);
  p2.setMZ(1001);
  s2.setRT(rts[k]);
  all_peaks.push_back(p2); all_spectra.push_back(s2);
  mt2.peaks.push_back(make_pair(&all_spectra.back(), &all_peaks.back()));
  // std::cout << "2_" << k /2 +1 << " :: " << mt2.peaks.back().first->getRT() << " : " << mt2.peaks.back().second->getIntensity() << std::endl;
}

mt1.updateMaximum();
mts.push_back(mt1);

mt2.updateMaximum();
mts.push_back(mt2);

// fix base line to 0 since we have no baseline here
mts.baseline = 0.0;

mts.max_trace = 0;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
// setup fitter

Param p;
p.setValue("max_iteration", 500);

EGHTraceFitter egh_trace_fitter;
egh_trace_fitter.setParameters(p);
egh_trace_fitter.fit(mts);

double expected_sigma = 1.5;
double expected_H = 10.0;
double expected_x0 = 680.1;
double expected_tau = 0.0;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

EGHTraceFitter* ptr = nullptr;
EGHTraceFitter* nullPointer = nullptr;
START_SECTION(EGHTraceFitter())
{
  ptr = new EGHTraceFitter();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~EGHTraceFitter())
{
	delete ptr;
}
END_SECTION

START_SECTION((EGHTraceFitter(const EGHTraceFitter& other)))
{
  EGHTraceFitter egh1(egh_trace_fitter);

  TEST_EQUAL(egh1.getCenter(),egh_trace_fitter.getCenter())
  TEST_EQUAL(egh1.getHeight(),egh_trace_fitter.getHeight())
  TEST_EQUAL(egh1.getLowerRTBound(),egh_trace_fitter.getLowerRTBound())
  TEST_EQUAL(egh1.getUpperRTBound(), egh_trace_fitter.getUpperRTBound())
}
END_SECTION

START_SECTION((EGHTraceFitter& operator=(const EGHTraceFitter& source)))
{
  EGHTraceFitter egh1;

  egh1 = egh_trace_fitter;

  TEST_EQUAL(egh1.getCenter(),egh_trace_fitter.getCenter())
  TEST_EQUAL(egh1.getHeight(),egh_trace_fitter.getHeight())
  TEST_EQUAL(egh1.getLowerRTBound(),egh_trace_fitter.getLowerRTBound())
  TEST_EQUAL(egh1.getUpperRTBound(), egh_trace_fitter.getUpperRTBound())
}
END_SECTION

START_SECTION((void fit(FeatureFinderAlgorithmPickedHelperStructs::MassTraces& traces)))
{
  // fit was already done before
  TEST_REAL_SIMILAR(egh_trace_fitter.getCenter(), expected_x0)
  TEST_REAL_SIMILAR(egh_trace_fitter.getHeight(), expected_H)
  EGHTraceFitter weighted_fitter;
  Param params = weighted_fitter.getDefaults();
  params.setValue("weighted", "true");
  weighted_fitter.setParameters(params);
  weighted_fitter.fit(mts);
  TEST_REAL_SIMILAR(weighted_fitter.getCenter(), expected_x0)
  TEST_REAL_SIMILAR(weighted_fitter.getHeight(), expected_H)
  mts[0].theoretical_int = 0.4;
  mts[1].theoretical_int = 0.6;
  weighted_fitter.fit(mts);
  TEST_REAL_SIMILAR(weighted_fitter.getCenter(), expected_x0)
  TEST_REAL_SIMILAR(weighted_fitter.getHeight(), 6.0825)
}
END_SECTION

START_SECTION((double getLowerRTBound() const))
{
  TEST_REAL_SIMILAR(egh_trace_fitter.getLowerRTBound(), expected_x0 - 2.5 * expected_sigma)
}
END_SECTION

START_SECTION((double getUpperRTBound() const))
{
  TEST_REAL_SIMILAR(egh_trace_fitter.getUpperRTBound(), expected_x0 + 2.5 * expected_sigma)
}
END_SECTION

START_SECTION((double getHeight() const))
{
  TEST_REAL_SIMILAR(egh_trace_fitter.getHeight(), expected_H)
}
END_SECTION

START_SECTION((double getCenter() const))
{
  TEST_REAL_SIMILAR(egh_trace_fitter.getCenter(), expected_x0)
}
END_SECTION

START_SECTION((double getTau() const))
{
  TEST_REAL_SIMILAR(egh_trace_fitter.getTau(), expected_tau)
}
END_SECTION

START_SECTION((double getSigma() const))
{
  TEST_REAL_SIMILAR(egh_trace_fitter.getSigma(), expected_sigma)
}
END_SECTION

START_SECTION((double getValue(double rt) const))
{
  TEST_REAL_SIMILAR(egh_trace_fitter.getValue(expected_x0), expected_H)
}
END_SECTION

START_SECTION((double computeTheoretical(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace& trace, Size k)))
{
  FeatureFinderAlgorithmPickedHelperStructs::MassTrace mt;
  mt.theoretical_int = 0.8;

  Peak1D peak;
  peak.setIntensity(8.0);

  MSSpectrum spec;
  spec.setRT(expected_x0);

  mt.peaks.push_back(make_pair(&spec, &peak));

  // theoretical should be expected_H * theoretical_int at position expected_x0
  TEST_REAL_SIMILAR(egh_trace_fitter.computeTheoretical(mt, 0), mt.theoretical_int * expected_H)
}
END_SECTION

START_SECTION((bool checkMaximalRTSpan(const double max_rt_span)))
{
  // Maximum RT span in relation to extended area that the model is allowed to have
  // 5.0 * sigma_ > max_rt_span * region_rt_span_

  double region_rt_span = mt1.peaks[mt1.peaks.size() - 1].first->getRT() - mt1.peaks[0].first->getRT();
  double max_rt_span = 5.0 * expected_sigma/ region_rt_span;

  TEST_EQUAL(egh_trace_fitter.checkMaximalRTSpan(max_rt_span), false);
  max_rt_span -= 0.1; // accept only smaller regions
  TEST_EQUAL(egh_trace_fitter.checkMaximalRTSpan(max_rt_span), true);
}
END_SECTION

START_SECTION((bool checkMinimalRTSpan(const std::pair<double, double>& rt_bounds, const double min_rt_span)))
{
  // is
  // (rt_bounds.second-rt_bounds.first) < min_rt_span * 5.0 * sigma_;
  // Minimum RT span in relation to extended area that has to remain after model fitting.

  pair<double, double> rt_bounds = make_pair(0.0,4.0);
  double min_rt_span = 0.5;

  TEST_EQUAL(egh_trace_fitter.checkMinimalRTSpan(rt_bounds, min_rt_span), false)
  min_rt_span += 0.5;
  TEST_EQUAL(egh_trace_fitter.checkMinimalRTSpan(rt_bounds, min_rt_span), true)
}
END_SECTION

START_SECTION((virtual double getArea()))
{
  TEST_REAL_SIMILAR(egh_trace_fitter.getArea(), sqrt(2 * PI) * expected_sigma * expected_H)
}
END_SECTION

START_SECTION((virtual String getGnuplotFormula(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace& trace, const char function_name, const double baseline, const double rt_shift)))
{
  String formula = egh_trace_fitter.getGnuplotFormula(mts[0], 'f', 0.0, 0.0);
  // should look like -- f(x)= 0 + (((4.5 + 3.93096e-15 * (x - 680.1 )) > 0) ? 8 * exp(-1 * (x - 680.1)**2 / ( 4.5 + 3.93096e-15 * (x - 680.1 ))) : 0) --
  TEST_EQUAL(formula.hasPrefix("f(x)= 0 + ((("), true)
  TEST_EQUAL(formula.hasSubstring(" )) > 0) ? "), true)
  TEST_EQUAL(formula.hasSubstring(" * exp(-1 * ("), true)
  TEST_EQUAL(formula.hasSubstring(")**2 / ( "), true)
  TEST_EQUAL(formula.hasSuffix(" ))) : 0)"), true)
}
END_SECTION

START_SECTION((double getFWHM() const))
{
  TEST_REAL_SIMILAR(egh_trace_fitter.getFWHM(), 3.53223007592464)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



