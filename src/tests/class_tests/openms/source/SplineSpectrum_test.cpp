// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FILTERING/DATAREDUCTION/SplineSpectrum.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;

double Gauss1(double x)
{
    return exp(- pow(x-416.8, 2)/(2*0.15*0.15));
}

double Gauss2(double x)
{
    return exp(- pow(x-418.7, 2)/(2*0.15*0.15));
}

START_TEST(SplineSpectrum, "$Id$")

std::vector<double> mz;
std::vector<double> intensity;
for (int i=0; i < 11; ++i)
{
    mz.push_back(416.3 + 0.1*i);
    intensity.push_back(Gauss1(416.3+0.1*i));
}
for (int i=0; i < 11; ++i)
{
    mz.push_back(418.2 + 0.1*i);
    intensity.push_back(Gauss2(418.2+0.1*i));
}

MSSpectrum<Peak1D> spectrum;
Peak1D peak;
spectrum.setRT(1789.0714);
for (unsigned i=0; i < mz.size(); ++i)
{
    peak.setMZ(mz[i]);
    peak.setIntensity(intensity[i]);
    spectrum.push_back(peak);
}

SplineSpectrum* nullPointer = 0;
SplineSpectrum* ptr;

START_SECTION(SplineSpectrum(const std::vector<double>& mz, const std::vector<double>& intensity))
    SplineSpectrum spline(mz, intensity);
    TEST_REAL_SIMILAR(spline.getMzMin(), 416.3);
    ptr = new SplineSpectrum(mz, intensity);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

START_SECTION(SplineSpectrum(const std::vector<double>& mz, const std::vector<double>& intensity, double scaling))
    SplineSpectrum spline(mz, intensity, 0.7);
    TEST_REAL_SIMILAR(spline.getMzMin(), 416.3)
    ptr = new SplineSpectrum(mz, intensity, 0.7);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

START_SECTION(SplineSpectrum(MSSpectrum<Peak1D>& raw_spectrum))
	SplineSpectrum spline(spectrum);
    TEST_REAL_SIMILAR(spline.getMzMin(), 416.3)
    ptr = new SplineSpectrum(spectrum);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

START_SECTION(SplineSpectrum(MSSpectrum<Peak1D>& raw_spectrum, double scaling))
	SplineSpectrum spline(spectrum, 0.7);
    TEST_REAL_SIMILAR(spline.getMzMin(), 416.3)
    ptr = new SplineSpectrum(spectrum, 0.7);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

SplineSpectrum spectrum2(mz, intensity);

START_SECTION(double getMzMin() const)
  TEST_EQUAL(spectrum2.getMzMin(), 416.3);
END_SECTION

START_SECTION(double getMzMax() const)
  TEST_EQUAL(spectrum2.getMzMax(), 419.2);
END_SECTION

MSSpectrum<> empty_spec;
SplineSpectrum ss_empty(empty_spec);

START_SECTION(unsigned getSplineCount() const)
  TEST_EQUAL(spectrum2.getSplineCount(), 2)
  
  // this should be used before getNavigator()
  TEST_EQUAL(ss_empty.getSplineCount(), 0)
END_SECTION

START_SECTION(SplineSpectrum::Navigator getNavigator())
  // just to test if it can be called
  SplineSpectrum::Navigator nav = spectrum2.getNavigator();

  // test exception on empty spectrum
  TEST_EXCEPTION(Exception::InvalidSize, ss_empty.getNavigator())

END_SECTION

START_SECTION(double SplineSpectrum::Navigator::eval(double mz))
  // outside range of Gaussians
  TEST_EQUAL(spectrum2.getNavigator().eval(400.0), 0);
  TEST_EQUAL(spectrum2.getNavigator().eval(417.8), 0);
  TEST_EQUAL(spectrum2.getNavigator().eval(500.0), 0);
  // near the edge
  TEST_REAL_SIMILAR(spectrum2.getNavigator().eval(416.33), 0.007848195698809);  // expected 0.00738068453767004 differs by 6%
  // near the maximum
  TEST_REAL_SIMILAR(spectrum2.getNavigator().eval(416.81), 0.997572728799559);  // expected 0.99778024508561 differs by 0.02%
  // evaluation in first package, then search in last package
  SplineSpectrum::Navigator nav = spectrum2.getNavigator();
  TEST_REAL_SIMILAR(nav.eval(416.81), 0.997572728799559);
  TEST_REAL_SIMILAR(nav.eval(418.75), 0.944147611428987);
  // evaluation in last package, then search in first package
  SplineSpectrum::Navigator nav2 = spectrum2.getNavigator();
  TEST_REAL_SIMILAR(nav2.eval(418.75), 0.944147611428987);
  TEST_REAL_SIMILAR(nav2.eval(416.81), 0.997572728799559);
END_SECTION

START_SECTION(double SplineSpectrum::Navigator::getNextMz(double mz))
  // advancing within package
  TEST_EQUAL(spectrum2.getNavigator().getNextMz(417.0), 417.07);
  // advancing to next package
  TEST_EQUAL(spectrum2.getNavigator().getNextMz(417.29), 418.2);
  // advancing beyond range
  TEST_REAL_SIMILAR(spectrum2.getNavigator().getNextMz(500.0), 419.2);
END_SECTION

std::vector<double> mz_short;
std::vector<double> intensity_short;
mz_short.push_back(416.3);
mz_short.push_back(416.5);
intensity_short.push_back(6.43);
intensity_short.push_back(8.88);

START_SECTION(SplineSpectrum(const std::vector<double>& mz, const std::vector<double>& intensity))
  // For 4 or fewer data points, a single SplinePackage is always generated disregarding the m/z spacing.
  SplineSpectrum spline_2datapoints(mz_short, intensity_short);
  TEST_EQUAL(spline_2datapoints.getSplineCount(), 1);
END_SECTION

mz_short.push_back(616.5);
mz_short.push_back(616.7);
intensity_short.push_back(398.41);
intensity_short.push_back(8.02);

START_SECTION(SplineSpectrum(const std::vector<double>& mz, const std::vector<double>& intensity))
  // For 4 or fewer data points, a single SplinePackage is always generated disregarding the m/z spacing.
  SplineSpectrum spline_4datapoints(mz_short, intensity_short);
  TEST_EQUAL(spline_4datapoints.getSplineCount(), 1);
END_SECTION

mz_short.push_back(616.95);
intensity_short.push_back(7.71);

START_SECTION(SplineSpectrum(const std::vector<double>& mz, const std::vector<double>& intensity))
  // For 5 or more data points, the data points can be split into multiple SplinePackages.
  SplineSpectrum spline_5datapoints(mz_short, intensity_short);
  TEST_EQUAL(spline_5datapoints.getSplineCount(), 2);
END_SECTION

END_TEST
