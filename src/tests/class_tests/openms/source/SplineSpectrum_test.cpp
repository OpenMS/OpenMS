// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

MSSpectrum spectrum;
Peak1D peak;
spectrum.setRT(1789.0714);
for (size_t i=0; i < mz.size(); ++i)
{
    peak.setMZ(mz[i]);
    peak.setIntensity(intensity[i]);
    spectrum.push_back(peak);
}

SplineSpectrum* nullPointer = nullptr;
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

START_SECTION(SplineSpectrum(MSSpectrum& raw_spectrum))
	SplineSpectrum spline(spectrum);
    TEST_REAL_SIMILAR(spline.getMzMin(), 416.3)
    ptr = new SplineSpectrum(spectrum);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

START_SECTION(SplineSpectrum(MSSpectrum& raw_spectrum, double scaling))
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

START_SECTION(size_t getSplineCount() const)
  TEST_EQUAL(spectrum2.getSplineCount(), 2)
END_SECTION

START_SECTION(SplineSpectrum::Navigator getNavigator())
  // just to test if it can be called
  SplineSpectrum::Navigator nav = spectrum2.getNavigator();
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

// Each SplinePackage in a SplineSpectrum must contain two or more data points. If this is not the case, the interpolation might lead to unexpected results.
// In the example below, a single data point @ 407.5 is placed between two packages. It does not form a SplinePackage on its own, but is instead part of the second SplinePackage.
std::vector<double> mz3;
std::vector<double> intensity3;
for (size_t i=0; i<4; ++i)
{
    mz3.push_back(400+i*0.5);
    intensity3.push_back(10.0);
}
mz3.push_back(407.5);
intensity3.push_back(10.0);
for (size_t i=0; i<4; ++i)
{
    mz3.push_back(410+i*0.5);
    intensity3.push_back(10.0);
}
SplineSpectrum spectrum3(mz3, intensity3);

START_SECTION(double SplineSpectrum::Navigator::eval(double mz))
  TEST_EQUAL(spectrum3.getSplineCount(),2);
  TEST_EQUAL(spectrum3.getNavigator().eval(405),0);    // Zero as expected, since 405 is between packages.
  TEST_EQUAL(spectrum3.getNavigator().eval(408),10);    // One might expect zero, but 407.5 is part of the second package.
END_SECTION

std::vector<double> mz4;
std::vector<double> intensity4;
mz4.push_back(407.5);
intensity4.push_back(10.0);
START_SECTION(SplineSpectrum(const std::vector<double>& mz, const std::vector<double>& intensity))
  TEST_EXCEPTION(Exception::IllegalArgument, new SplineSpectrum(mz4,intensity4));
END_SECTION

END_TEST
