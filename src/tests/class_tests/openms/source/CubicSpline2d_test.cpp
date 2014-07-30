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

#include <OpenMS/MATH/MISC/CubicSpline2d.h>

using namespace OpenMS;

START_TEST(CubicSpline2d, "$Id$")

std::vector<double> mz;
mz.push_back(486.784);
mz.push_back(486.787);
mz.push_back(486.790);
mz.push_back(486.793);
mz.push_back(486.795);
mz.push_back(486.797);
mz.push_back(486.800);
mz.push_back(486.802);
mz.push_back(486.805);
mz.push_back(486.808);
mz.push_back(486.811);
std::vector<double> intensity;
intensity.push_back(0.0);
intensity.push_back(154683.17);
intensity.push_back(620386.5);
intensity.push_back(1701390.12);
intensity.push_back(2848879.25);
intensity.push_back(3564045.5);
intensity.push_back(2744585.7);
intensity.push_back(1605583.0);
intensity.push_back(518984.0);
intensity.push_back(59152.21);
intensity.push_back(0.0);

std::map<double,double> map;
for (Size i=0; i<mz.size(); ++i)
{
    map.insert(std::pair<double,double>(mz[i], intensity[i]));
}

CubicSpline2d sp1(mz, intensity);
CubicSpline2d sp2(map);

CubicSpline2d* nullPointer = 0;

START_SECTION(CubicSpline2d(const std::vector<double>& x, const std::vector<double>& y))
  CubicSpline2d* sp3 = new CubicSpline2d(mz, intensity);
  TEST_NOT_EQUAL(sp3, nullPointer)
END_SECTION

START_SECTION(CubicSpline2d(const std::map<double, double>& m))
  CubicSpline2d* sp4 = new CubicSpline2d(map);
  TEST_NOT_EQUAL(sp4, nullPointer)
END_SECTION

START_SECTION(double eval(double x))
  // near border of spline range
  TEST_REAL_SIMILAR(sp1.eval(486.785), 35203.124211885166);
  TEST_REAL_SIMILAR(sp2.eval(486.785), 35203.124211885166);
  // inside spline range
  TEST_REAL_SIMILAR(sp1.eval(486.794), 2270517.50463);
  TEST_REAL_SIMILAR(sp2.eval(486.794), 2270517.50463);
  // at the input nodes
  TEST_REAL_SIMILAR(sp1.eval(486.784), 0.0);
  TEST_REAL_SIMILAR(sp1.eval(486.790), 620386.5);
  TEST_REAL_SIMILAR(sp2.eval(486.808), 59152.21);
  TEST_REAL_SIMILAR(sp2.eval(486.811), 0.0);
END_SECTION

START_SECTION(double derivatives(double x, unsigned order))
  // near border of spline range
  TEST_REAL_SIMILAR(sp1.derivatives(486.785,1), 39292607.3251133)
  TEST_REAL_SIMILAR(sp1.derivatives(486.785,2), 12268449342.7831);
  TEST_REAL_SIMILAR(sp2.derivatives(486.785,1), 39292607.3251133);
  TEST_REAL_SIMILAR(sp2.derivatives(486.785,2), 12268449342.7831);
  // inside spline range
  TEST_REAL_SIMILAR(sp1.derivatives(486.794,1), 594354391.618931);
  TEST_REAL_SIMILAR(sp1.derivatives(486.794,2), 9234360709.46161);
  TEST_REAL_SIMILAR(sp2.derivatives(486.794,1), 594354391.618931);
  TEST_REAL_SIMILAR(sp2.derivatives(486.794,2), 9234360709.46161);
END_SECTION

END_TEST
