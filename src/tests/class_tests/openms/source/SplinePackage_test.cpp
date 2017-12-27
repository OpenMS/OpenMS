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

#include <OpenMS/FILTERING/DATAREDUCTION/SplinePackage.h>

using namespace OpenMS;

START_TEST(SplinePackage, "$Id$")

std::vector<double> mz;
mz.push_back(413.8);
mz.push_back(413.9);
mz.push_back(414.0);
mz.push_back(414.1);
mz.push_back(414.2);
std::vector<double> intensity;
intensity.push_back(0.0);
intensity.push_back(100.2);
intensity.push_back(20.3);
intensity.push_back(2000.4);
intensity.push_back(4.3);

std::vector<double> mz1;
mz1.push_back(413.9);
std::vector<double> intensity1;
intensity1.push_back(100.2);

std::vector<double> mz2;
mz2.push_back(413.8);
mz2.push_back(413.9);
std::vector<double> intensity2;
intensity2.push_back(0.0);
intensity2.push_back(100.2);

SplinePackage sp1(mz, intensity, 0.7);

SplinePackage* nullPointer = nullptr;

START_SECTION(SplinePackage(std::vector<double> mz, std::vector<double> intensity, double scaling))
  SplinePackage* sp2 = new SplinePackage(mz, intensity, 0.7);
  TEST_NOT_EQUAL(sp2, nullPointer)
END_SECTION

START_SECTION(getMzMin())
  TEST_EQUAL(sp1.getMzMin(), 413.8);
END_SECTION

START_SECTION(getMzMax())
  TEST_EQUAL(sp1.getMzMax(), 414.2);
END_SECTION

START_SECTION(getMzStepWidth())
  TEST_REAL_SIMILAR(sp1.getMzStepWidth(), 0.07);
END_SECTION

START_SECTION(isInPackage(double mz))
  TEST_EQUAL(sp1.isInPackage(414.05), true);
END_SECTION

START_SECTION(eval(double mz))
  TEST_REAL_SIMILAR(sp1.eval(414.05), 1134.08593750018);
END_SECTION

START_SECTION(SplinePackage(std::vector<double> mz, std::vector<double> intensity, double scaling))
  TEST_EXCEPTION(Exception::IllegalArgument, SplinePackage(mz1, intensity1, 0.7));
END_SECTION

START_SECTION(SplinePackage(std::vector<double> mz, std::vector<double> intensity, double scaling))
  SplinePackage* sp4 = new SplinePackage(mz2, intensity2, 0.7);
  TEST_NOT_EQUAL(sp4, nullPointer);
  TEST_REAL_SIMILAR((*sp4).eval(413.85), 50.1);
END_SECTION

END_TEST
