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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMasses.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>

using namespace OpenMS;

START_TEST(MultiplexIsotopicPeakPattern, "$Id$")

MultiplexDeltaMasses mass_shifts;
mass_shifts.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0,"no_label"));
mass_shifts.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(6.031817,"Arg6"));

MultiplexIsotopicPeakPattern* nullPointer = nullptr;
MultiplexIsotopicPeakPattern* ptr;

START_SECTION(MultiplexIsotopicPeakPattern(int c, int ppp, MultiplexDeltaMasses ms, int msi))
    MultiplexIsotopicPeakPattern pattern(2, 4, mass_shifts, 3);
    TEST_EQUAL(pattern.getCharge(), 2);
    ptr = new MultiplexIsotopicPeakPattern(2, 4, mass_shifts, 3);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

MultiplexIsotopicPeakPattern pattern(2, 4, mass_shifts, 3);

START_SECTION(int getCharge() const)
  TEST_EQUAL(pattern.getCharge(), 2);
END_SECTION

START_SECTION(int getPeaksPerPeptide() const)
  TEST_EQUAL(pattern.getPeaksPerPeptide(), 4);
END_SECTION

START_SECTION(std::vector<double> getMassShifts() const)
  TEST_EQUAL(pattern.getMassShifts().getDeltaMasses()[0].delta_mass, 0);
  TEST_EQUAL(pattern.getMassShifts().getDeltaMasses()[1].delta_mass, 6.031817);
END_SECTION

START_SECTION(int getMassShiftIndex() const)
  TEST_EQUAL(pattern.getMassShiftIndex(), 3);
END_SECTION

START_SECTION(unsigned getMassShiftCount() const)
  TEST_EQUAL(pattern.getMassShiftCount(), 2);
END_SECTION

START_SECTION(double getMassShiftAt(int i) const)
  TEST_EQUAL(pattern.getMassShiftAt(0), 0);
  TEST_EQUAL(pattern.getMassShiftAt(1), 6.031817);
END_SECTION

START_SECTION(double getMZShiftAt(int i) const)
  TEST_REAL_SIMILAR(pattern.getMZShiftAt(0), -0.501677);
  TEST_REAL_SIMILAR(pattern.getMZShiftAt(1), 0);
  TEST_REAL_SIMILAR(pattern.getMZShiftAt(2), 0.501677);
  TEST_REAL_SIMILAR(pattern.getMZShiftAt(3), 1.00335);
  TEST_REAL_SIMILAR(pattern.getMZShiftAt(4), 1.50503);
  TEST_REAL_SIMILAR(pattern.getMZShiftAt(5), 2.51423);
  TEST_REAL_SIMILAR(pattern.getMZShiftAt(6), 3.01591);
  TEST_REAL_SIMILAR(pattern.getMZShiftAt(7), 3.51759);
  TEST_REAL_SIMILAR(pattern.getMZShiftAt(8), 4.01926);
  TEST_REAL_SIMILAR(pattern.getMZShiftAt(9), 4.52094);
END_SECTION

START_SECTION(unsigned getMZShiftCount() const)
  TEST_EQUAL(pattern.getMZShiftCount(), 10);
END_SECTION

END_TEST
