// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexMassPatternList.h>

using namespace OpenMS;

START_TEST(MultiplexMassPatternList, "$Id$")

std::vector<double> mass_shifts;
mass_shifts.push_back(0);
mass_shifts.push_back(6.031817);

MultiplexMassPatternList* nullPointer = 0;
MultiplexMassPatternList* ptr;

START_SECTION(MultiplexMassPatternList(std::vector<double> ms))
    MultiplexMassPatternList pattern(mass_shifts);
    TEST_EQUAL(pattern.getMassShiftCount(), 2);
    ptr = new MultiplexMassPatternList(mass_shifts);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

MultiplexMassPatternList pattern(mass_shifts);

START_SECTION(void addMassShifts(double ms) const)
  pattern.addMassShift(12.063634);
  TEST_EQUAL(pattern.getMassShifts()[2], 12.063634);
END_SECTION

START_SECTION(std::vector<double> getMassShifts() const)
  TEST_EQUAL(pattern.getMassShifts()[0], 0);
  TEST_EQUAL(pattern.getMassShifts()[1], 6.031817);
END_SECTION

START_SECTION(unsigned getMassShiftCount() const)
  TEST_EQUAL(pattern.getMassShiftCount(), 3);
END_SECTION

START_SECTION(double getMassShiftAt(int i) const)
  TEST_EQUAL(pattern.getMassShiftAt(0), 0);
  TEST_EQUAL(pattern.getMassShiftAt(1), 6.031817);
  TEST_EQUAL(pattern.getMassShiftAt(2), 12.063634);
END_SECTION

END_TEST
