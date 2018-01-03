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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResultRaw.h>

using namespace OpenMS;

START_TEST(MultiplexFilterResultRaw, "$Id$")

std::vector<double> mz_shifts;
mz_shifts.push_back(0);
mz_shifts.push_back(0.501677);
mz_shifts.push_back(3.01591);
mz_shifts.push_back(3.51759);

std::vector<double> intensities;
intensities.push_back(1789.0714);
intensities.push_back(1492.1012);
intensities.push_back(333.1105);
intensities.push_back(325.0520);

MultiplexFilterResultRaw* nullPointer = nullptr;
MultiplexFilterResultRaw* ptr;

START_SECTION(MultiplexFilterResultRaw(double mz, std::vector<double> mz_shifts, std::vector<double> intensities))
    MultiplexFilterResultRaw result(817.0411, mz_shifts, intensities);
    TEST_EQUAL(result.getMZ(), 817.0411);
    ptr = new MultiplexFilterResultRaw(817.0411, mz_shifts, intensities);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

MultiplexFilterResultRaw result(817.0411, mz_shifts, intensities);

START_SECTION(double getMZ() const)
  TEST_EQUAL(result.getMZ(), 817.0411);
END_SECTION

START_SECTION(std::vector<double> getMZShifts() const)
  TEST_EQUAL(result.getMZShifts()[0], 0);
  TEST_EQUAL(result.getMZShifts()[1], 0.501677);
  TEST_EQUAL(result.getMZShifts()[2], 3.01591);
  TEST_EQUAL(result.getMZShifts()[3], 3.51759);
END_SECTION

START_SECTION(std::vector<double> getIntensities() const)
  TEST_EQUAL(result.getIntensities()[0], 1789.0714);
  TEST_EQUAL(result.getIntensities()[1], 1492.1012);
  TEST_EQUAL(result.getIntensities()[2], 333.1105);
  TEST_EQUAL(result.getIntensities()[3], 325.0520);
END_SECTION

END_TEST
