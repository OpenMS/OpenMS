// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteredPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteredMSExperiment.h>

using namespace OpenMS;

START_TEST(MultiplexFilteredMSExperiment, "$Id$")

MultiplexFilteredMSExperiment* nullPointer = nullptr;
MultiplexFilteredMSExperiment* ptr;

START_SECTION(MultiplexFilteredMSExperiment())
    MultiplexFilteredMSExperiment exp;
    TEST_EQUAL(exp.size(), 0);
    ptr = new MultiplexFilteredMSExperiment();
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

MultiplexFilteredMSExperiment exp;
MultiplexFilteredPeak peak(654.32, 2345.67, 24, 110);
exp.addPeak(peak);
size_t n;

START_SECTION(addPeak(const MultiplexFilteredPeak& peak))
  n = exp.size();
  MultiplexFilteredPeak peak_temp(655.32, 2346.67, 25, 111);
  exp.addPeak(peak_temp);
  TEST_EQUAL(exp.size(), n + 1);
END_SECTION

START_SECTION(MultiplexFilteredPeak getPeak(size_t i))
  MultiplexFilteredPeak peak = exp.getPeak(0);
  TEST_REAL_SIMILAR(peak.getMZ(), 654.32);
END_SECTION

START_SECTION(double getMZ(size_t i))
  TEST_REAL_SIMILAR(exp.getMZ(0), 654.32);
END_SECTION

START_SECTION(std::vector<double> getMZ())
  TEST_REAL_SIMILAR(exp.getMZ()[0], 654.32);
END_SECTION

START_SECTION(double getRT(size_t i))
  TEST_REAL_SIMILAR(exp.getRT(0), 2345.67);
END_SECTION

START_SECTION(std::vector<double> getRT())
  TEST_REAL_SIMILAR(exp.getRT()[0], 2345.67);
END_SECTION

START_SECTION(size_t size())
  TEST_EQUAL(exp.size(), 2);
END_SECTION

END_TEST
