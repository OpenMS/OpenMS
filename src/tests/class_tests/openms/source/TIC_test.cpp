// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Chris Bielow$
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/QC/TIC.h>

///////////////////////////

START_TEST(TIC, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
using namespace OpenMS;
using namespace std;

TIC* ptr = nullptr;
TIC* nullPointer = nullptr;
START_SECTION(TIC())
  ptr = new TIC();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~TIC())
  delete ptr;
END_SECTION

// Dummy peakmap
PeakMap exp;
exp.resize(4);
Peak1D p;

// MS spectrum at RT = 0
p.setMZ(5.0);
p.setIntensity(3);
exp[0].push_back(p);
p.setMZ(10.0);
p.setIntensity(5);
exp[0].push_back(p);
exp[0].setMSLevel(1);
exp[0].setRT(0);

// MS spectrum at RT = 2
p.setMZ(5.0);
p.setIntensity(2);
exp[1].push_back(p);
exp[1].setMSLevel(1);
exp[1].setRT(2);

// MSMS spectrum at RT = 2
p.setMZ(5.0);
p.setIntensity(0.5);
exp[2].push_back(p);
exp[2].setMSLevel(2);
exp[2].setRT(2);

// MS spectrum at RT = 5
p.setMZ(5.0);
p.setIntensity(2.0);
exp[3].push_back(p);
p.setMZ(10.0);
p.setIntensity(3.0);
exp[3].push_back(p);
p.setMZ(15.0);
p.setIntensity(4.0);
exp[3].push_back(p);
exp[3].setMSLevel(1);
exp[3].setRT(5);

exp.updateRanges();

START_SECTION(QCBase::Status requires() const)
  TIC tic;
  TEST_EQUAL((tic.requires() == QCBase::Status(QCBase::Requires::RAWMZML)),true);
END_SECTION

START_SECTION(void compute(const MSExperiment &exp, float bin_size) && vector<MSChromatogram> getResults() const)
{
  TIC tic;
  TEST_EQUAL(tic.getResults().empty(),true);
  tic.compute(exp); // no binning
  tic.compute(exp,2.0); // bin size smaller than highest RT
  tic.compute(exp,6.0); // bin size bigger than highest RT
  tic.compute(exp,-1.0); // negative bin size

  TEST_EQUAL(tic.getResults().size(),4);
  TEST_EQUAL(tic.getResults()[0].size(),3);
  TEST_EQUAL(tic.getResults()[0][0].getIntensity(),8);
  TEST_EQUAL(tic.getResults()[0][1].getIntensity(),2);
  TEST_EQUAL(tic.getResults()[0][2].getIntensity(),9);

  TEST_EQUAL(tic.getResults()[1].size(),4);
  TEST_EQUAL(tic.getResults()[1][0].getIntensity(),8);
  TEST_EQUAL(tic.getResults()[1][1].getIntensity(),2);
  // Intensity at RT = 5 in between new data points at 4.0 and 6.0
  TEST_EQUAL(tic.getResults()[1][2].getIntensity(),4.5);
  TEST_EQUAL(tic.getResults()[1][3].getIntensity(),4.5);

  TEST_EQUAL(tic.getResults()[2].size(),2);
  // Intensities at RT = 2 and RT = 5 in between new data points at 0.0 and 6.0
  TEST_REAL_SIMILAR(tic.getResults()[2][0].getIntensity(),8.0 + 2.0* 4.0/6.0 + 9 * 1.0/6.0);
  TEST_REAL_SIMILAR(tic.getResults()[2][1].getIntensity(),2.0* 2.0/6.0 + 9 * 5.0/6.0);

  // same as in tic.getResults()[0]
  TEST_EQUAL(tic.getResults()[3].size(),3);
  TEST_EQUAL(tic.getResults()[3][0].getIntensity(),8);
  TEST_EQUAL(tic.getResults()[3][1].getIntensity(),2);
  TEST_EQUAL(tic.getResults()[3][2].getIntensity(),9);
}
END_SECTION

START_SECTION(void clear())
  TIC tic;
  tic.compute(exp);
  tic.clear();
  TEST_EQUAL(tic.getResults().empty(),true);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
