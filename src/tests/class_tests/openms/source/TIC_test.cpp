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

TIC tic;
START_SECTION(const String& getName() const override)
  TEST_EQUAL(tic.getName(), "TIC")
END_SECTION

START_SECTION(Status requires() const override)
  TEST_EQUAL((tic.requires() == QCBase::Status(QCBase::Requires::RAWMZML)),true);
END_SECTION

START_SECTION(void compute(const MSExperiment &exp, float bin_size))
  // very simple test ATM, since the computation is simply exp.getTIC(bin_size);
  MSExperiment exp;
  exp.setSpectra( { MSSpectrum() });
  TIC tic;
  tic.compute(exp, 0);
  auto r = tic.getResults();
  TEST_EQUAL(r.size(), 1);
  ABORT_IF(r[0].size() != 1); // one intensity per input spectrum
  ABORT_IF(r[0][0].getIntensity() != 0); // empty spectrum
END_SECTION

START_SECTION(vector<MSChromatogram> getResults() const)
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(void clear())
  TIC tic;
  MSExperiment exp2;
  tic.compute(exp2);
  TEST_EQUAL(tic.getResults().empty(), false);
  tic.clear();
  TEST_EQUAL(tic.getResults().empty(), true);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
