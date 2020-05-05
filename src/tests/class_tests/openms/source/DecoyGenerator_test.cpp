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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/DecoyGenerator.h>

///////////////////////////

START_TEST(DecoyGenerator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

DecoyGenerator* dg = nullptr;
DecoyGenerator* nullPointer = nullptr;

START_SECTION((DecoyGenerator()))
{
  dg = new DecoyGenerator();
  TEST_NOT_EQUAL(dg, nullPointer)
}
END_SECTION

START_SECTION((~DecoyGenerator()))
{
  delete dg;
}
END_SECTION

dg = new DecoyGenerator();
dg->setSeed(4711);

START_SECTION((AASequence reverseProtein(const AASequence& protein)))
  TEST_EQUAL(dg->reverseProtein(AASequence::fromString("PRTEINE")).toString(), "ENIETRP")
END_SECTION

START_SECTION((AASequence reversePeptide(const AASequence& protein, const String& protease)))
  TEST_EQUAL(dg->reversePeptides(AASequence::fromString("TESTPEPTIDE"), "Trypsin").toString(),"EDITPEPTSET")
  TEST_EQUAL(dg->reversePeptides(AASequence::fromString("TESTRPEPTRIDE"), "Trypsin/P").toString(),"TSETRTPEPREDI")
  TEST_EQUAL(dg->reversePeptides(AASequence::fromString("TESTRPEPTRIDE"), "Trypsin").toString(),"TPEPRTSETREDI")
END_SECTION

START_SECTION((AASequence shufflePeptides(const AASequence& aas, const String& protease, const int max_atempts, int seed)))
  TEST_EQUAL(dg->shufflePeptides(AASequence::fromString("TESTPEPTIDE"), "Trypsin").toString(),"PIDPETTSEET")
  TEST_EQUAL(dg->shufflePeptides(AASequence::fromString("TESTRPEPTRIDE"), "Trypsin/P").toString(),"ETTSRTPEPREID")
  TEST_EQUAL(dg->shufflePeptides(AASequence::fromString("TESTRPEPTRIDE"), "Trypsin").toString(), "ETPSERTTPREID")
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
