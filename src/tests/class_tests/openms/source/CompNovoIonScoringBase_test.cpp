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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIonScoringBase.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(CompNovoIonScoringBase, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(CompNovoIonScoringBase())
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION((CompNovoIonScoringBase(const CompNovoIonScoringBase &source)))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual ~CompNovoIonScoringBase()))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION((double scoreIsotopes(const PeakSpectrum &CID_spec, PeakSpectrum::ConstIterator it, Size charge)))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION((CompNovoIonScoringBase& operator=(const CompNovoIonScoringBase &source)))
{
	NOT_TESTABLE
}
END_SECTION

CompNovoIonScoringBase::IonScore * ptr = nullptr;
CompNovoIonScoringBase::IonScore * nullPointer = nullptr;
START_SECTION([CompNovoIonScoringBase::IonScore] IonScore())
	ptr=new CompNovoIonScoringBase::IonScore();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION([CompNovoIonScoringBase::IonScore] IonScore(const IonScore &rhs))
	CompNovoIonScoringBase::IonScore ion_score;
	ion_score.s_bion=5.0;
	TEST_EQUAL(CompNovoIonScoringBase::IonScore(ion_score).s_bion, 5.0)
END_SECTION

START_SECTION([CompNovoIonScoringBase::IonScore] IonScore& operator=(const IonScore &rhs))
	CompNovoIonScoringBase::IonScore ion_score, copy;
	ion_score.s_bion=5.0;
	copy = ion_score;
	TEST_EQUAL(copy.s_bion, ion_score.s_bion)
END_SECTION

START_SECTION([CompNovoIonScoringBase::IonScore] virtual ~IonScore())
	delete ptr;
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



