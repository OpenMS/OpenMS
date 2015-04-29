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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/NASequence.h>
#include <iostream>
#include <OpenMS/SYSTEM/StopWatch.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(NASequence, "$Id$")

/////////////////////////////////////////////////////////////

NASequence empty = NASequence();
NASequence constructedempty = NASequence("");
START_SECTION(NASequence())
  TEST_EQUAL(empty, constructedempty)
END_SECTION


START_SECTION(NASequence(const String& rhs))
  NASequence seq = NASequence();
  seq.setSequence("AAA");
  NASequence seq2("AAA");
  TEST_EQUAL(seq,seq2)
END_SECTION

START_SECTION(NASequence getSequence())
{
  NASequence seq = NASequence("UUU");
  TEST_EQUAL(seq.getSequence()=="UUU",TRUE);
  TEST_NOT_EQUAL(seq.getSequence(),"UU");
  TEST_NOT_EQUAL(seq.getSequence(),"AAA");
}
END_SECTION

START_SECTION(NASequence setSequence(const String & s)){
    NASequence seq;
    seq.setSequence("GU");
    TEST_EQUAL(seq.getSequence(),"GU");
}
END_SECTION

START_SECTION(NASequence size()){
    NASequence seq;
    TEST_EQUAL(seq.size(),0);
    seq.setSequence("UGG");
    TEST_EQUAL(seq.size(),3);
}
END_SECTION


START_SECTION(NASequence getFormula(Residue::ResidueType type,Int charge)){
    NASequence seq= NASequence("GG");
    TEST_EQUAL(seq.getFormula(Residue::NTerminal, 1),EmpiricalFormula("C10H12N5O4")+EmpiricalFormula("C10H12N5O7P"));
    TEST_EQUAL(seq.getFormula(Residue::CTerminal, 1),EmpiricalFormula("C10H12N5O7P")+EmpiricalFormula("C10H12N5O7P"));
    TEST_EQUAL(seq.getFormula(Residue::Internal, 1),EmpiricalFormula("C10H10N5O6P")+EmpiricalFormula("C10H12N5O7P"));
    TEST_EQUAL(seq.getFormula(Residue::Full, 1),EmpiricalFormula("C10H10N5O6P")+EmpiricalFormula("C10H12N5O6"));

}
END_SECTION
//START_SECTION(NASequence getMonoWeight(Residue::ResidueType type,Int charge)){

//}
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
