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
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/Ribonucleotide.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Ribonucleotide, "$Id$")

/////////////////////////////////////////////////////////////

Ribonucleotide* e_ptr = nullptr;
Ribonucleotide* e_nullPointer = nullptr;
START_SECTION((Ribonucleotide()))
  e_ptr = new Ribonucleotide();
  TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((virtual ~Ribonucleotide()))
  delete e_ptr;
END_SECTION

Ribonucleotide test_ribo("test", "T");

START_SECTION(Ribonucleotide(const Ribonucleotide& ribonucleotide))
  Ribonucleotide copy(test_ribo);
  TEST_EQUAL(copy, test_ribo)
END_SECTION

START_SECTION(Ribonucleotide& operator=(const Ribonucleotide& ribonucleotide))
  Ribonucleotide copy;
  copy = test_ribo;
  TEST_EQUAL(copy, test_ribo);
END_SECTION

START_SECTION(bool operator==(const Ribonucleotide& ribonucleotide) const)
  TEST_EQUAL(Ribonucleotide(),
             Ribonucleotide("unknown ribonucleotide", ".", "", ".", EmpiricalFormula(), '.', 0.0, 0.0, Ribonucleotide::ANYWHERE, EmpiricalFormula("C5H10O5")));
END_SECTION

START_SECTION(String getName() const)
  Ribonucleotide empty = Ribonucleotide();
  TEST_EQUAL(empty.getName(), "unknown ribonucleotide");
END_SECTION

START_SECTION(String getCode() const)
  Ribonucleotide empty = Ribonucleotide();
  TEST_EQUAL(empty.getCode(), ".");
END_SECTION

START_SECTION(String getNewCode() const)
  Ribonucleotide empty = Ribonucleotide();
  TEST_EQUAL(empty.getNewCode(), "");
END_SECTION

START_SECTION(String getHTMLCode() const)
  Ribonucleotide empty = Ribonucleotide();
  TEST_EQUAL(empty.getHTMLCode(), ".");
END_SECTION

START_SECTION(EmpiricalFormula getFormula() const)
  Ribonucleotide empty = Ribonucleotide();
  TEST_EQUAL(empty.getFormula(), EmpiricalFormula());
END_SECTION

START_SECTION(char getOrigin() const)
  Ribonucleotide empty = Ribonucleotide();
  TEST_EQUAL(empty.getOrigin(), '.');
END_SECTION

START_SECTION(double getMonoMass() const)
  Ribonucleotide empty = Ribonucleotide();
  TEST_EQUAL(empty.getMonoMass(), 0.0);
END_SECTION

START_SECTION(double getAvgMass() const)
  Ribonucleotide empty = Ribonucleotide();
  TEST_EQUAL(empty.getAvgMass(), 0.0);
END_SECTION

START_SECTION(Ribonucleotide::TermSpecificity getTermSpecificity() const)
  Ribonucleotide empty = Ribonucleotide();
  TEST_EQUAL(empty.getTermSpecificity(), Ribonucleotide::ANYWHERE);
END_SECTION

START_SECTION(EmpiricalFormula getBaselossFormula() const)
  Ribonucleotide empty = Ribonucleotide();
  TEST_EQUAL(empty.getBaselossFormula(), EmpiricalFormula("C5H10O5"));
END_SECTION


START_SECTION(setName(string name))
  Ribonucleotide empty = Ribonucleotide();
  empty.setName("foo");
  TEST_EQUAL(empty.getName(), "foo");
END_SECTION

START_SECTION(setCode(string code))
  Ribonucleotide empty = Ribonucleotide();
  empty.setCode("x");
  TEST_EQUAL(empty.getCode(), "x");
END_SECTION

START_SECTION(setNewCode(String newCode))
  Ribonucleotide empty = Ribonucleotide();
  empty.setNewCode("y");
  TEST_EQUAL(empty.getNewCode(), "y");
END_SECTION

START_SECTION(setHTMLCode(String hmtlCode))
  Ribonucleotide empty = Ribonucleotide();
  empty.setHTMLCode("z");
  TEST_EQUAL(empty.getHTMLCode(), "z");
END_SECTION

START_SECTION(setFormula(EmpiricalFormula formula))
  Ribonucleotide empty = Ribonucleotide();
  empty.setFormula(EmpiricalFormula("H2O"));
  TEST_EQUAL(empty.getFormula(), EmpiricalFormula("H2O"));
END_SECTION

START_SECTION(setOrigin(char origin))
  Ribonucleotide empty = Ribonucleotide();
  empty.setOrigin('q');
  TEST_EQUAL(empty.getOrigin(), 'q');
END_SECTION

START_SECTION(setMonoMass(double mono_mass))
  Ribonucleotide empty = Ribonucleotide();
  empty.setMonoMass(2.0);
  TEST_EQUAL(empty.getMonoMass(), 2.0);
END_SECTION

START_SECTION(setAvgMass(double avg_mass))
  Ribonucleotide empty = Ribonucleotide();
  empty.setAvgMass(3.0);
  TEST_EQUAL(empty.getAvgMass(), 3.0);
END_SECTION

START_SECTION(setTermSpecificity(Ribonucleotide::TermSpecificity specificity))
  Ribonucleotide empty = Ribonucleotide();
  empty.setTermSpecificity(Ribonucleotide::FIVE_PRIME);
  TEST_EQUAL(empty.getTermSpecificity(), Ribonucleotide::FIVE_PRIME);
END_SECTION

START_SECTION(void setBaselossFormula(EmpiricalFormula formula))
  Ribonucleotide empty = Ribonucleotide();
  empty.setBaselossFormula(EmpiricalFormula("CO2"));
  TEST_EQUAL(empty.getBaselossFormula(), EmpiricalFormula("CO2"));
END_SECTION

START_SECTION(bool isModified() const)
  TEST_EQUAL(test_ribo.isModified(), true);
  test_ribo.setOrigin('T');
  TEST_EQUAL(test_ribo.isModified(), false);
  test_ribo.setCode("Tm");
  TEST_EQUAL(test_ribo.isModified(), true);
END_SECTION

END_TEST

