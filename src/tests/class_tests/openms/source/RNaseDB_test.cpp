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
// $Maintainer: Xiao Liang $
// $Authors: Xiao Liang, Chris Bielow $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/RNaseDB.h>
#include <OpenMS/CHEMISTRY/DigestionEnzymeRNA.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(RNaseDB, "$Id$")

/////////////////////////////////////////////////////////////

RNaseDB* ptr = nullptr;
RNaseDB* nullPointer = nullptr;
String t1_regex("(?<=G)");

START_SECTION([EXTRA] multithreaded example)
{

  int nr_iterations (1e2), test (0);
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int k = 1; k < nr_iterations + 1; k++)
  {
    auto p = RNaseDB::getInstance();
    int tmp (0);
    if (p->hasEnzyme("Trypsin"), true)
    {
      tmp++;
    }

#ifdef _OPENMP
#pragma omp critical (add_test)
#endif
    {
      test += tmp;
    }
  }
  TEST_EQUAL(test, nr_iterations)
}
END_SECTION

START_SECTION(RNaseDB* getInstance())
    ptr = RNaseDB::getInstance();
    TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(virtual ~RNaseDB())
    NOT_TESTABLE
END_SECTION

START_SECTION((bool hasEnzyme(const String &name) const))
    TEST_EQUAL(ptr->hasEnzyme("RNAse"), false)
    TEST_EQUAL(ptr->hasEnzyme("RNase_T1"), true)
END_SECTION

START_SECTION((const DigestionEnzymeProtein* getEnzyme(const String &name) const))
    TEST_EQUAL(ptr->getEnzyme("RNase_T1")->getName(), "RNase_T1")
END_SECTION

START_SECTION((bool hasRegEx(const String& cleavage_regex) const))
    TEST_EQUAL(ptr->hasRegEx("(?<=[P])(?!P)"), false)
    TEST_EQUAL(ptr->hasRegEx(t1_regex), true)
END_SECTION

START_SECTION((const DigestionEnzymeRNA* getEnzymeByRegEx(const String& cleavage_regex) const))
    TEST_EQUAL(ptr->getEnzymeByRegEx(t1_regex)->getName(), "RNase_T1")
END_SECTION

START_SECTION(bool hasEnzyme(const DigestionEnzymeProtein* enzyme) const)
  TEST_EQUAL(ptr->hasEnzyme(ptr->getEnzyme("RNase_T1")), true)
  DigestionEnzymeRNA myNewEnzyme("bla", "blubb");
  TEST_EQUAL(ptr->hasEnzyme(&myNewEnzyme), false);
END_SECTION

START_SECTION(ConstEnzymeIterator beginEnzyme() const)
    auto it = ptr->beginEnzyme();
    Size count(0);
    while (it != ptr->endEnzyme())
    {
      ++it;
      ++count;
    }
    TEST_EQUAL(count >= 3, true)
END_SECTION

START_SECTION(ConstEnzymeIterator endEnzyme() const)
    NOT_TESTABLE // tested above
END_SECTION

END_TEST
