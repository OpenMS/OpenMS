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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/FineIsotopePatternGenerator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FineIsotopePatternGenerator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FineIsotopePatternGenerator* ptr = nullptr;
FineIsotopePatternGenerator* nullPointer = nullptr;
START_SECTION((FineIsotopePatternGenerator()))
  ptr = new FineIsotopePatternGenerator();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~FineIsotopePatternGenerator()))
  delete ptr;
END_SECTION

START_SECTION(( IsotopeDistribution run(const EmpiricalFormula&) const ))
{
  EmpiricalFormula ef ("C6H12O6");
  double threshold = 1e-5;

  {
    FineIsotopePatternGenerator gen;
    IsotopeDistribution id = gen.run(ef);
    TEST_EQUAL(id.size(), 3)

    TEST_REAL_SIMILAR(id[0].getMZ(), 180.063)
    TEST_REAL_SIMILAR(id[0].getIntensity(), 0.922119)

    TEST_REAL_SIMILAR(id[2].getMZ(), 182.068 ) 
    TEST_REAL_SIMILAR(id[2].getIntensity(), 0.0113774 )
  }

  {
    FineIsotopePatternGenerator gen(threshold);
    IsotopeDistribution id = gen.run(ef);
    TEST_EQUAL(id.size(), 14)

    TEST_REAL_SIMILAR(id[0].getMZ(), 180.063)
    TEST_REAL_SIMILAR(id[0].getIntensity(), 0.922119)

    TEST_REAL_SIMILAR(id[4].getMZ(), 182.068 ) 
    TEST_REAL_SIMILAR(id[4].getIntensity(), 0.0113774 )

    TEST_REAL_SIMILAR(id[13].getMZ(), 184.07434277234)
    TEST_REAL_SIMILAR(id[13].getIntensity(), 2.02975552383577e-05)
  }

  {
    FineIsotopePatternGenerator gen(1e-12);
    IsotopeDistribution id = gen.run(ef);
    TEST_EQUAL(id.size(), 104)

    gen.setThreshold(1e-25);
    TEST_EQUAL(gen.run(EmpiricalFormula(ef)).size(), 635)

    gen.setThreshold(1e-50);
    TEST_EQUAL(gen.run(EmpiricalFormula(ef)).size(), 1885)

    gen.setThreshold(1e-100);
    TEST_EQUAL(gen.run(EmpiricalFormula(ef)).size(), 2548)

    gen.setThreshold(1e-1000);
    TEST_EQUAL(gen.run(EmpiricalFormula(ef)).size(), 2548)
  }

  // For a C100 molecule
  {
    FineIsotopePatternGenerator gen;
    gen.setThreshold(1e-2);
    IsotopeDistribution id = gen.run(EmpiricalFormula("C100"));
    TEST_EQUAL(id.size(), 6)

    gen.setThreshold(1e-5);
    TEST_EQUAL(gen.run(EmpiricalFormula("C100")).size(), 9)

    gen.setThreshold(1e-10);
    TEST_EQUAL(gen.run(EmpiricalFormula("C100")).size(), 14)

    gen.setThreshold(1e-20);
    TEST_EQUAL(gen.run(EmpiricalFormula("C100")).size(), 22)

#if 0
    gen.setThreshold(1e-100);
    TEST_EQUAL(gen.run(EmpiricalFormula("C100")).size(), 65)

    gen.setThreshold(1e-150);
    TEST_EQUAL(gen.run(EmpiricalFormula("C100")).size(), 86)

    gen.setThreshold(1e-250);
    TEST_EQUAL(gen.run(EmpiricalFormula("C100")).size(), 101)

    gen.setThreshold(1e-1000);
    TEST_EQUAL(gen.run(EmpiricalFormula("C100")).size(), 101)
#endif
  }
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
