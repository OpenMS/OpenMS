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
// $Maintainer: Mathias Walzer $
// $Authors: Volker Mosthaf, Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(ThresholdMower, "$Id$")

/////////////////////////////////////////////////////////////

ThresholdMower* e_ptr = nullptr;
ThresholdMower* e_nullPointer = nullptr;
START_SECTION((ThresholdMower()))
  e_ptr = new ThresholdMower;
  TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~ThresholdMower()))
  delete e_ptr;
END_SECTION

e_ptr = new ThresholdMower();

START_SECTION((ThresholdMower(const ThresholdMower& source)))
  ThresholdMower copy(*e_ptr);
  TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
  TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((ThresholdMower& operator=(const ThresholdMower& source)))
  ThresholdMower copy;
  copy = *e_ptr;
  TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
  TEST_EQUAL(copy.getName(), e_ptr->getName());
END_SECTION

START_SECTION((template<typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)))
  DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);
  
  TEST_EQUAL(spec.size(), 121)

  Param p(e_ptr->getParameters());
  p.setValue("threshold", 1.0);
  e_ptr->setParameters(p);

  e_ptr->filterSpectrum(spec);
  TEST_EQUAL(spec.size(), 121)

  p.setValue("threshold", 10.0);
  e_ptr->setParameters(p);

  e_ptr->filterSpectrum(spec);
  TEST_EQUAL(spec.size(), 14)
END_SECTION

START_SECTION((void filterPeakMap(PeakMap& exp)))
  DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

  PeakMap pm;
  pm.addSpectrum(spec);

  TEST_EQUAL(pm.begin()->size(), 121)

  Param p(e_ptr->getParameters());
  p.setValue("threshold", 1.0);
  e_ptr->setParameters(p);

  e_ptr->filterPeakMap(pm);
  TEST_EQUAL(pm.begin()->size(), 121)

  p.setValue("threshold", 10.0);
  e_ptr->setParameters(p);
  e_ptr->filterPeakMap(pm);
  TEST_EQUAL(pm.begin()->size(), 14)

END_SECTION

START_SECTION((void filterPeakSpectrum(PeakSpectrum& spectrum)))
  DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

  TEST_EQUAL(spec.size(), 121)

  Param p(e_ptr->getParameters());
  p.setValue("threshold", 1.0);
  e_ptr->setParameters(p);

  e_ptr->filterPeakSpectrum(spec);
  TEST_EQUAL(spec.size(), 121)

  p.setValue("threshold", 10.0);
  e_ptr->setParameters(p);
  e_ptr->filterPeakSpectrum(spec);
  TEST_EQUAL(spec.size(), 14)
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
