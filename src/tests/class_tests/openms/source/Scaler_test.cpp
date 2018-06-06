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

#include <OpenMS/FILTERING/TRANSFORMERS/Scaler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Scaler, "$Id$")

/////////////////////////////////////////////////////////////

TOLERANCE_ABSOLUTE(0.01)

Scaler* e_ptr = nullptr;
Scaler* e_nullPointer = nullptr;
START_SECTION((Scaler()))
	e_ptr = new Scaler;
	TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~Scaler()))
	delete e_ptr;
END_SECTION

e_ptr = new Scaler();

START_SECTION((Scaler(const Scaler& source)))
	Scaler copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((Scaler& operator = (const Scaler& source)))
	Scaler copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((template<typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

	e_ptr->filterSpectrum(spec);
	
	TEST_EQUAL(spec.size(), 121)

	spec.sortByIntensity();
	TEST_REAL_SIMILAR(spec.begin()->getIntensity(), 96)
	TEST_REAL_SIMILAR((spec.end()-1)->getIntensity(), 121)
	TEST_REAL_SIMILAR((spec.end()-1)->getPosition()[0], 136.077)
	
END_SECTION

START_SECTION((void filterPeakMap(PeakMap& exp)))
	DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

	PeakMap pm;
	pm.addSpectrum(spec);

  e_ptr->filterPeakMap(pm);

  TEST_EQUAL(pm.begin()->size(), 121)

  pm.begin()->sortByIntensity();
  TEST_REAL_SIMILAR(pm.begin()->begin()->getIntensity(), 96)
  TEST_REAL_SIMILAR((pm.begin()->end()-1)->getIntensity(), 121)
  TEST_REAL_SIMILAR((pm.begin()->end()-1)->getPosition()[0], 136.077)

END_SECTION

START_SECTION((void filterPeakSpectrum(PeakSpectrum& spectrum)))
  DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

  e_ptr->filterPeakSpectrum(spec);

  TEST_EQUAL(spec.size(), 121)

  spec.sortByIntensity();
  TEST_REAL_SIMILAR(spec.begin()->getIntensity(), 96)
  TEST_REAL_SIMILAR((spec.end()-1)->getIntensity(), 121)
  TEST_REAL_SIMILAR((spec.end()-1)->getPosition()[0], 136.077)
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
