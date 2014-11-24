// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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

#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(WindowMower, "$Id$")

/////////////////////////////////////////////////////////////

WindowMower* e_ptr = 0;
WindowMower* e_nullPointer = 0;
START_SECTION((WindowMower()))
	e_ptr = new WindowMower;
	TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~WindowMower()))
	delete e_ptr;
END_SECTION

e_ptr = new WindowMower();

START_SECTION((WindowMower(const WindowMower& source)))
	WindowMower copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((WindowMower& operator = (const WindowMower& source)))
	WindowMower copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((template<typename SpectrumType> void filterPeakSpectrumForTopNInSlidingWindow(SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);
	TEST_EQUAL(spec.size(), 121)

	Param p(e_ptr->getParameters());
	p.setValue("windowsize", 50.0); // default
	p.setValue("peakcount", 2);  // default
	p.setValue("movetype", "slide"); // default and not needed as we directly call sliding window function
	e_ptr->setParameters(p);
	
	e_ptr->filterPeakSpectrumForTopNInSlidingWindow(spec);
	
	TEST_EQUAL(spec.size(), 56)
	
END_SECTION

START_SECTION((template<typename SpectrumType> void filterPeakSpectrumForTopNInJumpingWindow(SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);
	TEST_EQUAL(spec.size(), 121)

	Param p(e_ptr->getParameters());
	p.setValue("windowsize", 50.0); // default
	p.setValue("peakcount", 2); // default
	p.setValue("movetype", "jump");  // actually not needed as we directly call jumping window function
	e_ptr->setParameters(p);
	e_ptr->filterPeakSpectrumForTopNInJumpingWindow(spec);
	TEST_EQUAL(spec.size(), 30)	
END_SECTION

START_SECTION((void filterPeakMap(PeakMap& exp)))
  DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

  PeakMap pm;
  pm.addSpectrum(spec);

  TEST_EQUAL(pm.begin()->size(), 121)

  Param p(e_ptr->getParameters());
  p.setValue("windowsize", 50.0); // default
  p.setValue("peakcount", 2);
  p.setValue("movetype", "slide"); // default
  e_ptr->setParameters(p);

  e_ptr->filterPeakMap(pm);

  TEST_EQUAL(pm.begin()->size(), 56)
END_SECTION

START_SECTION((void filterPeakSpectrum(PeakSpectrum& spectrum)))
	DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);
  TEST_EQUAL(spec.size(), 121)

  Param p(e_ptr->getParameters());
  p.setValue("windowsize", 50.0); // default
  p.setValue("peakcount", 2);
  p.setValue("movetype", "slide");

  e_ptr->setParameters(p);

  e_ptr->filterPeakSpectrum(spec);

  TEST_EQUAL(spec.size(), 56)
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
