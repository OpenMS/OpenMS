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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Normalizer, "$Id$")

/////////////////////////////////////////////////////////////

Normalizer* e_ptr = nullptr;
Normalizer* e_nullPointer = nullptr;

START_SECTION((Normalizer()))
	e_ptr = new Normalizer;
  TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~Normalizer()))
	delete e_ptr;
END_SECTION

e_ptr = new Normalizer();

START_SECTION((Normalizer(const Normalizer& source)))
	Normalizer copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((Normalizer& operator = (const Normalizer& source)))
	Normalizer copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION


DTAFile dta_file;
PeakSpectrum spec_ref;
dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec_ref);
spec_ref.sortByIntensity();

START_SECTION((template<typename SpectrumType> void filterSpectrum(SpectrumType& spectrum) const))
	PeakSpectrum spec = spec_ref;
	TEST_EQUAL(spec.rbegin()->getIntensity(), 46)
	e_ptr->filterSpectrum(spec);
	TEST_EQUAL(spec.rbegin()->getIntensity(), 1)

	Param p(e_ptr->getParameters());
	p.setValue("method", "to_TIC");
	e_ptr->setParameters(p);
	e_ptr->filterSpectrum(spec);

	double sum(0);
	for (PeakSpectrum::ConstIterator it = spec.begin(); it != spec.end(); ++it)
	{
		sum += it->getIntensity();
	}

	TEST_REAL_SIMILAR(sum, 1.0);	
END_SECTION
	
START_SECTION((void filterPeakMap(PeakMap& exp) const))
	delete e_ptr;
	e_ptr = new Normalizer();

  PeakSpectrum spec = spec_ref;

	PeakMap pm;
	pm.addSpectrum(spec);

  TEST_EQUAL(pm.begin()->rbegin()->getIntensity(), 46)

  e_ptr->filterPeakMap(pm);

  TEST_EQUAL(pm.begin()->rbegin()->getIntensity(), 1)

	Param p(e_ptr->getParameters());
	p.setValue("method", "to_TIC");
	e_ptr->setParameters(p);
  e_ptr->filterPeakMap(pm);

  double sum(0);
  for (PeakMap::SpectrumType::ConstIterator it = pm.begin()->begin(); it != pm.begin()->end(); ++it)
  {
    sum += it->getIntensity();
  }
  TEST_REAL_SIMILAR(sum, 1.0);	
END_SECTION

START_SECTION((void filterPeakSpectrum(PeakSpectrum& spectrum) const))
	delete e_ptr;
	e_ptr = new Normalizer();
  PeakSpectrum spec = spec_ref;
  e_ptr->filterPeakSpectrum(spec);
  TEST_EQUAL(spec.rbegin()->getIntensity(), 1)

	Param p(e_ptr->getParameters());
	p.setValue("method", "to_TIC");
	e_ptr->setParameters(p);
  e_ptr->filterPeakSpectrum(spec);

  double sum(0);
  for (PeakSpectrum::ConstIterator it = spec.begin(); it != spec.end(); ++it)
  {
    sum += it->getIntensity();
  }
  TEST_REAL_SIMILAR(sum, 1.0);	
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
