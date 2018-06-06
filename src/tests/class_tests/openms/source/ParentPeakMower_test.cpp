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

#include <OpenMS/FILTERING/TRANSFORMERS/ParentPeakMower.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(ParentPeakMower, "$Id$")

/////////////////////////////////////////////////////////////

ParentPeakMower* e_ptr = nullptr;
ParentPeakMower* e_nullPointer = nullptr;

START_SECTION((ParentPeakMower()))
	e_ptr = new ParentPeakMower;
  TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~ParentPeakMower()))
	delete e_ptr;
END_SECTION

e_ptr = new ParentPeakMower();

START_SECTION((ParentPeakMower(const ParentPeakMower& source)))
	ParentPeakMower copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((ParentPeakMower& operator = (const ParentPeakMower& source)))
	ParentPeakMower copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((template<typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);
	spec.setMSLevel(2);
	
	spec.sortByPosition();

	TEST_REAL_SIMILAR((spec.begin() + 40)->getIntensity(), 37.5)

	double window_size(2.0);
	Param p(e_ptr->getParameters());
	p.setValue("window_size", window_size);
	p.setValue("default_charge", 2);
	p.setValue("clean_all_charge_states", (short)1);
	p.setValue("set_to_zero", (short)1);
	e_ptr->setParameters(p);

	e_ptr->filterSpectrum(spec);
	double pre_1_pos(spec.getPrecursors()[0].getMZ() * spec.getPrecursors()[0].getCharge());
	for (Int z = 1; z != spec.getPrecursors()[0].getCharge(); ++z)
	{
		for (PeakSpectrum::ConstIterator it = spec.begin(); it != spec.end(); ++it)
		{	
			if (fabs(it->getPosition()[0] - pre_1_pos / double(z)) <= window_size)
			{
				TEST_REAL_SIMILAR(it->getIntensity(), 0.0);
			}

			// test if NH3 loss is correct removed
			if (fabs(it->getPosition()[0] - (pre_1_pos - 17.0) / double(z)) <= window_size)
			{
				TEST_REAL_SIMILAR(it->getIntensity(), 0.0);
			}

			if (fabs(it->getPosition()[0] - (pre_1_pos - 18.0) / double(z)) <= window_size)
			{
				TEST_REAL_SIMILAR(it->getIntensity(), 0.0);
			}
		}
	}
	
END_SECTION

START_SECTION((void filterPeakMap(PeakMap& exp)))
  DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

	PeakMap pm;
	pm.addSpectrum(spec);

  pm.begin()->setMSLevel(2);

  pm.begin()->sortByPosition();

  TEST_REAL_SIMILAR((pm.begin()->begin() + 40)->getIntensity(), 37.5)

  double window_size(2.0);
	Param p(e_ptr->getParameters());
  p.setValue("window_size", window_size);
  p.setValue("default_charge", 2);
  p.setValue("clean_all_charge_states", (short)1);
  p.setValue("set_to_zero", (short)1);
	e_ptr->setParameters(p);

  e_ptr->filterPeakMap(pm);
  double pre_1_pos(pm.begin()->getPrecursors()[0].getMZ() * pm.begin()->getPrecursors()[0].getCharge());
  for (Int z = 1; z != pm.begin()->getPrecursors()[0].getCharge(); ++z)
  {
    for (PeakMap::SpectrumType::ConstIterator it = pm.begin()->begin(); it != pm.begin()->end(); ++it)
    {
      if (fabs(it->getPosition()[0] - pre_1_pos / double(z)) <= window_size)
      {
        TEST_REAL_SIMILAR(it->getIntensity(), 0.0);
      }

      // test if NH3 loss is correct removed
      if (fabs(it->getPosition()[0] - (pre_1_pos - 17.0) / double(z)) <= window_size)
      {
        TEST_REAL_SIMILAR(it->getIntensity(), 0.0);
      }

      if (fabs(it->getPosition()[0] - (pre_1_pos - 18.0) / double(z)) <= window_size)
      {
        TEST_REAL_SIMILAR(it->getIntensity(), 0.0);
      }
    }
  }


END_SECTION

START_SECTION((void filterPeakSpectrum(PeakSpectrum& spectrum)))
  DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);
  spec.setMSLevel(2);

  spec.sortByPosition();

  TEST_REAL_SIMILAR((spec.begin() + 40)->getIntensity(), 37.5)

  double window_size(2.0);
	Param p(e_ptr->getParameters());
  p.setValue("window_size", window_size);
  p.setValue("default_charge", 2);
  p.setValue("clean_all_charge_states", (short)1);
  p.setValue("set_to_zero", (short)1);
	e_ptr->setParameters(p);

  e_ptr->filterPeakSpectrum(spec);
  double pre_1_pos(spec.getPrecursors()[0].getMZ() * spec.getPrecursors()[0].getCharge());
  for (Int z = 1; z != spec.getPrecursors()[0].getCharge(); ++z)
  {
    for (PeakSpectrum::ConstIterator it = spec.begin(); it != spec.end(); ++it)
    {
      if (fabs(it->getPosition()[0] - pre_1_pos / double(z)) <= window_size)
      {
        TEST_REAL_SIMILAR(it->getIntensity(), 0.0);
      }

      // test if NH3 loss is correct removed
      if (fabs(it->getPosition()[0] - (pre_1_pos - 17.0) / double(z)) <= window_size)
      {
        TEST_REAL_SIMILAR(it->getIntensity(), 0.0);
      }

      if (fabs(it->getPosition()[0] - (pre_1_pos - 18.0) / double(z)) <= window_size)
      {
        TEST_REAL_SIMILAR(it->getIntensity(), 0.0);
      }
    }
  }

	
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
