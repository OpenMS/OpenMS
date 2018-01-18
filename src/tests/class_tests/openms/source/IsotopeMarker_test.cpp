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
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FILTERING/TRANSFORMERS/IsotopeMarker.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <map>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(IsotopeMarker, "$Id$")

/////////////////////////////////////////////////////////////

IsotopeMarker* e_ptr = nullptr;
IsotopeMarker* e_nullPointer = nullptr;

START_SECTION((IsotopeMarker()))
	e_ptr = new IsotopeMarker;
  TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~IsotopeMarker()))
	delete e_ptr;
END_SECTION

e_ptr = new IsotopeMarker();

START_SECTION((IsotopeMarker(const IsotopeMarker& source)))
	IsotopeMarker copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((IsotopeMarker& operator=(const IsotopeMarker& source)))
	IsotopeMarker copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((template<typename SpectrumType> void apply(std::map<double, bool>& marked, SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

	map<double, bool> marked;
	e_ptr->apply(marked, spec);
	
	TEST_EQUAL(marked.size(), 48)

	Param iso_param = e_ptr->getParameters();
	iso_param.setValue("marks", 2);
	e_ptr->setParameters(iso_param);
	marked.clear();
	e_ptr->apply(marked, spec);
	TEST_EQUAL(marked.size(), 17)
END_SECTION

START_SECTION((static PeakMarker* create()))
	PeakMarker* pm = IsotopeMarker::create();
	IsotopeMarker im;
	TEST_EQUAL(pm->getParameters(), im.getParameters())
	TEST_EQUAL(pm->getName(), im.getName())
END_SECTION

START_SECTION((static const String getProductName()))
	TEST_EQUAL(e_ptr->getProductName(), "IsotopeMarker")
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
