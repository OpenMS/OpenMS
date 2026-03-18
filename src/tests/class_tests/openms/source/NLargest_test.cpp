// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/PROCESSING/FILTERING/NLargest.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(NLargest, "$Id$")

/////////////////////////////////////////////////////////////

NLargest* e_ptr = nullptr;
NLargest* e_nullPointer = nullptr;

START_SECTION((NLargest()))
	e_ptr = new NLargest;
  TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION(NLargest(UInt n))
	NLargest filter(10);
  TEST_EQUAL((UInt)filter.getParameters().getValue("n"), 10)
END_SECTION

START_SECTION((~NLargest()))
	delete e_ptr;
END_SECTION

e_ptr = new NLargest();

START_SECTION((NLargest(const NLargest& source)))
	NLargest copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((NLargest& operator=(const NLargest& source)))
	NLargest copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((template<typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);
	TEST_EQUAL(spec.size(), 121)

	Param p(e_ptr->getParameters());
	p.setValue("n", 10);
	e_ptr->setParameters(p);
	e_ptr->filterSpectrum(spec);
	TEST_EQUAL(spec.size(), 10)
END_SECTION

START_SECTION((void filterPeakMap(PeakMap& exp)))
	delete e_ptr;
	e_ptr = new NLargest();
  DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

  PeakMap pm;
  pm.addSpectrum(spec);

  TEST_EQUAL(pm.begin()->size(), 121)

  Param p(e_ptr->getParameters());
  p.setValue("n", 10);
  e_ptr->setParameters(p);
  e_ptr->filterPeakMap(pm);
  TEST_EQUAL(pm.begin()->size(), 10)
END_SECTION

START_SECTION((void filterPeakSpectrum(PeakSpectrum& spectrum)))
  delete e_ptr;
  e_ptr = new NLargest();
  DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);
  TEST_EQUAL(spec.size(), 121)
	
  Param p(e_ptr->getParameters());
  p.setValue("n", 10);
  e_ptr->setParameters(p);
  e_ptr->filterPeakSpectrum(spec);
  TEST_EQUAL(spec.size(), 10)

  PeakSpectrum s_da;
  s_da.getIntegerDataArrays().resize(1); 
  s_da.getStringDataArrays().resize(1);
  // create a "triangle" shape with apex at i=50 
  for (Size i = 0; i != 50; ++i)
  {
    s_da.push_back(Peak1D(i, i + 0.1)); 
    s_da.getIntegerDataArrays()[0].push_back(i); 
    s_da.getStringDataArrays()[0].push_back("up"); 
  }
  for (int i = 50; i != 100; ++i)
  {
    s_da.push_back(Peak1D(i, (100 - i) + 0.2)); 
    s_da.getIntegerDataArrays()[0].push_back(i); 
    s_da.getStringDataArrays()[0].push_back("down"); 
  }
  e_ptr->filterPeakSpectrum(s_da);

/*
int  mz DA_int DA_string
50.2 50  50      down
49.2 51  51      down
49.1 49  49      up
48.2 52  52      down
48.1 48  48      up
47.2 53  53      down
47.1 47  47      up
46.2 54  54      down
46.1 46  46      up
45.2 55  55      down
*/

  TEST_EQUAL(s_da.size(), 10)
  TEST_EQUAL(s_da[0].getIntensity(), 50.2)
  TEST_EQUAL(s_da[1].getIntensity(), 49.2)
  TEST_EQUAL(s_da[2].getIntensity(), 49.1)
  TEST_EQUAL(s_da.getIntegerDataArrays()[0][0], 50)
  TEST_EQUAL(s_da.getIntegerDataArrays()[0][1], 51)
  TEST_EQUAL(s_da.getIntegerDataArrays()[0][2], 49)
  TEST_EQUAL(s_da.getStringDataArrays()[0][0], "down")
  TEST_EQUAL(s_da.getStringDataArrays()[0][1], "down")
  TEST_EQUAL(s_da.getStringDataArrays()[0][2], "up")
  TEST_EQUAL(s_da[7].getIntensity(), 46.2)
  TEST_EQUAL(s_da[8].getIntensity(), 46.1)
  TEST_EQUAL(s_da[9].getIntensity(), 45.2)
  TEST_EQUAL(s_da.getIntegerDataArrays()[0][7], 54)
  TEST_EQUAL(s_da.getIntegerDataArrays()[0][8], 46)
  TEST_EQUAL(s_da.getIntegerDataArrays()[0][9], 55)
  TEST_EQUAL(s_da.getStringDataArrays()[0][7], "down")
  TEST_EQUAL(s_da.getStringDataArrays()[0][8], "up")
  TEST_EQUAL(s_da.getStringDataArrays()[0][9], "down")
  
  // debug code
  // for (Size i = 0; i != s_da.size(); ++i) 
  // {
  //  cout << "int:" << s_da[i].getIntensity() << " mz:" << s_da[i].getMZ()  << "\t" << s_da.getIntegerDataArrays()[0][i] << "\t" << s_da.getStringDataArrays()[0][i] << endl;
  // }
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
