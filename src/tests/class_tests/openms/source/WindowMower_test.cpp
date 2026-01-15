// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Volker Mosthaf, Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/PROCESSING/FILTERING/WindowMower.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(WindowMower, "$Id$")

/////////////////////////////////////////////////////////////

WindowMower* e_ptr = nullptr;
WindowMower* e_nullPointer = nullptr;
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


  // test data array handling
  PeakSpectrum s_da;
  // create a "triangle" shape with apex at i=50 
/*
 int  mz   DA_int  DA_string
  0.1 0    0       up
  1.1 1    1       up
  2.1 2    2       up
  ...
  47.1 47  47      up
  48.1 48  48      up
  49.1 49  49      up
  50.2 50  50      down
  49.2 51  51      down
  48.2 52  52      down
  ...
  3.2 97   97      down
  2.2 98   98      down
  1.2 99   99      down
*/
  p.setValue("movetype", "slide");
  e_ptr->setParameters(p);
  s_da.getIntegerDataArrays().resize(1); 
  s_da.getStringDataArrays().resize(1);
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

/* result: the 4 rows in the middle: (48,49) + (49,50) + (50, 51) = 48,49,50,51
int  mz DA_int DA_string
48.1 48  48      up        
49.1 49  49      up
50.2 50  50      down
49.2 51  51      down
*/

  TEST_EQUAL(s_da.size(), 4)
  TEST_EQUAL(s_da[0].getIntensity(), 48.1)
  TEST_EQUAL(s_da[1].getIntensity(), 49.1)
  TEST_EQUAL(s_da[2].getIntensity(), 50.2)
  TEST_EQUAL(s_da[3].getIntensity(), 49.2)
  TEST_EQUAL(s_da.getIntegerDataArrays()[0][0], 48)
  TEST_EQUAL(s_da.getIntegerDataArrays()[0][1], 49)
  TEST_EQUAL(s_da.getIntegerDataArrays()[0][2], 50)
  TEST_EQUAL(s_da.getIntegerDataArrays()[0][3], 51)
  TEST_EQUAL(s_da.getStringDataArrays()[0][0], "up")
  TEST_EQUAL(s_da.getStringDataArrays()[0][1], "up")
  TEST_EQUAL(s_da.getStringDataArrays()[0][2], "down")
  TEST_EQUAL(s_da.getStringDataArrays()[0][3], "down")

  p.setValue("movetype", "jump");
  e_ptr->setParameters(p);
  s_da.clear(true);
  s_da.getIntegerDataArrays().resize(1); 
  s_da.getStringDataArrays().resize(1);

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

/* result: first window from m/z 0 to 49 and second window from m/z 50 to 99
int  mz  DA_int DA_string
48.1  48     48      up
49.1  49     49      up
50.2  50     50      down
49.2  51     51      down
*/
  TEST_EQUAL(s_da.size(), 4)
  TEST_EQUAL(s_da[0].getIntensity(), 48.1)
  TEST_EQUAL(s_da[1].getIntensity(), 49.1)
  TEST_EQUAL(s_da[2].getIntensity(), 50.2)
  TEST_EQUAL(s_da[3].getIntensity(), 49.2)
  TEST_EQUAL(s_da.getIntegerDataArrays()[0][0], 48)
  TEST_EQUAL(s_da.getIntegerDataArrays()[0][1], 49)
  TEST_EQUAL(s_da.getIntegerDataArrays()[0][2], 50)
  TEST_EQUAL(s_da.getIntegerDataArrays()[0][3], 51)
  TEST_EQUAL(s_da.getStringDataArrays()[0][0], "up")
  TEST_EQUAL(s_da.getStringDataArrays()[0][1], "up")
  TEST_EQUAL(s_da.getStringDataArrays()[0][2], "down")
  TEST_EQUAL(s_da.getStringDataArrays()[0][3], "down")

  p.setValue("windowsize", 10.0); 
  e_ptr->setParameters(p);
  s_da.clear(true);
  s_da.getIntegerDataArrays().resize(1); 
  s_da.getStringDataArrays().resize(1);

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
int  mz  DA_int DA_string
8.1 8 8 up
9.1 9 9 up
18.1 18 18 up
19.1 19 19 up
28.1 28 28 up
29.1 29 29 up
38.1 38 38 up
39.1 39 39 up
48.1 48 48 up
49.1 49 49 up
50.2 50 50 down
49.2 51 51 down
40.2 60 60 down
39.2 61 61 down
30.2 70 70 down
29.2 71 71 down
20.2 80 80 down
19.2 81 81 down
10.2 90 90 down
9.2 91 91 down
*/
// note that the last window contains only one peak 
// because the peak fraction in window mower is 0.9
  TEST_EQUAL(s_da.size(), 20)


END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
