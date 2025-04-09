// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/PROCESSING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/KERNEL/Peak2D.h>

///////////////////////////

START_TEST(SavitzkyGolayFilter<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

SavitzkyGolayFilter* dsg_ptr = nullptr;
SavitzkyGolayFilter* dsg_nullPointer = nullptr;
START_SECTION((SavitzkyGolayFilter()))
  dsg_ptr = new SavitzkyGolayFilter;
  TEST_NOT_EQUAL(dsg_ptr, dsg_nullPointer)
END_SECTION

START_SECTION((virtual ~SavitzkyGolayFilter()))
  delete dsg_ptr;
END_SECTION

Param param;
param.setValue("polynomial_order",2);
param.setValue("frame_length",3);

START_SECTION((template < typename PeakType > void filter(MSSpectrum &spectrum)))
  MSSpectrum spectrum;
  spectrum.resize(5);
  MSSpectrum::Iterator it = spectrum.begin();
  for (int i=0; i<5; ++i, ++it)
  {
  	it->setIntensity(0.0f);
    if (i==2)
    {
      it->setIntensity(1.0f);
    }
  }

  SavitzkyGolayFilter sgolay;
	sgolay.setParameters(param);
  sgolay.filter(spectrum);
  
  it=spectrum.begin();
  TEST_REAL_SIMILAR(it->getIntensity(),0.0)
  ++it;
  TEST_REAL_SIMILAR(it->getIntensity(),0.0)
  ++it;
  TEST_REAL_SIMILAR(it->getIntensity(),1.0)
  ++it;
  TEST_REAL_SIMILAR(it->getIntensity(),0.0)
  ++it;
  TEST_REAL_SIMILAR(it->getIntensity(),0.0)
END_SECTION 


START_SECTION((template <typename PeakType> void filterExperiment(MSExperiment<PeakType>& map)))
	TOLERANCE_ABSOLUTE(0.01)

	param.setValue("frame_length",4);

	PeakMap exp;
  exp.resize(4);

  Peak1D p;
  for (int i=0; i<9; ++i)
  {
  	p.setIntensity(0.0f);
    if (i==3)
    {
  		p.setIntensity(1.0f);
    }
    if (i==4)
    {
  		p.setIntensity(0.8f);
    }
    if (i==5)
    {
  		p.setIntensity(1.2f);
    }
    exp[0].push_back(p);
    exp[1].push_back(p);
  }
  exp[2].push_back(p);

  SavitzkyGolayFilter sgolay;
	sgolay.setParameters(param);
  sgolay.filterExperiment(exp);

	TEST_EQUAL(exp.size(),4)
	TEST_EQUAL(exp[0].size(),9)
	TEST_EQUAL(exp[1].size(),9)
	TEST_EQUAL(exp[2].size(),1)
	TEST_EQUAL(exp[3].size(),0)

	TEST_REAL_SIMILAR(exp[0][0].getIntensity(),0.0)
	TEST_REAL_SIMILAR(exp[0][1].getIntensity(),0.0571429)
	TEST_REAL_SIMILAR(exp[0][2].getIntensity(),0.274286)
	TEST_REAL_SIMILAR(exp[0][3].getIntensity(),0.657143)
	TEST_REAL_SIMILAR(exp[0][4].getIntensity(),1.14286)
	TEST_REAL_SIMILAR(exp[0][5].getIntensity(),0.771429)
	TEST_REAL_SIMILAR(exp[0][6].getIntensity(),0.342857)
	TEST_REAL_SIMILAR(exp[0][7].getIntensity(),0.0914286)
	TEST_REAL_SIMILAR(exp[0][8].getIntensity(),0.0)

	TEST_REAL_SIMILAR(exp[1][0].getIntensity(),0.0)
	TEST_REAL_SIMILAR(exp[1][1].getIntensity(),0.0571429)
	TEST_REAL_SIMILAR(exp[1][2].getIntensity(),0.274286)
	TEST_REAL_SIMILAR(exp[1][3].getIntensity(),0.657143)
	TEST_REAL_SIMILAR(exp[1][4].getIntensity(),1.14286)
	TEST_REAL_SIMILAR(exp[1][5].getIntensity(),0.771429)
	TEST_REAL_SIMILAR(exp[1][6].getIntensity(),0.342857)
	TEST_REAL_SIMILAR(exp[1][7].getIntensity(),0.0914286)
	TEST_REAL_SIMILAR(exp[1][8].getIntensity(),0.0)

	TEST_REAL_SIMILAR(exp[2][0].getIntensity(),0.0)

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
