// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/DECHARGING/FeatureDecharger.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FeatureDecharger, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureDecharger* ptr = 0;
START_SECTION((FeatureDecharger()))
	ptr = new FeatureDecharger();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((virtual ~FeatureDecharger()))
	delete ptr;
END_SECTION

START_SECTION((FeatureDecharger(const FeatureDecharger &source)))
  FeatureDecharger fdc;

  Param p = fdc.getParameters();
  p.setValue("cluster_rt_mz_relation", 12345.0);
  fdc.setParameters(p);
  
  FeatureDecharger fdc_copy(fdc);
  
  Param p_copy = fdc_copy.getParameters();
  
  TEST_EQUAL((double) p_copy.getValue("cluster_rt_mz_relation"), 12345.0);
  
END_SECTION

START_SECTION((FeatureDecharger& operator=(const FeatureDecharger &source)))
  FeatureDecharger fdc;

  Param p = fdc.getParameters();
  p.setValue("cluster_rt_mz_relation", 12345.0);
  fdc.setParameters(p);
  
  FeatureDecharger fdc_copy = fdc;
  
  Param p_copy = fdc_copy.getParameters();
  
  TEST_EQUAL((double) p_copy.getValue("cluster_rt_mz_relation"), 12345.0);
  
  
END_SECTION

START_SECTION((void compute(FeatureMapType &map)))
	//tested in getFeatureMap()
	NOT_TESTABLE
END_SECTION

START_SECTION((const FeatureMapType& getFeatureMap() const))
  FeatureMap<> map;

  // load a feature map
  FeatureXMLFile().load("data/FeatureDecharger_TestData.featureXML",map);  

  FeatureDecharger fdc;
  fdc.compute(map);
  map = fdc.getFeatureMap();
  
  // combined feature
  TEST_REAL_SIMILAR(map[0].getRT(), 144.576);
  //TEST_REAL_SIMILAR(map[0].getMZ(), 1332.470);
  TEST_REAL_SIMILAR(map[0].getIntensity(), 20000);
  
  // bad feature - but resolved
  TEST_REAL_SIMILAR(map[2].getRT(), 151.897);
  TEST_REAL_SIMILAR(map[2].getMZ(), 789.812);
  TEST_REAL_SIMILAR(map[2].getIntensity(), 55761);

  TEST_REAL_SIMILAR(map[1].getRT(), 151.6);
  TEST_REAL_SIMILAR(map[1].getMZ(), 793.812);
  TEST_REAL_SIMILAR(map[1].getIntensity(), 55761);

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



