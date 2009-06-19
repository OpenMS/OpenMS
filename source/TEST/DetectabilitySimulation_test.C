// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/SIMULATION/DetectabilitySimulation.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DetectabilitySimulation, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DetectabilitySimulation* ptr = 0;
START_SECTION(DetectabilitySimulation())
{
	ptr = new DetectabilitySimulation();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(DetectabilitySimulation())
{
	delete ptr;
}
END_SECTION

START_SECTION((DetectabilitySimulation(const DetectabilitySimulation &source)))
{
  DetectabilitySimulation source;
  Param p = source.getParameters();
  p.setValue("min_detect",0.0);
  source.setParameters(p);
  
  DetectabilitySimulation target(source);
  TEST_EQUAL(source.getParameters(), target.getParameters())
}
END_SECTION

START_SECTION((virtual ~IonizationSimulation()))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((DetectabilitySimulation& operator=(const DetectabilitySimulation &source)))
{
  DetectabilitySimulation detect_sim1;
  DetectabilitySimulation detect_sim2(detect_sim1);
  
	Param p = detect_sim1.getParameters();
	p.setValue("min_detect", 0.0);
	detect_sim1.setParameters(p);
	TEST_NOT_EQUAL(detect_sim1.getParameters(),detect_sim2.getParameters());
	detect_sim2 = detect_sim1;
	TEST_EQUAL(detect_sim2.getParameters(),detect_sim2.getParameters());
}
END_SECTION

START_SECTION((void filterDetectability(FeatureMapSim &)))
{
  // test no detect
  DetectabilitySimulation detect_off;
  Param p = detect_off.getParameters();
  p.setValue("dt_simulation_on","false");
  p.setValue("min_detect", 0.9);
  detect_off.setParameters(p);
  
  FeatureMapSim no_detect_features;
  StringList peps = StringList::create("TVQMENQFVAFVDK,ACHKKKKHHACAC,AAAAHTKLRTTIPPEFG,RYCNHKTUIKL");
	for (StringList::const_iterator it=peps.begin(); it!=peps.end(); ++it)
	{
		Feature f;
		PeptideIdentification pep_id;
		pep_id.insertHit(PeptideHit(1.0, 1, 1, *it));
		f.getPeptideIdentifications().push_back(pep_id);
		f.setIntensity(10);
		no_detect_features.push_back(f);
	}
  
  detect_off.filterDetectability(no_detect_features);
  
  TEST_EQUAL(no_detect_features.size(), 4)
  for(Size i = 0 ; i < no_detect_features.size() ; ++i) 
  {
    TEST_EQUAL(no_detect_features[i].getMetaValue("detectibility"), 1.0)
  }
  
  // test svm
  DetectabilitySimulation detect_svm;
  Param svm_params = detect_svm.getParameters();
  svm_params.setValue("dt_simulation_on","true");
  svm_params.setValue("min_detect", 0.4);
  svm_params.setValue("dt_model_file",OPENMS_GET_TEST_DATA_PATH("DetectabilitySimulation.svm"));
  detect_svm.setParameters(svm_params);
  
  FeatureMapSim svm_features;
	for (StringList::const_iterator it=peps.begin(); it!=peps.end(); ++it)
	{
		Feature f;
		PeptideIdentification pep_id;
		pep_id.insertHit(PeptideHit(1.0, 1, 1, *it));
		f.getPeptideIdentifications().push_back(pep_id);
		f.setIntensity(10);
		svm_features.push_back(f);
	}

  detect_svm.filterDetectability(svm_features);
  
  TEST_EQUAL(svm_features.size(), 2)
  TEST_EQUAL(svm_features[0].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "TVQMENQFVAFVDK")
  TEST_REAL_SIMILAR(svm_features[0].getMetaValue("detectability"), 0.869237485950867)
  TEST_EQUAL(svm_features[1].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "AAAAHTKLRTTIPPEFG")             
  TEST_REAL_SIMILAR(svm_features[1].getMetaValue("detectability"), 0.723545391996237)
  
  /*
  for(FeatureMapSim::const_iterator it = svm_features.begin() ; it != svm_features.end();
      ++it)
  {
    std::cout << (*it).getPeptideIdentifications()[0].getHits()[0].getSequence().toString()  << " " << (*it).getMetaValue("detectibility") << std::endl;
  }
  */
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



