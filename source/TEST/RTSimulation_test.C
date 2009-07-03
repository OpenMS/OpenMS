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
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/SIMULATION/RTSimulation.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(RTSimulation, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

const unsigned long rnd_gen_seed = 1;
RTSimulation* ptr = 0;
START_SECTION((RTSimulation(const gsl_rng *random_generator)))
{
	ptr = new RTSimulation(NULL);
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~RTSimulation())
{
	delete ptr;
}
END_SECTION

START_SECTION((RTSimulation(const RTSimulation &source)))
{
  RTSimulation source(NULL);
  Param p = source.getParameters();
  p.setValue("total_gradient_time",4000.0);
  source.setParameters(p);
  
  RTSimulation target(source);
  TEST_EQUAL(source.getParameters(), target.getParameters())
  TEST_EQUAL(source.getGradientTime(), target.getGradientTime())
}
END_SECTION

START_SECTION((virtual ~RTSimulation()))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((RTSimulation& operator=(const RTSimulation &source)))
{
  RTSimulation source(NULL);
  RTSimulation target(source);
  
  Param p = source.getParameters();
  p.setValue("total_gradient_time",4000.0);
  source.setParameters(p);
  
  TEST_NOT_EQUAL(source.getParameters(), target.getParameters())
  target = source;
  TEST_EQUAL(source.getParameters(), target.getParameters())
}
END_SECTION

START_SECTION(( void predictRT(FeatureMapSim & features, MSSimExperiment & experiment) ))
{
  // init rng 
  gsl_rng* rnd_gen = gsl_rng_alloc (gsl_rng_taus);
  gsl_rng_set(rnd_gen, rnd_gen_seed);
  
  // no rt scan
  RTSimulation no_rt_sim(rnd_gen);
  Param p = no_rt_sim.getParameters();
  p.setValue("rt_column_on","false");
  p.setValue("total_gradient_time",4000.0);
  no_rt_sim.setParameters(p);  

  FeatureMapSim no_rt_features;
  StringList peps = StringList::create("TVQMENQFVAFVDK,ACHKKKKHHACAC,AAAAHTKLRTTIPPEFG,RYCNHKTUIKL");
	for (StringList::const_iterator it=peps.begin(); it!=peps.end(); ++it)
	{
		Feature f;
		PeptideIdentification pep_id;
		pep_id.insertHit(PeptideHit(1.0, 1, 1, *it));
		f.getPeptideIdentifications().push_back(pep_id);
		f.setIntensity(10);
		no_rt_features.push_back(f);
	}
  
	MSSimExperiment experiment_no_rt;
	no_rt_sim.predictRT(no_rt_features, experiment_no_rt);
  TEST_EQUAL(experiment_no_rt.size(), 1);
  for(FeatureMapSim::const_iterator fIt = no_rt_features.begin(); fIt != no_rt_features.end();
      ++fIt)
  {
    TEST_EQUAL((*fIt).getRT(), -1);
  }
  
  // rt svm
  // no rt scan
  RTSimulation svm_rt_sim(rnd_gen);
  Param svm_params = svm_rt_sim.getParameters();
  svm_params.setValue("rt_column_on","true");
  svm_params.setValue("total_gradient_time",4000.0);
  svm_params.setValue("rt_model_file",OPENMS_GET_TEST_DATA_PATH("RTSimulation.svm"));
  svm_params.setValue("rt_shift_mean", 0);
  svm_params.setValue("rt_shift_stddev", 50);
  
  svm_rt_sim.setParameters(svm_params);
  
  FeatureMapSim svm_rt_features;

	for (StringList::const_iterator it=peps.begin(); it!=peps.end(); ++it)
	{
		Feature f;
		PeptideIdentification pep_id;
		pep_id.insertHit(PeptideHit(1.0, 1, 1, *it));
		f.getPeptideIdentifications().push_back(pep_id);
		f.setIntensity(10);
		svm_rt_features.push_back(f);
	}  

  MSSimExperiment experiment_rt;  
  svm_rt_sim.predictRT(svm_rt_features, experiment_rt);
  
  TEST_EQUAL(svm_rt_features.size(), 4)
  
  // TODO: check why these are different now & test MSExperiment generation
  //TEST_REAL_SIMILAR(svm_rt_features[0].getRT(), 1597.44)
  TEST_EQUAL(svm_rt_features[0].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "TVQMENQFVAFVDK")

  //TEST_REAL_SIMILAR(svm_rt_features[1].getRT(), 406.6)
  TEST_EQUAL(svm_rt_features[1].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "ACHKKKKHHACAC")
  
  //TEST_REAL_SIMILAR(svm_rt_features[2].getRT(), 1151.26)
  TEST_EQUAL(svm_rt_features[2].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "AAAAHTKLRTTIPPEFG")
  
  //TEST_REAL_SIMILAR(svm_rt_features[3].getRT(), 832.56)
  TEST_EQUAL(svm_rt_features[3].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "RYCNHKTUIKL")

  /*
  for(FeatureMapSim::const_iterator fIt = svm_rt_features.begin(); fIt != svm_rt_features.end();
      ++fIt)
  {
    std::cout << (*fIt).getRT() << " " << (*fIt).getPeptideIdentifications()[0].getHits()[0].getSequence().toString() << std::endl;
  }
   */
}
END_SECTION

START_SECTION((void predictContaminantsRT(FeatureMapSim &)))
{
  // TODO
}
END_SECTION

START_SECTION((bool isRTColumnOn() const ))
{
  RTSimulation rt_sim(NULL);

  Param p = rt_sim.getParameters();
  p.setValue("rt_column_on","true");
  rt_sim.setParameters(p);  
  
  TEST_EQUAL(rt_sim.isRTColumnOn(), true);
  
  p.setValue("rt_column_on","false");
  rt_sim.setParameters(p);  
  
  TEST_EQUAL(rt_sim.isRTColumnOn(), false);
}
END_SECTION

START_SECTION((SimCoordinateType getGradientTime() const ))
{
  RTSimulation rt_sim(NULL);
  
  Param p = rt_sim.getParameters();
  p.setValue("total_gradient_time",1000.0);
  rt_sim.setParameters(p);  
  
  TEST_EQUAL(rt_sim.getGradientTime(), 1000.0);
  
  p.setValue("total_gradient_time",4000.0);
  rt_sim.setParameters(p);  
  
  TEST_EQUAL(rt_sim.getGradientTime(), 4000.0);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



