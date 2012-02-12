// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/SIMULATION/IonizationSimulation.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(IonizationSimulation, "$Id$")

// to avoid parallel random number issues
TOPPBase::setMaxNumberOfThreads(1);

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IonizationSimulation* ptr = 0;
IonizationSimulation* nullPointer = 0;
const unsigned long rnd_gen_seed = 1;
SimRandomNumberGenerator rnd_gen;

// init reproducible rnd_gen
rnd_gen.technical_rng = gsl_rng_alloc(gsl_rng_mt19937);
gsl_rng_set(rnd_gen.technical_rng, 0);
rnd_gen.biological_rng = gsl_rng_alloc(gsl_rng_mt19937);
gsl_rng_set(rnd_gen.biological_rng, 0);

START_SECTION(IonizationSimulation())
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((IonizationSimulation(const SimRandomNumberGenerator& )))
{
  ptr = new IonizationSimulation(rnd_gen);
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~IonizationSimulation())
{
	delete ptr;
}
END_SECTION

START_SECTION((IonizationSimulation(const IonizationSimulation &source)))
{
  IonizationSimulation source(rnd_gen);
  Param p = source.getParameters();
  p.setValue("ionization_type","MALDI");
  source.setParameters(p);
  
  IonizationSimulation target(source);
  TEST_EQUAL(source.getParameters(), target.getParameters())
  
}
END_SECTION

START_SECTION((IonizationSimulation& operator=(const IonizationSimulation &source)))
{
  IonizationSimulation ion_sim1(rnd_gen);
  IonizationSimulation ion_sim2(ion_sim1);
  
	Param p = ion_sim1.getParameters();
	p.setValue("ionization_type", "MALDI");
	ion_sim1.setParameters(p);
	TEST_NOT_EQUAL(ion_sim1.getParameters(),ion_sim2.getParameters());
	ion_sim2 = ion_sim1;
	TEST_EQUAL(ion_sim2.getParameters(),ion_sim2.getParameters());
}
END_SECTION

START_SECTION((void ionize(FeatureMapSim &features, ConsensusMap &charge_consensus, MSSimExperiment &experiment)))
{
  // init rng 
  SimRandomNumberGenerator rnd_gen;

  rnd_gen.biological_rng = gsl_rng_alloc (gsl_rng_taus);
  gsl_rng_set(rnd_gen.biological_rng, rnd_gen_seed);
  rnd_gen.technical_rng = gsl_rng_alloc (gsl_rng_taus);
  gsl_rng_set(rnd_gen.technical_rng, rnd_gen_seed);
  
  // testing ESI
  IonizationSimulation esi_sim(rnd_gen);
  Param esi_param = esi_sim.getParameters();
  esi_param.setValue("ionization_type","ESI");
  esi_param.setValue("esi:ionized_residues",StringList::create("Arg,Lys,His"));
  esi_param.setValue("esi:ionization_probability", 0.8);
  esi_param.setValue("esi:charge_impurity", StringList::create("H+:1,NH4+:0.2,Ca++:0.1"));
  esi_param.setValue("esi:max_impurity_set_size", 3);
  
  esi_sim.setParameters(esi_param);
  
  FeatureMapSim esi_features;
	ConsensusMap cm;
  StringList peps = StringList::create("TVQMENQFVAFVDK,ACHKKKKHHACAC,AAAAHTKLRTTIPPEFG,RYCNHKTUIKL");
	for (StringList::const_iterator it=peps.begin(); it!=peps.end(); ++it)
	{
		Feature f;
		PeptideIdentification pep_id;
		pep_id.insertHit(PeptideHit(1.0, 1, 1, *it));
		f.getPeptideIdentifications().push_back(pep_id);
		f.setIntensity(10);
		esi_features.push_back(f);
	}
  
	MSSimExperiment exp;
	MSSimExperiment::SpectrumType spec;
	exp.push_back(spec);

  esi_sim.ionize(esi_features, cm, exp);
    
  TEST_EQUAL(esi_features.size(), 22)
  ABORT_IF(esi_features.size()!=22)

  TEST_EQUAL(esi_features[0].getCharge(), 2)
  TEST_EQUAL(esi_features[0].getIntensity(), 6)

  TEST_EQUAL(esi_features[1].getCharge(), 2)
  TEST_EQUAL(esi_features[1].getIntensity(), 2)

  TEST_EQUAL(esi_features[2].getCharge(), 1)
  TEST_EQUAL(esi_features[2].getIntensity(), 2)

  TEST_EQUAL(esi_features[3].getCharge(), 5)
  TEST_EQUAL(esi_features[3].getIntensity(), 3)
  
  TEST_EQUAL(esi_features[4].getCharge(), 7)
  TEST_EQUAL(esi_features[4].getIntensity(), 1)

  TEST_EQUAL(esi_features[5].getCharge(), 7)
  TEST_EQUAL(esi_features[5].getIntensity(), 1)

  TEST_EQUAL(esi_features[6].getCharge(), 6)
  TEST_EQUAL(esi_features[6].getIntensity(), 1)

  TEST_EQUAL(esi_features[7].getCharge(), 6)
  TEST_EQUAL(esi_features[7].getIntensity(), 1)

  TEST_EQUAL(esi_features[8].getCharge(), 4)
  TEST_EQUAL(esi_features[8].getIntensity(), 2)

  TEST_EQUAL(esi_features[9].getCharge(), 3)
  TEST_EQUAL(esi_features[9].getIntensity(), 2)

  TEST_EQUAL(esi_features[10].getCharge(), 5)
  TEST_EQUAL(esi_features[10].getIntensity(), 1)

  TEST_EQUAL(esi_features[11].getCharge(), 5)
  TEST_EQUAL(esi_features[11].getIntensity(), 1)
  
  TEST_EQUAL(esi_features[12].getCharge(), 4)
  TEST_EQUAL(esi_features[12].getIntensity(), 1)

  TEST_EQUAL(esi_features[13].getCharge(), 3)
  TEST_EQUAL(esi_features[13].getIntensity(), 1)

  TEST_EQUAL(esi_features[14].getCharge(), 3)
  TEST_EQUAL(esi_features[14].getIntensity(), 1)

  TEST_EQUAL(esi_features[15].getCharge(), 2)
  TEST_EQUAL(esi_features[15].getIntensity(), 1)

  TEST_EQUAL(esi_features[16].getCharge(), 4)
  TEST_EQUAL(esi_features[16].getIntensity(), 3)

	TEST_EQUAL(esi_features[17].getCharge(), 6)
  TEST_EQUAL(esi_features[17].getIntensity(), 2)

  TEST_EQUAL(esi_features[18].getCharge(), 5)
  TEST_EQUAL(esi_features[18].getIntensity(), 2)

  TEST_EQUAL(esi_features[19].getCharge(), 6)
  TEST_EQUAL(esi_features[19].getIntensity(), 1)

  TEST_EQUAL(esi_features[20].getCharge(), 4)
  TEST_EQUAL(esi_features[20].getIntensity(), 1)

  TEST_EQUAL(esi_features[21].getCharge(), 4)
  TEST_EQUAL(esi_features[21].getIntensity(), 1)

  for(FeatureMapSim::const_iterator fmIt = esi_features.begin(); fmIt != esi_features.end();
      ++fmIt)
  {
    std::cout << (*fmIt).getCharge() << " " 
							<< (*fmIt).getIntensity() << " " 
							<< (*fmIt).getPeptideIdentifications()[0].getHits()[0].getSequence().toString() 
							<< " Adducts: " << (*fmIt).getMetaValue("charge_adducts")
							<< " Parent: " << (*fmIt).getMetaValue("parent_feature_number")
							<< std::endl;
  }
  
  
  SimRandomNumberGenerator rnd_gen_maldi;

  rnd_gen_maldi.biological_rng = gsl_rng_alloc (gsl_rng_taus);
  gsl_rng_set(rnd_gen_maldi.biological_rng, rnd_gen_seed);
  rnd_gen_maldi.technical_rng = gsl_rng_alloc (gsl_rng_taus);
  gsl_rng_set(rnd_gen_maldi.technical_rng, rnd_gen_seed);

  // testing MALDI
  IonizationSimulation maldi_sim(rnd_gen_maldi);
  Param maldi_param = maldi_sim.getParameters();
  maldi_param.setValue("ionization_type","MALDI");
  maldi_param.setValue("maldi:ionization_probabilities", DoubleList::create("0.9,0.1"));
  
  maldi_sim.setParameters(maldi_param);
  
  FeatureMapSim maldi_features;
	for (StringList::const_iterator it=peps.begin(); it!=peps.end(); ++it)
	{
		Feature f;
		PeptideIdentification pep_id;
		pep_id.insertHit(PeptideHit(1.0, 1, 1, *it));
		f.getPeptideIdentifications().push_back(pep_id);
		f.setIntensity(10);
    maldi_features.push_back(f);
	}
  
	MSSimExperiment expt;
	MSSimExperiment::SpectrumType spect;
	expt.push_back(spect);
	maldi_sim.ionize(maldi_features, cm, expt);

  TEST_EQUAL(maldi_features.size(), 7)

  TEST_EQUAL(maldi_features[0].getCharge(), 1)
  TEST_EQUAL(maldi_features[0].getIntensity(), 7)
  
  TEST_EQUAL(maldi_features[1].getCharge(), 2)
  TEST_EQUAL(maldi_features[1].getIntensity(), 3)

  TEST_EQUAL(maldi_features[2].getCharge(), 1)
  TEST_EQUAL(maldi_features[2].getIntensity(), 10)

  TEST_EQUAL(maldi_features[3].getCharge(), 1)
  TEST_EQUAL(maldi_features[3].getIntensity(), 9)

  TEST_EQUAL(maldi_features[4].getCharge(), 2)
  TEST_EQUAL(maldi_features[4].getIntensity(), 1)

  TEST_EQUAL(maldi_features[5].getCharge(), 1)
  TEST_EQUAL(maldi_features[5].getIntensity(), 9)

  TEST_EQUAL(maldi_features[6].getCharge(), 2)
  TEST_EQUAL(maldi_features[6].getIntensity(), 1)
  
 
  for(FeatureMapSim::const_iterator fmIt = maldi_features.begin(); fmIt != maldi_features.end();
      ++fmIt)
  {
    std::cout << (*fmIt).getCharge() << " " << (*fmIt).getIntensity() << " " << (*fmIt).getPeptideIdentifications()[0].getHits()[0].getSequence().toString() << std::endl;
  }
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



