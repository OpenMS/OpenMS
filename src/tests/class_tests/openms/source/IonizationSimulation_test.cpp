// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Chris Bielow$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/SIMULATION/IonizationSimulation.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
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
MutableSimRandomNumberGeneratorPtr rnd_gen (new SimRandomNumberGenerator);

// init reproducible rnd_gen
rnd_gen->initialize(false, false);

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
  MutableSimRandomNumberGeneratorPtr rnd_gen (new SimRandomNumberGenerator);
  rnd_gen->setBiologicalRngSeed(rnd_gen_seed);
  rnd_gen->setTechnicalRngSeed(rnd_gen_seed);

  // testing ESI
  IonizationSimulation esi_sim(rnd_gen);
  Param esi_param = esi_sim.getParameters();
  esi_param.setValue("ionization_type","ESI");
  esi_param.setValue("esi:ionized_residues",ListUtils::create<String>("Arg,Lys,His"));
  esi_param.setValue("esi:ionization_probability", 0.8);
  esi_param.setValue("esi:charge_impurity", ListUtils::create<String>("H+:1,NH4+:0.2,Ca++:0.1"));
  esi_param.setValue("esi:max_impurity_set_size", 3);

  esi_sim.setParameters(esi_param);

  FeatureMapSim esi_features;
  ConsensusMap cm;
  StringList peps = ListUtils::create<String>("TVQMENQFVAFVDK,ACHKKKKHHACAC,AAAAHTKLRTTIPPEFG,RYCNHKTUIKL");
  for (StringList::const_iterator it=peps.begin(); it!=peps.end(); ++it)
  {
    Feature f;
    PeptideIdentification pep_id;
    pep_id.insertHit(PeptideHit(1.0, 1, 1, AASequence(*it)));
    f.getPeptideIdentifications().push_back(pep_id);
    f.setIntensity(10);
    esi_features.push_back(f);
  }

  MSSimExperiment exp;
  MSSimExperiment::SpectrumType spec;
  exp.addSpectrum(spec);

  esi_sim.ionize(esi_features, cm, exp);

  TEST_EQUAL(esi_features.size(), 18)
  ABORT_IF(esi_features.size()!=18)

  TEST_EQUAL(esi_features[0].getCharge(), 2)
  TEST_EQUAL(esi_features[0].getIntensity(), 6)

  TEST_EQUAL(esi_features[1].getCharge(), 2)
  TEST_EQUAL(esi_features[1].getIntensity(), 2)

  TEST_EQUAL(esi_features[2].getCharge(), 3)
  TEST_EQUAL(esi_features[2].getIntensity(), 1)

  TEST_EQUAL(esi_features[3].getCharge(), 1)
  TEST_EQUAL(esi_features[3].getIntensity(), 1)

  TEST_EQUAL(esi_features[4].getCharge(), 7)
  TEST_EQUAL(esi_features[4].getIntensity(), 2)

  TEST_EQUAL(esi_features[5].getCharge(), 7)
  TEST_EQUAL(esi_features[5].getIntensity(), 2)

  TEST_EQUAL(esi_features[6].getCharge(), 6)
  TEST_EQUAL(esi_features[6].getIntensity(), 1)

  TEST_EQUAL(esi_features[7].getCharge(), 4)
  TEST_EQUAL(esi_features[7].getIntensity(), 3)

  TEST_EQUAL(esi_features[8].getCharge(), 3)
  TEST_EQUAL(esi_features[8].getIntensity(), 2)

  TEST_EQUAL(esi_features[9].getCharge(), 4)
  TEST_EQUAL(esi_features[9].getIntensity(), 1)

  TEST_EQUAL(esi_features[10].getCharge(), 4)
  TEST_EQUAL(esi_features[10].getIntensity(), 1)

  TEST_EQUAL(esi_features[11].getCharge(), 3)
  TEST_EQUAL(esi_features[11].getIntensity(), 1)

  TEST_EQUAL(esi_features[12].getCharge(), 3)
  TEST_EQUAL(esi_features[12].getIntensity(), 1)

  TEST_EQUAL(esi_features[13].getCharge(), 1)
  TEST_EQUAL(esi_features[13].getIntensity(), 1)

  TEST_EQUAL(esi_features[14].getCharge(), 4)
  TEST_EQUAL(esi_features[14].getIntensity(), 5)

  TEST_EQUAL(esi_features[15].getCharge(), 4)
  TEST_EQUAL(esi_features[15].getIntensity(), 3)

  TEST_EQUAL(esi_features[16].getCharge(), 3)
  TEST_EQUAL(esi_features[16].getIntensity(), 1)

	TEST_EQUAL(esi_features[17].getCharge(), 2)
  TEST_EQUAL(esi_features[17].getIntensity(), 1)

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


  MutableSimRandomNumberGeneratorPtr rnd_gen_maldi (new SimRandomNumberGenerator);
  rnd_gen_maldi->setBiologicalRngSeed(rnd_gen_seed);
  rnd_gen_maldi->setTechnicalRngSeed(rnd_gen_seed);

  // testing MALDI
  IonizationSimulation maldi_sim(rnd_gen_maldi);
  Param maldi_param = maldi_sim.getParameters();
  maldi_param.setValue("ionization_type","MALDI");
  maldi_param.setValue("maldi:ionization_probabilities", ListUtils::create<DoubleReal>("0.9,0.1"));

  maldi_sim.setParameters(maldi_param);

  FeatureMapSim maldi_features;
	for (StringList::const_iterator it=peps.begin(); it!=peps.end(); ++it)
	{
		Feature f;
		PeptideIdentification pep_id;
		pep_id.insertHit(PeptideHit(1.0, 1, 1, AASequence(*it)));
		f.getPeptideIdentifications().push_back(pep_id);
		f.setIntensity(10);
    maldi_features.push_back(f);
	}

	MSSimExperiment expt;
	MSSimExperiment::SpectrumType spect;
	expt.addSpectrum(spect);
	maldi_sim.ionize(maldi_features, cm, expt);

  TEST_EQUAL(maldi_features.size(), 6)

  TEST_EQUAL(maldi_features[0].getCharge(), 1)
  TEST_EQUAL(maldi_features[0].getIntensity(), 9)

  TEST_EQUAL(maldi_features[1].getCharge(), 2)
  TEST_EQUAL(maldi_features[1].getIntensity(), 1)

  TEST_EQUAL(maldi_features[2].getCharge(), 1)
  TEST_EQUAL(maldi_features[2].getIntensity(), 9)

  TEST_EQUAL(maldi_features[3].getCharge(), 2)
  TEST_EQUAL(maldi_features[3].getIntensity(), 1)

  TEST_EQUAL(maldi_features[4].getCharge(), 1)
  TEST_EQUAL(maldi_features[4].getIntensity(), 10)

  TEST_EQUAL(maldi_features[5].getCharge(), 1)
  TEST_EQUAL(maldi_features[5].getIntensity(), 10)

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



