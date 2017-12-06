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
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/SIMULATION/RTSimulation.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/ListUtils.h>

using namespace OpenMS;
using namespace std;

START_TEST(RTSimulation, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

const unsigned long rnd_gen_seed = 1;
RTSimulation* ptr = nullptr;
RTSimulation* nullPointer = nullptr;
SimTypes::MutableSimRandomNumberGeneratorPtr empty_rnd_gen (new SimTypes::SimRandomNumberGenerator);

START_SECTION((RTSimulation(const SimTypes::MutableSimRandomNumberGeneratorPtr random_generator)))
{
  ptr = new RTSimulation(empty_rnd_gen);
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~RTSimulation())
{
  delete ptr;
}
END_SECTION

START_SECTION((RTSimulation(const RTSimulation &source)))
{
  RTSimulation source(empty_rnd_gen);
  Param p = source.getParameters();
  p.setValue("total_gradient_time",4000.0);
  source.setParameters(p);

  RTSimulation target(source);
  TEST_EQUAL(source.getParameters(), target.getParameters())
  TEST_EQUAL(source.getGradientTime(), target.getGradientTime())
}
END_SECTION

START_SECTION((RTSimulation& operator=(const RTSimulation &source)))
{
  RTSimulation source(empty_rnd_gen);
  RTSimulation target(source);

  Param p = source.getParameters();
  p.setValue("total_gradient_time",4000.0);
  source.setParameters(p);

  TEST_NOT_EQUAL(source.getParameters(), target.getParameters())
  target = source;
  TEST_EQUAL(source.getParameters(), target.getParameters())
}
END_SECTION

START_SECTION(( void predictRT(SimTypes::FeatureMapSim & features) ))
{
  // is fully tested by the different EXTRA tests for HPLC w absolute, HPLC w relative, none HPLC (and hopefully soon CE)
  NOT_TESTABLE
}
END_SECTION

START_SECTION(([EXTRA] Prediction Test - HPLC with relative RTs))
{
  // init rng
  // init rng
  SimTypes::MutableSimRandomNumberGeneratorPtr rnd_gen (new SimTypes::SimRandomNumberGenerator);
  rnd_gen->setBiologicalRngSeed(rnd_gen_seed);
  rnd_gen->setTechnicalRngSeed(rnd_gen_seed);

  // rt svm
  RTSimulation svm_rt_sim(rnd_gen);
  Param svm_params = svm_rt_sim.getParameters();
  svm_params.setValue("rt_column","HPLC");
  svm_params.setValue("total_gradient_time",4000.0);
  svm_params.setValue("scan_window:min",0.0);
  svm_params.setValue("scan_window:max",4000.0);
  svm_params.setValue("HPLC:model_file",OPENMS_GET_TEST_DATA_PATH("RTSimulation.svm"));
  svm_params.setValue("auto_scale", "true");
  svm_params.setValue("variation:affine_offset", 0);
  svm_params.setValue("variation:feature_stddev", 0);

  svm_rt_sim.setParameters(svm_params);

  SimTypes::FeatureMapSim svm_rt_features;
  StringList peps = ListUtils::create<String>("TVQMENQFVAFVDK,ACHKKKKHHACAC,AAAAHTKLRTTIPPEFG,RYCNHKTUIKL");
  for (StringList::const_iterator it=peps.begin(); it!=peps.end(); ++it)
  {
    Feature f;
    PeptideIdentification pep_id;
    pep_id.insertHit(PeptideHit(1.0, 1, 1, AASequence::fromString(*it)));
    f.getPeptideIdentifications().push_back(pep_id);
    f.setIntensity(10);
    svm_rt_features.push_back(f);
  }

  SimTypes::MSSimExperiment experiment_rt;
  svm_rt_sim.predictRT(svm_rt_features);

  TEST_EQUAL(svm_rt_features.size(), 4)

  TEST_REAL_SIMILAR(svm_rt_features[0].getRT(), 234.247)
  TEST_EQUAL(svm_rt_features[0].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "TVQMENQFVAFVDK")

  TEST_REAL_SIMILAR(svm_rt_features[1].getRT(), 471.292)
  TEST_EQUAL(svm_rt_features[1].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "RYCNHKTUIKL")

  TEST_REAL_SIMILAR(svm_rt_features[2].getRT(), 934.046)
  TEST_EQUAL(svm_rt_features[2].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "AAAAHTKLRTTIPPEFG")

  TEST_REAL_SIMILAR(svm_rt_features[3].getRT(), 946.127)
  TEST_EQUAL(svm_rt_features[3].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "ACHKKKKHHACAC")

}
END_SECTION

START_SECTION((void createExperiment(SimTypes::MSSimExperiment & experiment)))
{
  // init rng
  SimTypes::MutableSimRandomNumberGeneratorPtr rnd_gen (new SimTypes::SimRandomNumberGenerator);
  rnd_gen->setBiologicalRngSeed(rnd_gen_seed);
  rnd_gen->setTechnicalRngSeed(rnd_gen_seed);

  // rt svm
  RTSimulation svm_rt_sim(rnd_gen);
  Param svm_params = svm_rt_sim.getParameters();
  svm_params.setValue("rt_column","HPLC");
  svm_params.setValue("total_gradient_time",4000.0);
  svm_params.setValue("scan_window:min",200.0);
  svm_params.setValue("scan_window:max",500.0);
  svm_params.setValue("sampling_rate",5.0);
  svm_params.setValue("HPLC:model_file",OPENMS_GET_TEST_DATA_PATH("RTSimulation.svm"));
  svm_params.setValue("auto_scale", "true");
  svm_params.setValue("variation:affine_offset", 0);
  svm_params.setValue("variation:feature_stddev", 0);

  svm_rt_sim.setParameters(svm_params);

  SimTypes::FeatureMapSim svm_rt_features;
  StringList peps = ListUtils::create<String>("TVQMENQFVAFVDK,RYCNHKTUIKL");
  for (StringList::const_iterator it=peps.begin(); it!=peps.end(); ++it)
  {
    Feature f;
    PeptideIdentification pep_id;
    pep_id.insertHit(PeptideHit(1.0, 1, 1, AASequence::fromString(*it)));
    f.getPeptideIdentifications().push_back(pep_id);
    f.setIntensity(10);
    svm_rt_features.push_back(f);
  }

  SimTypes::MSSimExperiment experiment_rt;
  svm_rt_sim.predictRT(svm_rt_features);
  svm_rt_sim.createExperiment(experiment_rt);

  TEST_EQUAL(svm_rt_features.size(), 2)

  TEST_REAL_SIMILAR(experiment_rt.getMinRT(), 200.0)
  TEST_REAL_SIMILAR(experiment_rt.getMaxRT(), 500.0)

  SimTypes::MSSimExperiment::ConstIterator it = experiment_rt.RTBegin(200.0);
  SimTypes::MSSimExperiment::CoordinateType current_rt = 200.0;
  SimTypes::MSSimExperiment::CoordinateType scan_intervall = 5.0;
  while(it != experiment_rt.RTEnd(500.0))
  {
    TEST_REAL_SIMILAR((*it).getRT(), current_rt)
    ++it;
    current_rt += scan_intervall;
  }
}
END_SECTION

START_SECTION(([EXTRA] Prediction Test - No RT column))
{
  // init rng
  SimTypes::MutableSimRandomNumberGeneratorPtr rnd_gen (new SimTypes::SimRandomNumberGenerator);
  rnd_gen->setBiologicalRngSeed(rnd_gen_seed);
  rnd_gen->setTechnicalRngSeed(rnd_gen_seed);

  // no rt scan
  RTSimulation no_rt_sim(rnd_gen);
  Param p = no_rt_sim.getParameters();
  p.setValue("rt_column","none");
  p.setValue("total_gradient_time",4000.0);
  no_rt_sim.setParameters(p);

  SimTypes::FeatureMapSim no_rt_features;
  StringList peps = ListUtils::create<String>("TVQMENQFVAFVDK,ACHKKKKHHACAC,AAAAHTKLRTTIPPEFG,RYCNHKTUIKL");
  for (StringList::const_iterator it=peps.begin(); it!=peps.end(); ++it)
  {
    Feature f;
    PeptideIdentification pep_id;
    pep_id.insertHit(PeptideHit(1.0, 1, 1, AASequence::fromString(*it)));
    f.getPeptideIdentifications().push_back(pep_id);
    f.setIntensity(10);
    no_rt_features.push_back(f);
  }

  SimTypes::MSSimExperiment experiment_no_rt;
  no_rt_sim.predictRT(no_rt_features);
  no_rt_sim.createExperiment(experiment_no_rt);
  TEST_EQUAL(experiment_no_rt.size(), 1);
  for(SimTypes::FeatureMapSim::const_iterator fIt = no_rt_features.begin(); fIt != no_rt_features.end();
      ++fIt)
  {
    TEST_EQUAL((*fIt).getRT(), -1);
  }
}
END_SECTION

START_SECTION(([EXTRA] Prediction Test - HPLC with absolute RTs))
{
  // init rng
  SimTypes::MutableSimRandomNumberGeneratorPtr rnd_gen (new SimTypes::SimRandomNumberGenerator);
  rnd_gen->setBiologicalRngSeed(rnd_gen_seed);
  rnd_gen->setTechnicalRngSeed(rnd_gen_seed);

  // absolute rt values
  // rt svm
  RTSimulation rt_sim(rnd_gen);
  Param abs_svm_params = rt_sim.getParameters();
  abs_svm_params.setValue("rt_column","HPLC");
  abs_svm_params.setValue("total_gradient_time",4000.0);
  abs_svm_params.setValue("scan_window:min",200.0);
  abs_svm_params.setValue("scan_window:max",3000.0);
  abs_svm_params.setValue("HPLC:model_file",OPENMS_GET_TEST_DATA_PATH("RTSimulation_absolut_rt.model"));
  abs_svm_params.setValue("auto_scale", "false");
  abs_svm_params.setValue("variation:affine_offset", 0);
  abs_svm_params.setValue("variation:feature_stddev", 0);

  rt_sim.setParameters(abs_svm_params);

  SimTypes::FeatureMapSim features;

  // 2070, 1470, 2310, 3150
  StringList abs_peps = ListUtils::create<String>("QEFEVMEDHAGTYGLGDR,KGHHEAEIKPLAQSHATK,STPTAEDVTAPLVDEGAPGK,LSLEFPSGYPYNAPTVK");

  for (StringList::const_iterator it=abs_peps.begin(); it!=abs_peps.end(); ++it)
  {
    Feature f;
    PeptideIdentification pep_id;
    pep_id.insertHit(PeptideHit(1.0, 1, 1, AASequence::fromString(*it)));
    f.getPeptideIdentifications().push_back(pep_id);
    f.setIntensity(10);
    features.push_back(f);
  }

  SimTypes::MSSimExperiment experiment_rt;
  rt_sim.predictRT(features);

  TEST_EQUAL(features.size(), 3)

  // KGHHEAEIKPLAQSHATK 1560.7
  TEST_REAL_SIMILAR(features[0].getRT(), 1560.7)
  TEST_EQUAL(features[0].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "KGHHEAEIKPLAQSHATK")

  // QEFEVMEDHAGTYGLGDR 2160.7
  TEST_REAL_SIMILAR(features[1].getRT(), 2160.7)
  TEST_EQUAL(features[1].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "QEFEVMEDHAGTYGLGDR")

  // STPTAEDVTAPLVDEGAPGK 2400.69
  TEST_REAL_SIMILAR(features[2].getRT(), 2400.69)
  TEST_EQUAL(features[2].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "STPTAEDVTAPLVDEGAPGK")
}
END_SECTION

START_SECTION(([EXTRA] Prediction Test - CE column))
{
  // TODO: check CE rts
}
END_SECTION

START_SECTION((void predictContaminantsRT(SimTypes::FeatureMapSim &)))
{
  // TODO
}
END_SECTION

START_SECTION((bool isRTColumnOn() const ))
{
  RTSimulation rt_sim(empty_rnd_gen);

  Param p = rt_sim.getParameters();
  p.setValue("rt_column","HPLC");
  rt_sim.setParameters(p);

  TEST_EQUAL(rt_sim.isRTColumnOn(), true);

  p.setValue("rt_column","none");
  rt_sim.setParameters(p);

  TEST_EQUAL(rt_sim.isRTColumnOn(), false);
}
END_SECTION

START_SECTION((SimTypes::SimCoordinateType getGradientTime() const ))
{
  RTSimulation rt_sim(empty_rnd_gen);

  Param p = rt_sim.getParameters();
  p.setValue("total_gradient_time",1000.0);
  rt_sim.setParameters(p);

  TEST_EQUAL(rt_sim.getGradientTime(), 1000.0);

  p.setValue("total_gradient_time",4000.0);
  rt_sim.setParameters(p);

  TEST_EQUAL(rt_sim.getGradientTime(), 4000.0);
}
END_SECTION

START_SECTION((void wrapSVM(std::vector<AASequence>& peptide_sequences,std::vector<double>& predicted_retention_times)))
{
  // this method is called by "predictRT" so we already test it
  NOT_TESTABLE
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
