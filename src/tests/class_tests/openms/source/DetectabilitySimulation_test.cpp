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
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/SIMULATION/DetectabilitySimulation.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/ListUtils.h>

using namespace OpenMS;
using namespace std;

START_TEST(DetectabilitySimulation, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DetectabilitySimulation* ptr = nullptr;
DetectabilitySimulation* nullPointer = nullptr;
START_SECTION(DetectabilitySimulation())
{
  ptr = new DetectabilitySimulation();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~DetectabilitySimulation())
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

START_SECTION((void filterDetectability(SimTypes::FeatureMapSim & features)))
{
  // test no detect
  DetectabilitySimulation detect_off;
  Param p = detect_off.getParameters();
  p.setValue("dt_simulation_on","false");
  p.setValue("min_detect", 0.9);
  detect_off.setParameters(p);

  SimTypes::FeatureMapSim no_detect_features;
  StringList peps = ListUtils::create<String>("TVQMENQFVAFVDK,ACHKKKKHHACAC,AAAAHTKLRTTIPPEFG,RYCNHKTUIKL");
  for (StringList::const_iterator it=peps.begin(); it!=peps.end(); ++it)
  {
    Feature f;
    PeptideIdentification pep_id;
    pep_id.insertHit(PeptideHit(1.0, 1, 1, AASequence::fromString(*it)));
    f.getPeptideIdentifications().push_back(pep_id);
    f.setIntensity(10);
    no_detect_features.push_back(f);
  }

  detect_off.filterDetectability(no_detect_features);

  TEST_EQUAL(no_detect_features.size(), 4)
  for(Size i = 0 ; i < no_detect_features.size() ; ++i)
  {
    TEST_EQUAL(no_detect_features[i].getMetaValue("detectability"), 1.0)
  }

  // test svm
  DetectabilitySimulation detect_svm;
  Param svm_params = detect_svm.getParameters();
  svm_params.setValue("dt_simulation_on","true");
  svm_params.setValue("min_detect", 0.4);
  svm_params.setValue("dt_model_file",OPENMS_GET_TEST_DATA_PATH("DetectabilitySimulation.svm"));
  detect_svm.setParameters(svm_params);

  SimTypes::FeatureMapSim svm_features;
  for (StringList::const_iterator it=peps.begin(); it!=peps.end(); ++it)
  {
    Feature f;
    PeptideIdentification pep_id;
    pep_id.insertHit(PeptideHit(1.0, 1, 1, AASequence::fromString(*it)));
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
  for(SimTypes::FeatureMapSim::const_iterator it = svm_features.begin() ; it != svm_features.end();
      ++it)
  {
    std::cout << (*it).getPeptideIdentifications()[0].getHits()[0].getSequence().toString()  << " " << (*it).getMetaValue("detectibility") << std::endl;
  }
  */
}
END_SECTION

START_SECTION((void predictDetectabilities(std::vector<String>& peptides_vector,std::vector<double>& labels, std::vector<double>& detectabilities)))
{
  // this method is called by "filterDetectability" so we already test it
  NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



