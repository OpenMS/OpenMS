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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/SIMULATION/LABELING/ITRAQLabeler.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/ListUtils.h>

using namespace OpenMS;
using namespace std;

START_TEST(ITRAQLabeler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ITRAQLabeler* ptr = nullptr;
ITRAQLabeler* null_ptr = nullptr;
START_SECTION(ITRAQLabeler())
{
	ptr = new ITRAQLabeler();
	TEST_NOT_EQUAL(ptr, null_ptr)
  TEST_EQUAL(ptr->getParameters().getValue("iTRAQ"), "4plex")
}
END_SECTION

START_SECTION(virtual ~ITRAQLabeler())
{
	delete ptr;
}
END_SECTION

START_SECTION((void preCheck(Param &param) const ))
{
  ITRAQLabeler i;
  Param p;
  p.setValue("RawTandemSignal:status", "MS^E");
  TEST_EXCEPTION(Exception::InvalidParameter, i.preCheck(p));

  p.setValue("RawTandemSignal:status", "precursor");
  i.preCheck(p); // should work
}
END_SECTION

START_SECTION((void setUpHook(SimTypes::FeatureMapSimVector &)))
{
  ITRAQLabeler i;
  // check for correct number of channels
  SimTypes::FeatureMapSimVector f_maps;
  f_maps.push_back(FeatureMap());
  i.setUpHook(f_maps);

  // add another map
  Param p = i.getParameters();
  p.setValue("channel_active_4plex", ListUtils::create<String>("114:myReference, 117:blabla"), "Four-plex only: Each channel that was used in the experiment and its description (114-117) in format <channel>:<name>, e.g. \"114:myref\",\"115:liver\".");
  i.setParameters(p);
  f_maps.push_back(FeatureMap());
  i.setUpHook(f_maps);

  // if no Exception until here, all is good

  NOT_TESTABLE
}
END_SECTION

START_SECTION((void postDigestHook(SimTypes::FeatureMapSimVector &)))
{
  ITRAQLabeler i;

  SimTypes::FeatureMapSimVector f_maps;
  FeatureMap fm1, fm2, fm3;

  // create peptide
  PeptideHit pep_hit(1.0, 1, 0, AASequence::fromString("AAHJK"));
  PeptideEvidence pe1;
  pe1.setProteinAccession("p1");
  pep_hit.setPeptideEvidences(vector<PeptideEvidence>(1, pe1));
  PeptideIdentification pep_id;
  pep_id.insertHit(pep_hit);
  // --
  PeptideHit pep_hit2(1.0, 1, 0, AASequence::fromString("EEEEPPPK"));
  PeptideEvidence pe2;
  pe2.setProteinAccession("p2");
  pep_hit2.setPeptideEvidences(vector<PeptideEvidence>(1, pe2));

  PeptideIdentification pep_id2;
  pep_id2.insertHit(pep_hit2);
  // --
  PeptideHit pep_hit3(1.0, 1, 0, AASequence::fromString("EEEEPPPK")); // same peptide as #2, but from different protein
  PeptideEvidence pe3;
  pe3.setProteinAccession("p3");
  pep_hit3.setPeptideEvidences(vector<PeptideEvidence>(1, pe3));
  PeptideIdentification pep_id3;
  pep_id3.insertHit(pep_hit3);

  // generate Feature
  Feature f1;
  f1.getPeptideIdentifications().push_back(pep_id);
  fm1.push_back(f1);
  fm2.push_back(f1);

  // generate Feature
  Feature f2;
  f2.getPeptideIdentifications().push_back(pep_id2);
  fm3.push_back(f2);

  // generate Feature
  Feature f3;
  f3.getPeptideIdentifications().push_back(pep_id3);
  fm3.push_back(f3);

  // merge
  f_maps.push_back(fm1);
  f_maps.push_back(fm2);
  f_maps.push_back(fm3);

  i.postDigestHook(f_maps);

  // one merged map
  TEST_EQUAL(f_maps.size(), 1)

  TEST_EQUAL(f_maps[0].size(), 2)

  set<String> protein_accessions1 = f_maps[0][0].getPeptideIdentifications()[0].getHits()[0].extractProteinAccessionsSet();
  TEST_EQUAL(protein_accessions1.size(), 1)

  set<String> protein_accessions2 = f_maps[0][1].getPeptideIdentifications()[0].getHits()[0].extractProteinAccessionsSet();
  TEST_EQUAL(protein_accessions2.size(), 2)
}
END_SECTION

START_SECTION((void postRTHook(SimTypes::FeatureMapSimVector &)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void postDetectabilityHook(SimTypes::FeatureMapSimVector &)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void postIonizationHook(SimTypes::FeatureMapSimVector &)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void postRawMSHook(SimTypes::FeatureMapSimVector &)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void postRawTandemMSHook(SimTypes::FeatureMapSimVector &, SimTypes::MSSimExperiment &)))
{
  ITRAQLabeler i;
  SimTypes::MutableSimRandomNumberGeneratorPtr rnd_gen (new SimTypes::SimRandomNumberGenerator);
  rnd_gen->initialize(false, false);
  i.setRnd(rnd_gen);

  SimTypes::FeatureMapSimVector f_maps;
  FeatureMap fm1;

  SimTypes::MSSimExperiment exp;
  MSSpectrum spec;
  IntList il;
  il.push_back(0);
  spec.setMetaValue("parent_feature_ids",  il);
  spec.setRT(600);
  spec.setMSLevel(2);
  exp.addSpectrum(spec);

  SimTypes::MSSimExperiment exp2=exp;

  std::vector<double> eb(4);
  DoubleList elution_bounds(eb);
  elution_bounds[0] = 100; elution_bounds[1] = 509.2; elution_bounds[2] = 120; elution_bounds[3] = 734.3;
  std::vector<double> ei(5, 0.5); // 50% elution profile
  DoubleList elution_ints(ei);
  Feature f;
  f.setMetaValue("elution_profile_bounds", elution_bounds);
  f.setMetaValue("elution_profile_intensities", elution_ints);
  f.setIntensity(100); // should result in 100 * 0.5 = 50 intensity
  f.setMZ(400);
  f.setRT(601);
  f.getConvexHull().addPoint(DPosition<2>(509.2, 398));
  f.getConvexHull().addPoint(DPosition<2>(734.3, 402));
  f.setMetaValue(i.getChannelIntensityName(0), 100);
  f.setMetaValue(i.getChannelIntensityName(1), 100);
  f.setMetaValue(i.getChannelIntensityName(2), 100);
  f.setMetaValue(i.getChannelIntensityName(3), 100);


  fm1.push_back(f);

  f_maps.push_back(fm1);
  Param p;
  p = i.getParameters();
  // no isotope skewing
  StringList iso = ListUtils::create<String>("114:0/0/100/0,115:0/0/0/0,116:0/0/0/0,117:0/100/0/0");
  p.setValue("isotope_correction_values_4plex", iso);
  StringList ch = ListUtils::create<String>("114:c1,115:c2,116:c3,117:c4");
  p.setValue("channel_active_4plex", ch);
  p.setValue("iTRAQ", "4plex");
  i.setParameters(p);
  i.postRawTandemMSHook(f_maps, exp);

  TEST_EQUAL(exp.size(), 1)
  Size count(0);
  double expected_val4[4] = {0, 100, 100, 0};
  for (MSSpectrum::const_iterator it=exp[0].begin(); it!=exp[0].end() && it->getMZ()<118.0; ++it)
  {
    TEST_REAL_SIMILAR(it->getIntensity(), expected_val4[count]);
    ++count;
}
  TEST_EQUAL(count, 4)
  exp=exp2;//revert

  // with isotope skewing
  iso = ListUtils::create<String>("113:0/0/100/0,"
                           "114:0/0/50 /0,"
                           "115:0/100/0/0,"
                           "116:0/0/100/0,"
                           "117:0/0/0/100,"
                           "118:0/0/100/0,"
                           "119:0/0/100/0,"
                           "121:0/100/0/0");
  p.setValue("isotope_correction_values_8plex", iso);
  ch = ListUtils::create<String>("113:ch0,114:c1,115:c2,116:c3,117:c4,118:c5,119:c6,121:c7");
  p.setValue("channel_active_8plex", ch);
  p.setValue("iTRAQ", "8plex");
  i.setParameters(p);

  f.setMetaValue(i.getChannelIntensityName(4), 100);
  f.setMetaValue(i.getChannelIntensityName(5), 100);
  f.setMetaValue(i.getChannelIntensityName(6), 100);
  f.setMetaValue(i.getChannelIntensityName(7), 100);
  fm1.clear();
  fm1.push_back(f);
  f_maps.clear();
  f_maps.push_back(fm1);

  i.postRawTandemMSHook(f_maps, exp);

  TEST_EQUAL(exp.size(), 1)
  count=0;
  double expected_val8[8] = {0, 125, 25, 0, 50, 0, 100, 0};
  MSSpectrum::const_iterator it=exp[0].begin();
  for (; it!=exp[0].end(); ++it)
  {
    TEST_REAL_SIMILAR(it->getIntensity(),  expected_val8[count]);
    ++count;
  }
  TEST_EQUAL(count, 8)
}
END_SECTION

START_SECTION((static BaseLabeler* create()))
{
  BaseLabeler* labeler = ITRAQLabeler::create();
  BaseLabeler* nullPointer = nullptr;
  TEST_NOT_EQUAL(labeler, nullPointer)
  delete labeler;
}
END_SECTION

START_SECTION((static const String getProductName()))
{
  ITRAQLabeler i;
  TEST_EQUAL(i.getProductName(), "itraq")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
