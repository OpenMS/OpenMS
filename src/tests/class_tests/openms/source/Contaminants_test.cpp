// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Authors: Dominik Schmitz, Chris Bielow$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/QC/Contaminants.h>
//#include <OpenMS/FORMAT/FeatureXMLFile.h>
//#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>

///////////////////////////

using namespace OpenMS;
using namespace std;


START_TEST(Contaminants, "$Id$")





FeatureMap fmap;
FeatureMap emptyFmap;
FeatureMap emptyFeaturesFmap;
Feature f;
Feature emptyf;
PeptideIdentification emptyId;
PeptideHit emptyHit;
//emptyHit.setSequence(AASequence::fromString(""));
//emptyId.setHits(std::vector<PeptideHit>());
//Emptyf.setPeptideIdentifications({emptyId});

emptyFeaturesFmap.push_back(emptyf);


fmap.getProteinIdentifications().resize(1);
DigestionEnzymeProtein noenzyme("unknown_enzyme",
                         "",
                         set<String>(),
                         "");
//set no digestion enzyme
fmap.getProteinIdentifications()[0].getSearchParameters().digestion_enzyme = noenzyme;
//set empty contaminants database
vector<FASTAFile::FASTAEntry> contaminantsFile;

    //fill the featureMap of features with set sequence and intensity
    {
      PeptideIdentification id;

      PeptideHit hit;

      hit.setSequence(AASequence::fromString("AAAAAAAAAAK"));
      id.setHits(std::vector<PeptideHit>(1, hit));
      f.setPeptideIdentifications({id});
      f.setIntensity(12.0);
      fmap.push_back(f);

      hit.setSequence(AASequence::fromString("R"));
      id.setHits(std::vector<PeptideHit>(1, hit));
      f.setPeptideIdentifications({id});
      f.setIntensity(8.0);
      fmap.push_back(f);

      hit.setSequence(AASequence::fromString("R"));
      id.setHits(std::vector<PeptideHit>(1, hit));
      f.setPeptideIdentifications({id});
      f.setIntensity(10.0);
      fmap.push_back(f);

      hit.setSequence(AASequence::fromString("AAAAAAAAAAKR"));
      id.setHits(std::vector<PeptideHit>(1, hit));
      f.setPeptideIdentifications({id});
      f.setIntensity(20.0);
      fmap.push_back(f);

      hit.setSequence(AASequence::fromString("QQQQQQQQQQ"));
      id.setHits(std::vector<PeptideHit>(1, hit));
      f.setPeptideIdentifications({id});
      f.setIntensity(10.0);
      fmap.push_back(f);
    }

//fill the unassigned peptideidentifications
std::vector<PeptideIdentification> ids2(3);
PeptideHit hit2;
hit2.setSequence(AASequence::fromString("AAAAAAAAAAK"));
ids2[0].setHits(std::vector<PeptideHit>(1, hit2));
hit2.setSequence(AASequence::fromString("RCCCCCCCCCCK"));
ids2[1].setHits(std::vector<PeptideHit>(1, hit2));
hit2.setSequence(AASequence::fromString("DDDDDDDDDD"));
ids2[2].setHits(std::vector<PeptideHit>(1, hit2));
fmap.setUnassignedPeptideIdentifications(ids2);


    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////


//check the constructor
Contaminants* ptr = nullptr;
Contaminants* nullPointer = nullptr;
START_SECTION(Contaminants())
  ptr = new Contaminants();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~Contaminants())
  delete ptr;
END_SECTION

START_SECTION((void compute(FeatureMap& features, const std::vector<FASTAFile::FASTAEntry>& contaminants)))
{
  Contaminants conts1, conts2, conts3, conts4, conts5, conts6;

  //test if it aborts when contaminant database is empty
  conts1.compute(fmap, contaminantsFile);
  std::vector<Contaminants::ContaminantsSummary> result1 = conts1.getResults();
  ABORT_IF(!result1.empty());

  //set contaminant database "contaminantsFile"
  FASTAFile::FASTAEntry contaminantsProtein("test_protein", "protein consists of only Alanine or Cytosine", "AAAAAAAAAAKRAAAAAAAAAAKRCCCCCCCCCCKRCCCCCCCCCC");
  contaminantsFile.push_back(contaminantsProtein);

  //test if it aborts when featureMap is empty
  conts2.compute(emptyFmap, contaminantsFile);
  std::vector<Contaminants::ContaminantsSummary> result2 = conts2.getResults();
  ABORT_IF(!result2.empty());

  //test if it aborts when features are empty
  conts6.compute(emptyFeaturesFmap, contaminantsFile);
  std::vector<Contaminants::ContaminantsSummary> result6 = conts6.getResults();
  ABORT_IF(!result6.empty());

  //test without given missed cleavages and without given enzyme
  conts3.compute(fmap, contaminantsFile);
  std::vector<Contaminants::ContaminantsSummary> result3 = conts3.getResults();
  ABORT_IF(result3.size() != 1);
  TEST_REAL_SIMILAR(result3[0].assigned_contaminants_ratio, 0.0);
  TEST_REAL_SIMILAR(result3[0].assigned_contaminants_intensity, 0.0);
  TEST_REAL_SIMILAR(result3[0].unassigned_contaminants_ratio, 0.0);
  TEST_REAL_SIMILAR(result3[0].all_contaminants_ratio, 0.0);

  //set digestion enzyme to trypsin
  fmap.getProteinIdentifications()[0].getSearchParameters().digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme("trypsin");

  //test without given missed cleavages but with set enzyme
  conts4.compute(fmap, contaminantsFile);
  std::vector<Contaminants::ContaminantsSummary> result4 = conts4.getResults();
  ABORT_IF(result4.size() != 1);
  TEST_REAL_SIMILAR(result4[0].assigned_contaminants_ratio, 3/5.0);
  TEST_REAL_SIMILAR(result4[0].assigned_contaminants_intensity, 1/2.0);
  TEST_REAL_SIMILAR(result4[0].unassigned_contaminants_ratio, 1/3.0);
  TEST_REAL_SIMILAR(result4[0].all_contaminants_ratio, 4/8.0);
  TEST_EQUAL(fmap[0].getPeptideIdentifications()[0].getMetaValue("is_contaminant"), 1);
  TEST_EQUAL(fmap[1].getPeptideIdentifications()[0].getMetaValue("is_contaminant"), 1);
  TEST_EQUAL(fmap[2].getPeptideIdentifications()[0].getMetaValue("is_contaminant"), 1);
  TEST_EQUAL(fmap[3].getPeptideIdentifications()[0].getMetaValue("is_contaminant"), 0);
  TEST_EQUAL(fmap[4].getPeptideIdentifications()[0].getMetaValue("is_contaminant"), 0);
  TEST_EQUAL(fmap.getUnassignedPeptideIdentifications()[0].getMetaValue("is_contaminant"), 1);
  TEST_EQUAL(fmap.getUnassignedPeptideIdentifications()[1].getMetaValue("is_contaminant"), 0);
  TEST_EQUAL(fmap.getUnassignedPeptideIdentifications()[2].getMetaValue("is_contaminant"), 0);

  //set missed cleavages to 1
  fmap.getProteinIdentifications()[0].getSearchParameters().missed_cleavages = 1;

  //test with set missed cleavages and set enzyme
  conts5.compute(fmap, contaminantsFile);
  std::vector<Contaminants::ContaminantsSummary> result5 = conts5.getResults();
  ABORT_IF(result5.size() != 1);
  TEST_REAL_SIMILAR(result5[0].assigned_contaminants_ratio, 4/5.0);
  TEST_REAL_SIMILAR(result5[0].assigned_contaminants_intensity, 5/6.0);
  TEST_REAL_SIMILAR(result5[0].unassigned_contaminants_ratio, 2/3.0);
  TEST_REAL_SIMILAR(result5[0].all_contaminants_ratio, 6/8.0);
  TEST_EQUAL(fmap[0].getPeptideIdentifications()[0].getMetaValue("is_contaminant"), 1);
  TEST_EQUAL(fmap[1].getPeptideIdentifications()[0].getMetaValue("is_contaminant"), 1);
  TEST_EQUAL(fmap[2].getPeptideIdentifications()[0].getMetaValue("is_contaminant"), 1);
  TEST_EQUAL(fmap[3].getPeptideIdentifications()[0].getMetaValue("is_contaminant"), 1);
  TEST_EQUAL(fmap[4].getPeptideIdentifications()[0].getMetaValue("is_contaminant"), 0);
  TEST_EQUAL(fmap.getUnassignedPeptideIdentifications()[0].getMetaValue("is_contaminant"), 1);
  TEST_EQUAL(fmap.getUnassignedPeptideIdentifications()[1].getMetaValue("is_contaminant"), 1);
  TEST_EQUAL(fmap.getUnassignedPeptideIdentifications()[2].getMetaValue("is_contaminant"), 0);

}
END_SECTION

START_SECTION(Status requires() const)
{
  Contaminants temp;
  TEST_EQUAL(temp.requires(), QCBase::Status(QCBase::Requires::POSTFDRFEAT) | QCBase::Requires::CONTAMINANTS);
}
END_SECTION
    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
END_TEST
