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
// $Maintainer: Chris Bielow $
// $Authors: Swenja Wagner, Patricia Scheil $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/QC/QCBase.h>
#include <OpenMS/QC/MissedCleavages.h>
#include <string>
#include <vector>

//////////////////////////

using namespace OpenMS;

START_TEST(MissedCleavages, "$Id$")
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


  //construct PeptideHits
  PeptideHit pep_hit_no_cut;
  pep_hit_no_cut.setSequence(AASequence::fromString("AAAAAAAAAAAAAK"));
  PeptideHit pep_hit_one_cut;
  pep_hit_one_cut.setSequence(AASequence::fromString("AAAAKAAAAAR"));
  PeptideHit pep_hit_three_cuts;
  pep_hit_three_cuts.setSequence(AASequence::fromString("AAAKAAARAAAKAAAR"));

  //construct vectors of PeptideHits
  std::vector<PeptideHit> pep_hits_0 = {pep_hit_no_cut};
  std::vector<PeptideHit> pep_hits_1 = {pep_hit_one_cut};
  std::vector<PeptideHit> pep_hits_3 = {pep_hit_three_cuts};
  std::vector<PeptideHit> pep_hits_empty = {};

  //construct PeptideIdentification with PeptideHits
  PeptideIdentification pep_id_0;
  pep_id_0.setHits(pep_hits_0);
  PeptideIdentification pep_id_1;
  pep_id_1.setHits(pep_hits_1);
  PeptideIdentification pep_id_empty;
  pep_id_empty.setHits(pep_hits_empty);
  PeptideIdentification pep_id_3;
  pep_id_3.setHits(pep_hits_3);

  //construct vectors of PeptideIdentifications
  std::vector<PeptideIdentification> pep_ids = {pep_id_0, pep_id_1, pep_id_empty};
  std::vector<PeptideIdentification> pep_ids_1 = {pep_id_1, pep_id_1};
  std::vector<PeptideIdentification> pep_ids_empty{};
  std::vector<PeptideIdentification> pep_ids_3 = {pep_id_3};

  //construct features with peptideIdentifications
  Feature feat_empty_pi;
  feat_empty_pi.setPeptideIdentifications(pep_ids_empty);
  Feature feat;
  feat.setPeptideIdentifications(pep_ids);
  Feature feat_empty;
  Feature feat_3;
  feat_3.setPeptideIdentifications(pep_ids_3);

  //FeatureMap
  FeatureMap feature_map;
  FeatureMap feature_map_3;
  FeatureMap feature_map_empty;
  FeatureMap feature_map_no_protein;
  FeatureMap feature_map_no_enzyme;

  //stores data in the FeatureMaps
  feature_map.push_back(feat_empty_pi);
  feature_map.push_back(feat);
  feature_map.push_back(feat_empty);
  feature_map.setUnassignedPeptideIdentifications(pep_ids_1);
  feature_map.getProteinIdentifications().resize(1);
  feature_map.getProteinIdentifications()[0].getSearchParameters().digestion_enzyme = *ProteaseDB::getInstance() -> getEnzyme("trypsin");
  feature_map.getProteinIdentifications()[0].getSearchParameters().missed_cleavages = 2;

  //FeatureMap with more missed cleavages than allowed
  feature_map_3.push_back(feat_3);
  feature_map_3.getProteinIdentifications().resize(1);
  feature_map_3.getProteinIdentifications()[0].getSearchParameters().digestion_enzyme.setName("trypsin");
  feature_map_3.getProteinIdentifications()[0].getSearchParameters().missed_cleavages = 2;

  //FeatureMap without ProteinIdentifications
  feature_map_no_protein.push_back(feat);

  //FeatureMap without given enzyme
  feature_map_no_enzyme.push_back(feat);
  ProteinIdentification prot_id;
  feature_map_no_enzyme.setProteinIdentifications({prot_id});

  //////////////////////////////////////////////////////////////////
  //start Section
  /////////////////////////////////////////////////////////////////

  MissedCleavages* ptr = nullptr;
  MissedCleavages* nulpt = nullptr;
  START_SECTION(MissedCleavges())
  {
  ptr = new MissedCleavages();
  TEST_NOT_EQUAL(ptr, nulpt)
  }
  END_SECTION


  START_SECTION(~MissedCleavages())
  {
  delete ptr;
  }
  END_SECTION


  //tests compute function
  START_SECTION(void compute(FeatureMap& fmap))
  {
    //test with valid input
    MissedCleavages mc;
    mc.compute(feature_map);
    auto result = mc.getResults();

    TEST_EQUAL(result[0].size(),2)
    TEST_EQUAL(result[0][0], 1)
    TEST_EQUAL(result[0][1], 3)

    std::vector<UInt32> frequ;

    //test if result is stored as MetaInformation in PeptidHits in FeatureMap
    auto lam = [&frequ](PeptideIdentification& pep_id)
    {
      if(pep_id.getHits().empty())
      {
        return;
      }
      if (pep_id.getHits()[0].metaValueExists("missed_cleavages"))
      {
        frequ.push_back(pep_id.getHits()[0].getMetaValue("missed_cleavages"));
      }
    };

    OpenMS::QCBase::iterateFeatureMap(feature_map, lam);
    TEST_EQUAL(frequ.size(), 4)
    TEST_EQUAL(frequ[0], 1)
    TEST_EQUAL(frequ[1], 1)
    TEST_EQUAL(frequ[2], 0)
    TEST_EQUAL(frequ[3], 1)

    //empty feature map
    MissedCleavages mc_empty;
    mc_empty.compute(feature_map_empty);
    auto result_empty = mc_empty.getResults();

    TEST_EQUAL(result_empty[0].empty(), true)

    //Missing informations in ProteinIdentifications
    //fmap.getProteinIdentifications().empty()
    MissedCleavages mc_no_protein;
    TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, mc_no_protein.compute(feature_map_no_protein), "Missing information in ProteinIdentifications.")

    //no given enzyme
    //enzyme == "unknown_enzyme"
    MissedCleavages mc_no_enzyme;
    TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, mc_no_enzyme.compute(feature_map_no_enzyme), "No digestion enzyme in FeatureMap detected. No computation possible.")

    //Number of missed cleavages is greater than the allowed maximum number of missed cleavages.
    MissedCleavages mc_3;
    mc_3.compute(feature_map_3);
    auto result_3 = mc_3.getResults();

    TEST_EQUAL(result_3[0].size(),1)
    TEST_EQUAL(result_3[0][3], 1)
  }
  END_SECTION
  

  MissedCleavages mc;

  START_SECTION(const String& getName() const override)
  {
    TEST_EQUAL(mc.getName(), "MissedCleavages")
  }
  END_SECTION
  

  START_SECTION(QCBase::Status requires() const override)
  {
    QCBase::Status stat = QCBase::Status() | QCBase::Requires::POSTFDRFEAT;
    TEST_EQUAL(stat, mc.requires())
  }
  END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST