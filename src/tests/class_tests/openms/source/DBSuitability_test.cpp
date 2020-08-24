// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Tom Waschischeck $
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/QC/DBSuitability.h>

#include <vector>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Suitability, "$Id$")

/////////////////////////////////////////////////////////////
////////////////////// CREATE DATA //////////////////////////
/////////////////////////////////////////////////////////////

PeptideEvidence decoy_protein("DECOY_PROT", 0, 0, 'A', 'A');
PeptideEvidence target_protein("DB_PROT", 0, 0, 'A', 'A');
PeptideEvidence novo_protein(Constants::UserParam::CONCAT_PEPTIDE, 0, 0, 'A', 'A');

// target-db hits with different q-values
PeptideHit target_db_hit1;
target_db_hit1.setSequence(AASequence::fromString("PEP"));
target_db_hit1.setPeptideEvidences({ target_protein });
target_db_hit1.setMetaValue("target_decoy", "target+decoy");
target_db_hit1.setMetaValue("MS:1002252", 0.8);
target_db_hit1.setScore(0.002);

PeptideHit target_db_hit2;
target_db_hit2.setSequence(AASequence::fromString("PEP"));
target_db_hit2.setPeptideEvidences({ target_protein });
target_db_hit2.setMetaValue("target_decoy", "target");
target_db_hit2.setMetaValue("MS:1002252", 0.8);
target_db_hit2.setScore(0.011);

// target-novo hits with different xcorr scores
PeptideHit target_novo_hit1;
target_novo_hit1.setSequence(AASequence::fromString("PEP"));
target_novo_hit1.setPeptideEvidences({ novo_protein });
target_novo_hit1.setMetaValue("target_decoy", "target");
target_novo_hit1.setMetaValue("MS:1002252", 0.85); // diff to db 0.05
target_novo_hit1.setScore(0.001);

PeptideHit target_novo_hit2;
target_novo_hit2.setSequence(AASequence::fromString("PEP"));
target_novo_hit2.setPeptideEvidences({ novo_protein });
target_novo_hit2.setMetaValue("target_decoy", "target");
target_novo_hit2.setMetaValue("MS:1002252", 0.93); // diff to db 0.13
target_novo_hit2.setScore(0.001);

// decoy hits with different xcorr scores -> resulting cut-offs: 0.15, 0.1, 0.05 (devided by 377.3 - weight of "PEP")
PeptideHit decoy1;
decoy1.setSequence(AASequence::fromString("PEP"));
decoy1.setPeptideEvidences({ decoy_protein });
decoy1.setMetaValue("target_decoy", "decoy");
decoy1.setMetaValue("MS:1002252", 0.7);
decoy1.setScore(1);

PeptideHit decoy2;
decoy2.setSequence(AASequence::fromString("PEP"));
decoy2.setPeptideEvidences({ decoy_protein });
decoy2.setMetaValue("target_decoy", "decoy");
decoy2.setMetaValue("MS:1002252", 0.6);
decoy2.setScore(1);

PeptideHit decoy3;
decoy3.setSequence(AASequence::fromString("PEP"));
decoy3.setMetaValue("target_decoy", "decoy");
decoy3.setMetaValue("MS:1002252", 0.55);
decoy3.setScore(1);

PeptideHit high_decoy;
high_decoy.setSequence(AASequence::fromString("PEP"));
high_decoy.setMetaValue("target_decoy", "decoy");
high_decoy.setMetaValue("MS:1002252", 0.55);
high_decoy.setScore(0);

// some error throwing hits

PeptideHit no_xcorr_hit;
no_xcorr_hit.setSequence(AASequence::fromString("PEP"));
no_xcorr_hit.setPeptideEvidences({ decoy_protein });
no_xcorr_hit.setMetaValue("target_decoy", "decoy");
no_xcorr_hit.setScore(1);

// build identifications
vector<PeptideIdentification> pep_ids;
vector<PeptideIdentification> top_decoy;
vector<PeptideIdentification> few_decoys;
vector<PeptideIdentification> no_xcorr_ids;
PeptideIdentification pep_id;
pep_id.setScoreType("some_score");
pep_id.setHigherScoreBetter(false);
pep_id.setHits({ target_novo_hit1, decoy1, decoy2 });
pep_ids.push_back(pep_id);
top_decoy.push_back(pep_id);

pep_id.setHits({ target_db_hit1, decoy1, decoy3 });
pep_ids.push_back(pep_id);
top_decoy.push_back(pep_id);

pep_id.setHits({ target_db_hit2 });
pep_ids.push_back(pep_id);
few_decoys.push_back(pep_id);

pep_id.setHits({ high_decoy, target_db_hit2 });
top_decoy.push_back(pep_id);

pep_id.setHits({ target_novo_hit1, target_db_hit1, decoy2, decoy3});
pep_ids.push_back(pep_id);
top_decoy.push_back(pep_id);
no_xcorr_ids.push_back(pep_id);

pep_id.setHits({ target_novo_hit2, target_db_hit1 });
pep_ids.push_back(pep_id);
top_decoy.push_back(pep_id);

pep_id.setHits({ no_xcorr_hit });
no_xcorr_ids.push_back(pep_id);

vector<PeptideIdentification> pep_ids_2(pep_ids);
vector<PeptideIdentification> pep_ids_3(pep_ids);

vector<PeptideIdentification> FDR_id;
pep_id.setScoreType("q-value");
pep_id.setHits({ decoy1 });
FDR_id.push_back(pep_id);

/////////////////////////////////////////////////////////////
///////////////////// START TESTING /////////////////////////
/////////////////////////////////////////////////////////////

DBSuitability* ptr = nullptr;
DBSuitability* nulpt = nullptr;
START_SECTION(DBSuitability())
{
  ptr = new DBSuitability();
  TEST_NOT_EQUAL(ptr, nulpt)
}
END_SECTION

START_SECTION(~DBSuitability())
{
  delete ptr;
}
END_SECTION

START_SECTION(void compute(vector<PeptideIdentification>& pep_ids))
{
  DBSuitability s;
  Param p;
  p.setValue("reranking_cutoff_percentile", 1.);
  s.setParameters(p);
  s.compute(pep_ids);

  p.setValue("reranking_cutoff_percentile", 1./3);
  p.setValue("FDR", 0.);
  s.setParameters(p);
  s.compute(pep_ids_2);
  s.compute(top_decoy);

  p.setValue("reranking_cutoff_percentile", 0.);
  s.setParameters(p);
  s.compute(pep_ids_3);
  vector<DBSuitability::SuitabilityData> d = s.getResults();
  DBSuitability::SuitabilityData data_fract_1 = d[0];
  DBSuitability::SuitabilityData data_fract_05 = d[1];
  DBSuitability::SuitabilityData data_decoy_top = d[2];
  DBSuitability::SuitabilityData data_small_percentile = d[3];
  TEST_REAL_SIMILAR(data_fract_1.cut_off, 0.00044);
  TEST_REAL_SIMILAR(data_fract_05.cut_off, 0.00029);
  TEST_REAL_SIMILAR(data_decoy_top.cut_off, 0.00029);
  TEST_REAL_SIMILAR(data_small_percentile.cut_off, 0.00014);
  TEST_EQUAL(data_fract_1.num_interest, 2);
  TEST_EQUAL(data_fract_05.num_interest, 2);
  TEST_EQUAL(data_decoy_top.num_interest, 0);
  TEST_EQUAL(data_small_percentile.num_interest, 2);
  TEST_EQUAL(data_fract_1.num_re_ranked, 2);
  TEST_EQUAL(data_fract_05.num_re_ranked, 1);
  TEST_EQUAL(data_decoy_top.num_re_ranked, 0);
  TEST_EQUAL(data_small_percentile.num_re_ranked, 0);
  TEST_EQUAL(data_fract_1.num_top_db, 4);
  TEST_EQUAL(data_fract_05.num_top_db, 3);
  TEST_EQUAL(data_decoy_top.num_top_db, 0);
  TEST_EQUAL(data_small_percentile.num_top_db, 2);
  TEST_EQUAL(data_fract_1.num_top_novo, 1);
  TEST_EQUAL(data_fract_05.num_top_novo, 2);
  TEST_EQUAL(data_decoy_top.num_top_novo, 0);
  TEST_EQUAL(data_small_percentile.num_top_novo, 3);
  TEST_REAL_SIMILAR(data_fract_1.suitability, 4./5);
  TEST_REAL_SIMILAR(data_fract_05.suitability, 3./5);
  TEST_REAL_SIMILAR(data_small_percentile.suitability, 2./5);
  TEST_EQUAL(data_decoy_top.suitability, DBL_MAX);

  TEST_EXCEPTION_WITH_MESSAGE(Exception::Precondition, s.compute(FDR_id), "q-value found at PeptideIdentifications. That is not allowed! Please make sure FDR did not run previously.");
  TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, s.compute(few_decoys), "Under 20 % of peptide identifications have two decoy hits. This is not enough for re-ranking. Use the 'no_rerank' flag to still compute a suitability score.");
  TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, s.compute(no_xcorr_ids), "No cross correlation score found at peptide hit. Only Comet search engine is supported right now.");
}
END_SECTION

START_SECTION(getResults())
{
  NOT_TESTABLE;
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST