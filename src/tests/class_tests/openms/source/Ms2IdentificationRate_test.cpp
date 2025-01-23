// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Swenja Wagner, Patricia Scheil $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/QC/Ms2IdentificationRate.h>
#include <vector>

//////////////////////////

using namespace OpenMS;


START_TEST(Ms2IdentificationRate, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


// construct PeptideHits
PeptideHit pep_hit1_t1;
pep_hit1_t1.setMetaValue("target_decoy", "target");
PeptideHit pep_hit1_t2;
pep_hit1_t2.setMetaValue("target_decoy", "target");
PeptideHit pep_hit2_d;
pep_hit2_d.setMetaValue("target_decoy", "decoy");
PeptideHit pep_hit_fdr;

// construct vectors of PeptideHits
std::vector<PeptideHit> pep_hits_target = {pep_hit1_t1, pep_hit1_t2};
std::vector<PeptideHit> pep_hits_decoy = {pep_hit2_d};
std::vector<PeptideHit> pep_hits_empty = {};
std::vector<PeptideHit> pep_hits_fdr = {pep_hit_fdr};

// construct Peptideidentification with PeptideHits
PeptideIdentification pep_id_target;
pep_id_target.setHits(pep_hits_target);
PeptideIdentification pep_id_decoy;
pep_id_decoy.setHits(pep_hits_decoy);
PeptideIdentification pep_id_empty;
pep_id_empty.setHits(pep_hits_empty);
PeptideIdentification pep_id_fdr;
pep_id_fdr.setHits(pep_hits_fdr);

std::vector<PeptideIdentification> pep_ids = {pep_id_target, pep_id_decoy, pep_id_empty};
std::vector<PeptideIdentification> pep_ids_empty {};
std::vector<PeptideIdentification> pep_ids_fdr = {pep_id_fdr};

std::vector<PeptideIdentification> two_target_ids = {pep_id_target, pep_id_target, pep_id_decoy, pep_id_empty};

// construct features with peptideIdentifications
Feature feat_empty_pi;
feat_empty_pi.setPeptideIdentifications(pep_ids_empty);
Feature feat_target;
feat_target.setPeptideIdentifications(pep_ids);
Feature feat_empty;
Feature feat_fdr;
feat_fdr.setPeptideIdentifications(pep_ids_fdr);

// construct FeatureMap
FeatureMap fmap;
fmap.push_back(feat_empty_pi);
fmap.push_back(feat_target);
fmap.push_back(feat_empty);

FeatureMap fmap_fdr;
fmap_fdr.push_back(feat_fdr);

FeatureMap fmap_empty;


fmap.setUnassignedPeptideIdentifications(pep_ids);

// construct MSSpectrum
MSSpectrum ms2;
ms2.setMSLevel(2);
MSSpectrum ms1;
ms1.setMSLevel(1);
std::vector<MSSpectrum> ms_spectra = {ms2, ms2, ms2, ms2, ms2, ms2, ms1};
std::vector<MSSpectrum> ms1_spectra = {ms1};
std::vector<MSSpectrum> ms2_2_spectra = {ms2};

// construct MSExperiment
MSExperiment ms_exp;
ms_exp.setSpectra(ms_spectra);

// construct MSExperiment without MS2 spectra
MSExperiment ms1_exp;
ms1_exp.setSpectra(ms1_spectra);

// construct MSExperiment with two MS2 spectra
MSExperiment ms2_2_exp;
ms2_2_exp.setSpectra(ms2_2_spectra);

// construct empty MSExperiment
MSExperiment ms_empty_exp;


//////////////////////////////////////////////////////////////////
// start Section
/////////////////////////////////////////////////////////////////

Ms2IdentificationRate* ptr = nullptr;
Ms2IdentificationRate* nulpt = nullptr;
START_SECTION(Ms2IdentificationRate())
{
  ptr = new Ms2IdentificationRate();
  TEST_NOT_EQUAL(ptr, nulpt)
}
END_SECTION

START_SECTION(~Ms2IdentificationRate())
{
  delete ptr;
}
END_SECTION


Ms2IdentificationRate ms2ir;
Ms2IdentificationRate ms2ir_fdr;
Ms2IdentificationRate ms2ir_force_fdr;
Ms2IdentificationRate ms2ir_ms1;
Ms2IdentificationRate ms2ir_ms2_2;
Ms2IdentificationRate ms2ir_empty_msexp;
Ms2IdentificationRate ms2ir_empty_fmap;

// tests compute function with FeatureMap
START_SECTION(void compute(FeatureMap const& feature_map, MSExperiment const& exp, bool force_index = false))
{
  // test with valid input
  ms2ir.compute(fmap, ms_exp);
  std::vector<Ms2IdentificationRate::IdentificationRateData> result;
  result = ms2ir.getResults();

  for (const auto& idrd : result)
  {
    TEST_EQUAL(idrd.num_peptide_identification, 2)
    TEST_EQUAL(idrd.num_ms2_spectra, 6)
    TEST_REAL_SIMILAR(idrd.identification_rate, 1. / 3)
  }

  // less ms2 spectra than identifictions
  TEST_EXCEPTION_WITH_MESSAGE(Exception::Precondition, ms2ir_ms2_2.compute(fmap, ms2_2_exp), "There are more Identifications than MS2 spectra. Please check your data.")

  // empty ms experiment
  TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, ms2ir_empty_msexp.compute(fmap, ms_empty_exp), "MSExperiment is empty")

  // empty feature map
  ms2ir_empty_fmap.compute(fmap_empty, ms_exp);
  std::vector<Ms2IdentificationRate::IdentificationRateData> result_empty_fmap;
  result_empty_fmap = ms2ir_empty_fmap.getResults();

  for (const auto& idrd_empty_fmap : result_empty_fmap)
  {
    TEST_EQUAL(idrd_empty_fmap.num_peptide_identification, 0)
    TEST_EQUAL(idrd_empty_fmap.num_ms2_spectra, 6)
    TEST_REAL_SIMILAR(idrd_empty_fmap.identification_rate, 0)
  }

  // no fdr
  TEST_EXCEPTION_WITH_MESSAGE(Exception::Precondition, ms2ir_fdr.compute(fmap_fdr, ms_exp), "No target/decoy annotation found. If you want to continue regardless use -MS2_id_rate:assume_all_target")

  // force no fdr
  ms2ir_force_fdr.compute(fmap_fdr, ms_exp, true);
  std::vector<Ms2IdentificationRate::IdentificationRateData> result_force_fdr;
  result_force_fdr = ms2ir_force_fdr.getResults();

  for (const auto& idrd_force_fdr : result_force_fdr)
  {
    TEST_EQUAL(idrd_force_fdr.num_peptide_identification, 1)
    TEST_EQUAL(idrd_force_fdr.num_ms2_spectra, 6)
    TEST_REAL_SIMILAR(idrd_force_fdr.identification_rate, 1. / 6)
  }

  // no ms2 spectra
  TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, ms2ir_ms1.compute(fmap, ms1_exp), "No MS2 spectra found")
}
END_SECTION

Ms2IdentificationRate id_rate;
Ms2IdentificationRate id_rate_one_ms2;
Ms2IdentificationRate id_rate_empty_exp;
Ms2IdentificationRate id_rate_no_ids;
Ms2IdentificationRate id_rate_no_index;
Ms2IdentificationRate id_rate_force_no_index;
Ms2IdentificationRate id_rate_ms1;

// tests compute function with PeptideIdentifications
START_SECTION(void compute(const std::vector<PeptideIdentification>& pep_ids, const MSExperiment& exp, bool force_index = false))
{
  // test with valid input
  id_rate.compute(pep_ids, ms_exp);
  std::vector<Ms2IdentificationRate::IdentificationRateData> result;
  result = id_rate.getResults();

  for (const auto& idrd : result)
  {
    TEST_EQUAL(idrd.num_peptide_identification, 1)
    TEST_EQUAL(idrd.num_ms2_spectra, 6)
    TEST_REAL_SIMILAR(idrd.identification_rate, 1. / 6)
  }

  // less ms2 spectra than identifictions
  TEST_EXCEPTION_WITH_MESSAGE(Exception::Precondition, id_rate_one_ms2.compute(two_target_ids, ms2_2_exp), "There are more Identifications than MS2 spectra. Please check your data.")

  // empty ms experiment
  TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, id_rate_empty_exp.compute(pep_ids, ms_empty_exp), "MSExperiment is empty")

  // empty feature map
  id_rate_no_ids.compute(pep_ids_empty, ms_exp);
  std::vector<Ms2IdentificationRate::IdentificationRateData> result_no_pep_ids;
  result_no_pep_ids = id_rate_no_ids.getResults();

  for (const auto& idrd_empty_fmap : result_no_pep_ids)
  {
    TEST_EQUAL(idrd_empty_fmap.num_peptide_identification, 0)
    TEST_EQUAL(idrd_empty_fmap.num_ms2_spectra, 6)
    TEST_REAL_SIMILAR(idrd_empty_fmap.identification_rate, 0)
  }

  // no fdr
  TEST_EXCEPTION_WITH_MESSAGE(Exception::Precondition, id_rate_no_index.compute(pep_ids_fdr, ms_exp),
                              "No target/decoy annotation found. If you want to continue regardless use -MS2_id_rate:assume_all_target")

  // force no fdr
  id_rate_force_no_index.compute(pep_ids_fdr, ms_exp, true);
  std::vector<Ms2IdentificationRate::IdentificationRateData> result_force_fdr;
  result_force_fdr = id_rate_force_no_index.getResults();

  for (const auto& idrd_force_fdr : result_force_fdr)
  {
    TEST_EQUAL(idrd_force_fdr.num_peptide_identification, 1)
    TEST_EQUAL(idrd_force_fdr.num_ms2_spectra, 6)
    TEST_REAL_SIMILAR(idrd_force_fdr.identification_rate, 1. / 6)
  }

  // no ms2 spectra
  TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, id_rate_ms1.compute(pep_ids, ms1_exp), "No MS2 spectra found")
}
END_SECTION


START_SECTION(const String& getName() const override) {TEST_EQUAL(ms2ir.getName(), "Ms2IdentificationRate")} END_SECTION


  START_SECTION(QCBase::Status requirements() const override)
{
  QCBase::Status stat = QCBase::Status() | QCBase::Requires::RAWMZML | QCBase::Requires::POSTFDRFEAT;
  TEST_EQUAL(stat == ms2ir.requirements(), true)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
