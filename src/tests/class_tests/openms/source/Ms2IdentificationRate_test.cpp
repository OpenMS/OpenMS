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


//construct PeptideHits
PeptideHit pep_hit1_t1;
pep_hit1_t1.setMetaValue("target_decoy", "target");
PeptideHit pep_hit1_t2;
pep_hit1_t2.setMetaValue("target_decoy", "target");
PeptideHit pep_hit2_d;
pep_hit2_d.setMetaValue("target_decoy", "decoy");
PeptideHit pep_hit_fdr;

//construct vectors of PeptideHits
std::vector<PeptideHit> pep_hits_target = {pep_hit1_t1, pep_hit1_t2};
std::vector<PeptideHit> pep_hits_decoy = {pep_hit2_d};
std::vector<PeptideHit> pep_hits_empty = {};
std::vector<PeptideHit> pep_hits_fdr = {pep_hit_fdr};

//construct Peptideidentification with PeptideHits
PeptideIdentification pep_id_target;
pep_id_target.setHits(pep_hits_target);
PeptideIdentification pep_id_decoy;
pep_id_decoy.setHits(pep_hits_decoy);
PeptideIdentification pep_id_empty;
pep_id_empty.setHits(pep_hits_empty);
PeptideIdentification pep_id_fdr;
pep_id_fdr.setHits(pep_hits_fdr);

std::vector<PeptideIdentification> pep_ids = {pep_id_target, pep_id_decoy, pep_id_empty};
std::vector<PeptideIdentification> pep_ids_empty{};
std::vector<PeptideIdentification> pep_ids_fdr = {pep_id_fdr};

//construct features with peptideIdentifications
Feature feat_empty_pi;
feat_empty_pi.setPeptideIdentifications(pep_ids_empty);
Feature feat_target;
feat_target.setPeptideIdentifications(pep_ids);
Feature feat_empty;
Feature feat_fdr;
feat_fdr.setPeptideIdentifications(pep_ids_fdr);

//construct FeatureMap
FeatureMap fmap;
fmap.push_back(feat_empty_pi);
fmap.push_back(feat_target);
fmap.push_back(feat_empty);

FeatureMap fmap_fdr;
fmap_fdr.push_back(feat_fdr);

FeatureMap fmap_empty;


fmap.setUnassignedPeptideIdentifications(pep_ids);

//construct MSSpectrum
MSSpectrum ms2;
ms2.setMSLevel(2);
MSSpectrum ms1;
ms1.setMSLevel(1);
std::vector<MSSpectrum> ms_spectra = {ms2, ms2, ms2, ms2, ms2, ms2, ms1};
std::vector<MSSpectrum> ms1_spectra = {ms1};
std::vector<MSSpectrum> ms2_2_spectra = {ms2};

//construct MSExperiment
MSExperiment ms_exp;
ms_exp.setSpectra(ms_spectra);

//construct MSExperiment without MS2 spectra
MSExperiment ms1_exp;
ms1_exp.setSpectra(ms1_spectra);

//construct MSExperiment with two MS2 spectra
MSExperiment ms2_2_exp;
ms2_2_exp.setSpectra(ms2_2_spectra);

//construct empty MSExperiment
MSExperiment ms_empty_exp;



//////////////////////////////////////////////////////////////////
//start Section
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

//tests compute function
START_SECTION(void compute(FeatureMap const & feature_map, MSExperiment const & exp, bool force_fdr = false))
{
  //test with valid input
  ms2ir.compute(fmap, ms_exp);
  std::vector<Ms2IdentificationRate::IdentificationRateData> result;
  result = ms2ir.getResults();

  for (const auto& idrd : result)
  {
    TEST_EQUAL(idrd.num_peptide_identification, 2)
    TEST_EQUAL(idrd.num_ms2_spectra, 6)
    TEST_REAL_SIMILAR(idrd.identification_rate, 1./3)
  }

  //less ms2 spectra than identifictions
  TEST_EXCEPTION_WITH_MESSAGE(Exception::Precondition, ms2ir_ms2_2.compute(fmap, ms2_2_exp), "There are more Identifications than MS2 spectra. Please check your data.")

  //empty ms experiment
  TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, ms2ir_empty_msexp.compute(fmap, ms_empty_exp), "MSExperiment is empty")

  //empty feature map
  ms2ir_empty_fmap.compute(fmap_empty, ms_exp);
  std::vector<Ms2IdentificationRate::IdentificationRateData> result_empty_fmap;
  result_empty_fmap = ms2ir_empty_fmap.getResults();

  for (const auto& idrd_empty_fmap : result_empty_fmap)
  {
    TEST_EQUAL(idrd_empty_fmap.num_peptide_identification, 0)
    TEST_EQUAL(idrd_empty_fmap.num_ms2_spectra, 6)
    TEST_REAL_SIMILAR(idrd_empty_fmap.identification_rate, 0)
  }

  //no fdr
  TEST_EXCEPTION_WITH_MESSAGE(Exception::Precondition, ms2ir_fdr.compute(fmap_fdr, ms_exp), "FDR was not made. If you want to continue without FDR use -MS2_id_rate:force_no_fdr")

  // force no fdr
  ms2ir_force_fdr.compute(fmap_fdr, ms_exp, true);
  std::vector<Ms2IdentificationRate::IdentificationRateData> result_force_fdr;
  result_force_fdr = ms2ir_force_fdr.getResults();

  for (const auto& idrd_force_fdr : result_force_fdr)
  {
    TEST_EQUAL(idrd_force_fdr.num_peptide_identification, 1)
    TEST_EQUAL(idrd_force_fdr.num_ms2_spectra, 6)
    TEST_REAL_SIMILAR(idrd_force_fdr.identification_rate, 1./6)
  }

  //no ms2 spectra
  TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, ms2ir_ms1.compute(fmap, ms1_exp), "No MS2 spectra found")
}
END_SECTION


START_SECTION(const String& getName() const override)
{
  TEST_EQUAL(ms2ir.getName(), "Ms2IdentificationRate")
}
END_SECTION


START_SECTION(QCBase::Status requires() const override)
{
  QCBase::Status stat = QCBase::Status() | QCBase::Requires::RAWMZML | QCBase::Requires::POSTFDRFEAT;
  TEST_EQUAL(stat, ms2ir.requires())
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST