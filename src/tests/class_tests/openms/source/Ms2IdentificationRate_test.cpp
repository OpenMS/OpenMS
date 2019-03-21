/ --------------------------------------------------------------------------
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

#include <include/OpenMS/QC/Ms2IdentificationRate.h>

#include <iostream>
//////////////////////////

using namespace OpenMS;



//START_TEST(Ms2IdentificationRate, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


//construct PeptideHits
PeptideHit pep_hit1_t1;
pep_hit_t1.setMetaValue("target_decoy", "target");
PeptideHit pep_hit1_t2;
pep_hit_t2.setMetaValue("target_decoy", "target");
PeptideHit pep_hit2_d;
pep_hit_t1.setMetaValue("target_decoy", "decoy");

//construct vectors of PeptideHits
std::vector<PeptideHit> pep_hits_target = {pep_hit1_t1, pep_hit1_t2};
std::vector<PeptideHit> pep_hits_decoy = {pep_hit2_d};
std::vector<PeptideHit> pep_hits_empty = {};

//construct Peptideidentification with PeptideHits
PeptideIdentification pep_id_target;
pep_id_target.setHit(pep_hits_target);
PeptideIdentification pep_id_decoy;
pep_id_target.setHit(pep_hits_decoy);
PeptideIdentification pep_id_empty;
pep_id_target.setHit(pep_hits_empty);

std::vector<PeptideIdentification> pep_ids = {pep_id_target, pep_id_decoy, pep_id_empty};

//construct features with peptideIdentifications
Feature feat_2_targets;
feat_2_targets.setPeptideIdentification(pep_ids);
feat_2_targets.setPeptideIdentification(pep_ids);
Feature feat_target;
feat_target.setPeptideIdentification(pep_ids);
Feature feat_empty;

std::vector<Feature> features = {feat_2_targets, feat_target, feat_empty};

//construct FeatureMap
FeatureMap(features) fmap;
fmap.setUnassignedPeptideIdentifications(pep_ids);


//construct MSSpectrum
MSSpectrum ms2_1;
ms2_1.setLevel(2);
MSSpectrum ms2_2;
ms2_2.setLevel(2);
MSSpectrum ms2_3;
ms2_3.setLevel(2);
MSSpectrum ms2_4;
ms2_4.setLevel(2);
MSSpectrum ms2_5;
ms2_5.setLevel(2);
MSSpectrum ms2_6;
ms2_6.setLevel(2);
MSSpectrum ms1;
ms1.setLevel(1);

std::vector<MSSpectrum> ms_spectra = {ms2_1, ms2_2, ms2_3, ms2_4, ms2_5, ms2_6, ms1};

//construct MSExperiment
MSExperiment exp;
exp.setSpectra(ms_spectra);

//tests compute function
Ms2IdentificationRate ms2ir;
ms2ir.compute(fmap, exp);
std::vector<IdentificationRateData> result;
result = ms2ir.getResult();

for (auto idrd : result)
{
  std::cout << "Number of Peptide IDs: " << idrd.num_peptide_identification << std::endl;
  std::cout << "Number of MS2 Spectra: " << idrd.num_ms2_spectra << std::endl;
  std::cout << "Rate: " << idrd.identification_rate << std::endl;
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
//END_TEST