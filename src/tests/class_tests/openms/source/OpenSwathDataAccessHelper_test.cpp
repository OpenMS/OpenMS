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
// $Maintainer: Hannes Roest $
// $Authors: Witold Wolski, Hannes Roest, $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/TraMLFile.h>

#include <boost/assign/std/vector.hpp>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(OpenSwathDataAccessHelper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

OpenSwathDataAccessHelper* ptr = nullptr;
OpenSwathDataAccessHelper* nullPointer = nullptr;

START_SECTION(OpenSwathDataAccessHelper())
{
  ptr = new OpenSwathDataAccessHelper();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~OpenSwathDataAccessHelper())
{
  delete ptr;
}
END_SECTION

START_SECTION(OpenSwathDataAccessHelper::convertToSpectrumPtr(sptr))
{

  MSSpectrum sptr,omsptr;
  Peak1D p1;
  p1.setIntensity(1.0f);
  p1.setMZ(2.0);

  Peak1D p2;
  p2.setIntensity(2.0f);
  p2.setMZ(10.0);

  Peak1D p3;
  p3.setIntensity(3.0f);
  p3.setMZ(30.0);

  TEST_STRING_EQUAL(sptr.getName(), "")
  sptr.setName("my_fancy_name");
  sptr.push_back(p1);
  sptr.push_back(p2);
  sptr.push_back(p3);
  OpenSwath::SpectrumPtr p = OpenSwathDataAccessHelper::convertToSpectrumPtr(sptr);
  OpenSwathDataAccessHelper::convertToOpenMSSpectrum(p,omsptr);
  TEST_REAL_SIMILAR(p->getMZArray()->data[0],2.0);
  TEST_REAL_SIMILAR(p->getMZArray()->data[1],10.0);
  TEST_REAL_SIMILAR(p->getMZArray()->data[2],30.0);


  TEST_REAL_SIMILAR(p->getIntensityArray()->data[0],1.0f);
  TEST_REAL_SIMILAR(p->getIntensityArray()->data[1],2.0f);
  TEST_REAL_SIMILAR(p->getIntensityArray()->data[2],3.0f);
}
END_SECTION

START_SECTION(OpenSwathDataAccessHelper::convertToOpenMSChromatogram(cptr, chromatogram))
{
  //void OpenSwathDataAccessHelper::convertToOpenMSChromatogram(OpenMS::MSChromatogram & chromatogram,
  //                                                          const OpenSwath::ChromatogramPtr cptr)
  OpenSwath::ChromatogramPtr cptr(new OpenSwath::Chromatogram());
  cptr->getTimeArray()->data.push_back(1.0);
  cptr->getTimeArray()->data.push_back(2.0);
  cptr->getTimeArray()->data.push_back(3.0);
  cptr->getTimeArray()->data.push_back(4.0);

  cptr->getIntensityArray()->data.push_back(4.0);
  cptr->getIntensityArray()->data.push_back(3.0);
  cptr->getIntensityArray()->data.push_back(2.0);
  cptr->getIntensityArray()->data.push_back(1.0);

  MSChromatogram chromatogram;
  OpenSwathDataAccessHelper::convertToOpenMSChromatogram(cptr, chromatogram);

  TEST_REAL_SIMILAR(chromatogram[0].getRT(),1.);
  TEST_REAL_SIMILAR(chromatogram[0].getIntensity(),4.);
  TEST_REAL_SIMILAR(chromatogram[1].getRT(),2.);
  TEST_REAL_SIMILAR(chromatogram[1].getIntensity(),3.);
  TEST_REAL_SIMILAR(chromatogram[2].getRT(),3.);
  TEST_REAL_SIMILAR(chromatogram[2].getIntensity(),2.);

}
END_SECTION

START_SECTION(convertToOpenMSSpectrum(spectrum,sptr))
{
  OpenSwath::SpectrumPtr cptr(new OpenSwath::Spectrum());
  cptr->getMZArray()->data.push_back(1.0);
  cptr->getMZArray()->data.push_back(2.0);
  cptr->getMZArray()->data.push_back(3.0);
  cptr->getMZArray()->data.push_back(4.0);

  cptr->getIntensityArray()->data.push_back(4.0);
  cptr->getIntensityArray()->data.push_back(3.0);
  cptr->getIntensityArray()->data.push_back(2.0);
  cptr->getIntensityArray()->data.push_back(1.0);

  MSSpectrum spectrum;
  OpenSwathDataAccessHelper::convertToOpenMSSpectrum(cptr, spectrum);

  TEST_REAL_SIMILAR(spectrum[0].getMZ(),1.);
  TEST_REAL_SIMILAR(spectrum[0].getIntensity(),4.);
  TEST_REAL_SIMILAR(spectrum[1].getMZ(),2.);
  TEST_REAL_SIMILAR(spectrum[1].getIntensity(),3.);
  TEST_REAL_SIMILAR(spectrum[2].getMZ(),3.);
  TEST_REAL_SIMILAR(spectrum[2].getIntensity(),2.);
}
END_SECTION

START_SECTION((void OpenSwathDataAccessHelper::convertTargetedExp(const OpenMS::TargetedExperiment & transition_exp_, OpenSwath::LightTargetedExperiment & transition_exp)))
{
  OpenMS::TargetedExperiment transition_exp_;
  OpenSwath::LightTargetedExperiment transition_exp;

  {
    TargetedExperiment::Peptide pep;
    OpenSwath::LightCompound comp;

    pep.setChargeState(8);
    pep.setDriftTime(0.6);

    // add a RT
    TargetedExperimentHelper::RetentionTime rt;
    rt.setRT(5.1);
    rt.retention_time_unit = TargetedExperimentHelper::RetentionTime::RTUnit::SECOND;
    rt.retention_time_type = TargetedExperimentHelper::RetentionTime::RTType::PREDICTED;
    pep.rts.push_back(rt);
    pep.id = "my_id";

    pep.setPeptideGroupLabel("group1");
    pep.protein_refs.push_back("pr1");
    pep.protein_refs.push_back("pr2");
    pep.sequence = "PEPTIDE";

    // add a modification
    TargetedExperimentHelper::Peptide::Modification m;
    m.mono_mass_delta = 123;
    m.location = 3;
    m.unimod_id = 5;
    pep.mods.push_back(m);

    transition_exp_.addPeptide(pep);

    ReactionMonitoringTransition transition;
    transition.setName("tr1");
    transition.setNativeID("tr1_nid");
    transition.setPeptideRef("my_id");
    transition.setLibraryIntensity(400.2);
    transition.setPrecursorMZ(501.2);
    transition.setProductMZ(301.2);
    auto p = transition.getProduct();
    p.setChargeState(4);
    transition.setProduct(p);
    transition.setDecoyTransitionType( ReactionMonitoringTransition::DecoyTransitionType::DECOY );
    transition.setDetectingTransition(false);
    transition.setQuantifyingTransition(true);
    transition.setIdentifyingTransition(true);

    transition_exp_.addTransition(transition);
  }

  {
    TargetedExperiment::Compound pep;
    OpenSwath::LightCompound comp;

    pep.setChargeState(8);
    pep.setDriftTime(0.6);

    // add a RT
    TargetedExperimentHelper::RetentionTime rt;
    rt.setRT(5.3);
    rt.retention_time_unit = TargetedExperimentHelper::RetentionTime::RTUnit::SECOND;
    rt.retention_time_type = TargetedExperimentHelper::RetentionTime::RTType::PREDICTED;
    pep.rts.push_back(rt);
    pep.id = "my_id";

    pep.theoretical_mass = 46.069;
    pep.molecular_formula = "C2H6O";
    pep.smiles_string = "CCO";
    pep.setMetaValue("CompoundName", "some_name");

    transition_exp_.addCompound(pep);
  }

  OpenSwathDataAccessHelper::convertTargetedExp(transition_exp_, transition_exp);

  TEST_EQUAL( transition_exp.getTransitions().size(), 1)
  auto tr = transition_exp.getTransitions()[0];

  TEST_EQUAL(tr.transition_name, "tr1_nid")
  TEST_EQUAL(tr.peptide_ref, "my_id")
  TEST_REAL_SIMILAR(tr.library_intensity, 400.2)
  TEST_REAL_SIMILAR(tr.precursor_mz, 501.2)
  TEST_REAL_SIMILAR(tr.product_mz, 301.2)
  TEST_EQUAL(tr.fragment_charge, 4)
  TEST_EQUAL(tr.decoy, true)
  TEST_EQUAL(tr.detecting_transition, false)
  TEST_EQUAL(tr.quantifying_transition, true)
  TEST_EQUAL(tr.identifying_transition, true)
}

END_SECTION

START_SECTION((static void convertTargetedCompound(const TargetedExperiment::Peptide& pep, OpenSwath::LightCompound& comp)))
{
  TargetedExperiment::Peptide pep;
  OpenSwath::LightCompound comp;

  pep.setChargeState(8);
  pep.setDriftTime(0.6);

  // add a RT
  TargetedExperimentHelper::RetentionTime rt;
  rt.setRT(5.1);
  rt.retention_time_unit = TargetedExperimentHelper::RetentionTime::RTUnit::SECOND;
  rt.retention_time_type = TargetedExperimentHelper::RetentionTime::RTType::PREDICTED;
  pep.rts.push_back(rt);
  pep.id = "my_id";

  pep.setPeptideGroupLabel("group1");
  pep.protein_refs.push_back("pr1");
  pep.protein_refs.push_back("pr2");
  pep.sequence = "PEPTIDE";

  // add a modification
  TargetedExperimentHelper::Peptide::Modification m;
  m.mono_mass_delta = 123;
  m.location = 3;
  m.unimod_id = 5;
  pep.mods.push_back(m);

  OpenSwathDataAccessHelper::convertTargetedCompound(pep, comp);

  TEST_REAL_SIMILAR(comp.getDriftTime(), 0.6);
  TEST_EQUAL(comp.getChargeState(), 8);
  TEST_REAL_SIMILAR(comp.rt, 5.1);
  TEST_EQUAL(comp.sequence, "PEPTIDE");
  TEST_EQUAL(comp.modifications.size(), 1);
  TEST_EQUAL(comp.protein_refs.size(), 2);
  TEST_EQUAL(comp.peptide_group_label, "group1");
  TEST_EQUAL(comp.id, "my_id");

  TEST_EQUAL(comp.modifications[0].location, 3);
  TEST_EQUAL(comp.modifications[0].unimod_id, 5);
}
END_SECTION

START_SECTION((static void convertTargetedCompound(const TargetedExperiment::Compound& compound, OpenSwath::LightCompound& comp)))
{
  TargetedExperiment::Compound pep;
  OpenSwath::LightCompound comp;

  pep.setChargeState(8);
  pep.setDriftTime(0.6);

  // add a RT
  TargetedExperimentHelper::RetentionTime rt;
  rt.setRT(5.3);
  rt.retention_time_unit = TargetedExperimentHelper::RetentionTime::RTUnit::SECOND;
  rt.retention_time_type = TargetedExperimentHelper::RetentionTime::RTType::PREDICTED;
  pep.rts.push_back(rt);
  pep.id = "my_id";

  pep.theoretical_mass = 46.069;
  pep.molecular_formula = "C2H6O";
  pep.smiles_string = "CCO";
  pep.setMetaValue("CompoundName", "some_name");

  OpenSwathDataAccessHelper::convertTargetedCompound(pep, comp);

  TEST_REAL_SIMILAR(comp.getDriftTime(), 0.6);
  TEST_EQUAL(comp.getChargeState(), 8);
  TEST_REAL_SIMILAR(comp.rt, 5.3);
  TEST_EQUAL(comp.sum_formula, "C2H6O");
  TEST_EQUAL(comp.compound_name, "some_name");
  TEST_EQUAL(comp.id, "my_id");
}
END_SECTION

START_SECTION((static void convertPeptideToAASequence(const OpenSwath::LightCompound & peptide, AASequence & aa_sequence)))
{
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



