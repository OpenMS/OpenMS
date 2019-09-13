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
// $Authors: Juliane Schmachtenberg, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/QC/TopNoverRT.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/METADATA/DataProcessing.h>

///////////////////////////

START_TEST(TopNoverRT, "$Id$");

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

TopNoverRT* ptr = nullptr;
TopNoverRT* nullPointer = nullptr;

START_SECTION(TopNoverRT())
{
  ptr = new TopNoverRT;
  TEST_NOT_EQUAL(ptr, nullPointer);
}
END_SECTION

START_SECTION(~TopNoverRT())
{
  delete ptr;
}
END_SECTION

TopNoverRT top;
START_SECTION(const String& getName() const override)
{
  TEST_EQUAL(top.getName(), "TopNoverRT")
}
END_SECTION

START_SECTION(QCBase::Status requires() const override)
{
  TEST_EQUAL(top.requires() == (QCBase::Status() | QCBase::Requires::RAWMZML | QCBase::Requires::POSTFDRFEAT), true);
}
END_SECTION

START_SECTION(compute(const MSExperiment& exp, FeatureMap& features, const QCBase::SpectraMap& map_to_spectrum))
{
  //Valid FeatureMap
  FeatureMap fmap;
  PeptideIdentification peptide_ID;
  vector<PeptideIdentification> identifications;
  vector<PeptideIdentification> unassignedIDs;
  Feature f1;
  f1.setMetaValue("FWHM", 32.21);
  peptide_ID.setMetaValue("spectrum_reference","XTandem::0");
  identifications.push_back(peptide_ID);
  peptide_ID.setMetaValue("spectrum_reference","XTandem::1");
  identifications.push_back(peptide_ID);
  f1.setPeptideIdentifications(identifications);
  identifications.clear();
  fmap.push_back(f1);
  peptide_ID.setMetaValue("spectrum_reference","XTandem::10");
  identifications.push_back(peptide_ID);
  peptide_ID.setMetaValue("spectrum_reference","XTandem::12");
  identifications.push_back(peptide_ID);
  f1.setPeptideIdentifications(identifications);
  fmap.push_back(f1);
  //unassigned PeptideHits
  peptide_ID.setMetaValue("spectrum_reference","XTandem::1.5");
  unassignedIDs.push_back(peptide_ID);
  peptide_ID.setMetaValue("spectrum_reference","XTandem::2.5");
  unassignedIDs.push_back(peptide_ID);
  fmap.setUnassignedPeptideIdentifications(unassignedIDs);

  //MSExperiment
  PeakMap exp;
  MSSpectrum spec;
  Peak1D p;
  Precursor pre;
  pre.setMZ(5.5);
  std::vector< MSSpectrum> spectra;
  spec.setPrecursors({ pre });

  spec.setMSLevel(2);
  spec.setRT(0);
  spec.setNativeID("XTandem::0");
  p.setIntensity(2);
  spec.push_back(p);
  p.setIntensity(1);
  spec.push_back(p);
  spectra.push_back(spec);
  spec.clear(false);

  spec.setMSLevel(1);
  spec.setRT(0.5);
  spec.setNativeID("XTandem::0.5");
  spectra.push_back(spec);
  spec.clear(false);

  spec.setMSLevel(2);
  spec.setRT(1);
  spec.setNativeID("XTandem::1");
  p.setIntensity(4);
  spec.push_back(p);
  p.setIntensity(2);
  spec.push_back(p);
  spectra.push_back(spec);
  spec.clear(false);

  spec.setRT(1.5);
  spec.setNativeID("XTandem::1.5");
  spectra.push_back(spec);

  spec.setRT(2.5);
  spec.setNativeID("XTandem::2.5");
  spectra.push_back(spec);

  spec.setMSLevel(1);
  spec.setRT(9);
  spec.setNativeID("XTandem::9");
  spectra.push_back(spec);

  spec.setMSLevel(2);
  spec.setRT(10);
  spec.setNativeID("XTandem::10");
  p.setIntensity(3);
  spec.push_back(p);
  p.setIntensity(6);
  spec.push_back(p);
  spectra.push_back(spec);
  spec.clear(false);

  spec.setRT(12);
  spec.setNativeID("XTandem::12");
  p.setIntensity(1);
  spec.push_back(p);
  p.setIntensity(9);
  spec.push_back(p);
  spectra.push_back(spec);
  spec.clear(false);

  //not identified
  spec.setRT(20);
  spec.setNativeID("XTandem::20");
  p.setIntensity(5);
  spec.push_back(p);
  p.setIntensity(7);
  spec.push_back(p);
  spectra.push_back(spec);

  exp.setSpectra(spectra);

  QCBase::SpectraMap map_to_spectrum(exp);

  TopNoverRT top;
  vector<PeptideIdentification> new_unassigned_pep_ids;
  new_unassigned_pep_ids = top.compute(exp, fmap, map_to_spectrum);

  //test features
  TEST_EQUAL(fmap[0].getPeptideIdentifications()[0].getMetaValue("ScanEventNumber"), 1);
  TEST_EQUAL(fmap[0].getPeptideIdentifications()[0].getMetaValue("identified"), 1);
  TEST_EQUAL(fmap[0].getPeptideIdentifications()[1].getMetaValue("ScanEventNumber"), 1);
  TEST_REAL_SIMILAR(fmap[0].getPeptideIdentifications()[1].getMetaValue("total_ion_count"), 6);
  TEST_REAL_SIMILAR(fmap[0].getPeptideIdentifications()[1].getMetaValue("base_peak_intensity"), 4);
  TEST_EQUAL(fmap[1].getPeptideIdentifications()[0].getMetaValue("ScanEventNumber"), 1);
  TEST_REAL_SIMILAR(fmap[1].getPeptideIdentifications()[1].getMetaValue("total_ion_count"), 10);
  TEST_REAL_SIMILAR(fmap[1].getPeptideIdentifications()[1].getMetaValue("base_peak_intensity"), 9);
  TEST_EQUAL(fmap[1].getPeptideIdentifications()[1].getMetaValue("ScanEventNumber"), 2);
  //test unassigned
  TEST_EQUAL(fmap.getUnassignedPeptideIdentifications()[0].getMetaValue("ScanEventNumber"), 2);
  TEST_EQUAL(fmap.getUnassignedPeptideIdentifications()[0].getMetaValue("identified"), 1);
  TEST_EQUAL(fmap.getUnassignedPeptideIdentifications()[1].getMetaValue("ScanEventNumber"), 3);
  TEST_REAL_SIMILAR(new_unassigned_pep_ids[0].getRT(), 20);
  TEST_EQUAL(new_unassigned_pep_ids[0].getMetaValue("ScanEventNumber"), 3);
  TEST_EQUAL(new_unassigned_pep_ids[0].getMetaValue("identified"), 0);
  TEST_REAL_SIMILAR(new_unassigned_pep_ids[0].getMetaValue("total_ion_count"), 12);
  TEST_REAL_SIMILAR(new_unassigned_pep_ids[0].getMetaValue("base_peak_intensity"), 7);
  TEST_REAL_SIMILAR(new_unassigned_pep_ids[0].getMZ(), 5.5);

  // empty FeatureMap
  FeatureMap fmap_empty{};
  new_unassigned_pep_ids = top.compute(exp, fmap_empty, map_to_spectrum);
  TEST_EQUAL(new_unassigned_pep_ids.size(), 7);

  // empty feature
  fmap_empty.clear();
  Feature feature_empty{};
  fmap_empty.push_back(feature_empty);
  TEST_EXCEPTION(Exception::MissingInformation, top.compute(exp, fmap_empty, map_to_spectrum));

  // empty PeptideIdentifications
  fmap_empty.clear();
  fmap_empty.push_back(f1); // need a feature with FWHM
  feature_empty.setPeptideIdentifications( {} );
  fmap_empty.setUnassignedPeptideIdentifications( {} );
  new_unassigned_pep_ids = top.compute(exp, fmap_empty, map_to_spectrum);
  TEST_EQUAL(new_unassigned_pep_ids.size(), 5);
  // empty MSExperiment
  PeakMap exp_empty{};
  TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, top.compute(exp_empty, fmap, map_to_spectrum), "The mzml file / MSExperiment is empty.\n");

  // test exception PepID without 'spectrum_reference'
  PeptideIdentification pep_no_spec_ref;
  fmap[1].setPeptideIdentifications({pep_no_spec_ref});
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter, top.compute(exp, fmap, map_to_spectrum), "No spectrum reference annotated at peptide identification!");
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

