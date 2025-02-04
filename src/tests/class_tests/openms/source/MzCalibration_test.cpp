// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow$
// $Authors: Juliane Schmachtenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
///////////////////////////
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/QC/MzCalibration.h>
///////////////////////////

START_TEST(MzCalibration, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
using namespace OpenMS;
using namespace std;

MzCalibration* ptr = nullptr;
MzCalibration* nullPointer = nullptr;
START_SECTION(MzCalibration())
ptr = new MzCalibration();
TEST_NOT_EQUAL(ptr, nullPointer);
END_SECTION

START_SECTION(~MzCalibration())
delete ptr;
END_SECTION

START_SECTION(QCBase::Status requirements() const override)
{
  MzCalibration mzCal;
  TEST_EQUAL(mzCal.requirements() == (QCBase::Status() | QCBase::Requires::POSTFDRFEAT), true);
}
END_SECTION

// PeakMap
PeakMap exp;
MSSpectrum spec;
Precursor pre;
std::vector<MSSpectrum> spectra;
pre.setMetaValue("mz_raw", 5);
spec.setMSLevel(2);
spec.setPrecursors({pre});
spec.setRT(0);
spec.setNativeID("XTandem::1");
spectra.push_back(spec);

pre.setMetaValue("mz_raw", 6);
spec.setPrecursors({pre});
spec.setRT(0.5);
spec.setNativeID("XTandem::2");
spectra.push_back(spec);

pre.setMetaValue("mz_raw", 7);
spec.setPrecursors({pre});
spec.setRT(1);
spec.setNativeID("XTandem::3");
spectra.push_back(spec);

exp.setSpectra(spectra);

MSExperiment exp_no_calibration = exp;

// adding processing info
DataProcessing p;
p.setProcessingActions({OpenMS::DataProcessing::CALIBRATION});
boost::shared_ptr<DataProcessing> p_(new DataProcessing(p));
for (Size i = 0; i < exp.size(); ++i)
{
  exp[i].getDataProcessing().push_back(p_);
}
for (Size i = 0; i < exp.getNrChromatograms(); ++i)
{
  exp.getChromatogram(i).getDataProcessing().push_back(p_);
}
QCBase::SpectraMap spectra_map(exp);

// FeatureMap
FeatureMap fmap_ref;
PeptideHit peptide_hit;
std::vector<PeptideHit> peptide_hits;
PeptideIdentification peptide_ID;
vector<PeptideIdentification> identifications;
vector<PeptideIdentification> unassignedIDs;
Feature feature1;
peptide_hit.setSequence(AASequence::fromString("AAAA"));
peptide_hit.setCharge(2);
peptide_hits.push_back(peptide_hit);
peptide_ID.setHits(peptide_hits);
peptide_ID.setRT(0);
peptide_ID.setMZ(5.5);
peptide_ID.setSpectrumReference( "XTandem::1");
identifications.push_back(peptide_ID);
peptide_hits.clear();
peptide_hit.setSequence(AASequence::fromString("WWWW"));
peptide_hit.setCharge(3);
peptide_hits.push_back(peptide_hit);
peptide_ID.setHits(peptide_hits);
peptide_ID.setRT(1);
peptide_ID.setSpectrumReference( "XTandem::3");
identifications.push_back(peptide_ID);
peptide_hits.clear();
feature1.setPeptideIdentifications(identifications);
fmap_ref.push_back(feature1);
// unassigned PeptideHits
peptide_hit.setSequence(AASequence::fromString("YYYY"));
peptide_hit.setCharge(2);
peptide_hits.push_back(peptide_hit);
peptide_ID.setHits(peptide_hits);
peptide_hits.clear();
peptide_ID.setRT(0.5);
peptide_ID.setSpectrumReference( "XTandem::2");
unassignedIDs.push_back(peptide_ID);
fmap_ref.setUnassignedPeptideIdentifications(unassignedIDs);
MzCalibration cal;
// tests compute function
START_SECTION(void compute(FeatureMap& features, const MSExperiment& exp, const QCBase::SpectraMap map_to_spectrum))
{
  FeatureMap fmap = fmap_ref;
  cal.compute(fmap, exp, spectra_map);

  // things that shouldn't change
  ABORT_IF(fmap.size() != 1);
  ABORT_IF(fmap[0].getPeptideIdentifications().size() != 2);
  ABORT_IF(fmap[0].getPeptideIdentifications()[0].getHits().size() != 1);
  ABORT_IF(fmap[0].getPeptideIdentifications()[1].getHits().size() != 1);
  ABORT_IF(fmap.getUnassignedPeptideIdentifications().size() != 1);
  ABORT_IF(fmap.getUnassignedPeptideIdentifications()[0].getHits().size() != 1);
  // things that should now be there
  for (const Feature& f : fmap)
  {
    for (const PeptideIdentification& pepID : f.getPeptideIdentifications())
    {
      ABORT_IF(!pepID.getHits()[0].metaValueExists("mz_raw"));
      ABORT_IF(!pepID.getHits()[0].metaValueExists("mz_ref"));
      ABORT_IF(!pepID.getHits()[0].metaValueExists("uncalibrated_mz_error_ppm"));
      ABORT_IF(!pepID.getHits()[0].metaValueExists("calibrated_mz_error_ppm"));
    }
  }
  for (const PeptideIdentification& upepID : fmap.getUnassignedPeptideIdentifications())
  {
    ABORT_IF(!upepID.getHits()[0].metaValueExists("mz_raw"));
    ABORT_IF(!upepID.getHits()[0].metaValueExists("mz_ref"));
  }

  // test with valid input
  TEST_REAL_SIMILAR(fmap[0].getPeptideIdentifications()[0].getHits()[0].getMetaValue("mz_raw"), 5);
  TEST_REAL_SIMILAR(fmap[0].getPeptideIdentifications()[1].getHits()[0].getMetaValue("mz_raw"), 7);
  // test unassigned
  TEST_REAL_SIMILAR(fmap.getUnassignedPeptideIdentifications()[0].getHits()[0].getMetaValue("mz_raw"), 6);

  // test refMZ
  double ref = AASequence::fromString("AAAA").getMonoWeight(OpenMS::Residue::Full, 2) / 2;
  TEST_REAL_SIMILAR(fmap[0].getPeptideIdentifications()[0].getHits()[0].getMetaValue("mz_ref"), ref);
  TEST_REAL_SIMILAR(fmap[0].getPeptideIdentifications()[1].getHits()[0].getMetaValue("mz_ref"), AASequence::fromString("WWWW").getMonoWeight(OpenMS::Residue::Full, 3) / 3);
  TEST_REAL_SIMILAR(fmap.getUnassignedPeptideIdentifications()[0].getHits()[0].getMetaValue("mz_ref"), AASequence::fromString("YYYY").getMonoWeight(OpenMS::Residue::Full, 2) / 2);

  // test  mz_error
  TEST_REAL_SIMILAR(fmap[0].getPeptideIdentifications()[0].getHits()[0].getMetaValue("uncalibrated_mz_error_ppm"), (5 - ref) / ref * 1000000);
  TEST_REAL_SIMILAR(fmap[0].getPeptideIdentifications()[0].getHits()[0].getMetaValue("calibrated_mz_error_ppm"), (5.5 - ref) / ref * 1000000);

  // test empty MSExperiment
  MSExperiment exp_empty {};
  QCBase::SpectraMap spectra_map_empty(exp_empty);
  fmap = fmap_ref; // reset FeatureMap
  cal.compute(fmap, exp_empty, spectra_map_empty);
  TEST_REAL_SIMILAR(fmap[0].getPeptideIdentifications()[0].getHits()[0].getMetaValue("uncalibrated_mz_error_ppm"), (5.5 - ref) / ref * 1000000);

  // test with exp where no calibration was performed
  fmap = fmap_ref;
  cal.compute(fmap, exp_no_calibration, spectra_map);
  TEST_REAL_SIMILAR(fmap[0].getPeptideIdentifications()[0].getHits()[0].getMetaValue("uncalibrated_mz_error_ppm"), (5.5 - ref) / ref * 1000000);

  // test empty FeatureMap
  FeatureMap fmap_empty {};
  cal.compute(fmap_empty, exp, spectra_map);
  TEST_EQUAL(fmap_empty.isMetaEmpty(), true);

  // test feature is empty
  Feature feature_empty {};
  fmap_empty.push_back(feature_empty);
  cal.compute(fmap_empty, exp, spectra_map);
  TEST_EQUAL(fmap_empty.isMetaEmpty(), true);

  // test empty PeptideIdentification
  fmap_empty.clear();
  PeptideIdentification peptide_ID_empty {};
  identifications.push_back(peptide_ID_empty);
  feature1.setPeptideIdentifications(identifications);
  fmap_empty.push_back(feature1);
  cal.compute(fmap_empty, exp, spectra_map);
  TEST_EQUAL(fmap_empty.isMetaEmpty(), true);

  // test empty hit
  fmap_empty.clear();
  peptide_ID.setHits(std::vector<PeptideHit> {});
  identifications.clear();
  identifications.push_back(peptide_ID);
  feature1.setPeptideIdentifications(identifications);
  fmap_empty.push_back(feature1);
  cal.compute(fmap_empty, exp, spectra_map);
  TEST_EQUAL(fmap_empty.isMetaEmpty(), true);

  // test wrong MS-Level exception
  fmap = fmap_ref; // reset FeatureMap
  fmap[0].getPeptideIdentifications()[0].setSpectrumReference( "XTandem::4");
  exp.getSpectra()[0].setNativeID("XTandem::4");
  exp.getSpectra()[0].setMSLevel(1);
  spectra_map.calculateMap(exp);
  TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, cal.compute(fmap, exp, spectra_map), "The matching spectrum of the mzML is not an MS2 Spectrum.");

  // test exception PepID without 'spectrum_reference'
  fmap = fmap_ref; // reset FeatureMap
  PeptideIdentification pep_no_spec_ref;
  PeptideHit dummy_hit;
  dummy_hit.setSequence(AASequence::fromString("MMMMM"));
  dummy_hit.setCharge(2);
  pep_no_spec_ref.setHits({dummy_hit});
  fmap[0].setPeptideIdentifications({pep_no_spec_ref});
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter, cal.compute(fmap, exp, spectra_map), "No spectrum reference annotated at peptide identification!");
}
END_SECTION

START_SECTION(const String& getName() const)
{
  TEST_EQUAL(cal.getName(), "MzCalibration");
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
