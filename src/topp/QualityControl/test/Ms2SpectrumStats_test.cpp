// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Juliane Schmachtenberg, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/NativeIdentificationTest.h>

///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include "Ms2SpectrumStats.h"

///////////////////////////

START_TEST(Ms2SpectrumStats, "$Id$");

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

Ms2SpectrumStats* ptr = nullptr;
Ms2SpectrumStats* nullPointer = nullptr;

START_SECTION(Ms2SpectrumStats())
{
  ptr = new Ms2SpectrumStats;
  TEST_NOT_EQUAL(ptr, nullPointer);
}
END_SECTION

START_SECTION(~Ms2SpectrumStats())
{
  delete ptr;
}
END_SECTION

Ms2SpectrumStats top;
START_SECTION(const std::string& getName() const override) {TEST_EQUAL(top.getName(), "Ms2SpectrumStats")} END_SECTION

  START_SECTION(QCBase::Status requirements() const override)
{
  TEST_EQUAL(top.requirements() == (QCBase::Status() | QCBase::Requires::RAWMZML | QCBase::Requires::POSTFDRFEAT), true);
}
END_SECTION

START_SECTION(compute(const MSExperiment& exp, FeatureMap& features, const QCBase::SpectraMap& map_to_spectrum))
{
  // Valid FeatureMap
  FeatureMap fmap;
  PeptideIdentification peptide_ID;
  PeptideIdentificationList identifications;
  PeptideIdentificationList unassignedIDs;
  Feature f1;
  peptide_ID.setSpectrumReference( "XTandem::0");
  identifications.push_back(peptide_ID);
  peptide_ID.setSpectrumReference( "XTandem::1");
  identifications.push_back(peptide_ID);
  MapPeptideIdentifications ids;
  ids.features[{0}] = identifications;
  identifications.clear();
  fmap.push_back(f1);
  peptide_ID.setSpectrumReference( "XTandem::10");
  identifications.push_back(peptide_ID);
  peptide_ID.setSpectrumReference( "XTandem::12");
  identifications.push_back(peptide_ID);
  ids.features[{1}] = identifications;
  fmap.push_back(f1);
  // unassigned PeptideHits
  peptide_ID.setSpectrumReference( "XTandem::1.5");
  unassignedIDs.push_back(peptide_ID);
  peptide_ID.setSpectrumReference( "XTandem::2.5");
  unassignedIDs.push_back(peptide_ID);
  ids.unassigned = unassignedIDs;

  // MSExperiment
  PeakMap exp;
  MSSpectrum spec;
  Peak1D p;
  Precursor pre;
  pre.setMZ(5.5);
  std::vector<MSSpectrum> spectra;
  spec.setPrecursors({pre});

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

  // not identified
  spec.setRT(20);
  spec.setNativeID("XTandem::20");
  p.setIntensity(5);
  spec.push_back(p);
  p.setIntensity(7);
  spec.push_back(p);
  spectra.push_back(spec);

  exp.setSpectra(spectra);

  QCBase::SpectraMap map_to_spectrum(exp);

  Ms2SpectrumStats top;
  std::vector<IdentificationData::Observation> new_unassigned_pep_ids;
  Internal::ClassTest::addSearchRun(ids);
  Internal::ClassTest::setIdentifications(fmap, ids);
  new_unassigned_pep_ids = top.compute(exp, fmap, map_to_spectrum);
  const auto computed = Internal::ClassTest::identificationsOf(fmap);

  // test features
  TEST_EQUAL(computed.at({0})[0].getMetaValue("ScanEventNumber"), 1);
  TEST_EQUAL(computed.at({0})[0].getMetaValue("identified"), 1);
  TEST_EQUAL(computed.at({0})[1].getMetaValue("ScanEventNumber"), 1);
  TEST_REAL_SIMILAR(computed.at({0})[1].getMetaValue("total_ion_count"), 6);
  TEST_REAL_SIMILAR(computed.at({0})[1].getMetaValue("base_peak_intensity"), 4);
  TEST_EQUAL(computed.at({1})[0].getMetaValue("ScanEventNumber"), 1);
  TEST_REAL_SIMILAR(computed.at({1})[1].getMetaValue("total_ion_count"), 10);
  TEST_REAL_SIMILAR(computed.at({1})[1].getMetaValue("base_peak_intensity"), 9);
  TEST_EQUAL(computed.at({1})[1].getMetaValue("ScanEventNumber"), 2);
  // test unassigned
  TEST_EQUAL(computed.unassigned[0].getMetaValue("ScanEventNumber"), 2);
  TEST_EQUAL(computed.unassigned[0].getMetaValue("identified"), 1);
  TEST_EQUAL(computed.unassigned[1].getMetaValue("ScanEventNumber"), 3);
  TEST_REAL_SIMILAR(*new_unassigned_pep_ids[0].rt, 20);
  TEST_EQUAL(new_unassigned_pep_ids[0].getMetaValue("ScanEventNumber"), 3);
  TEST_EQUAL(new_unassigned_pep_ids[0].getMetaValue("identified"), 0);
  TEST_REAL_SIMILAR(new_unassigned_pep_ids[0].getMetaValue("total_ion_count"), 12);
  TEST_REAL_SIMILAR(new_unassigned_pep_ids[0].getMetaValue("base_peak_intensity"), 7);
  TEST_REAL_SIMILAR(*new_unassigned_pep_ids[0].mz, 5.5);

  // empty FeatureMap
  FeatureMap fmap_empty {};
  new_unassigned_pep_ids = top.compute(exp, fmap_empty, map_to_spectrum);
  TEST_EQUAL(new_unassigned_pep_ids.size(), 7);

  // empty PeptideIdentifications
  fmap_empty.clear();
  fmap_empty.push_back(f1); // need some non-empty feature
  MapPeptideIdentifications empty_ids;
  empty_ids.features[{0}] = identifications; // those of the second feature above
  Internal::ClassTest::addSearchRun(empty_ids);
  Internal::ClassTest::setIdentifications(fmap_empty, empty_ids);
  new_unassigned_pep_ids = top.compute(exp, fmap_empty, map_to_spectrum);
  TEST_EQUAL(new_unassigned_pep_ids.size(), 5);
  // empty MSExperiment
  PeakMap exp_empty {};
  TEST_EXCEPTION(Exception::MissingInformation, top.compute(exp_empty, fmap, map_to_spectrum));

  // test exception PepID without 'spectrum_reference'
  PeptideIdentification pep_no_spec_ref;
  auto no_reference = Internal::ClassTest::identificationsOf(fmap);
  no_reference.features[{1}] = {pep_no_spec_ref};
  Internal::ClassTest::addSearchRun(no_reference);
  Internal::ClassTest::setIdentifications(fmap, no_reference);
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter, top.compute(exp, fmap, map_to_spectrum), "No spectrum reference annotated at peptide identification!");
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
