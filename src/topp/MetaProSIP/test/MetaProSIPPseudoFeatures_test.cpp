// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include "MetaProSIPPseudoFeatures.h"

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>

#include <optional>
#include <string>

using namespace OpenMS;

namespace
{
/// An identification of spectrum @p reference at @p rt, with a match of @p sequence (none if empty)
PeptideIdentification makeId(const std::string& reference, double rt, const std::string& sequence, int charge = 2)
{
  PeptideIdentification id;
  id.setIdentifier("run0");
  id.setSpectrumReference(reference);
  id.setRT(rt);
  id.setMZ(sequence.empty() ? 500.0 : AASequence::fromString(sequence).getMZ(charge == 0 ? 2 : charge));
  id.setScoreType("q-value");
  id.setHigherScoreBetter(false);
  if (!sequence.empty())
  {
    PeptideHit hit;
    hit.setSequence(AASequence::fromString(sequence));
    hit.setCharge(charge);
    hit.setScore(0.01);
    id.insertHit(hit);
  }
  return id;
}

/// A feature of PEPTIDEK (spectrum at 100 s), and unassigned identifications: SAMPLER (200 s), one without charge
/// (250 s) and one without matches (300 s)
FeatureMap makeFeatures()
{
  FeatureMap features;
  Feature feature;
  feature.setRT(100.0);
  feature.setMZ(AASequence::fromString("PEPTIDEK").getMZ(2));
  feature.setCharge(2);
  feature.setMetaValue("feature_type", "feature");
  features.push_back(feature);
  MapPeptideIdentifications ids;
  ProteinIdentification protein;
  protein.setIdentifier("run0");
  protein.setPrimaryMSRunPath({"run0.mzML"});
  ids.proteins = {protein};
  ids.features[{0}] = {makeId("scan=1", 100.0, "PEPTIDEK")};
  ids.unassigned = {makeId("scan=3", 200.0, "SAMPLER"), makeId("scan=4", 250.0, "SAMPLER", 0), makeId("scan=5", 300.0, "")};
  IdentificationDataConverter::importMapIdentifications(features, ids);
  return features;
}

/// MS2 spectra at 100, 150, 200, 250 and 300 s with a precursor (charge 3 at 150 s), and one at 350 s without
PeakMap makeSpectra()
{
  PeakMap ms2;
  for (const double rt : {100.0, 150.0, 200.0, 250.0, 300.0, 350.0})
  {
    MSSpectrum spectrum;
    spectrum.setMSLevel(2);
    spectrum.setRT(rt);
    if (rt != 350.0)
    {
      Precursor precursor;
      precursor.setMZ(400.0 + rt);
      precursor.setCharge(rt == 150.0 ? 3 : 2);
      spectrum.setPrecursors({precursor});
    }
    ms2.addSpectrum(spectrum);
  }
  ms2.updateRanges();
  return ms2;
}
} // namespace

START_TEST(MetaProSIPPseudoFeatures, "$Id$")

START_SECTION((static Size addUnassignedIdentifications(FeatureMap& features, const std::string& feature_type)))
{
  FeatureMap features = makeFeatures();
  TEST_EQUAL(features.getUnassignedIdentifications().size(), 3)
  // SAMPLER only: the identification without charge and the one without matches have no position to model
  TEST_EQUAL(MetaProSIPPseudoFeatures::addUnassignedIdentifications(features, "id"), 1)
  ABORT_IF(features.size() != 2)
  const Feature& pseudo = features[1];
  TEST_EQUAL(pseudo.getMetaValue("feature_type"), "id")
  TEST_REAL_SIMILAR(pseudo.getRT(), 200.0)
  TEST_REAL_SIMILAR(pseudo.getMZ(), AASequence::fromString("SAMPLER").getMZ(2))
  // it links the identification with its matches, which is no longer unassigned
  const auto linked = pseudo.getLinkedIdentifications(features.getIdentificationData());
  ABORT_IF(linked.size() != 1)
  TEST_EQUAL(linked[0].query->data_id, "scan=3")
  TEST_EQUAL(linked[0].matches.size(), 1)
  TEST_EQUAL(features.getUnassignedIdentifications().size(), 2)
  TEST_EQUAL(features[0].getLinkedIdentifications(features.getIdentificationData()).size(), 1)
}
END_SECTION

START_SECTION((static Size addUnidentifiedSpectra(FeatureMap& features, const PeakMap& ms2, const std::string& feature_type)))
{
  // the spectra of identifications with matches are taken: 100 s (feature), 200 s and 250 s (unassigned); the
  // identification without matches (300 s) claims none, and the spectrum at 350 s has no precursor
  FeatureMap features = makeFeatures();
  TEST_EQUAL(MetaProSIPPseudoFeatures::addUnidentifiedSpectra(features, makeSpectra(), "unidentified"), 2)
  ABORT_IF(features.size() != 3)
  TEST_REAL_SIMILAR(features[1].getRT(), 150.0)
  TEST_REAL_SIMILAR(features[1].getMZ(), 550.0)
  TEST_EQUAL(features[1].getCharge(), 3)
  TEST_EQUAL(features[1].getMetaValue("feature_type"), "unidentified")
  TEST_TRUE(features[1].getIDMatches().empty() && features[1].getIDQueries().empty())
  TEST_REAL_SIMILAR(features[2].getRT(), 300.0)
  TEST_EQUAL(features[2].getCharge(), 2)

  // after the pseudo features of unassigned identifications, the spectrum of SAMPLER (200 s) stays taken
  features = makeFeatures();
  MetaProSIPPseudoFeatures::addUnassignedIdentifications(features, "id");
  TEST_EQUAL(MetaProSIPPseudoFeatures::addUnidentifiedSpectra(features, makeSpectra(), "unidentified"), 2)
  TEST_EQUAL(features.size(), 4)

  // without identifications, every spectrum with a precursor
  FeatureMap empty;
  TEST_EQUAL(MetaProSIPPseudoFeatures::addUnidentifiedSpectra(empty, makeSpectra(), "unidentified"), 5)
}
END_SECTION

END_TEST
