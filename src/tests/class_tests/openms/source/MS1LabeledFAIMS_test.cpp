// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include "MS1LabeledFAIMS.h"

#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Constants.h>

using namespace OpenMS;

namespace
{
PeptideIdentification makeId(const std::string& reference, std::optional<double> cv = std::nullopt)
{
  PeptideIdentification id;
  id.setIdentifier("run0");
  id.setSpectrumReference(reference);
  id.setMZ(500.0);
  id.setRT(100.0);
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDEK"));
  hit.setCharge(2);
  id.insertHit(hit);
  if (cv) { id.setMetaValue(Constants::UserParam::FAIMS_CV, *cv); }
  return id;
}

MSExperiment makeSpectra()
{
  MSExperiment spectra;
  for (int i = 0; i < 4; ++i)
  {
    MSSpectrum spectrum;
    spectrum.setNativeID("controllerType=0 controllerNumber=1 scan=" + std::to_string(i + 1));
    spectrum.setRT(100.0 + i * 0.01);
    spectrum.setMSLevel(i % 2 == 0 ? 1 : 2);
    // The first MS2 inherits its preceding MS1's CV, as in IMDataConverter.
    if (i != 1)
    {
      spectrum.setDriftTimeUnit(DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
      spectrum.setDriftTime(i < 2 ? -45.0 : -65.0);
    }
    spectra.addSpectrum(spectrum);
  }
  return spectra;
}

ConsensusMap makeMap(Size run)
{
  ConsensusMap map;
  map.setExperimentType("labeled_MS1");
  for (Size channel = 0; channel < 2; ++channel)
  {
    auto& header = map.getColumnHeaders()[channel];
    header.filename = "run" + std::to_string(run) + ".mzML";
    header.label = channel == 0 ? "no_label" : "Lys8";
    header.setMetaValue("channel_id", channel);
  }
  ProteinIdentification protein;
  protein.setIdentifier("run" + std::to_string(run));
  protein.setPrimaryMSRunPath({"run" + std::to_string(run) + ".mzML"});
  map.setProteinIdentifications({protein});
  for (const double cv : {-45.0, -65.0})
  {
    ConsensusFeature feature;
    feature.setUniqueId();
    feature.setMetaValue(Constants::UserParam::FAIMS_CV, cv);
    feature.setRT(100.0);
    feature.setMZ(cv == -45.0 ? 500.0 : 500.0001);
    feature.setCharge(2);
    for (Size channel = 0; channel < 2; ++channel)
    {
      FeatureHandle handle;
      handle.setUniqueId();
      handle.setMapIndex(channel);
      handle.setMZ(feature.getMZ() + channel * 4.0071);
      handle.setRT(100.0);
      handle.setCharge(2);
      handle.setIntensity(-cv * (run + 1));
      feature.insert(handle);
    }
    feature.setIntensity(-cv * (run + 1));
    map.push_back(feature);
  }
  map.updateRanges();
  return map;
}
} // namespace

START_TEST(MS1LabeledFAIMS, "$Id$")

START_SECTION((static void annotateCompensationVoltages(const MSExperiment&, PeptideIdentificationList&)))
{
  const auto spectra = makeSpectra();
  PeptideIdentificationList ids {makeId(spectra[1].getNativeID()), makeId("scan=4"), makeId("index=3")};
  MS1LabeledFAIMS::annotateCompensationVoltages(spectra, ids);
  TEST_REAL_SIMILAR((double)ids[0].getMetaValue(Constants::UserParam::FAIMS_CV), -45.0)
  TEST_REAL_SIMILAR((double)ids[1].getMetaValue(Constants::UserParam::FAIMS_CV), -65.0)
  TEST_REAL_SIMILAR((double)ids[2].getMetaValue(Constants::UserParam::FAIMS_CV), -65.0)

  ids = {makeId("", -65.001)};
  ids[0].setRT(spectra[3].getRT());
  MS1LabeledFAIMS::annotateCompensationVoltages(spectra, ids);
  TEST_REAL_SIMILAR((double)ids[0].getMetaValue(Constants::UserParam::FAIMS_CV), -65.0)
  TEST_EQUAL(ids[0].getSpectrumReference(), spectra[3].getNativeID())
  ids = {makeId("")};
  TEST_EXCEPTION(Exception::MissingInformation, MS1LabeledFAIMS::annotateCompensationVoltages(spectra, ids))
  ids = {makeId("scan=999")};
  TEST_EXCEPTION(Exception::MissingInformation, MS1LabeledFAIMS::annotateCompensationVoltages(spectra, ids))
  ids = {makeId("scan=4", -45.0)};
  TEST_EXCEPTION(Exception::InvalidParameter, MS1LabeledFAIMS::annotateCompensationVoltages(spectra, ids))
  ids = {makeId("", -85.0)};
  TEST_EXCEPTION(Exception::InvalidParameter, MS1LabeledFAIMS::annotateCompensationVoltages(spectra, ids))

  MSExperiment single_cv = spectra;
  single_cv.resize(2);
  ids = {makeId("")};
  ids[0].setRT(spectra[1].getRT());
  MS1LabeledFAIMS::annotateCompensationVoltages(single_cv, ids);
  TEST_REAL_SIMILAR((double)ids[0].getMetaValue(Constants::UserParam::FAIMS_CV), -45.0)
  TEST_EQUAL(ids[0].getSpectrumReference(), spectra[1].getNativeID())

  auto coincident = spectra;
  coincident[3].setRT(coincident[1].getRT());
  ids = {makeId("", -45.0), makeId("", -65.0)};
  for (auto& id : ids)
  {
    id.setRT(coincident[1].getRT());
  }
  MS1LabeledFAIMS::annotateCompensationVoltages(coincident, ids);
  TEST_EQUAL(ids[0].getSpectrumReference(), coincident[1].getNativeID())
  TEST_EQUAL(ids[1].getSpectrumReference(), coincident[3].getNativeID())
  ids = {makeId("", -65.0)};
  ids[0].setRT(200.0);
  TEST_EXCEPTION(Exception::MissingInformation, MS1LabeledFAIMS::annotateCompensationVoltages(coincident, ids))

  ids = {makeId("")};
  MS1LabeledFAIMS::annotateCompensationVoltages(MSExperiment(), ids);
  TEST_FALSE(ids[0].metaValueExists(Constants::UserParam::FAIMS_CV))
}
END_SECTION

START_SECTION((static void annotate(IDMapper&, ConsensusMap&, const PeptideIdentificationList&, const std::vector<ProteinIdentification>&)))
{
  auto map = makeMap(0);
  const auto proteins = map.getProteinIdentifications();
  map.getProteinIdentifications().clear();
  IDMapper mapper;
  // The wrong CV is an exact m/z match. It must never receive the identification.
  PeptideIdentificationList ids {makeId("scan=4", -65.0), makeId("scan=9", -85.0)};
  MS1LabeledFAIMS::annotate(mapper, map, ids, proteins);
  TEST_EQUAL(map.size(), 2)
  TEST_EQUAL(map[0].getPeptideIdentifications().size(), 0)
  TEST_EQUAL(map[1].getPeptideIdentifications().size(), 1)
  TEST_EQUAL((int)map[1].getPeptideIdentifications()[0].getMetaValue("map_index"), 0)
  TEST_EQUAL(map.getUnassignedPeptideIdentifications().size(), 1)
  TEST_EQUAL(map.getUnassignedPeptideIdentifications()[0].getSpectrumReference(), "scan=9")
  TEST_EQUAL(map.getProteinIdentifications().size(), 1)
  TEST_EQUAL(map.getColumnHeaders().at(0).filename, "run0.mzML")
}
END_SECTION

START_SECTION((static void group(FeatureGroupingAlgorithmQT&, const std::vector<ConsensusMap>&, ConsensusMap&)))
{
  std::vector<ConsensusMap> maps {makeMap(0), makeMap(1)};
  // Deliberately make the wrong CV the closest m/z match in the second run.
  for (auto& feature : maps[1])
  {
    const double shift = (double)feature.getMetaValue(Constants::UserParam::FAIMS_CV) == -45.0 ? 0.0001 : -0.0001;
    feature.setMZ(feature.getMZ() + shift);
    ConsensusFeature::HandleSetType handles;
    for (auto handle : feature.getFeatures())
    {
      handle.setMZ(handle.getMZ() + shift);
      handles.insert(handle);
    }
    feature.setFeatures(std::move(handles));
  }
  for (auto& feature : maps[0])
  {
    auto id = makeId("scan=" + std::to_string(feature.getUniqueId()), (double)feature.getMetaValue(Constants::UserParam::FAIMS_CV));
    id.setMetaValue("map_index", 0);
    id.setMetaValue("old_map_index", 0);
    feature.getPeptideIdentifications().push_back(id);
  }
  // Run 1 has no IDs: matching must transfer each ID only within its CV.
  FeatureGroupingAlgorithmQT linker;
  Param linking = linker.getParameters();
  linking.setValue("use_identifications", "true");
  linker.setParameters(linking);
  ConsensusMap result;
  MS1LabeledFAIMS::group(linker, maps, result);
  TEST_EQUAL(result.size(), 2)
  TEST_EQUAL(result.getColumnHeaders().size(), 4)
  TEST_EQUAL(result.getColumnHeaders().at(0).filename, "run0.mzML")
  TEST_EQUAL(result.getColumnHeaders().at(2).filename, "run1.mzML")
  TEST_EQUAL(result.getProteinIdentifications().size(), 2)
  for (const auto& feature : result)
  {
    const double cv = feature.getMetaValue(Constants::UserParam::FAIMS_CV);
    TEST_EQUAL(feature.size(), 4)
    TEST_EQUAL(feature.getPeptideIdentifications().size(), 1)
    TEST_REAL_SIMILAR((double)feature.getPeptideIdentifications()[0].getMetaValue(Constants::UserParam::FAIMS_CV), cv)
    for (const auto& handle : feature.getFeatures())
    {
      TEST_REAL_SIMILAR(handle.getIntensity(), -cv * (handle.getMapIndex() / 2 + 1))
    }
  }
  // Missing CV in one run must not remove its columns or allow a cross-CV match.
  maps[1].resize(1);
  MS1LabeledFAIMS::group(linker, maps, result);
  TEST_EQUAL(result.size(), 2)
  TEST_EQUAL(result.getColumnHeaders().size(), 4)
  for (const auto& feature : result)
  {
    TEST_EQUAL(feature.size(), (double)feature.getMetaValue(Constants::UserParam::FAIMS_CV) == -45.0 ? 4 : 2)
  }
  // A CV with IDs but no detected features retains its IDs once, without adding sample columns.
  maps[0].getUnassignedPeptideIdentifications().push_back(makeId("scan=999", -85.0));
  MS1LabeledFAIMS::group(linker, maps, result);
  TEST_EQUAL(result.size(), 2)
  TEST_EQUAL(result.getColumnHeaders().size(), 4)
  TEST_EQUAL(result.getUnassignedPeptideIdentifications().size(), 1)
  TEST_EQUAL(result.getUnassignedPeptideIdentifications()[0].getSpectrumReference(), "scan=999")
  TEST_EQUAL(result.getProteinIdentifications().size(), 2)
}
END_SECTION

END_TEST
