// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/VISUAL/LayerDataIdent.h>

///////////////////////////

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/METADATA/ID/IdentificationDataAdapter.h>
#include <OpenMS/VISUAL/LayerDataFeature.h>
#include <OpenMS/VISUAL/LayerDataPeak.h>
#include <OpenMS/VISUAL/VISITORS/LayerStatistics.h>
#include <OpenMS/VISUAL/VISITORS/LayerStoreData.h>

#include <numeric>

using namespace OpenMS;
using namespace std;

namespace
{
  using ID = IdentificationData;

  /**
    Identifications of run "search" (with database sequence P1): scan=1 at RT 100, m/z 500 (PEPTIDE, meta value
    "quality" 1.5) and scan=2 at RT 200, m/z 600 (SAMPLER, "quality" 2.5)
  */
  ID makeIdentifications()
  {
    ID data;
    auto& run = data.addRun("search");
    ID::ScoreDefinition score;
    score.name = "q-value";
    score.higher_better = false;
    run.setPrimaryScore(run.addScore(score));
    ID::Database database;
    database.path = "proteins.fasta";
    const auto db = run.addDatabase(database);
    run.setDatabaseSequences(std::vector<ID::DatabaseSequence> {{db, "P1"}});
    const auto source = run.addSource({});
    const auto add = [&](const std::string& data_id, double rt, double mz, const std::string& sequence, double quality) {
      ID::Observation observation;
      observation.data_id = data_id;
      observation.rt = rt;
      observation.mz = mz;
      observation.setMetaValue("quality", quality);
      const auto query = run.addIdentification(source, observation);
      ID::MatchData match;
      match.representation = sequence;
      match.charge = 2;
      match.sequence_evidence.push_back({db, "P1", std::nullopt, std::nullopt, 0, 0});
      run.addMatch(query, match, {0.01});
    };
    add("scan=1", 100.0, 500.0, "PEPTIDE", 1.5);
    add("scan=2", 200.0, 600.0, "SAMPLER", 2.5);
    return data;
  }

  /// A layer with the identifications of makeIdentifications()
  std::unique_ptr<LayerDataIdent> makeLayer()
  {
    auto layer = std::make_unique<LayerDataIdent>();
    layer->setIdentificationData(makeIdentifications());
    return layer;
  }

  /// The identifications that @p data stores as idXML (exported as peptide and protein identifications)
  IdentificationDataAdapter::LegacyResult storedAs(const LayerStoreData& data)
  {
    std::string path;
    NEW_TMP_FILE_EXT(path, ".idXML")
    data.saveToFile(path, ProgressLogger::NONE);
    ID loaded;
    FileHandler().loadIdentifications(path, loaded);
    return IdentificationDataAdapter::toLegacy(loaded);
  }
}

START_TEST(LayerDataIdent, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((void setIdentificationData(IdentificationData&& identifications)))
{
  auto layer = makeLayer();
  TEST_EQUAL(layer->type, LayerDataBase::DT_IDENT)
  TEST_EQUAL(layer->getIdentificationData().getRuns().size(), 1)
  TEST_EQUAL(layer->getIdentificationData().getUnlinked({}, {}).size(), 2)
}
END_SECTION

START_SECTION((RangeAllType getRange() const))
{
  const auto range = makeLayer()->getRange();
  TEST_REAL_SIMILAR(range.getMinRT(), 100.0)
  TEST_REAL_SIMILAR(range.getMaxRT(), 200.0)
  TEST_REAL_SIMILAR(range.getMinMZ(), 500.0)
  TEST_REAL_SIMILAR(range.getMaxMZ(), 600.0)
  TEST_EQUAL(LayerDataIdent().getRange().RangeRT::isEmpty(), true)
}
END_SECTION

START_SECTION((std::unique_ptr<LayerStoreData> storeVisibleData(const RangeAllType& visible_range, const DataFilters& layer_filters) const))
{
  auto layer = makeLayer();
  LayerDataBase::RangeAllType visible;
  visible.extendRT(50.0);
  visible.extendRT(150.0);
  visible.extendMZ(400.0);
  visible.extendMZ(550.0);
  // the visible identification, with its run and proteins
  const auto stored = storedAs(*layer->storeVisibleData(visible, DataFilters()));
  ABORT_IF(stored.peptides.size() != 1)
  TEST_EQUAL(stored.peptides[0].getHits()[0].getSequence().toString(), "PEPTIDE")
  ABORT_IF(stored.proteins.size() != 1)
  TEST_EQUAL(stored.proteins[0].getHits().size(), 1)
  // the layer is unchanged
  TEST_EQUAL(layer->getIdentificationData().getUnlinked({}, {}).size(), 2)
}
END_SECTION

START_SECTION((std::unique_ptr<LayerStoreData> storeFullData() const))
{
  const auto stored = storedAs(*makeLayer()->storeFullData());
  TEST_EQUAL(stored.peptides.size(), 2)
  ABORT_IF(stored.proteins.size() != 1)
  TEST_EQUAL(stored.proteins[0].getHits().size(), 1)
}
END_SECTION

START_SECTION((std::unique_ptr<LayerStatistics> getStats() const))
{
  auto layer = makeLayer();
  const auto stats = layer->getStats();
  const auto it = stats->getRangeStatistics().find({RangeStatsSource::METAINFO, "quality"});
  ABORT_IF(it == stats->getRangeStatistics().end())
  const auto& quality = std::get<RangeStatsDouble>(it->second);
  TEST_EQUAL(quality.getCount(), 2)
  TEST_REAL_SIMILAR(quality.getMin(), 1.5)
  TEST_REAL_SIMILAR(quality.getMax(), 2.5)
  const auto histogram = stats->getDistribution({RangeStatsSource::METAINFO, "quality"}, 2);
  TEST_REAL_SIMILAR(std::accumulate(histogram.begin(), histogram.end(), 0.0), 2.0)
}
END_SECTION

START_SECTION((bool annotate(const IdentificationData& identifications)))
{
  // an identification layer is not annotated with identifications
  TEST_EQUAL(LayerDataIdent().annotate(makeIdentifications()), false)

  // spectra: an MS2 spectrum at each identification, and one without
  LayerDataPeak peaks;
  auto run = std::make_shared<AnnotatedMSRun>();
  for (const auto& [rt, mz] : std::vector<std::pair<double, double>> {{100.0, 500.0}, {200.0, 600.0}, {300.0, 700.0}})
  {
    MSSpectrum spectrum;
    spectrum.setRT(rt);
    spectrum.setMSLevel(2);
    Precursor precursor;
    precursor.setMZ(mz);
    spectrum.getPrecursors().push_back(precursor);
    run->getMSExperiment().addSpectrum(spectrum);
  }
  peaks.setPeakData(run);
  TEST_EQUAL(peaks.annotate(makeIdentifications()), true)
  TEST_EQUAL(run->getSpectrumIdentifications(0).size(), 1)
  TEST_EQUAL(run->getSpectrumIdentifications(1).size(), 1)
  TEST_EQUAL(run->getSpectrumIdentifications(2).size(), 0)
  // the identifications replace those of the spectra
  TEST_EQUAL(peaks.annotate(makeIdentifications()), true)
  TEST_EQUAL(run->getIdentificationData().getRuns().size(), 1)
  TEST_EQUAL(run->getSpectrumIdentifications(0).size(), 1)
  // identifications without RT cannot be mapped
  auto without_rt = makeIdentifications();
  without_rt.eraseIdentifications([](const ID::Run&, const ID::Identification& query) { return query.data_id == "scan=2"; });
  {
    ID::Observation observation;
    observation.mz = 700.0;
    auto& search = without_rt.getRun("search");
    search.addIdentification(search.getSources()[0].id, observation);
  }
  TEST_EXCEPTION(Exception::MissingInformation, peaks.annotate(without_rt))
  TEST_EQUAL(run->getSpectrumIdentifications(0).size(), 1)

  // features: the identification at a feature (of its charge) is linked by it
  LayerDataFeature features;
  Feature feature;
  feature.setRT(200.0);
  feature.setMZ(600.0);
  feature.setCharge(2);
  features.getFeatureMap()->push_back(feature);
  TEST_EQUAL(features.annotate(makeIdentifications()), true)
  const auto& map = *features.getFeatureMap();
  const auto linked = map[0].getLinkedIdentifications(map.getIdentificationData());
  ABORT_IF(linked.size() != 1)
  TEST_EQUAL(linked[0].query->data_id, "scan=2")
  TEST_EQUAL(map.getUnassignedIdentifications().size(), 1)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
