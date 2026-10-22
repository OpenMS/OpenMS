// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <MetaDataBrowser_test.h>

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/VISUAL/MetaDataBrowser.h>
#include <OpenMS/VISUAL/VISUALIZER/IdentificationRunVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/MetaInfoVisualizer.h>

#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QStackedWidget>
#include <QtWidgets/QTreeWidget>

using namespace OpenMS;

namespace
{
  using ID = IdentificationData;

  /**
    A consensus map with identification data: run "search" (q-values) with database sequences P1 and P2,
    identification scan=1 (PEPTIDEA/P1 0.01, PEPTIDEB/P2 0.05) that feature 0 links, and the unassigned scan=2
    (PEPTIDEC/P2 0.02, with a numeric meta value)
  */
  ConsensusMap makeMap()
  {
    ConsensusMap map;
    auto& run = map.getIdentificationData().addRun("search");
    ID::ScoreDefinition score;
    score.name = "q-value";
    score.higher_better = false;
    run.setPrimaryScore(run.addScore(score));
    ID::Database database;
    database.path = "proteins.fasta";
    const auto db = run.addDatabase(database);
    run.setDatabaseSequences(std::vector<ID::DatabaseSequence> {{db, "P1"}, {db, "P2"}});
    const auto source = run.addSource({});
    ID::Observation observation;
    observation.data_id = "scan=1";
    const auto first = run.addIdentification(source, observation);
    observation.data_id = "scan=2";
    const auto second = run.addIdentification(source, observation);
    const auto add = [&](ID::QueryId query, const std::string& sequence, const std::string& accession, double value) {
      ID::MatchData match;
      match.representation = sequence;
      match.charge = 2;
      match.sequence_evidence.push_back({db, accession, std::nullopt, std::nullopt, 0, 0});
      if (sequence == "PEPTIDEC") match.setMetaValue("mass error", 1.5);
      return run.addMatch(query, match, {value});
    };
    // (the worse match first: the browser shows the best first)
    const auto b = add(first, "PEPTIDEB", "P2", 0.05);
    const auto a = add(first, "PEPTIDEA", "P1", 0.01);
    add(second, "PEPTIDEC", "P2", 0.02);
    ConsensusFeature feature;
    feature.setUniqueId(1);
    feature.addIDMatch({run.getUuid(), a});
    feature.addIDMatch({run.getUuid(), b});
    map.push_back(feature);
    return map;
  }

  /// The labels of the items of the tree of @p browser, in tree order
  QStringList labels(const MetaDataBrowser& browser)
  {
    QStringList result;
    auto* tree = browser.findChild<QTreeWidget*>();
    for (QTreeWidgetItemIterator it(tree); *it; ++it) result << (*it)->text(0);
    return result;
  }

  /// The visualizer of the item labeled @p label (its first child, if @p child)
  QWidget* visualizerOf(const MetaDataBrowser& browser, const QString& label, bool child = false)
  {
    auto* tree = browser.findChild<QTreeWidget*>();
    const auto items = tree->findItems(label, Qt::MatchExactly | Qt::MatchRecursive, 0);
    if (items.isEmpty()) return nullptr;
    auto* item = child ? items.first()->child(items.first()->childCount() - 1) : items.first();
    return browser.findChild<QStackedWidget*>()->widget(item->text(1).toInt());
  }
}

void TestMetaDataBrowser::testIdentificationTree()
{
  auto map = makeMap();
  MetaDataBrowser browser(true);
  browser.add(map);
  const auto shown = labels(browser);
  // the run with its database sequences, and the unassigned identification with its match
  QVERIFY(shown.contains("Run search"));
  QVERIFY(shown.contains("Prot P1"));
  QVERIFY(shown.contains("Prot P2"));
  QVERIFY(shown.contains("Identification scan=2"));
  QVERIFY(shown.contains("PEPTIDEC (0.02)"));
  // the identification that the feature links is not unassigned
  QVERIFY(!shown.contains("Identification scan=1"));
}

void TestMetaDataBrowser::testFeatureIdentifications()
{
  auto map = makeMap();
  MetaDataBrowser browser(false);
  browser.add(map, 0);
  const auto shown = labels(browser);
  // the linked identification, its matches best first (q-values: lower is better)
  const auto identification = shown.indexOf("Identification scan=1");
  QVERIFY(identification >= 0);
  QVERIFY(shown.indexOf("PEPTIDEA (0.01)") > identification);
  QVERIFY(shown.indexOf("PEPTIDEB (0.05)") > shown.indexOf("PEPTIDEA (0.01)"));
  QVERIFY(!shown.contains("Identification scan=2"));
}

void TestMetaDataBrowser::testEditAndSave()
{
  auto map = makeMap();
  const auto& data = map.getIdentificationData();
  const auto uuid = data.getRuns()[0].getUuid();
  const auto before = data;

  MetaDataBrowser browser(true);
  browser.add(map);

  // the settings of the run
  auto* run = qobject_cast<IdentificationRunVisualizer*>(visualizerOf(browser, "Run search"));
  QVERIFY(run != nullptr);
  run->findChild<QLineEdit*>("software")->setText("Comet");

  // a new meta value of a match (the fields for a new entry are the last three line edits of its metadata view)
  auto* metadata = qobject_cast<MetaInfoVisualizer*>(visualizerOf(browser, "PEPTIDEC (0.02)", true));
  QVERIFY(metadata != nullptr);
  const auto fields = metadata->findChildren<QLineEdit*>();
  QVERIFY(fields.size() >= 3);
  fields[fields.size() - 3]->setText("note");
  fields[fields.size() - 1]->setText("checked");
  for (auto* button : metadata->findChildren<QPushButton*>())
  {
    if (button->text() == "Add") button->click();
  }

  QVERIFY(QMetaObject::invokeMethod(&browser, "saveAll_"));

  const auto& edited = data.getRuns()[0];
  QCOMPARE(edited.getUuid(), uuid);
  QCOMPARE(edited.getSettings().software, std::string("Comet"));
  const auto& match = edited.getSources()[0].identifications[1].getMatches()[0];
  QCOMPARE(match.getMetaValue("note").toString(), std::string("checked"));
  // an unchanged value keeps its type
  QCOMPARE(match.getMetaValue("mass error").valueType(), DataValue::DOUBLE_VALUE);
  // records without edits are unchanged, and the feature still links its matches
  QVERIFY(edited.getSources()[0].identifications[0] == before.getRuns()[0].getSources()[0].identifications[0]);
  QCOMPARE(map[0].getLinkedIdentifications(data).size(), Size(1));
  QCOMPARE(map[0].getLinkedIdentifications(data)[0].matches.size(), Size(2));
}

// expands to a simple main() method that runs all the private slots
QTEST_MAIN(TestMetaDataBrowser)
