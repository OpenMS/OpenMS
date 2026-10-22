// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <TOPPView_test.h>

#include <QTimer>
#include <QElapsedTimer>
#include <QTemporaryDir>
#include <QSplashScreen>
#include <OpenMS/VISUAL/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/EnhancedTabBar.h>

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/VISUAL/GUIProgressLoggerImpl.h>
#include <OpenMS/VISUAL/LayerDataIdent.h>
#include <OpenMS/VISUAL/LayerDataPeak.h>
#include <OpenMS/VISUAL/Plot2DCanvas.h>

#include <fstream>

using namespace OpenMS;

void TestTOPPView::scheduleModalWidget_(const QString& key_sequence, const QString& title, const int delay)
{
  modal_key_sequence_.enqueue(ScheduleInfo(key_sequence, title, delay));
  std::cerr << "scheduled for window " << title.toStdString() << "\n";
  if (modal_key_sequence_.size()==1) // only schedule if this is the first entry
  {
    QTimer::singleShot(delay, this, SLOT(simulateClick_()));
  }
}

void TestTOPPView::waitForModalWidget(const int max_wait, const std::string& line)
{
  // test if accumulated scheduled time is less than max_wait
  int min_required_time = 0;
  foreach(ScheduleInfo i, modal_key_sequence_)
  {
    min_required_time+=i.delay;
  }
  if (min_required_time > max_wait)
  {
    QFAIL (std::string("Test is bound to fail due to a time restriction in line " + line + ". Please rethink!").c_str());
  }


  QElapsedTimer t;
  t.start();
  while (!modal_key_sequence_.isEmpty () && max_wait > t.elapsed())
  {
    QTest::qWait(50);
  }

  if (!modal_key_sequence_.isEmpty ())
  {
    QWARN (std::string("Modal dialogs timed out in line " + line + ". The following tests will most likely fail.").c_str());
    modal_key_sequence_.clear();
  }

}

void TestTOPPView::simulateClick_()
{
    if (!modal_key_sequence_.isEmpty ())
    {
      ScheduleInfo entry = modal_key_sequence_.head();
      std::cerr << "processing entry: '" << entry.keys.toStdString() << "' with dialog title '" << entry.title.toStdString() << "'\n";

      // search for a window
      QWidget * dialog = QApplication::activeModalWidget();
      if (!dialog) dialog = QApplication::activePopupWidget();
      if (!dialog) dialog = QApplication::activeWindow();

      if (!dialog || (dialog->windowTitle() != entry.title))
      {
        std::cerr << "item not found rescheduling...\n";
        QTimer::singleShot(100, this, SLOT(simulateClick_()));
        return;
      }
      QTest::keyClicks(dialog,entry.keys,Qt::NoModifier,20);
      QTest::keyClick(dialog,static_cast<Qt::Key>(Qt::Key_Return),Qt::NoModifier,20);
      QApplication::processEvents();

      // remove from queue
      modal_key_sequence_.dequeue();
    }

    if (!modal_key_sequence_.isEmpty ())
    {
      std::cerr << "Q not empty... rescheduling...\n";
      QTimer::singleShot(modal_key_sequence_.head().delay, this, SLOT(simulateClick_()));
    }
}

void TestTOPPView::testGui()
{
  // inject the GUIProgressLoggerImpl to be used by OpenMS lib via an extern variable
  make_gui_progress_logger = 
    []() -> ProgressLogger::ProgressLoggerImpl* { return new GUIProgressLoggerImpl(); };

  TOPPViewBase tv(TOPPViewBase::TOOL_SCAN::SKIP_SCAN);
  tv.show();
  QApplication::processEvents();

  QTest::qWait(1000);

#if 1 //def __APPLE__ // MAC OS does not support entering a filename via keyboard in the file-open menu
    tv.addDataFile(File::getOpenMSDataPath() + "/examples/peakpicker_tutorial_1.mzML", false, false);
    QCOMPARE(tv.tab_bar_.tabText(tv.tab_bar_.currentIndex()), QString("peakpicker_tutorial_1 (1D)"));
#else
    scheduleModalWidget_("peakpicker_tutorial_1.mzML", "Open file(s)",1000);                 // Open File dialog
    scheduleModalWidget_("", "Open data options for peakpicker_tutorial_1.mzML",1000); // layer data options dialog
    // open file dialog
    QTest::keyClicks(&tv,"f", Qt::AltModifier);
    QApplication::processEvents();
    // before we open the File-Open Dialog, we need to schedule the planned keyboard input
    // as this dialog is modal and won't return.
    // launch the modal widget
    QTest::keyClicks(0,"e");
    QApplication::processEvents();
    waitForModalWidget(15000, __LINE__);
#endif

  // compare the name of the opened tab
  QCOMPARE(tv.tab_bar_.tabText(tv.tab_bar_.currentIndex()), QString("peakpicker_tutorial_1 (1D)"));
}

void TestTOPPView::testIdentifications()
{
  using ID = IdentificationData;
  // identifications at RT 100/200 and m/z 500/600, and (with @p without_rt) one without RT
  const auto makeIdentifications = [](bool without_rt) {
    ID data;
    auto& run = data.addRun("search");
    ID::ScoreDefinition score;
    score.name = "q-value";
    score.higher_better = false;
    run.setPrimaryScore(run.addScore(score));
    const auto source = run.addSource({});
    const auto add = [&](std::optional<double> rt, double mz, const std::string& sequence) {
      ID::Observation observation;
      observation.rt = rt;
      observation.mz = mz;
      ID::MatchData match;
      match.representation = sequence;
      match.charge = 2;
      run.addMatch(run.addIdentification(source, observation), match, {0.01});
    };
    add(100.0, 500.0, "PEPTIDE");
    add(200.0, 600.0, "SAMPLER");
    if (without_rt) add(std::nullopt, 700.0, "PEPTIDER");
    return data;
  };
  const auto count = [](const ID& data) { return data.getUnlinked({}, {}).size(); };
  QTemporaryDir dir;
  const std::string path = dir.filePath("identifications.idXML").toStdString();
  FileHandler().storeIdentifications(path, makeIdentifications(true), {FileTypes::IDXML});
  const std::string path_with_rt = dir.filePath("identifications_with_rt.idXML").toStdString();
  FileHandler().storeIdentifications(path_with_rt, makeIdentifications(false), {FileTypes::IDXML});

  TOPPViewBase tv(TOPPViewBase::TOOL_SCAN::SKIP_SCAN);
  tv.show();
  QApplication::processEvents();

  // the identification layer (without the identification without RT)
  QCOMPARE(tv.addDataFile(path, false, false), TOPPViewBase::LOAD_RESULT::OK);
  auto* layer = dynamic_cast<LayerDataIdent*>(tv.getCurrentLayer());
  QVERIFY(layer != nullptr);
  QCOMPARE(count(layer->getIdentificationData()), Size(2));
  QCOMPARE(layer->getRange().getMinRT(), 100.0);
  QCOMPARE(layer->getRange().getMaxMZ(), 600.0);
  // merged identifications keep their runs
  auto* canvas = qobject_cast<Plot2DCanvas*>(tv.getActiveCanvas());
  QVERIFY(canvas != nullptr);
  canvas->mergeIntoLayer(0, makeIdentifications(true)); // the only layer of the window
  QCOMPARE(layer->getIdentificationData().getRuns().size(), Size(2));
  QCOMPARE(count(layer->getIdentificationData()), Size(5));
  QCOMPARE(layer->getRange().getMinRT(), 100.0);

  // a peak layer with MS2 spectra at the identifications
  auto run = std::make_shared<AnnotatedMSRun>();
  for (const auto& [rt, mz] : std::vector<std::pair<double, double>> {{100.0, 500.0}, {200.0, 600.0}})
  {
    MSSpectrum spectrum;
    spectrum.setRT(rt);
    spectrum.setMSLevel(2);
    Precursor precursor;
    precursor.setMZ(mz);
    spectrum.getPrecursors().push_back(precursor);
    spectrum.push_back(Peak1D(mz, 100.0f));
    run->getMSExperiment().addSpectrum(spectrum);
  }
  run->getMSExperiment().updateRanges();
  LayerDataBase::FeatureMapSharedPtrType features(new FeatureMap());
  LayerDataBase::ConsensusMapSharedPtrType consensus(new ConsensusMap());
  LayerDataBase::ODExperimentSharedPtrType on_disc(new OnDiscMSExperiment());
  ID none;
  tv.addData(features, consensus, none, run, on_disc, LayerDataBase::DT_PEAK, false, false, true, "", "spectra");
  QVERIFY(dynamic_cast<LayerDataPeak*>(tv.getCurrentLayer()) != nullptr);

  // annotated with the identifications of a file
  QSplashScreen splash;
  tv.loadFiles({"!", path_with_rt}, &splash);
  QCOMPARE(run->getSpectrumIdentifications(0).size(), Size(1));
  QCOMPARE(run->getSpectrumIdentifications(1).size(), Size(1));
  QCOMPARE(run->getSpectrumIdentifications(0)[0].matches[0]->representation, std::string("PEPTIDE"));

  // identifications without RT cannot be mapped, and a file that cannot be loaded leaves the spectra as they are (the
  // errors go to the log)
  tv.loadFiles({"!", path}, &splash);
  QCOMPARE(run->getSpectrumIdentifications(0).size(), Size(1));
  const std::string broken = dir.filePath("broken.idXML").toStdString();
  {
    std::ofstream out(broken);
    out << "not an idXML file";
  }
  tv.loadFiles({"!", broken}, &splash);
  QCOMPARE(run->getSpectrumIdentifications(0).size(), Size(1));
  QCOMPARE(run->getSpectrumIdentifications(0)[0].matches[0]->representation, std::string("PEPTIDE"));
}

// expands to a simple main() method that runs all the test functions
QTEST_MAIN(TestTOPPView) 

