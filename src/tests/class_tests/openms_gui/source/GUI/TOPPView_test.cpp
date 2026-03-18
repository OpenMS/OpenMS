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
#include <OpenMS/VISUAL/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/EnhancedTabBar.h>

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/VISUAL/GUIProgressLoggerImpl.h>

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

void TestTOPPView::waitForModalWidget(const int max_wait, const String& line)
{
  // test if accumulated scheduled time is less than max_wait
  int min_required_time = 0;
  foreach(ScheduleInfo i, modal_key_sequence_)
  {
    min_required_time+=i.delay;
  }
  if (min_required_time > max_wait)
  {
    QFAIL ( String("Test is bound to fail due to a time restriction in line " + line + ". Please rethink!").c_str());
  }


  QElapsedTimer t;
  t.start();
  while (!modal_key_sequence_.isEmpty () && max_wait > t.elapsed())
  {
    QTest::qWait(50);
  }

  if (!modal_key_sequence_.isEmpty ())
  {
    QWARN ( String("Modal dialogs timed out in line " + line + ". The following tests will most likely fail.").c_str());
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

// expands to a simple main() method that runs all the test functions
QTEST_MAIN(TestTOPPView) 

