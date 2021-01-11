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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <TOPPView_test.h>

#include <QTimer>

#include <OpenMS/VISUAL/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/EnhancedTabBar.h>

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/Factory.h>
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


  QTime t;
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
  // register a GUI logger
  Factory<ProgressLogger::ProgressLoggerImpl>::registerProduct(GUIProgressLoggerImpl::getProductName(), &GUIProgressLoggerImpl::create);

  TOPPViewBase tv;
  tv.show();
  QApplication::processEvents();

  QTest::qWait(1000);

#if 1 //def __APPLE__ // MAC OS does not support entering a filename via keyboard in the file-open menu
    tv.addDataFile(File::getOpenMSDataPath() + "/examples/peakpicker_tutorial_1.mzML", false, false);
    QCOMPARE(tv.tab_bar_.tabText(tv.tab_bar_.currentIndex()), QString("peakpicker_tutorial_1"));
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
  QCOMPARE(tv.tab_bar_.tabText(tv.tab_bar_.currentIndex()), QString("peakpicker_tutorial_1"));
}

// expands to a simple main() method that runs all the test functions
QTEST_MAIN(TestTOPPView) 

