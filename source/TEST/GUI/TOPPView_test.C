// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <QtGui>
#include <QtTest/QtTest>
#include <QTimer>
#include <QQueue>

#include <OpenMS/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/EnhancedTabBar.h>


namespace OpenMS
{


/**

	@todo write a proper GUI base class for the scheduler below (Chris)

	@todo document schedule for winXP
*/

class TestTOPPView: public QObject
{
	Q_OBJECT

	/**
		@brief Store information on timed keyboard input events
		
		Store the time offset, the key sequence and the expected window title.
		
	*/
	struct ScheduleInfo
	{
		ScheduleInfo():delay(0){};
		ScheduleInfo(QString p_keys, QString p_title, int p_delay)
			:keys(p_keys),
			 title(p_title),
			 delay(p_delay)
		{}
			
		QString keys; //< key sequence
		QString title;//< expected window title
		int delay;    //< delay in ms when event is fired off
	};
	
public slots:
	/**
		@brief Slot that tries to process the current event queue for modal dialogs until its empty.
		
		Slot that tries to process the current event queue for modal dialogs until its empty.
		The function will repeatedly invoke itself until the queue is empty, to allow other
		incoming events (e.g. loading a file) to be processed in between two scheduled dialogs.
	*/
	void simulateClick_();

private slots:
	void testGui();
	
	private:
		/**
			@brief Schedule a keyboard input using a QTimer signal to direct input to a modal window.
			
			Modal windows have their own event queue and once launched will halt the
			execution of the test script until closed. This implies one cannot simply direct keyboard
			input to them. To do that we pre-schedule the input in the main event loop
			using a timer. The keyboard input sequence @p key_sequence and a subsequent 'return' key press
			is then issued when the time out occurs, given that the current window has the correct
			@p title! Otherwise the event is rescheduled until the title is correct.
			The @p delay then the timer pops is relative to the last timed events successfull completion.
		*/
		void scheduleModalWidget_(const QString& key_sequence, const QString& title, const int delay=100);
		
		/**
			@brief Waits until the scheduled event queue is emtpy
			
			Waits until the scheduled event queue is emtpy.
		*/
		void waitForModalWidget(const int max_wait, const String& line);
		
		/// event queue for modal/popup dialogs
		QQueue<ScheduleInfo> modal_key_sequence_;
};



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
			QTest::keyClicks(0,entry.keys,Qt::NoModifier,20);
			QTest::keyClick(0,Qt::Key_Return,Qt::NoModifier,20);
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
	TOPPViewBase tv;
	tv.show();
	QApplication::processEvents();
	
	QTest::qWait(1000);

#if 1 //def __APPLE__ // MAC OS does not support entering a filename via keyboard in the file-open menu
		tv.addDataFile(File::getOpenMSDataPath() + "/examples/peakpicker_tutorial_1.mzML", false, false);
		QCOMPARE(tv.tab_bar_->tabText(tv.tab_bar_->currentIndex()), QString("peakpicker_tutorial_1.mzML"));
#else
	scheduleModalWidget_("peakpicker_tutorial_1.mzML", "Open file(s)",1000);								 // Open File dialog
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
	QCOMPARE(tv.tab_bar_->tabText(tv.tab_bar_->currentIndex()), QString("peakpicker_tutorial_1.mzML"));

}

} // Namespace OpenMS

using namespace OpenMS; // this is required to avoid linker errors on Windows

QTEST_MAIN(TestTOPPView)		 // expands to a simple main() method that runs all the test functions
#include <../source/TEST/GUI/moc_TOPPView_test.cxx> // for Qt's introspection

