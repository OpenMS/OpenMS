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
#include <OpenMS/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/EnhancedTabBar.h>

using namespace OpenMS;


class OpenMS::TestTOPPView: public QObject
{
	Q_OBJECT

private slots:
	void testGui();
};


void TestTOPPView::testGui()
{
	TOPPViewBase tv;
	tv.show();
	// open file dialog
	QTest::keyClicks(&tv,"f", Qt::AltModifier);
	//QTest::qWait(1250);
	QWidget * dialog = QApplication::activeWindow ();
	// problem here: this will open a modal dialog which has
	// its own event loop and will not return...
	// Kitware's Paraview might have a solution: http://www.itk.org/Wiki/Testing_design
	if (0)
	{
		QTest::keyClicks(0,"e");
		// get file dialog
		// open peakpicker_tutorial_1.mzML
		//QTest::qWait(1250);
		std::cerr << "entering filename:\n";
		dialog = QApplication::activeWindow ();
		QTest::keyClicks(0,"peakpicker_tutorial_1.mzML");
		QTest::keyClick(0,Qt::Key_Return);
		QTest::qWait(1250);
		// get layer dialog
		// confirm
		QTest::keyClick(0,Qt::Key_Return);
		QTest::qWait(1250);
		QCOMPARE(tv.tab_bar_->tabText(tv.tab_bar_->currentIndex()), QString("peakpicker_tutorial_1.mzML"));
	}
	else
	{
		tv.addDataFile(File::getOpenMSDataPath() + "/examples/peakpicker_tutorial_1.mzML", false, false);
		QCOMPARE(tv.tab_bar_->tabText(tv.tab_bar_->currentIndex()), QString("peakpicker_tutorial_1.mzML"));
	}
}

QTEST_MAIN(TestTOPPView)		 // expands to a simple main() method that runs all the test functions
#include <../source/TEST/GUI/moc_TOPPView_test.cxx> // for Qt's introspection

