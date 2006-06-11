// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
//
// --------------------------------------------------------------------------
// $Id: ViewProteinModificationScenario.h,v 1.2 2006/03/28 08:03:28 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:$
// --------------------------------------------------------------------------

#ifndef VIEWPROTEINMODIFICATIONSCENARIO_H
#define VIEWPROTEINMODIFICATIONSCENARIO_H

#include <qvariant.h>
#include <qdialog.h>

class QVBoxLayout;
class QHBoxLayout;
class QGridLayout;
class QSqlDatabase;
class QSqlCursor;
class QSqlForm;
class QPushButton;
class QDataBrowser;
class QSpinBox;
class QLabel;
class QLineEdit;
class QTextEdit;

namespace OpenMS
  {

  class ViewProteinModificationScenario : public QDialog
    {
      Q_OBJECT

    public:
      ViewProteinModificationScenario( QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
      ~ViewProteinModificationScenario();

      QPushButton* pushButton1;
      QDataBrowser* dataBrowser2;
      QSpinBox* QSpinBoxProtein_modification_scenario_ID;
      QSpinBox* QSpinBoxProtein_ID;
      QLabel* labelModification_combination_ID;
      QLabel* textLabel1;
      QSpinBox* QSpinBoxModification_combination_ID;
      QLabel* labelProtein_modification_scenario_ID;
      QPushButton* PushButtonFirst;
      QPushButton* PushButtonPrev;
      QPushButton* PushButtonNext;
      QPushButton* PushButtonLast;
      QLineEdit* lineEdit2;
      QLineEdit* QLineEditOverall_modifications;
      QLabel* labelOverall_modifications;
      QTextEdit* textEdit1;
      QLabel* labelPartial_modifications;
      QLabel* labelProtein_ID;
      QLabel* textLabel2;
      QLineEdit* lineEdit3;

    public slots:
      virtual void polish();

    protected:
      QHBoxLayout* layout5;

    protected slots:
      virtual void languageChange();

    };
}

#endif // VIEWPROTEINMODIFICATIONSCENARIO_H
