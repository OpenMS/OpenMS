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
// $Maintainer:$
// --------------------------------------------------------------------------

#ifndef EDITPROTEIN_H
#define EDITPROTEIN_H

#include <qvariant.h>
#include <qdialog.h>

class QVBoxLayout;
class QHBoxLayout;
class QGridLayout;
class QSqlDatabase;
class QSqlCursor;
class QSqlForm;
class QDataBrowser;
class QLabel;
class QLineEdit;
class QSpinBox;
class QPushButton;

namespace OpenMS
  {

  class EditProtein : public QDialog
    {
      Q_OBJECT

    public:
      EditProtein( QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
      ~EditProtein();

      QDataBrowser* dataBrowser1;
      QLabel* labelIdentifier;
      QLineEdit* QLineEditAverage_mass;
      QLineEdit* QLineEditPdb_filename;
      QLineEdit* QLineEditIdentifier;
      QLabel* labelProtein_ID;
      QLabel* labelMono_mass;
      QLabel* labelFasta_filename;
      QSpinBox* QSpinBoxProtein_ID;
      QLineEdit* QLineEditSequence_oneletter;
      QLabel* labelSequence_oneletter;
      QLabel* labelPdb_filename;
      QLineEdit* QLineEditMono_mass;
      QLabel* labelNo_of_aminoacids;
      QLabel* labelAverage_mass;
      QLineEdit* QLineEditFasta_filename;
      QSpinBox* QSpinBoxNo_of_aminoacids;
      QPushButton* PushButtonFirst;
      QPushButton* PushButtonPrev;
      QPushButton* PushButtonNext;
      QPushButton* PushButtonLast;
      QPushButton* PushButtonInsert;
      QPushButton* PushButtonUpdate;
      QPushButton* PushButtonDelete;
      QPushButton* pushButton8;
      QPushButton* Help;

    public slots:
      virtual void polish();

      virtual void help();

    protected:
      QGridLayout* dataBrowser1Layout;
      QGridLayout* layout2;
      QHBoxLayout* layout3;
      QHBoxLayout* layout4;

    protected slots:
      virtual void languageChange();

    };
}
#endif // EDITPROTEIN_H
