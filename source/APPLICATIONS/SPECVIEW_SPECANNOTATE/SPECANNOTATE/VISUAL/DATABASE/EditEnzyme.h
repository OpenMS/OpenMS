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

#ifndef EDITENZYME_H
#define EDITENZYME_H

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
class QPushButton;

namespace OpenMS
  {

  class EditEnzyme : public QDialog
    {
      Q_OBJECT

    public:
      EditEnzyme( QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
      ~EditEnzyme();

      QDataBrowser* dataBrowser1;
      QLabel* labelTerminality;
      QLineEdit* QLineEditTerminality;
      QLineEdit* QLineEditEnzyme_name;
      QLabel* labelCleavage_sites;
      QLineEdit* QLineEditCleavage_sites;
      QLabel* labelEnzyme_name;
      QPushButton* PushButtonFirst;
      QPushButton* PushButtonPrev;
      QPushButton* PushButtonNext;
      QPushButton* PushButtonLast;
      QPushButton* PushButtonInsert;
      QPushButton* PushButtonUpdate;
      QPushButton* PushButtonDelete;
      QPushButton* Done;
      QPushButton* Help;

    public slots:
      virtual void polish();

      virtual void help();

    protected:
      QGridLayout* dataBrowser1Layout;
      QGridLayout* layout1;
      QHBoxLayout* layout2;
      QHBoxLayout* layout3;

    protected slots:
      virtual void languageChange();

    };
}

#endif // EDITENZYME_H
