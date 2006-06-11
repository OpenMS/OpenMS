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
// $Id: EditModification.h,v 1.2 2006/03/28 08:03:28 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:$
// --------------------------------------------------------------------------

#ifndef EDITMODIFICATION_H
#define EDITMODIFICATION_H

#include <qvariant.h>
#include <qdialog.h>

class QVBoxLayout;
class QHBoxLayout;
class QGridLayout;
class QSqlDatabase;
class QSqlCursor;
class QSqlForm;
class QDataBrowser;
class QLineEdit;
class QLabel;
class QTextEdit;
class QPushButton;


namespace OpenMS
  {

  class EditModification : public QDialog
    {
      Q_OBJECT

    public:
      EditModification( QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
      ~EditModification();

      QDataBrowser* dataBrowser1;
      QLineEdit* QLineEditMinus_average_mass;
      QLineEdit* QLineEditPlus_average_mass;
      QLineEdit* QLineEditPlus_formula;
      QLineEdit* QLineEditModification_name;
      QLineEdit* QLineEditModification_ID;
      QLabel* labelPlus_mono_mass;
      QLabel* labelMinus_average_mass;
      QLineEdit* QLineEditMinus_formula;
      QLabel* labelPlus_formula;
      QLineEdit* QLineEditModification_sites;
      QLineEdit* QLineEditMinus_mono_mass;
      QLabel* labelMinus_mono_mass;
      QLabel* labelModification_sites;
      QLineEdit* QLineEditPlus_mono_mass;
      QLabel* labelMinus_formula;
      QLabel* labelPlus_average_mass;
      QLabel* labelModification_name;
      QLabel* labelModification_ID;
      QLabel* textLabel1;
      QTextEdit* textEdit1;
      QPushButton* PushButtonFirst;
      QPushButton* PushButtonPrev;
      QPushButton* PushButtonNext;
      QPushButton* PushButtonLast;
      QPushButton* PushButtonInsert;
      QPushButton* PushButtonUpdate;
      QPushButton* PushButtonDelete;
      QPushButton* Help;
      QPushButton* Done;

    public slots:
      virtual void polish();

      virtual void help();

    protected:
      QVBoxLayout* EditModificationLayout;
      QVBoxLayout* dataBrowser1Layout;
      QGridLayout* layout1;
      QVBoxLayout* layout5;
      QHBoxLayout* layout2;
      QHBoxLayout* layout3;
      QHBoxLayout* layout4;

    protected slots:
      virtual void languageChange();

    };
}
#endif // EDITMODIFICATION_H
