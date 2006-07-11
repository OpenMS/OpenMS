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

#ifndef EDITMODIFICATIONXML_H
#define EDITMODIFICATIONXML_H

#include <qvariant.h>
#include <qdialog.h>
#include <string>
#include <OpenMS/FORMAT/Param.h>

class QVBoxLayout;
class QHBoxLayout;
class QGridLayout;
class QSpacerItem;
class QLineEdit;
class QLabel;
class QTextEdit;
class QPushButton;


namespace OpenMS
  {


  class EditModificationXML : public QDialog
    {
      Q_OBJECT

    public:
      EditModificationXML( QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
      ~EditModificationXML();

      QLineEdit* lineEdit9;
      QLineEdit* lineEdit5;
      QLineEdit* lineEdit4;
      QLineEdit* lineEdit6;
      QLabel* textLabel2;
      QLabel* textLabel9;
      QLabel* textLabel5;
      QLineEdit* lineEdit8;
      QLabel* textLabel7;
      QLineEdit* lineEdit2;
      QLabel* textLabel8;
      QLineEdit* lineEdit3;
      QLabel* textLabel4;
      QLabel* textLabel3;
      QLabel* textLabel6;
      QLabel* textLabel1;
      QLineEdit* lineEdit1;
      QLineEdit* lineEdit7;
      QLabel* textLabel10;
      QTextEdit* textEdit1;
      QPushButton* pushButton1;
      QPushButton* pushButton2;
      QPushButton* pushButton3;
      QPushButton* pushButton4;

    public slots:
      virtual void clear();
      virtual void savenclear();
      virtual void lookup();
      virtual void setParamFilename( std::string filename );

    protected:
      OpenMS::Param param_;
      std::string param_filename_;

      QVBoxLayout* layout5;
      QVBoxLayout* layout3;
      QGridLayout* layout2;
      QSpacerItem* spacer1_2;
      QSpacerItem* spacer1;
      QHBoxLayout* layout4;

    protected slots:
      virtual void languageChange();

    };
}
#endif // EDITMODIFICATIONXML_H
