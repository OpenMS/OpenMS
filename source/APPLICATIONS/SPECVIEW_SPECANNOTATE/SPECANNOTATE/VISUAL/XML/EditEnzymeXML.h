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

#ifndef EDITENZYMEXML_H
#define EDITENZYMEXML_H

#include <qvariant.h>
#include <qdialog.h>
#include <OpenMS/FORMAT/Param.h>
#include <string>

class QVBoxLayout;
class QHBoxLayout;
class QGridLayout;
class QSpacerItem;
class QLabel;
class QLineEdit;
class QPushButton;


namespace OpenMS
  {

  class EditEnzymeXML : public QDialog
    {
      Q_OBJECT

    public:
      EditEnzymeXML( QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
      ~EditEnzymeXML();

      QLabel* textLabel3;
      QLabel* textLabel2;
      QLineEdit* lineEdit3;
      QLabel* textLabel1;
      QLineEdit* lineEdit2;
      QLineEdit* lineEdit1;
      QPushButton* done;
      QPushButton* lookup;
      QPushButton* clear;
      QPushButton* save;

    public slots:
      virtual void setParamFilename( std::string filename );

    protected:
      QGridLayout* EditEnzymeXMLLayout;
      QGridLayout* layout1;
      QHBoxLayout* layout2;

    protected slots:
      virtual void languageChange();

    private:
      std::string param_filename_;
      OpenMS::Param param_;

    private slots:
      virtual void clear_();
      virtual void savenclear();
      virtual void lookup_();

    };
}
#endif // EDITENZYMEXML_H
