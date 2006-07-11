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
// $Maintainer:  $
// --------------------------------------------------------------------------


#ifndef SETTINGSDIALOG_H
#define SETTINGSDIALOG_H

//QT includes
#include <qvariant.h>
#include <qdialog.h>
#include <qmainwindow.h>

//OpenMS includes
#include <OpenMS/FORMAT/Param.h>


class QVBoxLayout;
class QHBoxLayout;
class QGridLayout;
class QPushButton;
class QGroupBox;
class QLabel;
class QLineEdit;


namespace OpenMS
  {

  class SettingsDialog : public QDialog
    {
      Q_OBJECT

    public:
      SettingsDialog( QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
      ~SettingsDialog();

      QPushButton* buttonHelp;
      QPushButton* Save;
      QPushButton* buttonCancel;
      QPushButton* buttonOk;
      QGroupBox* groupBox1;
      QLabel* textLabel1;
      QLabel* textLabel2;
      QLineEdit* lineEdit1;
      // QLineEdit* lineEdit2;
      QLineEdit* lineEdit3;
      QLabel* textLabel3;
      QGroupBox* groupBox2;
      QLabel* textLabel1_2;
      QLineEdit* lineEdit5;
      QLabel* textLabel1_3;
      QLabel* textLabel2_2;
      QLineEdit* lineEdit6;
      QLineEdit* lineEdit7;

    public slots:
      virtual void help();
      virtual void save();
      virtual void ok();
      void setParamFilename(std::string filename);


    protected:
      QGridLayout* SettingsDialogLayout;
      QHBoxLayout* layout5;
      QGridLayout* groupBox1Layout;
      QGridLayout* groupBox2Layout;

    protected slots:
      virtual void languageChange();

    private:
      void init();
      void actualizeParentSettings();

      OpenMS::Param main_param_;
      std::string param_filename_;

    };
}

#endif // SETTINGSDIALOG_H
