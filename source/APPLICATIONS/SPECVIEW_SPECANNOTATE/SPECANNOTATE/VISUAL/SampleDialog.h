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


#ifndef SAMPLEDIALOG_H
#define SAMPLEDIALOG_H

#include <qvariant.h>
#include <qpixmap.h>
#include <qwidget.h>
#include <qfiledialog.h>
#include <qmap.h>
#include <qstring.h>
#include <string>
#include <vector>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>

class QVBoxLayout;
class QHBoxLayout;
class QGridLayout;
class QSqlDatabase;
class QSqlCursor;
class QSqlForm;
class QGroupBox;
class QLineEdit;
class QPushButton;
class QButtonGroup;
class QRadioButton;
class QLabel;
class QTextEdit;
class QListBox;
class QListBoxItem;
class QComboBox;


namespace OpenMS
  {

  class SpecAnnotate;

  class SampleDialog : public QWidget
    {
      Q_OBJECT

    public:
      SampleDialog( QWidget* parent = 0, const char* name = 0, WFlags fl = 0 );
      ~SampleDialog();

      QGroupBox* groupBox5;
      QLineEdit* lineEdit6;
      QPushButton* pushButton10;
      QPushButton* pushButton2;
      QPushButton* pushButton4;
      QGroupBox* groupBox3;
      QButtonGroup* buttonGroup1;
      QRadioButton* radioButton1;
      QRadioButton* radioButton2;
      QButtonGroup* buttonGroup2;
      QRadioButton* radioButton3;
      QRadioButton* radioButton4;
      QLineEdit* lineEdit3;
      QPushButton* pushButton5_2;
      QLineEdit* lineEdit2;
      QPushButton* pushButton5;
      QLabel* textLabel2_2;
      QLabel* textLabel1_4;
      QGroupBox* groupBox4;
      QPushButton* pushButton7;
      QLabel* textLabel4;
      QLabel* textLabel5;
      QTextEdit* textEdit1;
      QListBox* listBox2;
      QGroupBox* groupBox1;
      QLabel* textLabel2;
      QComboBox* comboBox2;
      QLabel* textLabel1;
      QComboBox* comboBox1;
      QGroupBox* groupBox3_2;
      QLabel* textLabel3_2;
      QLabel* textLabel1_3;
      QLabel* textLabel1_2;
      QComboBox* comboBox3;
      QComboBox* comboBox4;
      QComboBox* comboBox5;
      QLabel* textLabel2_3;
      QLineEdit* lineEdit3_2;
      QPushButton* pushButton1;
      QPushButton* pushButton9;
      QPushButton* pushButton3;

      virtual QString browse( QFileDialog::Mode mode, QString filetype );

    public slots:
      virtual void browsePeakfile();
      virtual void browseOutputdir();
      virtual void quit();
      virtual void loadSample( QString filename );
      virtual void saveSample();
      virtual void saveAs();
      virtual void annotate();
      virtual void inputModifications();
      virtual QString getProtein();
      virtual int getProteinSize();
      virtual void insertPartialMod( QString mod );
      virtual void loadSampleNoDefault();
      virtual void importPeaklistFromFile();
      virtual void importPeaks();
      virtual void exportFiles();
      virtual void exportMetadata();

    protected:
      QVBoxLayout* SampleDialogLayout;
      QHBoxLayout* groupBox5Layout;
      QVBoxLayout* groupBox3Layout;
      QHBoxLayout* layout4;
      QVBoxLayout* buttonGroup1Layout;
      QVBoxLayout* buttonGroup2Layout;
      QGridLayout* layout3;
      QGridLayout* layout4_2;
      QGridLayout* groupBox4Layout;
      QGridLayout* groupBox1Layout;
      QGridLayout* groupBox3_2Layout;
      QHBoxLayout* layout3_2;

    protected slots:
      virtual void languageChange();
      void dbConnect();

    private:
      std::vector<OpenMS::Spectrum1DWidget::Spectrum1D::iterator> peaklist;
      QMap<QString,QString>* settings_;
      SpecAnnotate* pa_msa;
      QSqlDatabase *defaultDB;


      QPixmap image0;

      void init();
      virtual std::string toStlString( QString s );

    };

}

#endif // SAMPLEDIALOG_H
