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
// $Id: SampleDialog.C,v 1.4 2006/05/30 15:46:40 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:  $
// --------------------------------------------------------------------------


#include "SampleDialog.h"

#include "../config_specannotate.h"

#include <qvariant.h>
#include <qfile.h>
#include <qtextstream.h>
#include <qdir.h>
#include <qmessagebox.h>
#include <sstream>
#include <ext/hash_map>
#include <iostream>
#include <qsqldatabase.h>
#include <qsqlcursor.h>
#include <qsqlform.h>
#include <qsqlrecord.h>
#include <qpushbutton.h>
#include <qgroupbox.h>
#include <qlineedit.h>
#include <qbuttongroup.h>
#include <qradiobutton.h>
#include <qlabel.h>
#include <qtextedit.h>
#include <qlistbox.h>
#include <qcombobox.h>
#include <qlayout.h>
#include <qtooltip.h>
#include <qwhatsthis.h>
#include <qimage.h>
#include <qpixmap.h>

#include "Annotate.h"
#include "InputModifications.h"
#include "../FUNCTION/string_hash_stl_fixes.h"
#include "../SpectrumMDIWindowEnhanced.h"

#include <qsqldatabase.h>
#include <qmap.h>
#include <qstring.h>
#include <qsqlquery.h>
#include <qstatusbar.h>
#include <qfiledialog.h>
#include "SpecAnnotate.h"

#include "fileopen.xpm"
#include "filesave.xpm"


void SampleDialog::dbConnect()
{
#ifndef ANNOTATE_XML
  //! the default database connection, using the QT-MySQL-Driver QTDATABASEDRIVER
  defaultDB = QSqlDatabase::addDatabase( QTDATABASEDRIVER );

  if ( ! defaultDB )
    {
      qWarning( "Failed to connect to driver" );

      ((QStatusBar*)(pa_msa->statusBar()))->message( tr("Could not connect to Database"), 2000 );
    }
  defaultDB->setDatabaseName( DATABASE );
  defaultDB->setUserName( (*settings_)["db_username"] );
  defaultDB->setPassword( (*settings_)["db_password"] );
  defaultDB->setHostName( (*settings_)["db_host"] );
  if ( ! defaultDB->open() )
    {
      qWarning( "Failed to open database: DATABASE!" +
                defaultDB->lastError().driverText() );
      qWarning( defaultDB->lastError().databaseText() );

      pa_msa->statusBar()->message( tr("Could not connect to Database"), 2000 );
    }
#endif
}


void SampleDialog::init()
{
  //get settings from parent widget (rtti: does pa point to instance of SpecAnnotate?)
  QWidget* pa = this->parentWidget();
  if ((pa_msa = dynamic_cast<SpecAnnotate*>(pa)))
    {
      settings_ = pa_msa->getSettings();
    }
  else
    {
      exit(1);
    }

  //connect to database
  dbConnect();

#ifndef ANNOTATE_XML
  //defaultDB is used automatically!!!!
  QSqlQuery query1( "SELECT identifier FROM protein;" );
  while ( query1.next() )
    {
      comboBox1->insertItem( query1.value( 0 ).toString());
    }

  QSqlQuery query2( "SELECT enzyme_name FROM enzyme;" );
  while ( query2.next() )
    {
      comboBox2->insertItem( query2.value( 0 ).toString());
    }


  QSqlQuery query3( "SELECT modification_name FROM modification ORDER BY modification_ID;" );
  while ( query3.next() )
    {
      listBox2->insertItem( query3.value( 0 ).toString());
    }
#endif

  //insert possibility to select no enzyme
  comboBox2->insertItem("");

  //insert possible annotation methods
#ifndef ANNOTATE_XML

  comboBox3->insertItem("enumerate");
  comboBox3->insertItem("improved_enumerate");
#endif

  comboBox3->insertItem("peakwise_cormen");

  //insert known masstypes
  comboBox4->insertItem("average");
  comboBox4->insertItem("mono");

  //insert known peakfileformats
  comboBox5->insertItem("toll");
  comboBox5->insertItem("kerber");

  //load the default samplefile
  loadSample(((*settings_)["spl_path"] + "default.spl"));
}


QString SampleDialog::browse(QFileDialog::Mode mode, QString filetype)
{
  QString fn;
  if (mode == QFileDialog::ExistingFile)
    {
      if (filetype == "ini")
        {
          fn = QFileDialog::getOpenFileName( (*settings_)["spl_path"], QString::null, this);
        }
      else if (filetype == "peak")
        {
          fn = QFileDialog::getOpenFileName( (*settings_)["peakfiles_path"], QString::null, this);
        }
    }
  else if (mode == QFileDialog::DirectoryOnly)
    {
      fn = QFileDialog::getExistingDirectory( (*settings_)["output_path"], this);
    }
  else
    {
      fn = QString::null;
    }

  if ( !fn.isEmpty() )
    {
      return fn;
    }
  else
    {
      return QString::null;
    }
}


void SampleDialog::browsePeakfile()
{
  QString fn = browse(QFileDialog::ExistingFile, "peak");
  if ((!fn.isEmpty()) && (!fn.isNull()))
    {
      lineEdit2->clear();
      lineEdit2->insert(fn);
    }
}


void SampleDialog::browseOutputdir()
{
  QString fn = browse(QFileDialog::DirectoryOnly, "");
  if (	(!fn.isEmpty()) && (!fn.isNull()))
    {
      lineEdit3->clear();
      lineEdit3->insert(fn);
    }
}


void SampleDialog::quit()
{
  pa_msa->close();
}


void SampleDialog::loadSample(QString filename)
{
  QString fn;
  if (filename == QString::null)
    {
      fn = browse(QFileDialog::ExistingFile, "ini");
      if ((!fn.isEmpty()) && (!fn.isNull()))
        {
          lineEdit6->clear();
          lineEdit6->insert(fn);
        }
    }
  else
    {
      lineEdit6->clear();
      lineEdit6->insert(filename);
      fn = filename;
    }

  //handle ini.file "by hand" with qt means
  QFile file( fn );
  if ( file.open( IO_ReadOnly ) )
    {
      QTextStream stream( &file );
      QString line;
      while ( !stream.atEnd() )
        {
          line = stream.readLine();
          if (line == "[SampleContents]")
            {
              for (int i=0; i<2; i++)
                {
                  line = stream.readLine();
                  if (line.contains("enzyme="))
                    {
                      line.remove("enzyme=");
                      for (int i = 0; i < comboBox2->count(); i++)
                        {
                          if (comboBox2->text(i) == line)
                            {
                              comboBox2->setCurrentItem(i);
                            }
                        }
                    }
                  else if (line.contains("protein="))
                    {
                      line.remove("protein=");
                      for (int i = 0; i < comboBox1->count(); i++)
                        {
                          if (comboBox1->text(i) == line)
                            {
                              comboBox1->setCurrentItem(i);
                            }
                        }
                    }
                }
            }
          else if (line == "[InputOutput]")
            {
              for (int i=0; i<4; i++)
                {
                  line = stream.readLine();
                  if (line.contains("peakfile="))
                    {
                      line.remove("peakfile=");
                      lineEdit2->clear();
                      lineEdit2->insert(line);
                    }
                  else if (line.contains("outputdir="))
                    {
                      line.remove("outputdir=");
                      lineEdit3->clear();
                      lineEdit3->insert(line);
                    }
                  else if (line.contains("using_peakFile=true"))
                    {
                      radioButton2->setChecked(true);
                    }
                  else if (line.contains("using_peakFile=false"))
                    {
                      radioButton1->setChecked(true);
                    }
                  else if (line.contains("using_outputDir=true"))
                    {
                      radioButton4->setChecked(true);
                    }
                  else if (line.contains("using_outputDir=false"))
                    {
                      radioButton3->setChecked(true);
                    }
                }
            }
          else if (line == "[Parameters]")
            {
              for (int i=0; i<4; i++)
                {
                  line = stream.readLine();
                  if (line.contains("search_range="))
                    {
                      line.remove("search_range=");
                      lineEdit3_2->clear();
                      lineEdit3_2->insert(line);
                    }
                  else if (line.contains("peakfile_format="))
                    {
                      line.remove("peakfile_format=");
                      for (int i = 0; i < comboBox5->count(); i++)
                        {
                          if (comboBox5->text(i) == line)
                            {
                              comboBox5->setCurrentItem(i);
                            }
                        }
                    }
                  else if (line.contains("masstype="))
                    {
                      line.remove("masstype=");
                      for (int i = 0; i < comboBox4->count(); i++)
                        {
                          if (comboBox4->text(i) == line)
                            {
                              comboBox4->setCurrentItem(i);
                            }
                        }
                    }
                  else if (line.contains("annotation_method="))
                    {
                      line.remove("annotation_method=");
                      for (int i = 0; i < comboBox3->count(); i++)
                        {
                          if (comboBox3->text(i) == line)
                            {
                              comboBox3->setCurrentItem(i);
                            }
                        }
                    }
                }
            }
          else if (line == "[PartialModifications]")
            {
              line = stream.readLine();
              textEdit1->clear();
              textEdit1->insert(line);
            }
          else if (line == "[OverallModifications]")
            {
              listBox2->clearSelection();
              while (!stream.atEnd())
                {
                  line = stream.readLine();
                  for (uint i = 0; i < listBox2->count(); i++)
                    {
                      if (listBox2->text(i) == line)
                        {
                          listBox2->setSelected(i, true);
                        }
                    }
                }
            }
        }
      file.close();
    }
}


void SampleDialog::saveSample()
{
  QFile file(lineEdit6->text());
  if ( file.open( IO_WriteOnly ) )
    {
      QTextStream stream( &file );
      stream << "[SampleContents]" << endl;
      stream << ("protein=" + comboBox1->currentText()) << endl;
      stream << ("enzyme=" + comboBox2->currentText()) << endl;

      stream << endl << endl;

      stream << "[InputOutput]" << endl;
      stream << ("peakfile=" + lineEdit2->text()) << endl;
      stream << ("outputdir=" + lineEdit3->text()) << endl;
      if (radioButton1->isChecked())
        {
          stream << "using_peakFile=false" << endl;
        }
      else
        {
          stream << "using_peakFile=true" << endl;
        }
      if (radioButton3->isChecked())
        {
          stream << "using_outputDir=false" << endl;
        }
      else
        {
          stream << "using_outputDir=true" << endl;
        }

      stream << endl << endl;

      stream << "[Parameters]" << endl;
      stream << ("search_range=" + lineEdit3_2->text()) << endl;
      stream << ("peakfile_format=" + comboBox5->currentText()) << endl;
      stream << ("masstype=" + comboBox4->currentText()) << endl;
      stream << ("annotation_method=" + comboBox3->currentText()) << endl;

      stream << endl << endl;

      stream << "[PartialModifications]" << endl;
      stream << textEdit1->text() << endl;

      stream << endl << endl;

      stream << "[OverallModifications]";
      for (uint i = 0; i < listBox2->count(); i++)
        {
          if (listBox2->isSelected(i))
            {
              stream << endl;
              stream << listBox2->text(i);
            }
        }
      file.close();
    }

  pa_msa->statusBar()->message( tr("Sample " + lineEdit6->text() + " Saved!"), 2000 );
}


void SampleDialog::saveAs()
{
  QString fn = QFileDialog::getSaveFileName(QString::null, QString::null, this);
  if ((!fn.isEmpty()) && (!fn.isNull()))
    {
      lineEdit6->clear();
      lineEdit6->insert(fn);
    }

  saveSample();
}


void SampleDialog::annotate()
{
  //check if wrong combination of input/output options is chosen
  if ((radioButton2->isChecked()) && (radioButton3->isChecked()))
    {
      QMessageBox::information(this, "Wrong selection", "Reading peaks from file and storing annotations as metadata in spectrum cannot be selected together. \nPlease correct your selection!");
    }
  else
    {
      __gnu_cxx::hash_map<std::string, std::string> sample_data;
      std::vector<std::string> ov_mods;
      sample_data["protein"] =  toStlString(comboBox1->currentText());
      sample_data["enzyme"] =  toStlString(comboBox2->currentText());
      sample_data["peakfile"] = toStlString(lineEdit2->text());

      //if annotations should be exported as metadata: take "" as outputdir
      if (radioButton3->isChecked())
        {
          sample_data["outputdir"] = "";
        }
      else
        {
          sample_data["outputdir"] = toStlString(lineEdit3->text());
        }
      sample_data["search_range"] = toStlString(lineEdit3_2->text());
      sample_data["peakfile_format"] = toStlString(comboBox5->currentText());
      sample_data["masstype"] = toStlString(comboBox4->currentText());
      sample_data["annotation_method"] = toStlString(comboBox3->currentText());
      sample_data["partial_modification_string"] = toStlString(textEdit1->text());

      for (uint i = 0; i < listBox2->count(); i++)
        {
          if (listBox2->isSelected(i))
            {
              ov_mods.push_back(listBox2->text(i));
            }
        }

      //clear peaklist from peaks of previous runs
      peaklist.clear();

      //if corresponding Button is checked: get peaklist from spectrum
      if (radioButton1->isChecked())
        {
          peaklist  = SpectrumMDIWindowEnhanced::getInstance()->getActiveSpectrumSelectedPeaks();
          if (peaklist.empty())
            {
              QMessageBox::information(this, "Missing Peaks", "Reading peaks from spectrum not successful, no peaks selected. \nPlease select peaks first!");
              return;
            }
        }

      Annotate* annotate = new Annotate(this, tr("Annotating Sample " + lineEdit6->text() + "..."));
      annotate->show();
      annotate->run(sample_data, peaklist, ov_mods, settings_);

      //renew Database connection
      /*
      while(QSqlDatabase::contains(QSqlDatabase::defaultConnection))
        {
          QSqlDatabase::removeDatabase(QSqlDatabase::defaultConnection);
        }
      while(QSqlDatabase::contains("db_handle_"))
        {
          QSqlDatabase::removeDatabase("db_handle_");
        }

      dbConnect();
      */
    }
}


void SampleDialog::inputModifications()
{

  InputModifications* inputmod = new InputModifications(this, tr("Partial Modification Inpu"));
  inputmod->show();
}


QString SampleDialog::getProtein()
{
  return comboBox1->currentText();
}


int SampleDialog::getProteinSize()
{
  QSqlQuery query( "SELECT no_of_aminoacids FROM protein WHERE identifier = \"" +  comboBox1->currentText() + "\";" );
  int result = 0;
  if ( query.next() )
    {
      result = query.value( 0 ).toInt();
    }
  return result;
}


void SampleDialog::insertPartialMod(QString mod)
{
  textEdit1->clear();
  textEdit1->insert(mod);
}


std::string SampleDialog::toStlString( QString s )
{
  std::ostringstream ofst;
  ofst << s;
  return ofst.str();
}


void SampleDialog::loadSampleNoDefault()
{
  loadSample(QString::null);
}


void SampleDialog::importPeaklistFromFile()
{
  lineEdit2->setEnabled(true);
  pushButton5->setEnabled(true);
  textLabel2_3->setEnabled(true);
  comboBox5->setEnabled(true);
  textLabel1_4->setEnabled(true);
  radioButton3->setEnabled(false);
  if (radioButton3->isChecked())
    {
      radioButton4->toggle();
    }
}


void SampleDialog::importPeaks()
{
  lineEdit2->setEnabled(false);
  textLabel2_3->setEnabled(false);
  comboBox5->setEnabled(false);
  textLabel1_4->setEnabled(false);
  pushButton5->setEnabled(false);
  radioButton3->setEnabled(true);
}


void SampleDialog::exportFiles()
{
  lineEdit3->setEnabled(true);
  textLabel2_2->setEnabled(true);
  pushButton5_2->setEnabled(true);
}


void SampleDialog::exportMetadata()
{
  lineEdit3->setEnabled(false);
  textLabel2_2->setEnabled(false);
  pushButton5_2->setEnabled(false);
}




static const char* const image0_data[] =
  {
    "32 32 4 1",
    "a c #0000c0",
    "# c #404000",
    "b c #ffffc0",
    ". c #ffffff",
    "..#.......#.....................",
    "..#......#.#....................",
    ".........#......................",
    "..#.####.##.############a.......",
    "..#.#..#.#..#bb#........a.......",
    "..#.#..#.#..#bb#........a.......",
    "..#.#..#.#..####........a.......",
    ".........#..............a.......",
    ".........#..............a.......",
    ".........#..............a.......",
    ".........#..............a.......",
    "........................a.......",
    "........................a.......",
    "........................a..a....",
    "........................a..a....",
    "........................a..a....",
    "....a...................a..a....",
    "....a...................a..a....",
    "....a...................a..a....",
    "....a...................a..a....",
    "....a..............a....a..a....",
    "....a..............a....a..a....",
    "....a..............a....a..a....",
    "....a.a............a....a..a....",
    "....a.a............a....a..a....",
    ".a..a.a.......a....a....a..a....",
    ".a..a.a..a....a....a.a..a..a....",
    ".a..a.a..a....a....a.a..a..a....",
    ".a..a.a..a....a....a.a..a.aa...a",
    ".aa.a.aa.a.a..aaa..a.a.aa.aa...a",
    ".aa.a.aa.a.a..aaa.aaaa.aa.aa..aa",
    "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
  };


/*
 *  Constructs a SampleDialog as a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 */
SampleDialog::SampleDialog( QWidget* parent, const char* name, WFlags fl )
    : QWidget( parent, name, fl ),
    image0( (const char **) image0_data )
{
  if ( !name )
    setName( "SampleDialog" );
  setSizePolicy( QSizePolicy( (QSizePolicy::SizeType)5, (QSizePolicy::SizeType)5, 50, 50, sizePolicy().hasHeightForWidth() ) );
  setMinimumSize( QSize( 50, 50 ) );
  setSizeIncrement( QSize( 1, 1 ) );
  setBaseSize( QSize( 50, 50 ) );
  setFocusPolicy( QWidget::ClickFocus );
  SampleDialogLayout = new QVBoxLayout( this, 11, 6, "SampleDialogLayout");

  groupBox5 = new QGroupBox( this, "groupBox5" );
  groupBox5->setSizePolicy( QSizePolicy( (QSizePolicy::SizeType)5, (QSizePolicy::SizeType)5, 3, 0, groupBox5->sizePolicy().hasHeightForWidth() ) );
  groupBox5->setColumnLayout(0, Qt::Vertical );
  groupBox5->layout()->setSpacing( 6 );
  groupBox5->layout()->setMargin( 11 );
  groupBox5Layout = new QHBoxLayout( groupBox5->layout() );
  groupBox5Layout->setAlignment( Qt::AlignTop );

  lineEdit6 = new QLineEdit( groupBox5, "lineEdit6" );
  groupBox5Layout->addWidget( lineEdit6 );

  pushButton10 = new QPushButton( groupBox5, "pushButton10" );
  QFont pushButton10_font(  pushButton10->font() );
  pushButton10_font.setBold( TRUE );
  pushButton10->setIconSet(QIconSet( fileopen));
  pushButton10->setFont( pushButton10_font );
  groupBox5Layout->addWidget( pushButton10 );

  pushButton2 = new QPushButton( groupBox5, "pushButton2" );
  pushButton2->setIconSet(QIconSet( filesave));
  groupBox5Layout->addWidget( pushButton2 );

  pushButton4 = new QPushButton( groupBox5, "pushButton4" );
  groupBox5Layout->addWidget( pushButton4 );
  pushButton4->setIconSet(QIconSet( filesave));
  SampleDialogLayout->addWidget( groupBox5 );

  groupBox3 = new QGroupBox( this, "groupBox3" );
  groupBox3->setColumnLayout(0, Qt::Vertical );
  groupBox3->layout()->setSpacing( 6 );
  groupBox3->layout()->setMargin( 11 );
  groupBox3Layout = new QVBoxLayout( groupBox3->layout() );
  groupBox3Layout->setAlignment( Qt::AlignTop );

  layout4 = new QHBoxLayout( 0, 0, 6, "layout4");

  buttonGroup1 = new QButtonGroup( groupBox3, "buttonGroup1" );
  buttonGroup1->setColumnLayout(0, Qt::Vertical );
  buttonGroup1->layout()->setSpacing( 6 );
  buttonGroup1->layout()->setMargin( 11 );
  buttonGroup1Layout = new QVBoxLayout( buttonGroup1->layout() );
  buttonGroup1Layout->setAlignment( Qt::AlignTop );

  radioButton1 = new QRadioButton( buttonGroup1, "radioButton1" );
  buttonGroup1Layout->addWidget( radioButton1 );

  radioButton2 = new QRadioButton( buttonGroup1, "radioButton2" );
  buttonGroup1Layout->addWidget( radioButton2 );
  layout4->addWidget( buttonGroup1 );

  buttonGroup2 = new QButtonGroup( groupBox3, "buttonGroup2" );
  buttonGroup2->setColumnLayout(0, Qt::Vertical );
  buttonGroup2->layout()->setSpacing( 6 );
  buttonGroup2->layout()->setMargin( 11 );
  buttonGroup2Layout = new QVBoxLayout( buttonGroup2->layout() );
  buttonGroup2Layout->setAlignment( Qt::AlignTop );

  radioButton3 = new QRadioButton( buttonGroup2, "radioButton3" );
  buttonGroup2Layout->addWidget( radioButton3 );

  radioButton4 = new QRadioButton( buttonGroup2, "radioButton4" );
  buttonGroup2Layout->addWidget( radioButton4 );
  layout4->addWidget( buttonGroup2 );
  groupBox3Layout->addLayout( layout4 );

  layout3 = new QGridLayout( 0, 1, 1, 0, 6, "layout3");

  lineEdit3 = new QLineEdit( groupBox3, "lineEdit3" );

  layout3->addWidget( lineEdit3, 1, 1 );

  pushButton5_2 = new QPushButton( groupBox3, "pushButton5_2" );
  pushButton5_2->setIconSet(QIconSet( fileopen));

  layout3->addWidget( pushButton5_2, 1, 2 );

  lineEdit2 = new QLineEdit( groupBox3, "lineEdit2" );

  layout3->addWidget( lineEdit2, 0, 1 );

  pushButton5 = new QPushButton( groupBox3, "pushButton5" );
  pushButton5->setIconSet(QIconSet( fileopen));

  layout3->addWidget( pushButton5, 0, 2 );

  textLabel2_2 = new QLabel( groupBox3, "textLabel2_2" );

  layout3->addWidget( textLabel2_2, 1, 0 );

  textLabel1_4 = new QLabel( groupBox3, "textLabel1_4" );

  layout3->addWidget( textLabel1_4, 0, 0 );
  groupBox3Layout->addLayout( layout3 );
  SampleDialogLayout->addWidget( groupBox3 );

  layout4_2 = new QGridLayout( 0, 1, 1, 0, 6, "layout4_2");

  groupBox4 = new QGroupBox( this, "groupBox4" );
  groupBox4->setSizePolicy( QSizePolicy( (QSizePolicy::SizeType)5, (QSizePolicy::SizeType)5, 2, 0, groupBox4->sizePolicy().hasHeightForWidth() ) );
  groupBox4->setColumnLayout(0, Qt::Vertical );
  groupBox4->layout()->setSpacing( 6 );
  groupBox4->layout()->setMargin( 11 );
  groupBox4Layout = new QGridLayout( groupBox4->layout() );
  groupBox4Layout->setAlignment( Qt::AlignTop );

  pushButton7 = new QPushButton( groupBox4, "pushButton7" );

  groupBox4Layout->addMultiCellWidget( pushButton7, 0, 0, 1, 2 );

  textLabel4 = new QLabel( groupBox4, "textLabel4" );

  groupBox4Layout->addWidget( textLabel4, 0, 0 );

  textLabel5 = new QLabel( groupBox4, "textLabel5" );

  groupBox4Layout->addWidget( textLabel5, 2, 0 );

  textEdit1 = new QTextEdit( groupBox4, "textEdit1" );

  groupBox4Layout->addMultiCellWidget( textEdit1, 1, 1, 0, 2 );

  listBox2 = new QListBox( groupBox4, "listBox2" );
  listBox2->setSelectionMode( QListBox::Multi );

  groupBox4Layout->addMultiCellWidget( listBox2, 3, 3, 0, 2 );

  layout4_2->addMultiCellWidget( groupBox4, 0, 1, 1, 1 );

  groupBox1 = new QGroupBox( this, "groupBox1" );
  groupBox1->setColumnLayout(0, Qt::Vertical );
  groupBox1->layout()->setSpacing( 6 );
  groupBox1->layout()->setMargin( 11 );
  groupBox1Layout = new QGridLayout( groupBox1->layout() );
  groupBox1Layout->setAlignment( Qt::AlignTop );

  textLabel2 = new QLabel( groupBox1, "textLabel2" );

  groupBox1Layout->addWidget( textLabel2, 1, 0 );

  comboBox2 = new QComboBox( FALSE, groupBox1, "comboBox2" );

  groupBox1Layout->addWidget( comboBox2, 1, 1 );

  textLabel1 = new QLabel( groupBox1, "textLabel1" );

  groupBox1Layout->addWidget( textLabel1, 0, 0 );

  comboBox1 = new QComboBox( FALSE, groupBox1, "comboBox1" );

  groupBox1Layout->addWidget( comboBox1, 0, 1 );

  layout4_2->addWidget( groupBox1, 0, 0 );

  groupBox3_2 = new QGroupBox( this, "groupBox3_2" );
  groupBox3_2->setColumnLayout(0, Qt::Vertical );
  groupBox3_2->layout()->setSpacing( 6 );
  groupBox3_2->layout()->setMargin( 11 );
  groupBox3_2Layout = new QGridLayout( groupBox3_2->layout() );
  groupBox3_2Layout->setAlignment( Qt::AlignTop );

  textLabel3_2 = new QLabel( groupBox3_2, "textLabel3_2" );

  groupBox3_2Layout->addWidget( textLabel3_2, 2, 0 );

  textLabel1_3 = new QLabel( groupBox3_2, "textLabel1_3" );

  groupBox3_2Layout->addWidget( textLabel1_3, 3, 0 );

  textLabel1_2 = new QLabel( groupBox3_2, "textLabel1_2" );

  groupBox3_2Layout->addWidget( textLabel1_2, 0, 0 );

  comboBox3 = new QComboBox( FALSE, groupBox3_2, "comboBox3" );

  groupBox3_2Layout->addWidget( comboBox3, 3, 1 );

  comboBox4 = new QComboBox( FALSE, groupBox3_2, "comboBox4" );

  groupBox3_2Layout->addWidget( comboBox4, 2, 1 );

  comboBox5 = new QComboBox( FALSE, groupBox3_2, "comboBox5" );

  groupBox3_2Layout->addWidget( comboBox5, 1, 1 );

  textLabel2_3 = new QLabel( groupBox3_2, "textLabel2_3" );

  groupBox3_2Layout->addWidget( textLabel2_3, 1, 0 );

  lineEdit3_2 = new QLineEdit( groupBox3_2, "lineEdit3_2" );

  groupBox3_2Layout->addWidget( lineEdit3_2, 0, 1 );

  layout4_2->addWidget( groupBox3_2, 1, 0 );
  SampleDialogLayout->addLayout( layout4_2 );

  layout3_2 = new QHBoxLayout( 0, 0, 6, "layout3_2");

  pushButton1 = new QPushButton( this, "pushButton1" );
  layout3_2->addWidget( pushButton1 );

  pushButton9 = new QPushButton( this, "pushButton9" );
  layout3_2->addWidget( pushButton9 );
  QSpacerItem* spacer = new QSpacerItem( 370, 29, QSizePolicy::Expanding, QSizePolicy::Minimum );
  layout3_2->addItem( spacer );

  pushButton3 = new QPushButton( this, "pushButton3" );
  QFont pushButton3_font(  pushButton3->font() );
  pushButton3_font.setBold( TRUE );
  pushButton3->setFont( pushButton3_font );
  pushButton3->setDefault( TRUE );
  pushButton3->setIconSet( QIconSet( image0 ) );
  pushButton3->setFlat( FALSE );
  layout3_2->addWidget( pushButton3 );
  SampleDialogLayout->addLayout( layout3_2 );

  languageChange();
  resize( QSize(910, 600).expandedTo(minimumSizeHint()) );
  clearWState( WState_Polished );

  // signals and slots connections
  connect( pushButton1, SIGNAL( clicked() ), this, SLOT( close() ) );
  connect( pushButton5, SIGNAL( clicked() ), this, SLOT( browsePeakfile() ) );
  connect( pushButton5_2, SIGNAL( clicked() ), this, SLOT( browseOutputdir() ) );
  connect( pushButton9, SIGNAL( clicked() ), this, SLOT( quit() ) );
  connect( pushButton10, SIGNAL( clicked() ), this, SLOT( loadSampleNoDefault() ) );
  connect( pushButton2, SIGNAL( clicked() ), this, SLOT( saveSample() ) );
  connect( pushButton4, SIGNAL( clicked() ), this, SLOT( saveAs() ) );
  connect( pushButton3, SIGNAL( clicked() ), this, SLOT( annotate() ) );
  connect( pushButton7, SIGNAL( clicked() ), this, SLOT( inputModifications() ) );
  connect( radioButton1, SIGNAL( toggled(bool) ), this, SLOT( importPeaks() ) );
  connect( radioButton2, SIGNAL( toggled(bool) ), this, SLOT( importPeaklistFromFile() ) );
  connect( radioButton3, SIGNAL( toggled(bool) ), this, SLOT( exportMetadata() ) );
  connect( radioButton4, SIGNAL( toggled(bool) ), this, SLOT( exportFiles() ) );

  // tab order
  setTabOrder( lineEdit6, comboBox1 );
  setTabOrder( comboBox1, comboBox2 );
  setTabOrder( comboBox2, lineEdit2 );
  setTabOrder( lineEdit2, lineEdit3 );
  setTabOrder( lineEdit3, lineEdit3_2 );
  setTabOrder( lineEdit3_2, textEdit1 );
  setTabOrder( textEdit1, listBox2 );
  setTabOrder( listBox2, pushButton3 );
  setTabOrder( pushButton3, pushButton1 );
  setTabOrder( pushButton1, pushButton10 );
  setTabOrder( pushButton10, pushButton2 );
  setTabOrder( pushButton2, pushButton4 );
  setTabOrder( pushButton4, pushButton5_2 );
  setTabOrder( pushButton5_2, pushButton5 );
  setTabOrder( pushButton5, pushButton7 );


#ifdef ANNOTATE_XML

  comboBox2->setEnabled(false);
  comboBox1->setEditable(true);
  pushButton7->setEnabled(false);
  textLabel2->setEnabled(false);
  textLabel5->setEnabled(false);
  listBox2->setEnabled(false);
#endif


  init();
}

/*
 *  Destroys the object and frees any allocated resources
 */
SampleDialog::~SampleDialog()
{
  // no need to delete child widgets, Qt does it all for us
}

/*
 *  Sets the strings of the subwidgets using the current
 *  language.
 */
void SampleDialog::languageChange()
{
  setCaption( tr( "Sample Dialog" ) );
  groupBox5->setTitle( tr( "Sample File" ) );
  QToolTip::add
    ( groupBox5, tr( "All entries in this dialog can be saved in a file, the \"Sample File\"" ) );
  QWhatsThis::add
    ( groupBox5, tr( "All entries in this dialog can be saved in a file, the \"Sample File\"" ) );
  pushButton10->setText( tr( "Load" ) );
  pushButton2->setText( tr( "Save" ) );
  pushButton4->setText( tr( "Save As" ) );
  groupBox3->setTitle( tr( "Input and Output" ) );
  QToolTip::add
    ( groupBox3, QString::null );
  QWhatsThis::add
    ( groupBox3, QString::null );
  buttonGroup1->setTitle( tr( "Input: Peaklist" ) );
  QToolTip::add
    ( buttonGroup1, tr( "In this box the user can decide, whether selected peaks in currend active spectrum in TOPPView schould be annotated, or whether a peaklist should be read out of a file" ) );
  QWhatsThis::add
    ( buttonGroup1, tr( "In this box the user can decide, whether selected peaks in currend active spectrum in TOPPView schould be annotated, or whether a peaklist should be read out of a file" ) );
  radioButton1->setText( tr( "Import Selected Peaks from active Spectrum" ) );
  radioButton2->setText( tr( "Use Peaklist in File" ) );
  buttonGroup2->setTitle( tr( "Output: Annotations" ) );
  QToolTip::add
    ( buttonGroup2, tr( "In this box the user can decide, whether found annotations should be written in one file per peak, or if they should be returned to the active spectrum as metadata" ) );
  QWhatsThis::add
    ( buttonGroup2, tr( "In this box the user can decide, whether found annotations should be written in one file per peak, or if they should be returned to the active spectrum as metadata" ) );
  radioButton3->setText( tr( "Store Annotations as Metadata in Spectrum" ) );
  radioButton4->setText( tr( "Export Annotations into Files (in Output Directory)" ) );
  QToolTip::add
    ( lineEdit3, tr( "Here the user can be specify the directory in which the output files should be created" ) );
  QWhatsThis::add
    ( lineEdit3, tr( "Here the user can be specify the directory in which the output files should be created" ) );
  pushButton5_2->setText( tr( "Browse" ) );
  QToolTip::add
    ( pushButton5_2, tr( "Here the user can be specify the directory in which the output files should be created" ) );
  QWhatsThis::add
    ( pushButton5_2, tr( "Here the user can be specify the directory in which the output files should be created" ) );
  QToolTip::add
    ( lineEdit2, tr( "Here the user can specify the peakfile to use" ) );
  QWhatsThis::add
    ( lineEdit2, tr( "Here the user can specify the peakfile to use" ) );
  pushButton5->setText( tr( "Browse" ) );
  QToolTip::add
    ( pushButton5, tr( "Here the user can specify the peakfile to use" ) );
  QWhatsThis::add
    ( pushButton5, tr( "Here the user can specify the peakfile to use" ) );
  textLabel2_2->setText( tr( "Output Directory" ) );
  QToolTip::add
    ( textLabel2_2, tr( "Here the user can be specify the directory in which the output files should be created" ) );
  QWhatsThis::add
    ( textLabel2_2, tr( "Here the user can be specify the directory in which the output files should be created" ) );
  textLabel1_4->setText( tr( "Peaklist File" ) );
  QToolTip::add
    ( textLabel1_4, tr( "Here the user can specify the peakfile to use" ) );
  QWhatsThis::add
    ( textLabel1_4, tr( "Here the user can specify the peakfile to use" ) );
  groupBox4->setTitle( tr( "Modifications" ) );
  pushButton7->setText( tr( "Insert Graphically" ) );
  QToolTip::add
    ( pushButton7, tr( "Opens  a dialog that helps you with the input of partial modifications" ) );
  QWhatsThis::add
    ( pushButton7, tr( "Opens  a dialog that helps you with the input of partial modificationss" ) );
  textLabel4->setText( tr( "Partial Modifications" ) );
  QToolTip::add
    ( textLabel4, tr( "Partial Modifications are modifications of the type \"Position x,y and z possibly can be modified with modifications A, B, C, or D ..." ) );
  QWhatsThis::add
    ( textLabel4, tr( "Partial Modifications are modifications of the type \"Position x,y and z possibly can be modified with modifications A, B, C, or D ..." ) );
  textLabel5->setText( tr( "Overall Modifications" ) );
  QToolTip::add
    ( textLabel5, tr( "Overall Modifications are modifications of following type: All cysteines are alkylated..." ) );
  QWhatsThis::add
    ( textLabel5, tr( "Overall Modifications are modifications of following type: All cysteines are alkylated..." ) );
  QToolTip::add
    ( textEdit1, tr( "Partial Modifications are modifications of the type \"Position x,y and z possibly can be modified with modifications A, B, C, or D ..." ) );
  QWhatsThis::add
    ( textEdit1, tr( "Partial Modifications are modifications of the type \"Position x,y and z possibly can be modified with modifications A, B, C, or D ..." ) );
  QToolTip::add
    ( listBox2, tr( "Overall Modifications are modifications of following type: All cysteines are alkylated..." ) );
  QWhatsThis::add
    ( listBox2, tr( "Overall Modifications are modifications of following type: All cysteines are alkylated..." ) );
  groupBox1->setTitle( tr( "Sample Contents" ) );
  textLabel2->setText( tr( "Enzyme" ) );
  QToolTip::add
    ( textLabel2, tr( "What enzyme is used for calculation of theoretical annotations? The user can add new items into the combo boxes by updating the database for proteins and enzymes" ) );
  QWhatsThis::add
    ( textLabel2, tr( "What enzyme is used for calculation of theoretical annotations? The user can add new items into the combo boxes by updating the database for proteins and enzymes" ) );
  QToolTip::add
    ( comboBox2, tr( "What enzyme is used for calculation of theoretical annotations? The user can add new items into the combo boxes by updating the database for proteins and enzymes" ) );
  QWhatsThis::add
    ( comboBox2, tr( "What enzyme is used for calculation of theoretical annotations? The user can add new items into the combo boxes by updating the database for proteins and enzymes" ) );
  textLabel1->setText( tr( "Protein" ) );
  QToolTip::add
    ( textLabel1, tr( "What protein is used for calculation of theoretical annotations? The user can add new items into the combo boxes by updating the database for proteins and enzymes" ) );
  QWhatsThis::add
    ( textLabel1, tr( "What protein is used for calculation of theoretical annotations? The user can add new items into the combo boxes by updating the database for proteins and enzymes" ) );
  QToolTip::add
    ( comboBox1, tr( "What protein is used for calculation of theoretical annotations? The user can add new items into the combo boxes by updating the database for proteins and enzymes" ) );
  QWhatsThis::add
    ( comboBox1, tr( "What protein is used for calculation of theoretical annotations? The user can add new items into the combo boxes by updating the database for proteins and enzymes" ) );
  groupBox3_2->setTitle( tr( "Parameters" ) );
  QWhatsThis::add
    ( groupBox3_2, tr( "ss" ) );
  textLabel3_2->setText( tr( "Mass Type" ) );
  QToolTip::add
    ( textLabel3_2, tr( "Mass Type specifies whether average molecular or monoisotopic molecular masses should be used" ) );
  QWhatsThis::add
    ( textLabel3_2, tr( "Mass Type specifies whether average molecular or monoisotopic molecular masses should be used" ) );
  textLabel1_3->setText( tr( "Annot. Method" ) );
  QToolTip::add
    ( textLabel1_3, tr( "Annotation Method specifies how annot. are calculated. \"enumerate\" and \"improved_enumerate\" store annotations in the database. \"improved_enumerate\" is recommended. " ) );
  QWhatsThis::add
    ( textLabel1_3, tr( "Annotation Method specifies how annot. are calculated. \"enumerate\" and \"improved_enumerate\" store annotations in the database. \"improved_enumerate\" is recommended. " ) );
  textLabel1_2->setText( tr( "Search Range" ) );
  QToolTip::add
    ( textLabel1_2, tr( "Search Range speciefies + / - how many Daltons a theoretically calculated annotation may differ from peakvalue to be recognized as annotation for this peak" ) );
  QWhatsThis::add
    ( textLabel1_2, tr( "Search Range speciefies + / - how many Daltons a theoretically calculated annotation may differ from peakvalue to be recognized as annotation for this peak" ) );
  QToolTip::add
    ( comboBox3, tr( "Annotation Method specifies how annot. are calculated. \"enumerate\" and \"improved_enumerate\" store annotations in the database. \"improved_enumerate\" is recommended. " ) );
  QWhatsThis::add
    ( comboBox3, tr( "Annotation Method specifies how annot. are calculated. \"enumerate\" and \"improved_enumerate\" store annotations in the database. \"improved_enumerate\" is recommended. " ) );
  QToolTip::add
    ( comboBox4, tr( "Mass Type specifies whether average molecular or monoisotopic molecular masses should be used" ) );
  QWhatsThis::add
    ( comboBox4, tr( "Mass Type specifies whether average molecular or monoisotopic molecular masses should be used" ) );
  QToolTip::add
    ( comboBox5, tr( "Peakfile Format specifies in what format the peaks are stored in the peakfile" ) );
  QWhatsThis::add
    ( comboBox5, tr( "Peakfile Format specifies in what format the peaks are stored in the peakfile" ) );
  textLabel2_3->setText( tr( "Peakfile Format" ) );
  QToolTip::add
    ( textLabel2_3, tr( "Peakfile Format specifies in what format the peaks are stored in the peakfile" ) );
  QWhatsThis::add
    ( textLabel2_3, tr( "Peakfile Format specifies in what format the peaks are stored in the peakfile" ) );
  QToolTip::add
    ( lineEdit3_2, tr( "Search Range speciefies + / - how many Daltons a theoretically calculated annotation may differ from peakvalue to be recognized as annotation for this peak" ) );
  QWhatsThis::add
    ( lineEdit3_2, tr( "Search Range speciefies + / - how many Daltons a theoretically calculated annotation may differ from peakvalue to be recognized as annotation for this peak" ) );
  pushButton1->setText( tr( "Cancel" ) );
  QToolTip::add
    ( pushButton1, tr( "This Button closes actual dialog widget in this window. If last dialog widget is closed, new dialog widget is opened with \"New Sample\" in the \"Sample\" menu." ) );
  QWhatsThis::add
    ( pushButton1, tr( "This Button closes actual dialog widget in this window. If last dialog widget is closed, new dialog widget is opened with \"New Sample\" in the \"Sample\" menu." ) );
  pushButton9->setText( tr( "Quit SpecAnnotate" ) );
  QToolTip::add
    ( pushButton9, tr( "Closes the annotation part of TOPPView" ) );
  QWhatsThis::add
    ( pushButton9, tr( "Closes the annotation part of TOPPView" ) );
  pushButton3->setText( tr( "Annotate" ) );
  QToolTip::add
    ( pushButton3, tr( "Starts annotation of peaks" ) );
}

