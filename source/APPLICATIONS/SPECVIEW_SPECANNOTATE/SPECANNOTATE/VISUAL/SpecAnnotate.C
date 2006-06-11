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
// $Id: SpecAnnotate.C,v 1.6 2006/05/30 15:46:40 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:  $
// --------------------------------------------------------------------------



#include "SpecAnnotate.h"

#include "../config_specannotate.h"
#include <fstream>

//QT-Includes
#include <qimage.h>
#include <qpixmap.h>
#include <qtoolbar.h>
#include <qtoolbutton.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qstatusbar.h>
#include <qmessagebox.h>
#include <qapplication.h>
#include <qaccel.h>
#include <qtextstream.h>
#include <qpainter.h>
#include <qpaintdevicemetrics.h>
#include <qwhatsthis.h>
#include <qsqldatabase.h>
#include <qsettings.h>
#include <qdir.h>
#include <qmap.h>
#include <qstring.h>
#include <qsqlquery.h>
#include <qstringlist.h>
#include <qfiledialog.h>
#include <qinputdialog.h>

//My Icons
#include "dbmod.xpm"
#include "dbenz.xpm"
#include "dbami.xpm"
#include "dbprot.xpm"
#include "annotateIcon.xpm"
#include "fileopen.xpm"
#include "filesave.xpm"

//Includes of my User-Interfaces
#include "EditModification.h"
#include "EditEnzyme.h"
#include "EditAminoacid.h"
#include "EditProtein.h"
#include "ViewSequence.h"
#include "ViewDigestFragment.h"
#include "ViewModificationCombination.h"
#include "ViewProteinModificationScenario.h"
#include "ViewRealizedModification.h"
#include "ViewSample.h"
#include "ViewAnnotation.h"
#include "ViewModificationCombinationPositionless.h"
#include "ViewRealizedModificationPositionless.h"

#include "SampleDialog.h"
#include "SettingsDialog.h"

#include "./XML/EditEnzymeXML.h"
#include "./XML/EditModificationXML.h"

#include <OpenMS/DATASTRUCTURES/String.h>


using namespace OpenMS;
using namespace std;



SpecAnnotate::SpecAnnotate()
    : QMainWindow( 0, "SpecAnnotate", WDestructiveClose )
{

  // ************************ Popup Menus *******************************************
  // ***** Sample:
  QPopupMenu * file = new QPopupMenu( this );
  menuBar()->insertItem( tr("&Sample"), file );

  file->insertItem(tr("&New Sample"), this, SLOT(newSample()), CTRL+Key_N);

  QPixmap fileOpenIc = QPixmap(fileopen);

  file->insertItem(fileOpenIc, tr("&Load Sample"), this, SLOT(loadSample()), CTRL+Key_L);

  QPixmap fileSaveIc = QPixmap(filesave);

  file->insertItem(fileSaveIc, tr("&Save Sample"), this, SLOT(saveSample()), CTRL+Key_S);

  QPixmap annotateIc = QPixmap(annotateIcon);

  file->insertItem(annotateIc, tr("&Annotate Sample"), this, SLOT(annotate()), CTRL+Key_A);

  file->insertSeparator();

  file->insertItem(tr("S&ettings"), this, SLOT(openSettingsDialog()), CTRL+Key_E );

  file->insertSeparator();

  file->insertItem(tr("&Quit SpecAnnotate"), this, SLOT( close()), CTRL+Key_Q );


  // ***** Database:
  QPopupMenu * database = new QPopupMenu( this );
#ifndef ANNOTATE_XML

  menuBar()->insertItem( tr("&Database"), database );
#endif

#ifdef ANNOTATE_XML

  database->setEnabled(false);
#endif

  database->insertItem(tr("R&eset Database"), this, SLOT(resetDB()));
  database->insertItem(tr("Setup Database from File"), this, SLOT(setupDB()));
  database->insertItem(tr("Import Protein from FASTA File"), this , SLOT(insertProtIntoDB()));

  database->insertSeparator();
  database->insertSeparator();

  QPixmap dbmodIcon = QPixmap(dbmod);
  database->insertItem(dbmodIcon, tr("&Modification"), this, SLOT(dbMod()));

  QPixmap dbenzIcon = QPixmap(dbenz);
  database->insertItem(dbenzIcon, tr("Enz&yme"), this, SLOT(dbEnz()));

  QPixmap dbamiIcon = QPixmap(dbami);
  database->insertItem(dbamiIcon, tr("&Aminoacid"), this, SLOT(dbAmi()));

  QPixmap dbprotIcon = QPixmap(dbprot);
  database->insertItem(dbprotIcon, tr("&Protein"), this, SLOT(dbProt()));


  database->insertSeparator();
  database->insertSeparator();


  database->insertItem(tr("&Sequence"), this, SLOT(dbSeq()));

  database->insertItem(tr("Di&gest_Fragment"), this, SLOT(dbDigFrag()));

  database->insertSeparator();

  database->insertItem(tr("P&rotein_Modification_Scenario"), this, SLOT(dbProtModScen()));

  database->insertItem(tr("M&odification_Combination"), this, SLOT(dbModComb()));

  database->insertItem(tr("Reali&zed_Modification"), this, SLOT(dbRealMod()));

  database->insertItem(tr("Modification_&Combination_Positionless"), this, SLOT(dbModCombPosless()));

  database->insertItem(tr("Realized_Mo&dification_Positionless"), this, SLOT(dbRealModPosless()));

  database->insertSeparator();

  database->insertItem(tr("Samp&le"), this, SLOT(dbSample()));

  database->insertItem(tr("A&nnotation"), this, SLOT(dbAnnot()));


  // ***** XML:
  QPopupMenu * xmlmenu = new QPopupMenu( this );
  menuBar()->insertItem( tr("&XML"), xmlmenu );

  xmlmenu->insertItem(dbenzIcon, tr("Edit Enzyme"), this, SLOT(xmlEnzyme()));
  xmlmenu->insertItem(dbmodIcon, tr("Edit Modification"), this, SLOT(xmlModification()));
  xmlmenu->insertItem(dbprotIcon, tr("Import Protein from FASTA file"), this, SLOT(insertProtIntoXML()));

  // ***** Help:
  menuBar()->insertSeparator();
  QPopupMenu * help = new QPopupMenu( this );
  menuBar()->insertItem( tr("&Help"), help );
  help->insertItem( tr("&About"), this, SLOT(about()), Key_F1 );
  help->insertItem( tr("About &Qt"), this, SLOT(aboutQt()) );
  help->insertItem( tr("Usage without database"), this, SLOT( usageWithoutDB()));
  help->insertItem( tr("Three annotation methods"), this, SLOT( threeMethods()));
  help->insertSeparator();
  help->insertItem( tr("What's &This"), this, SLOT(whatsThis()), SHIFT+Key_F1 );
  // ************************ End Popup Menus ******************************************


  //load local settings
  settings_ = new QMap<QString,QString>();
  loadSettings();

  resize( 870, 625 );

  //start with creating a new sample
  newSample();

  statusBar()->message( tr("Ready"), 2000 );

}


SpecAnnotate::~SpecAnnotate()
{
  delete settings_;
}




void SpecAnnotate::about()
{
  QMessageBox::about( this, tr("About SpecAnnotate"),
                      tr("This is a small programm to annotate mass spectra of "
                         "modified, e.g. glycosylated, proteins or protein digests \n \n \n"
                         "Copyright (C) 2003-2006 by Andreas Hofmann\n"
                         "deepsun@bioinf.uni-sb.de\n \n"
                         "This program is free software; you can redistribute it and/or modify\n"
                         "it under the terms of the GNU General Public License as published by\n"
                         "the Free Software Foundation; either version 2 of the License, or\n"
                         "(at your option) any later version.\n \n"
                         "This program is distributed in the hope that it will be useful,\n"
                         "but WITHOUT ANY WARRANTY; without even the implied warranty of \n"
                         "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the \n"
                         "GNU General Public License for more details.\n"
                         "You should have received a copy of the GNU General Public License \n"
                         "along with this program; if not, write to the \n"
                         "Free Software Foundation, Inc.,\n"
                         "59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.\n"));
}


void SpecAnnotate::aboutQt()
{
  QMessageBox::aboutQt( this, tr("Qt Application") );
}


void SpecAnnotate::dbMod()
{

  //! check whether connection to database can be established or not
  if ( ! connectToDatabase() )
    {
      statusBar()->message( tr("Could not connect to Database"), 2000 );
    }

  EditModification*  mod_connect = new EditModification(this, "Connection to MySQL Database, Table: modification");
  mod_connect->show();
}


void SpecAnnotate::dbEnz()
{

  //! check whether connection to database can be established or not
  if ( ! connectToDatabase() )
    {
      statusBar()->message( tr("Could not connect to Database"), 2000 );
    }

  EditEnzyme*  enz_connect = new EditEnzyme(this, "Connection to MySQL Database, Table: enzyme");
  enz_connect->show();
}


void SpecAnnotate::dbAmi()
{

  //! check whether connection to database can be established or not
  if ( ! connectToDatabase() )
    {
      statusBar()->message( tr("Could not connect to Database"), 2000 );
    }

  EditAminoacid*  ami_connect = new EditAminoacid(this, "Connection to MySQL Database, Table: aminoacid");
  ami_connect->show();
}


void SpecAnnotate::dbProt()
{

  //! check whether connection to database can be established or not
  if ( ! connectToDatabase() )
    {
      statusBar()->message( tr("Could not connect to Database"), 2000 );
    }

  EditProtein*  prot_connect = new EditProtein(this, "Connection to MySQL Database, Table: protein");
  prot_connect->show();
}


void SpecAnnotate::dbSeq()
{

  //! check whether connection to database can be established or not
  if ( ! connectToDatabase() )
    {
      statusBar()->message( tr("Could not connect to Database"), 2000 );
    }

  ViewSequence* seq_view = new ViewSequence(this, "View of MySQL Database, Table: sequence");
  seq_view->show();
}


void SpecAnnotate::dbDigFrag()
{

  //! check whether connection to database can be established or not
  if ( ! connectToDatabase() )
    {
      statusBar()->message( tr("Could not connect to Database"), 2000 );
    }

  ViewDigestFragment* digfrag_view = new ViewDigestFragment(this, "View of MySQL Database, Table: digest_fragment");
  digfrag_view->show();
}


void SpecAnnotate::dbModComb()
{

  //! check whether connection to database can be established or not
  if ( ! connectToDatabase() )
    {
      statusBar()->message( tr("Could not connect to Database"), 2000 );
    }
  ViewModificationCombination* m_view = new ViewModificationCombination(this,
                                        "View of MySQL Database, Table: modification_combination");
  m_view->show();
}


void SpecAnnotate::dbModCombPosless()
{

  //! check whether connection to database can be established or not
  if ( ! connectToDatabase() )
    {
      statusBar()->message( tr("Could not connect to Database"), 2000 );
    }
  ViewModificationCombinationPositionless* mp_view = new ViewModificationCombinationPositionless(this,
      "View of MySQL Database, Table: modification_combination_positionless");
  mp_view->show();
}


void SpecAnnotate::dbProtModScen()
{

  //! check whether connection to database can be established or not
  if ( ! connectToDatabase() )
    {
      statusBar()->message( tr("Could not connect to Database"), 2000 );
    }
  ViewProteinModificationScenario* protmodscen_view = new ViewProteinModificationScenario(this,
      "View of MySQL Database, Table: protein_modification_scenario");
  protmodscen_view->show();
}


void SpecAnnotate::dbRealMod()
{

  //! check whether connection to database can be established or not
  if ( ! connectToDatabase() )
    {
      statusBar()->message( tr("Could not connect to Database"), 2000 );
    }
  ViewRealizedModification* realmod_view = new ViewRealizedModification(this,
      "View of MySQL Database, Table: realized_modfication");
  realmod_view->show();
}


void SpecAnnotate::dbRealModPosless()
{

  //! check whether connection to database can be established or not
  if ( ! connectToDatabase() )
    {
      statusBar()->message( tr("Could not connect to Database"), 2000 );
    }
  ViewRealizedModificationPositionless* realmodp_view = new ViewRealizedModificationPositionless(this,
      "View of MySQL Database, Table: realized_modfication_positionless");
  realmodp_view->show();
}


void SpecAnnotate::dbSample()
{

  //! check whether connection to database can be established or not
  if ( ! connectToDatabase() )
    {
      statusBar()->message( tr("Could not connect to Database"), 2000 );
    }
  ViewSample* sample_view = new ViewSample(this,"View of MySQL Database, Table: sample");
  sample_view->show();
}


void SpecAnnotate::dbAnnot()
{

  //! check whether connection to database can be established or not
  if ( ! connectToDatabase() )
    {
      statusBar()->message( tr("Could not connect to Database"), 2000 );
    }
  ViewAnnotation* annot_view = new ViewAnnotation(this,"View of MySQL Database, Table: annotation");
  annot_view->show();
}



bool SpecAnnotate::connectToDatabase()
{
  //! the default database connection, using the QTDesigner-MySQL-Driver QTDATABASEDRIVER
  QSqlDatabase *defaultDB = QSqlDatabase::addDatabase( QTDATABASEDRIVER );
  if ( ! defaultDB )
    {
      qWarning( "Failed to connect to driver" );
      return FALSE;
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
      return FALSE;
    }

  return TRUE;
}


void SpecAnnotate::newSample()
{
  SampleDialog* create_sample = new SampleDialog(this,"Sample Dialog: Creating new Sample...");

  create_sample->setFocus();
  setCentralWidget( create_sample );

  create_sample->show();

  statusBar()->message( tr("Opened Dialog for Creating a new Sample..."), 2000 );
}


void SpecAnnotate::loadSettings()
{
  //Load global settings for this application out of file
  main_param_.load(std::string(qApp->argv()[0])+".ini");

  //If no settings for SpecAnnotate Present yet, open a Dialog for entering them
  if ((std::string)(main_param_.getValue("Preferences:SpecAnnotate:present")) != "true")
    {
      openSettingsDialog();
    }

  //fill preferences out of file into \c *settings_
  if ((main_param_.getValue("Preferences:DB:Login").isEmpty()) || (main_param_.getValue("Preferences:DB:Host").isEmpty()))
    {
      QMessageBox::information(this, "Notification:", (std::string)("Please enter database settings before starting SpecAnnotate!"));
      close();
      return;
    }

#ifndef ANNOTATE_XML
  (*settings_)["db_username"] = (std::string)(main_param_.getValue("Preferences:DB:Login"));
  (*settings_)["db_host"] = (std::string)(main_param_.getValue("Preferences:DB:Host"));


  if (main_param_.getValue("DBPassword").isEmpty())
    {
      stringstream ss;
      ss << "Enter database password for user '" <<(*settings_)["db_username"] << "' at '"<< (*settings_)["db_host"];  
      bool ok;
      QString text = QInputDialog::getText("TOPPView Database Password", ss.str(), QLineEdit::Password,QString::null, &ok, this );
      if ( ok ) 
	{
	  main_param_.setValue("DBPassword",text.ascii());
	}
      else
	{
	  close();
	  return;
	}
    }
  
  (*settings_)["db_password"] = (std::string)(main_param_.getValue("DBPassword"));
#else
	(*settings_)["db_username"] = "";
  (*settings_)["db_host"] = "";
  (*settings_)["db_password"] = "";
#endif

  (*settings_)["spl_path"] = (std::string)(main_param_.getValue("Preferences:SpecAnnotate:spl_path"));
  (*settings_)["peakfiles_path"] = (std::string)(main_param_.getValue("Preferences:SpecAnnotate:peakfiles_path"));
  (*settings_)["output_path"] = (std::string)(main_param_.getValue("Preferences:SpecAnnotate:output_path"));
}


void SpecAnnotate::openSettingsDialog()
{
  SettingsDialog* settings_diag = new SettingsDialog(this,"Settings Dialog: Please Enter Settings for your System!");
  settings_diag->setParamFilename(std::string(qApp->argv()[0])+".ini");
  settings_diag->show();
}


QMap<QString,QString>* SpecAnnotate::getSettings()
{
  return settings_;
}


void SpecAnnotate::quit()
{
  qApp->closeAllWindows();
}


void SpecAnnotate::loadSample()
{

  SampleDialog* load_sample = new SampleDialog(this,"Sample Dialog: Creating new Sample...");

  load_sample->setFocus();
  setCentralWidget(load_sample);

  load_sample->loadSample(QString::null);

  load_sample->show();
}

void SpecAnnotate::saveSample()
{
  QWidget* cw = centralWidget();
  SampleDialog* sd;
  if ((sd = dynamic_cast<SampleDialog*>(cw)))  //runtime type identification (rtti): is cw of type SampleDialog* ???
    {
      sd->saveSample();
    }
  else
    {
      statusBar()->message( tr("Could not Save Sample."), 2000 );
    }

}


void SpecAnnotate::resetDB()
{
  //! check whether connection to database can be established or not
  if ( ! connectToDatabase() )
    {
      statusBar()->message( tr("Could not connect to Database"), 2000 );
    }

  QSqlQuery query;

  if (!query.exec("TRUNCATE `annotation`"))
    {
      QMessageBox::warning(this, "Warning!", "Could not execute SQL Query:\n \"TRUNCATE `annotation`\"");
    }

  if (!query.exec("TRUNCATE `modification_combination`"))
    {
      QMessageBox::warning(this, "Warning!", "Could not execute SQL Query:\n \"TRUNCATE `modification_combination`\"");
    }

  if (!query.exec("TRUNCATE `protein_modification_scenario`"))
    {
      QMessageBox::warning(this, "Warning!", "Could not execute SQL Query:\n \"TRUNCATE `protein_modification_scenario`\"");
    }

  if (!query.exec("TRUNCATE `realized_modification`"))
    {
      QMessageBox::warning(this, "Warning!", "Could not execute SQL Query:\n \"TRUNCATE `realized_modification`\"");
    }

  if (!query.exec("TRUNCATE `sample`"))
    {
      QMessageBox::warning(this, "Warning!", "Could not execute SQL Query:\n \"TRUNCATE `sample`\"");
    }

  if (!query.exec("TRUNCATE `modification_combination_positionless`"))
    {
      QMessageBox::warning(this, "Warning!", "Could not execute SQL Query:\n \"TRUNCATE `modification_combination_positionless`\"");
    }

  if (!query.exec("TRUNCATE `realized_modification_positionless`"))
    {
      QMessageBox::warning(this, "Warning!", "Could not execute SQL Query:\n \"TRUNCATE `realized_modification_positionless`\"");
    }

  //erverything ok:
  QMessageBox::information(this, "Notification:", "Database reset successful. Executed Queries:\n \"TRUNCATE `annotation`\" \n \"TRUNCATE `modification_combination`\" \n \"TRUNCATE `protein_modification_scenario`\" \n \"TRUNCATE `realized_modification`\" \n \"TRUNCATE `sample`\" \n \"TRUNCATE `modification_combination_positionless`\" \n \"TRUNCATE `realized_modification_positionless`\"");

}



void SpecAnnotate::annotate()
{
  //get central widget: should be of class SampleDialog
  SampleDialog* cw_sd;
  QWidget* cw = this->centralWidget();
  if ((cw_sd = dynamic_cast<SampleDialog*>(cw)))
    {
      cw_sd->annotate();
    }
  else
    {
      exit(1);
    }
}


void SpecAnnotate::setupDB()
{
  QString q_database_file;
  q_database_file = QFileDialog::getOpenFileName(QString::null, QString::null, this, "",  "Please specify a .sql file containing the database:");

  std::ostringstream ofstr;
  ofstr << q_database_file;
  std::string database_file = ofstr.str();

  bool ok;
  QString database = QInputDialog::getText("Please insert the name of an exisiting database on your server:","QT's database drivers need an exisiting database to log in, so that they can create new databases.\nNo changes will be made to the existing database you specify here.",QLineEdit::Normal,"mysql",&ok,this);

  //! the default database connection, using the QTDesigner-MySQL-Driver QTDATABASEDRIVER
  QSqlDatabase *defaultDB = QSqlDatabase::addDatabase( QTDATABASEDRIVER );
  if ( ! defaultDB )
    {
      qWarning( "Failed to connect to driver" );
    }
  defaultDB->setUserName( (*settings_)["db_username"] );
  defaultDB->setPassword( (*settings_)["db_password"] );
  defaultDB->setHostName( (*settings_)["db_host"] );
  defaultDB->setDatabaseName( database );
  if ( ! defaultDB->open() )
    {
      qWarning( "Failed to connect to SQL server." +
                defaultDB->lastError().driverText() );
      qWarning( defaultDB->lastError().databaseText() );
    }

  QSqlQuery query(defaultDB);

  std::ifstream infile;
  char line[10000];
  infile.open (database_file.c_str());
  while (infile.good())
    {
      infile.getline (line,10000, ';');
      if (strlen(line)!=0)
        query.exec((std::string)line);
    }
  infile.close();

  QMessageBox::information(this, "Notification:", (std::string)("The database has been successfully set up from file \n" + database_file + "!"
                           + "\nPlease restart the program!"));
}



void SpecAnnotate::xmlEnzyme()
{
  EditEnzymeXML* xml_enz = new EditEnzymeXML(this,"edit information on enzymes stored in XML file");
  xml_enz->setParamFilename(std::string(qApp->argv()[0])+".ini");
  xml_enz->show();
}


void SpecAnnotate::xmlModification()
{
  EditModificationXML* xml_mod = new EditModificationXML(this,"edit information on modifications stored in XML file");
  xml_mod->setParamFilename(std::string(qApp->argv()[0])+".ini");
  xml_mod->show();
}



void SpecAnnotate::insertProtIntoDB()
{
  QString q_protein_file;
  q_protein_file = QFileDialog::getOpenFileName(QString::null, QString::null, this, "",  "Please specify a file that contains the protein to be imported in FASTA format:");

  std::ostringstream ofstr;
  ofstr << q_protein_file;
  std::string protein_filename = ofstr.str();

  QString q_identifier;
  bool ok;
  q_identifier = QInputDialog::getText("Input identifier","Please specify an unique identifier for this protein in the database:",QLineEdit::Normal,"",&ok,this);

  std::ostringstream ofstr1;
  ofstr1 << q_identifier;
  std::string identifier = ofstr1.str();

  if ((identifier == "") || (protein_filename == ""))
    {
      QMessageBox::warning(this, "Notification:", (std::string)("No protein has been imported!"));
      return;
    }

  //! reading out of FASTA File
  std::ifstream infile(protein_filename.c_str());
  if (!infile)
    {
      throw OpenMS::Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Wrong filename", ("Could not open file " + protein_filename).c_str() );
    }

  // throw first line away: not interesting
  std::string line;
  getline(infile,line);

  //iterate rest of lines: sequence_oneletter
  string sequence_oneletter;
  for (std::getline(infile,line);!infile.eof();std::getline(infile,line))
    {
      std::istringstream ist(line);
      std::string tmp;
      ist >> tmp;

      sequence_oneletter += tmp;
    }

  //fill some information into database values are calculated taking given filename into account. already present values will be overwritten
  //! check whether connection to database can be established or not
  if ( ! connectToDatabase() )
    {
      statusBar()->message( tr("Could not connect to Database"), 2000 );
    }

  //insert sequence_oneletter into database
  QSqlQuery query;
  query.exec("INSERT INTO " + (string)PROTEIN_TABLE + " ( `identifier` , `fasta_filename` ) "
             + " VALUES ( '" + identifier + "', '" + protein_filename + "' )");

  query.exec("UPDATE " + (string)PROTEIN_TABLE + " SET `sequence_oneletter` = \"" + sequence_oneletter + "\" WHERE "
             + " `identifier` = \"" + identifier + "\"");

  //insert number of aminoacids
  query.exec("UPDATE " + (string)PROTEIN_TABLE + " SET `no_of_aminoacids` = \"" + String(sequence_oneletter.size())
             + "\" WHERE `identifier` = \"" + identifier + "\"");

  QMessageBox::information(this, "Notification:", (std::string)("The protein " + identifier + " has been successfully imported from file "
                           + protein_filename + "."));

}

void SpecAnnotate::insertProtIntoXML()
{
  QString q_protein_file;
  q_protein_file = QFileDialog::getOpenFileName(QString::null, QString::null, this, "",  "Please specify a file that contains the protein to be imported in FASTA format:");

  std::ostringstream ofstr;
  ofstr << q_protein_file;
  std::string protein_filename = ofstr.str();

  QString q_identifier;
  bool ok;
  q_identifier = QInputDialog::getText("Input identifier","Please specify an unique identifier for this protein in the XML file:",QLineEdit::Normal,"",&ok,this);

  std::ostringstream ofstr1;
  ofstr1 << q_identifier;
  std::string identifier = ofstr1.str();

  if ((identifier == "") || (protein_filename == ""))
    {
      QMessageBox::warning(this, "Notification:", (std::string)("No protein has been imported!"));
      return;
    }

  //! reading out of FASTA File
  std::ifstream infile(protein_filename.c_str());
  if (!infile)
    {
      throw OpenMS::Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Wrong filename", ("Could not open file " + protein_filename).c_str() );
    }

  // throw first line away: not interesting
  std::string line;
  getline(infile,line);

  //iterate rest of lines: sequence_oneletter
  string sequence_oneletter;
  for (std::getline(infile,line);!infile.eof();std::getline(infile,line))
    {
      std::istringstream ist(line);
      std::string tmp;
      ist >> tmp;

      sequence_oneletter += tmp;
    }

  //reload XML File, store info into it and save it
  loadSettings();
  main_param_.setValue(("Preferences:SpecAnnotate:Protein:" + identifier + ":sequence_oneletter"), sequence_oneletter );
  main_param_.setValue(("Preferences:SpecAnnotate:Protein:" + identifier + ":no_of_aminoacids"), String(sequence_oneletter.size()));
  main_param_.save(std::string(qApp->argv()[0])+".ini");

  QMessageBox::information(this, "Notification:", (std::string)("The protein " + identifier + " has been successfully imported from file "
                           + protein_filename + "."));

}


void SpecAnnotate::usageWithoutDB()
{

  QMessageBox::information(this, "Help:", (std::string)
                           ((std::string)"If you do not have a database server at hand, you nevertheless can use this program. But with limited functionality.\n" +
                            (std::string)"You just have to define the macro ANNOTATE_XML in file SPECANNOTATE/config_specannotate.h and compile TOPPView again.\n" +
                            (std::string)"Then the database functionality is disabled and you can only use the peakwise_cormen annotation method.\n" +
                            (std::string)"Graphical insertion is also disabled, as well as selection of the protein via a ComboBox. You have to manually insert a \n"+
                            (std::string)"partial modification string and the identifier of the protein. Overall modifications also do not work, but you can emulate\n"+
                            (std::string)"them by setting partial modifications at the appropriate positions.\n"+
                            (std::string)"Information about amino acids is already hardcoded into the file TOPPView.ini. Proteins can be imported by clicking the \n"+
                            (std::string)"corresponding entry in the XML menu. Modifications can be inserted via a dialog that also can be accessed in the XML menu.\n"+
                            (std::string)"In this menu you can also open a dialog for insertion of information about Enzymes, but this has no use yet, since digested \n"+
                            (std::string)"proteins cannot be annotated without a database yet.\n"+
                            (std::string)"If you bear these things in mind, then you can annotate a spectrum of an undigested protein with the peakwise_cormen method \n"+
                            (std::string)"even without a database! Isn't that great?? :-)"));

}


void SpecAnnotate::threeMethods()
{
  QMessageBox::information(this, "Help:", (std::string)
			   ((std::string)"In this program the user can choose between three different annotation methods:\n\n"
			    +(std::string)"The enumerate method:\n"
			    +(std::string)"This method is the most naive and straightforward of the three. It has been implemented, but\n"
			    +(std::string)"its use is not recommended due to the not satisfying performance. This method calculates for\n"
			    +(string)"each digest-fragment all possible modification combinations and stores these, together with their\n"
			    +(string)"masses, in a database. Peak positions are then searched in this database and hopefully annotations\n"
			    +(string)"are found, that yield the same mass as the peak, within the search range. The bad performance of \n"
			    +(string)"this method is a result of two reasons: Masses for each pair of fragment and modification_combination"
                            +(string)"\nare stored in the database. Compared to the following improved_enumerate method, the number of\n"
			    +(string)"calculated masses is way too big and a it's calculation is a lot too time consuming. Besides, in \n"
			    +(string)"determining modification combinations for a fragment this method takes modification positions into\n"
			    +(string)"account. That means that permutations (with respect to the positions) of the same set of\n"
			    +(string)"modifications (that of course yield the same mass) show up as different annotations, again \n"
			    +(string)"increasing the number of masses to be calculated. This also shows another disadvantage of this \n"
			    +(string)"method, many essentially equal annotations (since containing the same modifications, only on \n"
			    +(string)"different positions) show up, which makes finding really interesting annotations a lot more \n"
			    +(string)"difficult for the user.\n\n"
			    +(string)"The improved_enumerate method:\n"
			   +(string)"In this method only each possible modification combination, with its mass is stored in the database.\n"
		         +(string)"For each combination between peak, to be annotated, and digest fragment, to take as unmodified mass to\n"
          	    +(string)"begin with, a seperate database search is executed. this greatly reduces the numer of database entries\n"
		    +(string)"to be calculated. Also this method does not take the actual position of the modification into account,\n"
			    +(string)"and therefore does not enumerate permutations of the same solution.\n"
			    +(string)"This method yields best performances and is recommended.\n\n"
			    +(string)"The peakwise_cormen method:\n"
			    +(string)"This method does not store anything in the database. It starts, for each combination of a peak and \n"
		    +(string)"a digest fragment, an instance of the exact subset sum problem, as seen in the \"Cormen\" algorithms\n"
			    +(string)"textbook. Thus the name. For undigested proteins, this method also works completely without a\n"
			    +(string)"database. Read the \"Usage without database\" help topic for more information."));
}
