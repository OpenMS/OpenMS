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


#include "SettingsDialog.h"

#include <qvariant.h>
#include <qpushbutton.h>
#include <qgroupbox.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <qlayout.h>
#include <qtooltip.h>
#include <qwhatsthis.h>


#include <qsqldatabase.h>
#include <qmap.h>
#include <qstring.h>
#include <qsqlquery.h>
#include <qstatusbar.h>
#include <qfiledialog.h>
#include "SpecAnnotate.h"


#include <qmessagebox.h>
#include <qmap.h>
#include <qstring.h>
#include <qsettings.h>
#include <qdir.h>
#include "SpecAnnotate.h"
#include "SpectrumMDIWindowEnhanced.h"


using namespace OpenMS;


void SettingsDialog::setParamFilename(std::string filename)
{
  param_filename_ = filename;
}


void SettingsDialog::init()
{
  //rtti: does pa point to instance of SpecAnnotate?
  QWidget* pa = parentWidget();
  if (SpecAnnotate* pa_msa = dynamic_cast<SpecAnnotate*>(pa))
    {
      QMap<QString,QString>* settings;
      settings = pa_msa->getSettings();

      lineEdit1->setText((*settings)["db_username"]);
      //lineEdit2->setText((*settings)["db_password"]);
      lineEdit3->setText((*settings)["db_host"]);
      lineEdit5->setText((*settings)["spl_path"]);
      lineEdit6->setText((*settings)["peakfiles_path"]);
      lineEdit7->setText((*settings)["output_path"]);
    }
}


void SettingsDialog::help()
{
  QMessageBox::information(this, tr("Help: Settings Dialog"),
                           tr("Please insert Settings valid for your system!"), 1);

}


void SettingsDialog::save()
{
  //Load global settings for this application out of file
  main_param_.load(param_filename_);


  //qDebug(QString(param_filename_));


  //main_param_.setValue("Preferences:SpecAnnotate:DB:username", lineEdit1->text());
  //main_param_.setValue("Preferences:SpecAnnotate:DB:password", lineEdit2->text());
  //main_param_.setValue("Preferences:SpecAnnotate:DB:host", lineEdit3->text());

  main_param_.setValue("Preferences:SpecAnnotate:spl_path", lineEdit5->text());
  main_param_.setValue("Preferences:SpecAnnotate:peakfiles_path", lineEdit6->text());
  main_param_.setValue("Preferences:SpecAnnotate:output_path", lineEdit7->text());

  //just a flag to check, whether Preferences for SpecAnnotate are present
  main_param_.setValue("Preferences:SpecAnnotate:present", (std::string)("true"));

  //save back into file the actualized settings for this application
  main_param_.save(param_filename_);
  

  //actualize settings in parent widget
  actualizeParentSettings();
}


void SettingsDialog::ok()
{
  save();
  accept();
}


void SettingsDialog::actualizeParentSettings()
{
  //rtti: does pa point to instance of SpecAnnotate?
  QWidget* pa = parentWidget();
  if (SpecAnnotate* pa_msa = dynamic_cast<SpecAnnotate*>(pa))
    {
      QMap<QString,QString>* settings;
      settings = pa_msa->getSettings();

      (*settings)["db_username"] = lineEdit1->text();
      //(*settings)["db_password"] = lineEdit2->text();
      (*settings)["db_host"] = lineEdit3->text();
      (*settings)["spl_path"] = lineEdit5->text();
      (*settings)["peakfiles_path"] = lineEdit6->text();
      (*settings)["output_path"] = lineEdit7->text();
    }

  //actualize settings in running Instance of SpectrumMDIWindow(Enhanced)
  SpectrumMDIWindowEnhanced::getInstance()->loadPreferences();

}





/*
 *  Constructs a SettingsDialog as a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
SettingsDialog::SettingsDialog( QWidget* parent, const char* name, bool modal, WFlags fl )
    : QDialog( parent, name, modal, fl )
{
  if ( !name )
    setName( "SettingsDialog" );
  QFont f( font() );
  setFont( f );
  setSizeGripEnabled( TRUE );
  SettingsDialogLayout = new QGridLayout( this, 1, 1, 11, 6, "SettingsDialogLayout");

  layout5 = new QHBoxLayout( 0, 0, 6, "layout5");

  buttonHelp = new QPushButton( this, "buttonHelp" );
  buttonHelp->setAutoDefault( TRUE );
  layout5->addWidget( buttonHelp );
  QSpacerItem* spacer = new QSpacerItem( 180, 20, QSizePolicy::Expanding, QSizePolicy::Minimum );
  layout5->addItem( spacer );

  Save = new QPushButton( this, "Save" );
  layout5->addWidget( Save );

  buttonCancel = new QPushButton( this, "buttonCancel" );
  buttonCancel->setAutoDefault( TRUE );
  layout5->addWidget( buttonCancel );

  buttonOk = new QPushButton( this, "buttonOk" );
  QFont buttonOk_font(  buttonOk->font() );
  buttonOk_font.setBold( TRUE );
  buttonOk->setFont( buttonOk_font );
  buttonOk->setAutoDefault( TRUE );
  buttonOk->setDefault( TRUE );
  layout5->addWidget( buttonOk );

  SettingsDialogLayout->addMultiCellLayout( layout5, 1, 1, 0, 1 );

  groupBox1 = new QGroupBox( this, "groupBox1" );
  groupBox1->setColumnLayout(0, Qt::Vertical );
  groupBox1->layout()->setSpacing( 6 );
  groupBox1->layout()->setMargin( 11 );
  groupBox1Layout = new QGridLayout( groupBox1->layout() );
  groupBox1Layout->setAlignment( Qt::AlignTop );

  textLabel1 = new QLabel( groupBox1, "textLabel1" );

  groupBox1Layout->addMultiCellWidget( textLabel1, 0, 0, 0, 1 );

  textLabel2 = new QLabel( groupBox1, "textLabel2" );

  groupBox1Layout->addMultiCellWidget( textLabel2, 1, 1, 0, 1 );

  lineEdit1 = new QLineEdit( groupBox1, "lineEdit1" );
  
  lineEdit1->setEnabled(false);

  groupBox1Layout->addWidget( lineEdit1, 0, 2 );

  //lineEdit2 = new QLineEdit( groupBox1, "lineEdit2" );

  //groupBox1Layout->addWidget( lineEdit2, 1, 2 );

  lineEdit3 = new QLineEdit( groupBox1, "lineEdit3" );

  lineEdit3->setEnabled(false);

  groupBox1Layout->addWidget( lineEdit3, 1, 2 );

  textLabel3 = new QLabel( groupBox1, "textLabel3" );

  groupBox1Layout->addWidget( textLabel3, 1, 0 );
  QSpacerItem* spacer_2 = new QSpacerItem( 40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum );
  groupBox1Layout->addItem( spacer_2, 2, 1 );

  SettingsDialogLayout->addWidget( groupBox1, 0, 0 );

  groupBox2 = new QGroupBox( this, "groupBox2" );
  groupBox2->setSizePolicy( QSizePolicy( (QSizePolicy::SizeType)5, (QSizePolicy::SizeType)5, 2, 0, groupBox2->sizePolicy().hasHeightForWidth() ) );
  groupBox2->setColumnLayout(0, Qt::Vertical );
  groupBox2->layout()->setSpacing( 6 );
  groupBox2->layout()->setMargin( 11 );
  groupBox2Layout = new QGridLayout( groupBox2->layout() );
  groupBox2Layout->setAlignment( Qt::AlignTop );

  textLabel1_2 = new QLabel( groupBox2, "textLabel1_2" );

  groupBox2Layout->addWidget( textLabel1_2, 0, 0 );

  lineEdit5 = new QLineEdit( groupBox2, "lineEdit5" );
  lineEdit5->setMinimumSize( QSize( 150, 0 ) );

  groupBox2Layout->addWidget( lineEdit5, 0, 1 );

  textLabel1_3 = new QLabel( groupBox2, "textLabel1_3" );

  groupBox2Layout->addWidget( textLabel1_3, 1, 0 );

  textLabel2_2 = new QLabel( groupBox2, "textLabel2_2" );

  groupBox2Layout->addWidget( textLabel2_2, 2, 0 );

  lineEdit6 = new QLineEdit( groupBox2, "lineEdit6" );

  groupBox2Layout->addWidget( lineEdit6, 1, 1 );

  lineEdit7 = new QLineEdit( groupBox2, "lineEdit7" );

  groupBox2Layout->addWidget( lineEdit7, 2, 1 );

  SettingsDialogLayout->addWidget( groupBox2, 0, 1 );
  languageChange();
  resize( QSize(823, 175).expandedTo(minimumSizeHint()) );
  clearWState( WState_Polished );

  // signals and slots connections
  connect( buttonOk, SIGNAL( clicked() ), this, SLOT( ok() ) );
  connect( buttonCancel, SIGNAL( clicked() ), this, SLOT( reject() ) );
  connect( buttonHelp, SIGNAL( clicked() ), this, SLOT( help() ) );
  connect( Save, SIGNAL( clicked() ), this, SLOT( save() ) );
  init();
}

/*
 *  Destroys the object and frees any allocated resources
 */
SettingsDialog::~SettingsDialog()
{
  // no need to delete child widgets, Qt does it all for us
}

/*
 *  Sets the strings of the subwidgets using the current
 *  language.
 */
void SettingsDialog::languageChange()
{
  setCaption( tr( "Settings Dialog: Please enter Settings for your System!" ) );
  buttonHelp->setText( tr( "&Help" ) );
  buttonHelp->setAccel( QKeySequence( tr( "F1" ) ) );
  Save->setText( tr( "Save" ) );
  buttonCancel->setText( tr( "&Cancel" ) );
  buttonCancel->setAccel( QKeySequence( QString::null ) );
  buttonOk->setText( tr( "&OK" ) );
  buttonOk->setAccel( QKeySequence( QString::null ) );
  groupBox1->setTitle( tr( "Database (edit in Preferences of TOPPView)" ) );
  textLabel1->setText( tr( "Username" ) );
  //textLabel2->setText( tr( "Password" ) );
  textLabel3->setText( tr( "Host" ) );
  groupBox2->setTitle( tr( "Paths" ) );
  textLabel1_2->setText( tr( ".spl default" ) );
  textLabel1_3->setText( tr( "peakfiles" ) );
  textLabel2_2->setText( tr( "output" ) );
}

