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


#include "Annotate.h"

#include "../config_specannotate.h"

#include <qvariant.h>
#include <qdir.h>
#include <qtimer.h>
#include <qevent.h>
#include <qsqlquery.h>
#include <qsqlresult.h>
#include <qsqldatabase.h>
#include <qmessagebox.h>
#include <iostream>
#include <qpushbutton.h>
#include <qlabel.h>
#include <qtextbrowser.h>
#include <qlcdnumber.h>
#include <qlayout.h>
#include <qtooltip.h>
#include <qwhatsthis.h>
#include "CustomEvents.h"

using namespace OpenMS;

void Annotate::run(__gnu_cxx::hash_map<std::string, std::string> sample_data, std::vector<OpenMS::Spectrum1DWidget::Spectrum1D::iterator>& peaklist, std::vector<std::string> ov_mods, QMap<QString, QString> * settings )
{
  //fill member variable
  settings_ = settings;
  setCaption("Annotating...");

  // for measuring the time
  t.start();                               //for counting time since start
  timer_ID = startTimer(1000);  //starting first one-second-timer for showing elapsed time every second
  db_display_update_timer = startTimer(60000); //first one-minute timer for updating db info
  updateDBDisplay();

  //creating  a new Annotation-Thread (object is not deleted until window is closed)
  qathread = new AnnotateThread(sample_data, peaklist, ov_mods, (*settings)["db_username"] , (*settings)["db_password"] , (*settings)["db_host"], this);

  qathread->start();
}


void Annotate::addOutput(std::string s)
{
  textBrowser1->append(QString(s));
  update();
}


void Annotate::ready()
{
  killTimer(timer_ID);
  killTimer(db_display_update_timer);
  killTimers();
  textBrowser1->setContentsPos(0, 0);
  updateDBDisplay();

  //clean up database
  /*
  while(QSqlDatabase::contains(QSqlDatabase::defaultConnection))
    {
  QSqlDatabase::removeDatabase(QSqlDatabase::defaultConnection);
    }
  while(QSqlDatabase::contains("db_handle_"))
    {
  QSqlDatabase::removeDatabase("db_handle_");
    }
  */
}


void Annotate::abort()
{
  //terminate annotation-thread
  if (qathread->running())
    {
      qathread->terminate();
      qathread->wait();
      killTimer(timer_ID);
      killTimer(db_display_update_timer);
      textBrowser1->setContentsPos(0, 0);
      updateDBDisplay();
      QMessageBox::information(this, tr("Warning:"), tr("Annotation aborted by user!"));
    }
}


void Annotate::closeWindow()
{
  //clean up
  if (qathread->running())
    {
      abort();
    }
  close();
}


void Annotate::timerEvent( QTimerEvent * e )
{
  if (e->timerId() == timer_ID)
    {
      //update lcd display
      QTime elapsed(0,0,0,0);
      elapsed = elapsed.addMSecs(t.elapsed());
      lCDNumber1->display(elapsed.toString()); //elapsed.toString("hh:mm:ss"));
      update();

      timer_ID = startTimer(1000);

    }
  else if (e->timerId() == db_display_update_timer)
    {
      updateDBDisplay();
    }
}



void Annotate::dbConnect()
{
#ifndef ANNOTATE_XML

  defaultDB = QSqlDatabase::addDatabase( QTDATABASEDRIVER );

  if ( ! defaultDB )
    {
      qWarning( "Failed to connect to driver" );
      return;
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
      return;
    }
#endif
}



void Annotate::updateDBDisplay()
{
#ifndef ANNOTATE_XML

  //connect to database for updating the display
  dbConnect();

  //get instance of Database (defaultConnection)
  defaultDB = QSqlDatabase::database();

  //execute queries
  QSqlQuery annotations = defaultDB->exec("SELECT count(*) FROM annotation");
  QSqlQuery mod_comb = defaultDB->exec("SELECT count(*) FROM modification_combination");
  QSqlQuery real_mod = defaultDB->exec("SELECT count(*) FROM realized_modification");
  QSqlQuery mod_comb_posless = defaultDB->exec("SELECT count(*) FROM modification_combination_positionless");
  QSqlQuery real_mod_posless = defaultDB->exec("SELECT count(*) FROM realized_modification_positionless");

  //process them and update lcd dislplays
  if (annotations.isActive())
    {
      annotations.next(); //first state is BEFORE value that we want
      lCDNumber2->display(annotations.value(0).toString());
    }
  if (mod_comb.isActive())
    {
      mod_comb.next(); //first state is BEFORE value that we want
      lCDNumber3->display(mod_comb.value(0).toString());
    }
  if (real_mod.isActive())
    {
      real_mod.next(); //first state is BEFORE value that we want
      lCDNumber4->display(real_mod.value(0).toString());
    }
  if (mod_comb_posless.isActive())
    {
      mod_comb_posless.next(); //first state is BEFORE value that we want
      lCDNumber5->display(mod_comb_posless.value(0).toString());
    }
  if (real_mod_posless.isActive())
    {
      real_mod_posless.next(); //first state is BEFORE value that we want
      lCDNumber6->display(real_mod_posless.value(0).toString());
    }

  //schedule repaint
  update();

  //start new one-minute-timer for updating database info
  db_display_update_timer = startTimer(60000);
#endif
}


void Annotate::customEvent( QCustomEvent * e )
{
  if (e->type() == 65432)   //an OutputEvent
    {
      OutputEvent* ue = (OutputEvent*) e;
      addOutput(ue->output());
    }
  else if (e->type() == 65433)  //a FinishEvent
    {
      ready();
      QMessageBox::information(this, tr("Notification:"), tr("Annotation of Peaks finished!"));
    }
}



/*
 *  Constructs a Annotate as a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
Annotate::Annotate( QWidget* parent, const char* name, bool modal, WFlags fl )
    : QDialog( parent, name, modal, fl )
{
  if ( !name )
    setName( "Annotate" );
  AnnotateLayout = new QVBoxLayout( this, 11, 6, "AnnotateLayout");

  layout3 = new QVBoxLayout( 0, 0, 6, "layout3");

  textLabel2 = new QLabel( this, "textLabel2" );
  layout3->addWidget( textLabel2 );

  textBrowser1 = new QTextBrowser( this, "textBrowser1" );
  layout3->addWidget( textBrowser1 );
  AnnotateLayout->addLayout( layout3 );

  layout5 = new QHBoxLayout( 0, 0, 6, "layout5");

  textLabel1 = new QLabel( this, "textLabel1" );
  layout5->addWidget( textLabel1 );

  lCDNumber1 = new QLCDNumber( this, "lCDNumber1" );
  lCDNumber1->setNumDigits( 8 );
  lCDNumber1->setSegmentStyle( QLCDNumber::Outline );
  layout5->addWidget( lCDNumber1 );
  AnnotateLayout->addLayout( layout5 );

  layout4 = new QHBoxLayout( 0, 0, 6, "layout4");

  textLabel1_2 = new QLabel( this, "textLabel1_2" );
  layout4->addWidget( textLabel1_2 );

  lCDNumber2 = new QLCDNumber( this, "lCDNumber2" );
  lCDNumber2->setNumDigits( 15 );
  layout4->addWidget( lCDNumber2 );
  AnnotateLayout->addLayout( layout4 );

  layout5_2 = new QHBoxLayout( 0, 0, 6, "layout5_2");

  textLabel2_2 = new QLabel( this, "textLabel2_2" );
  layout5_2->addWidget( textLabel2_2 );

  lCDNumber3 = new QLCDNumber( this, "lCDNumber3" );
  lCDNumber3->setNumDigits( 15 );
  layout5_2->addWidget( lCDNumber3 );
  AnnotateLayout->addLayout( layout5_2 );

  layout6 = new QHBoxLayout( 0, 0, 6, "layout6");

  textLabel3 = new QLabel( this, "textLabel3" );
  layout6->addWidget( textLabel3 );

  lCDNumber4 = new QLCDNumber( this, "lCDNumber4" );
  lCDNumber4->setNumDigits( 15 );
  layout6->addWidget( lCDNumber4 );
  AnnotateLayout->addLayout( layout6 );

  layout8 = new QHBoxLayout( 0, 0, 6, "layout8");

  textLabel1_3 = new QLabel( this, "textLabel1_3" );
  layout8->addWidget( textLabel1_3 );

  lCDNumber5 = new QLCDNumber( this, "lCDNumber5" );
  lCDNumber5->setNumDigits( 15 );
  layout8->addWidget( lCDNumber5 );
  AnnotateLayout->addLayout( layout8 );

  layout9 = new QHBoxLayout( 0, 0, 6, "layout9");

  textLabel2_3 = new QLabel( this, "textLabel2_3" );
  layout9->addWidget( textLabel2_3 );

  lCDNumber6 = new QLCDNumber( this, "lCDNumber6" );
  lCDNumber6->setNumDigits( 15 );
  layout9->addWidget( lCDNumber6 );
  AnnotateLayout->addLayout( layout9 );

  layout9_2 = new QHBoxLayout( 0, 0, 6, "layout9_2");

  pushButton1 = new QPushButton( this, "pushButton1" );
  pushButton1->setAutoDefault( FALSE );
  layout9_2->addWidget( pushButton1 );

  pushButton3 = new QPushButton( this, "pushButton3" );
  pushButton3->setAutoDefault( FALSE );
  layout9_2->addWidget( pushButton3 );
  AnnotateLayout->addLayout( layout9_2 );
  languageChange();
  resize( QSize(811, 621).expandedTo(minimumSizeHint()) );
  clearWState( WState_Polished );

  // signals and slots connections
  connect( pushButton1, SIGNAL( clicked() ), this, SLOT( abort() ) );
  connect( pushButton3, SIGNAL( clicked() ), this, SLOT( closeWindow() ) );

}

/*
 *  Destroys the object and frees any allocated resources
 */
Annotate::~Annotate()
{
  // no need to delete child widgets, Qt does it all for us
}

/*
 *  Sets the strings of the subwidgets using the current
 *  language.
 */
void Annotate::languageChange()
{
  setCaption( tr( "Annotating Sample..." ) );
  textLabel2->setText( tr( "<b><font size=\"+1\">Annotation Progress:</font></b>" ) );
  textLabel1->setText( tr( "<font size=\"+1\"><b>Time since Start:</b></font>" ) );
  textLabel1_2->setText( tr( "<b><font size=\"+1\"># Entries in</font></b> 'annotation'" ) );
  textLabel2_2->setText( tr( "<b><font size=\"+1\"># Entries in</font></b> 'modification_combination'" ) );
  textLabel3->setText( tr( "<b><font size=\"+1\"># Entries in</font></b> 'realized_modification'" ) );
  textLabel1_3->setText( tr( "<b><font size=\"+1\"># Entries in</font></b> 'modification_combination_positionless" ) );
  textLabel2_3->setText( tr( "<b><font size=\"+1\"># Entries in</font></b> 'realized_modification_positionless" ) );
  pushButton1->setText( tr( "Abort" ) );
  pushButton3->setText( tr( "Close Window" ) );
}

