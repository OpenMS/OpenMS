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

#include "EditEnzyme.h"

#include <qvariant.h>
#include <qsqldatabase.h>
#include <qsqlcursor.h>
#include <qsqlform.h>
#include <qsqlrecord.h>
#include <qpushbutton.h>
#include <qdatabrowser.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <qlayout.h>
#include <qtooltip.h>
#include <qwhatsthis.h>


#include <qmessagebox.h>


using namespace OpenMS;


void EditEnzyme::help()
{
  QMessageBox::information(this, tr("Database Help: enzyme"),
                           tr("To get some information about the different entries of the table, just place the cursor above their names, or use the \"What's this?\" function! \n" "To add an entry into the database, first click \"Insert\", enter your data, then press \"Update\"!"), 1);
}



/*
 *  Constructs a EditEnzyme as a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
EditEnzyme::EditEnzyme( QWidget* parent, const char* name, bool modal, WFlags fl )
    : QDialog( parent, name, modal, fl )
{
  if ( !name )
    setName( "EditEnzyme" );

  dataBrowser1 = new QDataBrowser( this, "dataBrowser1" );
  dataBrowser1->setGeometry( QRect( 10, 10, 370, 177 ) );
  QStringList dataBrowser1_stringlist;
  dataBrowser1_stringlist << "enzyme_ID ASC";
  dataBrowser1->setSort( dataBrowser1_stringlist );
  dataBrowser1Layout = new QGridLayout( dataBrowser1, 1, 1, 11, 6, "dataBrowser1Layout");

  layout1 = new QGridLayout( 0, 1, 1, 0, 6, "layout1");

  labelTerminality = new QLabel( dataBrowser1, "labelTerminality" );

  layout1->addWidget( labelTerminality, 2, 0 );

  QLineEditTerminality = new QLineEdit( dataBrowser1, "QLineEditTerminality" );

  layout1->addWidget( QLineEditTerminality, 2, 1 );

  QLineEditEnzyme_name = new QLineEdit( dataBrowser1, "QLineEditEnzyme_name" );

  layout1->addWidget( QLineEditEnzyme_name, 0, 1 );

  labelCleavage_sites = new QLabel( dataBrowser1, "labelCleavage_sites" );

  layout1->addWidget( labelCleavage_sites, 1, 0 );

  QLineEditCleavage_sites = new QLineEdit( dataBrowser1, "QLineEditCleavage_sites" );

  layout1->addWidget( QLineEditCleavage_sites, 1, 1 );

  labelEnzyme_name = new QLabel( dataBrowser1, "labelEnzyme_name" );

  layout1->addWidget( labelEnzyme_name, 0, 0 );

  dataBrowser1Layout->addLayout( layout1, 0, 0 );

  layout2 = new QHBoxLayout( 0, 0, 6, "layout2");

  PushButtonFirst = new QPushButton( dataBrowser1, "PushButtonFirst" );
  layout2->addWidget( PushButtonFirst );

  PushButtonPrev = new QPushButton( dataBrowser1, "PushButtonPrev" );
  layout2->addWidget( PushButtonPrev );

  PushButtonNext = new QPushButton( dataBrowser1, "PushButtonNext" );
  layout2->addWidget( PushButtonNext );

  PushButtonLast = new QPushButton( dataBrowser1, "PushButtonLast" );
  layout2->addWidget( PushButtonLast );

  dataBrowser1Layout->addLayout( layout2, 1, 0 );

  layout3 = new QHBoxLayout( 0, 0, 6, "layout3");

  PushButtonInsert = new QPushButton( dataBrowser1, "PushButtonInsert" );
  layout3->addWidget( PushButtonInsert );

  PushButtonUpdate = new QPushButton( dataBrowser1, "PushButtonUpdate" );
  layout3->addWidget( PushButtonUpdate );

  PushButtonDelete = new QPushButton( dataBrowser1, "PushButtonDelete" );
  layout3->addWidget( PushButtonDelete );

  dataBrowser1Layout->addLayout( layout3, 2, 0 );

  Done = new QPushButton( this, "Done" );
  Done->setGeometry( QRect( 280, 200, 90, 28 ) );

  Help = new QPushButton( this, "Help" );
  Help->setGeometry( QRect( 20, 200, 90, 29 ) );

  QSqlForm* dataBrowser1Form =  new QSqlForm( this, "dataBrowser1Form" );
  dataBrowser1Form->insert( QLineEditTerminality, "terminality" );
  dataBrowser1Form->insert( QLineEditEnzyme_name, "enzyme_name" );
  dataBrowser1Form->insert( QLineEditCleavage_sites, "cleavage_sites" );
  dataBrowser1->setForm( dataBrowser1Form );
  languageChange();
  resize( QSize(383, 251).expandedTo(minimumSizeHint()) );
  clearWState( WState_Polished );

  // signals and slots connections
  connect( PushButtonFirst, SIGNAL( clicked() ), dataBrowser1, SLOT( first() ) );
  connect( PushButtonPrev, SIGNAL( clicked() ), dataBrowser1, SLOT( prev() ) );
  connect( PushButtonNext, SIGNAL( clicked() ), dataBrowser1, SLOT( next() ) );
  connect( PushButtonLast, SIGNAL( clicked() ), dataBrowser1, SLOT( last() ) );
  connect( PushButtonInsert, SIGNAL( clicked() ), dataBrowser1, SLOT( insert() ) );
  connect( PushButtonUpdate, SIGNAL( clicked() ), dataBrowser1, SLOT( update() ) );
  connect( PushButtonDelete, SIGNAL( clicked() ), dataBrowser1, SLOT( del() ) );
  connect( Done, SIGNAL( clicked() ), this, SLOT( close() ) );
  connect( Help, SIGNAL( clicked() ), this, SLOT( help() ) );

  // tab order
  setTabOrder( QLineEditEnzyme_name, QLineEditCleavage_sites );
  setTabOrder( QLineEditCleavage_sites, QLineEditTerminality );
  setTabOrder( QLineEditTerminality, PushButtonFirst );
  setTabOrder( PushButtonFirst, PushButtonPrev );
  setTabOrder( PushButtonPrev, PushButtonNext );
  setTabOrder( PushButtonNext, PushButtonLast );
  setTabOrder( PushButtonLast, PushButtonInsert );
  setTabOrder( PushButtonInsert, PushButtonUpdate );
  setTabOrder( PushButtonUpdate, PushButtonDelete );
  setTabOrder( PushButtonDelete, Done );
}

/*
 *  Destroys the object and frees any allocated resources
 */
EditEnzyme::~EditEnzyme()
{
  // no need to delete child widgets, Qt does it all for us
}

/*
 *  Widget polish.  Reimplemented to handle
 *  default data browser initialization
 */
void EditEnzyme::polish()
{
  if ( dataBrowser1 )
    {
      if ( !dataBrowser1->sqlCursor() )
        {
          QSqlCursor* cursor = new QSqlCursor( "enzyme" );
          dataBrowser1->setSqlCursor( cursor, TRUE );
          dataBrowser1->refresh();
          dataBrowser1->first();
        }
    }
  QDialog::polish();
}

/*
 *  Sets the strings of the subwidgets using the current
 *  language.
 */
void EditEnzyme::languageChange()
{
  setCaption( tr( "Connect to MySQL-Database, Table: enzyme" ) );
  labelTerminality->setText( tr( "Terminality" ) );
  QToolTip::add
    ( labelTerminality, tr( "Signifies whether the protease cuts C- or N-terminal" ) );
  QWhatsThis::add
    ( labelTerminality, tr( "Signifies whether the protease cuts C- or N-terminal" ) );
  labelCleavage_sites->setText( tr( "Cleavage-Sites" ) );
  QToolTip::add
    ( labelCleavage_sites, tr( "Specifies the aminoacids before or after which the protease cuts. They are given in one-letter-code consecutively (e.g. MA would mean Methionine and Alanine)" ) );
  QWhatsThis::add
    ( labelCleavage_sites, tr( "Specifies the aminoacids before or after which the protease cuts. They are given in one-letter-code consecutively (e.g. MA would mean Methionine and Alanine)" ) );
  labelEnzyme_name->setText( tr( "Enzyme Name" ) );
  QToolTip::add
    ( labelEnzyme_name, tr( "The name of the protease in question" ) );
  QWhatsThis::add
    ( labelEnzyme_name, tr( "The name of the protease in question" ) );
  PushButtonFirst->setText( tr( "|< &First" ) );
  PushButtonFirst->setAccel( QKeySequence( tr( "Alt+F" ) ) );
  PushButtonPrev->setText( tr( "<< &Prev" ) );
  PushButtonNext->setText( tr( "&Next >>" ) );
  PushButtonLast->setText( tr( "&Last >|" ) );
  PushButtonInsert->setText( tr( "&Insert" ) );
  PushButtonUpdate->setText( tr( "&Update" ) );
  PushButtonDelete->setText( tr( "&Delete" ) );
  Done->setText( tr( "Done" ) );
  Help->setText( tr( "Help" ) );
}

