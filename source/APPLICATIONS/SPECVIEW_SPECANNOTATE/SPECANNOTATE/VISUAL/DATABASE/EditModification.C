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
// $Id: EditModification.C,v 1.2 2006/03/28 08:03:28 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:$
// --------------------------------------------------------------------------

#include "EditModification.h"

#include <qvariant.h>
#include <qsqldatabase.h>
#include <qsqlcursor.h>
#include <qsqlform.h>
#include <qsqlrecord.h>
#include <qpushbutton.h>
#include <qdatabrowser.h>
#include <qlineedit.h>
#include <qlabel.h>
#include <qtextedit.h>
#include <qlayout.h>
#include <qtooltip.h>
#include <qwhatsthis.h>

#include <qmessagebox.h>


using namespace OpenMS;


void EditModification::help()
{
  QMessageBox::information(this, tr("Database Help: modification"),
                           tr("To get some information about the different entries of the table, just place the cursor above their names, or use the \"What's this?\" function! \n" "To add an entry into the database, first click \"Insert\", enter your data, then press \"Update\"!"), 1);
}


/*
 *  Constructs a EditModification as a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
EditModification::EditModification( QWidget* parent, const char* name, bool modal, WFlags fl )
    : QDialog( parent, name, modal, fl )
{
  if ( !name )
    setName( "EditModification" );
  EditModificationLayout = new QVBoxLayout( this, 11, 6, "EditModificationLayout");

  dataBrowser1 = new QDataBrowser( this, "dataBrowser1" );
  QStringList dataBrowser1_stringlist;
  dataBrowser1_stringlist << "modification_ID ASC";
  dataBrowser1->setSort( dataBrowser1_stringlist );
  dataBrowser1Layout = new QVBoxLayout( dataBrowser1, 11, 6, "dataBrowser1Layout");

  layout1 = new QGridLayout( 0, 1, 1, 0, 6, "layout1");

  QLineEditMinus_average_mass = new QLineEdit( dataBrowser1, "QLineEditMinus_average_mass" );
  QLineEditMinus_average_mass->setAlignment( int( QLineEdit::AlignRight ) );

  layout1->addWidget( QLineEditMinus_average_mass, 3, 3 );

  QLineEditPlus_average_mass = new QLineEdit( dataBrowser1, "QLineEditPlus_average_mass" );
  QLineEditPlus_average_mass->setAlignment( int( QLineEdit::AlignRight ) );

  layout1->addWidget( QLineEditPlus_average_mass, 2, 3 );

  QLineEditPlus_formula = new QLineEdit( dataBrowser1, "QLineEditPlus_formula" );

  layout1->addWidget( QLineEditPlus_formula, 2, 1 );

  QLineEditModification_name = new QLineEdit( dataBrowser1, "QLineEditModification_name" );

  layout1->addWidget( QLineEditModification_name, 1, 1 );

  QLineEditModification_ID = new QLineEdit( dataBrowser1, "QLineEditModification_ID" );

  layout1->addWidget( QLineEditModification_ID, 0, 1 );

  labelPlus_mono_mass = new QLabel( dataBrowser1, "labelPlus_mono_mass" );

  layout1->addWidget( labelPlus_mono_mass, 4, 0 );

  labelMinus_average_mass = new QLabel( dataBrowser1, "labelMinus_average_mass" );

  layout1->addWidget( labelMinus_average_mass, 3, 2 );

  QLineEditMinus_formula = new QLineEdit( dataBrowser1, "QLineEditMinus_formula" );

  layout1->addWidget( QLineEditMinus_formula, 3, 1 );

  labelPlus_formula = new QLabel( dataBrowser1, "labelPlus_formula" );

  layout1->addWidget( labelPlus_formula, 2, 0 );

  QLineEditModification_sites = new QLineEdit( dataBrowser1, "QLineEditModification_sites" );

  layout1->addWidget( QLineEditModification_sites, 4, 3 );

  QLineEditMinus_mono_mass = new QLineEdit( dataBrowser1, "QLineEditMinus_mono_mass" );
  QLineEditMinus_mono_mass->setAlignment( int( QLineEdit::AlignRight ) );

  layout1->addWidget( QLineEditMinus_mono_mass, 1, 3 );

  labelMinus_mono_mass = new QLabel( dataBrowser1, "labelMinus_mono_mass" );

  layout1->addWidget( labelMinus_mono_mass, 1, 2 );

  labelModification_sites = new QLabel( dataBrowser1, "labelModification_sites" );

  layout1->addWidget( labelModification_sites, 4, 2 );

  QLineEditPlus_mono_mass = new QLineEdit( dataBrowser1, "QLineEditPlus_mono_mass" );
  QLineEditPlus_mono_mass->setAlignment( int( QLineEdit::AlignRight ) );

  layout1->addWidget( QLineEditPlus_mono_mass, 4, 1 );

  labelMinus_formula = new QLabel( dataBrowser1, "labelMinus_formula" );

  layout1->addWidget( labelMinus_formula, 3, 0 );

  labelPlus_average_mass = new QLabel( dataBrowser1, "labelPlus_average_mass" );

  layout1->addWidget( labelPlus_average_mass, 2, 2 );

  labelModification_name = new QLabel( dataBrowser1, "labelModification_name" );

  layout1->addWidget( labelModification_name, 1, 0 );

  labelModification_ID = new QLabel( dataBrowser1, "labelModification_ID" );

  layout1->addWidget( labelModification_ID, 0, 0 );
  dataBrowser1Layout->addLayout( layout1 );

  layout5 = new QVBoxLayout( 0, 0, 6, "layout5");

  textLabel1 = new QLabel( dataBrowser1, "textLabel1" );
  layout5->addWidget( textLabel1 );

  textEdit1 = new QTextEdit( dataBrowser1, "textEdit1" );
  layout5->addWidget( textEdit1 );
  dataBrowser1Layout->addLayout( layout5 );

  layout2 = new QHBoxLayout( 0, 0, 6, "layout2");

  PushButtonFirst = new QPushButton( dataBrowser1, "PushButtonFirst" );
  layout2->addWidget( PushButtonFirst );

  PushButtonPrev = new QPushButton( dataBrowser1, "PushButtonPrev" );
  layout2->addWidget( PushButtonPrev );

  PushButtonNext = new QPushButton( dataBrowser1, "PushButtonNext" );
  layout2->addWidget( PushButtonNext );

  PushButtonLast = new QPushButton( dataBrowser1, "PushButtonLast" );
  layout2->addWidget( PushButtonLast );
  dataBrowser1Layout->addLayout( layout2 );

  layout3 = new QHBoxLayout( 0, 0, 6, "layout3");

  PushButtonInsert = new QPushButton( dataBrowser1, "PushButtonInsert" );
  layout3->addWidget( PushButtonInsert );

  PushButtonUpdate = new QPushButton( dataBrowser1, "PushButtonUpdate" );
  layout3->addWidget( PushButtonUpdate );

  PushButtonDelete = new QPushButton( dataBrowser1, "PushButtonDelete" );
  layout3->addWidget( PushButtonDelete );
  dataBrowser1Layout->addLayout( layout3 );
  EditModificationLayout->addWidget( dataBrowser1 );

  layout4 = new QHBoxLayout( 0, 0, 6, "layout4");

  Help = new QPushButton( this, "Help" );
  layout4->addWidget( Help );
  QSpacerItem* spacer = new QSpacerItem( 400, 20, QSizePolicy::Expanding, QSizePolicy::Minimum );
  layout4->addItem( spacer );

  Done = new QPushButton( this, "Done" );
  layout4->addWidget( Done );
  EditModificationLayout->addLayout( layout4 );

  QSqlForm* dataBrowser1Form =  new QSqlForm( this, "dataBrowser1Form" );
  dataBrowser1Form->insert( QLineEditMinus_average_mass, "minus_average_mass" );
  dataBrowser1Form->insert( QLineEditPlus_average_mass, "plus_average_mass" );
  dataBrowser1Form->insert( QLineEditPlus_formula, "plus_formula" );
  dataBrowser1Form->insert( QLineEditModification_name, "modification_name" );
  dataBrowser1Form->insert( QLineEditModification_ID, "modification_ID" );
  dataBrowser1Form->insert( QLineEditMinus_formula, "minus_formula" );
  dataBrowser1Form->insert( QLineEditModification_sites, "modification_sites" );
  dataBrowser1Form->insert( QLineEditMinus_mono_mass, "minus_mono_mass" );
  dataBrowser1Form->insert( QLineEditPlus_mono_mass, "plus_mono_mass" );
  dataBrowser1Form->insert( textEdit1, "note" );
  dataBrowser1->setForm( dataBrowser1Form );
  languageChange();
  resize( QSize(634, 466).expandedTo(minimumSizeHint()) );
  clearWState( WState_Polished );

  // signals and slots connections
  connect( Done, SIGNAL( clicked() ), this, SLOT( close() ) );
  connect( PushButtonDelete, SIGNAL( clicked() ), dataBrowser1, SLOT( del() ) );
  connect( PushButtonFirst, SIGNAL( clicked() ), dataBrowser1, SLOT( first() ) );
  connect( Help, SIGNAL( clicked() ), this, SLOT( help() ) );
  connect( PushButtonInsert, SIGNAL( clicked() ), dataBrowser1, SLOT( insert() ) );
  connect( PushButtonLast, SIGNAL( clicked() ), dataBrowser1, SLOT( last() ) );
  connect( PushButtonNext, SIGNAL( clicked() ), dataBrowser1, SLOT( next() ) );
  connect( PushButtonPrev, SIGNAL( clicked() ), dataBrowser1, SLOT( prev() ) );
  connect( PushButtonUpdate, SIGNAL( clicked() ), dataBrowser1, SLOT( update() ) );

  // tab order
  setTabOrder( QLineEditModification_name, QLineEditPlus_formula );
  setTabOrder( QLineEditPlus_formula, QLineEditMinus_formula );
  setTabOrder( QLineEditMinus_formula, QLineEditPlus_mono_mass );
  setTabOrder( QLineEditPlus_mono_mass, QLineEditMinus_mono_mass );
  setTabOrder( QLineEditMinus_mono_mass, QLineEditPlus_average_mass );
  setTabOrder( QLineEditPlus_average_mass, QLineEditMinus_average_mass );
  setTabOrder( QLineEditMinus_average_mass, QLineEditModification_sites );
  setTabOrder( QLineEditModification_sites, PushButtonFirst );
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
EditModification::~EditModification()
{
  // no need to delete child widgets, Qt does it all for us
}

/*
 *  Widget polish.  Reimplemented to handle
 *  default data browser initialization
 */
void EditModification::polish()
{
  if ( dataBrowser1 )
    {
      if ( !dataBrowser1->sqlCursor() )
        {
          QSqlCursor* cursor = new QSqlCursor( "modification" );
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
void EditModification::languageChange()
{
  setCaption( tr( "Connection to MySQL-Database, Table: modification" ) );
  labelPlus_mono_mass->setText( tr( "Plus-Mass (monoisotopic)" ) );
  QToolTip::add
    ( labelPlus_mono_mass, tr( "This is the net monoisotopic molecular weight that is added to the molecule if one residue is modified with this modification" ) );
  QWhatsThis::add
    ( labelPlus_mono_mass, tr( "This is the net monoisotopic molecular weight that is added to the molecule if one residue is modified with this modification" ) );
  labelMinus_average_mass->setText( tr( "Minus-Mass (average)" ) );
  QToolTip::add
    ( labelMinus_average_mass, tr( "This is the net average molecular weight that is subtracted from the molecule if one residue is modified with this modification" ) );
  QWhatsThis::add
    ( labelMinus_average_mass, tr( "This is the net average molecular weight that is subtracted from the molecule if one residue is modified with this modification" ) );
  labelPlus_formula->setText( tr( "Plus-Formula" ) );
  QToolTip::add
    ( labelPlus_formula, tr( "This is the net molecular formula that is added to the molecule if one residue is modified with this modification" ) );
  QWhatsThis::add
    ( labelPlus_formula, tr( "This is the net molecular formula that is added to the molecule if one residue is modified with this modification" ) );
  labelMinus_mono_mass->setText( tr( "Minus-Mass (monoisotopic)" ) );
  QToolTip::add
    ( labelMinus_mono_mass, tr( "This is the net monoisotopic molecular weight that is subtracted from the molecule if one residue is modified with this modification" ) );
  QWhatsThis::add
    ( labelMinus_mono_mass, tr( "This is the net monoisotopic molecular weight that is subtracted from the molecule if one residue is modified with this modification" ) );
  labelModification_sites->setText( tr( "Modification Sites" ) );
  QToolTip::add
    ( labelModification_sites, tr( "These are the types of residues that can be modified by this modification. They are given as one-letter-codes consecutively (e.g. \"MA\" would mean Methionine and Alanine)" ) );
  QWhatsThis::add
    ( labelModification_sites, tr( "These are the types of residues that can be modified by this modification. They are given as one-letter-codes consecutively (e.g. \"MA\" would mean Methionine and Alanine)" ) );
  labelMinus_formula->setText( tr( "Minus-Formula" ) );
  QToolTip::add
    ( labelMinus_formula, tr( "This is the net molecular formula that is subtracted from the molecule if one residue is modified with this modification" ) );
  QWhatsThis::add
    ( labelMinus_formula, tr( "This is the net molecular formula that is subtracted from the molecule if one residue is modified with this modification" ) );
  labelPlus_average_mass->setText( tr( "Plus-Mass (average)" ) );
  QToolTip::add
    ( labelPlus_average_mass, tr( "This is the net average molecular weight that is added to the molecule if one residue is modified with this modification" ) );
  QWhatsThis::add
    ( labelPlus_average_mass, tr( "This is the net average molecular weight that is added to the molecule if one residue is modified with this modification" ) );
  labelModification_name->setText( tr( "Modification Name" ) );
  QToolTip::add
    ( labelModification_name, tr( "An arbitrary name for this modification - \"functor\"" ) );
  QWhatsThis::add
    ( labelModification_name, tr( "An arbitrary name for this modification - \"functor\"" ) );
  labelModification_ID->setText( tr( "Modification ID" ) );
  QToolTip::add
    ( labelModification_ID, tr( "Database-ID for this modification - \"functor\"" ) );
  QWhatsThis::add
    ( labelModification_ID, tr( "Database-ID for this modification - \"functor\"" ) );
  textLabel1->setText( tr( "Note" ) );
  QToolTip::add
    ( textLabel1, tr( "Here you can store whatever additional information to actual modification you want to!" ) );
  PushButtonFirst->setText( tr( "|< &First" ) );
  PushButtonFirst->setAccel( QKeySequence( tr( "Alt+F" ) ) );
  PushButtonPrev->setText( tr( "<< &Prev" ) );
  PushButtonPrev->setAccel( QKeySequence( tr( "Alt+P" ) ) );
  PushButtonNext->setText( tr( "&Next >>" ) );
  PushButtonLast->setText( tr( "&Last >|" ) );
  PushButtonInsert->setText( tr( "&Insert" ) );
  PushButtonInsert->setAccel( QKeySequence( tr( "Alt+I" ) ) );
  PushButtonUpdate->setText( tr( "&Update" ) );
  PushButtonUpdate->setAccel( QKeySequence( tr( "Alt+U" ) ) );
  PushButtonDelete->setText( tr( "&Delete" ) );
  Help->setText( tr( "Help" ) );
  Done->setText( tr( "Done" ) );
  QToolTip::add
    ( Done, tr( "close this window" ) );
}

