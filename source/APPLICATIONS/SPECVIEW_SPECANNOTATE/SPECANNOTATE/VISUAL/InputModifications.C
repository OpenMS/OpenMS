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
// $Id: InputModifications.C,v 1.3 2006/03/28 08:03:27 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:  $
// --------------------------------------------------------------------------


#include "InputModifications.h"

#include "../config_specannotate.h"

#include <qvariant.h>
#include <qsqldatabase.h>
#include <qstatusbar.h>
#include <list>
#include <qpushbutton.h>
#include <qlabel.h>
#include <qtextbrowser.h>
#include <qlistbox.h>
#include <qlayout.h>
#include <qtooltip.h>
#include <qwhatsthis.h>


using namespace OpenMS;


void InputModifications::init()
{
  QWidget* pa = this->parentWidget();
  if ((sd = dynamic_cast<SampleDialog*>(pa)))
    {
      //get settings from parent widget (rtti: does pa point to instance of SpecAnnotate?)
      QWidget* pa_pa = sd->parentWidget();
      if ((msa = dynamic_cast<SpecAnnotate*>(pa_pa)))
        {
          settings_ = msa->getSettings();
        }
      else
        {
          exit(1);
        }

      //! the default database connection, using the QTDesigner-MySQL-Driver QTDATABASEDRIVER
      QSqlDatabase *defaultDB = QSqlDatabase::addDatabase( QTDATABASEDRIVER );
      if ( ! defaultDB )
        {
          qWarning( "Failed to connect to driver" );
          ((QStatusBar*)(msa->statusBar()))->message( tr("Could not connect to Database"), 2000 );
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
          msa->statusBar()->message( tr("Could not connect to Database"), 2000 );
        }

      //positions
      QSqlQuery q1( "SELECT protein_ID FROM protein WHERE identifier = \"" + sd->getProtein() + "\";" );
      QString prot_ID;
      if ( q1.next() )
        {
          prot_ID = q1.value( 0 ).toString();
        }

      listBox1->clear();
      textLabel1->setText("Positions in Protein " + sd->getProtein());
      QString insert;
      for (int i = 0; i < sd->getProteinSize(); i++)
        {
          insert.setNum(i);
          QString aa_ID;
          QSqlQuery q2( "SELECT aminoacid_ID FROM sequence WHERE protein_ID = " + prot_ID +
                        " AND s_position = " + insert + ";" );
          if ( q2.next() )
            {
              aa_ID = q2.value( 0 ).toString();
            }
          QSqlQuery q3( "SELECT three_letter_code FROM aminoacid WHERE aminoacid_ID = " + aa_ID + ";" );
          if ( q3.next() )
            {
              insert += " (";
              insert += q3.value( 0 ).toString();
              insert += ")";
            }


          listBox1->insertItem(insert);
        }

      //modifications
      listBox1_2->clear();
      QSqlQuery query( "SELECT modification_name FROM modification ORDER BY modification_ID;" );
      while ( query.next() )
        {
          listBox1_2->insertItem( query.value( 0 ).toString());
        }

    }
  else
    {
      exit(1);
    }
}



void InputModifications::done()
{
  //FUEGE TEXT DES VIEWS INS ENSPR FELD IN QSAMPLEDIALOG EIN
  QString mod = textBrowser1->text();
  mod.truncate(mod.length()-2);
  mod += "*";
  sd->insertPartialMod(mod);

  close();
}





void InputModifications::addGroup()
{
  QString group;
  std::list<int> int_group;

  //get all modifications of current selection (group)
  for (uint i = 0; i < listBox1_2->count(); i++)
    {
      if (listBox1_2->isSelected(i))
        {
          QSqlQuery query( "SELECT modification_ID FROM modification WHERE modification_name = \"" + listBox1_2->text(i) + "\";" );
          if ( query.next() )
            {
              int_group.push_back(query.value( 0 ).toInt());
            }
        }
    }
  //to be sure same modification groups are recognized as same: sort
  int_group.sort();
  group = "( ";
  QString mod;
  for (std::list<int>::iterator it = int_group.begin(); it != int_group.end(); it++)
    {
      if (! (it == int_group.begin()))
        {
          group += " , ";
        }
      mod.setNum(*it);
      group += mod;
    }
  group += " )";

  //get all positions of current selection (group)
  QString add_string; //this string will be actually added to the overall modification string
  bool is_first = true;
  for (uint i = 0; i < listBox1->count(); i++)
    {
      if (listBox1->isSelected(i))
        {
          {
            if (!is_first)
              {
                add_string += " ; ";
              }
            add_string += (listBox1->text(i));
            //remove three-letter-codes
            add_string.remove((add_string.length() - 6), 6);
            add_string += " ";
            add_string += group;

            is_first = false;
          }
        }
    }

  add_string += " ; ";

  textBrowser1->insert(add_string);

  //clear selection
  resetSelection();
}


void InputModifications::resetSelection()
{
  listBox1->clearSelection();
  listBox1_2->clearSelection();
}


void InputModifications::resetString()
{
  textBrowser1->clear();
}


/*
 *  Constructs a InputModifications as a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
InputModifications::InputModifications( QWidget* parent, const char* name, bool modal, WFlags fl )
    : QDialog( parent, name, modal, fl )
{
  if ( !name )
    setName( "InputModifications" );
  InputModificationsLayout = new QVBoxLayout( this, 11, 6, "InputModificationsLayout");

  layout7 = new QVBoxLayout( 0, 0, 6, "layout7");

  textLabel3 = new QLabel( this, "textLabel3" );
  layout7->addWidget( textLabel3 );

  textBrowser1 = new QTextBrowser( this, "textBrowser1" );
  textBrowser1->setSizePolicy( QSizePolicy( (QSizePolicy::SizeType)7, (QSizePolicy::SizeType)7, 0, 0, textBrowser1->sizePolicy().hasHeightForWidth() ) );
  textBrowser1->setTextFormat( QTextBrowser::PlainText );
  textBrowser1->setAutoFormatting( int( QTextBrowser::AutoNone ) );
  layout7->addWidget( textBrowser1 );

  layout4 = new QHBoxLayout( 0, 0, 6, "layout4");

  textLabel1 = new QLabel( this, "textLabel1" );
  layout4->addWidget( textLabel1 );

  textLabel2 = new QLabel( this, "textLabel2" );
  layout4->addWidget( textLabel2 );
  layout7->addLayout( layout4 );
  InputModificationsLayout->addLayout( layout7 );

  layout1 = new QHBoxLayout( 0, 0, 6, "layout1");

  listBox1 = new QListBox( this, "listBox1" );
  listBox1->setSizePolicy( QSizePolicy( (QSizePolicy::SizeType)7, (QSizePolicy::SizeType)7, 0, 2, listBox1->sizePolicy().hasHeightForWidth() ) );
  listBox1->setSelectionMode( QListBox::Multi );
  layout1->addWidget( listBox1 );

  listBox1_2 = new QListBox( this, "listBox1_2" );
  listBox1_2->setSizePolicy( QSizePolicy( (QSizePolicy::SizeType)7, (QSizePolicy::SizeType)7, 0, 2, listBox1_2->sizePolicy().hasHeightForWidth() ) );
  listBox1_2->setSelectionMode( QListBox::Multi );
  layout1->addWidget( listBox1_2 );
  InputModificationsLayout->addLayout( layout1 );

  pushButton1 = new QPushButton( this, "pushButton1" );
  pushButton1->setAutoDefault( FALSE );
  InputModificationsLayout->addWidget( pushButton1 );

  layout5 = new QGridLayout( 0, 1, 1, 0, 6, "layout5");

  pushButton1_2 = new QPushButton( this, "pushButton1_2" );
  pushButton1_2->setAutoDefault( FALSE );

  layout5->addMultiCellWidget( pushButton1_2, 0, 0, 0, 1 );

  pushButton4 = new QPushButton( this, "pushButton4" );

  layout5->addMultiCellWidget( pushButton4, 0, 0, 2, 3 );

  pushButton5 = new QPushButton( this, "pushButton5" );

  layout5->addWidget( pushButton5, 1, 0 );

  pushButton3 = new QPushButton( this, "pushButton3" );
  QFont pushButton3_font(  pushButton3->font() );
  pushButton3_font.setBold( TRUE );
  pushButton3->setFont( pushButton3_font );
  pushButton3->setAutoDefault( FALSE );
  pushButton3->setDefault( TRUE );

  layout5->addWidget( pushButton3, 1, 3 );
  QSpacerItem* spacer = new QSpacerItem( 320, 20, QSizePolicy::Expanding, QSizePolicy::Minimum );
  layout5->addMultiCell( spacer, 1, 1, 1, 2 );
  InputModificationsLayout->addLayout( layout5 );
  languageChange();
  resize( QSize(629, 667).expandedTo(minimumSizeHint()) );
  clearWState( WState_Polished );

  // signals and slots connections
  connect( pushButton3, SIGNAL( clicked() ), this, SLOT( done() ) );
  connect( pushButton1, SIGNAL( clicked() ), this, SLOT( addGroup() ) );
  connect( pushButton1_2, SIGNAL( clicked() ), this, SLOT( resetSelection() ) );
  connect( pushButton4, SIGNAL( clicked() ), this, SLOT( resetString() ) );
  connect( pushButton5, SIGNAL( clicked() ), this, SLOT( close() ) );
  init();
}

/*
 *  Destroys the object and frees any allocated resources
 */
InputModifications::~InputModifications()
{
  // no need to delete child widgets, Qt does it all for us
}

/*
 *  Sets the strings of the subwidgets using the current
 *  language.
 */
void InputModifications::languageChange()
{
  setCaption( tr( "Partial Modification Input" ) );
  textLabel3->setText( tr( "Partial Modification String" ) );
  QToolTip::add
    ( textLabel3, tr( "Here the string signifying the selected partial modifications is displayed" ) );
  QWhatsThis::add
    ( textLabel3, tr( "Here the string signifying the selected partial modifications is displayed" ) );
  QToolTip::add
    ( textBrowser1, tr( "Here the string signifying the selected partial modifications is displayed" ) );
  QWhatsThis::add
    ( textBrowser1, tr( "Here the string signifying the selected partial modifications is displayed" ) );
  textLabel1->setText( tr( "Positions" ) );
  QToolTip::add
    ( textLabel1, tr( "Here the user can select the positions that possibly can be modified by modifications in actual modification group" ) );
  QWhatsThis::add
    ( textLabel1, tr( "Here the user can select the positions that possibly can be modified by modifications in actual modification group" ) );
  textLabel2->setText( tr( "Modifications" ) );
  QToolTip::add
    ( textLabel2, tr( "Here the user can select the modifications that all possibly can be realized at each of the positions of actual modification group" ) );
  QWhatsThis::add
    ( textLabel2, tr( "Here the user can select the modifications that all possibly can be realized at each of the positions of actual modification group" ) );
  listBox1->clear();
  listBox1->insertItem( tr( "New Item" ) );
  QToolTip::add
    ( listBox1, tr( "Here the user can select the positions that possibly can be modified by modifications in actual modification group" ) );
  QWhatsThis::add
    ( listBox1, tr( "Here the user can select the positions that possibly can be modified by modifications in actual modification group" ) );
  listBox1_2->clear();
  listBox1_2->insertItem( tr( "New Item" ) );
  QToolTip::add
    ( listBox1_2, tr( "Here the user can select the modifications that all possibly can be realized at each of the positions of actual modification group" ) );
  QWhatsThis::add
    ( listBox1_2, tr( "Here the user can select the modifications that all possibly can be realized at each of the positions of actual modification group" ) );
  pushButton1->setText( tr( "Add Group" ) );
  QToolTip::add
    ( pushButton1, tr( "This button builds a modification group out of actually selected positions and modifications and adds it in correct format to modification string displayed above" ) );
  QWhatsThis::add
    ( pushButton1, tr( "This button builds a modification group out of actually selected positions and modifications and adds it in correct format to modification string displayed above" ) );
  pushButton1_2->setText( tr( "Reset Selections" ) );
  QToolTip::add
    ( pushButton1_2, tr( "This button clears selections made in fields \"Positions\" and \"Modifications\"" ) );
  QWhatsThis::add
    ( pushButton1_2, tr( "This button clears selections made in fields \"Positions\" and \"Modifications\"" ) );
  pushButton4->setText( tr( "Reset Partial Modification String" ) );
  QToolTip::add
    ( pushButton4, tr( "This button erases the modification string built so far" ) );
  QWhatsThis::add
    ( pushButton4, tr( "This button erases the modification string built so far" ) );
  pushButton5->setText( tr( "Cancel" ) );
  QToolTip::add
    ( pushButton5, tr( "This button closes the partial modification input dialog without altering entry in sample dialog" ) );
  QWhatsThis::add
    ( pushButton5, tr( "This button closes the partial modification input dialog without altering entry in sample dialog" ) );
  pushButton3->setText( tr( "Done" ) );
  QToolTip::add
    ( pushButton3, tr( "This button adds the end sign * to partial modification string displayed above, adds it to sample dialog and closes this window" ) );
  QWhatsThis::add
    ( pushButton3, tr( "This button adds the end sign * to partial modification string displayed above, adds it to sample dialog and closes this window" ) );
}

