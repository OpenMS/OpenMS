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

#include "EditEnzymeXML.h"

#include <qvariant.h>
#include <qpushbutton.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <qlayout.h>
#include <qtooltip.h>
#include <qwhatsthis.h>


using namespace OpenMS;


void EditEnzymeXML::clear_()
{
  lineEdit1->clear();
  lineEdit2->clear();
  lineEdit3->clear();
}


void EditEnzymeXML::savenclear()
{
  param_.load(param_filename_);
  param_.setValue((std::string)(("Preferences:SpecAnnotate:Enzyme:" + lineEdit1->text() + ":cleav_sites").ascii()),
                  (std::string)(lineEdit2->text().ascii()));
  param_.setValue((std::string)(("Preferences:SpecAnnotate:Enzyme:" + lineEdit1->text() + ":terminality").ascii()),
                  (std::string)(lineEdit3->text().ascii()));
  param_.save(param_filename_);
  clear_();
}


void EditEnzymeXML::lookup_()
{
  param_.load(param_filename_);
  lineEdit2->setText((std::string)(param_.getValue((std::string)(("Preferences:SpecAnnotate:Enzyme:"
                                   + lineEdit1->text() + ":cleav_sites").ascii()))));
  lineEdit3->setText((std::string)(param_.getValue((std::string)(("Preferences:SpecAnnotate:Enzyme:"
                                   + lineEdit1->text() + ":terminality").ascii()))));
}


void EditEnzymeXML::setParamFilename( std::string filename )
{
  param_filename_ = filename;
}



/*
 *  Constructs a EditEnzymeXML as a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
EditEnzymeXML::EditEnzymeXML( QWidget* parent, const char* name, bool modal, WFlags fl )
    : QDialog( parent, name, modal, fl )
{
  if ( !name )
    setName( "EditEnzymeXML" );
  EditEnzymeXMLLayout = new QGridLayout( this, 1, 1, 11, 6, "EditEnzymeXMLLayout");

  layout1 = new QGridLayout( 0, 1, 1, 0, 6, "layout1");

  textLabel3 = new QLabel( this, "textLabel3" );

  layout1->addWidget( textLabel3, 2, 0 );

  textLabel2 = new QLabel( this, "textLabel2" );

  layout1->addWidget( textLabel2, 1, 0 );

  lineEdit3 = new QLineEdit( this, "lineEdit3" );

  layout1->addWidget( lineEdit3, 2, 1 );

  textLabel1 = new QLabel( this, "textLabel1" );

  layout1->addWidget( textLabel1, 0, 0 );

  lineEdit2 = new QLineEdit( this, "lineEdit2" );

  layout1->addWidget( lineEdit2, 1, 1 );

  lineEdit1 = new QLineEdit( this, "lineEdit1" );

  layout1->addWidget( lineEdit1, 0, 1 );

  EditEnzymeXMLLayout->addLayout( layout1, 0, 0 );

  done = new QPushButton( this, "done" );

  EditEnzymeXMLLayout->addWidget( done, 2, 0 );

  layout2 = new QHBoxLayout( 0, 0, 6, "layout2");

  lookup = new QPushButton( this, "lookup" );
  layout2->addWidget( lookup );

  clear = new QPushButton( this, "clear" );
  layout2->addWidget( clear );

  save = new QPushButton( this, "save" );
  layout2->addWidget( save );

  EditEnzymeXMLLayout->addLayout( layout2, 1, 0 );
  languageChange();
  resize( QSize(416, 185).expandedTo(minimumSizeHint()) );
  clearWState( WState_Polished );

  // signals and slots connections
  connect( done, SIGNAL( clicked() ), this, SLOT( close() ) );
  connect( clear, SIGNAL( clicked() ), this, SLOT( clear_() ) );
  connect( lookup, SIGNAL( clicked() ), this, SLOT( lookup_() ) );
  connect( save, SIGNAL( clicked() ), this, SLOT( savenclear() ) );


  QToolTip::add
    ( lookup, tr( "Push here to look up information for an enzyme with given NAME in XML file." ) );
}

/*
 *  Destroys the object and frees any allocated resources
 */
EditEnzymeXML::~EditEnzymeXML()
{
  // no need to delete child widgets, Qt does it all for us
}

/*
 *  Sets the strings of the subwidgets using the current
 *  language.
 */
void EditEnzymeXML::languageChange()
{
  setCaption( tr( "Edit information about enzymes stored in XML file" ) );
  textLabel3->setText( tr( "terminality" ) );
  textLabel2->setText( tr( "cleavage sites" ) );
  textLabel1->setText( tr( "enzyme name" ) );
  done->setText( tr( "done" ) );
  lookup->setText( tr( "lookup" ) );
  clear->setText( tr( "clear" ) );
  save->setText( tr( "save 'n' clear" ) );
  save->setAccel( QKeySequence( QString::null ) );
}

