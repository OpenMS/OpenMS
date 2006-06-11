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
// $Id: EditModificationXML.C,v 1.2 2006/03/28 08:03:30 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:$
// --------------------------------------------------------------------------

#include "EditModificationXML.h"

#include <qvariant.h>
#include <qpushbutton.h>
#include <qlineedit.h>
#include <qlabel.h>
#include <qtextedit.h>
#include <qlayout.h>
#include <qtooltip.h>
#include <qwhatsthis.h>


using namespace OpenMS;



void EditModificationXML::clear()
{
  lineEdit1->clear();
  lineEdit2->clear();
  lineEdit3->clear();
  lineEdit4->clear();
  lineEdit5->clear();
  lineEdit6->clear();
  lineEdit7->clear();
  lineEdit8->clear();
  lineEdit9->clear();
  textEdit1->clear();
}


void EditModificationXML::savenclear()
{
  param_.load(param_filename_);
  param_.setValue((std::string)("Preferences:SpecAnnotate:Modification:" + lineEdit1->text() + ":name").ascii(),
                  (std::string)(lineEdit2->text().ascii()));
  param_.setValue((std::string)("Preferences:SpecAnnotate:Modification:" + lineEdit1->text() + ":plus_formula").ascii(),
                  (std::string)(lineEdit3->text().ascii()));
  param_.setValue((std::string)("Preferences:SpecAnnotate:Modification:" + lineEdit1->text() + ":minus_formula").ascii(),
                  (std::string)(lineEdit4->text().ascii()));
  param_.setValue((std::string)("Preferences:SpecAnnotate:Modification:" + lineEdit1->text() + ":plus_mono_mass").ascii(),
                  (std::string)(lineEdit5->text().ascii()));
  param_.setValue((std::string)("Preferences:SpecAnnotate:Modification:" + lineEdit1->text() + ":minus_mono_mass").ascii(),
                  (std::string)(lineEdit6->text().ascii()));
  param_.setValue((std::string)("Preferences:SpecAnnotate:Modification:" + lineEdit1->text() + ":plus_average_mass").ascii(),
                  (std::string)(lineEdit7->text().ascii()));
  param_.setValue((std::string)("Preferences:SpecAnnotate:Modification:" + lineEdit1->text() + ":minus_average_mass").ascii(),
                  (std::string)(lineEdit8->text().ascii()));
  param_.setValue((std::string)("Preferences:SpecAnnotate:Modification:" + lineEdit1->text() + ":modification_sites").ascii(),
                  (std::string)(lineEdit9->text().ascii()));
  param_.setValue((std::string)("Preferences:SpecAnnotate:Modification:" + lineEdit1->text() + ":note").ascii(),
                  (std::string)(textEdit1->text().ascii()));
  param_.save(param_filename_);
  clear();
}


void EditModificationXML::lookup()
{
  param_.load(param_filename_);
  lineEdit2->setText((std::string)(param_.getValue((std::string)("Preferences:SpecAnnotate:Modification:"
                                   + lineEdit1->text() + ":name").ascii())));
  lineEdit3->setText((std::string)(param_.getValue((std::string)("Preferences:SpecAnnotate:Modification:"
                                   + lineEdit1->text() + ":plus_formula").ascii())));
  lineEdit4->setText((std::string)(param_.getValue((std::string)("Preferences:SpecAnnotate:Modification:"
                                   + lineEdit1->text() + ":minus_formula").ascii())));
  lineEdit5->setText((std::string)(param_.getValue((std::string)("Preferences:SpecAnnotate:Modification:"
                                   + lineEdit1->text() + ":plus_mono_mass").ascii())));
  lineEdit6->setText((std::string)(param_.getValue((std::string)("Preferences:SpecAnnotate:Modification:"
                                   + lineEdit1->text() + ":minus_mono_mass").ascii())));
  lineEdit7->setText((std::string)(param_.getValue((std::string)("Preferences:SpecAnnotate:Modification:"
                                   + lineEdit1->text() + ":plus_average_mass").ascii())));
  lineEdit8->setText((std::string)(param_.getValue((std::string)("Preferences:SpecAnnotate:Modification:"
                                   + lineEdit1->text() + ":minus_average_mass").ascii())));
  lineEdit9->setText((std::string)(param_.getValue((std::string)("Preferences:SpecAnnotate:Modification:"
                                   + lineEdit1->text() + ":modification_sites").ascii())));
  textEdit1->setText((std::string)(param_.getValue((std::string)("Preferences:SpecAnnotate:Modification:"
                                   + lineEdit1->text() + ":note").ascii())));
}


void EditModificationXML::setParamFilename( std::string filename )
{
  param_filename_ = filename;
}



/*
 *  Constructs a EditModificationXML as a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
EditModificationXML::EditModificationXML( QWidget* parent, const char* name, bool modal, WFlags fl )
    : QDialog( parent, name, modal, fl )
{
  if ( !name )
    setName( "EditModificationXML" );

  QWidget* privateLayoutWidget = new QWidget( this, "layout5" );
  privateLayoutWidget->setGeometry( QRect( 20, 10, 582, 372 ) );
  layout5 = new QVBoxLayout( privateLayoutWidget, 11, 6, "layout5");

  layout3 = new QVBoxLayout( 0, 0, 6, "layout3");

  layout2 = new QGridLayout( 0, 1, 1, 0, 6, "layout2");

  lineEdit9 = new QLineEdit( privateLayoutWidget, "lineEdit9" );

  layout2->addWidget( lineEdit9, 4, 3 );

  lineEdit5 = new QLineEdit( privateLayoutWidget, "lineEdit5" );

  layout2->addWidget( lineEdit5, 4, 1 );

  lineEdit4 = new QLineEdit( privateLayoutWidget, "lineEdit4" );

  layout2->addWidget( lineEdit4, 3, 1 );

  lineEdit6 = new QLineEdit( privateLayoutWidget, "lineEdit6" );

  layout2->addWidget( lineEdit6, 1, 3 );

  textLabel2 = new QLabel( privateLayoutWidget, "textLabel2" );

  layout2->addWidget( textLabel2, 1, 0 );

  textLabel9 = new QLabel( privateLayoutWidget, "textLabel9" );

  layout2->addWidget( textLabel9, 4, 2 );

  textLabel5 = new QLabel( privateLayoutWidget, "textLabel5" );

  layout2->addWidget( textLabel5, 4, 0 );

  lineEdit8 = new QLineEdit( privateLayoutWidget, "lineEdit8" );

  layout2->addWidget( lineEdit8, 3, 3 );

  textLabel7 = new QLabel( privateLayoutWidget, "textLabel7" );

  layout2->addWidget( textLabel7, 2, 2 );

  lineEdit2 = new QLineEdit( privateLayoutWidget, "lineEdit2" );

  layout2->addWidget( lineEdit2, 1, 1 );

  textLabel8 = new QLabel( privateLayoutWidget, "textLabel8" );

  layout2->addWidget( textLabel8, 3, 2 );

  lineEdit3 = new QLineEdit( privateLayoutWidget, "lineEdit3" );

  layout2->addWidget( lineEdit3, 2, 1 );

  textLabel4 = new QLabel( privateLayoutWidget, "textLabel4" );

  layout2->addWidget( textLabel4, 3, 0 );

  textLabel3 = new QLabel( privateLayoutWidget, "textLabel3" );

  layout2->addWidget( textLabel3, 2, 0 );

  textLabel6 = new QLabel( privateLayoutWidget, "textLabel6" );

  layout2->addWidget( textLabel6, 1, 2 );
  spacer1_2 = new QSpacerItem( 20, 30, QSizePolicy::Minimum, QSizePolicy::Expanding );
  layout2->addItem( spacer1_2, 0, 2 );

  textLabel1 = new QLabel( privateLayoutWidget, "textLabel1" );

  layout2->addWidget( textLabel1, 0, 0 );
  spacer1 = new QSpacerItem( 20, 30, QSizePolicy::Minimum, QSizePolicy::Expanding );
  layout2->addItem( spacer1, 0, 3 );

  lineEdit1 = new QLineEdit( privateLayoutWidget, "lineEdit1" );

  layout2->addWidget( lineEdit1, 0, 1 );

  lineEdit7 = new QLineEdit( privateLayoutWidget, "lineEdit7" );

  layout2->addWidget( lineEdit7, 2, 3 );
  layout3->addLayout( layout2 );

  textLabel10 = new QLabel( privateLayoutWidget, "textLabel10" );
  layout3->addWidget( textLabel10 );

  textEdit1 = new QTextEdit( privateLayoutWidget, "textEdit1" );
  layout3->addWidget( textEdit1 );
  layout5->addLayout( layout3 );

  layout4 = new QHBoxLayout( 0, 0, 6, "layout4");

  pushButton1 = new QPushButton( privateLayoutWidget, "pushButton1" );
  layout4->addWidget( pushButton1 );

  pushButton2 = new QPushButton( privateLayoutWidget, "pushButton2" );
  layout4->addWidget( pushButton2 );

  pushButton3 = new QPushButton( privateLayoutWidget, "pushButton3" );
  layout4->addWidget( pushButton3 );
  layout5->addLayout( layout4 );

  pushButton4 = new QPushButton( privateLayoutWidget, "pushButton4" );
  layout5->addWidget( pushButton4 );
  languageChange();
  resize( QSize(634, 421).expandedTo(minimumSizeHint()) );
  clearWState( WState_Polished );

  // signals and slots connections
  connect( pushButton4, SIGNAL( clicked() ), this, SLOT( close() ) );
  connect( pushButton2, SIGNAL( clicked() ), this, SLOT( clear() ) );
  connect( pushButton1, SIGNAL( clicked() ), this, SLOT( lookup() ) );
  connect( pushButton3, SIGNAL( clicked() ), this, SLOT( savenclear() ) );

}

/*
 *  Destroys the object and frees any allocated resources
 */
EditModificationXML::~EditModificationXML()
{
  // no need to delete child widgets, Qt does it all for us
}

/*
 *  Sets the strings of the subwidgets using the current
 *  language.
 */
void EditModificationXML::languageChange()
{
  setCaption( tr( "Edit information about modifications stored in XML file" ) );
  textLabel2->setText( tr( "Modification Name" ) );
  textLabel9->setText( tr( "Modification Sites" ) );
  textLabel5->setText( tr( "Plus Mass (monoisotopic)" ) );
  textLabel7->setText( tr( "Plus Mass (average)" ) );
  textLabel8->setText( tr( "Minus Mass (average)" ) );
  textLabel4->setText( tr( "Minus Formula" ) );
  textLabel3->setText( tr( "Plus Formula" ) );
  textLabel6->setText( tr( "Minus Mass (monoisotopic)" ) );
  textLabel1->setText( tr( "Mofication ID" ) );
  textLabel10->setText( tr( "Note" ) );
  pushButton1->setText( tr( "lookup" ) );
  QToolTip::add
    ( pushButton1, tr( "Push here to look up information about modification with given ID" ) );
  pushButton2->setText( tr( "clear" ) );
  pushButton3->setText( tr( "save 'n' clear" ) );
  pushButton4->setText( tr( "done" ) );
}

