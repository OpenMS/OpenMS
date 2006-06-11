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
// $Id: EditProtein.C,v 1.2 2006/03/28 08:03:28 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:$
// --------------------------------------------------------------------------

#include "EditProtein.h"

#include <qvariant.h>
#include <qsqldatabase.h>
#include <qsqlcursor.h>
#include <qsqlform.h>
#include <qsqlrecord.h>
#include <qpushbutton.h>
#include <qdatabrowser.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <qspinbox.h>
#include <qlayout.h>
#include <qtooltip.h>
#include <qwhatsthis.h>


#include <qmessagebox.h>


using namespace OpenMS;


void EditProtein::help()
{
  QMessageBox::information(this, tr("Database Help: protein"),
                           tr("To add a protein to database one must enter an unique Identifier (e.g. s. th. build from pdb entry code) \n and either a pdb-file name or a fasta-file name. Both as absolute names including correct paths.\n The rest of the fields, including protein ID, are filled automatically upon instanciation of class \"MyProtein\""), 1);

}



/*
 *  Constructs a EditProtein as a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
EditProtein::EditProtein( QWidget* parent, const char* name, bool modal, WFlags fl )
    : QDialog( parent, name, modal, fl )
{
  if ( !name )
    setName( "EditProtein" );

  dataBrowser1 = new QDataBrowser( this, "dataBrowser1" );
  dataBrowser1->setGeometry( QRect( 10, 10, 930, 200 ) );
  QStringList dataBrowser1_stringlist;
  dataBrowser1_stringlist << "protein_ID ASC";
  dataBrowser1->setSort( dataBrowser1_stringlist );
  dataBrowser1Layout = new QGridLayout( dataBrowser1, 1, 1, 11, 6, "dataBrowser1Layout");

  layout2 = new QGridLayout( 0, 1, 1, 0, 6, "layout2");

  labelIdentifier = new QLabel( dataBrowser1, "labelIdentifier" );

  layout2->addWidget( labelIdentifier, 1, 0 );

  QLineEditAverage_mass = new QLineEdit( dataBrowser1, "QLineEditAverage_mass" );
  QLineEditAverage_mass->setAlignment( int( QLineEdit::AlignRight ) );

  layout2->addWidget( QLineEditAverage_mass, 2, 3 );

  QLineEditPdb_filename = new QLineEdit( dataBrowser1, "QLineEditPdb_filename" );

  layout2->addWidget( QLineEditPdb_filename, 2, 1 );

  QLineEditIdentifier = new QLineEdit( dataBrowser1, "QLineEditIdentifier" );

  layout2->addWidget( QLineEditIdentifier, 1, 1 );

  labelProtein_ID = new QLabel( dataBrowser1, "labelProtein_ID" );

  layout2->addWidget( labelProtein_ID, 0, 0 );

  labelMono_mass = new QLabel( dataBrowser1, "labelMono_mass" );

  layout2->addWidget( labelMono_mass, 1, 2 );

  labelFasta_filename = new QLabel( dataBrowser1, "labelFasta_filename" );

  layout2->addWidget( labelFasta_filename, 3, 0 );

  QSpinBoxProtein_ID = new QSpinBox( dataBrowser1, "QSpinBoxProtein_ID" );
  QSpinBoxProtein_ID->setMaxValue( 2147483647 );

  layout2->addWidget( QSpinBoxProtein_ID, 0, 1 );

  QLineEditSequence_oneletter = new QLineEdit( dataBrowser1, "QLineEditSequence_oneletter" );

  layout2->addWidget( QLineEditSequence_oneletter, 0, 3 );

  labelSequence_oneletter = new QLabel( dataBrowser1, "labelSequence_oneletter" );

  layout2->addWidget( labelSequence_oneletter, 0, 2 );

  labelPdb_filename = new QLabel( dataBrowser1, "labelPdb_filename" );

  layout2->addWidget( labelPdb_filename, 2, 0 );

  QLineEditMono_mass = new QLineEdit( dataBrowser1, "QLineEditMono_mass" );
  QLineEditMono_mass->setAlignment( int( QLineEdit::AlignRight ) );

  layout2->addWidget( QLineEditMono_mass, 1, 3 );

  labelNo_of_aminoacids = new QLabel( dataBrowser1, "labelNo_of_aminoacids" );

  layout2->addWidget( labelNo_of_aminoacids, 3, 2 );

  labelAverage_mass = new QLabel( dataBrowser1, "labelAverage_mass" );

  layout2->addWidget( labelAverage_mass, 2, 2 );

  QLineEditFasta_filename = new QLineEdit( dataBrowser1, "QLineEditFasta_filename" );

  layout2->addWidget( QLineEditFasta_filename, 3, 1 );

  QSpinBoxNo_of_aminoacids = new QSpinBox( dataBrowser1, "QSpinBoxNo_of_aminoacids" );
  QSpinBoxNo_of_aminoacids->setMaxValue( 2147483647 );

  layout2->addWidget( QSpinBoxNo_of_aminoacids, 3, 3 );

  dataBrowser1Layout->addLayout( layout2, 0, 0 );

  layout3 = new QHBoxLayout( 0, 0, 6, "layout3");

  PushButtonFirst = new QPushButton( dataBrowser1, "PushButtonFirst" );
  layout3->addWidget( PushButtonFirst );

  PushButtonPrev = new QPushButton( dataBrowser1, "PushButtonPrev" );
  layout3->addWidget( PushButtonPrev );

  PushButtonNext = new QPushButton( dataBrowser1, "PushButtonNext" );
  layout3->addWidget( PushButtonNext );

  PushButtonLast = new QPushButton( dataBrowser1, "PushButtonLast" );
  layout3->addWidget( PushButtonLast );

  dataBrowser1Layout->addLayout( layout3, 1, 0 );

  layout4 = new QHBoxLayout( 0, 0, 6, "layout4");

  PushButtonInsert = new QPushButton( dataBrowser1, "PushButtonInsert" );
  layout4->addWidget( PushButtonInsert );

  PushButtonUpdate = new QPushButton( dataBrowser1, "PushButtonUpdate" );
  layout4->addWidget( PushButtonUpdate );

  PushButtonDelete = new QPushButton( dataBrowser1, "PushButtonDelete" );
  layout4->addWidget( PushButtonDelete );

  dataBrowser1Layout->addLayout( layout4, 2, 0 );

  pushButton8 = new QPushButton( this, "pushButton8" );
  pushButton8->setGeometry( QRect( 810, 220, 116, 29 ) );

  Help = new QPushButton( this, "Help" );
  Help->setGeometry( QRect( 20, 220, 116, 29 ) );

  QSqlForm* dataBrowser1Form =  new QSqlForm( this, "dataBrowser1Form" );
  dataBrowser1Form->insert( QLineEditAverage_mass, "average_mass" );
  dataBrowser1Form->insert( QLineEditPdb_filename, "pdb_filename" );
  dataBrowser1Form->insert( QLineEditIdentifier, "identifier" );
  dataBrowser1Form->insert( QSpinBoxProtein_ID, "protein_ID" );
  dataBrowser1Form->insert( QLineEditSequence_oneletter, "sequence_oneletter" );
  dataBrowser1Form->insert( QLineEditMono_mass, "mono_mass" );
  dataBrowser1Form->insert( QLineEditFasta_filename, "fasta_filename" );
  dataBrowser1Form->insert( QSpinBoxNo_of_aminoacids, "no_of_aminoacids" );
  dataBrowser1->setForm( dataBrowser1Form );
  languageChange();
  resize( QSize(944, 270).expandedTo(minimumSizeHint()) );
  clearWState( WState_Polished );

  // signals and slots connections
  connect( PushButtonFirst, SIGNAL( clicked() ), dataBrowser1, SLOT( first() ) );
  connect( PushButtonPrev, SIGNAL( clicked() ), dataBrowser1, SLOT( prev() ) );
  connect( PushButtonNext, SIGNAL( clicked() ), dataBrowser1, SLOT( next() ) );
  connect( PushButtonLast, SIGNAL( clicked() ), dataBrowser1, SLOT( last() ) );
  connect( PushButtonInsert, SIGNAL( clicked() ), dataBrowser1, SLOT( insert() ) );
  connect( PushButtonUpdate, SIGNAL( clicked() ), dataBrowser1, SLOT( update() ) );
  connect( PushButtonDelete, SIGNAL( clicked() ), dataBrowser1, SLOT( del() ) );
  connect( pushButton8, SIGNAL( clicked() ), this, SLOT( close() ) );
  connect( Help, SIGNAL( clicked() ), this, SLOT( help() ) );
}

/*
 *  Destroys the object and frees any allocated resources
 */
EditProtein::~EditProtein()
{
  // no need to delete child widgets, Qt does it all for us
}

/*
 *  Widget polish.  Reimplemented to handle
 *  default data browser initialization
 */
void EditProtein::polish()
{
  if ( dataBrowser1 )
    {
      if ( !dataBrowser1->sqlCursor() )
        {
          QSqlCursor* cursor = new QSqlCursor( "protein" );
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
void EditProtein::languageChange()
{
  setCaption( tr( "Connect to MySQL Database, table: protein" ) );
  labelIdentifier->setText( tr( "Identifier" ) );
  QToolTip::add
    ( labelIdentifier, tr( "User defined identifier for this protein (must be unique!)" ) );
  QWhatsThis::add
    ( labelIdentifier, tr( "User defined identifier for this protein (must be unique!)" ) );
  QToolTip::add
    ( QLineEditAverage_mass, tr( "Average molecular mass of protein" ) );
  QWhatsThis::add
    ( QLineEditAverage_mass, tr( "Average molecular mass of protein" ) );
  QToolTip::add
    ( QLineEditPdb_filename, tr( "Filename (including path) of the pdb file of this protein" ) );
  QWhatsThis::add
    ( QLineEditPdb_filename, tr( "Filename (including path) of the pdb file of this protein" ) );
  QToolTip::add
    ( QLineEditIdentifier, tr( "User defined identifier for this protein (must be unique!)" ) );
  QWhatsThis::add
    ( QLineEditIdentifier, tr( "User defined identifier for this protein (must be unique!)" ) );
  labelProtein_ID->setText( tr( "Protein ID" ) );
  QToolTip::add
    ( labelProtein_ID, tr( "Unique database ID of protein" ) );
  QWhatsThis::add
    ( labelProtein_ID, tr( "Unique database ID of protein" ) );
  labelMono_mass->setText( tr( "Mono mass" ) );
  QToolTip::add
    ( labelMono_mass, tr( "Monoisotopic molecular mass of protein" ) );
  QWhatsThis::add
    ( labelMono_mass, tr( "Monoisotopic molecular mass of protein" ) );
  labelFasta_filename->setText( tr( "Fasta filename" ) );
  QToolTip::add
    ( labelFasta_filename, tr( "Filename (including path) of the fasta file of this protein (required for SpecAnnotate)" ) );
  QWhatsThis::add
    ( labelFasta_filename, tr( "Filename (including path) of the fasta file of this protein (required for SpecAnnotate)" ) );
  QToolTip::add
    ( QSpinBoxProtein_ID, tr( "Unique database ID of protein" ) );
  QWhatsThis::add
    ( QSpinBoxProtein_ID, tr( "Unique database ID of protein" ) );
  QToolTip::add
    ( QLineEditSequence_oneletter, tr( "Sequence of protein in one-letter-code" ) );
  QWhatsThis::add
    ( QLineEditSequence_oneletter, tr( "Sequence of protein in one-letter-code" ) );
  labelSequence_oneletter->setText( tr( "Sequence oneletter" ) );
  QToolTip::add
    ( labelSequence_oneletter, tr( "Sequence of protein in one-letter-code" ) );
  QWhatsThis::add
    ( labelSequence_oneletter, tr( "Sequence of protein in one-letter-code" ) );
  labelPdb_filename->setText( tr( "Pdb filename" ) );
  QToolTip::add
    ( labelPdb_filename, tr( "Filename (including path) of the pdb file of this protein" ) );
  QWhatsThis::add
    ( labelPdb_filename, tr( "Filename (including path) of the pdb file of this protein" ) );
  QToolTip::add
    ( QLineEditMono_mass, tr( "Monoisotopic molecular mass of protein" ) );
  QWhatsThis::add
    ( QLineEditMono_mass, tr( "Monoisotopic molecular mass of protein" ) );
  labelNo_of_aminoacids->setText( tr( "No of aminoacids" ) );
  QToolTip::add
    ( labelNo_of_aminoacids, tr( "Number od aminoacids (lenght) of protein" ) );
  QWhatsThis::add
    ( labelNo_of_aminoacids, tr( "Number od aminoacids (lenght) of protein" ) );
  labelAverage_mass->setText( tr( "Average mass" ) );
  QToolTip::add
    ( labelAverage_mass, tr( "Average molecular mass of protein" ) );
  QWhatsThis::add
    ( labelAverage_mass, tr( "Average molecular mass of protein" ) );
  QToolTip::add
    ( QLineEditFasta_filename, tr( "Filename (including path) of the fasta file of this protein (required for SpecAnnotate)" ) );
  QWhatsThis::add
    ( QLineEditFasta_filename, tr( "Filename (including path) of the fasta file of this protein (required for SpecAnnotate)" ) );
  QToolTip::add
    ( QSpinBoxNo_of_aminoacids, tr( "Number od aminoacids (lenght) of protein" ) );
  QWhatsThis::add
    ( QSpinBoxNo_of_aminoacids, tr( "Number od aminoacids (lenght) of protein" ) );
  PushButtonFirst->setText( tr( "|< &First" ) );
  PushButtonPrev->setText( tr( "<< &Prev" ) );
  PushButtonNext->setText( tr( "&Next >>" ) );
  PushButtonLast->setText( tr( "&Last >|" ) );
  PushButtonInsert->setText( tr( "&Insert" ) );
  PushButtonUpdate->setText( tr( "&Update" ) );
  PushButtonDelete->setText( tr( "&Delete" ) );
  pushButton8->setText( tr( "Done" ) );
  Help->setText( tr( "Help" ) );
}

