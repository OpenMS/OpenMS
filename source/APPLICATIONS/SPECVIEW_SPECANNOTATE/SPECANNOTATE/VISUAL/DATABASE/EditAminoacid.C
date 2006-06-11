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
// $Id: EditAminoacid.C,v 1.3 2006/05/30 15:46:40 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:$
// --------------------------------------------------------------------------

#include "EditAminoacid.h"

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


#include <OpenMS/FORMAT/Param.h>
#include <string>


using namespace OpenMS;


void EditAminoacid::help()
{
  QMessageBox::information(this, tr("Database Help: aminoacid"),
                           tr("Information of different positions of amino acids are stored in database: \n \"Middle\": amino acid is in a chain between two other amino acids \n \"N-term.\": amino acid is at N-terminal position \n \"C-term.\": amino acid is on C-terminal position \n \"Single\": amino acid is in no chain at all"), 1);

}



void EditAminoacid::insertIntoXML()
{
  std::string filename = "TOPPView.ini";

  Param param;
  param.load(filename);
  param.setValue((std::string)(("Preferences:SpecAnnotate:Aminoacid:" + QLineEditAminoacid_name->text() + ":one_letter_code").ascii()),
                 (std::string)(QLineEditOne_letter_code->text().ascii()));
  param.setValue((std::string)(("Preferences:SpecAnnotate:Aminoacid:" + QLineEditAminoacid_name->text() + ":three_letter_code").ascii()),
                 (std::string)(QLineEditThree_letter_code->text().ascii()));
  param.setValue((std::string)(("Preferences:SpecAnnotate:Aminoacid:" + QLineEditAminoacid_name->text() + ":middle_formula").ascii()),
                 (std::string)(QLineEditMiddle_formula->text().ascii()));
  param.setValue((std::string)(("Preferences:SpecAnnotate:Aminoacid:" + QLineEditAminoacid_name->text() + ":single_formula").ascii()),
                 (std::string)(QLineEditSingle_formula->text().ascii()));
  param.setValue((std::string)(("Preferences:SpecAnnotate:Aminoacid:" + QLineEditAminoacid_name->text() + ":c_term_formula").ascii()),
                 (std::string)(QLineEditC_term_formula->text().ascii()));
  param.setValue((std::string)(("Preferences:SpecAnnotate:Aminoacid:" + QLineEditAminoacid_name->text() + ":n_term_formula").ascii()),
                 (std::string)(QLineEditN_term_formula->text().ascii()));
  param.setValue((std::string)(("Preferences:SpecAnnotate:Aminoacid:" + QLineEditAminoacid_name->text() + ":single_mono_mass").ascii()),
                 (std::string)(QLineEditSingle_mono_mass->text().ascii()));
  param.setValue((std::string)(("Preferences:SpecAnnotate:Aminoacid:" + QLineEditAminoacid_name->text() + ":middle_mono_mass").ascii()),
                 (std::string)(QLineEditMiddle_mono_mass->text().ascii()));
  param.setValue((std::string)(("Preferences:SpecAnnotate:Aminoacid:" + QLineEditAminoacid_name->text() + ":c_term_mono_mass").ascii()),
                 (std::string)(QLineEditC_term_mono_mass->text().ascii()));
  param.setValue((std::string)(("Preferences:SpecAnnotate:Aminoacid:" + QLineEditAminoacid_name->text() + ":n_term_mono_mass").ascii()),
                 (std::string)(QLineEditN_term_mono_mass->text().ascii()));
  param.setValue((std::string)(("Preferences:SpecAnnotate:Aminoacid:" + QLineEditAminoacid_name->text() + ":single_average_mass").ascii()),
                 (std::string)(QLineEditSingle_average_mass->text().ascii()));
  param.setValue((std::string)(("Preferences:SpecAnnotate:Aminoacid:" + QLineEditAminoacid_name->text() + ":middle_average_mass").ascii()),
                 (std::string)(QLineEditMiddle_average_mass->text().ascii()));
  param.setValue((std::string)(("Preferences:SpecAnnotate:Aminoacid:" + QLineEditAminoacid_name->text() + ":c_term_average_mass").ascii()),
                 (std::string)(QLineEditC_term_average_mass->text().ascii()));
  param.setValue((std::string)(("Preferences:SpecAnnotate:Aminoacid:" + QLineEditAminoacid_name->text() + ":n_term_average_mass").ascii()),
                 (std::string)(QLineEditN_term_average_mass->text().ascii()));
  param.save(filename);


  dataBrowser3->next();


}



/*
 *  Constructs a EditAminoacid as a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
EditAminoacid::EditAminoacid( QWidget* parent, const char* name, bool modal, WFlags fl )
    : QDialog( parent, name, modal, fl )
{
  if ( !name )
    setName( "EditAminoacid" );

  dataBrowser3 = new QDataBrowser( this, "dataBrowser3" );
  dataBrowser3->setGeometry( QRect( 10, 10, 768, 235 ) );
  QStringList dataBrowser3_stringlist;
  dataBrowser3_stringlist << "aminoacid_ID ASC";
  dataBrowser3->setSort( dataBrowser3_stringlist );
  dataBrowser3Layout = new QGridLayout( dataBrowser3, 1, 1, 11, 6, "dataBrowser3Layout");

  layout7 = new QGridLayout( 0, 1, 1, 0, 6, "layout7");

  labelOne_letter_code = new QLabel( dataBrowser3, "labelOne_letter_code" );

  layout7->addWidget( labelOne_letter_code, 1, 0 );

  labelAminoacid_name = new QLabel( dataBrowser3, "labelAminoacid_name" );

  layout7->addWidget( labelAminoacid_name, 0, 0 );

  QLineEditC_term_average_mass = new QLineEdit( dataBrowser3, "QLineEditC_term_average_mass" );
  QLineEditC_term_average_mass->setAlignment( int( QLineEdit::AlignRight ) );

  layout7->addWidget( QLineEditC_term_average_mass, 3, 5 );

  labelMiddle_formula = new QLabel( dataBrowser3, "labelMiddle_formula" );

  layout7->addWidget( labelMiddle_formula, 3, 0 );

  QLineEditThree_letter_code = new QLineEdit( dataBrowser3, "QLineEditThree_letter_code" );

  layout7->addWidget( QLineEditThree_letter_code, 2, 1 );

  QLineEditAminoacid_name = new QLineEdit( dataBrowser3, "QLineEditAminoacid_name" );

  layout7->addWidget( QLineEditAminoacid_name, 0, 1 );

  labelSingle_average_mass = new QLabel( dataBrowser3, "labelSingle_average_mass" );

  layout7->addWidget( labelSingle_average_mass, 1, 4 );

  labelC_term_average_mass = new QLabel( dataBrowser3, "labelC_term_average_mass" );

  layout7->addWidget( labelC_term_average_mass, 3, 4 );

  labelC_term_formula = new QLabel( dataBrowser3, "labelC_term_formula" );

  layout7->addWidget( labelC_term_formula, 0, 2 );

  QLineEditMiddle_average_mass = new QLineEdit( dataBrowser3, "QLineEditMiddle_average_mass" );
  QLineEditMiddle_average_mass->setAlignment( int( QLineEdit::AlignRight ) );

  layout7->addWidget( QLineEditMiddle_average_mass, 2, 5 );

  QLineEditMiddle_mono_mass = new QLineEdit( dataBrowser3, "QLineEditMiddle_mono_mass" );
  QLineEditMiddle_mono_mass->setAlignment( int( QLineEdit::AlignRight ) );

  layout7->addWidget( QLineEditMiddle_mono_mass, 3, 3 );

  labelMiddle_mono_mass = new QLabel( dataBrowser3, "labelMiddle_mono_mass" );

  layout7->addWidget( labelMiddle_mono_mass, 3, 2 );

  QLineEditMiddle_formula = new QLineEdit( dataBrowser3, "QLineEditMiddle_formula" );

  layout7->addWidget( QLineEditMiddle_formula, 3, 1 );

  labelN_term_formula = new QLabel( dataBrowser3, "labelN_term_formula" );

  layout7->addWidget( labelN_term_formula, 1, 2 );

  QLineEditC_term_formula = new QLineEdit( dataBrowser3, "QLineEditC_term_formula" );

  layout7->addWidget( QLineEditC_term_formula, 0, 3 );

  QLineEditOne_letter_code = new QLineEdit( dataBrowser3, "QLineEditOne_letter_code" );

  layout7->addWidget( QLineEditOne_letter_code, 1, 1 );

  QLineEditSingle_average_mass = new QLineEdit( dataBrowser3, "QLineEditSingle_average_mass" );
  QLineEditSingle_average_mass->setAlignment( int( QLineEdit::AlignRight ) );

  layout7->addWidget( QLineEditSingle_average_mass, 1, 5 );

  QLineEditSingle_mono_mass = new QLineEdit( dataBrowser3, "QLineEditSingle_mono_mass" );
  QLineEditSingle_mono_mass->setAlignment( int( QLineEdit::AlignRight ) );

  layout7->addWidget( QLineEditSingle_mono_mass, 2, 3 );

  labelSingle_mono_mass = new QLabel( dataBrowser3, "labelSingle_mono_mass" );

  layout7->addWidget( labelSingle_mono_mass, 2, 2 );

  QLineEditN_term_mono_mass = new QLineEdit( dataBrowser3, "QLineEditN_term_mono_mass" );
  QLineEditN_term_mono_mass->setAlignment( int( QLineEdit::AlignRight ) );

  layout7->addWidget( QLineEditN_term_mono_mass, 0, 5 );

  labelN_term_mono_mass = new QLabel( dataBrowser3, "labelN_term_mono_mass" );

  layout7->addWidget( labelN_term_mono_mass, 0, 4 );

  QLineEditN_term_formula = new QLineEdit( dataBrowser3, "QLineEditN_term_formula" );

  layout7->addWidget( QLineEditN_term_formula, 1, 3 );

  labelN_term_average_mass = new QLabel( dataBrowser3, "labelN_term_average_mass" );

  layout7->addWidget( labelN_term_average_mass, 4, 4 );

  QLineEditSingle_formula = new QLineEdit( dataBrowser3, "QLineEditSingle_formula" );

  layout7->addWidget( QLineEditSingle_formula, 4, 1 );

  labelMiddle_average_mass = new QLabel( dataBrowser3, "labelMiddle_average_mass" );

  layout7->addWidget( labelMiddle_average_mass, 2, 4 );

  QLineEditC_term_mono_mass = new QLineEdit( dataBrowser3, "QLineEditC_term_mono_mass" );
  QLineEditC_term_mono_mass->setAlignment( int( QLineEdit::AlignRight ) );

  layout7->addWidget( QLineEditC_term_mono_mass, 4, 3 );

  labelC_term_mono_mass = new QLabel( dataBrowser3, "labelC_term_mono_mass" );

  layout7->addWidget( labelC_term_mono_mass, 4, 2 );

  labelSingle_formula = new QLabel( dataBrowser3, "labelSingle_formula" );

  layout7->addWidget( labelSingle_formula, 4, 0 );

  QLineEditN_term_average_mass = new QLineEdit( dataBrowser3, "QLineEditN_term_average_mass" );
  QLineEditN_term_average_mass->setAlignment( int( QLineEdit::AlignRight ) );

  layout7->addWidget( QLineEditN_term_average_mass, 4, 5 );

  labelThree_letter_code = new QLabel( dataBrowser3, "labelThree_letter_code" );

  layout7->addWidget( labelThree_letter_code, 2, 0 );

  dataBrowser3Layout->addLayout( layout7, 0, 0 );

  layout8 = new QHBoxLayout( 0, 0, 6, "layout8");

  PushButtonFirst_2 = new QPushButton( dataBrowser3, "PushButtonFirst_2" );
  layout8->addWidget( PushButtonFirst_2 );

  PushButtonPrev_2 = new QPushButton( dataBrowser3, "PushButtonPrev_2" );
  layout8->addWidget( PushButtonPrev_2 );

  PushButtonNext_2 = new QPushButton( dataBrowser3, "PushButtonNext_2" );
  layout8->addWidget( PushButtonNext_2 );

  PushButtonLast_2 = new QPushButton( dataBrowser3, "PushButtonLast_2" );
  layout8->addWidget( PushButtonLast_2 );

  dataBrowser3Layout->addLayout( layout8, 1, 0 );

  layout9 = new QHBoxLayout( 0, 0, 6, "layout9");

  PushButtonInsert_2 = new QPushButton( dataBrowser3, "PushButtonInsert_2" );
  layout9->addWidget( PushButtonInsert_2 );

  PushButtonUpdate_2 = new QPushButton( dataBrowser3, "PushButtonUpdate_2" );
  layout9->addWidget( PushButtonUpdate_2 );

  PushButtonDelete_2 = new QPushButton( dataBrowser3, "PushButtonDelete_2" );
  layout9->addWidget( PushButtonDelete_2 );

  dataBrowser3Layout->addLayout( layout9, 2, 0 );

  Done = new QPushButton( this, "Done" );
  Done->setGeometry( QRect( 670, 260, 97, 28 ) );

  Help = new QPushButton( this, "Help" );
  Help->setGeometry( QRect( 20, 260, 100, 29 ) );

  QSqlForm* dataBrowser3Form =  new QSqlForm( this, "dataBrowser3Form" );
  dataBrowser3Form->insert( QLineEditC_term_average_mass, "c_term_average_mass" );
  dataBrowser3Form->insert( QLineEditThree_letter_code, "three_letter_code" );
  dataBrowser3Form->insert( QLineEditAminoacid_name, "aminoacid_name" );
  dataBrowser3Form->insert( QLineEditMiddle_average_mass, "middle_average_mass" );
  dataBrowser3Form->insert( QLineEditMiddle_mono_mass, "middle_mono_mass" );
  dataBrowser3Form->insert( QLineEditMiddle_formula, "middle_formula" );
  dataBrowser3Form->insert( QLineEditC_term_formula, "c_term_formula" );
  dataBrowser3Form->insert( QLineEditOne_letter_code, "one_letter_code" );
  dataBrowser3Form->insert( QLineEditSingle_average_mass, "single_average_mass" );
  dataBrowser3Form->insert( QLineEditSingle_mono_mass, "single_mono_mass" );
  dataBrowser3Form->insert( QLineEditN_term_mono_mass, "n_term_mono_mass" );
  dataBrowser3Form->insert( QLineEditN_term_formula, "n_term_formula" );
  dataBrowser3Form->insert( QLineEditSingle_formula, "single_formula" );
  dataBrowser3Form->insert( QLineEditC_term_mono_mass, "c_term_mono_mass" );
  dataBrowser3Form->insert( QLineEditN_term_average_mass, "n_term_average_mass" );
  dataBrowser3->setForm( dataBrowser3Form );
  languageChange();
  resize( QSize(785, 313).expandedTo(minimumSizeHint()) );
  clearWState( WState_Polished );

  // signals and slots connections
  connect( Done, SIGNAL( clicked() ), this, SLOT( close() ) );
  connect( PushButtonFirst_2, SIGNAL( clicked() ), dataBrowser3, SLOT( first() ) );
  connect( PushButtonPrev_2, SIGNAL( clicked() ), dataBrowser3, SLOT( prev() ) );
  connect( PushButtonNext_2, SIGNAL( clicked() ), dataBrowser3 , SLOT( next() ) ); //, this , SLOT ( insertIntoXML() ) );
  connect( PushButtonLast_2, SIGNAL( clicked() ), dataBrowser3, SLOT( last() ) );
  connect( PushButtonInsert_2, SIGNAL( clicked() ), dataBrowser3, SLOT( insert() ) );
  connect( PushButtonUpdate_2, SIGNAL( clicked() ), dataBrowser3, SLOT( update() ) );
  connect( PushButtonDelete_2, SIGNAL( clicked() ), dataBrowser3, SLOT( del() ) );
  connect( Help, SIGNAL( clicked() ), this, SLOT( help() ) );

  // tab order
  setTabOrder( QLineEditAminoacid_name, QLineEditOne_letter_code );
  setTabOrder( QLineEditOne_letter_code, QLineEditThree_letter_code );
  setTabOrder( QLineEditThree_letter_code, QLineEditMiddle_formula );
  setTabOrder( QLineEditMiddle_formula, QLineEditSingle_formula );
  setTabOrder( QLineEditSingle_formula, QLineEditC_term_formula );
  setTabOrder( QLineEditC_term_formula, QLineEditN_term_formula );
  setTabOrder( QLineEditN_term_formula, QLineEditSingle_mono_mass );
  setTabOrder( QLineEditSingle_mono_mass, QLineEditMiddle_mono_mass );
  setTabOrder( QLineEditMiddle_mono_mass, QLineEditC_term_mono_mass );
  setTabOrder( QLineEditC_term_mono_mass, QLineEditN_term_mono_mass );
  setTabOrder( QLineEditN_term_mono_mass, QLineEditSingle_average_mass );
  setTabOrder( QLineEditSingle_average_mass, QLineEditMiddle_average_mass );
  setTabOrder( QLineEditMiddle_average_mass, QLineEditC_term_average_mass );
  setTabOrder( QLineEditC_term_average_mass, QLineEditN_term_average_mass );
  setTabOrder( QLineEditN_term_average_mass, PushButtonFirst_2 );
  setTabOrder( PushButtonFirst_2, PushButtonPrev_2 );
  setTabOrder( PushButtonPrev_2, PushButtonNext_2 );
  setTabOrder( PushButtonNext_2, PushButtonLast_2 );
  setTabOrder( PushButtonLast_2, PushButtonInsert_2 );
  setTabOrder( PushButtonInsert_2, PushButtonUpdate_2 );
  setTabOrder( PushButtonUpdate_2, PushButtonDelete_2 );
  setTabOrder( PushButtonDelete_2, Done );
}

/*
 *  Destroys the object and frees any allocated resources
 */
EditAminoacid::~EditAminoacid()
{
  // no need to delete child widgets, Qt does it all for us
}

/*
 *  Widget polish.  Reimplemented to handle
 *  default data browser initialization
 */
void EditAminoacid::polish()
{
  if ( dataBrowser3 )
    {
      if ( !dataBrowser3->sqlCursor() )
        {
          QSqlCursor* cursor = new QSqlCursor( "aminoacid" );
          dataBrowser3->setSqlCursor( cursor, TRUE );
          dataBrowser3->refresh();
          dataBrowser3->first();
        }
    }
  QDialog::polish();
}

/*
 *  Sets the strings of the subwidgets using the current
 *  language.
 */
void EditAminoacid::languageChange()
{
  setCaption( tr( "Connect to MySQL-Database, Table: aminoacid" ) );
  labelOne_letter_code->setText( tr( "One-Letter-Code" ) );
  QToolTip::add
    ( labelOne_letter_code, tr( "One-letter-code of amino acid" ) );
  QWhatsThis::add
    ( labelOne_letter_code, tr( "One-letter-code of amino acid" ) );
  labelAminoacid_name->setText( tr( "Aminoacid Name" ) );
  QToolTip::add
    ( labelAminoacid_name, tr( "Name of the amino acid" ) );
  QWhatsThis::add
    ( labelAminoacid_name, tr( "Name of the amino acid" ) );
  QToolTip::add
    ( QLineEditC_term_average_mass, tr( "Average molecular mass of corresponding formula" ) );
  QWhatsThis::add
    ( QLineEditC_term_average_mass, tr( "Average molecular mass of corresponding formula" ) );
  labelMiddle_formula->setText( tr( "\"Middle\" Formula" ) );
  QToolTip::add
    ( labelMiddle_formula, tr( "Sum formula of amino acid, if it is in a \"middle\" position, neither n-terminal nor c-terminal" ) );
  QWhatsThis::add
    ( labelMiddle_formula, tr( "Sum formula of amino acid, if it is in a \"middle\" position, neither n-terminal nor c-terminal" ) );
  QToolTip::add
    ( QLineEditThree_letter_code, tr( "Three-letter-code of amino acid" ) );
  QWhatsThis::add
    ( QLineEditThree_letter_code, tr( "Three-letter-code of amino acid" ) );
  QToolTip::add
    ( QLineEditAminoacid_name, tr( "Name of the amino acid" ) );
  QWhatsThis::add
    ( QLineEditAminoacid_name, tr( "Name of the amino acid" ) );
  labelSingle_average_mass->setText( tr( "\"Single\" average Mass" ) );
  QToolTip::add
    ( labelSingle_average_mass, tr( "Average molecular mass of corresponding formula" ) );
  QWhatsThis::add
    ( labelSingle_average_mass, tr( "Average molecular mass of corresponding formula" ) );
  labelC_term_average_mass->setText( tr( "\"C-term.\" average Mass" ) );
  QToolTip::add
    ( labelC_term_average_mass, tr( "Average molecular mass of corresponding formula" ) );
  QWhatsThis::add
    ( labelC_term_average_mass, tr( "Average molecular mass of corresponding formula" ) );
  labelC_term_formula->setText( tr( "\"C-term.\" Formula" ) );
  QToolTip::add
    ( labelC_term_formula, tr( "Sum formula of amino acid, if it is in c-terminal position" ) );
  QWhatsThis::add
    ( labelC_term_formula, tr( "Sum formula of amino acid, if it is in c-terminal position" ) );
  QToolTip::add
    ( QLineEditMiddle_average_mass, tr( "Average molecular mass of corresponding formula" ) );
  QWhatsThis::add
    ( QLineEditMiddle_average_mass, tr( "Average molecular mass of corresponding formula" ) );
  QToolTip::add
    ( QLineEditMiddle_mono_mass, tr( "Monoisotopic molecular weight of corresponding formula" ) );
  QWhatsThis::add
    ( QLineEditMiddle_mono_mass, tr( "Monoisotopic molecular weight of corresponding formula" ) );
  labelMiddle_mono_mass->setText( tr( "\"Middle\" mono. Mass" ) );
  QToolTip::add
    ( labelMiddle_mono_mass, tr( "Monoisotopic molecular weight of corresponding formula" ) );
  QWhatsThis::add
    ( labelMiddle_mono_mass, tr( "Monoisotopic molecular weight of corresponding formula" ) );
  QToolTip::add
    ( QLineEditMiddle_formula, tr( "Sum formula of amino acid, if it is in a \"middle\" position, neither n-terminal nor c-terminal" ) );
  QWhatsThis::add
    ( QLineEditMiddle_formula, tr( "Sum formula of amino acid, if it is in a \"middle\" position, neither n-terminal nor c-terminal" ) );
  labelN_term_formula->setText( tr( "\"N-term.\" Formula" ) );
  QToolTip::add
    ( labelN_term_formula, tr( "Sum formula of amino acid, if it is in n-terminal position" ) );
  QWhatsThis::add
    ( labelN_term_formula, tr( "Sum formula of amino acid, if it is in n-terminal position" ) );
  QToolTip::add
    ( QLineEditC_term_formula, tr( "Sum formula of amino acid, if it is in c-terminal position" ) );
  QWhatsThis::add
    ( QLineEditC_term_formula, tr( "Sum formula of amino acid, if it is in c-terminal position" ) );
  QToolTip::add
    ( QLineEditOne_letter_code, tr( "One-letter-code of amino acid" ) );
  QWhatsThis::add
    ( QLineEditOne_letter_code, tr( "One-letter-code of amino acid" ) );
  QToolTip::add
    ( QLineEditSingle_average_mass, tr( "Average molecular mass of corresponding formula" ) );
  QWhatsThis::add
    ( QLineEditSingle_average_mass, tr( "Average molecular mass of corresponding formula" ) );
  QToolTip::add
    ( QLineEditSingle_mono_mass, tr( "Monoisotopic molecular weight of corresponding formula" ) );
  QWhatsThis::add
    ( QLineEditSingle_mono_mass, tr( "Monoisotopic molecular weight of corresponding formula" ) );
  labelSingle_mono_mass->setText( tr( "\"Single\" mono. Mass" ) );
  QToolTip::add
    ( labelSingle_mono_mass, tr( "Monoisotopic molecular weight of corresponding formula" ) );
  QWhatsThis::add
    ( labelSingle_mono_mass, tr( "Monoisotopic molecular weight of corresponding formula" ) );
  QToolTip::add
    ( QLineEditN_term_mono_mass, tr( "Monoisotopic molecular weight of corresponding formula" ) );
  QWhatsThis::add
    ( QLineEditN_term_mono_mass, tr( "Monoisotopic molecular weight of corresponding formula" ) );
  labelN_term_mono_mass->setText( tr( "\"N-term.\" mono. Mass" ) );
  QToolTip::add
    ( labelN_term_mono_mass, tr( "Monoisotopic molecular mass of corresponding formula" ) );
  QWhatsThis::add
    ( labelN_term_mono_mass, tr( "Monoisotopic molecular mass of corresponding formula" ) );
  QToolTip::add
    ( QLineEditN_term_formula, tr( "Sum formula of amino acid, if it is in n-terminal position" ) );
  QWhatsThis::add
    ( QLineEditN_term_formula, tr( "Sum formula of amino acid, if it is in n-terminal position" ) );
  labelN_term_average_mass->setText( tr( "\"N-term.\" average Mass" ) );
  QToolTip::add
    ( labelN_term_average_mass, tr( "Average molecular mass of corresponding formula" ) );
  QWhatsThis::add
    ( labelN_term_average_mass, tr( "Average molecular mass of corresponding formula" ) );
  QToolTip::add
    ( QLineEditSingle_formula, tr( "Sum formula of amino acid, if it is a single amino acid, not incorporated in any peptide or protein" ) );
  QWhatsThis::add
    ( QLineEditSingle_formula, tr( "Sum formula of amino acid, if it is a single amino acid, not incorporated in any peptide or protein" ) );
  labelMiddle_average_mass->setText( tr( "\"Middle\" average Mass" ) );
  QToolTip::add
    ( labelMiddle_average_mass, tr( "Average molecular mass of corresponding formula" ) );
  QWhatsThis::add
    ( labelMiddle_average_mass, tr( "Average molecular mass of corresponding formula" ) );
  QToolTip::add
    ( QLineEditC_term_mono_mass, tr( "Monoisotopic molecular weight of corresponding formula" ) );
  QWhatsThis::add
    ( QLineEditC_term_mono_mass, tr( "Monoisotopic molecular weight of corresponding formula" ) );
  labelC_term_mono_mass->setText( tr( "\"C-term.\" mono. Mass" ) );
  QToolTip::add
    ( labelC_term_mono_mass, tr( "Monoisotopic molecular weight of corresponding formula" ) );
  QWhatsThis::add
    ( labelC_term_mono_mass, tr( "Monoisotopic molecular weight of corresponding formula" ) );
  labelSingle_formula->setText( tr( "\"Single\" Formula" ) );
  QToolTip::add
    ( labelSingle_formula, tr( "Sum formula of amino acid, if it is a single amino acid, not incorporated in any peptide or protein" ) );
  QWhatsThis::add
    ( labelSingle_formula, tr( "Sum formula of amino acid, if it is a single amino acid, not incorporated in any peptide or protein" ) );
  QToolTip::add
    ( QLineEditN_term_average_mass, tr( "Average molecular mass of corresponding formula" ) );
  QWhatsThis::add
    ( QLineEditN_term_average_mass, tr( "Average molecular mass of corresponding formula" ) );
  labelThree_letter_code->setText( tr( "Three-Letter-Code" ) );
  QToolTip::add
    ( labelThree_letter_code, tr( "Three-letter-code of amino acid" ) );
  QWhatsThis::add
    ( labelThree_letter_code, tr( "Three-letter-code of amino acid" ) );
  PushButtonFirst_2->setText( tr( "|< &First" ) );
  PushButtonPrev_2->setText( tr( "<< &Prev" ) );
  PushButtonNext_2->setText( tr( "&Next >>" ) );
  PushButtonLast_2->setText( tr( "&Last >|" ) );
  PushButtonInsert_2->setText( tr( "&Insert" ) );
  PushButtonUpdate_2->setText( tr( "&Update" ) );
  PushButtonDelete_2->setText( tr( "&Delete" ) );
  Done->setText( tr( "Done" ) );
  Help->setText( tr( "Help" ) );
}

