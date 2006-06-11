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
// $Id: ViewProteinModificationScenario.C,v 1.2 2006/03/28 08:03:28 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:$
// --------------------------------------------------------------------------

#include "ViewProteinModificationScenario.h"

#include <qvariant.h>
#include <qsqldatabase.h>
#include <qsqlcursor.h>
#include <qsqlform.h>
#include <qsqlrecord.h>
#include <qpushbutton.h>
#include <qdatabrowser.h>
#include <qspinbox.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <qtextedit.h>
#include <qlayout.h>
#include <qtooltip.h>
#include <qwhatsthis.h>


using namespace OpenMS;

/*
 *  Constructs a ViewProteinModificationScenario as a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
ViewProteinModificationScenario::ViewProteinModificationScenario( QWidget* parent, const char* name, bool modal, WFlags fl )
    : QDialog( parent, name, modal, fl )
{
  if ( !name )
    setName( "ViewProteinModificationScenario" );

  pushButton1 = new QPushButton( this, "pushButton1" );
  pushButton1->setGeometry( QRect( 300, 430, 116, 29 ) );

  dataBrowser2 = new QDataBrowser( this, "dataBrowser2" );
  dataBrowser2->setGeometry( QRect( 10, 10, 680, 410 ) );
  QStringList dataBrowser2_stringlist;
  dataBrowser2_stringlist << "protein_modification_scenario_ID ASC";
  dataBrowser2->setSort( dataBrowser2_stringlist );
  dataBrowser2->setAutoEdit( FALSE );

  QSpinBoxProtein_modification_scenario_ID = new QSpinBox( dataBrowser2, "QSpinBoxProtein_modification_scenario_ID" );
  QSpinBoxProtein_modification_scenario_ID->setGeometry( QRect( 214, 11, 90, 21 ) );
  QSpinBoxProtein_modification_scenario_ID->setMaxValue( 2147483647 );

  QSpinBoxProtein_ID = new QSpinBox( dataBrowser2, "QSpinBoxProtein_ID" );
  QSpinBoxProtein_ID->setGeometry( QRect( 214, 38, 90, 21 ) );
  QSpinBoxProtein_ID->setMaxValue( 2147483647 );

  labelModification_combination_ID = new QLabel( dataBrowser2, "labelModification_combination_ID" );
  labelModification_combination_ID->setGeometry( QRect( 310, 40, 174, 21 ) );

  textLabel1 = new QLabel( dataBrowser2, "textLabel1" );
  textLabel1->setGeometry( QRect( 310, 10, 250, 20 ) );

  QSpinBoxModification_combination_ID = new QSpinBox( dataBrowser2, "QSpinBoxModification_combination_ID" );
  QSpinBoxModification_combination_ID->setGeometry( QRect( 570, 40, 90, 21 ) );
  QSpinBoxModification_combination_ID->setMaxValue( 2147483647 );

  labelProtein_modification_scenario_ID = new QLabel( dataBrowser2, "labelProtein_modification_scenario_ID" );
  labelProtein_modification_scenario_ID->setGeometry( QRect( 11, 11, 197, 21 ) );

  QWidget* privateLayoutWidget = new QWidget( dataBrowser2, "layout5" );
  privateLayoutWidget->setGeometry( QRect( 110, 341, 473, 50 ) );
  layout5 = new QHBoxLayout( privateLayoutWidget, 11, 6, "layout5");

  PushButtonFirst = new QPushButton( privateLayoutWidget, "PushButtonFirst" );
  layout5->addWidget( PushButtonFirst );

  PushButtonPrev = new QPushButton( privateLayoutWidget, "PushButtonPrev" );
  layout5->addWidget( PushButtonPrev );

  PushButtonNext = new QPushButton( privateLayoutWidget, "PushButtonNext" );
  layout5->addWidget( PushButtonNext );

  PushButtonLast = new QPushButton( privateLayoutWidget, "PushButtonLast" );
  layout5->addWidget( PushButtonLast );

  lineEdit2 = new QLineEdit( dataBrowser2, "lineEdit2" );
  lineEdit2->setGeometry( QRect( 570, 10, 90, 21 ) );

  QLineEditOverall_modifications = new QLineEdit( dataBrowser2, "QLineEditOverall_modifications" );
  QLineEditOverall_modifications->setGeometry( QRect( 150, 300, 510, 21 ) );

  labelOverall_modifications = new QLabel( dataBrowser2, "labelOverall_modifications" );
  labelOverall_modifications->setGeometry( QRect( 10, 300, 140, 21 ) );

  textEdit1 = new QTextEdit( dataBrowser2, "textEdit1" );
  textEdit1->setGeometry( QRect( 150, 120, 510, 160 ) );

  labelPartial_modifications = new QLabel( dataBrowser2, "labelPartial_modifications" );
  labelPartial_modifications->setGeometry( QRect( 10, 120, 130, 21 ) );

  labelProtein_ID = new QLabel( dataBrowser2, "labelProtein_ID" );
  labelProtein_ID->setGeometry( QRect( 11, 38, 197, 21 ) );

  textLabel2 = new QLabel( dataBrowser2, "textLabel2" );
  textLabel2->setGeometry( QRect( 10, 80, 120, 20 ) );

  lineEdit3 = new QLineEdit( dataBrowser2, "lineEdit3" );
  lineEdit3->setGeometry( QRect( 150, 80, 510, 21 ) );

  QSqlForm* dataBrowser2Form =  new QSqlForm( this, "dataBrowser2Form" );
  dataBrowser2Form->insert( QSpinBoxProtein_modification_scenario_ID, "protein_modification_scenario_ID" );
  dataBrowser2Form->insert( QSpinBoxProtein_ID, "protein_ID" );
  dataBrowser2Form->insert( QSpinBoxModification_combination_ID, "modification_combination_ID" );
  dataBrowser2Form->insert( lineEdit2, "modification_combination_positionless_ID" );
  dataBrowser2Form->insert( QLineEditOverall_modifications, "overall_modifications" );
  dataBrowser2Form->insert( textEdit1, "partial_modifications" );
  dataBrowser2Form->insert( lineEdit3, "annotation_method" );
  dataBrowser2->setForm( dataBrowser2Form );
  languageChange();
  resize( QSize(694, 473).expandedTo(minimumSizeHint()) );
  clearWState( WState_Polished );

  // signals and slots connections
  connect( pushButton1, SIGNAL( clicked() ), this, SLOT( close() ) );
  connect( PushButtonFirst, SIGNAL( clicked() ), dataBrowser2, SLOT( first() ) );
  connect( dataBrowser2, SIGNAL( firstRecordAvailable( bool ) ), PushButtonFirst, SLOT( setEnabled(bool) ) );
  connect( PushButtonPrev, SIGNAL( clicked() ), dataBrowser2, SLOT( prev() ) );
  connect( dataBrowser2, SIGNAL( prevRecordAvailable( bool ) ), PushButtonPrev, SLOT( setEnabled(bool) ) );
  connect( PushButtonNext, SIGNAL( clicked() ), dataBrowser2, SLOT( next() ) );
  connect( dataBrowser2, SIGNAL( nextRecordAvailable( bool ) ), PushButtonNext, SLOT( setEnabled(bool) ) );
  connect( PushButtonLast, SIGNAL( clicked() ), dataBrowser2, SLOT( last() ) );
  connect( dataBrowser2, SIGNAL( lastRecordAvailable( bool ) ), PushButtonLast, SLOT( setEnabled(bool) ) );
}

/*
 *  Destroys the object and frees any allocated resources
 */
ViewProteinModificationScenario::~ViewProteinModificationScenario()
{
  // no need to delete child widgets, Qt does it all for us
}

/*
 *  Widget polish.  Reimplemented to handle
 *  default data browser initialization
 */
void ViewProteinModificationScenario::polish()
{
  if ( dataBrowser2 )
    {
      if ( !dataBrowser2->sqlCursor() )
        {
          QSqlCursor* cursor = new QSqlCursor( "protein_modification_scenario" );
          dataBrowser2->setSqlCursor( cursor, TRUE );
          dataBrowser2->refresh();
          dataBrowser2->first();
        }
    }
  QDialog::polish();
}

/*
 *  Sets the strings of the subwidgets using the current
 *  language.
 */
void ViewProteinModificationScenario::languageChange()
{
  setCaption( tr( "View of MySQL Database, Table: protein_modification_scenario" ) );
  pushButton1->setText( tr( "Done" ) );
  labelModification_combination_ID->setText( tr( "Modification_combination_ID" ) );
  textLabel1->setText( tr( "Modification_combination_positionless_ID" ) );
  labelProtein_modification_scenario_ID->setText( tr( "Protein_modification_scenario_ID" ) );
  PushButtonFirst->setText( tr( "|< &First" ) );
  PushButtonPrev->setText( tr( "<< &Prev" ) );
  PushButtonPrev->setAccel( QKeySequence( tr( "Alt+P" ) ) );
  PushButtonNext->setText( tr( "&Next >>" ) );
  PushButtonLast->setText( tr( "&Last >|" ) );
  labelOverall_modifications->setText( tr( "Overall_modifications" ) );
  labelPartial_modifications->setText( tr( "Partial_modifications" ) );
  labelProtein_ID->setText( tr( "Protein_ID" ) );
  textLabel2->setText( tr( "Annotation_method" ) );
}

