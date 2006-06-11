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
// $Id: ViewSample.C,v 1.2 2006/03/28 08:03:28 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:$
// --------------------------------------------------------------------------

#include "ViewSample.h"

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


using namespace OpenMS;


/*
 *  Constructs a ViewSample as a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
ViewSample::ViewSample( QWidget* parent, const char* name, bool modal, WFlags fl )
    : QDialog( parent, name, modal, fl )
{
  if ( !name )
    setName( "ViewSample" );

  dataBrowser1 = new QDataBrowser( this, "dataBrowser1" );
  dataBrowser1->setGeometry( QRect( 10, 10, 370, 163 ) );
  QStringList dataBrowser1_stringlist;
  dataBrowser1_stringlist << "sample_ID ASC";
  dataBrowser1->setSort( dataBrowser1_stringlist );
  dataBrowser1->setAutoEdit( FALSE );
  dataBrowser1Layout = new QGridLayout( dataBrowser1, 1, 1, 11, 6, "dataBrowser1Layout");

  layout1 = new QGridLayout( 0, 1, 1, 0, 6, "layout1");

  labelEnzyme_ID = new QLabel( dataBrowser1, "labelEnzyme_ID" );

  layout1->addWidget( labelEnzyme_ID, 0, 0 );

  QLineEditAnnotation_method = new QLineEdit( dataBrowser1, "QLineEditAnnotation_method" );

  layout1->addWidget( QLineEditAnnotation_method, 2, 1 );

  QSpinBoxEnzyme_ID = new QSpinBox( dataBrowser1, "QSpinBoxEnzyme_ID" );
  QSpinBoxEnzyme_ID->setMaxValue( 2147483647 );

  layout1->addWidget( QSpinBoxEnzyme_ID, 0, 1 );

  QSpinBoxSample_ID = new QSpinBox( dataBrowser1, "QSpinBoxSample_ID" );
  QSpinBoxSample_ID->setMaxValue( 2147483647 );

  layout1->addWidget( QSpinBoxSample_ID, 3, 1 );

  QSpinBoxProtein_modification_scenario_ID = new QSpinBox( dataBrowser1, "QSpinBoxProtein_modification_scenario_ID" );
  QSpinBoxProtein_modification_scenario_ID->setMaxValue( 2147483647 );

  layout1->addWidget( QSpinBoxProtein_modification_scenario_ID, 1, 1 );

  labelSample_ID = new QLabel( dataBrowser1, "labelSample_ID" );

  layout1->addWidget( labelSample_ID, 3, 0 );

  labelProtein_modification_scenario_ID = new QLabel( dataBrowser1, "labelProtein_modification_scenario_ID" );

  layout1->addWidget( labelProtein_modification_scenario_ID, 1, 0 );

  labelAnnotation_method = new QLabel( dataBrowser1, "labelAnnotation_method" );

  layout1->addWidget( labelAnnotation_method, 2, 0 );

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

  pushButton1 = new QPushButton( this, "pushButton1" );
  pushButton1->setGeometry( QRect( 140, 190, 116, 29 ) );

  QSqlForm* dataBrowser1Form =  new QSqlForm( this, "dataBrowser1Form" );
  dataBrowser1Form->insert( QLineEditAnnotation_method, "annotation_method" );
  dataBrowser1Form->insert( QSpinBoxEnzyme_ID, "enzyme_ID" );
  dataBrowser1Form->insert( QSpinBoxSample_ID, "sample_ID" );
  dataBrowser1Form->insert( QSpinBoxProtein_modification_scenario_ID, "protein_modification_scenario_ID" );
  dataBrowser1->setForm( dataBrowser1Form );
  languageChange();
  resize( QSize(387, 236).expandedTo(minimumSizeHint()) );
  clearWState( WState_Polished );

  // signals and slots connections
  connect( pushButton1, SIGNAL( clicked() ), this, SLOT( close() ) );
  connect( PushButtonFirst, SIGNAL( clicked() ), dataBrowser1, SLOT( first() ) );
  connect( dataBrowser1, SIGNAL( firstRecordAvailable( bool ) ), PushButtonFirst, SLOT( setEnabled(bool) ) );
  connect( PushButtonPrev, SIGNAL( clicked() ), dataBrowser1, SLOT( prev() ) );
  connect( dataBrowser1, SIGNAL( prevRecordAvailable( bool ) ), PushButtonPrev, SLOT( setEnabled(bool) ) );
  connect( PushButtonNext, SIGNAL( clicked() ), dataBrowser1, SLOT( next() ) );
  connect( dataBrowser1, SIGNAL( nextRecordAvailable( bool ) ), PushButtonNext, SLOT( setEnabled(bool) ) );
  connect( PushButtonLast, SIGNAL( clicked() ), dataBrowser1, SLOT( last() ) );
  connect( dataBrowser1, SIGNAL( lastRecordAvailable( bool ) ), PushButtonLast, SLOT( setEnabled(bool) ) );
}

/*
 *  Destroys the object and frees any allocated resources
 */
ViewSample::~ViewSample()
{
  // no need to delete child widgets, Qt does it all for us
}

/*
 *  Widget polish.  Reimplemented to handle
 *  default data browser initialization
 */
void ViewSample::polish()
{
  if ( dataBrowser1 )
    {
      if ( !dataBrowser1->sqlCursor() )
        {
          QSqlCursor* cursor = new QSqlCursor( "sample" );
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
void ViewSample::languageChange()
{
  setCaption( tr( "View of MySQL Database, Table: sample" ) );
  labelEnzyme_ID->setText( tr( "Enzyme_ID" ) );
  labelSample_ID->setText( tr( "Sample_ID" ) );
  labelProtein_modification_scenario_ID->setText( tr( "Protein_modification_scenario_ID" ) );
  labelAnnotation_method->setText( tr( "Annotation_method" ) );
  PushButtonFirst->setText( tr( "|< &First" ) );
  PushButtonPrev->setText( tr( "<< &Prev" ) );
  PushButtonNext->setText( tr( "&Next >>" ) );
  PushButtonLast->setText( tr( "&Last >|" ) );
  pushButton1->setText( tr( "Done" ) );
}

