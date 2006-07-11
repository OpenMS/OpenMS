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

#include "ViewRealizedModificationPositionless.h"

#include <qvariant.h>
#include <qsqldatabase.h>
#include <qsqlcursor.h>
#include <qdatatable.h>
#include <qpushbutton.h>
#include <qlayout.h>
#include <qtooltip.h>
#include <qwhatsthis.h>


using namespace OpenMS;

/*
 *  Constructs a ViewRealizedModificationPositionless as a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
ViewRealizedModificationPositionless::ViewRealizedModificationPositionless( QWidget* parent, const char* name, bool modal, WFlags fl )
    : QDialog( parent, name, modal, fl )
{
  if ( !name )
    setName( "ViewRealizedModificationPositionless" );

  dataTable1 = new QDataTable( this, "dataTable1" );
  dataTable1->addColumn( "modification_ID", tr( "Modification_ID" ) );
  dataTable1->addColumn( "no_of_occurrences", tr( "No_of_occurrences" ) );
  dataTable1->addColumn( "next_realized_modification_positionless_ID", tr( "Next_realized_modification_positionless_ID" ) );
  dataTable1->addColumn( "realized_modification_positionless_ID", tr( "Realized_modification_positionless_ID" ) );
  dataTable1->setGeometry( QRect( 10, 10, 460, 550 ) );
  dataTable1->setReadOnly( TRUE );
  dataTable1->setSorting( TRUE );
  QStringList dataTable1_stringlist;
  dataTable1_stringlist << "realized_modification_positionless_ID ASC";
  dataTable1->setSort( dataTable1_stringlist );

  pushButton1 = new QPushButton( this, "pushButton1" );
  pushButton1->setGeometry( QRect( 180, 580, 116, 29 ) );

  languageChange();
  resize( QSize(480, 630).expandedTo(minimumSizeHint()) );
  clearWState( WState_Polished );

  // signals and slots connections
  connect( pushButton1, SIGNAL( clicked() ), this, SLOT( close() ) );
}

/*
 *  Destroys the object and frees any allocated resources
 */
ViewRealizedModificationPositionless::~ViewRealizedModificationPositionless()
{
  // no need to delete child widgets, Qt does it all for us
}

/*
 *  Widget polish.  Reimplemented to handle
 *  default data table initialization
 */
void ViewRealizedModificationPositionless::polish()
{
  if ( dataTable1 )
    {
      QSqlCursor* cursor = dataTable1->sqlCursor();
      if ( !cursor )
        {
          cursor = new QSqlCursor( "realized_modification_positionless" );
          if ( dataTable1->isReadOnly() )
            cursor->setMode( QSqlCursor::ReadOnly );
          dataTable1->setSqlCursor( cursor, FALSE, TRUE );
        }
      if ( !cursor->isActive() )
        dataTable1->refresh( QDataTable::RefreshAll );
    }
  QDialog::polish();
}

/*
 *  Sets the strings of the subwidgets using the current
 *  language.
 */
void ViewRealizedModificationPositionless::languageChange()
{
  setCaption( tr( "MySQL Database view, table: realized_modification_positionless" ) );
  pushButton1->setText( tr( "Done" ) );
}

