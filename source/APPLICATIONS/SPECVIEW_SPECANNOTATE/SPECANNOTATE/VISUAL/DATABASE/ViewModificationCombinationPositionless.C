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
// $Id: ViewModificationCombinationPositionless.C,v 1.2 2006/03/28 08:03:28 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:$
// --------------------------------------------------------------------------

#include "ViewModificationCombinationPositionless.h"

#include <qvariant.h>
#include <qsqldatabase.h>
#include <qsqlcursor.h>
#include <qpushbutton.h>
#include <qdatatable.h>
#include <qlayout.h>
#include <qtooltip.h>
#include <qwhatsthis.h>


using namespace OpenMS;


/*
 *  Constructs a ViewModificationCombinationPositionless as a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
ViewModificationCombinationPositionless::ViewModificationCombinationPositionless( QWidget* parent, const char* name, bool modal, WFlags fl )
    : QDialog( parent, name, modal, fl )
{
  if ( !name )
    setName( "ViewModificationCombinationPositionless" );

  pushButton1 = new QPushButton( this, "pushButton1" );
  pushButton1->setGeometry( QRect( 200, 580, 116, 29 ) );

  dataTable1 = new QDataTable( this, "dataTable1" );
  dataTable1->addColumn( "first_realized_modification_positionless_ID", tr( "First_realized_modification_positionless_ID" ) );
  dataTable1->addColumn( "next_modification_combination_positionless_ID", tr( "Next_modification_combination_positionless_ID" ) );
  dataTable1->addColumn( "modification_combination_positionless_ID", tr( "Modification_combination_positionless_ID" ) );
  dataTable1->setGeometry( QRect( 10, 10, 490, 550 ) );
  dataTable1->setReadOnly( TRUE );
  QStringList dataTable1_stringlist;
  dataTable1_stringlist << "modification_combination_positionless_ID ASC";
  dataTable1->setSort( dataTable1_stringlist );

  languageChange();
  resize( QSize(510, 621).expandedTo(minimumSizeHint()) );
  clearWState( WState_Polished );

  // signals and slots connections
  connect( pushButton1, SIGNAL( clicked() ), this, SLOT( close() ) );
}

/*
 *  Destroys the object and frees any allocated resources
 */
ViewModificationCombinationPositionless::~ViewModificationCombinationPositionless()
{
  // no need to delete child widgets, Qt does it all for us
}

/*
 *  Widget polish.  Reimplemented to handle
 *  default data table initialization
 */
void ViewModificationCombinationPositionless::polish()
{
  if ( dataTable1 )
    {
      QSqlCursor* cursor = dataTable1->sqlCursor();
      if ( !cursor )
        {
          cursor = new QSqlCursor( "modification_combination_positionless" );
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
void ViewModificationCombinationPositionless::languageChange()
{
  setCaption( tr( "MySQL Database view, table: modification_combination_positionless" ) );
  pushButton1->setText( tr( "Done" ) );
}

