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
// $Id: ViewModificationCombination.C,v 1.2 2006/03/28 08:03:28 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:$
// --------------------------------------------------------------------------

#include "ViewModificationCombination.h"

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
 *  Constructs a ViewModificationCombination as a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
ViewModificationCombination::ViewModificationCombination( QWidget* parent, const char* name, bool modal, WFlags fl )
    : QDialog( parent, name, modal, fl )
{
  if ( !name )
    setName( "ViewModificationCombination" );

  pushButton1 = new QPushButton( this, "pushButton1" );
  pushButton1->setGeometry( QRect( 150, 580, 116, 29 ) );

  dataTable2 = new QDataTable( this, "dataTable2" );
  dataTable2->addColumn( "modification_combination_ID", tr( "Modification_combination_ID" ) );
  dataTable2->addColumn( "first_realized_modification_ID", tr( "First_realized_modification_ID" ) );
  dataTable2->addColumn( "next_modification_combination_ID", tr( "Next_modification_combination_ID" ) );
  dataTable2->setGeometry( QRect( 10, 10, 372, 550 ) );
  dataTable2->setReadOnly( TRUE );
  QStringList dataTable2_stringlist;
  dataTable2_stringlist << "modification_combination_ID ASC";
  dataTable2->setSort( dataTable2_stringlist );

  languageChange();
  resize( QSize(394, 624).expandedTo(minimumSizeHint()) );
  clearWState( WState_Polished );

  // signals and slots connections
  connect( pushButton1, SIGNAL( clicked() ), this, SLOT( close() ) );
}

/*
 *  Destroys the object and frees any allocated resources
 */
ViewModificationCombination::~ViewModificationCombination()
{
  // no need to delete child widgets, Qt does it all for us
}

/*
 *  Widget polish.  Reimplemented to handle
 *  default data table initialization
 */
void ViewModificationCombination::polish()
{
  if ( dataTable2 )
    {
      QSqlCursor* cursor = dataTable2->sqlCursor();
      if ( !cursor )
        {
          cursor = new QSqlCursor( "modification_combination" );
          if ( dataTable2->isReadOnly() )
            cursor->setMode( QSqlCursor::ReadOnly );
          dataTable2->setSqlCursor( cursor, FALSE, TRUE );
        }
      if ( !cursor->isActive() )
        dataTable2->refresh( QDataTable::RefreshAll );
    }
  QDialog::polish();
}

/*
 *  Sets the strings of the subwidgets using the current
 *  language.
 */
void ViewModificationCombination::languageChange()
{
  setCaption( tr( "View of MySQL Database, Table: modification_combination" ) );
  pushButton1->setText( tr( "Done" ) );
}

