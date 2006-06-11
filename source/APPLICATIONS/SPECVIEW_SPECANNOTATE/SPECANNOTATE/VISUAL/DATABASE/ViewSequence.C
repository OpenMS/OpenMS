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
// $Id: ViewSequence.C,v 1.2 2006/03/28 08:03:28 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:$
// --------------------------------------------------------------------------

#include "ViewSequence.h"

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
 *  Constructs a ViewSequence as a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
ViewSequence::ViewSequence( QWidget* parent, const char* name, bool modal, WFlags fl )
    : QDialog( parent, name, modal, fl )
{
  if ( !name )
    setName( "ViewSequence" );

  Done = new QPushButton( this, "Done" );
  Done->setGeometry( QRect( 150, 510, 116, 29 ) );

  dataTable1 = new QDataTable( this, "dataTable1" );
  dataTable1->addColumn( "protein_ID", tr( "Protein_ID" ) );
  dataTable1->addColumn( "s_position", tr( "S_position" ) );
  dataTable1->addColumn( "aminoacid_ID", tr( "Aminoacid_ID" ) );
  dataTable1->setGeometry( QRect( 10, 10, 351, 490 ) );
  dataTable1->setReadOnly( TRUE );
  QStringList dataTable1_stringlist;
  dataTable1_stringlist << "protein_ID ASC";
  dataTable1_stringlist << "s_position ASC";
  dataTable1->setSort( dataTable1_stringlist );

  languageChange();
  resize( QSize(374, 547).expandedTo(minimumSizeHint()) );
  clearWState( WState_Polished );

  // signals and slots connections
  connect( Done, SIGNAL( clicked() ), this, SLOT( close() ) );
}

/*
 *  Destroys the object and frees any allocated resources
 */
ViewSequence::~ViewSequence()
{
  // no need to delete child widgets, Qt does it all for us
}

/*
 *  Widget polish.  Reimplemented to handle
 *  default data table initialization
 */
void ViewSequence::polish()
{
  if ( dataTable1 )
    {
      QSqlCursor* cursor = dataTable1->sqlCursor();
      if ( !cursor )
        {
          cursor = new QSqlCursor( "sequence" );
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
void ViewSequence::languageChange()
{
  setCaption( tr( "View of MySQL Database, Table: sequence" ) );
  Done->setText( tr( "Done" ) );
}

