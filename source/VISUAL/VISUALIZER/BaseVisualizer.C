// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2005 -- Oliver Kohlbacher, Knut Reinert
//
//  this library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  this library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: stefan_heess   $
// --------------------------------------------------------------------------
#include <OpenMS/DATASTRUCTURES/String.h>

#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/DataTable.h>
#include <OpenMS/VISUAL/MSMetaDataExplorer.h>
#include <OpenMS/METADATA/Digestion.h>

//QT
#include <qlayout.h>
#include <qwidget.h>
#include <qaction.h>
#include <qcombobox.h>
#include <qfiledialog.h>
#include <qlabel.h> 
#include <qlineedit.h>
#include <qmenubar.h>
#include <qmessagebox.h>
#include <qpopupmenu.h>
#include <qsettings.h>
#include <qstatusbar.h>
#include <qapplication.h>
#include <qlistview.h>
#include <qtextedit.h>
#include <qhbox.h>
#include <qgroupbox.h>
#include <qpushbutton.h>
#include <iostream>
#include <vector>


//using namespace std;
using namespace OpenMS;
using namespace std;

//Constructor
BaseVisualizer::BaseVisualizer(bool editable, QWidget *parent, const char *name) 
	: DataTable(editable, parent, name)
{
  
}

String BaseVisualizer::getType()
{
	return type_;
}

void BaseVisualizer::store()
{

}

void BaseVisualizer::reject()
{

}

void BaseVisualizer::finishAdding_()
{
	addVSpacer();
	
	if(isEditable())
	{	
		addSeperator();
		addButton(undobutton_, "Undo");
		connect(undobutton_, SIGNAL(clicked()), this, SLOT(reject()) );
	}
}


