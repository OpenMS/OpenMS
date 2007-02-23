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

#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>

#include <QtGui/QPushButton>

using namespace std;

namespace OpenMS
{

BaseVisualizer::BaseVisualizer(bool editable, QWidget *parent) 
	: DataTable(editable, parent)
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

}
