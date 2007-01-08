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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/LayerItem.h>
#include <qlabel.h>
#include <qcheckbox.h>
#include <qcolor.h>
#include <qcursor.h>

#include <iostream>

using namespace std;

namespace OpenMS
{

	LayerItem::LayerItem(UnsignedInt index, const std::string& text, QWidget * parent, const char * name, WFlags fl):
		LayerItemTemplate(parent,name,fl),
		activated_(false)
	{
		index_ = index;
		text_ = text;
		label->setText(text.c_str());
		label->setBuddy(this);
	}

	void LayerItem::setIndex(UnsignedInt index)
	{
		index_ = index;
	}

	void LayerItem::mousePressEvent ( QMouseEvent* /*e*/ )
	{
		emit activated(index_);
	}
	
	void LayerItem::activate()
	{
		if (!activated_)
		{
			activated_ = true;
			setPaletteBackgroundColor(Qt::blue);
			label->setText((String("<font color=white><b>")+text_+ String("</b></font>")).c_str());
		}
	}
	
	void LayerItem::deactivate()
	{
		if (activated_)
		{
			activated_ = false;
			unsetPalette();
			label->setText(text_.c_str());
		}
	}
	
	bool LayerItem::isActivated()
	{
		return activated_;
	}
	
	void LayerItem::changeState(bool state)
	{
		checkbox->setChecked(state);
	}
	
	void LayerItem::contextMenuEvent( QContextMenuEvent* e )
	{
		e->accept();
		
		QPopupMenu* context_menu = new QPopupMenu(this);
		context_menu->insertItem("Preferences",0,0);
		if (index_!=0)
		{
			context_menu->insertItem("Delete",1,1);
		}
		
		int result = context_menu->exec( QCursor::pos() );
		if (result == 0)
		{
			cout << "TODO" << endl; //????
			emit preferencesRequest(index_);
		}
		else if (result == 1)
		{
			emit removeRequest(index_);
		}
	}
	
	UnsignedInt LayerItem::getIndex() const
	{
		return index_;
	}
	
	const String& LayerItem::getLabel() const
	{
		return text_;
	}

	void LayerItem::toggled(bool state)
	{
		emit stateChanged(index_, state);
	}


} //namespace OpenMS
