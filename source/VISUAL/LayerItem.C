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

	LayerItem::LayerItem( QWidget * parent, const char * name, WFlags fl):
		LayerItemTemplate(parent,name,fl),
		activated_(false)
	{
		label->setBuddy(this);
	}
	
	LayerItem::~LayerItem()
	{
	
	}
	
	void LayerItem::mousePressEvent ( QMouseEvent* /*e*/ )
	{
		activate();
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
	
	void LayerItem::setIndex(UnsignedInt index)
	{
		index_ = index;
	}
	
	void LayerItem::changeState(bool state)
	{
		checkbox->setChecked(state);
	}
	
	void LayerItem::changeLabel(string l)
	{
		text_ = l;
		label->setText(l.c_str());
	}
	
	void LayerItem::toggled(bool state)
	{
		emit stateChanged(index_, state);
	}
	
	void LayerItem::contextMenuEvent( QContextMenuEvent* e )
	{
		e->accept();
		
		QPopupMenu* context_menu = new QPopupMenu(this);
		Q_CHECK_PTR(context_menu);
		context_menu->insertItem("Preferences",0,0);
		if (index_!=0)
		{
			context_menu->insertItem("Delete",1,1);
		}
		//cout << "1" << endl;
		int result = context_menu->exec( QCursor::pos() );
		if (result == 0)
		{
			//cout << "1.1" << endl;
			emit preferencesRequest(index_);
		}
		else if (result ==1)
		{
			//cout << "1.2" << endl;
			emit removeRequest(index_);
		}
		//cout << "2" << endl;
		delete(context_menu);
		//cout << "3" << endl;
	}
	
	UnsignedInt  LayerItem::getIndex() const
	{
		return index_;
	}
	
	String  LayerItem::getLabel() const
	{
		return label->text().ascii();
	}

} //namespace OpenMS
