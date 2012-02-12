// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/EnhancedTabBar.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <QtGui/QMouseEvent>
#include <QtGui/QMenu>
#include <QtGui/QMessageBox>

#include <iostream>

using namespace std;

namespace OpenMS
{

	EnhancedTabBar::EnhancedTabBar( QWidget * parent) 
		: QTabBar(parent)
	{
		connect(this,SIGNAL(currentChanged(int)),this,SLOT(currentChanged_(int)));
		
		//set up drag-and-drop
		setAcceptDrops(true);
	}
	
	EnhancedTabBar::~EnhancedTabBar()
	{
		
	}

	void EnhancedTabBar::dragEnterEvent(QDragEnterEvent* e)
	{
		e->acceptProposedAction();
	}
	
	void EnhancedTabBar::dropEvent(QDropEvent* e)
	{
		int tab = tabAt_(e->pos());
		if (tab!=-1)
		{
			emit dropOnTab(e->mimeData(), e->source(), tabData(tab).toInt());
		}
		else
		{
		  emit dropOnWidget(e->mimeData(), e->source());
		}

		e->acceptProposedAction();
	}

	void EnhancedTabBar::contextMenuEvent(QContextMenuEvent* e)
	{
		int tab = tabAt_(e->pos());
		if (tab!=-1)
		{
			QMenu menu(this);
			menu.addAction("Close");
			if (menu.exec(e->globalPos()))
			{
				emit aboutToCloseId(tabData(tab).toInt());
				removeTab(tab);
			}
		}
	}

	void EnhancedTabBar::mouseDoubleClickEvent(QMouseEvent* e)
	{
		if ( e->button() != Qt::LeftButton ) 
		{
			e->ignore();
			return;
    }
		int tab = tabAt_(e->pos());
		if (tab!=-1)
		{
			emit aboutToCloseId(tabData(tab).toInt());
			removeTab(tab);
		}
	}

	int EnhancedTabBar::addTab(const String& text, int id)
	{
		int tab_index = QTabBar::addTab(text.c_str());
    setTabData(tab_index, id);
		
		return tab_index;
	}
	
	void EnhancedTabBar::removeId(int id)
	{
		for (int i=0; i<this->count(); ++i)
    {
    	if (tabData(i).toInt()==id)
    	{
    		removeTab(i);
    		break;
    	}
    }
	}

	void EnhancedTabBar::setCurrentId(int id)
	{
		for (int i=0; i<this->count(); ++i)
    {
    	if (tabData(i).toInt()==id)
    	{
    		setCurrentIndex(i);
    		break;
    	}
    }
	}

	void EnhancedTabBar::currentChanged_(int id)
	{
		emit currentIdChanged(tabData(id).toInt());
	}

	int EnhancedTabBar::tabAt_(const QPoint& pos)
	{
		int tab = -1;

    for (int i=0; i<this->count(); ++i)
    {
			if (tabRect(i).contains(pos))
			{
				tab = i;
				break;
			}
		}

		return tab;
	}

} //namespace OpenMS	

