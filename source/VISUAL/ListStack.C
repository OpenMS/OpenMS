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
// $Id: ListStack.C,v 1.6 2006/03/28 08:03:39 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------


// OpenMS includes
#include <OpenMS/VISUAL/ListStack.h>
#include <OpenMS/VISUAL/PreferencesManager.h>

//QT includes
#include <qlayout.h>
#include <qheader.h>

//STL
#include <iostream>

using namespace std;
using namespace OpenMS;

ListStack::ListStack( QWidget * parent, const char * name ): 
QWidget(parent,name), 
last_(0)
{
		//layout
		QGridLayout* layout = new QGridLayout(this,1,2);
		layout->setSpacing(4);
		layout->setMargin(6);
		
		//listview (left)
		list_ =  new QListView(this);
		list_->setSizePolicy(QSizePolicy::QSizePolicy::Preferred,QSizePolicy::Minimum);
		list_->addColumn("Name");
		list_->setSorting(-1);
		list_->setHScrollBarMode(QScrollView::AlwaysOff);
		list_->header()->hide();
		layout->addWidget(list_,0,0);
			
		
		//widget stack (right)
		stack_ = new EnhancedWidgetStack(this);
		layout->addWidget(stack_,0,1);
		layout->setColStretch(1,1);
		
		connect(list_,SIGNAL( selectionChanged(QListViewItem*) ),stack_,SLOT( raiseWidget(QListViewItem*) ));
}

ListStack::~ListStack()
{
	
}

void ListStack::addWidget(std::string name, QWidget* widget, void* creator, void* parent)
{
	QListViewItem* i;
	if (parent==0 || w_to_item_[parent]==0)
	{
		if (last_==0)
		{
			i = new QListViewItem(list_,name.c_str());
		}
		else
		{
			i = new QListViewItem(list_,last_,name.c_str());
		}
		last_ = i;
	}
	else
	{
		i = new QListViewItem(w_to_item_[parent],name.c_str());
	}
	if (parent != creator)
	{
		w_to_item_[creator] = i;	
	}
	//cout << "c:" << creator <<" i:" << w_to_item_[creator] << " p:" <<parent<< " i[p]:" <<w_to_item_[parent]<<endl;

	stack_->addWidget(widget,i);

	if ( ((PreferencesManager*)creator)->isActive())  //TODO: This is ugly! creator always PreferenceManager* ?
	{
		list_->clearSelection();
		list_->setSelected(i,true);
		((PreferencesManager*)creator)->setActive(false);
	}
}

void ListStack::expand()
{
	QListViewItemIterator it( list_ );
	while ( it.current() ) 
	{
	  QListViewItem *item = it.current();
	  list_->setOpen(item,true);
	  ++it;
	}	
}


QWidget* ListStack::activeWidget()
{
	return stack_->visibleWidget();
}









