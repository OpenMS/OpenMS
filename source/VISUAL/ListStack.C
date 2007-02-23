// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/VISUAL/ListStack.h>
#include <OpenMS/VISUAL/PreferencesManager.h>

//QT includes
#include <QtGui/QGridLayout>
#include <QtGui/QHeaderView>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

	ListStack::ListStack( QWidget * parent): 
	QWidget(parent), 
	last_(0)
	{
		//layout
		QGridLayout* layout = new QGridLayout(this);
		layout->setSpacing(4);
		layout->setMargin(6);
		
		//listview (left)
		tree_ =  new QTreeWidget(this);
		tree_->setSizePolicy(QSizePolicy::QSizePolicy::Preferred,QSizePolicy::Minimum);
		tree_->setColumnCount(1);
		tree_->setHeaderLabel("Name");
		tree_->setSortingEnabled(false);
		tree_->setHorizontalScrollBarPolicy ( Qt::ScrollBarAlwaysOff );
		tree_->header()->hide();
		layout->addWidget(tree_,0,0);
			
		
		//widget stack (right)
		stack_ = new QStackedWidget(this);
		layout->addWidget(stack_,0,1);
		layout->setColumnStretch(1,1);
		
		connect(tree_,SIGNAL( itemSelectionChanged() ),this,SLOT( raiseActiveWidget_() ));
	}
	
	ListStack::~ListStack()
	{
		
	}
	
	void ListStack::addWidget(std::string name, QWidget* widget, void* creator, bool highlight, void* parent)
	{
		QTreeWidgetItem* i;
		if (parent==0 || w_to_item_[parent]==0)
		{
			if (last_==0)
			{
				i = new QTreeWidgetItem(tree_);
				i->setText(0,name.c_str());
			}
			else
			{
				i = new QTreeWidgetItem(tree_,last_);
				i->setText(0,name.c_str());
			}
			last_ = i;
		}
		else
		{
			i = new QTreeWidgetItem(w_to_item_[parent]);
			i->setText(0,name.c_str());
		}
		if (parent != creator)
		{
			w_to_item_[creator] = i;	
		}
		
		item_to_index_[i] = stack_->addWidget(widget);
		
		if (highlight)
		{
			tree_->setCurrentItem(i);
		}
	}
	
	void ListStack::expand()
	{
		QTreeWidgetItemIterator it( tree_ );
		while ( *it ) 
		{
		  tree_->expandItem(*it);
		  ++it;
		}	
	}
	
	
	QWidget* ListStack::activeWidget()
	{
		return stack_->currentWidget();
	}
	
	void ListStack::raiseActiveWidget_()
	{
		QTreeWidgetItem* item = tree_->selectedItems().first();
		stack_->setCurrentIndex(item_to_index_[item]);
	}

} //namespace
