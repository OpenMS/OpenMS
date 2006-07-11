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
#include <OpenMS/VISUAL/LayerManager.h>
#include <OpenMS/VISUAL/LayerItem.h>

#include <qlineedit.h>

#include <iostream>

using namespace std;
using namespace OpenMS;

LayerManager::LayerManager( QWidget * parent, const char * name, WFlags fl):
	LayerManagerTemplate(parent,name,fl),
	items_(),
	activated_item_(-1)
{
	main_layout_ = new QVBoxLayout(this);
	layout_ = new QVBoxLayout(main_layout_,2);
	main_layout_->addStretch(5);
}

LayerManager::~LayerManager()
{
	
}

void LayerManager::itemActivated(int index)
{
	activate(index);
	emit activatedChanged(index);
}

void LayerManager::activate(int index)
{
	//no change => ignore
	if (index==activated_item_)
	{
		return;
	}
	
	//deactivate previous active item
	if (activated_item_!=-1 && activated_item_<=int(items_.size()))
	{
		items_[activated_item_]->deactivate();
	}
	
	//store current activated item
	activated_item_ = index;
	items_[activated_item_]->activate();
}


void LayerManager::setVisible(UnsignedInt i, bool b)
{
	items_[i]->changeState(b);
}

int LayerManager::addLayer( std::string label )
{
	LayerItem* li = new LayerItem(this);
	li->changeLabel(label);
	li->setIndex(items_.size());
	layout_->addWidget(li);
	items_.push_back(li);
	connect(li,SIGNAL(stateChanged(int, bool)),this,SLOT(itemVisibilityChanged(int, bool)));
	connect(li,SIGNAL(activated(int)),this,SLOT(itemActivated(int)));
	connect(li,SIGNAL(removeRequest(int)),this,SLOT(itemRemoveRequest(int)));
	return (items_.size()-1);
}

void LayerManager::reset()
{
	for(vector<LayerItem*>::iterator it = items_.begin(); it != items_.end() ; ++it)
	{
		delete(*it);
	}
	items_.clear();
	delete(layout_);
	delete(main_layout_);
	main_layout_ = new QVBoxLayout(this);
	layout_ = new QVBoxLayout(main_layout_);
	main_layout_->addStretch(5);
	activated_item_ = -1;
}

void LayerManager::itemVisibilityChanged(int index, bool b)
{
	emit visibilityChanged(index,b);
}

void LayerManager::itemRemoveRequest(int index)
{
	//delete layer item
	LayerItem* item = items_[index];
	layout_->remove(item);
	delete(item);
	items_.erase(items_.begin()+index);
	
	//update activated item
	activated_item_ = -1;
	
	//update layer indices
	for (UnsignedInt i=index; i<items_.size(); ++i)
	{
		items_[i]->setIndex(i);
	}

//	for (UnsignedInt i=0; i<items_.size(); ++i)
//	{
//		cout << items_[i]->getIndex() << " "<< items_[i]->getLabel()<<endl;
//	}
	
	emit removed(index);
}
