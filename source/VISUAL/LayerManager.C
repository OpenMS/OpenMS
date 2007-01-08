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

namespace OpenMS
{

	LayerManager::LayerManager( QWidget * parent, const char * name, WFlags fl):
		LayerManagerTemplate(parent,name,fl),
		activated_item_(-1),
		count_(0)
	{
	}
	
	LayerManager::~LayerManager()
	{
		reset();
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
		
		if (index>=0 && index<int(count_))
		{
			//deactivate previous active item
			if (activated_item_!=-1)
			{
				//find right item
				QLayoutIterator it = layout_->iterator();
				for (SignedInt i=0; i<activated_item_; ++i)
				{
					++it;
				}
				// deactivate
				((LayerItem*)(it.current()->widget()))->deactivate();
			}
			
			//store current activated item
			activated_item_ = index;
			
			//find right item
			QLayoutIterator it = layout_->iterator();
			for (SignedInt i=0; i<activated_item_; ++i)
			{
				++it;
			}
			// activate
			((LayerItem*)(it.current()->widget()))->activate();
		}
	}
	
	
	void LayerManager::setVisible(UnsignedInt index, bool b)
	{
		if (index<count_)
		{
			//find right item
			QLayoutIterator it = layout_->iterator();
			for (UnsignedInt i=0; i<index; ++i)
			{
				++it;
			}
			// set visibilty
			((LayerItem*)(it.current()->widget()))->changeState(b);
		}
	}
	
	int LayerManager::addLayer( std::string label )
	{
		LayerItem* li = new LayerItem(count_,label,this);
		layout_->addWidget(li);
		connect(li,SIGNAL(stateChanged(int, bool)),this,SLOT(itemVisibilityChanged(int, bool)));
		connect(li,SIGNAL(activated(int)),this,SLOT(itemActivated(int)));
		connect(li,SIGNAL(removeRequest(int)),this,SLOT(itemRemoveRequest(int)));
		connect(li,SIGNAL(preferencesRequest(int)),this,SLOT(itemPreferencesRequest(int)));
		activate(count_);
		return (count_++);
	}
	
	void LayerManager::reset()
	{
		QLayoutIterator it = layout_->iterator();
    while ( it.current() != 0 ) 
    {
    	delete((LayerItem*)(it.current()->widget()));
    }
    count_ = 0;
    activated_item_ = -1;
	}
	
	void LayerManager::itemVisibilityChanged(int index, bool b)
	{
		emit visibilityChanged(index,b);
	}
	
	void LayerManager::itemRemoveRequest(int index)
	{
		if (index>=0 && index<int(count_))
		{
			//find right item
			QLayoutIterator it = layout_->iterator();
			int i=0;
			for (; i<index; ++i)
			{
				++it;
			}
			// delete
			delete((LayerItem*)(it.current()->widget()));
			// update indices
			while ( it.current() != 0 ) 
			{
				((LayerItem*)(it.current()->widget()))->setIndex(i++);
				++it;
			}

			//update activated item
			activated_item_ = -1;
			activate(0);
			
			--count_;
			
			emit removed(index);
		}
	}
	
	void LayerManager::itemPreferencesRequest(int index)
	{
		emit showPreferences(index);
	}
} //namespace OpenMS
