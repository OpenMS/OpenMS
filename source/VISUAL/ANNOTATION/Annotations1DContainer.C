// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/ANNOTATION/Annotations1DContainer.h>

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DDistanceItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DTextItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DPeakItem.h>

#include <QtCore/QPoint>
#include <QtCore/QObject>
#include <QtCore/QRectF>
#include <QtGui/QPainter>

#include <iostream>

namespace OpenMS
{	

	Annotations1DContainer::Annotations1DContainer()
		: std::list<Annotation1DItem*>()
	{
	}

	Annotations1DContainer::Annotations1DContainer(const Annotations1DContainer& rhs)
		: std::list<Annotation1DItem*>()
	{
		//copy annotations
		Annotation1DItem* new_item = 0;
		for (ConstIterator it = rhs.begin(); it != rhs.end(); ++it)
		{
			const Annotation1DDistanceItem* distance_item = dynamic_cast<const Annotation1DDistanceItem*>(*it);
			if (distance_item)
			{
				new_item = new Annotation1DDistanceItem(*distance_item);
				push_back(new_item);
				continue;
			}
			const Annotation1DTextItem* text_item = dynamic_cast<const Annotation1DTextItem*>(*it);
			if (text_item)
			{
				new_item = new Annotation1DTextItem(*text_item);
				push_back(new_item);
				continue;
			}
			const Annotation1DPeakItem* peak_item = dynamic_cast<const Annotation1DPeakItem*>(*it);
			if (peak_item)
			{
				new_item = new Annotation1DPeakItem(*peak_item);
				push_back(new_item);
				continue;
			}
		}
	}

	Annotations1DContainer& Annotations1DContainer::operator= (const Annotations1DContainer& rhs)
	{	
		if (this != &rhs)
		{
			//delete existing annotations
 			for (Iterator it = begin(); it != end(); ++it)
 			{
 				delete *it;
 			}
 			//clear list
 			clear();
 			//copy annotations
			Annotation1DItem* new_item = 0;
			for (ConstIterator it = rhs.begin(); it != rhs.end(); ++it)
			{
				const Annotation1DDistanceItem* distance_item = dynamic_cast<const Annotation1DDistanceItem*>(*it);
				if (distance_item)
				{
					new_item = new Annotation1DDistanceItem(*distance_item);
					push_back(new_item);
					continue;
				}
				const Annotation1DTextItem* text_item = dynamic_cast<const Annotation1DTextItem*>(*it);
				if (text_item)
				{
					new_item = new Annotation1DTextItem(*text_item);
					push_back(new_item);
					continue;
				}
				const Annotation1DPeakItem* peak_item = dynamic_cast<const Annotation1DPeakItem*>(*it);
				if (peak_item)
				{
					new_item = new Annotation1DPeakItem(*peak_item);
					push_back(new_item);
					continue;
				}
			}
		}
		return *this;
	}

	Annotations1DContainer::~Annotations1DContainer()
	{
		for (Iterator it = begin(); it != end(); ++it)
		{
			delete *it;
		}
	}
	
	Annotation1DItem* Annotations1DContainer::getItemAt(const QPoint& pos) const
	{
		for (ConstIterator it = begin(); it != end(); ++it)
		{
			if ((*it)->boundingBox().contains(pos))
			{
				return *it;
			}
		}
		return 0;
	}
	
	void Annotations1DContainer::selectItemAt(const QPoint& pos)
	{
		Annotation1DItem* item = getItemAt(pos);
		if (item != 0)
		{
			item->setSelected(true);
		}
	}
	
	void Annotations1DContainer::deselectItemAt(const QPoint& pos)
	{
		Annotation1DItem* item = getItemAt(pos);
		if (item != 0)
		{
			item->setSelected(false);
		}
	}
	
	void Annotations1DContainer::selectAll()
	{
		for (Iterator it = begin(); it != end(); ++it)
		{
			(*it)->setSelected(true);
		}
	}
	
	void Annotations1DContainer::deselectAll()
	{
		for (Iterator it = begin(); it != end(); ++it)
		{
			(*it)->setSelected(false);
		}
	}
		
	void Annotations1DContainer::removeSelectedItems()
	{
		for (Iterator it = begin(); it != end();)
		{
			if ((*it)->isSelected())
			{
				delete *it;
				it = erase(it);
			}
			else
			{
				++it;
			}
		}
	}
	
}//Namespace
