// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Steffen Sass $
// $Authors: $
// --------------------------------------------------------------------------


#include <OpenMS/DATASTRUCTURES/HashGrid.h>

namespace OpenMS
{

HashGrid::HashGrid(DoubleReal rt_threshold_,DoubleReal mz_threshold_)
{
	rt_threshold=rt_threshold_;
	mz_threshold=mz_threshold_;
	grid_size_x=-1;
	grid_size_y=-1;
	number_of_elements=0;
}

HashGrid::~HashGrid() {
	for (std::map<std::pair<Int,Int>, std::list<GridElement*> >::iterator it=elements.begin();it!=elements.end();++it)
	{
		std::list<GridElement*>& elements=it->second;
		for (std::list<GridElement*>::iterator lit=elements.begin();lit!=elements.end();++lit)
		{
			delete(*lit);
		}
	}

}


void HashGrid::removeElement(GridElement* element_,Int x,Int y)
{
	std::list<GridElement*>& subsets = elements[std::make_pair(x,y)];
	Size previous_size=subsets.size();
	subsets.remove(element_);
	if (subsets.size() < previous_size && number_of_elements>0)
		--number_of_elements;
	if (subsets.empty())
		elements.erase(std::make_pair(x,y));
}

void HashGrid::removeElement(GridElement* element_)
{
	int x = element_->mz / mz_threshold;
	int y = element_->rt / rt_threshold;
	removeElement(element_,x,y);
}

void HashGrid::removeCell(GridCells::iterator loc)
{
	number_of_elements-=loc->second.size();
	elements.erase(loc);
}


void HashGrid::insert(GridElement* element_)
{
	int x = element_->mz / mz_threshold;
	if (x>grid_size_x)
		grid_size_x=x;
	int y = element_->rt / rt_threshold;
	if (y>grid_size_y)
		grid_size_y=y;

	elements[std::make_pair(x,y)].push_back(element_);
	++number_of_elements;
}

void HashGrid::consoleOut()
{
	for (std::map<std::pair<Int,Int>, std::list<GridElement*> >::iterator it=elements.begin();it!=elements.end();++it)
	{
		std::pair<Int,Int> coords=it->first;
		std::list<GridElement*> act_elements= it->second;
		if (it->second.size()>0)
			std::cout << coords.first << "/" << coords.second<< ": ";
		for (std::list<GridElement*>::iterator list_it=act_elements.begin();list_it!=act_elements.end();++list_it)
		{
			std::cout << (*list_it)->getID() << " | ";
		}
		std::cout << std::endl;

	}
	std::cout << std::endl;
}

int HashGrid::size()
{
	return elements.size();
}


DoubleReal HashGrid::getRTThreshold() const
{
	return rt_threshold;
}

DoubleReal HashGrid::getMZThreshold() const
{
	return mz_threshold;
}
Int HashGrid::getGridSizeX()
{
	return grid_size_x;
}

Int HashGrid::getGridSizeY()
{
	return grid_size_y;
}
Int  HashGrid::getNumberOfElements()
{
	return number_of_elements;
}

GridCells::iterator HashGrid::begin()
{
	return elements.begin();
}

GridCells::iterator HashGrid::end()
{
	return elements.end();
}

GridCells::iterator HashGrid::find(std::pair<Int,Int> loc)
{
	return elements.find(loc);
}

}
