// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse, Holger Plattfaut, Steffen Sass $
// --------------------------------------------------------------------------


#include <OpenMS/DATASTRUCTURES/HashGrid.h>

namespace OpenMS
{

  HashGrid::HashGrid(DoubleReal rt_threshold, DoubleReal mz_threshold)
    : rt_threshold_(rt_threshold), mz_threshold_(mz_threshold), grid_size_x_(-1), grid_size_y_(-1), number_of_elements_(0), elements_()
  {

  }

  HashGrid::HashGrid()
  {

  }

  HashGrid::~HashGrid()
  {

  }

  void HashGrid::removeElement(GridElement * const element, const Int x, const Int y)
  {
	  std::list<GridElement*>& subsets = elements_[std::make_pair(x,y)];
    Size previous_size = subsets.size();
	  subsets.remove(element);
	  if (subsets.size() < previous_size && number_of_elements_>0) --number_of_elements_;
    if (subsets.empty()) elements_.erase(std::make_pair(x, y));
  }

  void HashGrid::removeElement(GridElement * const element)
  {
	  int x = element->mz / mz_threshold_;
	  int y = element->rt / rt_threshold_;
    removeElement(element, x, y);
  }

  void HashGrid::removeCell(const GridCells::iterator loc)
  {
    number_of_elements_-= loc->second.size();
	  elements_.erase(loc);
  }

  void HashGrid::insert(GridElement * const element)
  {
	  int x = element->mz / mz_threshold_;
    if (x > grid_size_x_) grid_size_x_ = x;
	  int y = element->rt / rt_threshold_;
    if (y > grid_size_y_) grid_size_y_ = y;

    elements_[std::make_pair(x, y)].push_back(element);
	  ++number_of_elements_;
  }

  void HashGrid::consoleOut() const
  {
    for (std::map<std::pair<Int, Int>, std::list<GridElement*> >::const_iterator it = elements_.begin(); it != elements_.end(); ++it)
	  {
      std::pair<Int, Int> coords = it->first;
		  std::list<GridElement*> act_elements = it->second;
		  if (it->second.size()>0) std::cout << coords.first << "/" << coords.second<< ": ";
      for (std::list<GridElement*>::const_iterator list_it = act_elements.begin();list_it != act_elements.end(); ++list_it)
		  {
			  std::cout << (*list_it)->getID() << " | ";
		  }
		  std::cout << "\n";

	  }
	  std::cout << std::endl;
  }

  Size HashGrid::size() const
  {
	  return elements_.size();
  }

  DoubleReal HashGrid::getRTThreshold() const
  {
	  return rt_threshold_;
  }

  DoubleReal HashGrid::getMZThreshold() const
  {
	  return mz_threshold_;
  }

  Int HashGrid::getGridSizeX() const
  {
	  return grid_size_x_;
  }

  Int HashGrid::getGridSizeY() const
  {
	  return grid_size_y_;
  }

  Size HashGrid::getNumberOfElements() const
  {
	  return number_of_elements_;
  }

  GridCells::iterator HashGrid::begin()
  {
	  return elements_.begin();
  }

  GridCells::iterator HashGrid::end()
  {
	  return elements_.end();
  }

  GridCells::iterator HashGrid::find(const std::pair<Int, Int> loc)
  {
	  return elements_.find(loc);
  }

}
