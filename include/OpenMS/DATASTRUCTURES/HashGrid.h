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


#include<map>
#include<list>
#include<set>
#include<vector>
#include<iostream>

#include <OpenMS/DATASTRUCTURES/GridElement.h>

#ifndef OPENMS_DATASTRUCTURES_HASHGRID_H
#define OPENMS_DATASTRUCTURES_HASHGRID_H

namespace OpenMS
{

  typedef std::map<std::pair<Int,Int>, std::list<GridElement*> > GridCells;

  /**
		@brief A data structure, which allows the arrangement of data points with an RT and m/z value in a two-dimensional grid.

		The size of each grid cell is determined by two values, namely <i>rt_threshold</i> and <i>mz_threshold</i>.
		<i>rt_threshold</i> defines the height of a grid cell and <i>mz_threshold</i> the width.
		The data points are stored in specific grid cells and are accessible via geometric hashing.
		So the corresponding cell of each data point can be calculated by dividing the rt and m/z values by its corresponding threshold.

		@image html HashGrid.png

		@ingroup Datastructures
   */

class OPENMS_DLLAPI HashGrid {

  private:

    DoubleReal rt_threshold_;
    DoubleReal mz_threshold_;
    Int grid_size_x_;
    Int grid_size_y_;
    Size number_of_elements_;
    GridCells elements_;

  public:

  /**
   * @brief default constructor
   */
   HashGrid();

  /**
   * @brief detailed constructor
   * @param rt_threshold_ defines the height of each grid cell
   * @param mz_threshold_ defines the width of each grid cell
   */
   HashGrid(DoubleReal rt_threshold_, DoubleReal mz_threshold_);

  /**
   * @brief destructor
   */
   ~HashGrid();

  /**
   * @brief removes an element from the hash grid. The cell, in which the element may be contained, is specified:
   * @param element the element to remove
   * @param x x-value of the cell
   * @param y y-value of the cell
   */
   void removeElement(GridElement * const element, Int x, Int y);

  /**
   * @brief removes an element from the hash grid. The cell, in which the element may be contained, is calculated out of the RT and m/z values of the element:
   * @param element_ the element to remove
   */
   void removeElement(GridElement * const element);

  /**
   * @brief removes the cell at location loc from the hash grid :
   * @param loc location of the cell
   */
   void removeCell(GridCells::iterator loc);

  /**
   * @brief inserts a new element in the grid:
   * @param element the element to insert
   */
   void insert(GridElement * const element);

  /**
   * @brief writes the content of the grid to the console:
   */
   void consoleOut() const;

  /**
   * @brief gets the number of element holding cells
   */
   Size size() const;

  /**
   * @brief gets the height of the cells
   */
   DoubleReal getRTThreshold() const;

  /**
   * @brief gets the width of the cells
   */
   DoubleReal getMZThreshold() const;

  /**
   * @brief gets the number of grids in m/z-direction
   */
   Int getGridSizeX() const;

  /**
   * @brief gets the number of grids in RT-direction
   */
   Int getGridSizeY() const;

  /**
	 * @brief gets the number of elements in the grid
	 */
   Size getNumberOfElements() const;

	/**
	 * @brief gets an iterator to the beginning of the GridCells map
	 */
   GridCells::iterator begin();

	/**
	 * @brief gets an iterator just past the end of the GridCells map
	 */
   GridCells::iterator end();

	/**
	 * @brief gets an iterator to the cell at location loc:
	 * @param loc location of the cell
	 */
   GridCells::iterator find(std::pair<Int,Int> loc);

  };

}

#endif /* HASHGRID_H_ */
