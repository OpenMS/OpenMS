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


#ifndef OPENMS_DATASTRUCTURES_GRIDELEMENT_H
#define OPENMS_DATASTRUCTURES_GRIDELEMENT_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/BaseFeature.h>

namespace OpenMS
{

  /**
		@brief Base class of all elements, which can be stored in a HashGrid.

		All GridElements must have an m/z and RT value, as well as an unique identifier.
		@see HashGrid
		@ingroup Datastructures
	*/

  class OPENMS_DLLAPI GridElement {
  public :

  /**
   * @brief default constructor
   */
   GridElement();

  /**
   * @brief destructor
   */
   ~GridElement();

  /** @brief copy constructor
   * 	@param source  this GridElement will be copied
   */
   GridElement(const GridElement& copy);

  /**
   * @brief m/z value of the element
   */
   DoubleReal mz;

  /**
   *@brief RT value of the element
   */
   DoubleReal rt;

	/**
	 * @brief detailed constructor
	 * @param rt RT value of the element
	 * @param mz m/z value of the element
	 */
   GridElement(DoubleReal rt_, DoubleReal mz_);

  /**
   * @brief gets the id of the element
	*/  
	virtual Int getID() const = 0;

  };
}

#endif /* GRIDELEMENT_H_ */
