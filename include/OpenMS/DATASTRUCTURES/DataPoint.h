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
// $Maintainer: Steffen Sass $
// $Authors: $
// --------------------------------------------------------------------------


#ifndef OPENMS_DATASTRUCTURES_DATAPOINT_H
#define OPENMS_DATASTRUCTURES_DATAPOINT_H

#include <OpenMS/DATASTRUCTURES/GridElement.h>

namespace OpenMS
{

/**
 * @brief A single data point, which can be stored in a HashGrid.
 * @see HashGrid
 * @ingroup Datastructures
 */

class OPENMS_DLLAPI DataPoint : public GridElement {

public:
	static const Int DOUBLE_TRIPLE = 1;
	static const Int DOUBLE = 2;
	static const Int TRIPLE = 3;
	/**
	 * @brief intensity at RT and m/z
	 */
	std::vector<std::vector<DoubleReal> > intensities;
	
	/**
	 * @brief mass shifts [Da] used in the filter
	 */	
	std::vector<DoubleReal> mass_shifts;

	/**
	 * @brief charge of the cluster (i.e. peptide) which the data point is part of
	 */
	Int charge;

	/**
	 * @brief number of isotopes per peptide of the cluster
	 */
	Int isotopes_per_peptide;
	
	/**
	 * @brief quality of the cluster
	 */
	DoubleReal quality;

	/**
	 * @brief ID number of the cluster the data point belongs to
	 */
	Int cluster_id;

	/**
	 * @brief size of the cluster which the data point is part of
	 */
	Int cluster_size;
	
	/**
	 * @brief ID of the data point
	 */
	Int feature_id;
	
	/**
	 * @brief default constructor
	 */
	DataPoint();
	
	/**
	 * @brief copy constructor
	 * @param this DataPoint will be copied
	 */
	DataPoint(const DataPoint &copyin);
	
	/// destructor
	~DataPoint(){};
  DataPoint& operator = (const DataPoint &rhs);
  bool operator == (const DataPoint &rhs) const;
  bool operator != (const DataPoint &rhs) const;
  bool operator < (const DataPoint &rhs) const;
	/**
	 * @brief gets the ID of the data point
	 */
	Int getID() const;
};
}


#endif /* DATAPOINT_H_ */
