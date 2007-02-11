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

#ifndef OPENMS_VISUAL_LAYERDATA_H
#define OPENMS_VISUAL_LAYERDATA_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/DFeatureMap.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS 
{
	/**
		@brief Struct that stores the data for one layer
		
		@ingroup Visual
	*/
  struct LayerData
	{
		/**	@name Type definitions */
		//@{
		///Dataset types
		enum DataType
		{
			DT_PEAK,		      ///< Peak/Raw data
			DT_FEATURE,	      ///< Feature data
			DT_FEATURE_PAIR,	///< Feature pair data (treated like features, but a line is drawn between even and uneven indices)
			DT_UNKNOWN			  ///< Undefined data type indicating an error
		};
		/// Main data type (experiment)
		typedef MSExperiment<> ExperimentType;
		/// Main data type (features)
		typedef DFeatureMap<2> FeatureMapType;	
		//@}
		
		/// Equality operator
		bool operator==(const LayerData& rhs) const
		{
			//check type
			if (type!=rhs.type || 
					min_int!=rhs.min_int || 
					max_int!=rhs.max_int || 
					visible!=rhs.visible) 
			{
				return false;
			}
			if (type==DT_PEAK)
			{
				if (peaks!=rhs.peaks || reduced!=rhs.reduced) return false;
			}
			else
			{
				if (features!=rhs.features) return false;
			}
			return true;
		}
		
		/// Default constructor
		LayerData()
			: visible(true),
				type(DT_UNKNOWN),
				name(),
				min_int(0.0),
				max_int(std::numeric_limits<double>::max()),
				peaks(),
				reduced(),
				features()
		{	
		}
		
		/// if this layer is visible
		bool visible;
		/// data type (peak of feature data)
		DataType type;
		/// layer name
		String name;
		
		/// minimum displayed intensity
		double min_int;
		/// maximum displayed intensity
		double max_int;
		
		/// peak data
		ExperimentType peaks;
		/// peak data (reduced)
		ExperimentType reduced;
		/// feature data
		FeatureMapType features;

	};

	///Print the contents to a stream.
	std::ostream& operator << (std::ostream& os, const LayerData& rhs);

} //namespace

#endif
