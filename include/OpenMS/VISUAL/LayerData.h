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
#include <OpenMS/KERNEL/FeatureMap.h>
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

		///Flags that determine which information is shown.
		enum Flags
		{
			F_HULLS,      ///< Features: Convex hull
			F_NUMBERS,    ///< Feature: Number
			P_SURFACE,    ///< Peaks: Surface calculated by marching squares
			P_CONTOURS,   ///< Peaks: Contour lines calculated by marching squares
			P_PRECURSORS, ///< Peaks: Mark precursor peaks of MS/MS scans
			P_PROJECTIONS ///< Peaks: Show projections
		};

		/// Main data type (experiment)
		typedef MSExperiment<> ExperimentType;
		/// Main data type (features)
		typedef FeatureMap<> FeatureMapType;	
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
				features(),
				f1(false),
				f2(false),
				f3(false),
				f4(false)
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
		
		/// Flag one (Feature: convex hull, Peak: surface)
		bool f1;
		/// Flag two (Feature: numbers , Peak: contours)
		bool f2;
		/// Flag three (Feature: - , Peak: precursors)
		bool f3;
		/// Flag four (Feature: - , Peak: projections)
		bool f4;
		
		/// Parameters of the layer
		Param param;
	};

	///Print the contents to a stream.
	std::ostream& operator << (std::ostream& os, const LayerData& rhs);

} //namespace

#endif
