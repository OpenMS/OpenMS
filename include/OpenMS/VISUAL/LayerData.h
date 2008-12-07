// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/VISUAL/MultiGradient.h>
#include <OpenMS/VISUAL/Annotations1DContainer.h>

#include <OpenMS/FILTERING/DATAREDUCTION/DataFilters.h>

namespace OpenMS 
{
	/**
		@brief Struct that stores the data for one layer
		
		@ingroup SpectrumWidgets
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
			DT_CONSENSUS,     ///< Consensus feature data
			DT_UNKNOWN			  ///< Undefined data type indicating an error
		};

		///Flags that determine which information is shown.
		enum Flags
		{
			F_HULL,       ///< Features: Overall convex hull
			F_HULLS,      ///< Features: Convex hulls of single mass traces 
			F_NUMBERS,    ///< Features: Number
			P_PRECURSORS, ///< Peaks: Mark precursor peaks of MS/MS scans
			P_PROJECTIONS,///< Peaks: Show projections
			C_ELEMENTS    ///< Consensus features: Show elements
		};
				
		/// Main data type (experiment)
		typedef MSExperiment<> ExperimentType;
		/// Main data type (features)
		typedef FeatureMap<> FeatureMapType;
		/// Main data type (consensus features)
		typedef ConsensusMap ConsensusMapType;
		//@}
		
		/// Default constructor
		LayerData()
			: visible(true),
				flipped(false),
				type(DT_UNKNOWN),
				name(),
				filename(),
				peaks(),
				features(),
				consensus(),
				f1(false),
				f2(false),
				f3(false),
				param(),
				gradient(),
				filters(),
				annotations_1d(),
				modified(false)
		{
		}
		
		/// if this layer is visible
		bool visible;
		/// if this layer is flipped (1d mirror view)
		mutable bool flipped;
		/// data type (peak of feature data)
		DataType type;
		/// layer name
		String name;
		/// file name of the file the data comes from (if available)
		String filename;
		/// peak data
		ExperimentType peaks;
		/// feature data
		FeatureMapType features;
		/// consensus feature data
		ConsensusMapType consensus;
		
		/// Flag one (Features: convex hulls, Peak: precursors, Consensus: elements)
		bool f1;
		/// Flag two (Features: numbers, Peak: projections, Consensus: -)
		bool f2;
		/// Flag tree (Features: convex hull, Peak: -, Consensus: -)
		bool f3;
		
		///Layer parameters
		Param param;
		
		///Gradient for 2D and 3D views
		MultiGradient gradient;
		
		///Filters to apply before painting
		DataFilters filters;
				
		///Annotations for the 1D view
		mutable Annotations1DContainer annotations_1d;
		
		///Flag that indicates that the input data was modified since loading it
		bool modified;
	};

	///Print the contents to a stream.
	std::ostream& operator << (std::ostream& os, const LayerData& rhs);

} //namespace

#endif
