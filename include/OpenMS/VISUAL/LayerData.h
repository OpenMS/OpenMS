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
#include <OpenMS/VISUAL/MultiGradient.h>

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
			DT_FEATURE_PAIR,	///< Feature pair data (treated like features, but a line is drawn between even and uneven indices)
			DT_UNKNOWN			  ///< Undefined data type indicating an error
		};

		///Flags that determine which information is shown.
		enum Flags
		{
			F_HULL,       ///< Features: Overall convex hull
			F_HULLS,      ///< Features: Convex hulls of single mass traces 
			F_NUMBERS,    ///< Features: Number
			P_PRECURSORS, ///< Peaks: Mark precursor peaks of MS/MS scans
			P_PROJECTIONS ///< Peaks: Show projections
		};

		///Information to filter
		enum FilterType
		{
			INTENSITY,		///< Filter the intensity value
			QUALITY,		  ///< Filter the overall quality value
			CHARGE		    ///< Filter the charge value
		};
		///Filter operation
		enum FilterOperation
		{
			LESS_EQUAL,		///< Less than the value or equal to the value
			EQUAL,		    ///< Equal to the value
			GREATER_EQUAL ///< Greater than the value or equal to the value
		};		
		///Class that represents a filtering step combining FilterType, FilterOperation and a value
		struct DataFilter
		{
			FilterType type;
			FilterOperation op;
			DoubleReal value;
		};
		///Array of DataFilters
		typedef std::vector<DataFilter> Filters;
		
		/// Main data type (experiment)
		typedef MSExperiment<> ExperimentType;
		/// Main data type (features)
		typedef FeatureMap<> FeatureMapType;	
		//@}
		
		/// Default constructor
		LayerData()
			: visible(true),
				type(DT_UNKNOWN),
				name(),
				peaks(),
				reduced(),
				features(),
				f1(false),
				f2(false),
				f3(false),
				param(),
				gradient(),
				filters()
		{
		}
		
		///Returns if the @p peak fulfills the current filter criteria
		bool passesFilters(const ExperimentType::PeakType& peak) const
		{
			for (Filters::const_iterator it=filters.begin(); it!=filters.end(); ++it)
			{
				if (it->type==INTENSITY)
				{
					if (it->op==GREATER_EQUAL && peak.getIntensity()<it->value) return false;
					else if (it->op==LESS_EQUAL && peak.getIntensity()>it->value) return false;
					else if (it->op==EQUAL && peak.getIntensity()!=it->value) return false;
				}
			}
			return true;
		}
		
		///Returns if the @p feature fulfills the current filter criteria
		bool passesFilters(const FeatureMapType::FeatureType& feature) const
		{
			for (Filters::const_iterator it=filters.begin(); it!=filters.end(); ++it)
			{
				if (it->type==INTENSITY)
				{
					if (it->op==GREATER_EQUAL && feature.getIntensity()<it->value) return false;
					else if (it->op==LESS_EQUAL && feature.getIntensity()>it->value) return false;
					else if (it->op==EQUAL && feature.getIntensity()!=it->value) return false;
				}
				else if (it->type==QUALITY)
				{
					if (it->op==GREATER_EQUAL && feature.getOverallQuality()<it->value) return false;
					else if (it->op==LESS_EQUAL && feature.getOverallQuality()>it->value) return false;
					else if (it->op==EQUAL && feature.getOverallQuality()!=it->value) return false;
				}
				else if (it->type==CHARGE)
				{
					if (it->op==EQUAL && feature.getCharge()!=it->value) return false;
					else if (it->op==GREATER_EQUAL && feature.getCharge()<it->value) return false;
					else if (it->op==LESS_EQUAL && feature.getCharge()>it->value) return false;
				}
			}
			return true;
		}		
		
		/// if this layer is visible
		bool visible;
		/// data type (peak of feature data)
		DataType type;
		/// layer name
		String name;
		
		/// peak data
		ExperimentType peaks;
		/// peak data (reduced is used only in 3D more right now)
		ExperimentType reduced;
		/// feature data
		FeatureMapType features;
		
		/// Flag one (Features: convex hulls, Peak: precursors)
		bool f1;
		/// Flag two (Features: numbers, Peak: projections)
		bool f2;
		/// Flag tree (Features: convex hull, Peak: -)
		bool f3;
		
		///Layer parameters
		Param param;
		
		///Gradient for 2D and 3D views
		MultiGradient gradient;
		
		///Filters to apply when painting
		std::vector<DataFilter> filters;
	};

	///Print the contents to a stream.
	std::ostream& operator << (std::ostream& os, const LayerData& rhs);

} //namespace

#endif
