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

#ifndef OPENMS_FILTERING_DATAREDUCTION_DATAFILTERS_H
#define OPENMS_FILTERING_DATAREDUCTION_DATAFILTERS_H

#include <OpenMS/DATASTRUCTURES/String.h>

#include <vector>

namespace OpenMS 
{
	class Feature;
	
	/**
		@brief DataFilter array providing some convenience functions
		
		@todo Add filtering of metadata, write tests (Marc, Johannes)
		@todo think about speeding up the whole filtering (Marc)
	*/
	class DataFilters
	{
		public:
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
				GREATER_EQUAL,///< Greater than the value or equal to the value
				EQUAL,		    ///< Equal to the value
				LESS_EQUAL		///< Less than the value or equal to the value				
			};

			///Representation of a peak/feature filter combining FilterType, FilterOperation and a value
			struct DataFilter
			{
				///Default constructor
				DataFilter()
					: field(DataFilters::INTENSITY),
						op(DataFilters::GREATER_EQUAL),
						value(0.0)
				{	
				}
				///Field to filter
				FilterType field;
				///Filter operation
				FilterOperation op;
				///Value for comparison
				DoubleReal value;

				/// Returns a string representation of the filter
				String toString() const;
				
				/**
					@brief Parses @p filter and sets the filter properties accordingly
					
					This method accepts the format provided by toString().
				*/
				void fromString(const String& filter) throw (Exception::InvalidValue);
			};
						
			///Filter count
			UInt size() const;
			
			///Filter accessor
			const DataFilter& operator[](UInt index) const;
			
			///Adds a filter
			void add(const DataFilter& filter);
			
			///Removes the filter corresponding to @p index
			void remove(UInt index) throw (Exception::IndexOverflow);
			
			///Replaces the filter corresponding to @p index
			void replace(UInt index, const DataFilter& filter)  throw (Exception::IndexOverflow);
			
			///Removes all filters
			void clear();
			
			///Returns if the @p peak fulfills the current filter criteria
			template<class PeakType>
			bool passes(const PeakType& peak) const
			{
				for (std::vector<DataFilter>::const_iterator it=filters_.begin(); it!=filters_.end(); ++it)
				{
					if (it->field==INTENSITY)
					{
						if (it->op==GREATER_EQUAL && peak.getIntensity()<it->value) return false;
						else if (it->op==LESS_EQUAL && peak.getIntensity()>it->value) return false;
						else if (it->op==EQUAL && peak.getIntensity()!=it->value) return false;
					}
				}
				return true;
			}

			
			///Returns if the @p feature fulfills the current filter criteria
			bool passes(const Feature& feature) const;

		protected:
			///Array of DataFilters
			std::vector<DataFilter> filters_;
	};		

} //namespace

#endif
