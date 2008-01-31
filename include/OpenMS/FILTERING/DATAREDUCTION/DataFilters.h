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
#include <OpenMS/METADATA/MetaInfoInterface.h>

#include <vector>

namespace OpenMS 
{
	class Feature;
	
	/**
		@brief DataFilter array providing some convenience functions
		
		@todo Write tests (Johannes)
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
				CHARGE,				///< Filter the charge value
				META_DATA			///< Filter meta data
			};
			///Filter operation
			enum FilterOperation
			{
				GREATER_EQUAL,///< Greater than the value or equal to the value
				EQUAL,		    ///< Equal to the value
				LESS_EQUAL,		///< Less than the value or equal to the value
				EXISTS				///< Only for META_DATA filter type, tests if meta data exists
			};

			///Representation of a peak/feature filter combining FilterType, FilterOperation and a value
			struct DataFilter
			{
				///Default constructor
				DataFilter()
					: field(DataFilters::INTENSITY),
						op(DataFilters::GREATER_EQUAL),
						value(0.0),
						value_string(),
						meta_name(),
						value_is_numerical(false)
				{	
				}
				///Field to filter
				FilterType field;
				///Filter operation
				FilterOperation op;
				///Value for comparison
				DoubleReal value;
				///String value for comparison (for meta data)
				String value_string;
				///Name of the considered meta information
				String meta_name;
				///Bool value that is true, if the specified value is numerical, else false 
				bool value_is_numerical;
				
				/// Returns a string representation of the filter
				String toString() const;
				
				/**
					@brief Parses @p filter and sets the filter properties accordingly
					
					This method accepts the format provided by toString().
				*/
				void fromString(const String& filter) throw (Exception::InvalidValue);
				
				///Equality operator
				bool operator==(const DataFilter& rhs) const
				{
					return field==rhs.field && op==rhs.op && value==rhs.value;
				}
				///Inequality operator
				bool operator!=(const DataFilter& rhs) const
				{
					return !operator==(rhs);
				}
				
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
					else if (it->field==META_DATA)
					{
						const MetaInfoInterface& mii = static_cast<MetaInfoInterface>(peak);
						if(!metaPasses(mii,it)) return false;
					}
				}
				return true;
			}

			///Returns if the @p feature fulfills the current filter criteria
			bool passes(const Feature& feature) const;

		protected:
			///Array of DataFilters
			std::vector<DataFilter> filters_;
			
			///Returns if the @p meta_interface (a peak or feature) passes the filter behind iterator @p it 
			inline bool metaPasses(const MetaInfoInterface& meta_interface, std::vector<DataFilter>::const_iterator it) const
			{
				UInt index = meta_interface.metaRegistry().getIndex(it->meta_name);
				if (!meta_interface.metaValueExists(index)) return false;
				else if (it->op!=EXISTS)
				{
					DataValue data_value = meta_interface.getMetaValue(index);
					if(!it->value_is_numerical)
					{
						if(data_value.valueType() != DataValue::STRVALUE) return false;
						else
						{
							// for string values, equality is the only valid operation (besides "exists", see above)
							if(it->op != EQUAL) return false;
							else if(it->value_string != data_value.toString()) return false;
						}	
					}
					else // value_is_numerical
					{
						if (data_value.valueType() == DataValue::STRVALUE || data_value.valueType() == DataValue::EMPTYVALUE) return false;
						else
						{
							if(it->op == EQUAL && (double)data_value != it->value) return false;
							else if(it->op == LESS_EQUAL && (double)data_value > it->value) return false;
							else if(it->op == GREATER_EQUAL && (double)data_value < it->value) return false;
						}
					}
				}
				return true;
			}
	};		

} //namespace

#endif
