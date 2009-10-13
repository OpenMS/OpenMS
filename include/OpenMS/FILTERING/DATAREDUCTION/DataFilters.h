// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FILTERING_DATAREDUCTION_DATAFILTERS_H
#define OPENMS_FILTERING_DATAREDUCTION_DATAFILTERS_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <iostream>

namespace OpenMS 
{
	class Feature;
	class ConsensusFeature;
	/**
		@brief DataFilter array providing some convenience functions
		
		@note For features the meta data filtering works on the MetaDataInterface of the Feature.
		For peaks it works on the FloatDataArrays definded in MSSpectrum.
	*/
	class OPENMS_DLLAPI DataFilters
	{
		public:
			DataFilters()
			 : filters_(),
			   meta_indices_(),
			   is_active_(false)
			{
			}
			
			///Information to filter
			enum FilterType
			{
				INTENSITY,		///< Filter the intensity value
				QUALITY,		  ///< Filter the overall quality value
				CHARGE,				///< Filter the charge value
				SIZE,         ///< Filter the number of subordinates/elements
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
			struct OPENMS_DLLAPI DataFilter
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
				///Bool value that indicates if the specified value is numerical
				bool value_is_numerical;
				
				/// Returns a string representation of the filter
				String toString() const;
				
				/**
					@brief Parses @p filter and sets the filter properties accordingly
					
					This method accepts the format provided by toString().
					
					@exception Exception::InvalidValue is thrown when the filter is not formatted properly
				*/
				void fromString(const String& filter);
				
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
			Size size() const;
			
			/**
			  @brief Filter accessor

			  @exception Exception::IndexOverflow is thrown for invalid indices
			*/
			const DataFilter& operator[](Size index) const;
			
			///Adds a filter
			void add(const DataFilter& filter);
			
			/**
			  @brief Removes the filter corresponding to @p index

			  @exception Exception::IndexOverflow is thrown for invalid indices
			*/
			void remove(Size index);
			
			/**
			  @brief Replaces the filter corresponding to @p index

			  @exception Exception::IndexOverflow is thrown for invalid indices
			*/
			void replace(Size index, const DataFilter& filter);
			
			///Removes all filters
			void clear();
			
			///Enables/disables the all the filters
			void setActive(bool is_active);
			
			/**
				@brief Returns if the filters are enabled
				
				They are automatically enabled when a filter is added and
				automatically disabled when the last filter is removed
			*/
			inline bool isActive() const
			{
				return is_active_;
			}
			
			///Returns if the @p feature fulfills the current filter criteria
			bool passes(const Feature& feature) const;

			///Returns if the @p consensus_feature fulfills the current filter criteria
			bool passes(const ConsensusFeature& consensus_feature) const;
			
			///Returns if the @p peak fulfills the current filter criteria
			template<class PeakType>
			inline bool passes(const MSSpectrum<PeakType>& spectrum, Size peak_index) const
			{
				if (!is_active_) return true;
				
				for (Size i = 0; i < filters_.size(); i++)
				{
					const DataFilters::DataFilter& filter = filters_[i];
					if (filter.field==INTENSITY)
					{
						switch (filter.op)
						{
							case GREATER_EQUAL:
								if (spectrum[peak_index].getIntensity()<filter.value) return false;
								break;
							case EQUAL:
								if (spectrum[peak_index].getIntensity()!=filter.value) return false;
								break;
							case LESS_EQUAL:
								if(spectrum[peak_index].getIntensity()>filter.value) return false;
								break;
							default:
								break;
						}
					}
					else if (filter.field==META_DATA)
					{
						const typename MSSpectrum<PeakType>::FloatDataArrays& f_arrays = spectrum.getFloatDataArrays();
						//find the right meta data array
						SignedSize f_index = -1;
						for (Size j=0; j<f_arrays.size(); ++j)
						{
							if (f_arrays[j].getName()==filter.meta_name)
							{
								f_index = j;
								break;
							}
						}
						//if it is present, compare it
						if (f_index!=-1)
						{
							if(filter.op == EQUAL && f_arrays[f_index][peak_index] != filter.value) return false;
							else if(filter.op == LESS_EQUAL && f_arrays[f_index][peak_index] > filter.value) return false;
							else if(filter.op == GREATER_EQUAL && f_arrays[f_index][peak_index] < filter.value) return false;
						}
						
						//if float array not found, search in integer arrays
						const typename MSSpectrum<PeakType>::IntegerDataArrays& i_arrays = spectrum.getIntegerDataArrays();
						//find the right meta data array
						SignedSize i_index = -1;
						for (Size j=0; j<i_arrays.size(); ++j)
						{
							if (i_arrays[j].getName()==filter.meta_name)
							{
								i_index = j;
								break;
							}
						}
						//if it is present, compare it
						if (i_index!=-1)
						{
							if(filter.op == EQUAL && i_arrays[i_index][peak_index] != filter.value) return false;
							else if(filter.op == LESS_EQUAL && i_arrays[i_index][peak_index] > filter.value) return false;
							else if(filter.op == GREATER_EQUAL && i_arrays[i_index][peak_index] < filter.value) return false;
						}
							
						//if it is not present, abort
						if (f_index==-1 && i_index==-1) return false;
					}
				}
				return true;
			}

		protected:
			///Array of DataFilters
			std::vector<DataFilter> filters_;
			///Vector of meta indices acting as index cache
			std::vector<Size> meta_indices_;

			///Determines if the filters are activated
			bool is_active_;
			
			///Returns if the meta value at @p index of @p meta_interface (a peak or feature) passes the @p filter
			inline bool metaPasses_(const MetaInfoInterface& meta_interface, const DataFilters::DataFilter& filter, Size index) const
			{
				if (!meta_interface.metaValueExists((UInt)index)) return false;
				else if (filter.op!=EXISTS)
				{
					const DataValue& data_value = meta_interface.getMetaValue((UInt)index);
					if(!filter.value_is_numerical)
					{
						if(data_value.valueType() != DataValue::STRING_VALUE) return false;
						else
						{
							// for string values, equality is the only valid operation (besides "exists", see above)
							if(filter.op != EQUAL) return false;
							else if(filter.value_string != data_value.toString()) return false;
						}	
					}
					else // value_is_numerical
					{
						if (data_value.valueType() == DataValue::STRING_VALUE || data_value.valueType() == DataValue::EMPTY_VALUE) return false;
						else
						{
							if(filter.op == EQUAL && (double)data_value != filter.value) return false;
							else if(filter.op == LESS_EQUAL && (double)data_value > filter.value) return false;
							else if(filter.op == GREATER_EQUAL && (double)data_value < filter.value) return false;
						}
					}
				}
				return true;
			}
	};		

} //namespace

#endif
