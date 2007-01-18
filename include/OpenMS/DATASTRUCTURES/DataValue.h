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

#ifndef OPENMS_DATASTRUCTURES_DATAVALUE_H
#define OPENMS_DATASTRUCTURES_DATAVALUE_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/Serialization.h>

namespace OpenMS
{

	/**
		@brief Class to hold a string or numeric value (integer, double, float, short integer, long integer)
		
		<UL>
		<LI> To choose one of these types, just use the apropriate constructor.
		<LI> Automatic conversion is supported and throws Exceptions in case of invalid ones.
		<LI> An empty objects is created with the default constructor.
		</UL>
		@ingroup Datastructures, Serialization
	*/
	class DataValue
	{

    public:

			/// Empty data value for comparisons
			static const DataValue EMPTY;
	
			/// Supported types for DataValue
			enum DataType {
				STRVALUE,       ///< String value
				INTVALUE,       ///< Integer value
				DOUVALUE,       ///< Double value
				FLOVALUE,       ///< Float value
				SHOVALUE,       ///< Short value
				LONVALUE,       ///< Long value
				EMPTYVALUE      ///< Empty value
				};
			/// default constructor
			DataValue();
			/// specific constructor for char* (converted to string)
			DataValue(char*);
			/// specific constructor for string
			DataValue(std::string);
			/// specific constructor for double
			DataValue(double);
			/// specific constructor for float
			DataValue(float);
			/// specific constructor for int
			DataValue(int);
			/// specific constructor for short
			DataValue(short);
			/// specific constructor for long
			DataValue(long);
	
			/// returns the type of value stored
			inline DataType valueType() const
			{
				return value_type_;
			}
	
			/// copy constructor
			DataValue(const DataValue&);
	
			/// destructor
			virtual ~DataValue();
	
			/// assignment operator
			DataValue& operator = (const DataValue&);
	
			/// conversion operator to string
			operator std::string() const  throw(Exception::ConversionError);
			/// conversion operator to double
			operator double() const  throw(Exception::ConversionError);
			/// conversion operator to float
			operator float() const throw(Exception::ConversionError);
			/// conversion operator to int
			operator int() const throw(Exception::ConversionError) ;
			/// conversion operator to unsigned int
			operator unsigned int() const throw(Exception::ConversionError);
			/// conversion operator to short
			operator short() const throw(Exception::ConversionError);
			/// conversion operator to long
			operator long() const throw(Exception::ConversionError);
	
			/// test if empty
			bool isEmpty() const;
	
			/// output stream operator
			friend std::ostream& operator<<(std::ostream&, const DataValue&);
	
			/// Equality comparator
			friend bool operator==(const DataValue&, const DataValue&);
			/// Equality comparator
			friend bool operator!=(const DataValue&, const DataValue&);
	
			/**
				@brief Convert DataValues to char*
				
				If the DataValue contains a string, a pointer to it's char* is returned.
				<BR>
				If the DataValue is empty, NULL is returned.
				<BR>
				Else ConversionError is thrown.
			*/	
			const char* toChar() const throw(Exception::ConversionError);
	
			/** 
				@brief Convert DataValues to String.
				
				Used to read out all types of data for writing them to file, so no exceptions are thrown.
				If you expect a DataValue of DataType STRVALUE, you better use the cast operator!
			*/
			std::string toString() const;

	  protected:
	  	/// Type of the currently stored value
			DataType value_type_;
	
			/// Space to store the data
			union
			{
 				int    int_;
				double dou_;
				float	 flo_;
				short  sho_;
				long   lon_;
				std::string* str_;			
			} 
			data_;

		///@name Serialization
		//@{
	 private:

		/// Serialization interface
		template<class Archive>
		void save(Archive & ar, const unsigned int /* version */ ) const
		{
			ar & boost::serialization::make_nvp("typ",value_type_);
			switch (value_type_)
			{
			case INTVALUE:
				ar & boost::serialization::make_nvp("int",data_.int_);
				break;
			case FLOVALUE:
				ar & boost::serialization::make_nvp("flo",data_.flo_);
				break;
			case DOUVALUE:
				ar & boost::serialization::make_nvp("dou",data_.dou_);
				break;
			case SHOVALUE:
				ar & boost::serialization::make_nvp("sho",data_.sho_);
				break;
			case LONVALUE:
				ar & boost::serialization::make_nvp("lon",data_.lon_);
				break;
			case STRVALUE:
				{
					const std::string & data_str_(*data_.str_);
					ar & boost::serialization::make_nvp("str",data_str_);
				}
				break;
			case EMPTYVALUE:
			default:
				break;
			}
		}

		/// Serialization interface
		template<class Archive>
		void load(Archive & ar, const unsigned int /* version */ )
		{
			ar & boost::serialization::make_nvp("typ",value_type_);
			switch (value_type_)
			{
			case INTVALUE:
				ar & boost::serialization::make_nvp("int",data_.int_);
				break;
			case FLOVALUE:
				ar & boost::serialization::make_nvp("flo",data_.flo_);
				break;
			case DOUVALUE:
				ar & boost::serialization::make_nvp("dou",data_.dou_);
				break;
			case SHOVALUE:
				ar & boost::serialization::make_nvp("sho",data_.sho_);
				break;
			case LONVALUE:
				ar & boost::serialization::make_nvp("lon",data_.lon_);
				break;
			case STRVALUE:
				{
					// TODO: avoid copying of string content
					std::string data_str_;
					ar & boost::serialization::make_nvp("str",data_str_);
					*this = DataValue(data_str_);
				}
				break;
			case EMPTYVALUE:
			default:
				break;
			}
		}

		BOOST_SERIALIZATION_SPLIT_MEMBER()

		//@}

		/// Serialization
		friend class boost::serialization::access;
		
	};
}

#endif // OPENMS_DATASTRUCTURES_DATAVALUE_H

