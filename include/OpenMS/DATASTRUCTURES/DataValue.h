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
#include <OpenMS/DATASTRUCTURES/String.h>
#include <QtCore/QString>

namespace OpenMS
{

	/**
		@brief Class to hold a string or numeric value (UInt, Int, Real, DoubleReal)
		
		<UL>
			<LI> To choose one of these types, just use the apropriate constructor.
			<LI> Automatic conversion is supported and throws Exceptions in case of invalid ones.
			<LI> An empty objects is created with the default constructor.
		</UL>
		
		@ingroup Datastructures
	*/
	class DataValue
	{

    public:

			/// Empty data value for comparisons
			static const DataValue EMPTY;
	
			/// Supported types for DataValue
			enum DataType {
				STRING_VALUE,  ///< String value
				INT_VALUE,     ///< UInt/Int value
				DOUBLE_VALUE,  ///< DoubleReal/Real value
				EMPTY_VALUE    ///< Empty value
				};
			/// default constructor
			DataValue();
			/// specific constructor for char* (converted to string)
			DataValue(const char*);
			/// specific constructor for String
			DataValue(const String&);
			/// specific constructor for DoubleReal
			DataValue(DoubleReal);
			/// specific constructor for Real
			DataValue(Real);
			/// specific constructor for Int
			DataValue(Int);
			/// specific constructor for UInt
			DataValue(UInt);
	
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
			operator unsigned int() const throw(Exception::ConversionError) ;
	
			/// test if the value is empty
			inline bool isEmpty() const
			{
				return value_type_ == EMPTY_VALUE;
			}
	
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
			String toString() const;

			/**
				@brief Convert DataValues to a QString. 
			
				This method does not throw an exception, if the data is not string data.
				The data is simply converted.
			*/
			QString toQString() const;

	  protected:
	  	/// Type of the currently stored value
			DataType value_type_;
	
			/// Space to store the data
			union
			{
 				PointerSizeInt int_;
				DoubleReal dou_;
				String* str_;			
			} data_;
	};
}

#endif // OPENMS_DATASTRUCTURES_DATAVALUE_H

