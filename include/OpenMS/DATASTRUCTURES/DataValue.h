// -*- Mode: C++; tab-width: 2; -*-
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
			
			///@name Coinstructors and destructors
			//@{
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
			/// copy constructor
			DataValue(const DataValue&);
			/// destructor
			virtual ~DataValue();
			//@}
			
			///@name cast operators
			///These methods are used when the DataType is known. 
			///If they are applied to a DataValue with the wrong DataType, an exception is thrown.
			//@{
			/// conversion operator to string
			operator std::string() const  throw(Exception::ConversionError);
			/// conversion operator to double
			operator DoubleReal() const  throw(Exception::ConversionError);
			/// conversion operator to float
			operator Real() const throw(Exception::ConversionError);
			/// conversion operator to int
			operator Int() const throw(Exception::ConversionError) ;
			/// conversion operator to unsigned int
			operator UInt() const throw(Exception::ConversionError) ;
			/**
				@brief Convert DataValues to char*
				
				If the DataValue contains a string, a pointer to it's char* is returned.
				If the DataValue is empty, NULL is returned.
			*/	
			const char* toChar() const throw(Exception::ConversionError);
			//@}
			
			///@name conversion operators
			///These methods can be used independent of the DataType. If you already know the DataType, you should use a cast operator!
			/// <BR>For conversion of string DataValues to numeric types, first use toString() and then the conversion methods of String.
			//@{
			///Conversion to String.
			String toString() const;
			///Conversuin to QString
			QString toQString() const;
			/**
				@brief Conversion to bool
			  
			  Converts the strings 'true' and 'false' to a bool.

			  @exception Exception::ConversionError is thrown for non-string parameters and string parameters with values other than 'true' and 'false'.
			*/
			bool toBool() const;
			//@}

			/// returns the type of value stored
			inline DataType valueType() const
			{
				return value_type_;
			}
	
			/// assignment operator
			DataValue& operator = (const DataValue&);
			
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

