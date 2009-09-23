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
// $Maintainer: Marc Sturm $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_DATAVALUE_H
#define OPENMS_DATASTRUCTURES_DATAVALUE_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/DATASTRUCTURES/IntList.h>
#include <OpenMS/DATASTRUCTURES/DoubleList.h>
#include <QtCore/QString>

namespace OpenMS
{

	/**
		@brief Class to hold strings, numeric values, lists of strings and lists of numeric values.
		
		- To choose one of these types, just use the apropriate constructor.
		- Automatic conversion is supported and throws Exceptions in case of invalid conversions.
		- An empty object is created with the default constructor.

		@ingroup Datastructures
	*/
	class OPENMS_DLLAPI DataValue
	{

		public:

			/// Empty data value for comparisons
			static const DataValue EMPTY;

			/// Supported types for DataValue
			enum DataType
			{
				STRING_VALUE,     ///< string value
				INT_VALUE,        ///< integer value
				DOUBLE_VALUE,     ///< double value
				STRING_LIST,      ///< string list
				INT_LIST,					///< integer list
				DOUBLE_LIST,			///< double list
				EMPTY_VALUE     	///< empty value
			};

			/// @name Constructors and destructors
			//@{
			/// default constructor
			DataValue();
			/// specific constructor for char* (converted to string)
			DataValue(const char*);
			/// specific constructor for std::string values
			DataValue(const std::string&);
			/// specific constructor for string values
			DataValue(const String&);
			/// specific constructor for QString values
			DataValue(const QString&);
			/// specific constructor for string lists
			DataValue(const StringList&);
			/// specific constructor for integer lists
			DataValue(const IntList&);
			/// specific constructor for double lists
			DataValue(const DoubleList&);
			/// specific constructor for long double values (note: the implementation uses DoubleReal)
			DataValue(long double);
			/// specific constructor for double values (note: the implementation uses DoubleReal)
			DataValue(double);
			/// specific constructor for float values (note: the implementation uses DoubleReal)
			DataValue(float);
			/// specific constructor for short int values (note: the implementation uses SignedSize)
			DataValue(short int);
			/// specific constructor for unsigned short int values (note: the implementation uses SignedSize)
			DataValue(unsigned short int);
      /// specific constructor for int values (note: the implementation uses SignedSize)
      DataValue(int);
      /// specific constructor for unsigned int values (note: the implementation uses SignedSize)
      DataValue(unsigned);
      /// specific constructor for long int values (note: the implementation uses SignedSize)
      DataValue(long int);
      /// specific constructor for unsigned long int values (note: the implementation uses SignedSize)
      DataValue(unsigned long);
      /// specific constructor for long long int values (note: the implementation uses SignedSize)
      DataValue(long long);
      /// specific constructor for unsigned long long int values (note: the implementation uses SignedSize)
      DataValue(unsigned long long);
			/// copy constructor
			DataValue(const DataValue&);
			/// destructor
			virtual ~DataValue();
			//@}

			///@name cast operators
			///These methods are used when the DataType is known.
			///If they are applied to a DataValue with the wrong DataType, an exception is thrown.
			//@{
			/**
				@brief conversion operator to string

				@exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
			*/
			operator std::string() const ;
			/**
				@brief conversion operator to string list

				@exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
			*/
			operator StringList() const ;
			/**
				@brief conversion operator to integer list

				@exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
			*/
			operator IntList() const ;
			/**
				@brief conversion operator to double list
	
				@exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
			*/
			operator DoubleList() const ;
			/**
				@brief conversion operator to long double

				Note: The implementation uses typedef DoubleReal (as opposed to float, double, long double.)

				@exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
			*/
			operator long double() const ;
			/**
				@brief conversion operator to double
			
				Note: The implementation uses typedef DoubleReal (as opposed to float, double, long double.)

				@exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
			*/
			operator double() const ;
			/**
				@brief conversion operator to float

				Note: The implementation uses typedef DoubleReal (as opposed to float, double, long double.)

				@exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
			*/
			operator float() const;
			/**
        @brief conversion operator to short int

        Note: The implementation uses typedef SignedSize.

        @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
      */
      operator short int () const;
      /**
        @brief conversion operator to unsigned short int

        Note: The implementation uses typedef SignedSize.

        @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
      */
      operator unsigned short int () const;
      /**
				@brief conversion operator to int

        Note: The implementation uses typedef SignedSize.

				@exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
			*/
      operator int () const;
			/**
				@brief conversion operator to unsigned int

        Note: The implementation uses typedef SignedSize.

				@exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
			*/
      operator unsigned int () const;
			/**
				@brief conversion operator to long int

        Note: The implementation uses typedef SignedSize.

				@exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
			*/
			operator long int () const;
			/**
				@brief conversion operator to unsigned long int

        Note: The implementation uses typedef SignedSize.

				@exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
			*/
			operator unsigned long int () const;
			/**
				@brief conversion operator to long long

        Note: The implementation uses typedef SignedSize.

				@exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
			*/
			operator long long () const;
			/**
				@brief conversion operator to unsigned long long

        Note: The implementation uses typedef SignedSize.

				@exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
			*/
			operator unsigned long long () const;
			/**
				@brief Convert DataValues to char*

				If the DataValue contains a string, a pointer to it's char* is returned.
				If the DataValue is empty, NULL is returned.
			*/
			const char* toChar() const;
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
			friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream&, const DataValue&);

			/// Equality comparator
			friend OPENMS_DLLAPI bool operator==(const DataValue&, const DataValue&);
			/// Equality comparator
			friend OPENMS_DLLAPI bool operator!=(const DataValue&, const DataValue&);

		protected:
			
			/// Type of the currently stored value
			DataType value_type_;

			/// Space to store the data
			union
			{
 				SignedSize ssize_;
				DoubleReal dou_;
				String* str_;
				StringList* str_list_;
				IntList* int_list_;
				DoubleList* dou_list_;
			} data_;
	};
}

#endif // OPENMS_DATASTRUCTURES_DATAVALUE_H
