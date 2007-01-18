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
// $Id: DataValue.sip,v 1.1.1.1 2005/04/04 13:03:04 marc_sturm Exp $
// $Author: marc_sturm $
// --------------------------------------------------------------------------
//

class DataValue
{
%TypeHeaderCode
	#include <OpenMS/KERNEL/DataValue.h>
	using namespace OpenMS;
%End
	public:

	enum DataType {STRVALUE, INTVALUE, DOUVALUE, FLOVALUE, SHOVALUE, LONVALUE, EMPTY};
		
	class InvalidConversion	
		: Exception::Base
	{
		public:
		InvalidConversion(const char* /file/, int /line/, const std::string& /s/)	throw();
		~InvalidConversion() throw();
	};

	class EmptyDataValue 
		: Exception::Base
	{
		public:	
		EmptyDataValue(const char* /file/, int /line/)	throw();
		~EmptyDataValue()	throw();
	};

	static void setNumericNull(long);   
	static void unsetNumericNull();
	static void setWriteMode(long);	
	static void unsetWriteMode();
		
	DataValue();
		
	DataValue(const char*);
	DataValue(std::string);
	DataValue(double);
	// DataValue(float);
	// DataValue(int);
	// DataValue(short);
	// DataValue(long);
		
	DataValue(const DataValue&);
	virtual ~DataValue();

	//	DataValue& operator = (const DataValue&);
	
	DataType valueType() const;

	/*	
	// conversion operators
	operator string() const;
	operator double() const;
	operator float() const;
	operator int() const;
	operator unsigned int() const;
	operator short() const;
	operator long() const;
	*/

		
	bool isEmpty() const;
		
	//		friend std::ostream& operator<<(std::ostream&, const DataValue&);
	// ???? __str__/__repr__

	bool operator == (const DataValue&);
	bool operator != (const DataValue&);
		
	char* toChar();
	std::string toString() const;
};
