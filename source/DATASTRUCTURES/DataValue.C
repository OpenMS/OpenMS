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

#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/CONCEPT/Types.h>
#include <iostream>
#include <sstream>
#include <math.h>
using namespace std;

namespace OpenMS
{

	const DataValue DataValue::EMPTY;
	
	// default ctor
	DataValue::DataValue() : value_type_(EMPTYVALUE) 
	{
		
	}
	
	// destructor
	DataValue::~DataValue()
	{
		if (value_type_ == STRVALUE)
		{
			delete (data_.str_);
		}
	}
	
	// Is DataValue object empty?
	bool DataValue::isEmpty() const 
	{
		return value_type_ == EMPTYVALUE; 
	}
	
	//-------------------------------------------------------------------
	//    ctor for all supported types a DataValue object can hold
	//--------------------------------------------------------------------
	DataValue::DataValue(double p) : value_type_(DOUVALUE)
	{ 
		data_.dou_ = p;
	}
	
	DataValue::DataValue(float p) : value_type_(FLOVALUE)
	{ 
		data_.flo_ = p;
	}
	
	DataValue::DataValue(int p) : value_type_(INTVALUE)
	{
		data_.int_ = p;
	}
	
	DataValue::DataValue(short p) : value_type_(SHOVALUE)
	{ 
		data_.sho_ = p;
	}

	DataValue::DataValue(long p) : value_type_(LONVALUE)
	{ 
		data_.lon_ = p;
	}
	
	DataValue::DataValue(char* p)	:	value_type_(STRVALUE)
	{ 
		data_.str_ = new string(p);
	}
	
	DataValue::DataValue(std::string p): value_type_(STRVALUE)
	{ 
		data_.str_ = new string(p);
	}
	
	//--------------------------------------------------------------------
	//                       copy constructor
	//--------------------------------------------------------------------
	DataValue::DataValue(const DataValue& p): value_type_(p.value_type_), data_(p.data_)
	{
		if (value_type_==STRVALUE)
		{
			data_.str_ = new string(*(p.data_.str_));
		}
	}
	
	//--------------------------------------------------------------------
	//                      assignment operator
	//--------------------------------------------------------------------
	DataValue& DataValue::operator = (const DataValue& p)
	{
		// Check for self-assignment
		if (this==&p) return *this;
		
		// handle string pointers
		if (p.value_type_ != STRVALUE && value_type_!=STRVALUE)
		{
			data_ = p.data_;
		}
		else if (p.value_type_ == STRVALUE && value_type_==STRVALUE)
		{
			*(data_.str_) = *(p.data_.str_);
		}
		else if (p.value_type_ != STRVALUE && value_type_==STRVALUE)
		{
			delete(data_.str_);
			data_ = p.data_;
		}		
		else if (p.value_type_ == STRVALUE && value_type_!=STRVALUE)
		{
			data_.str_ = new string(*(p.data_.str_));
		}
		
		// copy type
		value_type_		= p.value_type_;
		
		return *this;
	}
	
	//---------------------------------------------------------------------------
	//                      Conversion operator
	//
	// No conversion if DataValue object empty 
	//    -> ConversionError-Exception or
  //    -> return numeric_null_ if DataValue of numerical type is in write mode
	//       (in case that empty values can't be written to file format)
	// No conversion if desired type does not fit DataValue type 
	//    -> ConversionError-Exception
	//----------------------------------------------------------------------------
	DataValue::operator double() const throw(Exception::ConversionError)
	{
		if (value_type_ == EMPTYVALUE)
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert DataValue::EMPTY to double");
		}
		else if (value_type_ == STRVALUE) 
		{
		  return atof(data_.str_->c_str());
		}
		else if (value_type_ == FLOVALUE) 
		{
		  return double(data_.flo_);
		}
		else if (value_type_ == INTVALUE) 
		{
		  return double(data_.int_);
		}
		else if (value_type_ != DOUVALUE) 
		{
		  throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert DataValue to double");
		}
		return data_.dou_; 
	}
	
	DataValue::operator float() const throw(Exception::ConversionError)
	{
		if (value_type_ == EMPTYVALUE)
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert DataValue::EMPTY to float");
		}
		else if (value_type_ == STRVALUE) 
		{
		  return atof(data_.str_->c_str());
		}
		else if (value_type_ == DOUVALUE) 
		{
		  return float(data_.dou_);
		}
		else if (value_type_ == INTVALUE) 
		{
		  return float(data_.int_);
		}
		else if (value_type_ != FLOVALUE)
		{
		  throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Could not convert DataValue to float");
		}

		return data_.flo_; 
	}
	
	DataValue::operator int() const throw(Exception::ConversionError)
	{
		if (value_type_ == EMPTYVALUE)	
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert DataValue::EMPTY to int");
		}
		else if (value_type_ == STRVALUE) 
		{
		  return atoi(data_.str_->c_str());
		}
		else if (value_type_ == FLOVALUE) 
		{
		  return int(data_.flo_);
		}
		else if (value_type_ == DOUVALUE) 
		{
		  return int(data_.dou_);
		}
		else if(value_type_ != INTVALUE) 
		{
		  throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert DataValue to int");
		}
		return data_.int_;
	}

	DataValue::operator unsigned int() const throw(Exception::ConversionError)
	{
		if (value_type_ == EMPTYVALUE)	
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert DataValue::EMPTY to unsigend int");
		}
		else if (value_type_ == STRVALUE) 
		{
		  return abs(atoi(data_.str_->c_str()));
		}
		else if (value_type_ == FLOVALUE) 
		{
		  return (unsigned int)fabs(data_.flo_);
		}
		else if (value_type_ == DOUVALUE) 
		{
		  return (unsigned int)fabs(data_.dou_);
		}
		else if(value_type_ != INTVALUE)
		{
		  throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Could not convert DataValue to unsigned int");
		}

		return UnsignedInt( abs(data_.int_));
	}

	DataValue::operator long() const throw(Exception::ConversionError)
	{
		if (value_type_ == EMPTYVALUE)	
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert DataValue::EMPTY to long");
		}
		else if (value_type_ == STRVALUE) 
		{
		  return atoi(data_.str_->c_str());
		}
		else if (value_type_ == FLOVALUE) 
		{
		  return long(data_.flo_);
		}
		else if (value_type_ == DOUVALUE) 
		{
		  return long(data_.dou_);
		}
		else if (value_type_ == INTVALUE) 
		{
		  return long(data_.int_);
		}
		else if(value_type_ != LONVALUE)
		{
		  throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Could not convert DataValue to long");
		}

		return data_.lon_;
	}
	
	DataValue::operator std::string() const throw(Exception::ConversionError)
	{
		if(value_type_ != STRVALUE)
		{
		  throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Could not convert DataValue to string");
		}
		return *(data_.str_);
	}

	DataValue::operator short() const throw(Exception::ConversionError)
	{
		if (value_type_ == EMPTYVALUE)	
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert DataValue::EMPTY to short");
		}
		else if (value_type_ == STRVALUE) 
		{
		  return atoi(data_.str_->c_str());
		}
		else if (value_type_ == FLOVALUE) 
		{
		  return short(data_.flo_);
		}
		else if (value_type_ == DOUVALUE) 
		{
		  return short(data_.dou_);
		}
		else if (value_type_ == INTVALUE) 
		{
		  return short(data_.int_);
		}
		else if (value_type_ != SHOVALUE)
		{
		  throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Could not convert DataValue to short");
		}

		return data_.sho_;
	}
	
	// Convert DataValues to char*
	const char* DataValue::toChar() const throw(Exception::ConversionError)
	{
		switch(value_type_) 
		{
			case DataValue::STRVALUE: return const_cast<const char*>( data_.str_->c_str() );
			case DataValue::EMPTYVALUE: return NULL;
			default: throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Could not convert DataValue to char*");
		};
	}

	// Convert DataValues to String
	string DataValue::toString() const
	{
		stringstream ss;
		ss.precision(7);
		switch(value_type_) 
		{
			case DataValue::EMPTYVALUE: break;
			case DataValue::STRVALUE: return *(data_.str_); break;
			case DataValue::INTVALUE: ss << data_.int_ ; break;
			case DataValue::DOUVALUE: ss << data_.dou_ ; break;
			case DataValue::FLOVALUE: ss << data_.flo_ ; break;
			case DataValue::SHOVALUE: ss << data_.sho_ ; break;
			case DataValue::LONVALUE: ss << data_.lon_ ; break;
		};
		return ss.str();
	}
	
	// ----------------- Comperator ----------------------
	
	bool operator==(const DataValue& a, const  DataValue& b)
	{
		if (a.value_type_ == b.value_type_)
		{
			switch(a.value_type_) 
			{
				case DataValue::EMPTYVALUE: return true;
	  		case DataValue::STRVALUE: return *(a.data_.str_) == *(b.data_.str_);
				case DataValue::INTVALUE: return a.data_.int_ == b.data_.int_;
			  case DataValue::DOUVALUE: return fabs(a.data_.dou_ - b.data_.dou_)<1e-6;
	  		case DataValue::FLOVALUE: return fabs(a.data_.flo_ - b.data_.flo_)<1e-6;
				case DataValue::SHOVALUE: return a.data_.sho_ == b.data_.sho_;
				case DataValue::LONVALUE: return a.data_.lon_ == b.data_.lon_;
			};
		}
		return false;
	}

	bool operator!=(const DataValue& a, const DataValue& b)
	{
		return !(a==b);
	}
	
	// ----------------- Output operator ----------------------
	
	std::ostream& operator<<(std::ostream& os, const DataValue& p)
	{
		switch(p.value_type_) 
		{
			case DataValue::STRVALUE: os << *(p.data_.str_); break;
			case DataValue::INTVALUE: os << p.data_.int_; break;
			case DataValue::DOUVALUE: os << p.data_.dou_; break;
			case DataValue::FLOVALUE: os << p.data_.flo_; break;
			case DataValue::SHOVALUE: os << p.data_.sho_; break;
			case DataValue::LONVALUE: os << p.data_.lon_; break;
			case DataValue::EMPTYVALUE: break;
		};
		return os;
	}
	
} //namespace
