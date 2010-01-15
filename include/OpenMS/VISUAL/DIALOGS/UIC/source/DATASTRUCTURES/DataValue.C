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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/CONCEPT/Types.h>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>

using namespace std;

namespace OpenMS
{

	const DataValue DataValue::EMPTY;
	
	// default ctor
	DataValue::DataValue() : value_type_(EMPTY_VALUE) 
	{
		
	}
	
	// destructor
	DataValue::~DataValue()
	{
		if (value_type_ == STRING_VALUE)
		{
			delete (data_.str_);
		}
		else if (value_type_ == STRING_LIST)
		{
			delete (data_.str_list_);
		}
		else if(value_type_ == INT_LIST)
		{
			delete (data_.int_list_);
		}
		else if(value_type_ == DOUBLE_LIST)
		{
			delete (data_.dou_list_);
		}
	}
	
	//-------------------------------------------------------------------
	//    ctor for all supported types a DataValue object can hold
	//--------------------------------------------------------------------
	DataValue::DataValue(long double p) : value_type_(DOUBLE_VALUE)
	{ 
		data_.dou_ = p;
	}
	
	DataValue::DataValue(double p) : value_type_(DOUBLE_VALUE)
	{ 
		data_.dou_ = p;
	}
	
	DataValue::DataValue(float p) : value_type_(DOUBLE_VALUE)
	{ 
		data_.dou_ = p;
	}
	
	DataValue::DataValue(short int p) : value_type_(INT_VALUE)
	{
		data_.ssize_ = p;
	}

	DataValue::DataValue(unsigned short int p) : value_type_(INT_VALUE)
	{
		data_.ssize_ = p;
	}

  DataValue::DataValue(int p) : value_type_(INT_VALUE)
  {
    data_.ssize_ = p;
  }

  DataValue::DataValue(unsigned int p) : value_type_(INT_VALUE)
  {
    data_.ssize_ = p;
  }

  DataValue::DataValue(long int p) : value_type_(INT_VALUE)
  {
    data_.ssize_ = p;
  }

  DataValue::DataValue(unsigned long int p) : value_type_(INT_VALUE)
	{
		data_.ssize_ = p;
	}

  DataValue::DataValue(long long p) : value_type_(INT_VALUE)
	{
		data_.ssize_ = p;
	}

  DataValue::DataValue(unsigned long long p) : value_type_(INT_VALUE)
	{
		data_.ssize_ = p;
	}

	DataValue::DataValue(const char* p)	:	value_type_(STRING_VALUE)
	{ 
		data_.str_ = new String(p);
	}

	DataValue::DataValue(const string& p): value_type_(STRING_VALUE)
	{ 
		data_.str_ = new String(p);
	}

	DataValue::DataValue(const QString& p): value_type_(STRING_VALUE)
	{ 
		data_.str_ = new String(p);
	}
	
	DataValue::DataValue(const String& p): value_type_(STRING_VALUE)
	{ 
		data_.str_ = new String(p);
	}

	DataValue::DataValue(const StringList& p): value_type_(STRING_LIST)
	{ 
		data_.str_list_ = new StringList(p);
	}
	
	DataValue::DataValue(const IntList& p): value_type_(INT_LIST)
	{ 
		data_.int_list_ = new IntList(p);
	}
	
	DataValue::DataValue(const DoubleList& p): value_type_(DOUBLE_LIST)
	{ 
		data_.dou_list_ = new DoubleList(p);
	}

	
	//--------------------------------------------------------------------
	//                       copy constructor
	//--------------------------------------------------------------------
	DataValue::DataValue(const DataValue& p): value_type_(p.value_type_), data_(p.data_)
	{
		if (value_type_==STRING_VALUE)
		{
			data_.str_ = new String(*(p.data_.str_));
		}
		else if (value_type_==STRING_LIST)
		{
			data_.str_list_ = new StringList(*(p.data_.str_list_));
		}
		else if (value_type_==INT_LIST)
		{
			data_.int_list_ = new IntList(*(p.data_.int_list_));
		}
		else if (value_type_==DOUBLE_LIST)
		{
			data_.dou_list_ = new DoubleList(*(p.data_.dou_list_));
		}
	}
	
	//--------------------------------------------------------------------
	//                      assignment operator
	//--------------------------------------------------------------------
	DataValue& DataValue::operator = (const DataValue& p)
	{
		// Check for self-assignment
		if (this==&p) return *this;

		// clean up
		if (value_type_==STRING_LIST)
		{
			delete(data_.str_list_);
		}
		else if (value_type_==STRING_VALUE)
		{
			delete(data_.str_);
		}
		else if (value_type_==INT_LIST)
		{
			delete(data_.int_list_);
		}
		else if (value_type_==DOUBLE_LIST)
		{
			delete(data_.dou_list_);
		}
		
		// assign
		if (p.value_type_ == STRING_LIST)
		{
			data_.str_list_ = new StringList(*(p.data_.str_list_));
		}
		else if (p.value_type_ == STRING_VALUE)
		{
			data_.str_ = new String(*(p.data_.str_));
		}
		else if (p.value_type_ == INT_LIST)
		{
			data_.int_list_ = new IntList(*(p.data_.int_list_));
		}
		else if (p.value_type_ == DOUBLE_LIST)
		{
			data_.dou_list_ = new DoubleList(*(p.data_.dou_list_));
		}
		else
		{
			data_ = p.data_;
		}
		
		// copy type
		value_type_		= p.value_type_;
		
		return *this;
	}
	
	//---------------------------------------------------------------------------
	//                      Conversion operators
	//----------------------------------------------------------------------------
	DataValue::operator long double() const
	{
		if (value_type_ == EMPTY_VALUE)
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert DataValue::EMPTY to long double");
		}
		else if (value_type_ == INT_VALUE) 
		{
		  return (long double)(data_.ssize_);
		}
		return data_.dou_; 
	}
	
	DataValue::operator double() const
	{
		if (value_type_ == EMPTY_VALUE)
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert DataValue::EMPTY to double");
		}
		else if (value_type_ == INT_VALUE) 
		{
		  return double(data_.ssize_);
		}
		return data_.dou_; 
	}
	
	DataValue::operator float() const
	{
		if (value_type_ == EMPTY_VALUE)
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert DataValue::EMPTY to float");
		}
		else if (value_type_ == INT_VALUE) 
		{
		  return float(data_.ssize_);
		}
		return data_.dou_; 
	}
	
	DataValue::operator short int() const
	{
		if (value_type_ != INT_VALUE)	
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert non-integer DataValue to short int");
		}
		return data_.ssize_;
	}

	DataValue::operator unsigned short int() const
	{
		if (value_type_ != INT_VALUE)
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert non-integer DataValue to UInt");
		}
		if (data_.ssize_ < 0.0)
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert negative integer DataValue to unsigned short int");
		}
		return data_.ssize_;
	}

  DataValue::operator int() const
  {
    if (value_type_ != INT_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert non-integer DataValue to int");
    }
    return data_.ssize_;
  }

  DataValue::operator unsigned int() const
  {
    if (value_type_ != INT_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert non-integer DataValue to unsigned int");
    }
    if (data_.ssize_ < 0.0)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert negative integer DataValue to unsigned int");
    }
    return data_.ssize_;
  }

  DataValue::operator long int() const
  {
    if (value_type_ != INT_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert non-integer DataValue to long int");
    }
    return data_.ssize_;
  }

  DataValue::operator unsigned long int() const
	{
		if (value_type_ != INT_VALUE)	
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert non-integer DataValue to unsigned long int");
		}
		if (data_.ssize_ < 0.0)	
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert negative integer DataValue to unsigned long int");
		}
		return data_.ssize_;
	}
	
  DataValue::operator long long() const
  {
    if (value_type_ != INT_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert non-integer DataValue to Int");
    }
    return data_.ssize_;
  }

  DataValue::operator unsigned long long() const
	{
		if (value_type_ != INT_VALUE)	
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert non-integer DataValue to UInt");
		}
		if (data_.ssize_ < 0.0)	
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert negative integer DataValue to UInt");
		}
		return data_.ssize_;
	}
	
	DataValue::operator std::string() const
	{
		if(value_type_ != STRING_VALUE)
		{
		  throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Could not convert non-string DataValue to string");
		}
		return *(data_.str_);
	}

	DataValue::operator StringList() const
	{
		if(value_type_ != STRING_LIST)
		{
		  throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Could not convert non-StringList DataValue to StringList");
		}
		return *(data_.str_list_);
	}

	DataValue::operator IntList() const
	{
		if(value_type_ != INT_LIST)
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Could not convert non-IntList DataValue to IntList");
		}
		return *(data_.int_list_);
	}

	DataValue::operator DoubleList() const
	{
		if(value_type_ != DOUBLE_LIST)
		{
		  throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Could not convert non-DoubleList DataValue to DoubleList");
		}
		return *(data_.dou_list_);
	}
	
	// Convert DataValues to char*
	const char* DataValue::toChar() const
	{
		switch(value_type_) 
		{
			case DataValue::STRING_VALUE: return const_cast<const char*>( data_.str_->c_str() );
			case DataValue::EMPTY_VALUE: return NULL;
			default: throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Could not convert DataValue to char*");
		};
	}

	// Convert DataValues to String
	String DataValue::toString() const
	{
		stringstream ss;
		switch(value_type_) 
		{
			case DataValue::EMPTY_VALUE: break;
			case DataValue::STRING_VALUE: return *(data_.str_); break;
			case DataValue::STRING_LIST: ss << *(data_.str_list_) ; break;
			case DataValue::INT_LIST: ss << *(data_.int_list_) ; break;
			case DataValue::DOUBLE_LIST: ss << *(data_.dou_list_) ; break;
			case DataValue::INT_VALUE: ss << data_.ssize_ ; break;
			case DataValue::DOUBLE_VALUE: ss << precisionWrapper(data_.dou_) ; break;
			default: throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Could not convert DataValue to String");
		};
		return ss.str();
	}

	QString DataValue::toQString() const
	{
		QString result;
		switch(value_type_) 
		{
			case DataValue::EMPTY_VALUE: break;
			case DataValue::STRING_VALUE: result = QString::fromStdString(*(data_.str_)); break;
			case DataValue::STRING_LIST: result = QString::fromStdString(this->toString()) ; break;
			case DataValue::INT_LIST: result = QString::fromStdString(this->toString()) ; break;
			case DataValue::DOUBLE_LIST: result = QString::fromStdString(this->toString()) ; break;
			case DataValue::INT_VALUE: result.setNum(data_.ssize_); break;
			case DataValue::DOUBLE_VALUE: result.setNum(data_.dou_,'f'); break;
			default: throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Could not convert DataValue to QString");
		};
		return result;
	}

	bool DataValue::toBool() const
	{
		if (value_type_ != STRING_VALUE)	
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert non-string DataValue to bool.");
		}
		else if (*(data_.str_)!="true" &&  *(data_.str_)!="false")	
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Could not convert '") + *(data_.str_) + "' to bool. Valid stings are 'true' and 'false'.");
		}

		return (*(data_.str_)=="true");
	}


	// ----------------- Comparator ----------------------
	
	bool operator==(const DataValue& a, const  DataValue& b)
	{
		if (a.value_type_ == b.value_type_)
		{
			switch(a.value_type_) 
			{
				case DataValue::EMPTY_VALUE: return b.value_type_== DataValue::EMPTY_VALUE;
	  		case DataValue::STRING_VALUE: return *(a.data_.str_) == *(b.data_.str_);
	  		case DataValue::STRING_LIST: return *(a.data_.str_list_) == *(b.data_.str_list_);
	  		case DataValue::INT_LIST: return *(a.data_.int_list_) == *(b.data_.int_list_);
	  		case DataValue::DOUBLE_LIST: return *(a.data_.dou_list_)==*(b.data_.dou_list_);
				case DataValue::INT_VALUE: return a.data_.ssize_ == b.data_.ssize_;
			  case DataValue::DOUBLE_VALUE: return fabs(a.data_.dou_ - b.data_.dou_)<1e-6;
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
			case DataValue::STRING_VALUE: os << *(p.data_.str_); break;
			case DataValue::STRING_LIST: os << *(p.data_.str_list_); break;
			case DataValue::INT_LIST: os << *(p.data_.int_list_);break;
			case DataValue::DOUBLE_LIST: os << *(p.data_.dou_list_);break;
			case DataValue::INT_VALUE: os << p.data_.ssize_; break;
			case DataValue::DOUBLE_VALUE: os << precisionWrapper(p.data_.dou_); break;
			case DataValue::EMPTY_VALUE: break;
		};
		return os;
	}
	
} //namespace
