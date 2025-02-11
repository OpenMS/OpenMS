// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/DataValue.h>

#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
#include <OpenMS/DATASTRUCTURES/ParamValue.h>

#include <QtCore/QString>

#include <sstream>

using namespace std;

namespace OpenMS
{

  const DataValue DataValue::EMPTY;

  // default ctor
  DataValue::DataValue() :
    value_type_(EMPTY_VALUE),
    unit_type_(OTHER),
    unit_(-1)
  {
  }

  // destructor
  DataValue::~DataValue()
  {
    clear_();
  }

  const std::string DataValue::NamesOfDataType[] = {
    "String", 
    "Int", 
    "Double", 
    "StringList", 
    "IntList", 
    "DoubleList", 
    "Empty"
    };

  //-------------------------------------------------------------------
  //    ctor for all supported types a DataValue object can hold
  //--------------------------------------------------------------------
  DataValue::DataValue(long double p) :
    value_type_(DOUBLE_VALUE), unit_type_(OTHER), unit_(-1)
  {
    data_.dou_ = p;
  }

  DataValue::DataValue(double p) :
    value_type_(DOUBLE_VALUE), unit_type_(OTHER), unit_(-1)
  {
    data_.dou_ = p;
  }

  DataValue::DataValue(float p) :
    value_type_(DOUBLE_VALUE), unit_type_(OTHER), unit_(-1)
  {
    data_.dou_ = p;
  }

  DataValue::DataValue(short int p) :
    value_type_(INT_VALUE), unit_type_(OTHER), unit_(-1)
  {
    data_.ssize_ = p;
  }

  DataValue::DataValue(unsigned short int p) :
    value_type_(INT_VALUE), unit_type_(OTHER), unit_(-1)
  {
    data_.ssize_ = p;
  }

  DataValue::DataValue(int p) :
    value_type_(INT_VALUE), unit_type_(OTHER), unit_(-1)
  {
    data_.ssize_ = p;
  }

  DataValue::DataValue(unsigned int p) :
    value_type_(INT_VALUE), unit_type_(OTHER), unit_(-1)
  {
    data_.ssize_ = p;
  }

  DataValue::DataValue(long int p) :
    value_type_(INT_VALUE), unit_type_(OTHER), unit_(-1)
  {
    data_.ssize_ = p;
  }

  DataValue::DataValue(unsigned long int p) :
    value_type_(INT_VALUE), unit_type_(OTHER), unit_(-1)
  {
    data_.ssize_ = p;
  }

  DataValue::DataValue(long long p) :
    value_type_(INT_VALUE), unit_type_(OTHER), unit_(-1)
  {
    data_.ssize_ = p;
  }

  DataValue::DataValue(unsigned long long p) :
    value_type_(INT_VALUE), unit_type_(OTHER), unit_(-1)
  {
    data_.ssize_ = p;
  }

  DataValue::DataValue(const char* p) :
    value_type_(STRING_VALUE), unit_type_(OTHER), unit_(-1)
  {
    data_.str_ = new String(p);
  }

  DataValue::DataValue(const string& p) :
    value_type_(STRING_VALUE), unit_type_(OTHER), unit_(-1)
  {
    data_.str_ = new String(p);
  }

  DataValue::DataValue(const QString& p) :
    value_type_(STRING_VALUE), unit_type_(OTHER), unit_(-1)
  {
    data_.str_ = new String(p);
  }

  DataValue::DataValue(const String& p) :
    value_type_(STRING_VALUE), unit_type_(OTHER), unit_(-1)
  {
    data_.str_ = new String(p);
  }

  DataValue::DataValue(const StringList& p) :
    value_type_(STRING_LIST), unit_type_(OTHER), unit_(-1)
  {
    data_.str_list_ = new StringList(p);
  }

  DataValue::DataValue(const IntList& p) :
    value_type_(INT_LIST), unit_type_(OTHER), unit_(-1)
  {
    data_.int_list_ = new IntList(p);
  }

  DataValue::DataValue(const DoubleList& p) :
    value_type_(DOUBLE_LIST), unit_type_(OTHER), unit_(-1)
  {
    data_.dou_list_ = new DoubleList(p);
  }

  DataValue::DataValue(const ParamValue& p) :
    unit_type_(OTHER), unit_(-1)
  {
    switch (p.valueType()) 
    {
    case ParamValue::EMPTY_VALUE:
        value_type_ = EMPTY_VALUE;
    break;
    case ParamValue::INT_VALUE:
        value_type_ = INT_VALUE;
        data_.ssize_ = p;
    break;
    case ParamValue::DOUBLE_VALUE:
        value_type_ = DOUBLE_VALUE;
        data_.dou_ = p;
    break;
    case ParamValue::STRING_VALUE:
        value_type_ = STRING_VALUE;
        data_.str_ = new String(p.toChar());
    break;
    case ParamValue::INT_LIST:
        value_type_ = INT_LIST;
        data_.int_list_ = new IntList(p.toIntVector());
    break;
    case ParamValue::DOUBLE_LIST:
        value_type_ = DOUBLE_LIST;
        data_.dou_list_ = new DoubleList(p.toDoubleVector());
    break;
    case ParamValue::STRING_LIST:
        value_type_ = STRING_LIST;
        data_.str_list_ = new StringList(ListUtils::toStringList<std::string>(p));
    break;
    }
  }

  //--------------------------------------------------------------------
  //                   copy and move constructors
  //--------------------------------------------------------------------
  DataValue::DataValue(const DataValue& p) :
    value_type_(p.value_type_),
    unit_type_(p.unit_type_),
    unit_(p.unit_),
    data_(p.data_)
  {
    if (value_type_ == STRING_VALUE)
    {
      data_.str_ = new String(*(p.data_.str_));
    }
    else if (value_type_ == STRING_LIST)
    {
      data_.str_list_ = new StringList(*(p.data_.str_list_));
    }
    else if (value_type_ == INT_LIST)
    {
      data_.int_list_ = new IntList(*(p.data_.int_list_));
    }
    else if (value_type_ == DOUBLE_LIST)
    {
      data_.dou_list_ = new DoubleList(*(p.data_.dou_list_));
    }
  }

  DataValue::DataValue(DataValue&& rhs) noexcept :
    value_type_(std::move(rhs.value_type_)),
    unit_type_(std::move(rhs.unit_type_)),
    unit_(std::move(rhs.unit_)),
    data_(std::move(rhs.data_))
  {
    // clean up rhs, take ownership of data_
    // NOTE: value_type_ == EMPTY_VALUE implies data_ is empty and can be reset
    rhs.value_type_ = EMPTY_VALUE;
    rhs.unit_type_ = OTHER;
    rhs.unit_ = -1;
  }

  void DataValue::clear_() noexcept
  {
    if (value_type_ == STRING_LIST)
    {
      delete(data_.str_list_);
    }
    else if (value_type_ == STRING_VALUE)
    {
      delete(data_.str_);
    }
    else if (value_type_ == INT_LIST)
    {
      delete(data_.int_list_);
    }
    else if (value_type_ == DOUBLE_LIST)
    {
      delete(data_.dou_list_);
    }

    value_type_ = EMPTY_VALUE;
    unit_type_ = OTHER;
    unit_ = -1;
  }

  //--------------------------------------------------------------------
  //                    copy and move assignment operators
  //--------------------------------------------------------------------
  DataValue& DataValue::operator=(const DataValue& p)
  {
    // Check for self-assignment
    if (this == &p)
    {
      return *this;
    }

    // clean up
    clear_();

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
    value_type_ = p.value_type_;
    unit_type_ = p.unit_type_;
    unit_ = p.unit_;

    return *this;
  }

  /// Move assignment operator
  DataValue& DataValue::operator=(DataValue&& rhs) noexcept
  {
    // Check for self-assignment
    if (this == &rhs)
    {
      return *this;
    }

    // clean up *this
    clear_();

    // assign values to *this
    data_ = rhs.data_;
    value_type_ = rhs.value_type_;
    unit_type_ = rhs.unit_type_;
    unit_ = rhs.unit_;

    // clean up rhs 
    rhs.value_type_ = EMPTY_VALUE;
    rhs.unit_type_ = OTHER;
    rhs.unit_ = -1;

    return *this;
  }

  //--------------------------------------------------------------------
  //                assignment conversion operator
  //--------------------------------------------------------------------

  DataValue& DataValue::operator=(const char* arg)
  {
    clear_();
    data_.str_ = new String(arg);
    value_type_ = STRING_VALUE;
    return *this;
  }

  DataValue& DataValue::operator=(const std::string& arg)
  {
    clear_();
    data_.str_ = new String(arg);
    value_type_ = STRING_VALUE;
    return *this;
  }

  DataValue& DataValue::operator=(const String& arg)
  {
    clear_();
    data_.str_ = new String(arg);
    value_type_ = STRING_VALUE;
    return *this;
  }

  DataValue& DataValue::operator=(const QString& arg)
  {
    clear_();
    data_.str_ = new String(arg);
    value_type_ = STRING_VALUE;
    return *this;
  }

  DataValue& DataValue::operator=(const StringList& arg)
  {
    clear_();
    data_.str_list_ = new StringList(arg);
    value_type_ = STRING_LIST;
    return *this;
  }

  DataValue& DataValue::operator=(const IntList& arg)
  {
    clear_();
    data_.int_list_ = new IntList(arg);
    value_type_ = INT_LIST;
    return *this;
  }

  DataValue& DataValue::operator=(const DoubleList& arg)
  {
    clear_();
    data_.dou_list_ = new DoubleList(arg);
    value_type_ = DOUBLE_LIST;
    return *this;
  }

  DataValue& DataValue::operator=(const long double arg)
  {
    clear_();
    data_.dou_ = arg;
    value_type_ = DOUBLE_VALUE;
    return *this;
  }

  DataValue& DataValue::operator=(const double arg)
  {
    clear_();
    data_.dou_ = arg;
    value_type_ = DOUBLE_VALUE;
    return *this;
  }

  DataValue& DataValue::operator=(const float arg)
  {
    clear_();
    data_.dou_ = arg;
    value_type_ = DOUBLE_VALUE;
    return *this;
  }

  DataValue& DataValue::operator=(const short int arg)
  {
    clear_();
    data_.ssize_ = arg;
    value_type_ = INT_VALUE;
    return *this;
  }

  DataValue& DataValue::operator=(const unsigned short int arg)
  {
    clear_();
    data_.ssize_ = arg;
    value_type_ = INT_VALUE;
    return *this;
  }

  DataValue& DataValue::operator=(const int arg)
  {
    clear_();
    data_.ssize_ = arg;
    value_type_ = INT_VALUE;
    return *this;
  }

  DataValue& DataValue::operator=(const unsigned arg)
  {
    clear_();
    data_.ssize_ = arg;
    value_type_ = INT_VALUE;
    return *this;
  }

  DataValue& DataValue::operator=(const long int arg)
  {
    clear_();
    data_.ssize_ = arg;
    value_type_ = INT_VALUE;
    return *this;
  }

  DataValue& DataValue::operator=(const unsigned long arg)
  {
    clear_();
    data_.ssize_ = arg;
    value_type_ = INT_VALUE;
    return *this;
  }

  DataValue& DataValue::operator=(const long long arg)
  {
    clear_();
    data_.ssize_ = arg;
    value_type_ = INT_VALUE;
    return *this;
  }

  DataValue& DataValue::operator=(const unsigned long long arg)
  {
    clear_();
    data_.ssize_ = arg;
    value_type_ = INT_VALUE;
    return *this;
  }

  //---------------------------------------------------------------------------
  //                      Conversion operators
  //----------------------------------------------------------------------------
  DataValue::operator long double() const
  {
    if (value_type_ == EMPTY_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert DataValue::EMPTY to long double");
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
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert DataValue::EMPTY to double");
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
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert DataValue::EMPTY to float");
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
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
        "Could not convert non-integer DataValue of type '" + NamesOfDataType[value_type_] + "' and value '" + this->toString(true) + "' to short int");
    }
    return data_.ssize_;
  }

  DataValue::operator unsigned short int() const
  {
    if (value_type_ != INT_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
        "Could not convert non-integer DataValue of type '" + NamesOfDataType[value_type_] + "' and value '" + this->toString(true) + "' to unsigned short int");
    }
    if (data_.ssize_ < 0.0)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert negative integer DataValue to unsigned short int");
    }
    return data_.ssize_;
  }

  DataValue::operator int() const
  {
    if (value_type_ != INT_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
        "Could not convert non-integer DataValue of type '" + NamesOfDataType[value_type_] + "' and value '" + this->toString(true) + "' to int");
    }
    return data_.ssize_;
  }

  DataValue::operator unsigned int() const
  {
    if (value_type_ != INT_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
        "Could not convert non-integer DataValue of type '" + NamesOfDataType[value_type_] + "' and value '" + this->toString(true) + "' to unsigned int");
    }
    if (data_.ssize_ < 0.0)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert negative integer DataValue to unsigned int");
    }
    return data_.ssize_;
  }

  DataValue::operator long int() const
  {
    if (value_type_ != INT_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
        "Could not convert non-integer DataValue of type '" + NamesOfDataType[value_type_] + "' and value '" + this->toString(true) + "' to long int");
    }
    return data_.ssize_;
  }

  DataValue::operator unsigned long int() const
  {
    if (value_type_ != INT_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
        "Could not convert non-integer DataValue of type '" + NamesOfDataType[value_type_] + "' and value '" + this->toString(true) + "' to unsigned long int");
    }
    if (data_.ssize_ < 0.0)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert negative integer DataValue to unsigned long int");
    }
    return data_.ssize_;
  }

  DataValue::operator long long() const
  {
    if (value_type_ != INT_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
        "Could not convert non-integer DataValue of type '" + NamesOfDataType[value_type_] + "' and value '" + this->toString(true) + "' to long");
    }
    return data_.ssize_;
  }

  DataValue::operator unsigned long long() const
  {
    if (value_type_ != INT_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
        "Could not convert non-integer DataValue of type '" + NamesOfDataType[value_type_] + "' and value '" + this->toString(true) + "' to unsigned long");
    }
    if (data_.ssize_ < 0.0)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert negative integer DataValue to UInt");
    }
    return data_.ssize_;
  }

  DataValue::operator ParamValue() const
  {
    switch (value_type_)
    {
      case EMPTY_VALUE:
        return ParamValue();
      case INT_VALUE:
        return ParamValue(int(*this));
      case DOUBLE_VALUE:
        return ParamValue(double(*this));
      case STRING_VALUE:
        return ParamValue(std::string(*this));
      case INT_LIST:
        return ParamValue(this->toIntList());
      case DOUBLE_LIST:
        return ParamValue(this->toDoubleList());
      case STRING_LIST:
        {
          // DataValue uses OpenMS::String while ParamValue uses std:string.
          // Therefore the StringList isn't castable.
          vector<std::string> v;
          for (const String& s : this->toStringList())
          {
            v.push_back(s);
          }
          return ParamValue(v);
        }
      default:
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Type of DataValue is unknown!");
    }
    throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Type of DataValue is unknown!");
  }

  DataValue::operator std::string() const
  {
    if (value_type_ != STRING_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
        "Could not convert non-string DataValue of type '" + NamesOfDataType[value_type_] + "' and value '" + this->toString(true) + "' to string");
    }
    return *(data_.str_);
  }

  DataValue::operator StringList() const
  {
    return this->toStringList();
  }

  DataValue::operator IntList() const
  {
    return this->toIntList();
  }

  DataValue::operator DoubleList() const
  {
    return this->toDoubleList();
  }

  // Convert DataValues to char*
  const char* DataValue::toChar() const
  {
    switch (value_type_)
    {
    case DataValue::STRING_VALUE: 
      return const_cast<const char*>(data_.str_->c_str());

    case DataValue::EMPTY_VALUE: 
      return nullptr;

    default: 
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
        "Could not convert DataValue of type '" + NamesOfDataType[value_type_] + "' and value '" + this->toString(true) + "' to char*");
    }
  }

  StringList DataValue::toStringList() const
  {
    if (value_type_ != STRING_LIST)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Could not convert non-StringList DataValue of type '" + NamesOfDataType[value_type_] + "' and value '" + this->toString(true) + "' to StringList");
    }
    return *(data_.str_list_);
  }

  IntList DataValue::toIntList() const
  {
    if (value_type_ != INT_LIST)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Could not convert non-IntList DataValue of type '" + NamesOfDataType[value_type_] + "' and value '" + this->toString(true) + "' to IntList");
    }
    return *(data_.int_list_);
  }

  DoubleList DataValue::toDoubleList() const
  {
    if (value_type_ != DOUBLE_LIST)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
        "Could not convert non-DoubleList DataValue of type '" + NamesOfDataType[value_type_] + "' and value '" + this->toString(true) + "' to DoubleList");
    }
    return *(data_.dou_list_);
  }

  // Convert DataValues to String
  String DataValue::toString(bool full_precision) const
  {
    std::stringstream ss;
    switch (value_type_)
    {
      case DataValue::EMPTY_VALUE: 
        break;
      case DataValue::STRING_VALUE: 
        return *(data_.str_);
      case DataValue::STRING_LIST: ss << *(data_.str_list_); 
        break;
      case DataValue::INT_LIST: ss << *(data_.int_list_); 
        break;
      case DataValue::DOUBLE_LIST: 
        if (full_precision) 
        {
          ss << *(data_.dou_list_);
        }
        else 
        {
          ss << VecLowPrecision<double>(*(data_.dou_list_));
        }
        break;

      case DataValue::INT_VALUE: 
        return String(data_.ssize_);
      case DataValue::DOUBLE_VALUE: 
        return String(data_.dou_, full_precision);

      default:
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
          "Could not convert DataValue of type '" + NamesOfDataType[value_type_] + "' to String");
    }
    return ss.str();
  }

  QString DataValue::toQString() const
  {
    return toString(true).toQString();
  }

  bool DataValue::toBool() const
  {
    if (value_type_ != STRING_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
        "Could not convert non-string DataValue of type '" + NamesOfDataType[value_type_] + "' and value '" + this->toString(true) + "' to bool");
    }
    else if (*(data_.str_) != "true" &&  *(data_.str_) != "false")
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
        String("Could not convert non-string DataValue of type '") + NamesOfDataType[value_type_] + 
        "' and value '" + *(data_.str_) + "' to bool. Valid stings are 'true' and 'false'.");
    }

    return *(data_.str_) == "true";
  }

  // ----------------- Comparator ----------------------

  bool operator==(const DataValue& a, const  DataValue& b)
  {
    if (a.value_type_ == b.value_type_ && a.unit_type_ == b.unit_type_ && a.unit_ == b.unit_)
    {
      switch (a.value_type_)
      {
      case DataValue::EMPTY_VALUE: return b.value_type_ == DataValue::EMPTY_VALUE;

      case DataValue::STRING_VALUE: return *(a.data_.str_) == *(b.data_.str_);

      case DataValue::STRING_LIST: return *(a.data_.str_list_) == *(b.data_.str_list_);

      case DataValue::INT_LIST: return *(a.data_.int_list_) == *(b.data_.int_list_);

      case DataValue::DOUBLE_LIST: return *(a.data_.dou_list_) == *(b.data_.dou_list_);

      case DataValue::INT_VALUE: return a.data_.ssize_ == b.data_.ssize_;

      case DataValue::DOUBLE_VALUE: return fabs(a.data_.dou_ - b.data_.dou_) < 1e-6;

      default: throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Type of DataValue is unkown!");    
      }
    }
    return false;
  }

  bool operator<(const DataValue& a, const  DataValue& b)
  {
    if (a.value_type_ == b.value_type_)
    {
      switch (a.value_type_)
      {
      case DataValue::EMPTY_VALUE: return false;

      case DataValue::STRING_VALUE: return *(a.data_.str_) < *(b.data_.str_);

      case DataValue::STRING_LIST: return a.data_.str_list_->size() < b.data_.str_list_->size();

      case DataValue::INT_LIST: return a.data_.int_list_->size() < b.data_.int_list_->size();

      case DataValue::DOUBLE_LIST: return a.data_.dou_list_->size() < b.data_.dou_list_->size();

      case DataValue::INT_VALUE: return a.data_.ssize_ < b.data_.ssize_;

      case DataValue::DOUBLE_VALUE: return a.data_.dou_ < b.data_.dou_;

      default: throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Type of DataValue is unkown!");    
      }
    }
    return false;
  }

  bool operator>(const DataValue& a, const  DataValue& b)
  {
    if (a.value_type_ == b.value_type_)
    {
      switch (a.value_type_)
      {
      case DataValue::EMPTY_VALUE: return false;

      case DataValue::STRING_VALUE: return *(a.data_.str_) > *(b.data_.str_);

      case DataValue::STRING_LIST: return a.data_.str_list_->size() > b.data_.str_list_->size();

      case DataValue::INT_LIST: return a.data_.int_list_->size() > b.data_.int_list_->size();

      case DataValue::DOUBLE_LIST: return a.data_.dou_list_->size() > b.data_.dou_list_->size();

      case DataValue::INT_VALUE: return a.data_.ssize_ > b.data_.ssize_;

      case DataValue::DOUBLE_VALUE: return a.data_.dou_ > b.data_.dou_;

      default: throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Type of DataValue is unkown!");    

      }
    }
    return false;
  }

  bool operator!=(const DataValue& a, const DataValue& b)
  {
    return !(a == b);
  }

  // ----------------- Output operator ----------------------

  /// for doubles or lists of doubles, you get full precision. Use DataValue::toString(false) if you only need low precision
  std::ostream& operator<<(std::ostream& os, const DataValue& p)
  {
    switch (p.value_type_)
    {
    case DataValue::STRING_VALUE: os << *(p.data_.str_); break;

    case DataValue::STRING_LIST: os << *(p.data_.str_list_); break;

    case DataValue::INT_LIST: os << *(p.data_.int_list_); break;

    case DataValue::DOUBLE_LIST: os << *(p.data_.dou_list_); break;

    case DataValue::INT_VALUE: os << String(p.data_.ssize_); break; // using our String conversion (faster than os)

    case DataValue::DOUBLE_VALUE: os << String(p.data_.dou_); break; // using our String conversion (faster than os)

    case DataValue::EMPTY_VALUE: break;

    default: throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Type of DataValue is unkown!");
    }
    return os;
  }

  // ----------------- Unit methods ----------------------

  const int32_t& DataValue::getUnit() const
  {
    return unit_;
  }

  void DataValue::setUnit(const int32_t& unit)
  {
    unit_ = unit;
  }

} //namespace
