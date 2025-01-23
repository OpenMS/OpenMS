// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Authors: Ruben Grünberg $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/ParamValue.h>

#include <OpenMS/DATASTRUCTURES/String.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <sstream>
#include <cmath>

namespace OpenMS
{
  const ParamValue ParamValue::EMPTY;

  // default ctor
  ParamValue::ParamValue() :
    value_type_(EMPTY_VALUE)
  {
  }

  // destructor
  ParamValue::~ParamValue()
  {
    clear_();
  }

  //-------------------------------------------------------------------
  //    ctor for all supported types a ParamValue object can hold
  //--------------------------------------------------------------------
  ParamValue::ParamValue(long double p) :
    value_type_(DOUBLE_VALUE)
  {
    data_.dou_ = p;
  }

  ParamValue::ParamValue(double p) :
    value_type_(DOUBLE_VALUE)
  {
    data_.dou_ = p;
  }

  ParamValue::ParamValue(float p) :
    value_type_(DOUBLE_VALUE)
  {
    data_.dou_ = p;
  }

  ParamValue::ParamValue(short int p) :
    value_type_(INT_VALUE)
  {
    data_.ssize_ = p;
  }

  ParamValue::ParamValue(unsigned short int p) :
    value_type_(INT_VALUE)
  {
    data_.ssize_ = p;
  }

  ParamValue::ParamValue(int p) :
    value_type_(INT_VALUE)
  {
    data_.ssize_ = p;
  }

  ParamValue::ParamValue(unsigned int p) :
    value_type_(INT_VALUE)
  {
    data_.ssize_ = p;
  }

  ParamValue::ParamValue(long int p) :
    value_type_(INT_VALUE)
  {
    data_.ssize_ = p;
  }

  ParamValue::ParamValue(unsigned long int p) :
    value_type_(INT_VALUE)
  {
    data_.ssize_ = p;
  }

  ParamValue::ParamValue(long long p) :
    value_type_(INT_VALUE)
  {
    data_.ssize_ = p;
  }

  ParamValue::ParamValue(unsigned long long p) :
    value_type_(INT_VALUE)
  {
    data_.ssize_ = p;
  }

  ParamValue::ParamValue(const char* p) :
    value_type_(STRING_VALUE)
  {
    data_.str_ = new std::string(p);
  }

  ParamValue::ParamValue(const std::string& p) :
    value_type_(STRING_VALUE)
  {
    data_.str_ = new std::string(p);
  }

  ParamValue::ParamValue(const std::vector<std::string>& p) :
    value_type_(STRING_LIST)
  {
    data_.str_list_ = new std::vector<std::string>(p);
  }

  ParamValue::ParamValue(const std::vector<int>& p) :
    value_type_(INT_LIST)
  {
    data_.int_list_ = new std::vector<int>(p);
  }

  ParamValue::ParamValue(const std::vector<double>& p) :
    value_type_(DOUBLE_LIST)
  {
    data_.dou_list_ = new std::vector<double>(p);
  }

  //--------------------------------------------------------------------
  //                   copy and move constructors
  //--------------------------------------------------------------------
  ParamValue::ParamValue(const ParamValue& p) :
    value_type_(p.value_type_)
  {
    switch (value_type_) 
    {
    case STRING_VALUE:
        data_.str_ = new std::string(*p.data_.str_);
    break;
    case STRING_LIST:
        data_.str_list_ = new std::vector<std::string>(*p.data_.str_list_);
    break;
    case INT_LIST:
        data_.int_list_ = new std::vector<int>(*p.data_.int_list_);
    break;
    case DOUBLE_LIST:
        data_.dou_list_ = new std::vector<double>(*p.data_.dou_list_);
    break;
    default:
        data_ = p.data_;
    break;
    }
  }

  ParamValue::ParamValue(ParamValue&& rhs) noexcept :
    value_type_(std::move(rhs.value_type_)),
    data_(std::move(rhs.data_))
  {
    // clean up rhs, take ownership of data_
    // NOTE: value_type_ == EMPTY_VALUE implies data_ is empty and can be reset
    rhs.value_type_ = EMPTY_VALUE;
  }

  void ParamValue::clear_() noexcept
  {
    switch (value_type_) 
    {
      case STRING_VALUE:
        delete data_.str_;
      break;
      case STRING_LIST:
        delete data_.str_list_;
      break;
      case INT_LIST:
        delete data_.int_list_;
      break;
      case DOUBLE_LIST:
        delete data_.dou_list_;
      break;
      default:
      break;
    }
    value_type_ = EMPTY_VALUE;
  }

  //--------------------------------------------------------------------
  //                    copy and move assignment operators
  //--------------------------------------------------------------------
  ParamValue& ParamValue::operator=(const ParamValue& p)
  {
    // Check for self-assignment
    if (this == &p)
    {
      return *this;
    }

    // clean up
    clear_();

    // assign
    switch (p.value_type_) 
    {
    case STRING_VALUE:
        data_.str_ = new std::string(*p.data_.str_);
        break;
    case STRING_LIST:
        data_.str_list_ = new std::vector<std::string>(*p.data_.str_list_);
        break;
    case INT_LIST:
        data_.int_list_ = new std::vector<int>(*p.data_.int_list_);
        break;
    case DOUBLE_LIST:
        data_.dou_list_ = new std::vector<double>(*p.data_.dou_list_);
        break;
    default:
        data_ = p.data_;
        break;
    }

    // copy type
    value_type_ = p.value_type_;

    return *this;
  }

  /// Move assignment operator
  ParamValue& ParamValue::operator=(ParamValue&& rhs) noexcept
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
    // clean up rhs 
    rhs.value_type_ = EMPTY_VALUE;

    return *this;
  }

  //--------------------------------------------------------------------
  //                assignment conversion operator
  //--------------------------------------------------------------------

  ParamValue& ParamValue::operator=(const char* arg)
  {
    clear_();
    data_.str_ = new std::string(arg);
    value_type_ = STRING_VALUE;
    return *this;
  }

  ParamValue& ParamValue::operator=(const std::string& arg)
  {
    clear_();
    data_.str_ = new std::string(arg);
    value_type_ = STRING_VALUE;
    return *this;
  }

  ParamValue& ParamValue::operator=(const std::vector<std::string>& arg)
  {
    clear_();
    data_.str_list_ = new std::vector<std::string>(arg);
    value_type_ = STRING_LIST;
    return *this;
  }

  ParamValue& ParamValue::operator=(const std::vector<int>& arg)
  {
    clear_();
    data_.int_list_ = new std::vector<int>(arg);
    value_type_ = INT_LIST;
    return *this;
  }

  ParamValue& ParamValue::operator=(const std::vector<double>& arg)
  {
    clear_();
    data_.dou_list_ = new std::vector<double>(arg);
    value_type_ = DOUBLE_LIST;
    return *this;
  }

  ParamValue& ParamValue::operator=(const long double arg)
  {
    clear_();
    data_.dou_ = arg;
    value_type_ = DOUBLE_VALUE;
    return *this;
  }

  ParamValue& ParamValue::operator=(const double arg)
  {
    clear_();
    data_.dou_ = arg;
    value_type_ = DOUBLE_VALUE;
    return *this;
  }

  ParamValue& ParamValue::operator=(const float arg)
  {
    clear_();
    data_.dou_ = arg;
    value_type_ = DOUBLE_VALUE;
    return *this;
  }

  ParamValue& ParamValue::operator=(const short int arg)
  {
    clear_();
    data_.ssize_ = arg;
    value_type_ = INT_VALUE;
    return *this;
  }

  ParamValue& ParamValue::operator=(const unsigned short int arg)
  {
    clear_();
    data_.ssize_ = arg;
    value_type_ = INT_VALUE;
    return *this;
  }

  ParamValue& ParamValue::operator=(const int arg)
  {
    clear_();
    data_.ssize_ = arg;
    value_type_ = INT_VALUE;
    return *this;
  }

  ParamValue& ParamValue::operator=(const unsigned arg)
  {
    clear_();
    data_.ssize_ = arg;
    value_type_ = INT_VALUE;
    return *this;
  }

  ParamValue& ParamValue::operator=(const long int arg)
  {
    clear_();
    data_.ssize_ = arg;
    value_type_ = INT_VALUE;
    return *this;
  }

  ParamValue& ParamValue::operator=(const unsigned long arg)
  {
    clear_();
    data_.ssize_ = arg;
    value_type_ = INT_VALUE;
    return *this;
  }

  ParamValue& ParamValue::operator=(const long long arg)
  {
    clear_();
    data_.ssize_ = arg;
    value_type_ = INT_VALUE;
    return *this;
  }

  ParamValue& ParamValue::operator=(const unsigned long long arg)
  {
    clear_();
    data_.ssize_ = arg;
    value_type_ = INT_VALUE;
    return *this;
  }

  //---------------------------------------------------------------------------
  //                      Conversion operators
  //----------------------------------------------------------------------------
  ParamValue::operator long double() const
  {
    if (value_type_ == EMPTY_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert ParamValue::EMPTY to long double");
    }
    else if (value_type_ == INT_VALUE)
    {
      return (long double)(data_.ssize_);
    }
    return data_.dou_;
  }

  ParamValue::operator double() const
  {
    if (value_type_ == EMPTY_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert ParamValue::EMPTY to double");
    }
    else if (value_type_ == INT_VALUE)
    {
      return double(data_.ssize_);
    }
    return data_.dou_;
  }

  ParamValue::operator float() const
  {
    if (value_type_ == EMPTY_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert ParamValue::EMPTY to float");
    }
    else if (value_type_ == INT_VALUE)
    {
      return float(data_.ssize_);
    }
    return data_.dou_;
  }

  ParamValue::operator short int() const
  {
    if (value_type_ != INT_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert non-integer ParamValue to short int");
    }
    return data_.ssize_;
  }

  ParamValue::operator unsigned short int() const
  {
    if (value_type_ != INT_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert non-integer ParamValue to UInt");
    }
    if (data_.ssize_ < 0.0)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert negative integer ParamValue to unsigned short int");
    }
    return data_.ssize_;
  }

  ParamValue::operator int() const
  {
    if (value_type_ != INT_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert non-integer ParamValue to int");
    }
    return data_.ssize_;
  }

  ParamValue::operator unsigned int() const
  {
    if (value_type_ != INT_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert non-integer ParamValue to unsigned int");
    }
    if (data_.ssize_ < 0.0)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert negative integer ParamValue to unsigned int");
    }
    return data_.ssize_;
  }

  ParamValue::operator long int() const
  {
    if (value_type_ != INT_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert non-integer ParamValue to long int");
    }
    return data_.ssize_;
  }

  ParamValue::operator unsigned long int() const
  {
    if (value_type_ != INT_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert non-integer ParamValue to unsigned long int");
    }
    if (data_.ssize_ < 0.0)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert negative integer ParamValue to unsigned long int");
    }
    return data_.ssize_;
  }

  ParamValue::operator long long() const
  {
    if (value_type_ != INT_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert non-integer ParamValue to Int");
    }
    return data_.ssize_;
  }

  ParamValue::operator unsigned long long() const
  {
    if (value_type_ != INT_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert non-integer ParamValue to UInt");
    }
    if (data_.ssize_ < 0.0)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert negative integer ParamValue to UInt");
    }
    return data_.ssize_;
  }

  ParamValue::operator std::string() const
  {
    if (value_type_ != STRING_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert non-string ParamValue to string");
    }
    return *(data_.str_);
  }

  ParamValue::operator std::vector<std::string>() const
  {
    return this->toStringVector();
  }

  ParamValue::operator std::vector<int>() const
  {
    return this->toIntVector();
  }

  ParamValue::operator std::vector<double>() const
  {
    return this->toDoubleVector();
  }

  // Convert ParamValues to char*
  const char* ParamValue::toChar() const
  {
    switch (value_type_)
    {
    case STRING_VALUE:
        return data_.str_->c_str();
    break;
    case EMPTY_VALUE:
        return nullptr;
    break;
    default:
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert non-string ParamValue to char*");
    break;
    }
  }

  std::string ParamValue::toString(bool full_precision) const 
  {
    std::string str;
    switch (value_type_) 
    {
    case EMPTY_VALUE:
        return "";
    break;
    case STRING_VALUE:
        return *data_.str_;
    break;
    case INT_VALUE:
        return std::to_string(data_.ssize_);
    break;
    case DOUBLE_VALUE:
    {
        return doubleToString(data_.dou_, full_precision);
    }
    break;
    case STRING_LIST:
        str = "[";
        if (!data_.str_list_->empty()) 
        {
            for (std::vector<std::string>::const_iterator it = data_.str_list_->begin();
                 it != data_.str_list_->end() - 1; ++it) 
            {
              str += *it + ", ";
            }
            str += data_.str_list_->back();
        }
        str += "]";
    break;
    case INT_LIST:
        str = "[";
        if (!data_.int_list_->empty()) {
            for (std::vector<int>::const_iterator it = data_.int_list_->begin();
                 it != data_.int_list_->end() - 1; ++it)
            {
                str += std::to_string(*it) + ", ";
            }
            str += std::to_string(data_.int_list_->back());
        }
        str += "]";
    break;
    case DOUBLE_LIST:
        str = "[";
        if (!data_.dou_list_->empty()) {
            for (std::vector<double>::const_iterator it = data_.dou_list_->begin();
                 it != data_.dou_list_->end() - 1; ++it) {
                str += doubleToString(*it, full_precision) + ", ";
            }
            str += doubleToString(data_.dou_list_->back(), full_precision);
        }
        str +=  "]";
    break;
    default:
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert ParamValue to String");
    break;
    }
    return str;
  }

  std::vector<std::string> ParamValue::toStringVector() const 
  {
    if (value_type_ != STRING_LIST)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert non-std::vector<std::string> ParamValue to std::vector<std::string>");
    }
    return *(data_.str_list_);
  }

  std::vector<int> ParamValue::toIntVector() const
  {
    if (value_type_ != INT_LIST)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert non-std::vector<int> ParamValue to std::vector<int>");
    }
    return *(data_.int_list_);
  }

  std::vector<double> ParamValue::toDoubleVector() const {
    if (value_type_ != DOUBLE_LIST)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert non-std::vector<double> ParamValue to std::vector<double>");
    }
    return *(data_.dou_list_);
  }

  bool ParamValue::toBool() const
  {
    if (value_type_ != STRING_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert non-string ParamValue to bool.");
    }
    else if (!(*(data_.str_) == "true" || *(data_.str_) == "false"))
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not convert '" + *(data_.str_) + "' to bool. Valid stings are 'true' and 'false'.");
    }

    return *(data_.str_) == "true";
  }

  // ----------------- Comparator ----------------------

  bool operator==(const ParamValue& a, const  ParamValue& b)
  {
    if (a.value_type_ == b.value_type_)
    {
      switch (a.value_type_)
      {
      case ParamValue::EMPTY_VALUE:
          return true;
      break;
      case ParamValue::STRING_VALUE:
          return *(a.data_.str_) == *(b.data_.str_);
      break;
      case ParamValue::STRING_LIST:
          return *(a.data_.str_list_) == *(b.data_.str_list_);
      break;
      case ParamValue::INT_LIST:
          return *(a.data_.int_list_) == *(b.data_.int_list_);
      break;
      case ParamValue::DOUBLE_LIST:
          return *(a.data_.dou_list_) == *(b.data_.dou_list_);
      break;
      case ParamValue::INT_VALUE:
          return a.data_.ssize_ == b.data_.ssize_;
      break;
      case ParamValue::DOUBLE_VALUE:
          return a.data_.dou_ == b.data_.dou_;
          //return std::fabs(a.data_.dou_ - b.data_.dou_) < 1e-6; This would add an include for <cmath>
      break;
      }
    }
    return false;
  }

  bool operator<(const ParamValue& a, const  ParamValue& b)
  {
    if (a.value_type_ == b.value_type_)
    {
      switch (a.value_type_)
      {
      case ParamValue::EMPTY_VALUE:
          return false;
      break;
      case ParamValue::STRING_VALUE:
          return *(a.data_.str_) < *(b.data_.str_);
      break;
      case ParamValue::STRING_LIST:
          return a.data_.str_list_->size() < b.data_.str_list_->size();
      break;
      case ParamValue::INT_LIST:
          return a.data_.int_list_->size() < b.data_.int_list_->size();
      break;
      case ParamValue::DOUBLE_LIST:
          return a.data_.dou_list_->size() < b.data_.dou_list_->size();
      break;
      case ParamValue::INT_VALUE:
          return a.data_.ssize_ < b.data_.ssize_;
      break;
      case ParamValue::DOUBLE_VALUE:
          return a.data_.dou_ < b.data_.dou_;
      break;
      }
    }
    return false;
  }

  bool operator>(const ParamValue& a, const  ParamValue& b)
  {
    if (a.value_type_ == b.value_type_)
    {
      switch (a.value_type_)
      {
      case ParamValue::EMPTY_VALUE:
          return false;
      break;
      case ParamValue::STRING_VALUE:
          return *(a.data_.str_) > *(b.data_.str_);
      break;
      case ParamValue::STRING_LIST:
          return a.data_.str_list_->size() > b.data_.str_list_->size();
      break;
      case ParamValue::INT_LIST:
          return a.data_.int_list_->size() > b.data_.int_list_->size();
      break;
      case ParamValue::DOUBLE_LIST:
          return a.data_.dou_list_->size() > b.data_.dou_list_->size();
      break;
      case ParamValue::INT_VALUE:
          return a.data_.ssize_ > b.data_.ssize_;
      break;
      case ParamValue::DOUBLE_VALUE:
          return a.data_.dou_ > b.data_.dou_;
      break;
      }
    }
    return false;
  }

  bool operator!=(const ParamValue& a, const ParamValue& b)
  {
    return !(a == b);
  }

  // ----------------- Output operator ----------------------

  /// for doubles or lists of doubles, you get full precision. Use ParamValue::toString(false) if you only need low precision
  std::ostream& operator<<(std::ostream& os, const ParamValue& p)
  {
    switch (p.value_type_)
    {
    case ParamValue::STRING_VALUE:
        os << *(p.data_.str_);
    break;
    case ParamValue::STRING_LIST:
        os << "[";
        if (!p.data_.str_list_->empty())
        {
          for (auto it = p.data_.str_list_->begin(), end = p.data_.str_list_->end() - 1; it != end; ++it) 
          {
            os << *it << ", ";
          }
          os << p.data_.str_list_->back();
        }
        os << "]";
    break;
    case ParamValue::INT_LIST:
        os << "[";
        if (!p.data_.int_list_->empty())
        {
          for (auto it = p.data_.int_list_->begin(), end = p.data_.int_list_->end() - 1; it != end; ++it) 
          {
            os << *it << ", ";
          }
          os << p.data_.int_list_->back();
        }
        os << "]";
    break;
    case ParamValue::DOUBLE_LIST:
        os << "[";
        if (!p.data_.dou_list_->empty()) 
        {
          for (auto it = p.data_.dou_list_->begin(), end = p.data_.dou_list_->end() - 1; it != end; ++it) 
          {
            os << *it << ", ";
          }
          os << p.data_.dou_list_->back();
        }
        os << "]";
    break;
    case ParamValue::INT_VALUE:
        os << p.data_.ssize_;
    break;
    case ParamValue::DOUBLE_VALUE:
        os << p.data_.dou_;
    break;
    case ParamValue::EMPTY_VALUE:
    break;
    }
    return os;
  }

  std::string ParamValue::doubleToString(double value, bool full_precision)
  {
    return String(value, full_precision);
  }

} //namespace
