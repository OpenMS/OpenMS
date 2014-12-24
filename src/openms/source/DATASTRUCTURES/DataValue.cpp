// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>

#include <OpenMS/CONCEPT/PrecisionWrapper.h>

#include <OpenMS/config.h>

#include <QtCore/QString>

#include <cstddef>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>

using namespace std;

namespace OpenMS
{

  const DataValue DataValue::EMPTY;

  // default ctor
  DataValue::DataValue() :
    value_type_(EMPTY_VALUE), unit_("")
  {
  }

  // destructor
  DataValue::~DataValue()
  {
    clear_();
  }

  //-------------------------------------------------------------------
  //    ctor for all supported types a DataValue object can hold
  //--------------------------------------------------------------------
  DataValue::DataValue(long double p) :
    value_type_(DOUBLE_VALUE), unit_("")
  {
    data_.dou_ = p;
  }

  DataValue::DataValue(double p) :
    value_type_(DOUBLE_VALUE), unit_("")
  {
    data_.dou_ = p;
  }

  DataValue::DataValue(float p) :
    value_type_(DOUBLE_VALUE), unit_("")
  {
    data_.dou_ = p;
  }

  DataValue::DataValue(short int p) :
    value_type_(INT_VALUE), unit_("")
  {
    data_.ssize_ = p;
  }

  DataValue::DataValue(unsigned short int p) :
    value_type_(INT_VALUE), unit_("")
  {
    data_.ssize_ = p;
  }

  DataValue::DataValue(int p) :
    value_type_(INT_VALUE), unit_("")
  {
    data_.ssize_ = p;
  }

  DataValue::DataValue(unsigned int p) :
    value_type_(INT_VALUE), unit_("")
  {
    data_.ssize_ = p;
  }

  DataValue::DataValue(long int p) :
    value_type_(INT_VALUE), unit_("")
  {
    data_.ssize_ = p;
  }

  DataValue::DataValue(unsigned long int p) :
    value_type_(INT_VALUE), unit_("")
  {
    data_.ssize_ = p;
  }

  DataValue::DataValue(long long p) :
    value_type_(INT_VALUE), unit_("")
  {
    data_.ssize_ = p;
  }

  DataValue::DataValue(unsigned long long p) :
    value_type_(INT_VALUE), unit_("")
  {
    data_.ssize_ = p;
  }

  DataValue::DataValue(const char* p) :
    value_type_(STRING_VALUE), unit_("")
  {
    data_.str_ = new String(p);
  }

  DataValue::DataValue(const string& p) :
    value_type_(STRING_VALUE), unit_("")
  {
    data_.str_ = new String(p);
  }

  DataValue::DataValue(const QString& p) :
    value_type_(STRING_VALUE), unit_("")
  {
    data_.str_ = new String(p);
  }

  DataValue::DataValue(const String& p) :
    value_type_(STRING_VALUE), unit_("")
  {
    data_.str_ = new String(p);
  }

  DataValue::DataValue(const StringList& p) :
    value_type_(STRING_LIST), unit_("")
  {
    data_.str_list_ = new StringList(p);
  }

  DataValue::DataValue(const IntList& p) :
    value_type_(INT_LIST), unit_("")
  {
    data_.int_list_ = new IntList(p);
  }

  DataValue::DataValue(const DoubleList& p) :
    value_type_(DOUBLE_LIST), unit_("")
  {
    data_.dou_list_ = new DoubleList(p);
  }

  //--------------------------------------------------------------------
  //                       copy constructor
  //--------------------------------------------------------------------
  DataValue::DataValue(const DataValue& p) :
    value_type_(p.value_type_), data_(p.data_)
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

    if (p.hasUnit())
    {
      unit_ = p.unit_;
    }
  }

  void DataValue::clear_()
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
    unit_ = "";
  }

  //--------------------------------------------------------------------
  //                      assignment operator
  //--------------------------------------------------------------------
  DataValue& DataValue::operator=(const DataValue& p)
  {
    // Check for self-assignment
    if (this == &p)
      return *this;

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
    value_type_     = p.value_type_;

    // copy unit if necessary
    if (p.hasUnit())
    {
      unit_ = p.unit_;
    }

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
    if (value_type_ != STRING_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert non-string DataValue to string");
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
    case DataValue::STRING_VALUE: return const_cast<const char*>(data_.str_->c_str());

    case DataValue::EMPTY_VALUE: return NULL;

    default: throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert DataValue to char*");
    }
  }

  StringList DataValue::toStringList() const
  {
    if (value_type_ != STRING_LIST)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert non-StringList DataValue to StringList");
    }
    return *(data_.str_list_);
  }

  IntList DataValue::toIntList() const
  {
    if (value_type_ != INT_LIST)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert non-IntList DataValue to IntList");
    }
    return *(data_.int_list_);
  }

  DoubleList DataValue::toDoubleList() const
  {
    if (value_type_ != DOUBLE_LIST)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert non-DoubleList DataValue to DoubleList");
    }
    return *(data_.dou_list_);
  }

  // Convert DataValues to String
  String DataValue::toString() const
  {
    stringstream ss;
    switch (value_type_)
    {
    case DataValue::EMPTY_VALUE: break;

    case DataValue::STRING_VALUE: return *(data_.str_);

    case DataValue::STRING_LIST: ss << *(data_.str_list_); break;

    case DataValue::INT_LIST: ss << *(data_.int_list_); break;

    case DataValue::DOUBLE_LIST: ss << *(data_.dou_list_); break;

    case DataValue::INT_VALUE: ss << data_.ssize_; break;

    case DataValue::DOUBLE_VALUE: ss << precisionWrapper(data_.dou_); break;

    default: throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert DataValue to String");
    }
    return ss.str();
  }

  QString DataValue::toQString() const
  {
    QString result;
    switch (value_type_)
    {
    case DataValue::EMPTY_VALUE: break;

    case DataValue::STRING_VALUE: result = QString::fromStdString(*(data_.str_)); break;

    case DataValue::STRING_LIST: result = QString::fromStdString(this->toString()); break;

    case DataValue::INT_LIST: result = QString::fromStdString(this->toString()); break;

    case DataValue::DOUBLE_LIST: result = QString::fromStdString(this->toString()); break;

    case DataValue::INT_VALUE: result.setNum(data_.ssize_); break;

    case DataValue::DOUBLE_VALUE: result.setNum(data_.dou_, 'f'); break;

    default: throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert DataValue to QString");
    }
    return result;
  }

  bool DataValue::toBool() const
  {
    if (value_type_ != STRING_VALUE)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not convert non-string DataValue to bool.");
    }
    else if (*(data_.str_) != "true" &&  *(data_.str_) != "false")
    {
      throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Could not convert '") + *(data_.str_) + "' to bool. Valid stings are 'true' and 'false'.");
    }

    return *(data_.str_) == "true";
  }

  // ----------------- Comparator ----------------------

  bool operator==(const DataValue& a, const  DataValue& b)
  {
    if (a.value_type_ == b.value_type_)
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
      }
    }
    return false;
  }

  bool operator!=(const DataValue& a, const DataValue& b)
  {
    return !(a == b);
  }

  // ----------------- Output operator ----------------------

  std::ostream& operator<<(std::ostream& os, const DataValue& p)
  {
    switch (p.value_type_)
    {
    case DataValue::STRING_VALUE: os << *(p.data_.str_); break;

    case DataValue::STRING_LIST: os << *(p.data_.str_list_); break;

    case DataValue::INT_LIST: os << *(p.data_.int_list_); break;

    case DataValue::DOUBLE_LIST: os << *(p.data_.dou_list_); break;

    case DataValue::INT_VALUE: os << p.data_.ssize_; break;

    case DataValue::DOUBLE_VALUE: os << precisionWrapper(p.data_.dou_); break;

    case DataValue::EMPTY_VALUE: break;
    }
    return os;
  }

  // ----------------- Unit methods ----------------------

  const String& DataValue::getUnit() const
  {
    return unit_;
  }

  void DataValue::setUnit(const OpenMS::String& unit)
  {
    unit_ = unit;
  }

} //namespace
