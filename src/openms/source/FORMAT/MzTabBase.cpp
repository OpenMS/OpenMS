// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Timo Sachsenberg, Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzTabBase.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <cassert>

namespace OpenMS
{

  bool MzTabParameterList::isNull() const
  {
    return parameters_.empty();
  }

  void MzTabParameterList::setNull(bool b)
  {
    if (b) { parameters_.clear(); }
  }

  std::string MzTabParameterList::toCellString() const
  {
    if (isNull())
    {
      return "null";
    }
    else
    {
      std::string ret;
      for (std::vector<MzTabParameter>::const_iterator it = parameters_.begin(); it != parameters_.end(); ++it)
      {
        if (it != parameters_.begin())
        {
          ret += "|";
        }
        ret += it->toCellString();
      }
      return ret;
    }
  }

  void MzTabParameterList::fromCellString(const std::string& s)
  {
    std::string lower = s;
    StringUtils::trim(StringUtils::toLower(lower));

    if (lower == "null")
    {
      setNull(true);
    }
    else
    {
      std::vector<std::string> fields;
      StringUtils::split(s, "|", fields);
      for (Size i = 0; i != fields.size(); ++i)
      {
        MzTabParameter p;
        lower = fields[i];
        StringUtils::trim(StringUtils::toLower(lower));
        if (lower == "null")
        {
          throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,StringUtils::toStr("MzTabParameter in MzTabParameterList must not be null '") + s);
        }
        p.fromCellString(fields[i]);
        parameters_.push_back(p);
      }
    }
  }

  std::vector<MzTabParameter> MzTabParameterList::get() const
  {
    return parameters_;
  }

  void MzTabParameterList::set(const std::vector<MzTabParameter>& parameters)
  {
    parameters_ = parameters;
  }

  MzTabStringList::MzTabStringList() :
      sep_('|')
  {
  }

  void MzTabStringList::setSeparator(char sep)
  {
    sep_ = sep;
  }

  bool MzTabStringList::isNull() const
  {
    return entries_.empty();
  }

  void MzTabStringList::setNull(bool b)
  {
    if (b)
    {
      entries_.clear();
    }
  }

  std::string MzTabStringList::toCellString() const
  {
    if (isNull())
    {
      return "null";
    }
    else
    {
      std::string ret;
      for (std::vector<MzTabString>::const_iterator it = entries_.begin(); it != entries_.end(); ++it)
      {
        if (it != entries_.begin())
        {
          ret += sep_;
        }
        ret += it->toCellString();
      }
      return ret;
    }
  }

  void MzTabStringList::fromCellString(const std::string& s)
  {
    std::string lower = s;
    StringUtils::trim(StringUtils::toLower(lower));

    if (lower == "null")
    {
      setNull(true);
    }
    else
    {
      std::vector<std::string> fields;
      StringUtils::split(s, sep_, fields);
      for (Size i = 0; i != fields.size(); ++i)
      {
        MzTabString ts;
        ts.fromCellString(fields[i]);
        entries_.push_back(ts);
      }
    }
  }

  std::vector<MzTabString> MzTabStringList::get() const
  {
    return entries_;
  }

  void MzTabStringList::set(const std::vector<MzTabString>& entries)
  {
    entries_ = entries;
  }

  MzTabSpectraRef::MzTabSpectraRef() :
      ms_run_(0)
  {
  }

  bool MzTabSpectraRef::isNull() const
  {
    return (ms_run_ < 1) || (spec_ref_.empty());
  }

  void MzTabSpectraRef::setNull(bool b)
  {
    if (b)
    {
      ms_run_ = 0;
      spec_ref_.clear();
    }
  }

  void MzTabSpectraRef::setMSFile(Size index)
  {
    assert(index >= 1);
    if (index >= 1)
    {
      ms_run_ = index;
    }
  }

  void MzTabSpectraRef::setSpecRef(const std::string& spec_ref)
  {
    assert(!spec_ref.empty());
    if (!spec_ref.empty())
    {
      spec_ref_ = spec_ref;
    }
    else
    {
      OPENMS_LOG_WARN << "Spectrum reference not set." << std::endl;
    }
  }

  std::string MzTabSpectraRef::getSpecRef() const
  {
    assert(!isNull());
    return spec_ref_;
  }

  Size MzTabSpectraRef::getMSFile() const
  {
    assert(!isNull());
    return ms_run_;
  }

  void MzTabSpectraRef::setSpecRefFile(const std::string& spec_ref)
  {
    assert(!spec_ref.empty());
    if (!spec_ref.empty())
    {
      spec_ref_ = spec_ref;
    }
  }

  std::string MzTabSpectraRef::toCellString() const
  {
    if (isNull())
    {
      return "null";
    }
    else
    {
      return StringUtils::toStr("ms_run[") + StringUtils::toStr(ms_run_) + "]:" + spec_ref_;
    }
  }

  void MzTabSpectraRef::fromCellString(const std::string& s)
  {
    std::string lower = s;
    StringUtils::trim(StringUtils::toLower(lower));
    if (lower == "null")
    {
      setNull(true);
    }
    else
    {
      std::vector<std::string> fields;
      StringUtils::split(s, ":", fields);
      if (fields.size() != 2)
      {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,StringUtils::toStr("Can not convert to MzTabSpectraRef from '") + s + "'");
      }

      spec_ref_ = fields[1];
      { std::string f0 = fields[0]; StringUtils::remove(f0, ']'); ms_run_ = (Size)(StringUtils::toInt32(f0)); }
    }
  }
  MzTabParameter::MzTabParameter()
      : CV_label_(""),
        accession_(""),
        name_(""),
        value_("")
  {
  }

  bool MzTabParameter::isNull() const
  {
    return CV_label_.empty() && accession_.empty() && name_.empty() && value_.empty();
  }

  void MzTabParameter::setNull(bool b)
  {
    if (b)
    {
      CV_label_.clear();
      accession_.clear();
      name_.clear();
      value_.clear();
    }
  }

  void MzTabParameter::setCVLabel(const std::string& CV_label)
  {
    CV_label_ = CV_label;
  }

  void MzTabParameter::setAccession(const std::string& accession)
  {
    accession_ = accession;
  }

  void MzTabParameter::setName(const std::string& name)
  {
    name_ = name;
  }

  void MzTabParameter::setValue(const std::string& value)
  {
    value_ = value;
  }

  std::string MzTabParameter::getCVLabel() const
  {
    assert(!isNull());
    return CV_label_;
  }

  std::string MzTabParameter::getAccession() const
  {
    assert(!isNull());
    return accession_;
  }

  std::string MzTabParameter::getName() const
  {
    assert(!isNull());
    return name_;
  }

  std::string MzTabParameter::getValue() const
  {
    assert(!isNull());
    return value_;
  }

  std::string MzTabParameter::toCellString() const
  {
    if (isNull())
    {
      return "null";
    }
    else
    {
      std::string ret = "[";
      ret += CV_label_ + ", ";
      ret += accession_ + ", ";

      if (StringUtils::hasSubstring(name_, ", "))
      {
        ret +="\"" + name_ + "\""; // quote name if it contains a ","
      }
      else
      {
        ret += name_;
      }

      ret +=StringUtils::toStr(", ");

      if (StringUtils::hasSubstring(value_, ", "))
      {
        ret +="\"" + value_ + "\""; // quote value if it contains a ","
      }
      else
      {
        ret += value_;
      }

      ret += "]";
      return ret;
    }
  }

  void MzTabParameter::fromCellString(const std::string& s)
  {
    std::string lower = s;
    StringUtils::trim(StringUtils::toLower(lower));
    if (lower == "null")
    {
      setNull(true);
    }
    else
    {
      StringList fields;
      std::string field;
      bool in_quotes = false;
      for (std::string::const_iterator sit = s.begin(); sit != s.end(); ++sit)
      {
        if (*sit == '\"') // start or end of quotes
        {
          in_quotes = !in_quotes;
        }
        else if (*sit == ',') // , encountered
        {
          if (in_quotes) // case 1: , in quote
          {
            field += ','; // add , (no split)
          }
          else // split at , if not in quotes
          {
            fields.push_back(StringUtils::trim(field));
            field.clear();
          }
        }
        else if (*sit != '[' && *sit != ']')
        {
          // skip leading ws
          if (*sit == ' ' && field.empty())
          {
            continue;
          }
          field += *sit;
        }
      }

      fields.push_back(StringUtils::trim(field));

      if (fields.size() != 4)
      {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,StringUtils::toStr("Could not convert String '") + s + "' to MzTabParameter");
      }

      CV_label_ = fields[0];
      accession_ = fields[1];
      name_ = fields[2];
      value_ = fields[3];
    }
  }

  MzTabString::MzTabString(const std::string& s)
  {
    set(s);
  }

  void MzTabString::set(const std::string& value)
  {
    std::string lower = value;
    StringUtils::trim(StringUtils::toLower(lower));
    if (lower == "null")
    {
      setNull(true);
    }
    else
    {
      value_ = value;
      StringUtils::trim(value_);
    }
  }

  std::string MzTabString::get() const
  {
    return value_;
  }

  bool MzTabString::isNull() const
  {
    return value_.empty();
  }

  void MzTabString::setNull(bool b)
  {
    if (b)
    {
      value_.clear();
    }
  }

  MzTabString::MzTabString()
      : value_()
  {
  }

  std::string MzTabString::toCellString() const
  {
    if (isNull())
    {
      return "null";
    }
    else
    {
      return value_;
    }
  }

  void MzTabString::fromCellString(const std::string& s)
  {
    set(s);
  }

  MzTabBoolean::MzTabBoolean(bool v)
  {
    set((int)v);
  }

  MzTabBoolean::MzTabBoolean()
      : value_(-1)
  {
  }

  void MzTabBoolean::set(const bool& value)
  {
    value_ = (int)value;
  }

  Int MzTabBoolean::get() const
  {
    return value_;
  }

  bool MzTabBoolean::isNull() const
  {
    return value_ < 0;
  }

  void MzTabBoolean::setNull(bool b)
  {
    if (!b)
      value_ = -1;
    else
      value_ = 0;
  }

  std::string MzTabBoolean::toCellString() const
  {
    if (isNull())
    {
      return "null";
    }
    else
    {
      if (value_)
      {
        return "1";
      }
      else
      {
        return "0";
      }
    }
  }

  void MzTabBoolean::fromCellString(const std::string& s)
  {
    std::string lower = s;
    StringUtils::trim(StringUtils::toLower(lower));
    if (lower == "null")
    {
      setNull(true);
    }
    else
    {
      if (s == "0")
      {
        set(false);
      }
      else if (s == "1")
      {
        set(true);
      }
      else
      {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,StringUtils::toStr("Could not convert String '") + s + "' to MzTabBoolean");
      }
    }
  }

  bool MzTabIntegerList::isNull() const
  {
    return entries_.empty();
  }

  void MzTabIntegerList::setNull(bool b)
  {
    if (b)
    {
      entries_.clear();
    }
  }

  std::string MzTabIntegerList::toCellString() const
  {
    if (isNull())
    {
      return "null";
    }
    else
    {
      std::string ret;
      for (std::vector<MzTabInteger>::const_iterator it = entries_.begin(); it != entries_.end(); ++it)
      {
        if (it != entries_.begin())
        {
          ret += ",";
        }
        ret += it->toCellString();
      }
      return ret;
    }
  }

  void MzTabIntegerList::fromCellString(const std::string& s)
  {
    std::string lower = s;
    StringUtils::trim(StringUtils::toLower(lower));
    if (lower == "null")
    {
      setNull(true);
    }
    else
    {
      std::vector<std::string> fields;
      StringUtils::split(s, ",", fields);
      for (Size i = 0; i != fields.size(); ++i)
      {
        MzTabInteger ds;
        ds.fromCellString(fields[i]);
        entries_.push_back(ds);
      }
    }
  }

  std::vector<MzTabInteger> MzTabIntegerList::get() const
  {
    return entries_;
  }

  void MzTabIntegerList::set(const std::vector<MzTabInteger>& entries)
  {
    entries_ = entries;
  }

  MzTabInteger::MzTabInteger(const int v)
  {
    set(v);
  }

  MzTabInteger::MzTabInteger()
      : value_(0), state_(MZTAB_CELLSTATE_NULL)
  {
  }

  void MzTabInteger::set(const Int& value)
  {
    state_ = MZTAB_CELLSTATE_DEFAULT;
    value_ = value;
  }

  Int MzTabInteger::get() const
  {
    if (state_ == MZTAB_CELLSTATE_DEFAULT)
    {
      return value_;
    }
    else
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,StringUtils::toStr("Trying to extract MzTab Integer value from non-integer valued cell. Did you check the cell state before querying the value?"));
    }
  }

  std::string MzTabInteger::toCellString() const
  {
    switch (state_)
    {
      case MZTAB_CELLSTATE_NULL:
        return "null";

      case MZTAB_CELLSTATE_NAN:
        return "NaN";

      case MZTAB_CELLSTATE_INF:
        return "Inf";

      case MZTAB_CELLSTATE_DEFAULT:
      default:
        return StringUtils::toStr(value_);
    }
  }

  void MzTabInteger::fromCellString(const std::string& s)
  {
    std::string lower = s;
    StringUtils::trim(StringUtils::toLower(lower));
    if (lower == "null") { setNull(true); }
    else if (lower == "nan") { setNaN(); }
    else if (lower == "inf") { setInf(); }
    else // default case
    {
      // some mzTab files from external sources contain floating point numbers in integer columns
      auto val = StringUtils::toDouble(lower);
      if (val != (Int)val) // check if the value is actually an integer (e.g. 4.0)
      {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,StringUtils::toStr("Could not convert String '") + s + "' to MzTabInteger");
      }
      set((Int)val);
    }
  }

  bool MzTabInteger::isNull() const
  {
    return state_ == MZTAB_CELLSTATE_NULL;
  }

  void MzTabInteger::setNull(bool b)
  {
    state_ = b ? MZTAB_CELLSTATE_NULL : MZTAB_CELLSTATE_DEFAULT;
  }

  bool MzTabInteger::isNaN() const
  {
    return state_ == MZTAB_CELLSTATE_NAN;
  }

  void MzTabInteger::setNaN()
  {
    state_ = MZTAB_CELLSTATE_NAN;
  }

  bool MzTabInteger::isInf() const
  {
    return state_ == MZTAB_CELLSTATE_INF;
  }

  void MzTabInteger::setInf()
  {
    state_ = MZTAB_CELLSTATE_INF;
  }

  bool MzTabDouble::isNull() const
  {
    return state_ == MZTAB_CELLSTATE_NULL;
  }

  void MzTabDouble::setNull(bool b)
  {
    state_ = b ? MZTAB_CELLSTATE_NULL : MZTAB_CELLSTATE_DEFAULT;
  }

  bool MzTabDouble::isNaN() const
  {
    return state_ == MZTAB_CELLSTATE_NAN;
  }

  void MzTabDouble::setNaN()
  {
    state_ = MZTAB_CELLSTATE_NAN;
  }

  bool MzTabDouble::isInf() const
  {
    return state_ == MZTAB_CELLSTATE_INF;
  }

  void MzTabDouble::setInf()
  {
    state_ = MZTAB_CELLSTATE_INF;
  }

  MzTabDouble::MzTabDouble()
      : value_(0.0), state_(MZTAB_CELLSTATE_NULL)
  {
  }

  MzTabDouble::MzTabDouble(const double v)
  {
    set(v);
  }

  void MzTabDouble::set(const double& value)
  {
    state_ = MZTAB_CELLSTATE_DEFAULT;
    value_ = value;
  }

  double MzTabDouble::get() const
  {
    if (state_ != MZTAB_CELLSTATE_DEFAULT)
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,StringUtils::toStr("Trying to extract MzTab Double value from non-double valued cell. Did you check the cell state before querying the value?"));
    }

    return value_;
  }

  std::string MzTabDouble::toCellString() const
  {
    switch (state_)
    {
      case MZTAB_CELLSTATE_NULL:
        return "null";

      case MZTAB_CELLSTATE_NAN:
        return "NaN";

      case MZTAB_CELLSTATE_INF:
        return "Inf";

      case MZTAB_CELLSTATE_DEFAULT:
      default:
        return StringUtils::toStr(value_);
    }
  }

  void MzTabDouble::fromCellString(const std::string& s)
  {
    std::string lower = s;
    StringUtils::trim(StringUtils::toLower(lower));
    if (lower == "null")
    {
      setNull(true);
    }
    else if (lower == "nan")
    {
      setNaN();
    }
    else if (lower == "inf")
    {
      setInf();
    }
    else // default case
    {
      set(StringUtils::toDouble(lower));
    }
  }

  bool MzTabDouble::operator<(const MzTabDouble& rhs) const
  {
    return this->value_ < rhs.value_;
  }

  bool MzTabDouble::operator==(const MzTabDouble& rhs) const
  {
    return this->value_ == rhs.value_;
  }

  bool MzTabDoubleList::isNull() const
  {
    return entries_.empty();
  }

  void MzTabDoubleList::setNull(bool b)
  {
    if (b)
    {
      entries_.clear();
    }
  }

  std::string MzTabDoubleList::toCellString() const
  {
    if (isNull())
    {
      return "null";
    }
    else
    {
      std::string ret;
      for (std::vector<MzTabDouble>::const_iterator it = entries_.begin(); it != entries_.end(); ++it)
      {
        if (it != entries_.begin())
        {
          ret += "|";
        }
        ret += it->toCellString();
      }
      return ret;
    }
  }

  void MzTabDoubleList::fromCellString(const std::string& s)
  {
    std::string lower = s;
    StringUtils::trim(StringUtils::toLower(lower));
    if (lower == "null")
    {
      setNull(true);
    }
    else
    {
      std::vector<std::string> fields;
      StringUtils::split(s, "|", fields);
      for (Size i = 0; i != fields.size(); ++i)
      {
        MzTabDouble ds;
        ds.fromCellString(fields[i]);
        entries_.push_back(ds);
      }
    }
  }

  std::vector<MzTabDouble> MzTabDoubleList::get() const
  {
    return entries_;
  }

  void MzTabDoubleList::set(const std::vector<MzTabDouble>& entries)
  {
    entries_ = entries;
  }

} // namespace OpenMS