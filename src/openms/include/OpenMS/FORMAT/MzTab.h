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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MZTAB_H
#define OPENMS_FORMAT_MZTAB_H

#include <OpenMS/FORMAT/SVOutStream.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <map>
#include <vector>
#include <list>
#include <algorithm>
#include <OpenMS/KERNEL/StandardTypes.h>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wnon-virtual-dtor"

namespace OpenMS
{
/**
      @brief Data model of MzTab files.
      Please see the official MzTab specification at https://code.google.com/p/mztab/

      @ingroup FileIO
  */

// MzTab supports null, NaN, Inf for cells with Integer or Double values. MzTabCellType explicitly defines the state of the cell for these types.
enum MzTabCellStateType
{
  MZTAB_CELLSTATE_DEFAULT,
  MZTAB_CELLSTATE_NULL,
  MZTAB_CELLSTATE_NAN,
  MZTAB_CELLSTATE_INF,
  SIZE_OF_MZTAB_CELLTYPE
};

// basic interface for all MzTab datatypes (can be null; are converted from and to cell string)
class MzTabNullAbleInterface
{
  public:
    virtual ~MzTabNullAbleInterface() {}
    virtual bool isNull() const = 0;
    virtual void setNull(bool b) = 0;
    virtual String toCellString() const = 0;
    virtual void fromCellString(const String&) = 0;
};

// interface for NaN- and Inf- able datatypes (Double and Integer in MzTab). These are as well null-able
class MzTabNullNaNAndInfAbleInterface :
    public MzTabNullAbleInterface
{
  public:
    virtual ~MzTabNullNaNAndInfAbleInterface() {}
    virtual bool isNaN() const = 0;
    virtual void setNaN() = 0;
    virtual bool isInf() const = 0;
    virtual void setInf() = 0;
};

// base class for atomic, non-container types (Double, Int)
class MzTabNullAbleBase :
    public MzTabNullAbleInterface
{
  public:
    MzTabNullAbleBase() :
      null_(true)
    {
    }

    virtual ~MzTabNullAbleBase()
    {
    }

    bool isNull() const
    {
      return null_;
    }

    void setNull(bool b)
    {
      null_ = b;
    }

  protected:
    bool null_;
};

// base class for the atomic non-container like MzTab data types (Double, Int)
class MzTabNullNaNAndInfAbleBase :
    public MzTabNullNaNAndInfAbleInterface
{
  public:
    MzTabNullNaNAndInfAbleBase() :
      state_(MZTAB_CELLSTATE_NULL)
    {
    }

    virtual ~MzTabNullNaNAndInfAbleBase() {}

    bool isNull() const
    {
      return state_ == MZTAB_CELLSTATE_NULL;
    }

    void setNull(bool b)
    {
      state_ = b ? MZTAB_CELLSTATE_NULL : MZTAB_CELLSTATE_DEFAULT;
    }

    bool isNaN() const
    {
      return state_ == MZTAB_CELLSTATE_NAN;
    }

    void setNaN()
    {
      state_ = MZTAB_CELLSTATE_NAN;
    }

    bool isInf() const
    {
      return state_ == MZTAB_CELLSTATE_INF;
    }

    void setInf()
    {
      state_ = MZTAB_CELLSTATE_INF;
    }

  protected:
    MzTabCellStateType state_;
};

class MzTabDouble :
    public MzTabNullNaNAndInfAbleBase
{
  public:
    MzTabDouble()
      : value_(0.0)
    {
    }

    virtual ~MzTabDouble() {}

    void set(const double& value)
    {
      state_ = MZTAB_CELLSTATE_DEFAULT;
      value_ = value;
    }
    explicit MzTabDouble(const double v)
    {
      set(v);
    }

    double get() const
    {
      if (state_ != MZTAB_CELLSTATE_DEFAULT)
      {
        throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Trying to extract MzTab Double value from non-double valued cell. Did you check the cell state before querying the value?"));
      }

      return value_;
    }

    String toCellString() const
    {
      switch (state_)
      {
      case MZTAB_CELLSTATE_NULL:
        return String("null");

      case MZTAB_CELLSTATE_NAN:
        return String("NaN");

      case MZTAB_CELLSTATE_INF:
        return String("Inf");

      case MZTAB_CELLSTATE_DEFAULT:
      default:
        return String(value_);
      }
    }

    void fromCellString(const String& s)
    {
      String lower = s;
      lower.toLower().trim();
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
      else   // default case
      {
        set(lower.toDouble());
      }
    }

  protected:
    double value_;
};

class MzTabDoubleList :
    public MzTabNullAbleBase
{
  public:
    MzTabDoubleList()
    {
    }

    virtual ~MzTabDoubleList() {}

    bool isNull() const
    {
      return entries_.empty();
    }

    void setNull(bool b)
    {
      if (b)
      {
        entries_.clear();
      }
    }

    String toCellString() const
    {
      if (isNull())
      {
        return "null";
      }
      else
      {
        String ret;
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

    void fromCellString(const String& s)
    {
      String lower = s;
      lower.toLower().trim();
      if (lower == "null")
      {
        setNull(true);
      }
      else
      {
        String ss = s;
        std::vector<String> fields;
        ss.split("|", fields);
        for (Size i = 0; i != fields.size(); ++i)
        {
          MzTabDouble ds;
          ds.fromCellString(fields[i]);
          entries_.push_back(ds);
        }
      }
    }

    std::vector<MzTabDouble> get() const
    {
      return entries_;
    }

    void set(const std::vector<MzTabDouble>& entries)
    {
      entries_ = entries;
    }

  protected:
    std::vector<MzTabDouble> entries_;
};

class MzTabInteger :
    public MzTabNullNaNAndInfAbleBase
{
  public:

    virtual ~MzTabInteger() {}

    void set(const Int& value)
    {
      state_ = MZTAB_CELLSTATE_DEFAULT;
      value_ = value;
    }

    Int get() const
    {
      if (state_ == MZTAB_CELLSTATE_DEFAULT)
      {
        return value_;
      }
      else
      {
        throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Trying to extract MzTab Integer value from non-integer valued cell. Did you check the cell state before querying the value?"));
      }
    }

    String toCellString() const
    {
      switch (state_)
      {
      case MZTAB_CELLSTATE_NULL:
        return String("null");

      case MZTAB_CELLSTATE_NAN:
        return String("NaN");

      case MZTAB_CELLSTATE_INF:
        return String("Inf");

      case MZTAB_CELLSTATE_DEFAULT:
      default:
        return String(value_);
      }
    }

    void fromCellString(const String& s)
    {
      String lower = s;
      lower.toLower().trim();
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
      else   // default case
      {
        set(lower.toInt());
      }
    }

  protected:
    Int value_;
};

class MzTabIntegerList :
    public MzTabNullAbleBase
{
  public:
    MzTabIntegerList()
    {
    }

    bool isNull() const
    {
      return entries_.empty();
    }

    void setNull(bool b)
    {
      if (b)
      {
        entries_.clear();
      }
    }

    String toCellString() const
    {
      if (isNull())
      {
        return "null";
      }
      else
      {
        String ret;
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

    void fromCellString(const String& s)
    {
      String lower = s;
      lower.toLower().trim();
      if (lower == "null")
      {
        setNull(true);
      }
      else
      {
        String ss = s;
        std::vector<String> fields;
        ss.split(",", fields);
        for (Size i = 0; i != fields.size(); ++i)
        {
          MzTabInteger ds;
          ds.fromCellString(fields[i]);
          entries_.push_back(ds);
        }
      }
    }

    std::vector<MzTabInteger> get() const
    {
      return entries_;
    }

    void set(const std::vector<MzTabInteger>& entries)
    {
      entries_ = entries;
    }

  protected:
    std::vector<MzTabInteger> entries_;
};

class MzTabBoolean :
    public MzTabNullAbleBase
{
  public:
    MzTabBoolean()
      :value_(false)
    {

    }

    explicit MzTabBoolean(bool v)
    {
      set(v);
    }

    virtual ~MzTabBoolean() {}

    void set(const bool& value)
    {
      setNull(false);
      value_ = value;
    }

    Int get() const
    {
      return value_;
    }

    String toCellString() const
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

    void fromCellString(const String& s)
    {
      String lower = s;
      lower.toLower().trim();
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
          throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Could not convert String '") + s + "' to MzTabBoolean");
        }
      }
    }

  protected:
    bool value_;
};

class MzTabString :
    public MzTabNullAbleInterface
{
  public:

    MzTabString()
      :value_()
    {
    }

    explicit MzTabString(const String& s)
    {
      set(s);
    }

    virtual ~MzTabString() {}

    void set(const String& value)
    {
      String lower = value;
      lower.toLower().trim();
      if (lower == "null")
      {
        setNull(true);
      }
      else
      {
        value_ = value;
      }
    }

    String get() const
    {
      return value_;
    }

    bool isNull() const
    {
      return value_.empty();
    }

    void setNull(bool b)
    {
      if (b)
      {
        value_.clear();
      }
    }

    String toCellString() const
    {
      if (isNull())
      {
        return String("null");
      }
      else
      {
        return value_;
      }
    }

    void fromCellString(const String& s)
    {
      set(s);
    }

  protected:
    String value_;
};

class MzTabParameter :
    public MzTabNullAbleInterface
{
  public:

    virtual ~MzTabParameter() {}

    bool isNull() const
    {
      return CV_label_.empty() && accession_.empty() && name_.empty() && value_.empty();
    }

    void setNull(bool b)
    {
      if (b)
      {
        CV_label_.clear();
        accession_.clear();
        name_.clear();
        value_.clear();
      }
    }

    void setCVLabel(const String& CV_label)
    {
      CV_label_ = CV_label;
    }

    void setAccession(const String& accession)
    {
      accession_ = accession;
    }

    void setName(const String& name)
    {
      name_ = name;
    }

    void setValue(const String& value)
    {
      value_ = value;
    }

    String getCVLabel() const
    {
      assert(!isNull());
      return CV_label_;
    }

    String getAccession() const
    {
      assert(!isNull());
      return accession_;
    }

    String getName() const
    {
      assert(!isNull());
      return name_;
    }

    String getValue() const
    {
      assert(!isNull());
      return value_;
    }

    String toCellString() const
    {
      if (isNull())
      {
        return "null";
      }
      else
      {
        String ret = "[";
        ret += CV_label_ + ", ";
        ret += accession_ + ", ";

        if (name_.hasSubstring(", "))
        {
          ret += String("\"") + name_ + String("\""); // quote name if it contains a ","
        } else
        {
          ret += name_;
        }

        ret += String(", ");

        if (value_.hasSubstring(", "))
        {
          ret += String("\"") + value_ + String("\""); // quote value if it contains a ","
        } else
        {
          ret += value_;
        }

        ret += "]";
        return ret;
      }
    }

    void fromCellString(const String& s)
    {
      String lower = s;
      lower.toLower().trim();
      if (lower == "null")
      {
        setNull(true);
      }
      else
      {
        StringList fields;
        String field;
        bool in_quotes = false;
        String::const_iterator quote_start = s.begin();
        for (String::const_iterator sit = s.begin(); sit != s.end(); ++sit)
        {
          if (*sit == '\"')  // start or end of quotes
          {
            in_quotes = !in_quotes;
          } else if (*sit == ',')   // , encountered
          {
            if (in_quotes)  // case 1: , in quote
            {
              field += ','; // add , (no split)
            } else // split at , if not in quotes
            {
              fields.push_back(field.trim());
              field.clear();
            }
          } else if (*sit != '[' && *sit != ']')
          {
            // skip leading ws
            if (*sit == ' ' && field.empty())
            {
              continue;
            }
            field += *sit;
          }
        }

        fields.push_back(field.trim());

        if (fields.size() != 4)
        {
          throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Could not convert String '") + s + "' to MzTabParameter");
        }

        CV_label_ = fields[0];
        accession_ = fields[1];
        name_ = fields[2];
        value_ = fields[3];
      }
    }

  protected:
    String CV_label_;
    String accession_;
    String name_;
    String value_;
};

class MzTabParameterList :
    public MzTabNullAbleInterface
{
  public:

    virtual ~MzTabParameterList() {}

    bool isNull() const
    {
      return parameters_.empty();
    }

    void setNull(bool b)
    {
      if (b)
      {
        parameters_.clear();
      }
    }

    String toCellString() const
    {
      if (isNull())
      {
        return "null";
      }
      else
      {
        String ret;
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

    void fromCellString(const String& s)
    {
      String lower = s;
      lower.toLower().trim();

      if (lower == "null")
      {
        setNull(true);
      }
      else
      {
        String ss = s;
        std::vector<String> fields;
        ss.split("|", fields);
        for (Size i = 0; i != fields.size(); ++i)
        {
          MzTabParameter p;
          lower = fields[i];
          lower.toLower().trim();
          if (lower == "null")
          {
            throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("MzTabParameter in MzTabParameterList must not be null '") + s);
          }
          p.fromCellString(fields[i]);
          parameters_.push_back(p);
        }
      }
    }

    std::vector<MzTabParameter> get() const
    {
      return parameters_;
    }

    void set(const std::vector<MzTabParameter>& parameters)
    {
      parameters_ = parameters;
    }

  protected:
    std::vector<MzTabParameter> parameters_;
};

class MzTabStringList :
    public MzTabNullAbleInterface
{
  public:
    MzTabStringList() :
      sep_('|')
    {
    }

    virtual ~MzTabStringList() {}

    // needed for e.g. ambiguity_members and GO accessions as these use ',' as separator while the others use '|'
    void setSeparator(char sep)
    {
      sep_ = sep;
    }

    bool isNull() const
    {
      return entries_.empty();
    }

    void setNull(bool b)
    {
      if (b)
      {
        entries_.clear();
      }
    }

    String toCellString() const
    {
      if (isNull())
      {
        return "null";
      }
      else
      {
        String ret;
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

    void fromCellString(const String& s)
    {
      String lower = s;
      lower.toLower().trim();

      if (lower == "null")
      {
        setNull(true);
      }
      else
      {
        String ss = s;
        std::vector<String> fields;
        ss.split(sep_, fields);
        for (Size i = 0; i != fields.size(); ++i)
        {
          MzTabString ts;
          ts.fromCellString(fields[i]);
          entries_.push_back(ts);
        }
      }
    }

    std::vector<MzTabString> get() const
    {
      return entries_;
    }

    void set(const std::vector<MzTabString>& entries)
    {
      entries_ = entries;
    }

  protected:
    std::vector<MzTabString> entries_;
    char sep_;
};

struct MzTabModification :
    public MzTabNullAbleInterface
{
  public:

    virtual ~MzTabModification() {}

    bool isNull() const
    {
      return pos_param_pairs_.empty() && mod_identifier_.isNull();
    }

    void setNull(bool b)
    {
      if (b)
      {
        pos_param_pairs_.clear();
        mod_identifier_.setNull(true);
      }
    }

    // set (potentially ambiguous) position(s) with associated parameter (might be null if not set)
    void setPositionsAndParameters(const std::vector<std::pair<Size, MzTabParameter> >& ppp)
    {
      pos_param_pairs_ = ppp;
    }

    std::vector<std::pair<Size, MzTabParameter> > getPositionsAndParameters() const
    {
      return pos_param_pairs_;
    }

    void setModificationIdentifier(const MzTabString& mod_id)
    {
      mod_identifier_ = mod_id;
    }

    MzTabString getModOrSubstIdentifier() const
    {
      assert(!isNull());
      return mod_identifier_;
    }

    String toCellString() const
    {
      if (isNull())
      {
        return String("null");
      }
      else
      {
        String pos_param_string;

        for (Size i = 0; i != pos_param_pairs_.size(); ++i)
        {
          pos_param_string += pos_param_pairs_[i].first;

          // attach MzTabParameter if available
          if (!pos_param_pairs_[i].second.isNull())
          {
            pos_param_string += pos_param_pairs_[i].second.toCellString();
          }

          // add | as separator (except for last one)
          if (i < pos_param_pairs_.size() - 1)
          {
            pos_param_string += String("|");
          }
        }

        // quick sanity check
        if (mod_identifier_.isNull())
        {
          throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Modification or Substitution identifier MUST NOT be null or empty in MzTabModification"));
        }

        String res;
        // only add '-' if we have position information
        if (!pos_param_string.empty())
        {
          res = pos_param_string + "-" + mod_identifier_.toCellString();
        }
        else
        {
          res = mod_identifier_.toCellString();
        }
        return res;
      }
    }

    void fromCellString(const String& s)
    {
      String lower = s;
      lower.toLower().trim();
      if (lower == "null")
      {
        setNull(true);
      }
      else
      {
        if (!lower.hasSubstring("-")) // no positions? simply use s as mod identifier
        {
          mod_identifier_.set(String(s).trim());
        }
        else
        {
          String ss = s;
          ss.trim();
          std::vector<String> fields;
          ss.split("-", fields);

          if (fields.size() != 2)
          {
            throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Can't convert to MzTabModification from '") + s);
          }
          mod_identifier_.fromCellString(fields[1].trim());

          std::vector<String> position_fields;
          fields[0].split("|", position_fields);

          for (Size i = 0; i != position_fields.size(); ++i)
          {
            Size spos = position_fields[i].find_first_of("[");

            if (spos == std::string::npos) // only position information and no parameter
            {
              pos_param_pairs_.push_back(std::make_pair(position_fields[i].toInt(), MzTabParameter()));
            }
            else
            {
              // extract position part
              Int pos = String(position_fields[i].begin(), position_fields[i].begin() + spos).toInt();

              // extract [,,,] part
              MzTabParameter param;
              param.fromCellString(position_fields[i].substr(spos));
              pos_param_pairs_.push_back(std::make_pair(pos, param));
            }
          }
        }
      }
    }

  protected:
    std::vector<std::pair<Size, MzTabParameter> > pos_param_pairs_;
    MzTabString mod_identifier_;
};

class MzTabModificationList :
    public MzTabNullAbleBase
{
  public:
    virtual ~MzTabModificationList() {}

    bool isNull() const
    {
      return entries_.empty();
    }

    void setNull(bool b)
    {
      if (b)
      {
        entries_.clear();
      }
    }

    String toCellString() const
    {
      if (isNull())
      {
        return "null";
      }
      else
      {
        String ret;
        for (std::vector<MzTabModification>::const_iterator it = entries_.begin(); it != entries_.end(); ++it)
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

    void fromCellString(const String& s)
    {
      String lower = s;
      lower.toLower().trim();
      if (lower == "null")
      {
        setNull(true);
      }
      else
      {
        String ss = s;
        std::vector<String> fields;

        if (!ss.hasSubstring("[")) // no parameters
        {
          ss.split(",", fields);
          for (Size i = 0; i != fields.size(); ++i)
          {
            MzTabModification ms;
            ms.fromCellString(fields[i]);
            entries_.push_back(ms);
          }
        }
        else
        {
          // example string: 3|4[a,b,,v]|8[,,"blabla, [bla]",v],1|2|3[a,b,,v]-mod:123
          // we don't want to split at the , inside of [ ]  MzTabParameter brackets.
          // Additionally,  and we don't want to recognise quoted brackets inside the MzTabParameter where they can occur in quoted text (see example string)
          bool in_param_bracket = false;
          bool in_quotes = false;

          for (Size pos = 0; pos != ss.size(); ++pos)
          {
            // param_bracket state
            if (ss[pos] == '[' && !in_quotes)
            {
              in_param_bracket = true;
              continue;
            }

            if (ss[pos] == ']' && !in_quotes)
            {
              in_param_bracket = false;
              continue;
            }

            // quote state
            if (ss[pos] == '\"')
            {
              in_quotes = !in_quotes;
              continue;
            }

            // comma in param bracket
            if (ss[pos] == ',' && !in_quotes && in_param_bracket)
            {
              ss[pos] = ((char)007); // use ASCII bell as temporary separator
              continue;
            }
          }

          // now the split at comma is save
          ss.split(",", fields);

          for (Size i = 0; i != fields.size(); ++i)
          {
            fields[i].substitute(((char)007), ','); // resubstitute comma after split
            MzTabModification ms;
            ms.fromCellString(fields[i]);
            entries_.push_back(ms);
          }
        }
      }
    }

    std::vector<MzTabModification> get() const
    {
      return entries_;
    }

    void set(const std::vector<MzTabModification>& entries)
    {
      entries_ = entries;
    }

  protected:
    std::vector<MzTabModification> entries_;

};

class MzTabSpectraRef :
    public MzTabNullAbleInterface
{
  public:
    MzTabSpectraRef() :
      ms_run_(0)
    {
    }

    virtual ~MzTabSpectraRef() {}

    bool isNull() const
    {
      return (ms_run_ < 1) || (spec_ref_.empty());
    }

    void setNull(bool b)
    {
      if (b)
      {
        ms_run_ = 0;
        spec_ref_.clear();
      }
    }

    void setMSFile(Size index)
    {
      assert(index >= 1);
      if (index >= 1)
      {
        ms_run_ = index;
      }
    }

    void setSpecRef(String spec_ref)
    {
      assert(!spec_ref.empty());
      if (!spec_ref.empty())
      {
        spec_ref_ = spec_ref;
      }
    }

    String getSpecRef() const
    {
      assert(!isNull());
      return spec_ref_;
    }

    Size getMSFile() const
    {
      assert(!isNull());
      return ms_run_;
    }

    void setSpecRefFile(const String& spec_ref)
    {
      assert(!spec_ref.empty());
      if (!spec_ref.empty())
      {
        spec_ref_ = spec_ref;
      }
    }

    String toCellString() const
    {
      if (isNull())
      {
        return String("null");
      }
      else
      {
        return String("ms_run[") + String(ms_run_) + "]:" + spec_ref_;
      }
    }

    void fromCellString(const String& s)
    {
      String lower = s;
      lower.toLower().trim();
      if (lower == "null")
      {
        setNull(true);
      }
      else
      {
        String ss = s;
        std::vector<String> fields;
        ss.split(":", fields);
        if (fields.size() != 2)
        {
          throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Can not convert to MzTabSpectraRef from '") + s);
        }

        spec_ref_ = fields[1];
        ms_run_ = (Size)(fields[0].substitute("ms_run[", "").remove(']').toInt());
      }
    }

  protected:
    Size ms_run_; // number is specified in the meta data section.
    String spec_ref_;
};

// MTD

struct MzTabSampleMetaData
{
    MzTabString description;
    std::map<Size, MzTabParameter> species;
    std::map<Size, MzTabParameter> tissue;
    std::map<Size, MzTabParameter> cell_type;
    std::map<Size, MzTabParameter> disease;
    std::map<Size, MzTabParameter> custom;
};

struct MzTabSoftwareMetaData
{
    MzTabParameter software;
    std::map<Size, MzTabString> setting;
};

struct MzTabModificationMetaData
{
    MzTabParameter modification;
    MzTabString site;
    MzTabString position;
};

struct MzTabAssayMetaData
{
    MzTabParameter quantification_reagent;
    std::map<Size, MzTabModificationMetaData> quantification_mod;
    MzTabString sample_ref;
    MzTabString ms_run_ref;
};

struct MzTabCVMetaData
{
    MzTabString label;
    MzTabString full_name;
    MzTabString version;
    MzTabString url;
};

struct MzTabInstrumentMetaData
{
    MzTabParameter name;
    MzTabParameter source;
    std::map<Size, MzTabParameter> analyzer;
    MzTabParameter detector;
};

struct MzTabContactMetaData
{
    MzTabString name;
    MzTabString affiliation;
    MzTabString email;
};

struct MzTabMSRunMetaData
{
    MzTabParameter format;
    MzTabString location;
    MzTabParameter id_format;
    MzTabParameterList fragmentation_method;
};

struct MzTabStudyVariableMetaData
{
    MzTabIntegerList assay_refs;
    MzTabIntegerList sample_refs;
    MzTabString description;
};

// all meta data of a mzTab file. Please refer to specification for documentation.
struct MzTabMetaData
{
    MzTabMetaData()
    {
      mz_tab_version.fromCellString(String("1.0.0"));
    }

    MzTabString mz_tab_version;
    MzTabString mz_tab_mode;
    MzTabString mz_tab_type;
    MzTabString mz_tab_id;
    MzTabString title;
    MzTabString description;

    std::map<Size, MzTabParameter> protein_search_engine_score;
    std::map<Size, MzTabParameter> peptide_search_engine_score;
    std::map<Size, MzTabParameter> psm_search_engine_score;
    std::map<Size, MzTabParameter> smallmolecule_search_engine_score;

    std::map<Size, MzTabParameterList> sample_processing;

    std::map<Size, MzTabInstrumentMetaData> instrument;

    std::map<Size, MzTabSoftwareMetaData> software;

    MzTabParameterList false_discovery_rate;

    std::map<Size, MzTabString> publication;

    std::map<Size, MzTabContactMetaData> contact;

    std::map<Size, MzTabString> uri;

    std::map<Size, MzTabModificationMetaData> fixed_mod;

    std::map<Size, MzTabModificationMetaData> variable_mod;

    MzTabParameter quantification_method;

    MzTabParameter protein_quantification_unit;
    MzTabParameter peptide_quantification_unit;
    MzTabParameter small_molecule_quantification_unit;

    std::map<Size, MzTabMSRunMetaData> ms_run;

    std::map<Size, MzTabParameter> custom;

    std::map<Size, MzTabSampleMetaData> sample;

    std::map<Size, MzTabAssayMetaData> assay;

    std::map<Size, MzTabStudyVariableMetaData> study_variable;

    std::map<Size, MzTabCVMetaData> cv;

    std::vector<String> colunit_protein;
    std::vector<String> colunit_peptide;
    std::vector<String> colunit_psm;
    std::vector<String> colunit_small_molecule;
};

typedef std::pair<String, MzTabString> MzTabOptionalColumnEntry; //  column name (not null able), value (null able)

// PRT - Protein section (Table based)
struct MzTabProteinSectionRow
{
    MzTabProteinSectionRow()
    {
      // use "," as list separator because "|" can be used for go terms and protein accessions
      go_terms.setSeparator(',');
      ambiguity_members.setSeparator(',');
    }

    MzTabString accession; // The protein’s accession.
    MzTabString description; // Human readable description (i.e. the name)
    MzTabInteger taxid; // NEWT taxonomy for the species.
    MzTabString species; // Human readable name of the species
    MzTabString database; // Name of the protein database.
    MzTabString database_version; // String Version of the protein database.
    MzTabParameterList search_engine; // Search engine(s) identifying the protein.
    std::map<Size, MzTabDouble>  best_search_engine_score;  // best_search_engine_score[1-n]
    std::map<Size, std::map<Size, MzTabDouble> > search_engine_score_ms_run; // search_engine_score[index1]_ms_run[index2]
    MzTabInteger reliability;
    std::map<Size, MzTabInteger> num_psms_ms_run;
    std::map<Size, MzTabInteger> num_peptides_distinct_ms_run;
    std::map<Size, MzTabInteger> num_peptides_unique_ms_run;
    MzTabStringList ambiguity_members; // Alternative protein identifications.
    MzTabModificationList modifications; // Modifications identified in the protein.
    MzTabString uri; // Location of the protein’s source entry.
    MzTabStringList go_terms; // List of GO terms for the protein.
    MzTabDouble protein_coverage; // (0-1) Amount of protein sequence identified.
    std::map<Size, MzTabDouble> protein_abundance_assay;
    std::map<Size, MzTabDouble> protein_abundance_study_variable;
    std::map<Size, MzTabDouble> protein_abundance_stdev_study_variable;
    std::map<Size, MzTabDouble> protein_abundance_std_error_study_variable;
    std::vector<MzTabOptionalColumnEntry> opt_; // Optional Columns must start with “opt_”
};

// PEP - Peptide section (Table based)
struct MzTabPeptideSectionRow
{
    MzTabString sequence; // The peptide’s sequence.
    MzTabString accession; // The protein’s accession.
    MzTabBoolean unique; // 0=false, 1=true, null else: Peptide is unique for the protein.
    MzTabString database; // Name of the sequence database.
    MzTabString database_version; // Version (and optionally # of entries).
    MzTabParameterList search_engine; // Search engine(s) that identified the peptide.
    std::map<Size, MzTabDouble> best_search_engine_score; // Search engine(s) score(s) for the peptide.
    std::map<Size, std::map<Size, MzTabDouble> > search_engine_score_ms_run;
    MzTabInteger reliability; // (1-3) 0=null Identification reliability for the peptide.
    MzTabModificationList modifications; // Modifications identified in the peptide.
    MzTabDoubleList retention_time; // Time points in seconds. Semantics may vary.
    MzTabDoubleList retention_time_window;
    MzTabDouble charge; // Precursor ion’s charge.
    MzTabDouble mass_to_charge; // Precursor ion’s m/z.
    MzTabString uri; // Location of the PSMs source entry.
    MzTabSpectraRef spectra_ref; // Spectra identifying the peptide.
    std::map<Size, MzTabDouble> peptide_abundance_assay;
    std::map<Size, MzTabDouble> peptide_abundance_study_variable;
    std::map<Size, MzTabDouble> peptide_abundance_stdev_study_variable;
    std::map<Size, MzTabDouble> peptide_abundance_std_error_study_variable;
    std::vector<MzTabOptionalColumnEntry> opt_; // Optional columns must start with “opt_”.
};

// PSM - PSM section (Table based)
struct MzTabPSMSectionRow
{
    MzTabString sequence; // The peptide’s sequence.
    MzTabInteger PSM_ID;
    MzTabString accession; // The protein’s accession.
    MzTabBoolean unique; // 0=false, 1=true, null else: Peptide is unique for the protein.
    MzTabString database; // Name of the sequence database.
    MzTabString database_version; // Version (and optionally # of entries).
    MzTabParameterList search_engine; // Search engine(s) that identified the peptide.
    std::map<Size, MzTabDouble> search_engine_score; // Search engine(s) score(s) for the peptide.
    MzTabInteger reliability; // (1-3) 0=null Identification reliability for the peptide.
    MzTabModificationList modifications; // Modifications identified in the peptide.
    MzTabDoubleList retention_time; // Time points in seconds. Semantics may vary.
    MzTabDouble charge; // Precursor ion’s charge.
    MzTabDouble exp_mass_to_charge; // Precursor ion’s m/z.
    MzTabDouble calc_mass_to_charge;
    MzTabString uri; // Location of the PSM’s source entry.
    MzTabSpectraRef spectra_ref; // Spectra identifying the peptide.
    MzTabString pre;
    MzTabString post;
    MzTabString start;
    MzTabString end;
    std::vector<MzTabOptionalColumnEntry> opt_; // Optional columns must start with “opt_”.
};

// SML Small molecule section (table based)
struct MzTabSmallMoleculeSectionRow
{
    MzTabStringList identifier; // The small molecule’s identifier.
    MzTabString chemical_formula; // Chemical formula of the identified compound.
    MzTabString smiles; // Molecular structure in SMILES format.
    MzTabString inchi_key; // InChi Key of the identified compound.
    MzTabString description; // Human readable description (i.e. the name)
    MzTabDouble exp_mass_to_charge; // Precursor ion’s m/z.
    MzTabDouble calc_mass_to_charge; // Precursor ion’s m/z.
    MzTabDouble charge; // Precursor ion’s charge.
    MzTabDoubleList retention_time; // Time points in seconds. Semantics may vary.
    MzTabInteger taxid; // NEWT taxonomy for the species.
    MzTabString species; // Human readable name of the species
    MzTabString database; // Name of the used database.
    MzTabString database_version; // String Version of the database (and optionally # of compounds).
    MzTabInteger reliability; // (1-3) The identification reliability.
    MzTabString uri; // The source entry’s location.
    MzTabSpectraRef spectra_ref; // Spectra identifying the small molecule.
    MzTabParameterList search_engine; // Search engine(s) identifying the small molecule.
    std::map<Size, MzTabDouble> best_search_engine_score; // Search engine(s) identifications score(s).
    std::map<Size, std::map<Size, MzTabDouble> > search_engine_score_ms_run;
    MzTabString modifications; // Modifications identified on the small molecule.
    std::map<Size, MzTabDouble> smallmolecule_abundance_assay;
    std::map<Size, MzTabDouble> smallmolecule_abundance_study_variable;
    std::map<Size, MzTabDouble> smallmolecule_abundance_stdev_study_variable;
    std::map<Size, MzTabDouble> smallmolecule_abundance_std_error_study_variable;
    std::vector<MzTabOptionalColumnEntry> opt_; // Optional columns must start with “opt_”.
};

typedef std::vector<MzTabProteinSectionRow> MzTabProteinSectionRows;
typedef std::vector<MzTabPeptideSectionRow> MzTabPeptideSectionRows;
typedef std::vector<MzTabPSMSectionRow> MzTabPSMSectionRows;
typedef std::vector<MzTabSmallMoleculeSectionRow> MzTabSmallMoleculeSectionRows;

/**
      @brief Data model of MzTab files.
      Please see the official MzTab specification at https://code.google.com/p/mztab/

      @ingroup FileIO
 */
class OPENMS_DLLAPI MzTab
{
  public:
    /// Default constructor
    MzTab() {}

    /// Destructor
    ~MzTab() {}

    const MzTabMetaData& getMetaData() const
    {
      return meta_data_;
    }

    void setMetaData(const MzTabMetaData& md)
    {
      meta_data_ = md;
    }

    const MzTabProteinSectionRows& getProteinSectionRows() const
    {
      return protein_data_;
    }

    void setProteinSectionRows(const MzTabProteinSectionRows& psd)
    {
      protein_data_ = psd;
    }

    const MzTabPeptideSectionRows& getPeptideSectionRows() const
    {
      return peptide_data_;
    }

    void setPeptideSectionRows(const MzTabPeptideSectionRows& psd)
    {
      peptide_data_ = psd;
    }

    const MzTabPSMSectionRows& getPSMSectionRows() const
    {
      return psm_data_;
    }

    void setPSMSectionRows(const MzTabPSMSectionRows& psd)
    {
      psm_data_ = psd;
    }

    void setCommentRows(const std::map<Size, String> & com)
    {
      comment_rows_ = com;
    }

    void setEmptyRows(const std::vector<Size> & empty)
    {
      empty_rows_ = empty;
    }

    const std::vector<Size>& getEmptyRows() const
    {
      return empty_rows_;
    }

    const std::map<Size, String>& getCommentRows() const
    {
      return comment_rows_;
    }

    const MzTabSmallMoleculeSectionRows& getSmallMoleculeSectionRows() const
    {
      return small_molecule_data_;
    }

    void setSmallMoleculeSectionRows(const MzTabSmallMoleculeSectionRows& smsd)
    {
      small_molecule_data_ = smsd;
    }

    // Extract opt_ (custom, optional column names)
    std::vector<String> getProteinOptionalColumnNames() const
    {
      std::vector<String> names;
      if (!protein_data_.empty())
      {
        const std::vector<MzTabOptionalColumnEntry>& opt_ = protein_data_[0].opt_;
        for (std::vector<MzTabOptionalColumnEntry>::const_iterator it = opt_.begin(); it != opt_.end(); ++it)
        {
          names.push_back(it->first);
        }
      }
      return names;
    }

    // Extract opt_ (custom, optional column names)
    std::vector<String> getPeptideOptionalColumnNames() const
    {
      std::vector<String> names;
      if (!peptide_data_.empty())
      {
        const std::vector<MzTabOptionalColumnEntry>& opt_ = peptide_data_[0].opt_;
        for (std::vector<MzTabOptionalColumnEntry>::const_iterator it = opt_.begin(); it != opt_.end(); ++it)
        {
          names.push_back(it->first);
        }
      }
      return names;
    }

    // Extract opt_ (custom, optional column names)
    std::vector<String> getPSMOptionalColumnNames() const
    {
      std::vector<String> names;
      if (!psm_data_.empty())
      {
        const std::vector<MzTabOptionalColumnEntry>& opt_ = psm_data_[0].opt_;
        for (std::vector<MzTabOptionalColumnEntry>::const_iterator it = opt_.begin(); it != opt_.end(); ++it)
        {
          names.push_back(it->first);
        }
      }
      return names;
    }

    // Extract opt_ (custom, optional column names)
    std::vector<String> getSmallMoleculeOptionalColumnNames() const
    {
      std::vector<String> names;
      if (!small_molecule_data_.empty())
      {
        const std::vector<MzTabOptionalColumnEntry>& opt_ = small_molecule_data_[0].opt_;
        for (std::vector<MzTabOptionalColumnEntry>::const_iterator it = opt_.begin(); it != opt_.end(); ++it)
        {
          names.push_back(it->first);
        }
      }
      return names;
    }

  protected:
    MzTabMetaData meta_data_;
    MzTabProteinSectionRows protein_data_;
    MzTabPeptideSectionRows peptide_data_;
    MzTabPSMSectionRows psm_data_;
    MzTabSmallMoleculeSectionRows small_molecule_data_;
    std::vector<Size> empty_rows_;  // index of empty rows
    std::map<Size, String> comment_rows_; // comments
};

} // namespace OpenMS

#pragma clang diagnostic pop

#endif // OPENMS_FORMAT_MZTAB_H
