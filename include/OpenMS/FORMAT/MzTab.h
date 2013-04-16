// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
#include <OpenMS/DATASTRUCTURES/StringList.h>

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
    virtual bool isNull() const = 0;
    virtual void setNull(bool b) = 0;
    virtual String toCellString() const = 0;
    virtual void fromCellString(String) = 0;
  };

  // interface for NaN- and Inf- able datatypes (Double and Integer in MzTab). These are as well null-able
  class MzTabNullNaNAndInfAbleInterface : public MzTabNullAbleInterface
  {
  public:
    virtual bool isNaN() const = 0;
    virtual void setNaN() = 0;
    virtual bool isInf() const = 0;
    virtual void setInf() = 0;
  };

  // base class for the atomic non-container like MzTab data types (Double, Int)
  class MzTabNullAbleBase : public MzTabNullAbleInterface
  {
  public:
    MzTabNullAbleBase():
        null_(true)
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
  class MzTabNullNaNAndInfAbleBase : public MzTabNullNaNAndInfAbleInterface
  {
  public:
    MzTabNullNaNAndInfAbleBase():
      state_(MZTAB_CELLSTATE_NULL)
    {
    }
        
    bool isNull() const
    {
      return (state_ == MZTAB_CELLSTATE_NULL);
    }

    void setNull(bool b)
    {
      state_ = b ? MZTAB_CELLSTATE_NULL : MZTAB_CELLSTATE_DEFAULT;
    }

    bool isNaN() const
    {
      return (state_ == MZTAB_CELLSTATE_NAN);
    }

    void setNaN()
    {
      state_ = MZTAB_CELLSTATE_NAN;
    }

    bool isInf() const
    {
      return (state_ == MZTAB_CELLSTATE_INF);
    }

    void setInf()
    {
      state_ = MZTAB_CELLSTATE_INF;
    }

  protected:
    MzTabCellStateType state_;
  };

  class MzTabDouble : public MzTabNullNaNAndInfAbleBase
  {
  public:
    void set(const DoubleReal& value)
    {
      state_ = MZTAB_CELLSTATE_DEFAULT;
      value_ = value;
    }

    DoubleReal get() const
    {
      if (state_ == MZTAB_CELLSTATE_DEFAULT)
      {
        return value_;
      } else
      {
        throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Trying to extract MzTab Double value from non-double valued cell. Did you check the cell state before querying the value?"));
        return 0;
      }
    }

    String toCellString() const
    {
      switch(state_)
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

    void fromCellString(String s)
    {
      s.trim();
      if (s.toUpper() == "null")
      {
        setNull(true);
      } else if (s.toUpper() == "NaN")
      {
        setNaN();
      } else if (s.toUpper() == "Inf")
      {
        setInf();
      } else // default case
      {
        value_ = s.toDouble();
      }
    }

  protected:
    DoubleReal value_;
  };

  class MzTabDoubleList : public MzTabNullAbleBase
  {
  public:
    MzTabDoubleList()
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
      } else
      {
        String ret;
        for (std::vector<MzTabDouble>::const_iterator it = entries_.begin(); it != entries_.end(); ++it)
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

    void fromCellString(String s)
    {
      s.trim();
      if (s.toUpper() == "null")
      {
        setNull(true);
      } else
      {
        std::vector<String> fields;
        s.split(",", fields);
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

  class MzTabInteger : public MzTabNullNaNAndInfAbleBase
  {
  public:
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
      } else
      {
        throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Trying to extract MzTab Integer value from non-integer valued cell. Did you check the cell state before querying the value?"));
        return 0;
      }
    }

    String toCellString() const
    {
      switch(state_)
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

    void fromCellString(String s)
    {
      s.trim();
      if (s.toUpper() == "null")
      {
        setNull(true);
      } else if (s.toUpper() == "NaN")
      {
        setNaN();
      } else if (s.toUpper() == "Inf")
      {
        setInf();
      } else // default case
      {
        value_ = s.toInt();
      }
    }

  protected:
    Int value_;
  };

  class MzTabBoolean : public MzTabNullAbleBase
  {
  public:
    void set(const bool& value)
    {
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
      } else
      {
        if (value_)
        {
          return "1";
        } else
        {
          return "0";
        }
      }
    }

    void fromCellString(String s)
    {
      s.trim();
      if (s.toUpper() == "null")
      {
        setNull(true);
      } else
      {
        if (s == "0")
        {
          value_ = false;
        } else if (s == "1")
        {
          value_ = true;
        } else
        {
          throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Could not convert String '") +s + "' to MzTabBoolean");
        }
      }
    }

  protected:
    bool value_;
  };

  class MzTabString : public MzTabNullAbleInterface
  {
  public:
    void set(const String& value)
    {
      value_ = value;
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
        return "null";
      } else
      {
        return value_;
      }
    }

    void fromCellString(String s)
    {
      s.trim();
      if (s.toUpper() == "null")
      {
        setNull(true);
      } else
      {
        value_ = s;
      }
    }

  protected:
    String value_;
  };

  class MzTabParameter : public MzTabNullAbleInterface
  {
  public:
    bool isNull() const
    {
      return (CV_label_.empty() && accession_.empty() && name_.empty() && value_.empty());
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
      } else
      {
        String ret = "[";
        ret += CV_label_ + ",";
        ret += accession_ + ",";
        ret += name_ + ",";
        ret += value_;
        ret += "]";
        return ret;
      }
    }

    void fromCellString(String s)
    {
      s.trim();
      if (s.toUpper() == "null")
      {
        setNull(true);
      } else
      {
        std::vector<String> fields;
        s.split(',', fields);
        if (fields.size() != 4)
        {
          throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Could not convert String '") + s + "' to MzTabParameter");
        } else
        {
          CV_label_ = fields[0];
          accession_ = fields[1];
          name_ = fields[2];
          value_ = fields[3];
        }
      }
    }

  protected:
    String CV_label_;
    String accession_;
    String name_;
    String value_;
  };

  class MzTabParameterList : public MzTabNullAbleInterface
  {
  public:
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
      } else
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

    void fromCellString(String s)
    {
      s.trim();
      if (s.toUpper() == "null")
      {
        setNull(true);
      } else
      {
        std::vector<String> fields;
        s.split('|', fields);
        for (Size i = 0; i != fields.size(); ++i)
        {
          MzTabParameter p;
          if (fields[i].toUpper() == "null")
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

  class MzTabStringList : public MzTabNullAbleInterface
  {
  public:
    MzTabStringList():
        sep_('|')
    {
    }

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
      } else
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

    void fromCellString(String s)
    {
      s.trim();
      if (s.toUpper() == "null")
      {
        setNull(true);
      } else
      {
        std::vector<String> fields;
        s.split(sep_, fields);
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

  struct MzTabModification : public MzTabNullAbleInterface
  {
  public:
    MzTabModification():
        position_(-1),
        reliability_score_(-1)
    {
    }

    bool isNull() const
    {
      return ((position_ == -1) && (reliability_score_ == -1) && (mod_identifier_.empty()));
    }

    void setNull(bool b)
    {
      if (b)
      {
        position_ = -1;
        reliability_score_ = -1;
        mod_identifier_.clear();
      }
    }

    void setPosition(Int index)
    {
      position_ = index;
    }

    void setReliabilityScore(DoubleReal score)
    {
      reliability_score_ = score;
    }

    void setModIdentifier(String mod_id)
    {
      mod_identifier_ = mod_id;
    }

    Int getPosition() const
    {
      assert(!isNull());
      return position_;
    }

    DoubleReal getReliabilityScore() const
    {
      assert(!isNull());
      return reliability_score_;
    }

    String getModIdentifier() const
    {
      assert(!isNull());
      return mod_identifier_;
    }

    String toCellString() const
    {
      if (isNull())
      {
        return String("null");
      } else
      {
        String position_string;
        if (position_  >= 0)
        {
          position_string = String(position_);
        }

        String reliablility_string;
        if (reliability_score_  >= -1e-10)
        {
          reliablility_string = String("[") + String(reliability_score_) + String("]");
        }

        String res;
        if (!(position_string.empty() && reliablility_string.empty()))
        {
          res = position_string + reliablility_string + "-" + mod_identifier_;
        } else
        {
          res = mod_identifier_;
        }
        return res;
      }
    }

    void fromCellString(String s)
    {
      s.trim();
      if (s.toUpper() == "null")
      {
        setNull(true);
      } else
      {
        if (!s.hasSubstring("-"))  // no position or reliability fields? simply use s as mod identifier
        {
          position_ = -1;
          reliability_score_ = -1;
          mod_identifier_ = s;
        } else
        {
          std::vector<String> fields;
          s.split("-", fields);
          if (fields.size() != 2)
          {
            throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Can't convert to MzTabModification from '") + s);
          }
          mod_identifier_ = fields[1];

          Size spos = fields[0].find_first_of("[");
          if (spos == std::string::npos)  // only position information
          {
            position_ = fields[0].toInt();
          } else
          {
            String pos_string(fields[0].begin(), fields[0].begin() + spos);
            String rel_string(fields[0].begin() + spos, fields[0].end());

            if (pos_string.empty())
            {
              position_ = -1;
            } else
            {
              position_ = pos_string.toInt();
            }

            if (rel_string.empty())
            {
              reliability_score_ = -1;
            } else
            {
              rel_string.remove('[').remove(']');
              reliability_score_ = rel_string.toDouble();
            }
          }
        }
      }
    }

  protected:
    Int position_;
    DoubleReal reliability_score_;
    String mod_identifier_;
  };

  class MzTabModificationList : public MzTabNullAbleBase
  {
  public:
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
      } else
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

    void fromCellString(String s)
    {
      s.trim();
      if (s.toUpper() == "null")
      {
        setNull(true);
      } else
      {
        std::vector<String> fields;
        s.split(",", fields);
        for (Size i = 0; i != fields.size(); ++i)
        {
          MzTabModification ms;
          ms.fromCellString(fields[i]);
          entries_.push_back(ms);
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

  class MzTabSpectraRef : public MzTabNullAbleInterface
  {
  public:
    MzTabSpectraRef():
        ms_file_(0)
    {
    }

    bool isNull() const
    {
      return ((ms_file_ < 1) || (spec_ref_.empty()));
    }

    void setNull(bool b)
    {
      if (b)
      {
        ms_file_ = 0;
        spec_ref_.clear();
      }
    }

    void setMSFile(Size index)
    {
      assert(index >= 1);
      if (index >= 1)
      {
        ms_file_ = index;
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
      return ms_file_;
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
      } else
      {
        return String("ms_file[") + String(ms_file_) + "]:" + spec_ref_;
      }
    }

    void fromCellString(String s)
    {
      s.trim();
      if (s.toUpper() == "null")
      {
        setNull(true);
      } else
      {
        std::vector<String> fields;
        s.split(':', fields);
        if (fields.size() != 2)
        {
          throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Can not convert to MzTabSpectraRef from '") + s);
        }

        spec_ref_ = fields[1];
        ms_file_= (Size)(fields[0].substitute("ms_file[", "").remove(']').toInt());
      }
    }

  protected:
    Size ms_file_; // number is specified in the meta data section.
    String spec_ref_;
  };

  // MTD - Metadata section (Key-value)

  // all meta data belonging to a (potentially empty) sub unit id
  struct MzTabSubIdMetaData
  {
    // ranges denote multiplicity as specified in the specification document
    String sub_id; // 0..1 empty string if no sub samples have been recorded and SUB_ID is optional
    std::vector<MzTabParameter> species; // 0..* Species of the unit / subsample.
    std::vector<MzTabParameter> tissue; // 0..* Tissue of the unit / subsample.
    std::vector<MzTabParameter> cell_type; // 0..* Parameter  Cell type of the unit / subsample.
    std::vector<MzTabParameter> disease; // 0..* Disease state of the unit / subsample.
    std::vector<String> description; // 0..* Description of the subsample.
    std::vector<MzTabParameter> quantification_reagent; // 0..* Quantification reagent used to label the subsample.
    std::vector<MzTabParameter> custom; // 0..* Additional parameters for the subsample.
  };

  // all meta data belonging to one unit id
  struct MzTabUnitIdMetaData
  {
    MzTabUnitIdMetaData():
        title("null"),
        description("null")
    {
    }

    //String unit_id; // the unit id not null able!
    String title; // 0..1 The unit’s title
    String description; // 0..1
    std::vector<MzTabParameterList> sample_processing; // 0..* Description of the sample processing.
    std::vector<MzTabParameter> instrument_name; // 0..* The instrument’s name
    std::vector<MzTabParameter> instrument_source; // 0..* The instrument’s source
    std::vector<MzTabParameter> instrument_analyzer; // 0..* The instrument’s analyzer
    std::vector<MzTabParameter> instrument_detector; // 0..* The instrument’s detector
    std::vector<MzTabParameter> software; // 0..* Analysis software used in the order it was used.
    std::vector<std::vector<String> > software_setting; // 0..* A sotware setting used. This field MAY occur multiple times for a single software (=same index). 
    std::vector<MzTabParameterList> false_discovery_rate; // 0..* False discovery rate(s)for the experiment.
    std::vector<MzTabStringList> publication; // 0..* Publication ids (pubmed / doi).
    std::vector<MzTabString> contact_name; // 0..* Contact name.
    std::vector<MzTabString> contact_affiliation; // 0..* Contact affiliation.
    std::vector<MzTabString> contact_email; // 0..* Contact’s e-mail address.
    std::vector<String> uri; // 0..* Points to the unit’s source data.
    MzTabParameterList mod; // 0..1 Modifications reported in the unit.
    MzTabParameter quantification_method; // 0..1 Quantification method used.
    MzTabParameter protein_quantification_unit; // 0..1 Unit of protein quantification results.
    MzTabParameter peptide_quantification_unit; // 0..1 Unit of peptide quantification results.
    MzTabParameter small_molecule_quantification_unit; // 0..1 Unit of small molecule quantification results.
    std::vector<MzTabParameter> ms_file_format; // // 0..* Data format of the external MS data file.
    std::vector<MzTabParameter> ms_file_location; // 0..* Location of the external MS data file.
    std::vector<MzTabParameter> ms_file_id_format; // 0..* Identifier format of the external MS data file.
    std::vector<MzTabParameter> custom; // 0..*  Additional parameters.
    std::vector<MzTabSubIdMetaData> sub_id_data; // can contain none, one or multiple sub ids

    // Units: The format of the value has to be {column name}={Parameter defining the unit}
    // This field MUST NOT be used to define a unit for quantification columns. 
    std::vector<String> colunit_protein; // 0..* Defines the used unit for a column in the protein section.
    std::vector<String> colunit_peptide; // 0..* Defines the used unit for a column in the peptide section.
    std::vector<String> colunit_small_molecule; // 0..* Defines the used unit for a column in the small molecule section.
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
    //String unit_id; // The unit’s id. not null able!
    MzTabString description; // Human readable description (i.e. the name)
    MzTabInteger taxid; // NEWT taxonomy for the species.
    MzTabString species; // Human readable name of the species
    MzTabString database; // Name of the protein database.
    MzTabString database_version; // String Version of the protein database.
    MzTabParameterList search_engine; // Search engine(s) identifying the protein.
    MzTabParameterList search_engine_score; // Search engine(s) reliability score(s).
    MzTabInteger reliability; // (1-3) Identification reliability.
    MzTabInteger num_peptides; // Number of PSMs assigned to the protein.
    MzTabInteger num_peptides_distinct; // Distinct (sequence + modifications) # of peptides.
    MzTabInteger num_peptides_unambiguous; // Distinct number of unambiguous peptides.
    MzTabStringList ambiguity_members; // Alternative protein identifications.
    MzTabModificationList modifications; // Modifications identified in the protein.
    MzTabString uri; // Location of the protein’s source entry.
    MzTabStringList go_terms; // List of GO terms for the protein.
    MzTabDouble protein_coverage; // (0-1) Amount of protein sequence identified.
    std::vector<DoubleReal> protein_abundance_sub; // Protein abundance in the subsample.
    std::vector<DoubleReal> protein_abundance_stdev_sub; // Standard deviation of the protein abundance.
    std::vector<DoubleReal> protein_abundance_std_error_sub; // Standard error of the protein abundance.
    std::vector<MzTabOptionalColumnEntry> opt_; // Optional Columns must start with “opt_”
  };

  // PEP - Peptide section (Table based)
  struct MzTabPeptideSectionRow
  {
    MzTabString sequence; // The peptide’s sequence.
    MzTabString accession; // The protein’s accession.
    //String unit_id; // The unit’s id.
    MzTabBoolean unique; // 0=false, 1=true, null else: Peptide is unique for the protein.
    MzTabString database; // Name of the sequence database.
    MzTabString database_version; // Version (and optionally # of entries).
    MzTabParameterList search_engine; // Search engine(s) that identified the peptide.
    MzTabParameterList search_engine_score; // Search engine(s) score(s) for the peptide.
    MzTabInteger reliability; // (1-3) 0=null Identification reliability for the peptide.
    MzTabModificationList modifications; // Modifications identified in the peptide.
    MzTabDoubleList retention_time; // Time points in seconds. Semantics may vary.
    MzTabDouble charge; // Precursor ion’s charge.
    MzTabDouble mass_to_charge; // Precursor ion’s m/z.
    MzTabString uri; // Location of the PSM’s source entry.
    MzTabSpectraRef spectra_ref; // Spectra identifying the peptide.
    std::vector<DoubleReal> peptide_abundance_sub; // Peptide abundance in the subsample;
    std::vector<DoubleReal> peptide_abundance_stdev_sub; // Peptide abundance standard deviation.
    std::vector<DoubleReal> peptide_abundance_std_error_sub; // Peptide abundance standard error.
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
    MzTabDouble mass_to_charge; // Precursor ion’s m/z.
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
    MzTabParameterList search_engine_score; // Search engine(s) identifications score(s).
    MzTabModificationList modifications; // Modifications identified on the small molecule.
    std::vector<MzTabDouble> smallmolecule_abundance_sub; // Abundance in the subsample;
    std::vector<MzTabDouble> smallmolecule_abundance_stdev_sub; // Standard deviation of the abundance.
    std::vector<MzTabDouble> smallmolecule_abundance_std_error_sub; // Standard errpr of the abundance.
    std::vector<MzTabOptionalColumnEntry> opt_; // Optional columns must start with “opt_”.
  };

  typedef std::vector<MzTabProteinSectionRow> MzTabProteinSectionRows;

  typedef std::vector<MzTabPeptideSectionRow> MzTabPeptideSectionRows;

  typedef std::vector<MzTabSmallMoleculeSectionRow> MzTabSmallMoleculeSectionRows;

  typedef std::map<String, MzTabUnitIdMetaData> MzTabMetaData;

  typedef std::map<String, MzTabProteinSectionRows> MzTabProteinSectionData;

  typedef std::map<String, MzTabPeptideSectionRows> MzTabPeptideSectionData;

  typedef std::map<String, MzTabSmallMoleculeSectionRows> MzTabSmallMoleculeSectionData;

  /**
      @brief Data model of MzTab files.
      Please see the official MzTab specification at https://code.google.com/p/mztab/

      @ingroup FileIO
 */  
  class OPENMS_DLLAPI MzTab
  {
  public:
    /// Default constructor
    MzTab();

    /// Destructor
    ~MzTab();

    const MzTabMetaData& getMetaData() const
    {
      return map_unitid_to_meta_data_;
    }

    const MzTabProteinSectionData& getProteinSectionData() const
    {
      return map_unitid_to_protein_data_;
    }

    const MzTabPeptideSectionData& getPeptideSectionData() const
    {
      return map_unitid_to_peptide_data_;
    }

    const MzTabSmallMoleculeSectionData& getSmallMoleculeSectionData() const
    {
      return map_unitid_to_small_molecule_data_;
    }

  protected:
    MzTabMetaData map_unitid_to_meta_data_;
    MzTabProteinSectionData map_unitid_to_protein_data_;
    MzTabPeptideSectionData map_unitid_to_peptide_data_;
    MzTabSmallMoleculeSectionData map_unitid_to_small_molecule_data_;
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_MZTAB_H
