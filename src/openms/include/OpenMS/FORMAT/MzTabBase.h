// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Oliver Alka $
// $Authors: Timo Sachsenberg, Oliver Alka $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>


#include <map>
#include <vector>
#include <list>
#include <algorithm>


namespace OpenMS
{
  /**
      @brief Base functionality to for MzTab data models

      @ingroup FileIO
  */

  /// MzTab supports null, NaN, Inf for cells with Integer or Double values. MzTabCellType explicitly defines the state of the cell for these types.
  enum MzTabCellStateType
  {
    MZTAB_CELLSTATE_DEFAULT,
    MZTAB_CELLSTATE_NULL,
    MZTAB_CELLSTATE_NAN,
    MZTAB_CELLSTATE_INF,
    SIZE_OF_MZTAB_CELLTYPE
  };

  class OPENMS_DLLAPI MzTabDouble
  {
public:
    MzTabDouble();

    explicit MzTabDouble(const double v);

    void set(const double& value);

    double get() const;

    String toCellString() const;

    void fromCellString(const String& s);

    bool isNull() const;

    void setNull(bool b);

    bool isNaN() const;

    void setNaN();

    bool isInf() const;

    void setInf();

    ~MzTabDouble() = default;

    bool operator<(const MzTabDouble& rhs) const;

    bool operator==(const MzTabDouble& rhs) const;

protected:
    double value_;
    MzTabCellStateType state_;
  };

  class OPENMS_DLLAPI MzTabDoubleList
  {
public:
    MzTabDoubleList() = default;

    bool isNull() const;

    void setNull(bool b);

    String toCellString() const;

    void fromCellString(const String& s);

    std::vector<MzTabDouble> get() const;

    void set(const std::vector<MzTabDouble>& entries);

    ~MzTabDoubleList() = default;
protected:
    std::vector<MzTabDouble> entries_;
  };

  class OPENMS_DLLAPI MzTabInteger
  {
public:
    MzTabInteger();

    explicit MzTabInteger(const int v);

    void set(const Int& value);

    Int get() const;

    String toCellString() const;

    void fromCellString(const String& s);

    bool isNull() const;

    void setNull(bool b);

    bool isNaN() const;

    void setNaN();

    bool isInf() const;

    void setInf();

    ~MzTabInteger() = default;
protected:
    Int value_;
    MzTabCellStateType state_;
  };

  class OPENMS_DLLAPI MzTabIntegerList
  {
public:
    MzTabIntegerList() = default;

    bool isNull() const;

    void setNull(bool b);

    String toCellString() const;

    void fromCellString(const String& s);

    std::vector<MzTabInteger> get() const;

    void set(const std::vector<MzTabInteger>& entries);

    ~MzTabIntegerList() = default;
protected:
    std::vector<MzTabInteger> entries_;
  };

  class OPENMS_DLLAPI MzTabBoolean
  {
public:
    MzTabBoolean();

    bool isNull() const;

    void setNull(bool b);

    explicit MzTabBoolean(bool v);

    void set(const bool& value);

    Int get() const;

    String toCellString() const;

    void fromCellString(const String& s);

    ~MzTabBoolean() = default;
protected:
    int value_;
  };

  class OPENMS_DLLAPI MzTabString
  {
public:
    MzTabString();

    explicit MzTabString(const String& s);

    bool isNull() const;

    void setNull(bool b);

    void set(const String& value);

    String get() const;

    String toCellString() const;

    void fromCellString(const String& s);

    ~MzTabString() = default;
protected:
    String value_;
  };

  typedef std::pair<String, MzTabString> MzTabOptionalColumnEntry; //<  column name (not null able), value (null able)

  class OPENMS_DLLAPI MzTabParameter
  {
  public:
    MzTabParameter();

    bool isNull() const;

    void setNull(bool b);

    void setCVLabel(const String& CV_label);

    void setAccession(const String& accession);

    void setName(const String& name);

    void setValue(const String& value);

    String getCVLabel() const;

    String getAccession() const;

    String getName() const;

    String getValue() const;

    String toCellString() const;

    void fromCellString(const String& s);

    ~MzTabParameter() = default;
  protected:
    String CV_label_;
    String accession_;
    String name_;
    String value_;
  };

  class OPENMS_DLLAPI MzTabParameterList
  {
  public:
    MzTabParameterList() = default;

    bool isNull() const;

    void setNull(bool b);

    String toCellString() const;

    void fromCellString(const String& s);

    std::vector<MzTabParameter> get() const;

    void set(const std::vector<MzTabParameter>& parameters);

    ~MzTabParameterList() = default;

  protected:
    std::vector<MzTabParameter> parameters_;
  };

  class OPENMS_DLLAPI MzTabStringList
  {
  public:
    MzTabStringList();

    bool isNull() const;

    void setNull(bool b);

    /// needed for e.g. ambiguity_members and GO accessions as these use ',' as separator while the others use '|'
    void setSeparator(char sep);

    String toCellString() const;

    void fromCellString(const String& s);

    std::vector<MzTabString> get() const;

    void set(const std::vector<MzTabString>& entries);

    ~MzTabStringList() = default;
  protected:
    std::vector<MzTabString> entries_;
    char sep_;
  };

  class OPENMS_DLLAPI MzTabSpectraRef
  {
  public:
    MzTabSpectraRef();

    bool isNull() const;

    void setNull(bool b);

    void setMSFile(Size index);

    void setSpecRef(const String& spec_ref);

    String getSpecRef() const;

    Size getMSFile() const;

    void setSpecRefFile(const String& spec_ref);

    String toCellString() const;

    void fromCellString(const String& s);

    ~MzTabSpectraRef() = default;
  protected:
    Size ms_run_; //< number is specified in the meta data section.
    String spec_ref_;
  };

  // MTD
  struct OPENMS_DLLAPI MzTabSoftwareMetaData
  {
    MzTabParameter software;
    std::map<Size, MzTabString> setting;
  };

  struct OPENMS_DLLAPI MzTabSampleMetaData
  {
    MzTabString description;
    std::map<Size, MzTabParameter> species;
    std::map<Size, MzTabParameter> tissue;
    std::map<Size, MzTabParameter> cell_type;
    std::map<Size, MzTabParameter> disease;
    std::map<Size, MzTabParameter> custom;
  };

  struct OPENMS_DLLAPI MzTabCVMetaData
  {
    MzTabString label;
    MzTabString full_name;
    MzTabString version;
    MzTabString url;
  };

  struct OPENMS_DLLAPI MzTabInstrumentMetaData
  {
    MzTabParameter name;
    MzTabParameter source;
    std::map<Size, MzTabParameter> analyzer;
    MzTabParameter detector;
  };

  struct OPENMS_DLLAPI MzTabContactMetaData
  {
    MzTabString name;
    MzTabString affiliation;
    MzTabString email;
  };

  class OPENMS_DLLAPI MzTabBase
  {
  public:
    MzTabBase() = default;
    virtual ~MzTabBase() = default;

  protected:
    /// Helper function for "get...OptionalColumnNames" functions
    template <typename SectionRows>
    std::vector<String> getOptionalColumnNames_(const SectionRows& rows) const
    {
      // vector is used to preserve the column order
      std::vector<String> names;
      if (!rows.empty())
      {
        for (typename SectionRows::const_iterator it = rows.begin(); it != rows.end(); ++it)
        {
          for (auto it_opt = it->opt_.cbegin(); it_opt != it->opt_.cend(); ++it_opt)
          {
            if (std::find(names.begin(), names.end(), it_opt->first) == names.end())
            {
              names.push_back(it_opt->first);
            }
          }
        }
      }
      return names;
    }
  };
} // namespace OpenMS
