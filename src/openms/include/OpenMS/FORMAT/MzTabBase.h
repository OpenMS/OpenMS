// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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

    std::string toCellString() const;

    void fromCellString(const std::string& s);

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

    std::string toCellString() const;

    void fromCellString(const std::string& s);

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

    std::string toCellString() const;

    void fromCellString(const std::string& s);

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

    std::string toCellString() const;

    void fromCellString(const std::string& s);

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

    std::string toCellString() const;

    void fromCellString(const std::string& s);

    ~MzTabBoolean() = default;
protected:
    int value_;
  };

  class OPENMS_DLLAPI MzTabString
  {
public:
    MzTabString();

    explicit MzTabString(const std::string& s);

    bool isNull() const;

    void setNull(bool b);

    void set(const std::string& value);

    std::string get() const;

    std::string toCellString() const;

    void fromCellString(const std::string& s);

    ~MzTabString() = default;
protected:
    std::string value_;
  };

  typedef std::pair<std::string, MzTabString> MzTabOptionalColumnEntry; //<  column name (not null able), value (null able)

  class OPENMS_DLLAPI MzTabParameter
  {
  public:
    MzTabParameter();

    bool isNull() const;

    void setNull(bool b);

    void setCVLabel(const std::string& CV_label);

    void setAccession(const std::string& accession);

    void setName(const std::string& name);

    void setValue(const std::string& value);

    std::string getCVLabel() const;

    std::string getAccession() const;

    std::string getName() const;

    std::string getValue() const;

    std::string toCellString() const;

    void fromCellString(const std::string& s);

    ~MzTabParameter() = default;
  protected:
    std::string CV_label_;
    std::string accession_;
    std::string name_;
    std::string value_;
  };

  class OPENMS_DLLAPI MzTabParameterList
  {
  public:
    MzTabParameterList() = default;

    bool isNull() const;

    void setNull(bool b);

    std::string toCellString() const;

    void fromCellString(const std::string& s);

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

    std::string toCellString() const;

    void fromCellString(const std::string& s);

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

    void setSpecRef(const std::string& spec_ref);

    std::string getSpecRef() const;

    Size getMSFile() const;

    void setSpecRefFile(const std::string& spec_ref);

    std::string toCellString() const;

    void fromCellString(const std::string& s);

    ~MzTabSpectraRef() = default;
  protected:
    Size ms_run_; //< number is specified in the meta data section.
    std::string spec_ref_;
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
    std::vector<std::string> getOptionalColumnNames_(const SectionRows& rows) const
    {
      // vector is used to preserve the column order
      std::vector<std::string> names;
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
      return names;
    }
  };
} // namespace OpenMS
