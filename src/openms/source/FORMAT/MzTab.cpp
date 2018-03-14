// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#include <OpenMS/FORMAT/MzTab.h>

#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/METADATA/MetaInfoInterfaceUtils.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/CONCEPT/VersionInfo.h>

using namespace std;

namespace OpenMS
{

  MzTabParameterList::~MzTabParameterList()
  {
  }

  bool MzTabParameterList::isNull() const
  {
    return parameters_.empty();
  }

  void MzTabParameterList::setNull(bool b)
  {
    if (b) { parameters_.clear(); }
  }

  String MzTabParameterList::toCellString() const
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

  void MzTabParameterList::fromCellString(const String& s)
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
          throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("MzTabParameter in MzTabParameterList must not be null '") + s);
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

  MzTabStringList::~MzTabStringList()
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

  String MzTabStringList::toCellString() const
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

  void MzTabStringList::fromCellString(const String& s)
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

  std::vector<MzTabString> MzTabStringList::get() const
  {
    return entries_;
  }

  void MzTabStringList::set(const std::vector<MzTabString>& entries)
  {
    entries_ = entries;
  }

  MzTabModification::MzTabModification()
  {
  }

  MzTabModification::~MzTabModification()
  {
  }

  bool MzTabModification::isNull() const
  {
    return pos_param_pairs_.empty() && mod_identifier_.isNull();
  }

  void MzTabModification::setNull(bool b)
  {
    if (b)
    {
      pos_param_pairs_.clear();
      mod_identifier_.setNull(true);
    }
  }

  void MzTabModification::setPositionsAndParameters(const std::vector<std::pair<Size, MzTabParameter> >& ppp)
  {
    pos_param_pairs_ = ppp;
  }

  std::vector<std::pair<Size, MzTabParameter> > MzTabModification::getPositionsAndParameters() const
  {
    return pos_param_pairs_;
  }

  void MzTabModification::setModificationIdentifier(const MzTabString& mod_id)
  {
    mod_identifier_ = mod_id;
  }

  MzTabString MzTabModification::getModOrSubstIdentifier() const
  {
    assert(!isNull());
    return mod_identifier_;
  }

  String MzTabModification::toCellString() const
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
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Modification or Substitution identifier MUST NOT be null or empty in MzTabModification"));
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

  void MzTabModification::fromCellString(const String& s)
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
          throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Can't convert to MzTabModification from '") + s);
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

  MzTabModificationList::~MzTabModificationList()
  {

  }

  bool MzTabModificationList::isNull() const
  {
    return entries_.empty();
  }

  void MzTabModificationList::setNull(bool b)
  {
    if (b)
    {
      entries_.clear();
    }
  }

  String MzTabModificationList::toCellString() const
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

  void MzTabModificationList::fromCellString(const String& s)
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

  std::vector<MzTabModification> MzTabModificationList::get() const
  {
    return entries_;
  }

  void MzTabModificationList::set(const std::vector<MzTabModification>& entries)
  {
    entries_ = entries;
  }

  MzTabSpectraRef::MzTabSpectraRef() :
    ms_run_(0)
  {
  }

  MzTabSpectraRef::~MzTabSpectraRef()
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

  void MzTabSpectraRef::setSpecRef(String spec_ref)
  {
    assert(!spec_ref.empty());
    if (!spec_ref.empty())
    {
      spec_ref_ = spec_ref;
    }
  }

  String MzTabSpectraRef::getSpecRef() const
  {
    assert(!isNull());
    return spec_ref_;
  }

  Size MzTabSpectraRef::getMSFile() const
  {
    assert(!isNull());
    return ms_run_;
  }

  void MzTabSpectraRef::setSpecRefFile(const String& spec_ref)
  {
    assert(!spec_ref.empty());
    if (!spec_ref.empty())
    {
      spec_ref_ = spec_ref;
    }
  }

  String MzTabSpectraRef::toCellString() const
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

  void MzTabSpectraRef::fromCellString(const String& s)
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
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Can not convert to MzTabSpectraRef from '") + s);
      }

      spec_ref_ = fields[1];
      ms_run_ = (Size)(fields[0].substitute("ms_run[", "").remove(']').toInt());
    }
  }

  MzTabProteinSectionRow::MzTabProteinSectionRow()
  {
    // use "," as list separator because "|" can be used for go terms and protein accessions
    go_terms.setSeparator(',');
    ambiguity_members.setSeparator(',');
  }

  MzTabMetaData::MzTabMetaData()
  {
    mz_tab_version.fromCellString(String("1.0.0"));
  }

  MzTab::MzTab()
  {

  }

  MzTab::~MzTab()
  {

  }

  const MzTabMetaData& MzTab::getMetaData() const
  {
    return meta_data_;
  }

  void MzTab::setMetaData(const MzTabMetaData& md)
  {
    meta_data_ = md;
  }

  const MzTabProteinSectionRows& MzTab::getProteinSectionRows() const
  {
    return protein_data_;
  }

  void MzTab::setProteinSectionRows(const MzTabProteinSectionRows& psd)
  {
    protein_data_ = psd;
  }

  const MzTabPeptideSectionRows& MzTab::getPeptideSectionRows() const
  {
    return peptide_data_;
  }

  void MzTab::setPeptideSectionRows(const MzTabPeptideSectionRows& psd)
  {
    peptide_data_ = psd;
  }

  const MzTabPSMSectionRows& MzTab::getPSMSectionRows() const
  {
    return psm_data_;
  }

  void MzTab::setPSMSectionRows(const MzTabPSMSectionRows& psd)
  {
    psm_data_ = psd;
  }

  void MzTab::setCommentRows(const std::map<Size, String>& com)
  {
    comment_rows_ = com;
  }

  void MzTab::setEmptyRows(const std::vector<Size>& empty)
  {
    empty_rows_ = empty;
  }

  const std::vector<Size>& MzTab::getEmptyRows() const
  {
    return empty_rows_;
  }

  const std::map<Size, String>& MzTab::getCommentRows() const
  {
    return comment_rows_;
  }

  const MzTabSmallMoleculeSectionRows& MzTab::getSmallMoleculeSectionRows() const
  {
    return small_molecule_data_;
  }

  void MzTab::setSmallMoleculeSectionRows(const MzTabSmallMoleculeSectionRows& smsd)
  {
    small_molecule_data_ = smsd;
  }

  std::vector<String> MzTab::getProteinOptionalColumnNames() const
  {
    // vector is used to preserve the column order
    std::vector<String> names;
    if (!protein_data_.empty())
    {
      for (MzTabProteinSectionRows::const_iterator it = protein_data_.begin(); it != protein_data_.end(); ++it)
      {
        for (std::vector<MzTabOptionalColumnEntry>::const_iterator it_opt = it->opt_.begin(); it_opt != it->opt_.end(); ++it_opt)
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

  std::vector<String> MzTab::getPeptideOptionalColumnNames() const
  {
    // vector is used to preserve the column order
    std::vector<String> names;
    if (!peptide_data_.empty())
    {
      for (MzTabPeptideSectionRows::const_iterator it = peptide_data_.begin(); it != peptide_data_.end(); ++it)
      {
        for (std::vector<MzTabOptionalColumnEntry>::const_iterator it_opt = it->opt_.begin(); it_opt != it->opt_.end(); ++it_opt)
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

  std::vector<String> MzTab::getPSMOptionalColumnNames() const
  {
    // vector is used to preserve the column order
    std::vector<String> names;
    if (!psm_data_.empty())
    {
      for (MzTabPSMSectionRows::const_iterator it = psm_data_.begin(); it != psm_data_.end(); ++it)
      {
        for (std::vector<MzTabOptionalColumnEntry>::const_iterator it_opt = it->opt_.begin(); it_opt != it->opt_.end(); ++it_opt)
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

  std::vector<String> MzTab::getSmallMoleculeOptionalColumnNames() const
  {
    // vector is used to preserve the column order
    std::vector<String> names;
    if (!small_molecule_data_.empty())
    {
      for (MzTabSmallMoleculeSectionRows::const_iterator it = small_molecule_data_.begin(); it != small_molecule_data_.end(); ++it)
      {
        for (std::vector<MzTabOptionalColumnEntry>::const_iterator it_opt = it->opt_.begin(); it_opt != it->opt_.end(); ++it_opt)
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

  MzTabParameter::MzTabParameter()
    : CV_label_(""),
    accession_(""),
    name_(""),
    value_("")
  {

  }

  MzTabParameter::~MzTabParameter()
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

  void MzTabParameter::setCVLabel(const String& CV_label)
  {
    CV_label_ = CV_label;
  }

  void MzTabParameter::setAccession(const String& accession)
  {
    accession_ = accession;
  }

  void MzTabParameter::setName(const String& name)
  {
    name_ = name;
  }

  void MzTabParameter::setValue(const String& value)
  {
    value_ = value;
  }

  String MzTabParameter::getCVLabel() const
  {
    assert(!isNull());
    return CV_label_;
  }

  String MzTabParameter::getAccession() const
  {
    assert(!isNull());
    return accession_;
  }

  String MzTabParameter::getName() const
  {
    assert(!isNull());
    return name_;
  }

  String MzTabParameter::getValue() const
  {
    assert(!isNull());
    return value_;
  }

  String MzTabParameter::toCellString() const
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
      }
      else
      {
        ret += name_;
      }

      ret += String(", ");

      if (value_.hasSubstring(", "))
      {
        ret += String("\"") + value_ + String("\""); // quote value if it contains a ","
      }
      else
      {
        ret += value_;
      }

      ret += "]";
      return ret;
    }
  }

  void MzTabParameter::fromCellString(const String& s)
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
      for (String::const_iterator sit = s.begin(); sit != s.end(); ++sit)
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
            fields.push_back(field.trim());
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

      fields.push_back(field.trim());

      if (fields.size() != 4)
      {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Could not convert String '") + s + "' to MzTabParameter");
      }

      CV_label_ = fields[0];
      accession_ = fields[1];
      name_ = fields[2];
      value_ = fields[3];
    }
  }

  MzTabString::MzTabString(const String& s)
  {
    set(s);
  }

  MzTabString::~MzTabString()
  {

  }

  void MzTabString::set(const String& value)
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
      value_.trim();
    }
  }

  String MzTabString::get() const
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

  String MzTabString::toCellString() const
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

  void MzTabString::fromCellString(const String& s)
  {
    set(s);
  }

  MzTabBoolean::MzTabBoolean(bool v)
  {
    set(v);
  }

  MzTabBoolean::~MzTabBoolean()
  {

  }

  MzTabBoolean::MzTabBoolean()
    : value_(false)
  {

  }

  void MzTabBoolean::set(const bool& value)
  {
    setNull(false);
    value_ = value;
  }

  Int MzTabBoolean::get() const
  {
    return value_;
  }

  String MzTabBoolean::toCellString() const
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

  void MzTabBoolean::fromCellString(const String& s)
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
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Could not convert String '") + s + "' to MzTabBoolean");
      }
    }
  }

  MzTabIntegerList::MzTabIntegerList()
  {
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

  String MzTabIntegerList::toCellString() const
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

  void MzTabIntegerList::fromCellString(const String& s)
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

  MzTabInteger::~MzTabInteger()
  {

  }

  MzTabInteger::MzTabInteger()
    : value_(0)
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
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Trying to extract MzTab Integer value from non-integer valued cell. Did you check the cell state before querying the value?"));
    }
  }

  String MzTabInteger::toCellString() const
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

  void MzTabInteger::fromCellString(const String& s)
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
    else // default case
    {
      set(lower.toInt());
    }
  }

  MzTabNullAbleBase::MzTabNullAbleBase() :
    null_(true)
  {
  }

  MzTabNullAbleBase::~MzTabNullAbleBase()
  {
  }

  bool MzTabNullAbleBase::isNull() const
  {
    return null_;
  }

  void MzTabNullAbleBase::setNull(bool b)
  {
    null_ = b;
  }

  MzTabNullNaNAndInfAbleBase::MzTabNullNaNAndInfAbleBase() :
    state_(MZTAB_CELLSTATE_NULL)
  {
  }

  MzTabNullNaNAndInfAbleBase::~MzTabNullNaNAndInfAbleBase()
  {
  }

  bool MzTabNullNaNAndInfAbleBase::isNull() const
  {
    return state_ == MZTAB_CELLSTATE_NULL;
  }

  void MzTabNullNaNAndInfAbleBase::setNull(bool b)
  {
    state_ = b ? MZTAB_CELLSTATE_NULL : MZTAB_CELLSTATE_DEFAULT;
  }

  bool MzTabNullNaNAndInfAbleBase::isNaN() const
  {
    return state_ == MZTAB_CELLSTATE_NAN;
  }

  void MzTabNullNaNAndInfAbleBase::setNaN()
  {
    state_ = MZTAB_CELLSTATE_NAN;
  }

  bool MzTabNullNaNAndInfAbleBase::isInf() const
  {
    return state_ == MZTAB_CELLSTATE_INF;
  }

  void MzTabNullNaNAndInfAbleBase::setInf()
  {
    state_ = MZTAB_CELLSTATE_INF;
  }

  MzTabDouble::MzTabDouble()
    : value_(0.0)
  {
  }

  MzTabDouble::MzTabDouble(const double v)
  {
    set(v);
  }

  MzTabDouble::~MzTabDouble()
  {

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
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Trying to extract MzTab Double value from non-double valued cell. Did you check the cell state before querying the value?"));
    }

    return value_;
  }

  String MzTabDouble::toCellString() const
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

  void MzTabDouble::fromCellString(const String& s)
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
    else // default case
    {
      set(lower.toDouble());
    }
  }

  MzTabDoubleList::MzTabDoubleList()
  {
  }

  MzTabDoubleList::~MzTabDoubleList()
  {

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

  String MzTabDoubleList::toCellString() const
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

  void MzTabDoubleList::fromCellString(const String& s)
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

  std::vector<MzTabDouble> MzTabDoubleList::get() const
  {
    return entries_;
  }

  void MzTabDoubleList::set(const std::vector<MzTabDouble>& entries)
  {
    entries_ = entries;
  }

  MzTabNullAbleInterface::~MzTabNullAbleInterface()
  {
  }

  MzTabNullNaNAndInfAbleInterface::~MzTabNullNaNAndInfAbleInterface()
  {
  }

  void MzTab::addPepEvidenceToRows(const vector<PeptideEvidence>& peptide_evidences, MzTabPSMSectionRow& row, MzTabPSMSectionRows& rows)
  {
    if (!peptide_evidences.empty())
    {
      for (Size i = 0; i != peptide_evidences.size(); ++i)
      {
        // get AABefore and AAAfter as well as start and end for all pep evidences

        // pre/post
        // from spec: Amino acid preceding the peptide (coming from the PSM) in the protein
        // sequence. If unknown “null” MUST be used, if the peptide is N-terminal “-“
        // MUST be used.
        if (peptide_evidences[i].getAABefore() == PeptideEvidence::UNKNOWN_AA)
        {
          row.pre = MzTabString("null");
        }
        else if (peptide_evidences[i].getAABefore() == PeptideEvidence::N_TERMINAL_AA)
        {
          row.pre = MzTabString("-");
        }
        else
        {
          row.pre = MzTabString(String(peptide_evidences[i].getAABefore()));
        }

        if (peptide_evidences[i].getAAAfter() == PeptideEvidence::UNKNOWN_AA)
        {
          row.post = MzTabString("null");
        }
        else if (peptide_evidences[i].getAAAfter() == PeptideEvidence::C_TERMINAL_AA)
        {
          row.post = MzTabString("-");
        }
        else
        {
          row.post = MzTabString(String(peptide_evidences[i].getAAAfter()));
        }

        // start/end
        if (peptide_evidences[i].getStart() == PeptideEvidence::UNKNOWN_POSITION)
        {
          row.start = MzTabString("null");
        }
        else
        {
          row.start = MzTabString(String(peptide_evidences[i].getStart() + 1)); // counting in mzTab starts at 1
        }

        if (peptide_evidences[i].getEnd() == PeptideEvidence::UNKNOWN_POSITION)
        {
          row.end = MzTabString("null");
        }
        else
        {
          row.end = MzTabString(String(peptide_evidences[i].getEnd() + 1)); // counting in mzTab starts at 1
        }

        row.accession = MzTabString(peptide_evidences[i].getProteinAccession());

        rows.push_back(row);
      }
    }
    else
    { // report without pep evidence information
      row.pre = MzTabString("null");
      row.post = MzTabString("null");
      row.start = MzTabString("null");
      row.end = MzTabString("null");
      rows.push_back(row);
    }
  }

  void MzTab::addMetaInfoToOptionalColumns(
    const set<String>& keys, 
    vector<MzTabOptionalColumnEntry>& opt, 
    const String id, 
    const MetaInfoInterface meta)
  {
    for (String const & key : keys)
    {
      MzTabOptionalColumnEntry opt_entry;
      opt_entry.first = String("opt_") + id + String("_") + String(key).substitute(' ','_');
      if (meta.metaValueExists(key))
      {
        opt_entry.second = MzTabString(meta.getMetaValue(key).toString().substitute(' ','_'));
      } // otherwise it is default ("null")
      opt.push_back(opt_entry);
    }
  }

  map<Size, MzTabModificationMetaData> MzTab::generateMzTabStringFromModifications(const vector<String>& mods)
  {
    map<Size, MzTabModificationMetaData> mods_mztab;
    Size index(1);
    for (String const & s : mods)
    {
      MzTabModificationMetaData mod;
      MzTabParameter mp;
      mp.setCVLabel("UNIMOD");
      ModificationsDB* mod_db = ModificationsDB::getInstance();
      // MzTab standard is to just report Unimod accession.
      ResidueModification m = mod_db->getModification(s);
      String unimod_accession = m.getUniModAccession();
      mp.setAccession(unimod_accession.toUpper());
      mp.setName(m.getId());
      mod.modification = mp;

      if (m.getTermSpecificity() == ResidueModification::C_TERM)
      {
        mod.position = MzTabString("Any C-term");
      }
      else if (m.getTermSpecificity() == ResidueModification::N_TERM)
      {
        mod.position = MzTabString("Any N-term");
      }
      else if (m.getTermSpecificity() == ResidueModification::ANYWHERE)
      {
        mod.position = MzTabString("Anywhere");
      }
      else if (m.getTermSpecificity() == ResidueModification::PROTEIN_C_TERM)
      {
        mod.position = MzTabString("Protein C-term");
      }
      else if (m.getTermSpecificity() == ResidueModification::PROTEIN_N_TERM)
      {
        mod.position = MzTabString("Protein N-term");
      }

      mod.site = MzTabString(m.getOrigin());
      mods_mztab[index] = mod;
      ++index;
    }
    return mods_mztab;
  }

  map<Size, MzTabModificationMetaData> MzTab::generateMzTabStringFromVariableModifications(const vector<String>& mods)
  {
    if (mods.empty())
    {
      map<Size, MzTabModificationMetaData> mods_mztab;
      MzTabModificationMetaData mod_mtd;
      mod_mtd.modification.fromCellString("[MS, MS:1002454, No variable modifications searched, ]");
      mods_mztab.insert(make_pair(1, mod_mtd));
      return mods_mztab;
    }
    else
    {
      return generateMzTabStringFromModifications(mods);
    }
  }

  map<Size, MzTabModificationMetaData> MzTab::generateMzTabStringFromFixedModifications(const vector<String>& mods)
  {
    if (mods.empty())
    {
      map<Size, MzTabModificationMetaData> mods_mztab;
      MzTabModificationMetaData mod_mtd;
      mod_mtd.modification.fromCellString("[MS, MS:1002453, No fixed modifications searched, ]");
      mods_mztab.insert(make_pair(1, mod_mtd));
      return mods_mztab;
    }
    else
    {
      return generateMzTabStringFromModifications(mods);
    }
  }

  MzTab MzTab::exportFeatureMapToMzTab(
    const FeatureMap & feature_map, 
    const String & filename)
  {
    LOG_INFO << "exporting feature map: \"" << filename << "\" to mzTab: " << std::endl;
    MzTab mztab;
    MzTabMetaData meta_data;

    vector<ProteinIdentification> prot_ids = feature_map.getProteinIdentifications();        
    vector<String> var_mods, fixed_mods;
    MzTabString db, db_version;
    if (!prot_ids.empty())
    {
      ProteinIdentification::SearchParameters sp = prot_ids[0].getSearchParameters();
      var_mods = sp.variable_modifications;
      fixed_mods = sp.fixed_modifications;
      db = sp.db.empty() ? MzTabString() : MzTabString(sp.db);
      db_version = sp.db_version.empty() ? MzTabString() : MzTabString(sp.db_version);
    }

    meta_data.variable_mod = generateMzTabStringFromVariableModifications(var_mods);
    meta_data.fixed_mod = generateMzTabStringFromFixedModifications(fixed_mods);


    // mandatory meta values
    meta_data.mz_tab_type = MzTabString("Quantification");
    meta_data.mz_tab_mode = MzTabString("Summary");
    meta_data.description = MzTabString("OpenMS export from featureXML");

    MzTabMSRunMetaData ms_run;
    StringList spectra_data;
    feature_map.getPrimaryMSRunPath(spectra_data);
    ms_run.location = spectra_data.empty() ? MzTabString("null") : MzTabString(spectra_data[0]);
    meta_data.ms_run[1] = ms_run;
    meta_data.uri[1] = MzTabString(filename);
    meta_data.psm_search_engine_score[1] = MzTabParameter(); // TODO: we currently only support psm search engine scores annotated to the identification run
    meta_data.peptide_search_engine_score[1] = MzTabParameter();

    mztab.setMetaData(meta_data);

    // pre-analyze data for occuring meta values at feature and peptide hit level
    // these are used to build optional columns containing the meta values in internal data structures

    set<String> feature_user_value_keys;
    set<String> peptide_hit_user_value_keys;
    for (Size i = 0; i < feature_map.size(); ++i)
    {
      const Feature& f = feature_map[i];
      vector<String> keys;
      f.getKeys(keys); //TODO: why not just return it?
      for (String & s : keys) 
      { 
        if (s.has(' '))
         {
           s.substitute(' ', '_');
         }
      }

      feature_user_value_keys.insert(keys.begin(), keys.end());

      const vector<PeptideIdentification>& pep_ids = f.getPeptideIdentifications();
      for (PeptideIdentification const & pep_id : pep_ids)
      {
        for (PeptideHit const & hit : pep_id.getHits())
        {
          vector<String> ph_keys;
          hit.getKeys(ph_keys);
          for (String & s : ph_keys) 
          { 
            if (s.has(' '))
            {
              s.substitute(' ', '_');
            }
          } 
          peptide_hit_user_value_keys.insert(ph_keys.begin(), ph_keys.end());
        }
      }
    }

    MzTabPeptideSectionRows rows;
    for (Size i = 0; i < feature_map.size(); ++i)
    {
      MzTabPeptideSectionRow row;
      const Feature& f = feature_map[i];
      row.mass_to_charge = MzTabDouble(f.getMZ());
      MzTabDoubleList rt_list;
      vector<MzTabDouble> rts;
      rts.push_back(MzTabDouble(f.getRT()));
      rt_list.set(rts);
      row.retention_time = rt_list;

      // set rt window if a bounding box has been set
      vector<MzTabDouble> window;
      if (f.getConvexHull().getBoundingBox() != DBoundingBox<2>())
      {
        window.push_back(MzTabDouble(f.getConvexHull().getBoundingBox().minX()));
        window.push_back(MzTabDouble(f.getConvexHull().getBoundingBox().maxX()));
      }

      MzTabDoubleList rt_window;
      rt_window.set(window);
      row.retention_time_window = rt_window;
      row.charge = MzTabInteger(f.getCharge());
      row.peptide_abundance_stdev_study_variable[1];
      row.peptide_abundance_std_error_study_variable[1];
      row.peptide_abundance_study_variable[1] = MzTabDouble(f.getIntensity());
      row.best_search_engine_score[1] = MzTabDouble();
      row.search_engine_score_ms_run[1][1] = MzTabDouble();

      // create opt_ column for peptide sequence containing modification
      MzTabOptionalColumnEntry opt_global_modified_sequence;
      opt_global_modified_sequence.first = String("opt_global_modified_sequence");
      row.opt_.push_back(opt_global_modified_sequence);

      // create and fill opt_ columns for feature (peptide) user values
      addMetaInfoToOptionalColumns(feature_user_value_keys, row.opt_, String("global"), f);

      const vector<PeptideIdentification>& pep_ids = f.getPeptideIdentifications();
      if (pep_ids.empty())
      {
        rows.push_back(row);
        continue;
      }

      // TODO: here we assume that all have the same score type etc.
      vector<PeptideHit> all_hits;
      for (vector<PeptideIdentification>::const_iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
      {
        all_hits.insert(all_hits.end(), it->getHits().begin(), it->getHits().end());
      }

      if (all_hits.empty())
      {
        rows.push_back(row);
        continue;
      }

      // create new peptide id object to assist in sorting
      PeptideIdentification new_pep_id = pep_ids[0];
      new_pep_id.setHits(all_hits);
      new_pep_id.assignRanks();

      const PeptideHit& best_ph = new_pep_id.getHits()[0];
      const AASequence& aas = best_ph.getSequence();
      row.sequence = MzTabString(aas.toUnmodifiedString());

      row.modifications = extractModificationListFromAASequence(aas, fixed_mods);

      const set<String>& accessions = best_ph.extractProteinAccessionsSet();
      const vector<PeptideEvidence> peptide_evidences = best_ph.getPeptideEvidences();

      row.unique = accessions.size() == 1 ? MzTabBoolean(true) : MzTabBoolean(false);
      // select accession of first peptide_evidence as representative ("leading") accession
      row.accession = peptide_evidences.empty() ? MzTabString("null") : MzTabString(peptide_evidences[0].getProteinAccession());
      row.best_search_engine_score[1] = MzTabDouble(best_ph.getScore());
      row.search_engine_score_ms_run[1][1] = MzTabDouble(best_ph.getScore());

      // find opt_global_modified_sequence in opt_ and set it to the OpenMS amino acid string (easier human readable than unimod accessions)
      for (Size i = 0; i != row.opt_.size(); ++i)
      {
        MzTabOptionalColumnEntry& opt_entry = row.opt_[i];

        if (opt_entry.first == String("opt_global_modified_sequence"))
        {
          opt_entry.second = MzTabString(aas.toString());
        }
      }

      // create and fill opt_ columns for psm (PeptideHit) user values
      addMetaInfoToOptionalColumns(peptide_hit_user_value_keys, row.opt_, String("global"), best_ph);

      rows.push_back(row);
    }
    mztab.setPeptideSectionRows(rows);
    return mztab;
  }

  MzTab MzTab::exportIdentificationsToMzTab(
    const vector<ProteinIdentification>& prot_ids, 
    const vector<PeptideIdentification>& peptide_ids, 
    const String& filename)
  {
    LOG_INFO << "exporting identifications: \"" << filename << "\" to mzTab: " << std::endl;
    vector<PeptideIdentification> pep_ids = peptide_ids;
    MzTab mztab;
    MzTabMetaData meta_data;
    vector<String> var_mods, fixed_mods;
    MzTabString db, db_version;

    // search engine and version -> MS runs index
    map<pair<String, String>, set<Size>> search_engine_to_runs;
    map<Size, vector<pair<String, String>>> run_to_search_engines;

    // helper to map between peptide identifications and MS run
    map<size_t, size_t> map_pep_idx_2_run;

    // used to report quantitative study variables
    Size quant_study_variables(0);

    if (!prot_ids.empty())
    {
      // Check if abundances are annotated to the ind. protein groups
      // if so, we will output the abundances as in a quantification file
      for (auto & p : prot_ids[0].getIndistinguishableProteins())
      {
        if (p.getFloatDataArrays().empty() 
          || p.getFloatDataArrays()[0].getName() != "abundances")
        {
          quant_study_variables = 0;
          break;
        }
        quant_study_variables = p.getFloatDataArrays()[0].size();
      } 

      // TODO: maybe use something different to determine if it is inference data
      bool has_inference_data = prot_ids[0].getSearchEngine() == "Fido" ? true : false; 

      MzTabParameter protein_score_type;
      protein_score_type.fromCellString("[,," + prot_ids[0].getSearchEngine() + ",]"); // TODO: check if we need one for every run (should not be redundant!)
      meta_data.protein_search_engine_score[1] = protein_score_type; // TODO add meta value to ProteinIdentification

      // collect variable and fixed modifications from different runs
      for (auto const & pid : prot_ids)
      {
        const ProteinIdentification::SearchParameters & sp = pid.getSearchParameters();
        var_mods.insert(std::end(var_mods), std::begin(sp.variable_modifications), std::end(sp.variable_modifications));
        fixed_mods.insert(std::end(fixed_mods), std::begin(sp.fixed_modifications), std::end(sp.fixed_modifications));
      }
      // make mods unique
      std::sort(var_mods.begin(), var_mods.end());
      auto v_it = std::unique(var_mods.begin(), var_mods.end()); 
      var_mods.resize(std::distance(var_mods.begin(), v_it));
      std::sort(fixed_mods.begin(), fixed_mods.end());
      auto f_it = std::unique(fixed_mods.begin(), fixed_mods.end()); 
      fixed_mods.resize(std::distance(fixed_mods.begin(), f_it));

      // TODO: check if standard should provide run level info
      const ProteinIdentification::SearchParameters & sp = prot_ids[0].getSearchParameters();
      db = sp.db.empty() ? MzTabString() : MzTabString(prot_ids[0].getSearchParameters().db);
      db_version = sp.db_version.empty() ? MzTabString() : MzTabString(prot_ids[0].getSearchParameters().db_version);

      // MS runs of a peptide identification object is stored in 
      // the protein identification object with the same "identifier".
      // Thus, we build a map from psm_idx->run_index (aka index of PeptideHit -> run index)
      // TODO: move to method
      {
        map<String, size_t> map_id_to_run;

        // first: map run identifier to run index
        size_t run_index(1);
        for (auto it = prot_ids.begin(); it != prot_ids.end(); ++it)
        {
          // First entry might be the inference result without (single) associated ms_run. We skip it.
          if (has_inference_data && it == prot_ids.begin()) { continue; }
          map_id_to_run[it->getIdentifier()] = run_index;
          ++run_index;
        }

        // second: map peptide index to run index
        size_t psm_idx(0);
        for (auto it = peptide_ids.begin(); it != peptide_ids.end(); ++it, ++psm_idx)
        {          
          size_t run_idx = map_id_to_run[it->getIdentifier()];
          map_pep_idx_2_run[psm_idx] = run_idx;
        }
      }

      // Determine search engines used in the different MS runs. TODO: move to method
      {
        size_t run_index(1);
        for (auto it = prot_ids.begin(); it != prot_ids.end(); ++it)
        {
          // First entry might be the inference result without (single) associated ms_run. We skip it.
          if (has_inference_data && it == prot_ids.begin()) { continue; }

          size_t hit_index = std::distance(prot_ids.begin(), it);
          const String & search_engine_name = prot_ids[hit_index].getSearchEngine();
          const String & search_engine_version = prot_ids[hit_index].getSearchEngineVersion();
          search_engine_to_runs[make_pair(search_engine_name, search_engine_version)].insert(run_index);
          run_to_search_engines[run_index].push_back(make_pair(search_engine_name, search_engine_version));
          ++run_index;
        }
      }

      Size psm_search_engine_index(1);
      for (auto const & se : search_engine_to_runs) // loop over (unique) search engine names
      {
        MzTabParameter psm_score_type;
        const pair<String, String>& name_version = se.first;
        psm_score_type.fromCellString("[,," + name_version.first + "," + name_version.second + "]");
        meta_data.psm_search_engine_score[psm_search_engine_index] = psm_score_type;
        meta_data.peptide_search_engine_score[psm_search_engine_index] = psm_score_type; // same score type for peptides
      }

      // TODO: sp.digestion_enzyme
      // TODO: sp.missed_cleavages

      ////////////////////////////////////////////////////////////////
      // generate protein section
      MzTabProteinSectionRows protein_rows;

      Size current_run_index(1);
      for (auto it = prot_ids.begin(); it != prot_ids.end(); ++it, ++current_run_index)
      {
        std::vector<ProteinIdentification::ProteinGroup> protein_groups;
        std::vector<ProteinHit> protein_hits; 

        // We only report quantitative data for indistinguishable groups (which may be composed of single proteins).
        // We skip the more extensive reporting of single proteins and general groups with complex shared peptide relations. 
        if (quant_study_variables == 0)
        {
          protein_hits = it->getHits();        
          protein_groups = it->getProteinGroups();
        }

        const std::vector<ProteinIdentification::ProteinGroup>& indist_groups = it->getIndistinguishableProteins();

        MzTabMSRunMetaData ms_run;
        StringList ms_run_in_data;
        it->getPrimaryMSRunPath(ms_run_in_data);
        ms_run.location = ms_run_in_data.empty() ? MzTabString("null") : MzTabString(ms_run_in_data[0]);
        meta_data.ms_run[current_run_index - 1] = ms_run;

        // TODO: add processing information that this file has been exported from "filename"

        // pre-analyze data for occuring meta values at protein hit level
        // these are used to build optional columns containing the meta values in internal data structures
        set<String> protein_hit_user_value_keys = 
          MetaInfoInterfaceUtils::findCommonMetaKeys<vector<ProteinHit>, set<String> >(protein_hits.begin(), protein_hits.end(), 100.0);

        // column headers may not contain spaces
        for (String s : protein_hit_user_value_keys) 
        { 
          if (s.has(' '))
          {
            protein_hit_user_value_keys.erase(s);
            s.substitute(' ', '_');
            protein_hit_user_value_keys.insert(s);
          }
        } 

        // we do not want descriptions twice
        protein_hit_user_value_keys.erase("Description");

        for (Size i = 0; i != protein_hits.size(); ++i)
        {
          const ProteinHit& hit = protein_hits[i];

          MzTabProteinSectionRow protein_row;

          protein_row.accession = MzTabString(hit.getAccession());
          protein_row.description = MzTabString(hit.getDescription()); 
       // protein_row.taxid = hit.getTaxonomyID(); // TODO maybe add as meta value to protein hit NEWT taxonomy for the species.
       // MzTabString species = hit.getSpecies(); // Human readable name of the species
          protein_row.database = db; // Name of the protein database.
          protein_row.database_version = db_version; // String Version of the protein database.
          protein_row.best_search_engine_score[1] = MzTabDouble(hit.getScore());
       // MzTabParameterList search_engine; // Search engine(s) identifying the protein.
       // std::map<Size, MzTabDouble>  best_search_engine_score; // best_search_engine_score[1-n]
       // std::map<Size, std::map<Size, MzTabDouble> > search_engine_score_ms_run; // search_engine_score[index1]_ms_run[index2]
       // MzTabInteger reliability;
       // std::map<Size, MzTabInteger> num_psms_ms_run;
       // std::map<Size, MzTabInteger> num_peptides_distinct_ms_run;
       // std::map<Size, MzTabInteger> num_peptides_unique_ms_run;
       // MzTabModificationList modifications; // Modifications identified in the protein. TODO: mandatory
       // MzTabString uri; // Location of the protein’s source entry.
       // MzTabStringList go_terms; // List of GO terms for the protein.
          double coverage = hit.getCoverage();
          protein_row.protein_coverage = coverage >= 0 ? MzTabDouble(coverage) : MzTabDouble(); // (0-1) Amount of protein sequence identified.
       // std::vector<MzTabOptionalColumnEntry> opt_; // Optional Columns must start with “opt_”

          // create and fill opt_ columns for protein hit user values
          addMetaInfoToOptionalColumns(protein_hit_user_value_keys, protein_row.opt_, String("global"), hit);

          // optional column for protein groups
          MzTabOptionalColumnEntry opt_column_entry;
          opt_column_entry.first = "opt_global_protein_group_type";
          opt_column_entry.second = MzTabString("single_protein");
          protein_row.opt_.push_back(opt_column_entry);
           
          protein_rows.push_back(protein_row);
        }

        /////////////////////////////////////////////////////////////
        // reporting of general protein groups (not supported for quant data)
        for (Size i = 0; i != protein_groups.size(); ++i)
        {
          const ProteinIdentification::ProteinGroup& group = protein_groups[i];
          MzTabProteinSectionRow protein_row;
          protein_row.database = db; // Name of the protein database.
          protein_row.database_version = db_version; // String Version of the protein database.
            
          MzTabStringList ambiguity_members;
          ambiguity_members.setSeparator(',');
          vector<MzTabString> entries;
          for (Size j = 0; j != group.accessions.size() ; ++j)
          {
            // set accession and description to first element of group
            if (j == 0)
            {
              protein_row.accession = MzTabString(group.accessions[j]);
              // protein_row.description  // TODO: how to set description? information not contained in group
            }
            entries.push_back(MzTabString(group.accessions[j]));
          }
          ambiguity_members.set(entries);
          protein_row.ambiguity_members = ambiguity_members; // Alternative protein identifications.
          protein_row.best_search_engine_score[1] = MzTabDouble(group.probability);

          double coverage = group.coverage; // TODO: create getter / setter for coverage
          if (coverage >= 0) { protein_row.protein_coverage = MzTabDouble(coverage); }

          MzTabOptionalColumnEntry opt_column_entry;
          opt_column_entry.first = "opt_global_protein_group_type";
          opt_column_entry.second = MzTabString("protein_group");
          protein_row.opt_.push_back(opt_column_entry);
          protein_rows.push_back(protein_row);
        }

        /////////////////////////////////////////////////////////////
        // reporting of protein groups composed of indistinguishable proteins
        for (Size i = 0; i != indist_groups.size(); ++i)
        {
          const ProteinIdentification::ProteinGroup& group = indist_groups[i];
          MzTabProteinSectionRow protein_row;
          protein_row.database = db; // Name of the protein database.
          protein_row.database_version = db_version; // String Version of the protein database.

          // column: accession and ambiguity_members            
          MzTabStringList ambiguity_members;
          ambiguity_members.setSeparator(',');
          vector<MzTabString> entries;
          for (Size j = 0; j != group.accessions.size() ; ++j)
          {
            // set accession and description to first element of group
            if (j == 0) { protein_row.accession = MzTabString(group.accessions[j]); }
            entries.push_back(MzTabString(group.accessions[j]));
          }
          ambiguity_members.set(entries);
          protein_row.ambiguity_members = ambiguity_members; // set of indistinguishable proteins
          double coverage = group.coverage;
          if (coverage >= 0) { protein_row.protein_coverage = MzTabDouble(coverage); }

          MzTabOptionalColumnEntry opt_column_entry;
          opt_column_entry.first = "opt_global_protein_group_type";
          opt_column_entry.second = MzTabString("indistinguishable_group");
          protein_row.opt_.push_back(opt_column_entry);
          protein_row.best_search_engine_score[1] = MzTabDouble(group.probability);
          
          // Store quantitative value attached to abundances in study variables
          if (group.getFloatDataArrays().size() == 1 
            && group.getFloatDataArrays()[0].getName() == "abundances")
          {
            const ProteinIdentification::ProteinGroup::FloatDataArray & fa = group.getFloatDataArrays()[0];
            Size i(1);
            for (float f : fa)
            {
              protein_row.protein_abundance_assay[i] = MzTabDouble(f); // assay has same information as SV (without design)
              protein_row.protein_abundance_study_variable[i] = MzTabDouble(f);
              protein_row.protein_abundance_stdev_study_variable[i] = MzTabDouble();
              protein_row.protein_abundance_std_error_study_variable[i] = MzTabDouble();
              ++i;
            }
          } 
          protein_rows.push_back(protein_row);
        }
      }
      mztab.setProteinSectionRows(protein_rows);
    }
    // end protein groups

    ////////////////////////////////////////////////////
    // PSMs

    // mandatory meta values
    if (quant_study_variables == 0)
    {
      meta_data.mz_tab_type = MzTabString("Identification");
    }
    else
    {
      meta_data.mz_tab_type = MzTabString("Quantification");
    }

    meta_data.mz_tab_mode = MzTabString("Summary");
    meta_data.description = MzTabString("OpenMS export from idXML");

    meta_data.variable_mod = generateMzTabStringFromModifications(var_mods);
    meta_data.fixed_mod = generateMzTabStringFromModifications(fixed_mods);

    MzTabSoftwareMetaData sw;
    sw.software.fromCellString("[MS,MS:1000752,TOPP software," + VersionInfo::getVersion() + "]");
    meta_data.software[1] = sw;

    MzTabPSMSectionRows rows;
    Size psm_id(0);
    for (auto it = pep_ids.begin(); it != pep_ids.end(); ++it, ++psm_id)
    {
      // skip empty peptide identification objects
      if (it->getHits().empty()) { continue; }

      // sort by rank
      it->assignRanks();

      MzTabPSMSectionRow row;

      // link to spectrum in MS run
      String spectrum_nativeID = it->getMetaValue("spectrum_reference").toString();
      size_t run_index = map_pep_idx_2_run[psm_id];
      MzTabSpectraRef spec_ref;
      row.spectra_ref.setMSFile(run_index);
      row.spectra_ref.setSpecRef(spectrum_nativeID);

      // only consider best peptide hit for export
      const PeptideHit& best_ph = it->getHits()[0];
      const AASequence& aas = best_ph.getSequence();
      row.sequence = MzTabString(aas.toUnmodifiedString());

      // extract all modifications in the current sequence for reporting. In contrast to peptide and protein section all modifications are reported.
      row.modifications = extractModificationListFromAASequence(aas);

      row.PSM_ID = MzTabInteger(psm_id);
      row.database = db;
      row.database_version = db_version;
      MzTabParameterList search_engines;

      if (run_to_search_engines[run_index].size() != 1)
      {
        throw; // multiple search engines not supported
      }

      pair<String, String> name_version = *run_to_search_engines[run_index].begin();
      search_engines.fromCellString("[,," + name_version.first + "," + name_version.second + "]");
      row.search_engine = search_engines;

      row.search_engine_score[1] = MzTabDouble(best_ph.getScore());

      vector<MzTabDouble> rts_vector;
      rts_vector.push_back(MzTabDouble(it->getRT()));

      MzTabDoubleList rts;
      rts.set(rts_vector);
      row.retention_time = rts;
      row.charge = MzTabInteger(best_ph.getCharge());
      row.exp_mass_to_charge = MzTabDouble(it->getMZ());
      row.calc_mass_to_charge = best_ph.getCharge() != 0 ? MzTabDouble(aas.getMonoWeight(Residue::Full, best_ph.getCharge()) / best_ph.getCharge()) : MzTabDouble();

      // add opt_global_modified_sequence in opt_ and set it to the OpenMS amino acid string (easier human readable than unimod accessions)
      MzTabOptionalColumnEntry opt_entry;
      opt_entry.first = String("opt_global_modified_sequence");
      opt_entry.second = MzTabString(aas.toString());
      row.opt_.push_back(opt_entry);

      // meta data on PSMs
      vector<String> ph_keys;
      best_ph.getKeys(ph_keys);
      for (String & s : ph_keys)      
      { 
        if (s.has(' '))
        {
          s.substitute(' ', '_');
        }
      } 
      set<String> ph_key_set(ph_keys.begin(), ph_keys.end());
      addMetaInfoToOptionalColumns(ph_key_set, row.opt_, String("global"), best_ph);

      // meta data on peptide identifications
      vector<String> pid_keys;
      it->getKeys(pid_keys);
      for (String & s : pid_keys)      
      { 
        if (s.has(' '))
        {
          s.substitute(' ', '_');
        }
      } 
      set<String> pid_key_set(pid_keys.begin(), pid_keys.end());
      addMetaInfoToOptionalColumns(pid_key_set, row.opt_, String("global"), *it);

      // TODO Think about if the uniqueness can be determined by # of peptide evidences
      // b/c this would only differ when evidences come from different DBs
      const set<String>& accessions = best_ph.extractProteinAccessionsSet();
      row.unique = accessions.size() == 1 ? MzTabBoolean(true) : MzTabBoolean(false);

      // create row for every PeptideEvidence entry (mapping to a protein)
      const vector<PeptideEvidence> peptide_evidences = best_ph.getPeptideEvidences();

      // pass common row entries and create rows for all peptide evidences
      addPepEvidenceToRows(peptide_evidences, row, rows);
    }

    mztab.setMetaData(meta_data);
    mztab.setPSMSectionRows(rows);

    return mztab;
  }

  MzTabModificationList MzTab::extractModificationListFromAASequence(const AASequence& aas, const vector<String>& fixed_mods)
  {
    MzTabModificationList mod_list;
    vector<MzTabModification> mods;

    if (aas.isModified())
    {
      if (aas.hasNTerminalModification())
      {
        MzTabModification mod;
        const ResidueModification& res_mod = *(aas.getNTerminalModification());
        if (std::find(fixed_mods.begin(), fixed_mods.end(), res_mod.getId()) == fixed_mods.end())
        {
          String unimod = res_mod.getUniModAccession();
          MzTabString unimod_accession = MzTabString(unimod.toUpper());
          vector<std::pair<Size, MzTabParameter> > pos;
          pos.push_back(make_pair(0, MzTabParameter()));
          mod.setModificationIdentifier(unimod_accession);
          mod.setPositionsAndParameters(pos);
          mods.push_back(mod);
        }
      }

      for (Size ai = 0; ai != aas.size(); ++ai)
      {
        if (aas[ai].isModified())
        {
          MzTabModification mod;
          const ResidueModification& res_mod = *(aas[ai].getModification());
          if (std::find(fixed_mods.begin(), fixed_mods.end(), res_mod.getId()) == fixed_mods.end())
          {
            // MzTab standard is to just report Unimod accession.
            String unimod = res_mod.getUniModAccession();
            MzTabString unimod_accession = MzTabString(unimod.toUpper());
            vector<std::pair<Size, MzTabParameter> > pos;
            pos.push_back(make_pair(ai + 1, MzTabParameter()));
            mod.setPositionsAndParameters(pos);
            mod.setModificationIdentifier(unimod_accession);
            mods.push_back(mod);
          }
        }
      }

      if (aas.hasCTerminalModification())
      {
        MzTabModification mod;
        const ResidueModification& res_mod = *(aas.getCTerminalModification());
        if (std::find(fixed_mods.begin(), fixed_mods.end(), res_mod.getId()) == fixed_mods.end())
        {
          String unimod = res_mod.getUniModAccession();
          MzTabString unimod_accession = MzTabString(unimod.toUpper());
          vector<std::pair<Size, MzTabParameter> > pos;
          pos.push_back(make_pair(aas.size() + 1, MzTabParameter()));
          mod.setPositionsAndParameters(pos);
          mod.setModificationIdentifier(unimod_accession);
          mods.push_back(mod);
        }
      }
    }
    mod_list.set(mods);
    return mod_list;
  }

  MzTab MzTab::exportConsensusMapToMzTab(
    const ConsensusMap& consensus_map, 
    const String& filename, 
    const bool export_unidentified_features,
    const bool export_unassigned_ids)
  {
    LOG_INFO << "exporting consensus map: \"" << filename << "\" to mzTab: " << std::endl;
    vector<ProteinIdentification> prot_ids = consensus_map.getProteinIdentifications();

    // extract mapped IDs (TODO: there should be a helper function)
    vector<PeptideIdentification> pep_ids;
    for (Size i = 0; i < consensus_map.size(); ++i)
    {
      const ConsensusFeature& c = consensus_map[i];
      const vector<PeptideIdentification>& p = c.getPeptideIdentifications();
      pep_ids.insert(std::end(pep_ids), std::begin(p), std::end(p));
    }

    // used to export PSMs of unassigned peptide identifications
    if (export_unassigned_ids)
    {
      const vector<PeptideIdentification>& up = consensus_map.getUnassignedPeptideIdentifications();
      pep_ids.insert(std::end(pep_ids), std::begin(up), std::end(up));
    }

    ///////////////////////////////////////////////////////////////////////
    // Export protein/-group quantifications (stored as meta value in protein IDs)
    MzTab mztab = exportIdentificationsToMzTab(prot_ids, pep_ids, filename);

    // determine number of channels
    // TODO: check if/how this works with fractions and multiplexed experiments
    Size n_study_variables = consensus_map.getFileDescriptions().size();

    // collect variable and fixed modifications from different runs
    vector<String> var_mods, fixed_mods;
    for (auto const & pid : prot_ids)
    {
      const ProteinIdentification::SearchParameters & sp = pid.getSearchParameters();
      var_mods.insert(std::end(var_mods), std::begin(sp.variable_modifications), std::end(sp.variable_modifications));
      fixed_mods.insert(std::end(fixed_mods), std::begin(sp.fixed_modifications), std::end(sp.fixed_modifications));
    }
    // make mods unique
    std::sort(var_mods.begin(), var_mods.end());
    auto v_it = std::unique(var_mods.begin(), var_mods.end()); 
    var_mods.resize(std::distance(var_mods.begin(), v_it));
    std::sort(fixed_mods.begin(), fixed_mods.end());
    auto f_it = std::unique(fixed_mods.begin(), fixed_mods.end()); 
    fixed_mods.resize(std::distance(fixed_mods.begin(), f_it));

    ///////////////////////////////////////////////////////////////////////
    // MetaData section
    
    MzTabMetaData meta_data = mztab.getMetaData();

    // add some mandatory meta values
    meta_data.mz_tab_type = MzTabString("Quantification");
    meta_data.mz_tab_mode = MzTabString("Summary");
    meta_data.description = MzTabString("OpenMS export from consensusXML");

    MzTabParameter quantification_method;
    quantification_method.fromCellString("[MS,MS:1001834,LC-MS label-free quantitation analysis,]");
    meta_data.quantification_method = quantification_method;
    MzTabParameter protein_quantification_unit;
    protein_quantification_unit.fromCellString("[,,Abundance,]"); // TODO: add better term to obo
    meta_data.protein_quantification_unit = protein_quantification_unit;
    MzTabParameter peptide_quantification_unit;
    peptide_quantification_unit.fromCellString("[,,Abundance,]");
    meta_data.peptide_quantification_unit = peptide_quantification_unit;
 
    StringList ms_runs;
    meta_data.ms_run.clear();
    consensus_map.getPrimaryMSRunPath(ms_runs);
    for (Size i = 0; i != ms_runs.size(); ++i)
    {
      MzTabMSRunMetaData ms_run;
      ms_run.location = MzTabString(ms_runs[i]);
      MzTabParameter format_param;
      // format_param.fromCellString("[MS,MS:1000768,Thermo nativeID format]"); //TODO: must be read from mzML and stored as meta information
      // ms_run.format = format_param;
      meta_data.ms_run[i + 1] = ms_run;

      // assay meta data
      MzTabAssayMetaData assay;    
      MzTabParameter quantification_reagent;
      quantification_reagent.fromCellString("[MS,MS:1002038,unlabeled sample,]");  //TODO: add support for other labels
      assay.quantification_reagent = quantification_reagent;
      MzTabString ms_run_ref;
      ms_run_ref.fromCellString(String("ms_run[") + String(i+1)  + "]");
      assay.ms_run_ref = ms_run_ref;
      meta_data.assay[ i + 1] = assay;

      // study variable meta data
      MzTabStudyVariableMetaData sv;
      MzTabString sv_assay_refs;
      sv_assay_refs.fromCellString("assay[" + String(i+1)  + "]");
      sv.assay_refs = sv_assay_refs;

      MzTabString sv_description;
      sv_description.fromCellString("no description given"); // TODO: read from design file
      sv.description = sv_description;      
      meta_data.study_variable[ i + 1] = sv;
    }

    mztab.setMetaData(meta_data);

    // Pre-analyze data for re-occurring meta values at consensus feature and peptide hit level.
    // These are stored in optional columns.
    set<String> consensus_feature_user_value_keys;
    set<String> peptide_hit_user_value_keys;
    for (ConsensusFeature const & c : consensus_map)
    {
      vector<String> keys;
      c.getKeys(keys);
      for (String & s : keys) 
      { 
        if (s.has(' '))
        {
          s.substitute(' ', '_');
        }
      } 
      
      consensus_feature_user_value_keys.insert(keys.begin(), keys.end());

      const vector<PeptideIdentification> & pep_ids = c.getPeptideIdentifications();
      for (auto const & pep_id : pep_ids)
      {
        for (auto const & hit : pep_id.getHits())
        {
          vector<String> ph_keys;
          hit.getKeys(ph_keys);
          for (String & s : ph_keys) 
          { 
            if (s.has(' '))
            {
              s.substitute(' ', '_');
            }
          }
          peptide_hit_user_value_keys.insert(ph_keys.begin(), ph_keys.end());
        }
      }
    }

    MzTabPeptideSectionRows rows;

    for (ConsensusFeature const & c : consensus_map)
    {
      MzTabPeptideSectionRow row;

      // create opt_ column for peptide sequence containing modification
      MzTabOptionalColumnEntry opt_global_modified_sequence;
      opt_global_modified_sequence.first = String("opt_global_modified_sequence");
      row.opt_.push_back(opt_global_modified_sequence);

      // create opt_ columns for consensus feature (peptide) user values
      for (String const & key : consensus_feature_user_value_keys)
      {
        MzTabOptionalColumnEntry opt_entry;
        opt_entry.first = String("opt_global_") + key;
        if (c.metaValueExists(key))
        {
          opt_entry.second = MzTabString(c.getMetaValue(key).toString());
        } // otherwise it is default ("null")
        row.opt_.push_back(opt_entry);
      }

      // create opt_ columns for psm (PeptideHit) user values
      for (String const & key : peptide_hit_user_value_keys)
      {
        MzTabOptionalColumnEntry opt_entry;
        opt_entry.first = String("opt_global_") + key;
        // leave value empty as we have to fill it with the value from the best peptide hit
        row.opt_.push_back(opt_entry);
      }

      row.mass_to_charge = MzTabDouble(c.getMZ());
      MzTabDoubleList rt_list;
      vector<MzTabDouble> rts;
      rts.push_back(MzTabDouble(c.getRT()));
      rt_list.set(rts);
      row.retention_time = rt_list;
      MzTabDoubleList rt_window;
      row.retention_time_window = rt_window;
      row.charge = MzTabInteger(c.getCharge());
      row.best_search_engine_score[1] = MzTabDouble();

      // initialize columns
      for (Size study_variable = 1; study_variable <= n_study_variables; ++study_variable)
      {
        row.peptide_abundance_stdev_study_variable[study_variable] = MzTabDouble();
        row.peptide_abundance_std_error_study_variable[study_variable] = MzTabDouble();
        row.peptide_abundance_study_variable[study_variable] = MzTabDouble();
      }

      for (Size ms_run = 1; ms_run <= ms_runs.size(); ++ms_run)
      {
        row.search_engine_score_ms_run[1][ms_run] = MzTabDouble();
      }

      ConsensusFeature::HandleSetType fs = c.getFeatures();
      for (auto fit = fs.begin(); fit != fs.end(); ++fit)
      {
        Size study_variable = fit->getMapIndex() + 1;
        row.peptide_abundance_stdev_study_variable[study_variable];
        row.peptide_abundance_std_error_study_variable[study_variable];
        row.peptide_abundance_study_variable[study_variable] = MzTabDouble(fit->getIntensity());
      }

      vector<PeptideIdentification> pep_ids = c.getPeptideIdentifications();
      if (!pep_ids.empty())
      {
        if (pep_ids.size() != 1)
        {
          throw OpenMS::Exception::IllegalArgument(
              __FILE__
            , __LINE__
            , __FUNCTION__
            , "Consensus features may contain at most one identification. Run IDConflictResolver first to remove ambiguities!");
        }

        pep_ids[0].assignRanks();
        const PeptideHit& best_ph = pep_ids[0].getHits()[0];
        const AASequence& aas = best_ph.getSequence();
        row.sequence = MzTabString(aas.toUnmodifiedString());

        // annotate variable modifications (no fixed ones)
        row.modifications = extractModificationListFromAASequence(aas, fixed_mods);

        const set<String>& accessions = best_ph.extractProteinAccessionsSet();
        const vector<PeptideEvidence> peptide_evidences = best_ph.getPeptideEvidences();

        row.unique = accessions.size() == 1 ? MzTabBoolean(true) : MzTabBoolean(false);
        // select accession of first peptide_evidence as representative ("leading") accession
        row.accession = peptide_evidences.empty() ? MzTabString("null") : MzTabString(peptide_evidences[0].getProteinAccession());

        row.best_search_engine_score[1] = MzTabDouble(best_ph.getScore());

        // TODO: support run level scores - for now we assume we got the same score from every ms run
        for (Size ms_run = 1; ms_run <= ms_runs.size(); ++ms_run)
        {
          row.search_engine_score_ms_run[1][ms_run] = MzTabDouble(best_ph.getScore());
        }

        // fill opt_ columns

        // find opt_global_modified_sequence in opt_ and set it to the OpenMS amino acid string (easier human readable than unimod accessions)
        for (Size i = 0; i != row.opt_.size(); ++i)
        {
          MzTabOptionalColumnEntry& opt_entry = row.opt_[i];

          if (opt_entry.first == String("opt_global_modified_sequence"))
          {
            opt_entry.second = MzTabString(aas.toString());
          }
        }

        // fill opt_ column of psm
        vector<String> ph_keys;
        best_ph.getKeys(ph_keys);
        for (String & s : ph_keys) 
        { 
          if (s.has(' '))
          {
            s.substitute(' ', '_');
          }
        }

        for (Size k = 0; k != ph_keys.size(); ++k)
        {
          const String& key = ph_keys[k];

          // find matching entry in opt_ (TODO: speed this up)
          for (Size i = 0; i != row.opt_.size(); ++i)
          {
            MzTabOptionalColumnEntry& opt_entry = row.opt_[i];

            if (opt_entry.first == String("opt_global_") + key)
            {
              opt_entry.second = MzTabString(best_ph.getMetaValue(key).toString());
            }
          }
        }
      }

      // skip export of unidentified feature TODO: move logic a bit up to me more efficient
      if (export_unidentified_features == false 
        && row.accession.isNull()) 
      { 
        continue; 
      } 
      rows.push_back(row);
    }

    mztab.setPeptideSectionRows(rows);
    return mztab;
  }
}

