// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
    if (b)
    {
      parameters_.clear();
    }
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

}
