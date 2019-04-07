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
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/METADATA/MetaInfoInterfaceUtils.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>

#include <tuple>

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
    if (!spec_ref.empty())
    {
      spec_ref_ = spec_ref;
    }
    else
    {
      LOG_WARN << "Spectrum reference not set." << endl;
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

  const MzTabSmallMoleculeSectionRows& MzTab::getSmallMoleculeSectionRows() const
  {
    return small_molecule_data_;
  }

  void MzTab::setSmallMoleculeSectionRows(const MzTabSmallMoleculeSectionRows& smsd)
  {
    small_molecule_data_ = smsd;
  }

  const MzTabNucleicAcidSectionRows& MzTab::getNucleicAcidSectionRows() const
  {
    return nucleic_acid_data_;
  }

  void MzTab::setNucleicAcidSectionRows(const MzTabNucleicAcidSectionRows& nasd)
  {
    nucleic_acid_data_ = nasd;
  }

  const MzTabOligonucleotideSectionRows& MzTab::getOligonucleotideSectionRows() const
  {
    return oligonucleotide_data_;
  }

  void MzTab::setOligonucleotideSectionRows(const MzTabOligonucleotideSectionRows& onsd)
  {
    oligonucleotide_data_ = onsd;
  }

  const MzTabOSMSectionRows& MzTab::getOSMSectionRows() const
  {
    return osm_data_;
  }

  void MzTab::setOSMSectionRows(const MzTabOSMSectionRows& osd)
  {
    osm_data_ = osd;
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

  std::vector<String> MzTab::getProteinOptionalColumnNames() const
  {
    return getOptionalColumnNames_(protein_data_);
  }

  std::vector<String> MzTab::getPeptideOptionalColumnNames() const
  {
    return getOptionalColumnNames_(peptide_data_);
  }

  std::vector<String> MzTab::getPSMOptionalColumnNames() const
  {
    return getOptionalColumnNames_(psm_data_);
  }

  std::vector<String> MzTab::getSmallMoleculeOptionalColumnNames() const
  {
    return getOptionalColumnNames_(small_molecule_data_);
  }

  std::vector<String> MzTab::getNucleicAcidOptionalColumnNames() const
  {
    return getOptionalColumnNames_(nucleic_acid_data_);
  }

  std::vector<String> MzTab::getOligonucleotideOptionalColumnNames() const
  {
    return getOptionalColumnNames_(oligonucleotide_data_);
  }

  std::vector<String> MzTab::getOSMOptionalColumnNames() const
  {
    return getOptionalColumnNames_(osm_data_);
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
      const ResidueModification* m = mod_db->getModification(s);
      String unimod_accession = m->getUniModAccession();
      mp.setAccession(unimod_accession.toUpper());
      mp.setName(m->getId());
      mod.modification = mp;

      if (m->getTermSpecificity() == ResidueModification::C_TERM)
      {
        mod.position = MzTabString("Any C-term");
      }
      else if (m->getTermSpecificity() == ResidueModification::N_TERM)
      {
        mod.position = MzTabString("Any N-term");
      }
      else if (m->getTermSpecificity() == ResidueModification::ANYWHERE)
      {
        mod.position = MzTabString("Anywhere");
      }
      else if (m->getTermSpecificity() == ResidueModification::PROTEIN_C_TERM)
      {
        mod.position = MzTabString("Protein C-term");
      }
      else if (m->getTermSpecificity() == ResidueModification::PROTEIN_N_TERM)
      {
        mod.position = MzTabString("Protein N-term");
      }

      mod.site = MzTabString(m->getOrigin());
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

    if (!spectra_data.empty())
    {
      // prepend file:// if not there yet
      String m = spectra_data[0];
      if (!m.hasPrefix("file://")) {m = String("file://") + m; }
      ms_run.location = MzTabString(m);
    }
    else
    {
      ms_run.location = MzTabString("null");
    }

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
      rts.emplace_back(MzTabDouble(f.getRT()));
      rt_list.set(rts);
      row.retention_time = rt_list;

      // set rt window if a bounding box has been set
      vector<MzTabDouble> window;
      if (f.getConvexHull().getBoundingBox() != DBoundingBox<2>())
      {
        window.emplace_back(MzTabDouble(f.getConvexHull().getBoundingBox().minX()));
        window.emplace_back(MzTabDouble(f.getConvexHull().getBoundingBox().maxX()));
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
    const String& filename,
    bool first_run_inference_only)
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

    // old/secondary/overwritten search engines and versions.
    // TODO we could potentially make a map too, but our mzTabs currently do not support
    //  associating a PSM with multiple SEs in the metadata section. (we write them as opt_ cols)
    vector<pair<String, String>> secondary_search_engines;
    vector<vector<pair<String, String>>> secondary_search_engines_settings;

    // helper to map between peptide identifications and MS run
    //map<size_t, size_t> map_pep_idx_2_run;
    map<pair<size_t,size_t>,size_t> map_run_fileidx_2_msfileidx;
    map<String, size_t> idrun_2_run_index;

    // used to report quantitative study variables
    Size quant_study_variables(0);

    if (!prot_ids.empty())
    {
      // Check if abundances are annotated to the ind. protein groups
      // if so, we will output the abundances as in a quantification file
      // TODO: we currently assume groups are only in the first run, if at all
      //  if we add a field to an ProtIDRun to specify to which condition it belongs,
      //  a vector of ProtIDRuns can potentially hold multiple groupings with quants
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

      // TODO: use a different identifier to determine if it is inference data (check other places!)
      bool skip_first_run = prot_ids[0].hasInferenceData() && first_run_inference_only;
      if (skip_first_run)
      {
        LOG_DEBUG << "MzTab: Inference data provided. Considering first run only for inference data." << std::endl;
      }

      MzTabParameter protein_score_type;
      protein_score_type.fromCellString("[,," + prot_ids[0].getInferenceEngine() + "," + prot_ids[0].getInferenceEngineVersion() + "]"); // TODO: check if we need one for every run (should not be redundant!)
      meta_data.protein_search_engine_score[1] = protein_score_type; // TODO add meta value to ProteinIdentification

      // maps the filenames in the runs to the MSRun index in the metadata section
      //TODO I think this is a temporary map and could be freed after usage
      map<String, size_t> msfilename_2_msfileindex;


      // collect variable and fixed modifications plus file origins from different runs
      size_t current_ms_run_index(1);
      size_t current_idrun_index(0);
      for (auto const & pid : prot_ids)
      {
        idrun_2_run_index[pid.getIdentifier()] = current_idrun_index;
        const ProteinIdentification::SearchParameters & sp = pid.getSearchParameters();
        var_mods.insert(std::end(var_mods), std::begin(sp.variable_modifications), std::end(sp.variable_modifications));
        fixed_mods.insert(std::end(fixed_mods), std::begin(sp.fixed_modifications), std::end(sp.fixed_modifications));


        StringList ms_run_in_data;
        pid.getPrimaryMSRunPath(ms_run_in_data);

        if (!ms_run_in_data.empty())
        {
          // prepend file:// if not there yet
          for (const String& s : ms_run_in_data)
          {

            String m;
            if (!s.hasPrefix("file://"))
            {
               m = String("file://") + s;
            }
            else
            {
               m = s;
            }
            // use the string without file: prefix for the map
            const auto& msfileidxpair_success = msfilename_2_msfileindex.emplace(s, current_ms_run_index);
            if (msfileidxpair_success.second) // newly inserted
            {
              MzTabMSRunMetaData ms_run;
              ms_run.location = MzTabString(m); // use the string with file: prefix
              meta_data.ms_run[current_ms_run_index] = ms_run;
              current_ms_run_index++;
            }
          }
        }
        else
        {
          MzTabMSRunMetaData ms_run;
          ms_run.location = MzTabString("null");
          // next line is a hack. In case we would ever have some idXML where some runs are annotated
          // and others are not. If a run is not annotated use its index as a String key.
          msfilename_2_msfileindex.emplace(String(current_idrun_index), current_ms_run_index);
          meta_data.ms_run[current_ms_run_index] = ms_run;
          current_ms_run_index++;
        }
        current_idrun_index++;
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
      //The following "rescoring" algorithms will overwrite search engine names and scores but not the settings.
      // Settings and old searchengines should then be stored in the ProteinIDRun in Metavalues (see PercolatorAdapter)
      // They are parsed later as "secondary search engines".
      if (prot_ids[0].getSearchEngine() != "Percolator" && prot_ids[0].getSearchEngine() != "IDPosteriorErrorProbability")
      {
        MzTabSoftwareMetaData sesoftwaremd;
        MzTabParameter sesoftware;
        sesoftware.fromCellString("[,," + prot_ids[0].getSearchEngine() + "," + prot_ids[0].getSearchEngineVersion() + "]");
        sesoftwaremd.software = sesoftware;
        sesoftwaremd.setting[1] = sp.db.empty() ? MzTabString() : MzTabString("db:"+sp.db);
        sesoftwaremd.setting[2] = sp.db_version.empty() ? MzTabString() : MzTabString("db_version:"+sp.db_version);
        sesoftwaremd.setting[3] = sp.taxonomy.empty() ? MzTabString() : MzTabString("taxonomy:"+sp.taxonomy);
        sesoftwaremd.setting[4] = MzTabString("fragment_mass_tolerance:"+String(sp.fragment_mass_tolerance));
        sesoftwaremd.setting[5] = MzTabString("fragment_mass_tolerance_unit:" + String(sp.fragment_mass_tolerance_ppm ? "ppm" : "Da"));
        sesoftwaremd.setting[6] = MzTabString("precursor_mass_tolerance:"+String(sp.precursor_mass_tolerance));
        sesoftwaremd.setting[7] = MzTabString("precursor_mass_tolerance_unit:" + String(sp.precursor_mass_tolerance_ppm ? "ppm" : "Da"));
        sesoftwaremd.setting[8] = MzTabString(String("enzyme:") + sp.digestion_enzyme.getName());
        meta_data.software[1] = sesoftwaremd;
      }

      // MS runs of a peptide identification object is stored in
      // the protein identification object with the same "identifier".
      // Thus, we build a map from psm_idx->run_index (aka index of PeptideHit -> run index)
      // TODO: move to method
      {

        //TODO the following assumes that every file occurs max. once in all runs
        size_t run_index(0);
        for (const auto& run : prot_ids)
        {
          // First entry might be the inference result without (single) associated ms_run. We skip it.
          if (skip_first_run && run_index == 1) { continue; }

          StringList files;
          run.getPrimaryMSRunPath(files);
          if (!files.empty())
          {
            size_t file_index(0);
            for (const String& file : files)
            {
              map_run_fileidx_2_msfileidx[{run_index,file_index}] = msfilename_2_msfileindex[file];
              file_index++;
            }
          }
          else
          {
            map_run_fileidx_2_msfileidx[{run_index,0}] = msfilename_2_msfileindex[String(run_index)];
          }
          run_index++;
        }
      }

      // Determine search engines used in the different MS runs. TODO: move to method
      {
        size_t run_index(0);
        for (auto it = prot_ids.begin(); it != prot_ids.end(); ++it)
        {
          // First entry might be the inference result without (single) associated ms_run. We skip it.
          if (skip_first_run && it == prot_ids.begin())
          {
            continue;
          }

          size_t hit_index = std::distance(prot_ids.begin(), it);
          const String &search_engine_name = prot_ids[hit_index].getSearchEngine();
          const String &search_engine_version = prot_ids[hit_index].getSearchEngineVersion();
          search_engine_to_runs[make_pair(search_engine_name, search_engine_version)].insert(run_index);
          run_to_search_engines[run_index].push_back(make_pair(search_engine_name, search_engine_version));
          vector<String> mvkeys;
          const ProteinIdentification::SearchParameters& sp = prot_ids[hit_index].getSearchParameters();
          sp.getKeys(mvkeys);
          if (prot_ids[hit_index].metaValueExists("percolator"))
          {
            secondary_search_engines.emplace_back(make_pair("Percolator", sp.getMetaValue("percolator")));
          }
          for (const String & mvkey : mvkeys)
          {
            // currently only supported for SEs overwritten by PercolatorAdapter
            // we also do not support different settings but assume that they were run with the same settings
            // and that these settings were taken over into the ones for PercolatorAdapter
            if (mvkey.hasPrefix("SE:"))
            {
              secondary_search_engines.emplace_back(make_pair(mvkey.substr(3), sp.getMetaValue(mvkey)));
            }
          }
          secondary_search_engines_settings.resize(secondary_search_engines.size());
          for (const String & mvkey : mvkeys)
          {
            for (size_t s = 0; s < secondary_search_engines.size(); ++s)
            {
              if (mvkey.hasPrefix(secondary_search_engines[s].first))
              {
                secondary_search_engines_settings[s].emplace_back(make_pair(mvkey.substr(secondary_search_engines[s].first.size()+1), sp.getMetaValue(mvkey)));
              }
            }
          }
          ++run_index;
        }
      }
      Size sw_idx(meta_data.software.size());
      Size cnt(0);
      for (auto const & se : secondary_search_engines)
      {
        MzTabSoftwareMetaData sesoftwaremd;
        MzTabParameter sesoftware;
        sesoftware.fromCellString("[,," + se.first + "," + se.second + "]");
        sesoftwaremd.software = sesoftware;
        Size cnt2(1);
        for (auto const & sesetting : secondary_search_engines_settings[cnt])
        {
          sesoftwaremd.setting[cnt2] = MzTabString(sesetting.first + ":" + sesetting.second);
          cnt2++;
        }
        meta_data.software[sw_idx] = sesoftwaremd;
        sw_idx++;
        cnt++;
      }

      //TODO actually we would need to go through all MetaValues in the PeptideHits to get all scores.
      // and also: Why do we go through all runs here, while above for the search engine settings we only take the first.
      // nonsense.
      Size psm_search_engine_index(1);
      for (auto const & se : search_engine_to_runs) // loop over (unique) search engine names
      {
        MzTabParameter psm_score_type;
        const pair<String, String>& name_version = se.first;
        psm_score_type.fromCellString("[,," + name_version.first + "," + name_version.second + "]");
        meta_data.psm_search_engine_score[psm_search_engine_index] = psm_score_type;
        meta_data.peptide_search_engine_score[psm_search_engine_index] = psm_score_type; // same score type for peptides
        psm_search_engine_index++;
      }

      // map (indist.)protein groups to their protein hits (by index).
      map<Size, set<Size>> ind2prot; // indistinguishable protein groups
      map<Size, set<Size>> pg2prot; // general protein groups

      const std::vector<ProteinHit> proteins = prot_ids.front().getHits();

      // map indistinguishable groups to the contained proteins
      const std::vector<ProteinIdentification::ProteinGroup>& indist_groups = prot_ids.front().getIndistinguishableProteins();
      Size ind_idx{0};
      for (const ProteinIdentification::ProteinGroup & p : indist_groups)
      {
        for (const String & a : p.accessions)
        {
          // find protein corresponding to accession stored in group
          auto it = std::find_if(proteins.begin(), proteins.end(), [&a](const ProteinHit & ph)
            {
              return ph.getAccession() == a;
            }
          );
          if (it == proteins.end()) { continue; };
          Size protein_index = std::distance(proteins.begin(), it);
          ind2prot[ind_idx].insert(protein_index);
        }
        ++ind_idx;
      }

      // map general protein groups to the contained proteins
      const std::vector<ProteinIdentification::ProteinGroup>& protein_groups = prot_ids.front().getProteinGroups();
      Size pg_idx{0};
      for (const ProteinIdentification::ProteinGroup & p : protein_groups)
      {
        for (const String & a : p.accessions)
        {
          // find protein corresponding to accession stored in group
          auto it = std::find_if(proteins.begin(), proteins.end(), [&a](const ProteinHit & ph)
            {
              return ph.getAccession() == a;
            }
          );
          if (it == proteins.end()) { continue; }
          Size protein_index = std::distance(proteins.begin(), it);
          pg2prot[pg_idx].insert(protein_index);
        }
      }
      ++pg_idx;
        
      ////////////////////////////////////////////////////////////////
      // generate protein section

      MzTabProteinSectionRows protein_rows;

      for (auto it = prot_ids.begin(); it != prot_ids.end(); ++it)
      {
        const std::vector<ProteinHit>& protein_hits = it->getHits(); 
        const std::vector<ProteinIdentification::ProteinGroup>& indist_groups = it->getIndistinguishableProteins();

        // TODO: add processing information that this file has been exported from "filename"

        // pre-analyze data for occurring meta values at protein hit level
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


        // We only report quantitative data for indistinguishable groups (which may be composed of single proteins).
        // We skip the more extensive reporting of general groups with complex shared peptide relations. 
        std::vector<ProteinIdentification::ProteinGroup> protein_groups;
        if (quant_study_variables == 0)
        {
          protein_groups = it->getProteinGroups();
        }

        /*
        * protein_hits are supposed to contain all inferred proteins (single proteins and part of groups)
        * indist_groups define the indistinguishable groups and reference proteins in protein_hits
        * protein_groups define general protein groups and reference proteins in protein_hits
        */
       if (!skip_first_run)
       {
              
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
          MzTabModificationList modifications; // Modifications identified in the protein.
          const std::set<pair<Size, ResidueModification>>& leader_mods = hit.getModifications();
          for (auto const & m : leader_mods)
          {
            MzTabModification mztab_mod;
            String unimod = m.second.getUniModAccession();
            MzTabString unimod_accession = MzTabString(unimod.toUpper());
            mztab_mod.setModificationIdentifier(unimod_accession);
            vector<std::pair<Size, MzTabParameter> > pos;
            pos.emplace_back(make_pair(m.first, MzTabParameter())); // position, parameter pair (e.g. FLR)
            mztab_mod.setPositionsAndParameters(pos);
          }
          protein_row.modifications = modifications;

       // MzTabString uri; // Location of the protein’s source entry.
       // MzTabStringList go_terms; // List of GO terms for the protein.
          double coverage = hit.getCoverage() / 100.0; // convert percent to fraction
          protein_row.coverage = coverage >= 0 ? MzTabDouble(coverage) : MzTabDouble(); // (0-1) Amount of protein sequence identified.
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
            entries.emplace_back(MzTabString(group.accessions[j]));
          }
          ambiguity_members.set(entries);
          protein_row.ambiguity_members = ambiguity_members; // Alternative protein identifications.
          protein_row.best_search_engine_score[1] = MzTabDouble(group.probability);
      
          protein_row.coverage = MzTabDouble(); 

          MzTabOptionalColumnEntry opt_column_entry;
          opt_column_entry.first = "opt_global_protein_group_type";
          opt_column_entry.second = MzTabString("protein_group");
          protein_row.opt_.push_back(opt_column_entry);
          protein_rows.push_back(protein_row);
        }
       }

        /////////////////////////////////////////////////////////////
        // reporting of protein groups composed of indistinguishable proteins
        for (Size i = 0; i != indist_groups.size(); ++i)
        {
          const ProteinIdentification::ProteinGroup& group = indist_groups[i];

          // get references (indices) into proteins vector
          const set<Size> & protein_hits_idx = ind2prot[i];

          // determine group leader
          const ProteinHit& leader_protein = protein_hits[*protein_hits_idx.begin()];

          MzTabProteinSectionRow protein_row;
          protein_row.database = db; // Name of the protein database.
          protein_row.database_version = db_version; // String Version of the protein database.

          // column: accession and ambiguity_members            
          MzTabStringList ambiguity_members;
          ambiguity_members.setSeparator(',');
          vector<MzTabString> entries;

          // set accession and description to first element of group
          protein_row.accession = MzTabString(leader_protein.getAccession());
          
          // TODO: check with standard if it is important to also place leader at first position 
          //       (because order in set and vector may differ)
          for (Size j = 0; j != group.accessions.size() ; ++j)
          {
            entries.emplace_back(MzTabString(group.accessions[j]));
          }
          ambiguity_members.set(entries);
          protein_row.ambiguity_members = ambiguity_members; // set of indistinguishable proteins

          // annotate if group contains only one or multiple proteins
          MzTabOptionalColumnEntry opt_column_entry;
          opt_column_entry.first = "opt_global_protein_group_type";

          // TODO: we could count the number of targets or set it to target if at least one target is inside the group
          // we will always call them "indistinguishable_proteins" to differentiate between e.g.
          // protein scores based on grouping or on single proteins
          opt_column_entry.second = MzTabString("indistinguishable_proteins");
          protein_row.opt_.push_back(opt_column_entry);

          // column: coverage
          // calculate mean coverage from individual protein coverages
          double coverage{0};
          for (const Size & prot_idx : protein_hits_idx)
          {
            coverage += (1.0 / (double)protein_hits_idx.size()) * 0.01 * protein_hits[prot_idx].getCoverage();
          }
          if (coverage >= 0) { protein_row.coverage = MzTabDouble(coverage); }
                    
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

          // add protein description of first (leader) protein
          protein_row.description = MzTabString(leader_protein.getDescription());
          protein_row.taxid = (leader_protein.metaValueExists("TaxID")) ? 
            MzTabInteger(static_cast<int>(leader_protein.getMetaValue("TaxID"))) :
            MzTabInteger();

          protein_row.species = (leader_protein.metaValueExists("Species")) ? 
            MzTabString(leader_protein.getMetaValue("Species")) :
            MzTabString();

          protein_row.uri = (leader_protein.metaValueExists("URI")) ? 
            MzTabString(leader_protein.getMetaValue("URI")) :
            MzTabString();

          if (leader_protein.metaValueExists("GO"))
          { 
            StringList sl = leader_protein.getMetaValue("GO");
            String s;
            s.concatenate(sl.begin(), sl.end(), ",");
            protein_row.go_terms.fromCellString(s);
          }

          protein_row.best_search_engine_score[1] = MzTabDouble(group.probability); // TODO: group probability or search engine score?

          protein_row.reliability = MzTabInteger();

          MzTabParameterList search_engine; // Search engine(s) identifying the protein.
          protein_row.search_engine = search_engine;

          MzTabModificationList modifications; // Modifications identified in the protein.
          const std::set<pair<Size, ResidueModification>>& leader_mods = leader_protein.getModifications();
          for (auto const & m : leader_mods)
          {
            MzTabModification mztab_mod;
            String unimod = m.second.getUniModAccession();
            MzTabString unimod_accession = MzTabString(unimod.toUpper());
            mztab_mod.setModificationIdentifier(unimod_accession);
            vector<std::pair<Size, MzTabParameter> > pos;

            // mzTab position is one-based, internal is 0-based so we need to +1
            pos.emplace_back(make_pair(m.first + 1, MzTabParameter())); // position, parameter pair (e.g. FLR)
            mztab_mod.setPositionsAndParameters(pos);
            vector<MzTabModification> mztab_mods(1, mztab_mod);
            modifications.set(mztab_mods);
          }
          protein_row.modifications = modifications;

          if (leader_protein.metaValueExists("num_psms_ms_run"))
          {
            const IntList& il = leader_protein.getMetaValue("num_psms_ms_run");
            for (Size i = 0; i != il.size(); ++i)
            {
              protein_row.num_psms_ms_run[i+1] = MzTabInteger(il[i]);
            }
          }

          if (leader_protein.metaValueExists("num_peptides_distinct_ms_run"))
          {
            const IntList& il = leader_protein.getMetaValue("num_peptides_distinct_ms_run");
            for (Size i = 0; i != il.size(); ++i)
            {
              protein_row.num_peptides_distinct_ms_run[i+1] = MzTabInteger(il[i]);
            }
          }

          if (leader_protein.metaValueExists("num_peptides_unique_ms_run"))
          {
            const IntList& il = leader_protein.getMetaValue("num_peptides_unique_ms_run");
            for (Size i = 0; i != il.size(); ++i)
            {
              protein_row.num_peptides_unique_ms_run[i+1] = MzTabInteger(il[i]);
            }
          }

/*
TODO:
Not sure how to handle these:
       // std::map<Size, MzTabDouble>  best_search_engine_score; // best_search_engine_score[1-n]
       // std::map<Size, std::map<Size, MzTabDouble> > search_engine_score_ms_run; // search_engine_score[index1]_ms_run[index2]
*/

          // Add protein(group) row to MzTab 
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
    meta_data.software[std::max<size_t>(1u, meta_data.software.size())] = sw;

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
      size_t run_index = idrun_2_run_index[it->getIdentifier()];
      StringList filenames;
      prot_ids[run_index].getPrimaryMSRunPath(filenames);
      size_t msfile_index(0);
      if (filenames.size() <= 1) //either none or only one file for this run
      {
        msfile_index = map_run_fileidx_2_msfileidx[{run_index, 0}];
      }
      else
      {
        if (it->metaValueExists("map_index"))
        {
          msfile_index = map_run_fileidx_2_msfileidx[{run_index, it->getMetaValue("map_index")}];
        }
        else
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Multiple files in a run, but no map_index in PeptideIdentification found.");
        }
      }

      MzTabSpectraRef spec_ref;
      row.spectra_ref.setMSFile(msfile_index);
      if (spectrum_nativeID.empty())
      {
        LOG_WARN << "spectrum_reference not set in ID with precursor (RT, m/z) " << it->getRT() << ", " << it->getMZ() << endl;
      }
      else
      {
        row.spectra_ref.setSpecRef(spectrum_nativeID);
      }
      
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
        throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, run_to_search_engines[run_index].size()); // multiple search engines not supported yet
      }

      pair<String, String> name_version = *run_to_search_engines[run_index].begin();
      search_engines.fromCellString("[,," + name_version.first + "," + name_version.second + "]");
      row.search_engine = search_engines;

      row.search_engine_score[1] = MzTabDouble(best_ph.getScore());

      vector<MzTabDouble> rts_vector;
      rts_vector.emplace_back(MzTabDouble(it->getRT()));

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
      const vector<PeptideEvidence>& peptide_evidences = best_ph.getPeptideEvidences();

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


/*
 -
 -      // mandatory meta values
 -      meta_data.mz_tab_type = MzTabString("Quantification");
 -      meta_data.mz_tab_mode = MzTabString("Summary");
 -      meta_data.description = MzTabString("Export from consensusXML");
 -
 -      // For consensusXML we export a "Summary Quantification" file. This means we don't need to report feature quantification values at the assay level
 -      // but only at the study variable variable level.
 -
 -      meta_data.variable_mod = generateMzTabStringFromModifications(var_mods);
 -      meta_data.fixed_mod = generateMzTabStringFromModifications(fixed_mods);
 -      meta_data.peptide_search_engine_score[1] = MzTabParameter();
 -      meta_data.psm_search_engine_score[1] = MzTabParameter(); // TODO insert search engine information
 -
 -      StringList ms_runs;
 -      consensus_map.getPrimaryMSRunPath(ms_runs); 
 -
 -      // condense consecutive unique MS runs to get the different MS files
 -      auto it = std::unique(ms_runs.begin(), ms_runs.end());
 -      ms_runs.resize(std::distance(ms_runs.begin(), it)); 
 -
 -      // set run meta data
 -      Size run_index{1};
 -      for (auto const & m : ms_runs)
 -      {
 -        MzTabMSRunMetaData mztab_run_metadata;
 -        mztab_run_metadata.format.fromCellString("[MS,MS:1000584,mzML file,]");
 -        mztab_run_metadata.id_format.fromCellString("[MS,MS:1001530,mzML unique identifier,]");
 -        mztab_run_metadata.location = MzTabString(m);
 -        meta_data.ms_run[run_index] = mztab_run_metadata;
 -        LOG_DEBUG << "Adding MS run for file: " << m << endl;
 -        ++run_index;
 -      }
 -
 -      mztab.setMetaData(meta_data);
 -
 -      // pre-analyze data for occurring meta values at consensus feature and peptide hit level
 -      // these are used to build optional columns containing the meta values in internal data structures
 -      set<String> consensus_feature_user_value_keys;
 */


  MzTab MzTab::exportConsensusMapToMzTab(
    const ConsensusMap& consensus_map, 
    const String& filename, 
    const bool export_unidentified_features,
    const bool export_unassigned_ids,
    String title)
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
    // In this case, the first run is only for inference, get peptide info from the rest of the runs.
    MzTab mztab = exportIdentificationsToMzTab(prot_ids, pep_ids, filename, true);


    // determine number of samples
    ExperimentalDesign ed = ExperimentalDesign::fromConsensusMap(consensus_map);  
    Size n_study_variables = ed.getNumberOfSamples();

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
    meta_data.title = MzTabString(title);

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

    // condense consecutive unique MS runs to get the different MS files
    auto it = std::unique(ms_runs.begin(), ms_runs.end());
    ms_runs.resize(std::distance(ms_runs.begin(), it)); 

    // set run meta data
    Size run_index{1};
    for (String m : ms_runs)
    {
      MzTabMSRunMetaData mztab_run_metadata;
      mztab_run_metadata.format.fromCellString("[MS,MS:1000584,mzML file,]");
      mztab_run_metadata.id_format.fromCellString("[MS,MS:1001530,mzML unique identifier,]");

      // prepend file:// if not there yet
      if (!m.hasPrefix("file://")) {m = String("file://") + m; }

      mztab_run_metadata.location = MzTabString(m);
      meta_data.ms_run[run_index] = mztab_run_metadata;
      LOG_DEBUG << "Adding MS run for file: " << m << endl;
      ++run_index;
    }

    // assay index (and sample index) must be unique numbers 1..n
    // fraction_group + label define the quant. values of an assay
    auto pl2fg = ed.getPathLabelToFractionGroupMapping(false);

    const String & experiment_type = consensus_map.getExperimentType();    

    // assay meta data
    for (auto const & c : consensus_map.getColumnHeaders())
    {
      Size assay_index{1};
      
      MzTabAssayMetaData assay;    
      MzTabParameter quantification_reagent;
      if (experiment_type == "label-free")
      {
        quantification_reagent.fromCellString("[MS,MS:1002038,unlabeled sample,]");
        auto pl = make_pair(c.second.filename, 1); // TODO: only label-free here -> adapt to multiplexed
        assay_index = pl2fg[pl];
      } 
      else if (experiment_type == "labeled_MS1")
      {
        // TODO: check if there are appropriate CV terms
        quantification_reagent.fromCellString("[MS,MS:XXXXXX,MS1 labeled sample," + c.second.label + "]");
        Size label{1};
        if (c.second.metaValueExists("channel_id"))
        {
          label = static_cast<unsigned int>(c.second.getMetaValue("channel_id")) + 1;
        }
        assay_index = label;
      } 
      else if (experiment_type == "labeled_MS2")
      {
        // TODO: check if there are appropriate CV terms
        quantification_reagent.fromCellString("[MS,MS:XXXXXX,MS2 labeled sample," + c.second.label + "]");
        Size label{1};
        if (c.second.metaValueExists("channel_id"))
        {
          label = static_cast<unsigned int>(c.second.getMetaValue("channel_id")) + 1;
        }
        assay_index = label;
      }

      // look up run index by filename
      auto it = find_if(meta_data.ms_run.begin(), meta_data.ms_run.end(), 
        [&c] (const pair<Size, MzTabMSRunMetaData>& m) { return m.second.location.toCellString().hasSuffix(c.second.filename); } );
      Size run_index = it->first;

      meta_data.assay[assay_index].quantification_reagent = quantification_reagent;
      meta_data.assay[assay_index].ms_run_ref.push_back(run_index);
      
      // study variable meta data
      MzTabString sv_description;
      meta_data.study_variable[assay_index].description.fromCellString("no description given");
      IntList al;
      al.push_back(assay_index);
      meta_data.study_variable[assay_index].assay_refs = al;
    }

    mztab.setMetaData(meta_data);

    // optional meta value columns
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
      LOG_DEBUG << "Initializing study variables:" << n_study_variables << endl;
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
        Size study_variable{1};
        const int index = fit->getMapIndex();
        const ConsensusMap::ColumnHeader& ch = consensus_map.getColumnHeaders().at(index);

        if (experiment_type == "label-free")
        {
          // convert from column index to study variable index
          auto pl = make_pair(ch.filename, 1);
          study_variable = pl2fg[pl];
        } 
        else if (experiment_type == "labeled_MS1")
        {
          Size label{1};
          if (ch.metaValueExists("channel_id"))
          {
            label = static_cast<unsigned int>(ch.getMetaValue("channel_id")) + 1;
          }
          study_variable = label;
        } 
        else if (experiment_type == "labeled_MS2")
        {
          Size label{1};
          if (ch.metaValueExists("channel_id"))
          {
            label = static_cast<unsigned int>(ch.getMetaValue("channel_id")) + 1;
          }
          study_variable = label;
        }

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

