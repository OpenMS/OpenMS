// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg  $
// $Authors: Marc Sturm, Andreas Bertsch, Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ControlledVocabulary.h>

#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/SYSTEM/File.h>

#include <iostream>
#include <fstream>
#include <map>

using namespace std;

namespace OpenMS
{

  ControlledVocabulary::CVTerm::CVTerm() :
    name(),
    id(),
    parents(),
    children(),
    obsolete(false),
    description(),
    synonyms(),
    unparsed(),
    xref_type(XRefType::NONE),
    xref_binary()
  {
  }

  ControlledVocabulary::CVTerm::CVTerm(const CVTerm& rhs) = default;

  ControlledVocabulary::CVTerm& ControlledVocabulary::CVTerm::operator=(const CVTerm& rhs)
  {
    if (this != &rhs)
    {
      name = rhs.name;
      id = rhs.id;
      parents = rhs.parents;
      children = rhs.children;
      obsolete = rhs.obsolete;
      description = rhs.description;
      synonyms = rhs.synonyms;
      unparsed = rhs.unparsed;
      xref_type = rhs.xref_type;
      xref_binary = rhs.xref_binary;
      units = rhs.units;
    }
    return *this;
  }

  std::string ControlledVocabulary::CVTerm::getXRefTypeName(XRefType type)
  {
    switch (type)
    {
    case XRefType::XSD_STRING: return "xsd:string";

    case XRefType::XSD_INTEGER: return "xsd:integer";

    case XRefType::XSD_DECIMAL: return "xsd:decimal";

    case XRefType::XSD_NEGATIVE_INTEGER: return "xsd:negativeInteger";

    case XRefType::XSD_POSITIVE_INTEGER: return "xsd:positiveInteger";

    case XRefType::XSD_NON_NEGATIVE_INTEGER: return "xsd:nonNegativeInteger";

    case XRefType::XSD_NON_POSITIVE_INTEGER: return "xsd:nonPositiveInteger";

    case XRefType::XSD_BOOLEAN: return "xsd:boolean";

    case XRefType::XSD_DATE: return "xsd:date";

    case XRefType::XSD_ANYURI: return "xsd:anyURI";

    default: return "none";
    }
  }

//  bool ControlledVocabulary::CVTerm::isSearchEngineSpecificScore()
//  { //maybe unsafe?
//    if (this->parents.find("MS:1001143")!=this->parents.end()) return true;
//    return false;
//  }

  bool ControlledVocabulary::CVTerm::isHigherBetterScore(ControlledVocabulary::CVTerm term)
  {
//      for (StringList::const_iterator unp = this->unparsed.begin(); unp != this->unparsed.end(); ++unp)
//      {
//        if (unp->hasPrefix("relationship: has_order MS:1002108")) return true;
//      }
//      return false;
      //most scores are higher better, but most entries in CV for these are not annotated -> default is true
      for (StringList::const_iterator unp = term.unparsed.begin(); unp != term.unparsed.end(); ++unp)
      {
        if (StringUtils::hasPrefix(*unp, "relationship: has_order MS:1002109")) return false;
      }
      return true;
  }

  std::string ControlledVocabulary::CVTerm::toXMLString(const std::string& ref, const std::string& value) const
  {
    std::string s =  "<cvParam accession=\"" + id + "\" cvRef=\"" + ref + "\" name=\"" + Internal::XMLHandler::writeXMLEscape(name);
    if (!value.empty())
    {
      s += "\" value=\"" + Internal::XMLHandler::writeXMLEscape(StringUtils::toStr(value));
    }
    s +=  "\"/>";
    return s;
    //~ TODO: handle unknown cvparams in ControlledVocabulary to get same formatting but more userdefined interface
  }

  std::string ControlledVocabulary::CVTerm::toXMLString(const std::string& ref, const OpenMS::DataValue& value) const
  {
    std::string s =  "<cvParam accession=\"" + id + "\" cvRef=\"" + ref + "\" name=\"" + Internal::XMLHandler::writeXMLEscape(name);
    if (!value.isEmpty())
    {
      s += "\" value=\"" + Internal::XMLHandler::writeXMLEscape(StringUtils::toStr(value));
    }
    if (value.hasUnit())
    {
      std::string un = *(this->units.begin());
      s += "\" unitAccession=\"" + un + "\" unitCvRef=\"" + StringUtils::prefix(un, 2);
      // TODO: Currently we do not store the unit name in the CVTerm, only the
      // accession number (we would need the ControlledVocabulary to look up
      // the unit CVTerm).
      // "\" unitName=\"" + unit.name
    }
    s +=  "\"/>";
    return s;
  }

  ControlledVocabulary::ControlledVocabulary() = default;

  ControlledVocabulary::~ControlledVocabulary() = default;

  void ControlledVocabulary::loadFromOBO(const std::string& name, const std::string& filename)
  {
    bool in_term = false;
    name_ = name;

    ifstream is(filename.c_str());
    if (!is)
    {
      if (!File::exists(filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      else if (!File::readable(filename))
      {
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      else
      {
        throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
    }

    std::string line, line_wo_spaces;
    CVTerm term;

    //parse file
    while (getline(is, line, '\n'))
    {
      StringUtils::trim(line);
      line_wo_spaces = line;
      StringUtils::removeWhitespaces(line_wo_spaces);

      //do nothing for empty lines
      if (line.empty())
      {
        continue;
      }

      if (StringUtils::hasPrefix(line_wo_spaces, "data-version:"))
      {
        version_ = StringUtils::trimmed(StringUtils::substr(line, line.find(':') + 1));
      }
      if (StringUtils::hasPrefix(line_wo_spaces, "default-namespace:"))
      {
        label_ = StringUtils::trimmed(StringUtils::substr(line, line.find(':') + 1));
      }
      if (StringUtils::hasPrefix(line_wo_spaces, "remark:URL:"))
      {
        // Find the position of "http://" or "https://"
        size_t httpPos = line.find("http://");
        size_t httpsPos = line.find("https://");

        // Determine the starting position of the URL
        if (httpPos != std::string::npos) 
        {
          url_ = StringUtils::trimmed(StringUtils::substr(line, httpPos));
        } else if (httpsPos != std::string::npos) 
        {
          url_ = StringUtils::trimmed(StringUtils::substr(line, httpsPos));
        } else 
        {
          // No URL found
          std::cerr << "No URL found in the line." << std::endl;
        }
      }

      //********************************************************************************
      //stanza line
      if (line_wo_spaces[0] == '[')
      {
        //[term] stanza
        if (StringUtils::toLower(line_wo_spaces) == "[term]") //new term
        {
          in_term = true;
          if (!term.id.empty()) //store last term
          {
            terms_[term.id] = term;
          }

          //clear temporary term members
          term = CVTerm();
        }
        // other stanza => not in a term
        else
        {
          in_term = false;
        }
      }
      //********************************************************************************
      //data line
      else if (in_term)
      {
        if (StringUtils::hasPrefix(line_wo_spaces, "id:"))
        {
          term.id = StringUtils::trimmed(StringUtils::substr(line, line.find(':') + 1));
        }
        else if (StringUtils::hasPrefix(line_wo_spaces, "name:"))
        {
          term.name = StringUtils::trimmed(StringUtils::substr(line, line.find(':') + 1));
        }
        else if (StringUtils::hasPrefix(line_wo_spaces, "is_a:"))
        {
          if (StringUtils::has(line, '!'))
          {
            std::string parent_id = StringUtils::trimmed(StringUtils::prefix(StringUtils::substr(line, line.find(':') + 1), '!'));
            term.parents.insert(parent_id);

            //check if the parent term name is correct
            std::string parent_name = StringUtils::trimmed(StringUtils::suffix(line, '!'));
            if (!checkName_(parent_id, parent_name))
              cerr << "Warning: while loading term '" << term.id << "' of CV '" << name_ << "': parent term name '" << parent_name << "' and id '" << parent_id << "' differ." << "\n";
          }
          else
          {
            term.parents.insert(StringUtils::trimmed(StringUtils::substr(line, line.find(':') + 1)));
          }
        }
        // brenda tissue special relationships, DRV (derived and part of)
        else if (StringUtils::hasPrefix(line_wo_spaces, "relationship:DRV") && name == "brenda")
        {
          if (StringUtils::has(line, '!'))
          {
            // e.g. relationship: DRV BTO:0000142 ! brain
            std::string parent_id = StringUtils::prefix(StringUtils::substr(line, line.find("DRV") + 4), ':') + ":" + StringUtils::trimmed(StringUtils::prefix(StringUtils::suffix(line, ':'), '!'));
            term.parents.insert(parent_id);

            //check if the parent term name is correct
            std::string parent_name = StringUtils::trimmed(StringUtils::suffix(line, '!'));
            if (!checkName_(parent_id, parent_name))
              cerr << "Warning: while loading term '" << term.id << "' of CV '" << name_ << "': DRV relationship term name '" << parent_name << "' and id '" << parent_id << "' differ." << "\n";
          }
          else
          {
            // e.g. relationship: DRV BTO:0000142
            term.parents.insert(StringUtils::prefix(StringUtils::substr(line, line.find("DRV") + 4), ':') + ":" + StringUtils::trimmed(StringUtils::suffix(line, ':')));
          }
        }
        else if (StringUtils::hasPrefix(line_wo_spaces, "relationship:part_of") && name == "brenda")
        {
          if (StringUtils::has(line, '!'))
          {
            std::string parent_id = StringUtils::prefix(StringUtils::substr(line, line.find("part_of") + 8), ':') + ":" + StringUtils::trimmed(StringUtils::prefix(StringUtils::suffix(line, ':'), '!'));
            term.parents.insert(parent_id);

            //check if the parent term name is correct
            std::string parent_name = StringUtils::trimmed(StringUtils::suffix(line, '!'));
            if (!checkName_(parent_id, parent_name))
            {
              cerr << "Warning: while loading term '" << term.id << "' of CV '" << name_ << "': part_of relationship term name '" << parent_name << "' and id '" << parent_id << "' differ." << "\n";
            }
          }
          else
          {
            term.parents.insert(StringUtils::prefix(StringUtils::substr(line, line.find("part_of") + 8), ':') + ":" + StringUtils::trimmed(StringUtils::suffix(line, ':')));
          }
        }
        else if (StringUtils::hasPrefix(line_wo_spaces, "relationship:has_units"))
        {
          if (StringUtils::has(line, '!'))
          {
            std::string unit_id = StringUtils::prefix(StringUtils::substr(line, line.find("has_units") + 10), ':') + ":" + StringUtils::trimmed(StringUtils::prefix(StringUtils::suffix(line, ':'), '!'));
            term.units.insert(unit_id);

            //check if the parent term name is correct
            std::string unit_name = StringUtils::trimmed(StringUtils::suffix(line, '!'));
            if (!checkName_(unit_id, unit_name))
            {
              cerr << "Warning: while loading term '" << term.id << "' of CV '" << name_ << "': has_units relationship term name '" << unit_name << "' and id '" << unit_id << "' differ." << "\n";
            }
          }
          else
          {
            term.units.insert(StringUtils::prefix(StringUtils::substr(line, line.find("has_units") + 10), ':') + ":" + StringUtils::trimmed(StringUtils::suffix(line, ':')));
          }
        }
        else if (StringUtils::hasPrefix(line_wo_spaces, "def:"))
        {
          std::string description = StringUtils::substr(line, line.find('"') + 1);
          StringUtils::trim(description);
          description = StringUtils::substr(description, 0, description.find('"'));
          StringUtils::trim(description);
          term.description = description;
        }
        else if (StringUtils::hasPrefix(line_wo_spaces, "synonym:"))
        {
          std::string synonym = StringUtils::substr(line, line.find('"') + 1);
          StringUtils::trim(synonym);
          synonym = StringUtils::substr(synonym, 0, synonym.find('"'));
          StringUtils::trim(synonym);
          term.synonyms.push_back(synonym);
        }
        else if (line_wo_spaces == "is_obsolete:true")
        {
          term.obsolete = true;
        }
        else if (StringUtils::hasPrefix(line_wo_spaces, "xref:value-type") 
          || StringUtils::hasPrefix(line_wo_spaces, "xref_analog:value-type")
        )
        {
          StringUtils::remove(line_wo_spaces, '\\');
          if (StringUtils::hasSubstring(line_wo_spaces, "value-type:xsd:string"))
          {
            term.xref_type = CVTerm::XRefType::XSD_STRING;
            continue;
          }
          if (StringUtils::hasSubstring(line_wo_spaces, "value-type:xsd:integer") || StringUtils::hasSubstring(line_wo_spaces, "value-type:xsd:int"))
          {
            term.xref_type = CVTerm::XRefType::XSD_INTEGER;
            continue;
          }
          if (StringUtils::hasSubstring(line_wo_spaces, "value-type:xsd:decimal") ||
              StringUtils::hasSubstring(line_wo_spaces, "value-type:xsd:float") ||
              StringUtils::hasSubstring(line_wo_spaces, "value-type:xsd:double"))
          {
            term.xref_type = CVTerm::XRefType::XSD_DECIMAL;
            continue;
          }
          if (StringUtils::hasSubstring(line_wo_spaces, "value-type:xsd:negativeInteger"))
          {
            term.xref_type = CVTerm::XRefType::XSD_NEGATIVE_INTEGER;
            continue;
          }
          if (StringUtils::hasSubstring(line_wo_spaces, "value-type:xsd:positiveInteger"))
          {
            term.xref_type = CVTerm::XRefType::XSD_POSITIVE_INTEGER;
            continue;
          }
          if (StringUtils::hasSubstring(line_wo_spaces, "value-type:xsd:nonNegativeInteger"))
          {
            term.xref_type = CVTerm::XRefType::XSD_NON_NEGATIVE_INTEGER;
            continue;
          }
          if (StringUtils::hasSubstring(line_wo_spaces, "value-type:xsd:nonPositiveInteger"))
          {
            term.xref_type = CVTerm::XRefType::XSD_NON_POSITIVE_INTEGER;
            continue;
          }
          if (StringUtils::hasSubstring(line_wo_spaces, "value-type:xsd:boolean") || StringUtils::hasSubstring(line_wo_spaces, "value-type:xsd:bool"))
          {
            term.xref_type = CVTerm::XRefType::XSD_BOOLEAN;
            continue;
          }
          if (StringUtils::hasSubstring(line_wo_spaces, "value-type:xsd:date"))
          {
            term.xref_type = CVTerm::XRefType::XSD_DATE;
            continue;
          }
          if (StringUtils::hasSubstring(line_wo_spaces, "value-type:xsd:anyURI"))
          {
            term.xref_type = CVTerm::XRefType::XSD_ANYURI;
            continue;
          }
          cerr << "ControlledVocabulary: OBOFile: unknown xsd type: " << line_wo_spaces << ", ignoring" << "\n";
        }
        else if (StringUtils::hasPrefix(line_wo_spaces, "relationship:has_value_type")) // since newer obo type in relationship instead of xref
        {
          if (StringUtils::hasSubstring(line_wo_spaces, "xsd:string"))
          {
            term.xref_type = CVTerm::XRefType::XSD_STRING;
            continue;
          }
          if (StringUtils::hasSubstring(line_wo_spaces, "xsd:integer") 
          || StringUtils::hasSubstring(line_wo_spaces, "xsd:int"))
          {
            term.xref_type = CVTerm::XRefType::XSD_INTEGER;
            continue;
          }
          if (StringUtils::hasSubstring(line_wo_spaces, "xsd:decimal") ||
              StringUtils::hasSubstring(line_wo_spaces, "xsd:float") ||
              StringUtils::hasSubstring(line_wo_spaces, "xsd:double"))
          {
            term.xref_type = CVTerm::XRefType::XSD_DECIMAL;
            continue;
          }
          if (StringUtils::hasSubstring(line_wo_spaces, "xsd:negativeInteger"))
          {
            term.xref_type = CVTerm::XRefType::XSD_NEGATIVE_INTEGER;
            continue;
          }
          if (StringUtils::hasSubstring(line_wo_spaces, "xsd:positiveInteger"))
          {
            term.xref_type = CVTerm::XRefType::XSD_POSITIVE_INTEGER;
            continue;
          }
          if (StringUtils::hasSubstring(line_wo_spaces, "xsd:nonNegativeInteger"))
          {
            term.xref_type = CVTerm::XRefType::XSD_NON_NEGATIVE_INTEGER;
            continue;
          }
          if (StringUtils::hasSubstring(line_wo_spaces, "xsd:nonPositiveInteger"))
          {
            term.xref_type = CVTerm::XRefType::XSD_NON_POSITIVE_INTEGER;
            continue;
          }
          if (StringUtils::hasSubstring(line_wo_spaces, "xsd:boolean") 
          || StringUtils::hasSubstring(line_wo_spaces, "xsd:bool"))
          {
            term.xref_type = CVTerm::XRefType::XSD_BOOLEAN;
            continue;
          }
          if (StringUtils::hasSubstring(line_wo_spaces, "xsd:date"))
          {
            term.xref_type = CVTerm::XRefType::XSD_DATE;
            continue;
          }
          if (StringUtils::hasSubstring(line_wo_spaces, "xsd:anyURI"))
          {
            term.xref_type = CVTerm::XRefType::XSD_ANYURI;
            continue;
          }
          if (
            StringUtils::hasSubstring(line_wo_spaces, "MS:1002711") ||
            StringUtils::hasSubstring(line_wo_spaces, "MS:1002712") ||
            StringUtils::hasSubstring(line_wo_spaces, "MS:1002713")
          )
          {
            term.xref_type = CVTerm::XRefType::XSD_STRING; // store list as string
            continue;
          }
          cerr << "ControlledVocabulary: OBOFile: unknown xsd type: " << line_wo_spaces << ", ignoring" << "\n";
        }       
        else if (StringUtils::hasPrefix(line_wo_spaces, "xref:binary-data-type") || StringUtils::hasPrefix(line_wo_spaces, "xref_analog:binary-data-type"))
        {
          StringUtils::remove(line_wo_spaces, '\\');
          //remove description (if present)
          // according to rev1165 of the cv comments are here quoted, see http://psidev.cvs.sourceforge.net/viewvc/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo?revision=1.165&view=markup
          if (StringUtils::has(line_wo_spaces, '\"'))
          {
            line_wo_spaces = StringUtils::substr(line_wo_spaces, 0, line_wo_spaces.find('\"'));
          }
          //trim prefix
          line_wo_spaces = StringUtils::substr(line_wo_spaces, 22);
          //trim just to be sure
          StringUtils::trim(line_wo_spaces);
          term.xref_binary.push_back(line_wo_spaces);
        }
        else if (!line.empty())
        {
          term.unparsed.push_back(line);
        }
      }
    }

    if (!term.id.empty()) //store last term
    {
      terms_[term.id] = term;
    }

    // now build all child terms
    for (auto it = terms_.begin(); it != terms_.end(); ++it)
    {
      //cerr << it->first << "\n";
      for (auto pit = it->second.parents.begin(); pit != it->second.parents.end(); ++pit)
      {
        //cerr << "Parent: " << *pit << "\n";
        terms_[*pit].children.insert(it->first);
      }

      auto mit = namesToIds_.find(it->second.name);
      if (mit == namesToIds_.end())
      {
        namesToIds_.insert(pair<std::string, std::string>(it->second.name, it->first));
      }
      else
      {
        //~ TODO that case would be bad do something
        std::string s = it->second.name + it->second.description;
        namesToIds_.insert(pair<std::string, std::string>(s, it->first));
      }
    }
  }

  const ControlledVocabulary::CVTerm& ControlledVocabulary::getTerm(const std::string& id) const
  {
    if (const auto it = terms_.find(id); it == terms_.end())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid CV identifier!", id);
    }
    else
    {
      return it->second;
    }
  }

  const std::map<std::string, ControlledVocabulary::CVTerm>& ControlledVocabulary::getTerms() const
  {
    return terms_;
  }

  void ControlledVocabulary::getAllChildTerms(set<std::string>& terms, const std::string& parent) const
  {
    //cerr << "Parent: " << parent << "\n";
    for (const auto& child : getTerm(parent).children)
    {
      terms.insert(child);
      //TODO: This is not safe for cyclic graphs. Are they allowed in CVs?
      getAllChildTerms(terms, child);
    }
  }

  void ControlledVocabulary::addAllChildTerms(set<std::string>& terms, const std::string& parent) const
  {
    getTerm(parent); // validate before mutating caller-owned output (throws if parent is unknown)
    terms.insert(parent);
    getAllChildTerms(terms, parent);
  }

  const ControlledVocabulary::CVTerm& ControlledVocabulary::getTermByName(const std::string& name, const std::string& desc) const
  {
    //slow, but Vocabulary is very finite and this method will be called only a few times during write of a ML file using a CV
    auto it = namesToIds_.find(name);
    if (it == namesToIds_.end())
    {
      if (!desc.empty())
      {
        it = namesToIds_.find(StringUtils::toStr(name + desc));
        if (it == namesToIds_.end())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid CV name!", name);
        }
      }
      else
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid CV name!", name);
      }
    }

    return terms_.at(it->second);
  }

  bool ControlledVocabulary::exists(const std::string& id) const
  {
    return terms_.contains(id);
  }

  const ControlledVocabulary::CVTerm* ControlledVocabulary::checkAndGetTermByName(const std::string& name) const
  {
    const auto it = namesToIds_.find(name);
    if (it == namesToIds_.end()) return nullptr;
    return &terms_.at(it->second);
  }

  bool ControlledVocabulary::hasTermWithName(const std::string& name) const
  {
    const auto it = namesToIds_.find(name);
    return it != namesToIds_.end();
  }

  bool ControlledVocabulary::isChildOf(const std::string& child, const std::string& parent) const
  {
    // cout << "CHECK child:" << child << " parent: " << parent << "\n";
    const CVTerm& ch = getTerm(child);

    for (const auto & it : ch.parents)
    {
      // cout << "Parent: " << it << "\n";

      // check if it is a direct parent
      if (it == parent)
      {
        return true;
      }
      // check if it is an indirect parent
      else if (isChildOf(it, parent))
      {
        return true;
      }
    }

    return false;
  }

  std::ostream& operator<<(std::ostream& os, const ControlledVocabulary& cv)
  {
    for (const auto & it : cv.terms_)
    {
      os << "[Term]\n";
      os << "id: '" << it.second.id << "'\n";
      os << "name: '" << it.second.name <<  "'\n";
      for (const auto & parent_term : it.second.parents)
      {
        cout << "is_a: '" << parent_term <<  "'\n";
      }
    }
    return os;
  }

  const std::string& ControlledVocabulary::name() const
  {
    return name_;
  }

  const std::string& ControlledVocabulary::label() const
  {
    return label_;
  }

  const std::string& ControlledVocabulary::version() const
  {
    return version_;
  }

  const std::string& ControlledVocabulary::url() const
  {
    return url_;
  }

  const ControlledVocabulary& ControlledVocabulary::getPSIMSCV()
  {
    static const ControlledVocabulary cv = []() {
      ControlledVocabulary cv;
      cv.loadFromOBO("MS", File::find("/CV/psi-ms.obo"));
      cv.loadFromOBO("PATO", File::find("/CV/quality.obo"));
      cv.loadFromOBO("UO", File::find("/CV/unit.obo"));
      cv.loadFromOBO("BTO", File::find("/CV/brenda.obo"));
      cv.loadFromOBO("GO", File::find("/CV/goslim_goa.obo"));
      return cv;
    }();
    return cv;
  }

  bool ControlledVocabulary::checkName_(const std::string& id, const std::string& name, bool ignore_case) const
  {
    if (!exists(id))
    {
      return true; //what?!
    }
    std::string parent_name = name;
    std::string real_parent_name = getTerm(id).name;

    if (ignore_case)
    {
      StringUtils::toLower(parent_name);
      StringUtils::toLower(real_parent_name);
    }

    return real_parent_name == parent_name;
  }

} // namespace OpenMS
