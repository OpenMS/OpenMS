// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
    xref_type(NONE),
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

  String ControlledVocabulary::CVTerm::getXRefTypeName(XRefType type)
  {
    switch (type)
    {
    case XSD_STRING: return "xsd:string";

    case XSD_INTEGER: return "xsd:integer";

    case XSD_DECIMAL: return "xsd:decimal";

    case XSD_NEGATIVE_INTEGER: return "xsd:negativeInteger";

    case XSD_POSITIVE_INTEGER: return "xsd:positiveInteger";

    case XSD_NON_NEGATIVE_INTEGER: return "xsd:nonNegativeInteger";

    case XSD_NON_POSITIVE_INTEGER: return "xsd:nonPositiveInteger";

    case XSD_BOOLEAN: return "xsd:boolean";

    case XSD_DATE: return "xsd:date";

    case XSD_ANYURI: return "xsd:anyURI";

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
        if (unp->hasPrefix("relationship: has_order MS:1002109")) return false;
      }
      return true;
  }

  String ControlledVocabulary::CVTerm::toXMLString(const OpenMS::String& ref, const String& value) const
  {
    String s =  "<cvParam accession=\"" + id + "\" cvRef=\"" + ref + "\" name=\"" + Internal::XMLHandler::writeXMLEscape(name);
    if (!value.empty())
    {
      s += "\" value=\"" + Internal::XMLHandler::writeXMLEscape(value);
    }
    s +=  "\"/>";
    return s;
    //~ TODO: handle unknown cvparams in ControlledVocabulary to get same formatting but more userdefined interface
  }

  String ControlledVocabulary::CVTerm::toXMLString(const OpenMS::String& ref, const OpenMS::DataValue& value) const
  {
    String s =  "<cvParam accession=\"" + id + "\" cvRef=\"" + ref + "\" name=\"" + Internal::XMLHandler::writeXMLEscape(name);
    if (!value.isEmpty())
    {
      s += "\" value=\"" + Internal::XMLHandler::writeXMLEscape(value);
    }
    if (value.hasUnit())
    {
      String un = *(this->units.begin());
      s += "\" unitAccession=\"" + un + "\" unitCvRef=\"" + un.prefix(2);
      // TODO: Currently we do not store the unit name in the CVTerm, only the
      // accession number (we would need the ControlledVocabulary to look up
      // the unit CVTerm).
      // "\" unitName=\"" + unit.name
    }
    s +=  "\"/>";
    return s;
  }

  ControlledVocabulary::ControlledVocabulary() :
    terms_(),
    name_("")
  {

  }

  ControlledVocabulary::~ControlledVocabulary() = default;

  void ControlledVocabulary::loadFromOBO(const String& name, const String& filename)
  {
    bool in_term = false;
    name_ = name;

    ifstream is(filename.c_str());
    if (!is)
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    String line, line_wo_spaces;
    CVTerm term;

    //parse file
    while (getline(is, line, '\n'))
    {
      line.trim();
      line_wo_spaces = line;
      line_wo_spaces.removeWhitespaces();

      //do nothing for empty lines
      if (line.empty())
      {
        continue;
      }

      if (line_wo_spaces.hasPrefix("data-version:"))
      {
        version_ = line.substr(line.find(':') + 1).trim();
      }
      if (line_wo_spaces.hasPrefix("default-namespace:"))
      {
        label_ = line.substr(line.find(':') + 1).trim();
      }
      if (line_wo_spaces.hasPrefix("remark:URL:"))
      {
        // Find the position of "http://" or "https://"
        size_t httpPos = line.find("http://");
        size_t httpsPos = line.find("https://");

        // Determine the starting position of the URL
        if (httpPos != std::string::npos) 
        {
          url_ = line.substr(httpPos).trim();
        } else if (httpsPos != std::string::npos) 
        {
          url_ = line.substr(httpsPos).trim();
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
        if (line_wo_spaces.toLower() == "[term]") //new term
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
        if (line_wo_spaces.hasPrefix("id:"))
        {
          term.id = line.substr(line.find(':') + 1).trim();
        }
        else if (line_wo_spaces.hasPrefix("name:"))
        {
          term.name = line.substr(line.find(':') + 1).trim();
        }
        else if (line_wo_spaces.hasPrefix("is_a:"))
        {
          if (line.has('!'))
          {
            String parent_id = line.substr(line.find(':') + 1).prefix('!').trim();
            term.parents.insert(parent_id);

            //check if the parent term name is correct
            String parent_name = line.suffix('!').trim();
            if (!checkName_(parent_id, parent_name))
              cerr << "Warning: while loading term '" << term.id << "' of CV '" << name_ << "': parent term name '" << parent_name << "' and id '" << parent_id << "' differ." << "\n";
          }
          else
          {
            term.parents.insert(line.substr(line.find(':') + 1).trim());
          }
        }
        // brenda tissue special relationships, DRV (derived and part of)
        else if (line_wo_spaces.hasPrefix("relationship:DRV") && name == "brenda")
        {
          if (line.has('!'))
          {
            // e.g. relationship: DRV BTO:0000142 ! brain
            String parent_id = line.substr(line.find("DRV") + 4).prefix(':') + ":" + line.suffix(':').prefix('!').trim();
            term.parents.insert(parent_id);

            //check if the parent term name is correct
            String parent_name = line.suffix('!').trim();
            if (!checkName_(parent_id, parent_name))
              cerr << "Warning: while loading term '" << term.id << "' of CV '" << name_ << "': DRV relationship term name '" << parent_name << "' and id '" << parent_id << "' differ." << "\n";
          }
          else
          {
            // e.g. relationship: DRV BTO:0000142
            term.parents.insert(line.substr(line.find("DRV") + 4).prefix(':') + ":" + line.suffix(':').trim());
          }
        }
        else if (line_wo_spaces.hasPrefix("relationship:part_of") && name == "brenda")
        {
          if (line.has('!'))
          {
            String parent_id = line.substr(line.find("part_of") + 8).prefix(':') + ":" + line.suffix(':').prefix('!').trim();
            term.parents.insert(parent_id);

            //check if the parent term name is correct
            String parent_name = line.suffix('!').trim();
            if (!checkName_(parent_id, parent_name))
            {
              cerr << "Warning: while loading term '" << term.id << "' of CV '" << name_ << "': part_of relationship term name '" << parent_name << "' and id '" << parent_id << "' differ." << "\n";
            }
          }
          else
          {
            term.parents.insert(line.substr(line.find("part_of") + 8).prefix(':') + ":" + line.suffix(':').trim());
          }
        }
        else if (line_wo_spaces.hasPrefix("relationship:has_units"))
        {
          if (line.has('!'))
          {
            String unit_id = line.substr(line.find("has_units") + 10).prefix(':') + ":" + line.suffix(':').prefix('!').trim();
            term.units.insert(unit_id);

            //check if the parent term name is correct
            String unit_name = line.suffix('!').trim();
            if (!checkName_(unit_id, unit_name))
            {
              cerr << "Warning: while loading term '" << term.id << "' of CV '" << name_ << "': has_units relationship term name '" << unit_name << "' and id '" << unit_id << "' differ." << "\n";
            }
          }
          else
          {
            term.units.insert(line.substr(line.find("has_units") + 10).prefix(':') + ":" + line.suffix(':').trim());
          }
        }
        else if (line_wo_spaces.hasPrefix("def:"))
        {
          String description = line.substr(line.find('"') + 1);
          description.trim();
          description = description.substr(0, description.find('"'));
          description.trim();
          term.description = description;
        }
        else if (line_wo_spaces.hasPrefix("synonym:"))
        {
          String synonym = line.substr(line.find('"') + 1);
          synonym.trim();
          synonym = synonym.substr(0, synonym.find('"'));
          synonym.trim();
          term.synonyms.push_back(synonym);
        }
        else if (line_wo_spaces == "is_obsolete:true")
        {
          term.obsolete = true;
        }
        else if (line_wo_spaces.hasPrefix("xref:value-type") 
          || line_wo_spaces.hasPrefix("xref_analog:value-type")
        )
        {
          line_wo_spaces.remove('\\');
          if (line_wo_spaces.hasSubstring("value-type:xsd:string"))
          {
            term.xref_type = CVTerm::XSD_STRING;
            continue;
          }
          if (line_wo_spaces.hasSubstring("value-type:xsd:integer") || line_wo_spaces.hasSubstring("value-type:xsd:int"))
          {
            term.xref_type = CVTerm::XSD_INTEGER;
            continue;
          }
          if (line_wo_spaces.hasSubstring("value-type:xsd:decimal") ||
              line_wo_spaces.hasSubstring("value-type:xsd:float") ||
              line_wo_spaces.hasSubstring("value-type:xsd:double"))
          {
            term.xref_type = CVTerm::XSD_DECIMAL;
            continue;
          }
          if (line_wo_spaces.hasSubstring("value-type:xsd:negativeInteger"))
          {
            term.xref_type = CVTerm::XSD_NEGATIVE_INTEGER;
            continue;
          }
          if (line_wo_spaces.hasSubstring("value-type:xsd:positiveInteger"))
          {
            term.xref_type = CVTerm::XSD_POSITIVE_INTEGER;
            continue;
          }
          if (line_wo_spaces.hasSubstring("value-type:xsd:nonNegativeInteger"))
          {
            term.xref_type = CVTerm::XSD_NON_NEGATIVE_INTEGER;
            continue;
          }
          if (line_wo_spaces.hasSubstring("value-type:xsd:nonPositiveInteger"))
          {
            term.xref_type = CVTerm::XSD_NON_POSITIVE_INTEGER;
            continue;
          }
          if (line_wo_spaces.hasSubstring("value-type:xsd:boolean") || line_wo_spaces.hasSubstring("value-type:xsd:bool"))
          {
            term.xref_type = CVTerm::XSD_BOOLEAN;
            continue;
          }
          if (line_wo_spaces.hasSubstring("value-type:xsd:date"))
          {
            term.xref_type = CVTerm::XSD_DATE;
            continue;
          }
          if (line_wo_spaces.hasSubstring("value-type:xsd:anyURI"))
          {
            term.xref_type = CVTerm::XSD_ANYURI;
            continue;
          }
          cerr << "ControlledVocabulary: OBOFile: unknown xsd type: " << line_wo_spaces << ", ignoring" << "\n";
        }
        else if (line_wo_spaces.hasPrefix("relationship:has_value_type")) // since newer obo type in relationship instead of xref
        {
          if (line_wo_spaces.hasSubstring("xsd:string"))
          {
            term.xref_type = CVTerm::XSD_STRING;
            continue;
          }
          if (line_wo_spaces.hasSubstring("xsd:integer") 
          || line_wo_spaces.hasSubstring("xsd:int"))
          {
            term.xref_type = CVTerm::XSD_INTEGER;
            continue;
          }
          if (line_wo_spaces.hasSubstring("xsd:decimal") ||
              line_wo_spaces.hasSubstring("xsd:float") ||
              line_wo_spaces.hasSubstring("xsd:double"))
          {
            term.xref_type = CVTerm::XSD_DECIMAL;
            continue;
          }
          if (line_wo_spaces.hasSubstring("xsd:negativeInteger"))
          {
            term.xref_type = CVTerm::XSD_NEGATIVE_INTEGER;
            continue;
          }
          if (line_wo_spaces.hasSubstring("xsd:positiveInteger"))
          {
            term.xref_type = CVTerm::XSD_POSITIVE_INTEGER;
            continue;
          }
          if (line_wo_spaces.hasSubstring("xsd:nonNegativeInteger"))
          {
            term.xref_type = CVTerm::XSD_NON_NEGATIVE_INTEGER;
            continue;
          }
          if (line_wo_spaces.hasSubstring("xsd:nonPositiveInteger"))
          {
            term.xref_type = CVTerm::XSD_NON_POSITIVE_INTEGER;
            continue;
          }
          if (line_wo_spaces.hasSubstring("xsd:boolean") 
          || line_wo_spaces.hasSubstring("xsd:bool"))
          {
            term.xref_type = CVTerm::XSD_BOOLEAN;
            continue;
          }
          if (line_wo_spaces.hasSubstring("xsd:date"))
          {
            term.xref_type = CVTerm::XSD_DATE;
            continue;
          }
          if (line_wo_spaces.hasSubstring("xsd:anyURI"))
          {
            term.xref_type = CVTerm::XSD_ANYURI;
            continue;
          }
          if (
            line_wo_spaces.hasSubstring("MS:1002711") ||
            line_wo_spaces.hasSubstring("MS:1002712") ||
            line_wo_spaces.hasSubstring("MS:1002713")
          )
          {
            term.xref_type = CVTerm::XSD_STRING; // store list as string
            continue;
          }
          cerr << "ControlledVocabulary: OBOFile: unknown xsd type: " << line_wo_spaces << ", ignoring" << "\n";
        }       
        else if (line_wo_spaces.hasPrefix("xref:binary-data-type") || line_wo_spaces.hasPrefix("xref_analog:binary-data-type"))
        {
          line_wo_spaces.remove('\\');
          //remove description (if present)
          // according to rev1165 of the cv comments are here quoted, see http://psidev.cvs.sourceforge.net/viewvc/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo?revision=1.165&view=markup
          if (line_wo_spaces.has('\"'))
          {
            line_wo_spaces = line_wo_spaces.substr(0, line_wo_spaces.find('\"'));
          }
          //trim prefix
          line_wo_spaces = line_wo_spaces.substr(22);
          //trim just to be sure
          line_wo_spaces.trim();
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
    for (std::map<String, CVTerm>::iterator it = terms_.begin(); it != terms_.end(); ++it)
    {
      //cerr << it->first << "\n";
      for (set<String>::const_iterator pit = it->second.parents.begin(); pit != it->second.parents.end(); ++pit)
      {
        //cerr << "Parent: " << *pit << "\n";
        terms_[*pit].children.insert(it->first);
      }

      std::map<String, String>::iterator mit = namesToIds_.find(it->second.name);
      if (mit == namesToIds_.end())
      {
        namesToIds_.insert(pair<String, String>(it->second.name, it->first));
      }
      else
      {
        //~ TODO that case would be bad do something
        String s = it->second.name + it->second.description;
        namesToIds_.insert(pair<String, String>(s, it->first));
      }
    }
  }

  const ControlledVocabulary::CVTerm& ControlledVocabulary::getTerm(const String& id) const
  {
    std::map<String, CVTerm>::const_iterator it = terms_.find(id);
    if (it == terms_.end())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid CV identifier!", id);
    }
    return it->second;
  }

  const std::map<String, ControlledVocabulary::CVTerm>& ControlledVocabulary::getTerms() const
  {
    return terms_;
  }

  void ControlledVocabulary::getAllChildTerms(set<String>& terms, const String& parent) const
  {
    //cerr << "Parent: " << parent << "\n";
    for (const auto& child : getTerm(parent).children)
    {
      terms.insert(child);
      //TODO: This is not safe for cyclic graphs. Are they allowed in CVs?
      getAllChildTerms(terms, child);
    }
  }

  const ControlledVocabulary::CVTerm& ControlledVocabulary::getTermByName(const String& name, const String& desc) const
  {
    //slow, but Vocabulary is very finite and this method will be called only a few times during write of a ML file using a CV
    std::map<String, String>::const_iterator it = namesToIds_.find(name);
    if (it == namesToIds_.end())
    {
      if (!desc.empty())
      {
        it = namesToIds_.find(String(name + desc));
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

  bool ControlledVocabulary::exists(const String& id) const
  {
    return terms_.find(id) != terms_.end();
  }

  const ControlledVocabulary::CVTerm* ControlledVocabulary::checkAndGetTermByName(const OpenMS::String& name) const
  {
    std::map<String, String>::const_iterator it = namesToIds_.find(name);
    if (it == namesToIds_.end()) return nullptr;
    return &terms_.at(it->second);
  }

  bool ControlledVocabulary::hasTermWithName(const OpenMS::String& name) const
  {
    std::map<String, String>::const_iterator it = namesToIds_.find(name);
    return it != namesToIds_.end();
  }

  bool ControlledVocabulary::isChildOf(const String& child, const String& parent) const
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

  const String& ControlledVocabulary::name() const
  {
    return name_;
  }

  const String& ControlledVocabulary::label() const
  {
    return label_;
  }

  const String& ControlledVocabulary::version() const
  {
    return version_;
  }

  const String& ControlledVocabulary::url() const
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

  bool ControlledVocabulary::checkName_(const String& id, const String& name, bool ignore_case) const
  {
    if (!exists(id))
    {
      return true; //what?!
    }
    String parent_name = name;
    String real_parent_name = getTerm(id).name;

    if (ignore_case)
    {
      parent_name.toLower();
      real_parent_name.toLower();
    }

    return real_parent_name == parent_name;
  }

} // namespace OpenMS
