// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/VALIDATORS/SemanticValidator.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/DATASTRUCTURES/CVMappingRule.h>
#include <OpenMS/DATASTRUCTURES/CVMappingTerm.h>
#include <map>

using namespace xercesc;
using namespace std;

namespace OpenMS::Internal
{

    SemanticValidator::SemanticValidator(const CVMappings& mapping, const ControlledVocabulary& cv) :
      XMLHandler("", 0),
      XMLFile(),
      mapping_(mapping),
      cv_(cv),
      open_tags_(),
      cv_tag_("cvParam"),
      accession_att_("accession"),
      name_att_("name"),
      value_att_("value"),
      unit_accession_att_("unitAccession"),
      unit_name_att_("unitName"),
      check_term_value_types_(true),
      check_units_(false)
    {
      //order rules by element
      for (Size r = 0; r < mapping_.getMappingRules().size(); ++r)
      {
        rules_[mapping_.getMappingRules()[r].getElementPath()].push_back(mapping_.getMappingRules()[r]);
      }
    }

    SemanticValidator::~SemanticValidator()
    = default;

    void SemanticValidator::setTag(const String& tag)
    {
      cv_tag_ = tag;
    }

    void SemanticValidator::setAccessionAttribute(const String& accession)
    {
      accession_att_ = accession;
    }

    void SemanticValidator::setNameAttribute(const String& name)
    {
      name_att_ = name;
    }

    void SemanticValidator::setValueAttribute(const String& value)
    {
      value_att_ = value;
    }

    void SemanticValidator::setCheckTermValueTypes(bool check)
    {
      check_term_value_types_ = check;
    }

    void SemanticValidator::setCheckUnits(bool check)
    {
      check_units_ = check;
    }

    void SemanticValidator::setUnitAccessionAttribute(const String& accession)
    {
      unit_accession_att_ = accession;
    }

    void SemanticValidator::setUnitNameAttribute(const String& name)
    {
      unit_name_att_ = name;
    }

    bool SemanticValidator::validate(const String& filename, StringList& errors, StringList& warnings)
    {
      //TODO Check if all required CVs are loaded => exception if not

      //try to open file
      if (!File::exists(filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }

      //initialize
      errors_.clear();
      warnings_.clear();

      //parse
      file_ = filename;
      parse_(filename, this);

      //set output
      errors = errors_;
      warnings = warnings_;

      return errors_.empty();
    }

    void SemanticValidator::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const Attributes& attributes)
    {
      String tag = sm_.convert(qname);
      String path = getPath_() + "/" + cv_tag_ + "/@" + accession_att_;
      open_tags_.push_back(tag);

      if (tag == cv_tag_)
      {
        //extract accession, name and value
        CVTerm parsed_term;
        getCVTerm_(attributes, parsed_term);

        //check if the term is unknown
        if (!cv_.exists(parsed_term.accession))
        {
          warnings_.push_back(String("Unknown CV term: '") + parsed_term.accession + " - " + parsed_term.name + "' at element '" + getPath_(1) + "'");
          return;
        }

        //check if the term is obsolete
        if (cv_.getTerm(parsed_term.accession).obsolete)
        {
          warnings_.push_back(String("Obsolete CV term: '") + parsed_term.accession + " - " + parsed_term.name + "' at element '" + getPath_(1) + "'");
        }

        //actual handling of the term
        handleTerm_(path, parsed_term);
      }
    }

    void SemanticValidator::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
    {
      String tag = sm_.convert(qname);
      String path = getPath_() + "/" + cv_tag_ + "/@" + accession_att_;

      //look up rules and fulfilled rules/terms
      vector<CVMappingRule>& rules = rules_[path];
      std::map<String, std::map<String, UInt> >& fulfilled = fulfilled_[path]; //(rule ID => term ID => term count)

      //check how often each term appeared
      for (Size r = 0; r < rules.size(); ++r)
      {
        for (Size t = 0; t < rules[r].getCVTerms().size(); ++t)
        {
          if (!rules[r].getCVTerms()[t].getIsRepeatable() && fulfilled[rules[r].getIdentifier()][rules[r].getCVTerms()[t].getAccession()] > 1)
          {
            errors_.push_back(String("Violated mapping rule '") + rules[r].getIdentifier() + "' number of term repeats at element '" + getPath_() + "'");
          }
        }
      }

      //check if all required rules are fulfilled
      for (Size r = 0; r < rules.size(); ++r)
      {
        //Count the number of distinct matched terms
        Size terms_count = rules[r].getCVTerms().size();
        UInt match_count = 0;
        for (Size t = 0; t < terms_count; ++t)
        {
          if (fulfilled[rules[r].getIdentifier()][rules[r].getCVTerms()[t].getAccession()] >= 1)
          {
            ++match_count;
          }
        }

        //MUST / AND - all terms must be matched
        if (rules[r].getRequirementLevel() == CVMappingRule::MUST && rules[r].getCombinationsLogic() == CVMappingRule::AND)
        {
          if (match_count != terms_count)
          {
            errors_.push_back(String("Violated mapping rule '") + rules[r].getIdentifier() + "' at element '" + getPath_() + "', " + String(terms_count) + " term(s) should be present, " + String(match_count) + " found!");
          }
        }
        //MUST / OR - at least one terms must be matched
        else if (rules[r].getRequirementLevel() == CVMappingRule::MUST && rules[r].getCombinationsLogic() == CVMappingRule::OR)
        {
          if (match_count == 0)
          {
            errors_.push_back(String("Violated mapping rule '") + rules[r].getIdentifier() + "' at element '" + getPath_() + "', at least one term must be present!");
          }
        }
        //MUST / XOR - exactly one term must be matched
        else if (rules[r].getRequirementLevel() == CVMappingRule::MUST && rules[r].getCombinationsLogic() == CVMappingRule::XOR)
        {
          if (match_count != 1)
          {
            errors_.push_back(String("Violated mapping rule '") + rules[r].getIdentifier() + "' at element '" + getPath_() + "' exactly one of the allowed terms must be used!");
          }
        }
        //MAY(SHOULD) / AND - none or all terms must be matched
        else if (rules[r].getRequirementLevel() != CVMappingRule::SHOULD && rules[r].getCombinationsLogic() == CVMappingRule::AND)
        {
          if (match_count != 0 && match_count != terms_count)
          {
            errors_.push_back(String("Violated mapping rule '") + rules[r].getIdentifier() + "' at element '" + getPath_() + "', if any, all terms must be given!");
          }
        }
        //MAY(SHOULD) / XOR - zero or one terms must be matched
        else if (rules[r].getRequirementLevel() != CVMappingRule::SHOULD && rules[r].getCombinationsLogic() == CVMappingRule::XOR)
        {
          if (match_count > 1)
          {
            errors_.push_back(String("Violated mapping rule '") + rules[r].getIdentifier() + "' at element '" + getPath_() + "', if any, only exactly one of the allowed terms can be used!");
          }
        }
        //MAY(SHOULD) / OR - none to all terms must be matched => always true
      }

      //clear fulfilled rules
      fulfilled_.erase(path);

      open_tags_.pop_back();
    }

    void SemanticValidator::characters(const XMLCh* const /*chars*/, const XMLSize_t /*length*/)
    {
      //nothing to do here
    }

    String SemanticValidator::getPath_(UInt remove_from_end) const
    {
      String path;
      path.concatenate(open_tags_.begin(), open_tags_.end() - remove_from_end, "/");
      path = String("/") + path;
      return path;
    }

    void SemanticValidator::getCVTerm_(const Attributes& attributes, CVTerm& parsed_term)
    {
      parsed_term.accession = attributeAsString_(attributes, accession_att_.c_str());
      parsed_term.name = attributeAsString_(attributes, name_att_.c_str());
      parsed_term.has_value = optionalAttributeAsString_(parsed_term.value, attributes, value_att_.c_str());
      if (check_units_)
      {
        parsed_term.has_unit_accession = optionalAttributeAsString_(parsed_term.unit_accession, attributes, unit_accession_att_.c_str());
        parsed_term.has_unit_name = optionalAttributeAsString_(parsed_term.unit_name, attributes, unit_name_att_.c_str());
      }
      else
      {
        parsed_term.has_unit_accession = false;
        parsed_term.has_unit_name = false;
      }
    }

    //~ void SemanticValidator::makeCVTerm_(const ControlledVocabulary::CVTerm lc, CVTerm & parsed_term)
    //~ {
    //~ parsed_term.accession = lc.id;
    //~ parsed_term.name = lc.name;
    //~ //no value in ControlledVocabulary::CVTerm
    //~ //no units either yet
    //~ {
    //~ parsed_term.has_unit_accession = false;
    //~ parsed_term.has_unit_name = false;
    //~ }
    //~ }

    //reimplemented to
    // - ignore values (not known)
    // - allow more names (upper-lower-case + spaces)
    void SemanticValidator::handleTerm_(const String& path, const CVTerm& parsed_term)
    {
      //check if the term is allowed in this element
      //and if there is a mapping rule for this element
      //Also store fulfilled rule term counts - this count is used to check of the MUST/MAY and AND/OR/XOR is fulfilled
      bool allowed = false;
      bool rule_found = false;
      vector<CVMappingRule>& rules = rules_[path];
      for (Size r = 0; r < rules.size(); ++r) //go thru all rules
      {
        rule_found = true;
        for (Size t = 0; t < rules[r].getCVTerms().size(); ++t) //go thru all terms
        {
          const CVMappingTerm& term = rules[r].getCVTerms()[t];
          if (term.getUseTerm() && term.getAccession() == parsed_term.accession) //check if the term itself is allowed
          {
            allowed = true;
            fulfilled_[path][rules[r].getIdentifier()][term.getAccession()]++;
            break;
          }

          UInt& counter = fulfilled_[path][rules[r].getIdentifier()][term.getAccession()];
          auto searcher = [&allowed, &counter, &parsed_term] (const String& child)
          {
            if (child == parsed_term.accession)
            {
              allowed = true;
              ++counter;
              return true;
            }
            return false;
          };
          if (term.getAllowChildren() && cv_.iterateAllChildren(term.getAccession(), searcher)) //check if the term's children are allowed
          {
            break;
          }
        }
      }

      // check units
      if (check_units_ && cv_.exists(parsed_term.accession))
      {
        ControlledVocabulary::CVTerm term = cv_.getTerm(parsed_term.accession);
        // check if the cv term has units
        if (!term.units.empty())
        {
          if (!parsed_term.has_unit_accession)
          {
            errors_.push_back(String("CV term must have a unit: " + parsed_term.accession + " - " + parsed_term.name));
          }
          else
          {
            // check if the accession is ok
            if (cv_.exists(parsed_term.unit_accession))
            {
              // check whether this unit is allowed within the cv term
              if (term.units.find(parsed_term.unit_accession) == term.units.end())
              {
                // last chance, a child term of the units was used
                set<String> child_terms;

                bool found_unit(false);
                auto lambda = [&parsed_term, &found_unit] (const String& child)
                {
                  if (child == parsed_term.accession)
                  {
                    found_unit = true;
                    return true;
                  }
                  return false;
                };
                for (set<String>::const_iterator it = term.units.begin(); it != term.units.end(); ++it)
                {
                  if (cv_.iterateAllChildren(*it, lambda)) break;
                }

                if (!found_unit)
                {
                  errors_.push_back(String("Unit CV term not allowed: " + parsed_term.unit_accession + " - " + parsed_term.unit_name + " of term " + parsed_term.accession + " - " + parsed_term.name));
                }
              }
            }
            else
            {
              errors_.push_back(String("Unit CV term not found: " + parsed_term.unit_accession + " - " + parsed_term.unit_name + " of term " + parsed_term.accession + " - " + parsed_term.name));
            }
          }
        }
        else
        {
          // check whether unit was used
          if (parsed_term.has_unit_accession || parsed_term.has_unit_name)
          {
            warnings_.push_back(String("Unit CV term used, but not allowed: " + parsed_term.unit_accession + " - " + parsed_term.unit_name + " of term " + parsed_term.accession + " - " + parsed_term.name));
          }
        }
      }

      if (!rule_found) //No rule found
      {
        warnings_.push_back(String("No mapping rule found for element '") + getPath_(1) + "'");
      }
      else if (!allowed) //if rule found and not allowed
      {
        errors_.push_back(String("CV term used in invalid element: '") + parsed_term.accession + " - " + parsed_term.name + "' at element '" + getPath_(1) + "'");
      }

      //check if term accession and term name match
      if (cv_.exists(parsed_term.accession))
      {
        String parsed_name = parsed_term.name;
        parsed_name.trim();
        String correct_name = cv_.getTerm(parsed_term.accession).name;
        correct_name.trim();
        if (parsed_name != correct_name)
        {
          errors_.push_back(String("Name of CV term not correct: '") + parsed_term.accession + " - " + parsed_name + "' should be '" + correct_name + "'");
        }
      }

      if (check_term_value_types_) //check values
      {
        ControlledVocabulary::CVTerm::XRefType type = cv_.getTerm(parsed_term.accession).xref_type;

        // get value, if it exists
        if (parsed_term.has_value && (!parsed_term.value.empty() || type == ControlledVocabulary::CVTerm::XSD_STRING))
        {
          String value = parsed_term.value;
          if (type == ControlledVocabulary::CVTerm::NONE)
          {
            //Quality CV does not state value type :(
            if (!parsed_term.accession.hasPrefix("PATO:"))
            {
              errors_.push_back(String("Value of CV term not allowed: '") + parsed_term.accession + " - " + parsed_term.name + "' value='" + value + "' at element '" + getPath_(1) + "'");
            }
          }
          else if (type == ControlledVocabulary::CVTerm::XSD_STRING)
          {
            // nothing to check
          }
          else if (type == ControlledVocabulary::CVTerm::XSD_INTEGER)
          {
            try
            {
              parsed_term.value.toInt();
            }
            catch (Exception::ConversionError& /*e*/)
            {
              errors_.push_back(String("Value-type of CV term wrong, should be xsd:integer: '") + parsed_term.accession + " - " + parsed_term.name + "' value='" + value + "' at element '" + getPath_(1) + "'");
            }
          }
          else if (type == ControlledVocabulary::CVTerm::XSD_DECIMAL)
          {
            try
            {
              value.toDouble();
            }
            catch (Exception::ConversionError& /*e*/)
            {
              errors_.push_back(String("Value-type of CV term wrong, should be xsd:decimal: '") + parsed_term.accession + " - " + parsed_term.name + "' value='" + value + "' at element '" + getPath_(1) + "'");
            }
          }
          else if (type == ControlledVocabulary::CVTerm::XSD_NEGATIVE_INTEGER)
          {
            try
            {
              int int_value = value.toInt();
              if (int_value >= 0)
              {
                throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "should be negative");
              }
            }
            catch (Exception::ConversionError& /*e*/)
            {
              errors_.push_back(String("Value-type of CV term wrong, should be xsd:negativeInteger: '") + parsed_term.accession + " - " + parsed_term.name + "' value='" + value + "' at element '" + getPath_(1) + "'");
            }
          }
          else if (type == ControlledVocabulary::CVTerm::XSD_POSITIVE_INTEGER)
          {
            try
            {
              int int_value = value.toInt();
              if (int_value <= 0)
              {
                throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "should be positive");
              }
            }
            catch (Exception::ConversionError& /*e*/)
            {
              errors_.push_back(String("Value-type of CV term wrong, should be xsd:positiveInteger: '") + parsed_term.accession + " - " + parsed_term.name + "' value='" + value + "' at element '" + getPath_(1) + "'");
            }
          }
          else if (type == ControlledVocabulary::CVTerm::XSD_NON_NEGATIVE_INTEGER)
          {
            try
            {
              int int_value = value.toInt();
              if (int_value < 0)
              {
                throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "should not be negative");
              }
            }
            catch (Exception::ConversionError& /*e*/)
            {
              errors_.push_back(String("Value-type of CV term wrong, should be xsd:nonNegativeInteger: '") + parsed_term.accession + " - " + parsed_term.name + "' value='" + value + "' at element '" + getPath_(1) + "'");
            }
          }
          else if (type == ControlledVocabulary::CVTerm::XSD_NON_POSITIVE_INTEGER)
          {
            try
            {
              int int_value = value.toInt();
              if (int_value > 0)
              {
                throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "should not be positive");
              }
            }
            catch (Exception::ConversionError& /*e*/)
            {
              errors_.push_back(String("Value-type of CV term wrong, should be xsd:nonPositiveInteger: '") + parsed_term.accession + " - " + parsed_term.name + "' value='" + value + "' at element '" + getPath_(1) + "'");
            }
          }
          else if (type == ControlledVocabulary::CVTerm::XSD_BOOLEAN)
          {
            String value_copy = value;
            value_copy.trim();
            value_copy.toLower();
            if (value_copy != "1" && value_copy  != "0" && value_copy != "true" && value_copy != "false")
            {
              errors_.push_back(String("Value-type of CV term wrong, should be xsd:boolean: '") + parsed_term.accession + " - " + parsed_term.name + "' value='" + value + "' at element '" + getPath_(1) + "'");
            }
          }
          else if (type == ControlledVocabulary::CVTerm::XSD_DATE)
          {
            try
            {
              DateTime tmp;
              tmp.set(value);
            }
            catch (Exception::ParseError&)
            {
              errors_.push_back(String("Value-type of CV term wrong, should be xsd:date: '") + parsed_term.accession + " - " + parsed_term.name + "' value='" + value + "' at element '" + getPath_(1) + "'");
            }
          }
          else if (type == ControlledVocabulary::CVTerm::XSD_ANYURI)
          { // according to RFC 2396 this is there must be a colon (looked only 2 minutes on it)
            if (!value.has(':'))
            {
              errors_.push_back(String("Value-type of CV term wrong, should be xsd:anyURI (at least a colon is needed): '") + parsed_term.accession + " - " + parsed_term.name + "' value=" + value + "' at element '" + getPath_(1) + "'");
            }
          }
          else
          {
            errors_.push_back(String("Value-type unknown (type #" + String(type) + "): '") + parsed_term.accession + " - " + parsed_term.name + "' value='" + value + "' at element '" + getPath_(1) + "'");
          }
        }
        else if (cv_.getTerm(parsed_term.accession).xref_type != ControlledVocabulary::CVTerm::NONE)
        {
          errors_.push_back(String("Value-type required, but not given (" + ControlledVocabulary::CVTerm::getXRefTypeName(cv_.getTerm(parsed_term.accession).xref_type) + "): '") + parsed_term.accession + " - " + parsed_term.name + "' value='" + parsed_term.value + "' at element '" + getPath_(1) + "'");
        }
      }
    }

    bool SemanticValidator::locateTerm(const String& path, const CVTerm& parsed_term) const
    {
      //check if the term is allowed in this element
      //and if there is a mapping rule for this element
      //Also store fulfilled rule term counts - this count is used to check of the MUST/MAY and AND/OR/XOR is fulfilled
      const vector<CVMappingRule>& rules = rules_.at(path);
      for (Size r = 0; r < rules.size(); ++r) //go thru all rules
      {
        for (Size t = 0; t < rules[r].getCVTerms().size(); ++t) //go thru all terms
        {
          const CVMappingTerm& term = rules[r].getCVTerms()[t];
          if (term.getUseTerm() && term.getAccession() == parsed_term.accession) //check if the term itself is allowed
          {
            return true;
          }
          auto searcher = [&parsed_term] (const String& child) { return child == parsed_term.accession; };
          if (term.getAllowChildren() && cv_.iterateAllChildren(term.getAccession(), searcher))
          {
            return true;
          }
        }
      }
      return false;
    }
    
} // namespace OpenMS // namespace Internal
