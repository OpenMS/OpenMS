// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/VALIDATORS/MzDataValidator.h>
#include <OpenMS/DATASTRUCTURES/CVMappingTerm.h>
#include <OpenMS/DATASTRUCTURES/CVMappingRule.h>

#include <OpenMS/FORMAT/ControlledVocabulary.h>

using namespace xercesc;
using namespace std;

namespace OpenMS::Internal
{

    MzDataValidator::MzDataValidator(const CVMappings& mapping, const ControlledVocabulary& cv) :
      SemanticValidator(mapping, cv)
    {
      setCheckUnits(true);
    }

    MzDataValidator::~MzDataValidator()
    = default;

    void MzDataValidator::handleTerm_(const String& path, const CVTerm& parsed_term)
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
            ++fulfilled_[path][rules[r].getIdentifier()][term.getAccession()];
            break;
          }
          if (term.getAllowChildren()) //check if the term's children are allowed
          {
            auto searcher = [&parsed_term] (const String& child)
            {
              return child == parsed_term.accession;
            };

            if (cv_.iterateAllChildren(term.getAccession(), searcher))
            {
              allowed = true;
              ++fulfilled_[path][rules[r].getIdentifier()][term.getAccession()];
              break;
            }
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
                auto lambda = [&parsed_term] (const String& child)
                {
                  return child == parsed_term.unit_accession;
                };

                bool found_unit(false);
                for (set<String>::const_iterator it = term.units.begin(); it != term.units.end(); ++it)
                {
                  if (cv_.iterateAllChildren(*it, lambda))
                  {
                    found_unit = true;
                    break;
                  }
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

        //be a bit more soft: ignore upper-lower case
        parsed_name.toLower();
        correct_name.toLower();

        //be a bit more soft: ignore spaces
        parsed_name.removeWhitespaces();
        correct_name.removeWhitespaces();

        if (parsed_name != correct_name)
        {
          errors_.push_back(String("Name of CV term not correct: '") + parsed_term.accession + " - " + parsed_name + "' should be '" + correct_name + "'");
        }
      }
    }

} // namespace OpenMS // namespace Internal
